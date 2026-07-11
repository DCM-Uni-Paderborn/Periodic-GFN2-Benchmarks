#!/usr/bin/env python3
"""Follow a converged Klimes-23 SCC solution through nearby cell volumes."""

from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import time
from pathlib import Path

import numpy as np

from run_literature_eos_benchmarks import (
    CP2K_DEFAULT,
    Job,
    cp2k_input,
    load_structures,
    parse_cp2k_output,
    scale_structure,
)


REPO = Path(__file__).resolve().parents[1]


def parse_factors(value: str) -> list[float]:
    factors = [float(field) for field in value.split(",")]
    if not factors or any(factor <= 0 for factor in factors):
        raise argparse.ArgumentTypeError("volume factors must be positive")
    return factors


def restart_file(root: Path) -> Path | None:
    files = list(root.glob("*-RESTART.kp"))
    return max(files, key=lambda path: path.stat().st_mtime) if files else None


def result_file(root: Path) -> Path:
    return root / "result.json"


def load_completed_step(root: Path) -> tuple[dict[str, object], Path] | None:
    path = result_file(root)
    restart = restart_file(root)
    if not path.is_file() or restart is None:
        return None
    try:
        result = json.loads(path.read_text())
    except json.JSONDecodeError:
        return None
    if not result.get("ok"):
        return None
    return result, restart


def publish_step(job: Job, source: Path, result: dict[str, object]) -> None:
    job.root.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source / "cp2k.inp", job.root / "cp2k.inp")
    shutil.copy2(source / "cp2k.out", job.root / "cp2k.out")
    result_file(job.root).write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")


def run_step(
    args: argparse.Namespace,
    structure,
    factor: float,
    previous_factor: float,
    previous_restart: Path,
    index: int,
) -> tuple[dict[str, object], Path]:
    root = REPO / "Klimes-Solids23/results/continuation" / args.method / args.system / f"v_{factor:.4f}"
    root.mkdir(parents=True, exist_ok=True)
    completed = load_completed_step(root)
    if completed is not None:
        print(f"CONTINUE SKIP {args.method} {args.system} v={factor:.4f}", flush=True)
        return completed

    job = Job("klimes", structure, args.method, factor, args.cp2k, args.threads, True)
    text = cp2k_input(job, args.damping, args.max_scf, args.memory, "TBLITE")
    text = text.replace(
        "  &DFT\n",
        f"  &DFT\n    WFN_RESTART_FILE_NAME {previous_restart.resolve()}\n",
        1,
    )
    text = text.replace("SCF_GUESS MOPAC", "SCF_GUESS RESTART")
    text = text.replace("&RESTART OFF", "&RESTART ON")
    (root / "cp2k.inp").write_text(text)
    output = root / "cp2k.out"
    if output.exists():
        output.unlink()
    env = os.environ.copy()
    env.update(
        {
            "OMP_NUM_THREADS": str(args.threads),
            "OMP_PROC_BIND": "false",
            "OPENBLAS_NUM_THREADS": "1",
            "MKL_NUM_THREADS": "1",
            "VECLIB_MAXIMUM_THREADS": "1",
        }
    )
    started = time.time()
    with (root / "launcher.log").open("w") as launcher:
        proc = subprocess.run(
            [str(args.cp2k), "-i", "cp2k.inp", "-o", "cp2k.out"],
            cwd=root,
            env=env,
            stdout=launcher,
            stderr=subprocess.STDOUT,
            check=False,
        )
    ok, energy, reason = parse_cp2k_output(output)
    current_restart = restart_file(root)
    ok = bool(ok and proc.returncode == 0 and current_restart is not None)
    if current_restart is None and reason == "ok":
        reason = "CP2K did not write a k-point restart"
    cell, _ = scale_structure(structure, factor)
    result: dict[str, object] = {
        "ok": ok,
        "reason": "ok" if ok else reason,
        "return_code": proc.returncode,
        "dataset": "klimes",
        "system": structure.system,
        "label": structure.label,
        "category": structure.category,
        "method": args.method,
        "volume_factor": factor,
        "volume_A3": abs(float(np.linalg.det(cell))),
        "energy_hartree": energy,
        "natoms": structure.natoms,
        "kmesh": list(structure.kmesh),
        "kpoint_path": "CP2K native Bloch",
        "full_grid": True,
        "initial_guess": "volume_continuation",
        "seed_volume_factor": previous_factor,
        "restart_source": str(previous_restart.resolve()),
        "continuation_index": index,
        "scc_mixer": "TBLITE",
        "damping": args.damping,
        "tblite_mixer_damping": args.damping,
        "tblite_mixer_iterations": args.max_scf,
        "tblite_mixer_memory": args.memory,
        "elapsed_seconds": time.time() - started,
    }
    result_file(root).write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    state = "OK" if ok else "FAIL"
    print(f"CONTINUE {state} {args.method} {args.system} v={factor:.4f} ({result['elapsed_seconds']:.1f}s)", flush=True)
    if not ok or current_restart is None:
        raise RuntimeError(f"continuation failed at v={factor:.4f}: {reason}")
    return result, current_restart


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--system", required=True)
    parser.add_argument("--method", required=True, choices=("GFN1", "GFN2"))
    parser.add_argument("--factors", required=True, type=parse_factors)
    parser.add_argument("--seed-factor", required=True, type=float)
    parser.add_argument("--seed-restart", required=True, type=Path)
    parser.add_argument("--publish-factors", type=parse_factors, default=[])
    parser.add_argument("--cp2k", type=Path, default=CP2K_DEFAULT)
    parser.add_argument("--damping", type=float, default=0.20)
    parser.add_argument("--max-scf", type=int, default=750)
    parser.add_argument("--memory", type=int, default=100)
    parser.add_argument("--threads", type=int, default=1)
    args = parser.parse_args()

    if not args.cp2k.is_file():
        parser.error(f"CP2K binary not found: {args.cp2k}")
    if not args.seed_restart.is_file():
        parser.error(f"seed restart not found: {args.seed_restart}")
    structures = [structure for structure in load_structures("klimes") if structure.system == args.system]
    if len(structures) != 1:
        parser.error(f"unknown Klimes-23 system: {args.system}")
    publish = {round(factor, 8) for factor in args.publish_factors}
    previous_factor = args.seed_factor
    previous_restart = args.seed_restart
    for index, factor in enumerate(args.factors, start=1):
        result, previous_restart = run_step(
            args,
            structures[0],
            factor,
            previous_factor,
            previous_restart,
            index,
        )
        if round(factor, 8) in publish:
            job = Job("klimes", structures[0], args.method, factor, args.cp2k, args.threads, True)
            source = REPO / "Klimes-Solids23/results/continuation" / args.method / args.system / f"v_{factor:.4f}"
            publish_step(job, source, result)
            print(f"PUBLISH {args.method} {args.system} v={factor:.4f}", flush=True)
        previous_factor = factor
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
