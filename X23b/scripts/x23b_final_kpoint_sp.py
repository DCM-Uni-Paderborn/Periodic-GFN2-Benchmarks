#!/usr/bin/env python3
"""Run native-Bloch single points on converged X23b k222 cell-opt geometries."""

from __future__ import annotations

import argparse
import csv
import fcntl
import json
import math
import os
import re
import shutil
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import x23b_kpoint_cellopt as cellopt


FIELDS = (
    "method",
    "system",
    "target_mesh",
    "source_mesh",
    "program_ended",
    "source_energy_hartree",
    "target_energy_hartree",
    "source_lattice_energy_kJmol",
    "target_lattice_energy_kJmol",
    "delta_target_minus_source_kJmol",
    "target_error_kJmol",
    "source_restart",
    "output",
)


def variant(mesh: int) -> str:
    return f"k{mesh}{mesh}{mesh}_sp_on_k222"


def gamma_centered_shift(mesh: int) -> float:
    return 0.0 if mesh % 2 else (mesh - 1) / (2.0 * mesh)


def final_restart(run_dir: Path) -> Path:
    candidates = list(run_dir.glob("*-1.restart"))
    if not candidates:
        raise ValueError(f"final restart not found in {run_dir}")
    return max(candidates, key=lambda path: path.stat().st_mtime_ns)


def restart_to_single_point(source: Path, project: str, mesh: int) -> str:
    lines = source.read_text().splitlines()
    global_start = next(index for index, line in enumerate(lines) if line.strip().upper() == "&GLOBAL")
    lines = lines[global_start:]

    motion_start, motion_end = cellopt.section_bounds(lines, "MOTION")
    del lines[motion_start : motion_end + 1]

    replacements = (
        (re.compile(r"^(\s*PROJECT(?:_NAME)?\s+).*$", re.I), rf'\1"{project}"'),
        (re.compile(r"^(\s*RUN_TYPE\s+).*$", re.I), r"\1ENERGY"),
    )
    for pattern, replacement in replacements:
        for index, line in enumerate(lines):
            if pattern.match(line):
                lines[index] = pattern.sub(replacement, line, count=1)
                break
        else:
            raise ValueError(f"required GLOBAL keyword missing in {source}: {pattern.pattern}")

    k_start, k_end = cellopt.section_bounds(lines, "KPOINTS")
    shift = gamma_centered_shift(mesh)
    for index in range(k_start, k_end + 1):
        if lines[index].strip().upper().startswith("SCHEME"):
            indent = re.match(r"\s*", lines[index]).group(0)
            lines[index] = (
                f"{indent}SCHEME MACDONALD {mesh} {mesh} {mesh} "
                f"{shift:.12g} {shift:.12g} {shift:.12g}"
            )
            break
    else:
        raise ValueError(f"KPOINTS SCHEME missing in {source}")
    return "\n".join(lines) + "\n"


def manifest_path(output_root: Path) -> Path:
    return output_root / "manifest.csv"


def prepare(args: argparse.Namespace) -> None:
    output_root = args.output_root.resolve()
    output_root.parent.mkdir(parents=True, exist_ok=True)
    lock_path = output_root.parent / f".{output_root.name}.manifest.lock"
    with lock_path.open("w") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        _prepare(args, output_root)


def _prepare(args: argparse.Namespace, output_root: Path) -> None:
    if args.clean and output_root.exists():
        shutil.rmtree(output_root)
    methods = (args.method,) if args.method else cellopt.METHODS
    selected_systems = set(args.system) if args.system else None
    manifest: dict[tuple[str, str], dict[str, str]] = {}
    path = manifest_path(output_root)
    if path.is_file() and not args.clean:
        with path.open(newline="") as handle:
            for row in csv.DictReader(handle):
                row.setdefault("mesh", "3")
                manifest[(row["method"], row["system"])] = row

    prepared = 0
    for system_data in cellopt.systems():
        system = str(system_data["id"])
        if selected_systems is not None and system not in selected_systems:
            continue
        for method in methods:
            source_run_dir = args.cellopt_root.resolve() / method / system / cellopt.VARIANT
            source_output = source_run_dir / "cp2k.out"
            if not cellopt.cp2k_completed(source_output):
                raise ValueError(f"converged k222 optimization required: {method}/{system}")
            source_restart = final_restart(source_run_dir)
            target_variant = variant(args.mesh)
            project = f"{system}_{method}_{target_variant}".replace("-", "_")
            run_dir = output_root / method / system / target_variant
            run_dir.mkdir(parents=True, exist_ok=True)
            input_path = run_dir / f"{project}.inp"
            input_path.write_text(restart_to_single_point(source_restart, project, args.mesh))
            manifest[(method, system)] = {
                "method": method,
                "system": system,
                "mesh": str(args.mesh),
                "source_run_dir": str(source_run_dir),
                "source_restart": str(source_restart),
                "input": str(input_path),
                "run_dir": str(run_dir),
            }
            prepared += 1

    rows = [manifest[key] for key in sorted(manifest)]
    output_root.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=("method", "system", "mesh", "source_run_dir", "source_restart", "input", "run_dir"),
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)
    print(f"Prepared {prepared} k{args.mesh}{args.mesh}{args.mesh} single points; manifest has {len(rows)} entries")


def output_ok(path: Path) -> bool:
    if not path.is_file():
        return False
    text = path.read_text(errors="ignore")
    return "PROGRAM ENDED" in text and "ENERGY| Total FORCE_EVAL" in text


def run_one(input_path: Path, cp2k: Path, threads: int, force: bool) -> tuple[Path, int, str]:
    run_dir = input_path.parent
    output = run_dir / "cp2k.out"
    with input_path.open() as input_lock:
        try:
            fcntl.flock(input_lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            return input_path, 0, "BUSY"
        if not force and output_ok(output):
            return input_path, 0, "SKIP"
        if force:
            for path in run_dir.iterdir():
                if path != input_path:
                    path.unlink() if path.is_file() else shutil.rmtree(path)
        env = os.environ.copy()
        env.update(
            {
                "OMP_NUM_THREADS": str(threads),
                "OPENBLAS_NUM_THREADS": str(threads),
                "MKL_NUM_THREADS": str(threads),
            }
        )
        code = subprocess.run(
            [str(cp2k), "-i", input_path.name, "-o", output.name],
            cwd=run_dir,
            env=env,
            check=False,
        ).returncode
        action = "CONVERGED" if code == 0 and output_ok(output) else "FAILED"
        return input_path, code, action


def run(args: argparse.Namespace) -> None:
    inputs = sorted(args.output_root.resolve().glob(f"*/**/{variant(args.mesh)}/*.inp"))
    if args.method:
        inputs = [path for path in inputs if args.method in path.parts]
    if args.system:
        wanted = set(args.system)
        inputs = [path for path in inputs if any(system in path.parts for system in wanted)]
    if not inputs:
        raise ValueError("no prepared single-point inputs selected")
    failed = []
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        futures = {
            pool.submit(run_one, path, args.cp2k.resolve(), args.threads_per_job, args.force): path
            for path in inputs
        }
        for future in as_completed(futures):
            input_path, code, action = future.result()
            relative = input_path.relative_to(args.output_root.resolve())
            print(f"{action:9s} {relative} rc={code}", flush=True)
            if code != 0 or action == "FAILED":
                failed.append(relative)
    if failed:
        raise SystemExit(f"{len(failed)} single-point jobs failed")


def finite(value: float | None) -> str:
    return "" if value is None or not math.isfinite(value) else f"{value:.12f}"


def collect(args: argparse.Namespace) -> None:
    output_root = args.output_root.resolve()
    with manifest_path(output_root).open(newline="") as handle:
        manifest = list(csv.DictReader(handle))
    metadata = {str(row["id"]): row for row in cellopt.systems()}
    gas = cellopt.load_molecule_rows(args.molecule_run_root.resolve())
    rows = []
    for entry in manifest:
        method, system = entry["method"], entry["system"]
        target_mesh = int(entry["mesh"])
        k222_output = Path(entry["source_run_dir"]) / "cp2k.out"
        k333_output = Path(entry["run_dir"]) / "cp2k.out"
        k222_text = k222_output.read_text(errors="ignore")
        k333_text = k333_output.read_text(errors="ignore") if k333_output.is_file() else ""
        pattern = r"^\s*ENERGY\| Total FORCE_EVAL .*?([-+0-9.Ee]+)\s*$"
        e222 = cellopt.last_float(k222_text, pattern)
        e333 = cellopt.last_float(k333_text, pattern)
        n_molecules = int(metadata[system]["molecules_per_cell"])
        gas_energy = float(gas[(method, system)]["gas_energy_hartree"])
        ref = float(metadata[system]["ref_energy"])
        lattice222 = None if e222 is None else (gas_energy - e222 / n_molecules) * cellopt.HARTREE_TO_KJMOL
        lattice333 = None if e333 is None else (gas_energy - e333 / n_molecules) * cellopt.HARTREE_TO_KJMOL
        rows.append(
            {
                "method": method,
                "system": system,
                "target_mesh": f"k{target_mesh}{target_mesh}{target_mesh}",
                "source_mesh": "k222_cellopt",
                "program_ended": "PROGRAM ENDED" in k333_text,
                "source_energy_hartree": finite(e222),
                "target_energy_hartree": finite(e333),
                "source_lattice_energy_kJmol": finite(lattice222),
                "target_lattice_energy_kJmol": finite(lattice333),
                "delta_target_minus_source_kJmol": finite(
                    None if lattice222 is None or lattice333 is None else lattice333 - lattice222
                ),
                "target_error_kJmol": finite(None if lattice333 is None else lattice333 - ref),
                "source_restart": entry["source_restart"],
                "output": str(k333_output),
            }
        )
    args.csv.parent.mkdir(parents=True, exist_ok=True)
    with args.csv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    for target_mesh in sorted({row["target_mesh"] for row in rows}):
        for method in cellopt.METHODS:
            deltas = [
                abs(float(row["delta_target_minus_source_kJmol"]))
                for row in rows
                if row["method"] == method
                and row["target_mesh"] == target_mesh
                and row["delta_target_minus_source_kJmol"]
            ]
            if deltas:
                print(
                    f"{method} {target_mesh}: {len(deltas)} complete, "
                    f"mean |target-k222|={sum(deltas) / len(deltas):.6f}, max={max(deltas):.6f} kJ/mol"
                )
            else:
                print(f"{method} {target_mesh}: 0 complete")


def main() -> None:
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="command", required=True)
    choices = sorted(str(row["id"]) for row in cellopt.systems())

    prepare_parser = subparsers.add_parser("prepare")
    prepare_parser.add_argument("--cellopt-root", type=Path, required=True)
    prepare_parser.add_argument("--output-root", type=Path, required=True)
    prepare_parser.add_argument("--mesh", type=int, default=3)
    prepare_parser.add_argument("--method", choices=cellopt.METHODS)
    prepare_parser.add_argument("--system", action="append", choices=choices)
    prepare_parser.add_argument("--clean", action="store_true")
    prepare_parser.set_defaults(function=prepare)

    run_parser = subparsers.add_parser("run")
    run_parser.add_argument("--output-root", type=Path, required=True)
    run_parser.add_argument("--cp2k", type=Path, required=True)
    run_parser.add_argument("--jobs", type=int, default=4)
    run_parser.add_argument("--threads-per-job", type=int, default=1)
    run_parser.add_argument("--mesh", type=int, default=3)
    run_parser.add_argument("--method", choices=cellopt.METHODS)
    run_parser.add_argument("--system", action="append", choices=choices)
    run_parser.add_argument("--force", action="store_true")
    run_parser.set_defaults(function=run)

    collect_parser = subparsers.add_parser("collect")
    collect_parser.add_argument("--output-root", type=Path, required=True)
    collect_parser.add_argument("--molecule-run-root", type=Path, required=True)
    collect_parser.add_argument("--csv", type=Path, required=True)
    collect_parser.set_defaults(function=collect)

    args = parser.parse_args()
    args.function(args)


if __name__ == "__main__":
    main()
