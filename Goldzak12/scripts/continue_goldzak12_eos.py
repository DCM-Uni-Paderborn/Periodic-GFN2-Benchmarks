#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import os
import shutil
import zipfile
from pathlib import Path

import numpy as np

import run_goldzak12_benchmark as base
import run_goldzak12_eos_benchmark as eos


SCC_RESTART_KEYS = ("tblite0_qat", "tblite0_qsh", "tblite0_dpat", "tblite0_qpat")


def extrapolate_scc_restart(
    previous: Path,
    current: Path,
    previous_scale: float,
    current_scale: float,
    target_scale: float,
    output: Path,
) -> None:
    denominator = current_scale - previous_scale
    if abs(denominator) < 1.0e-12:
        raise ValueError("Cannot extrapolate SCC restart from identical scales")
    factor = (target_scale - current_scale) / denominator
    with (
        np.load(previous) as old_data,
        np.load(current) as current_data,
        zipfile.ZipFile(current) as current_archive,
    ):
        payloads: dict[str, bytes] = {}
        for key in current_data.files:
            payload = current_archive.read(f"{key}.npy")
            if key in SCC_RESTART_KEYS:
                value = np.array(current_data[key], order="F", copy=True)
                value += factor * (current_data[key] - old_data[key])
                major = payload[6]
                header_size = 2 if major == 1 else 4
                header_length = int.from_bytes(payload[8 : 8 + header_size], "little")
                data_offset = 8 + header_size + header_length
                raw = np.asfortranarray(value).tobytes(order="F")
                if len(raw) != len(payload) - data_offset:
                    raise ValueError(f"Unexpected NPY payload size for {key}")
                payload = payload[:data_offset] + raw
            payloads[key] = payload
    output.parent.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(output, mode="w", compression=zipfile.ZIP_STORED, allowZip64=False) as archive:
        for key, payload in payloads.items():
            info = zipfile.ZipInfo(f"{key}.npy", date_time=(1980, 1, 1, 0, 0, 0))
            info.compress_type = zipfile.ZIP_STORED
            archive.writestr(info, payload)


def continuation_scale_tag(scale: float) -> str:
    return f"s{scale:.5f}".replace(".", "p")


def continuation_input(
    ref: base.Reference,
    method: str,
    mesh: str,
    scale: float,
    project: str,
    restart: Path | None,
    scc_restart: Path | None,
    scc_restart_out: Path,
    mixer: str,
    max_scf: int,
    damping: float,
    memory: int,
    alpha: float,
    disable_diis: bool,
    omega0: float,
    min_weight: float,
    max_weight: float,
    weight_factor: float,
) -> str:
    text = base.solid_input(ref, method, "ENERGY", mesh, ref.a_exp * scale, project)
    text = text.replace("        &RESTART OFF", "        &RESTART ON", 1)
    text = text.replace("      MAX_SCF 300", f"      MAX_SCF {max_scf}", 1)
    text = text.replace("        ALPHA 0.2", f"        ALPHA {alpha:.8f}", 1)
    scc_lines = f"          SCC_RESTART_WRITE_FILE_NAME {scc_restart_out.name}\n"
    if scc_restart is not None:
        restart_name = os.path.relpath(scc_restart, start=scc_restart_out.parent)
        scc_lines = f"          SCC_RESTART_FILE_NAME {restart_name}\n" + scc_lines
    text = text.replace(
        "          ACCURACY 0.05\n        &END TBLITE",
        "          ACCURACY 0.05\n" + scc_lines + "        &END TBLITE",
        1,
    )
    if mixer == "tblite":
        marker = "        &END TBLITE\n      &END XTB"
        replacement = (
            "        &END TBLITE\n"
            "        SCC_MIXER TBLITE\n"
            "        &TBLITE_MIXER\n"
            f"          ITERATIONS {max_scf}\n"
            f"          MEMORY {memory}\n"
            f"          DAMPING {damping:.8f}\n"
            f"          OMEGA0 {omega0:.8f}\n"
            f"          MIN_WEIGHT {min_weight:.8f}\n"
            f"          MAX_WEIGHT {max_weight:.8f}\n"
            f"          WEIGHT_FACTOR {weight_factor:.8f}\n"
            "        &END TBLITE_MIXER\n"
            "      &END XTB"
        )
        text = text.replace(marker, replacement, 1)
    elif mixer == "cp2k":
        text = text.replace("        GFN_TYPE TBLITE", "        GFN_TYPE TBLITE\n        SCC_MIXER CP2K", 1)
        if disable_diis:
            text = text.replace("      EPS_SCF 1.0E-9", "      EPS_SCF 1.0E-9\n      EPS_DIIS 1.0E-10", 1)
    if restart is not None:
        text = text.replace(
            "  &DFT\n",
            f"  &DFT\n    WFN_RESTART_FILE_NAME {restart.resolve()}\n",
            1,
        )
        text = text.replace("      SCF_GUESS MOPAC", "      SCF_GUESS RESTART", 1)
    return text


def run_point(
    ref: base.Reference,
    method: str,
    mesh: str,
    scale: float,
    restart: Path | None,
    scc_restart: Path | None,
    previous_scc_restart: Path | None,
    previous_scale: float | None,
    restart_scale: float,
    extrapolate_scc: bool,
    cp2k: Path,
    threads: int,
    variant: str,
    mixer: str,
    max_scf: int,
    damping: float,
    memory: int,
    alpha: float,
    disable_diis: bool,
    omega0: float,
    min_weight: float,
    max_weight: float,
    weight_factor: float,
) -> tuple[Path, Path, Path, Path]:
    tag = continuation_scale_tag(scale)
    project = f"{ref.solid}_{method}_eos_cont_{variant}_{mesh}_{tag}"
    run_dir = base.ROOT / "runs" / f"continuation_{variant}" / method / ref.solid / mesh / tag
    inp = run_dir / f"{project}.inp"
    out = run_dir / f"{project}.out"
    scc_restart_out = run_dir / f"{project}-SCC-RESTART.npz"
    scc_restart_input = scc_restart
    if extrapolate_scc and scc_restart is not None and previous_scc_restart is not None:
        if previous_scale is None:
            raise ValueError("previous_scale is required for SCC extrapolation")
        scc_restart_input = run_dir / f"{project}-SCC-PREDICTED.npz"
        extrapolate_scc_restart(
            previous_scc_restart,
            scc_restart,
            previous_scale,
            restart_scale,
            scale,
            scc_restart_input,
        )
    base.write_file(
        inp,
        continuation_input(
            ref,
            method,
            mesh,
            scale,
            project,
            restart,
            scc_restart_input,
            scc_restart_out,
            mixer,
            max_scf,
            damping,
            memory,
            alpha,
            disable_diis,
            omega0,
            min_weight,
            max_weight,
            weight_factor,
        ),
    )
    code = base.run_cp2k(cp2k, inp, out, threads)
    restart_out = run_dir / f"{project}-RESTART.kp"
    if code != 0 or not base.output_ok(out) or not restart_out.exists() or not scc_restart_out.exists():
        raise RuntimeError(f"Continuation failed at scale {scale:.5f}; inspect {out}")
    energy = base.parse_energy(out)
    print(f"ok {ref.solid} {method} {mesh} scale={scale:.5f} energy={energy:.12f}", flush=True)
    return inp, out, restart_out, scc_restart_out


def promote(
    ref: base.Reference,
    method: str,
    mesh: str,
    scale: float,
    predecessor_scale: float,
    inp: Path,
    out: Path,
) -> None:
    project = eos.eos_project(ref.solid, method, mesh, scale)
    run_dir = base.ROOT / "runs" / "eos" / method / ref.solid / mesh / eos.scale_tag(scale)
    target_out = run_dir / f"{project}.out"
    run_dir.mkdir(parents=True, exist_ok=True)
    shutil.copy2(out, target_out)
    eos.strategy_path(target_out).write_text(
        json.dumps(
            {
                "strategy": "volume_continuation",
                "completed": True,
                "predecessor_scale": predecessor_scale,
                "continuation_input": str(inp.relative_to(base.ROOT)),
            },
            indent=2,
        )
        + "\n"
    )


def main() -> int:
    parser = argparse.ArgumentParser(description="Continue an LC10 EOS branch through neighboring volumes.")
    parser.add_argument("--solid", required=True)
    parser.add_argument("--method", choices=base.METHODS, required=True)
    parser.add_argument("--mesh", default="k444")
    parser.add_argument("--start-scale", type=float, required=True)
    parser.add_argument("--start-restart", type=Path)
    parser.add_argument("--start-scc-restart", type=Path)
    parser.add_argument("--previous-scc-restart", type=Path)
    parser.add_argument("--previous-scale", type=float)
    parser.add_argument("--scale", type=float, action="append", required=True)
    parser.add_argument("--cp2k", type=Path, default=base.DEFAULT_CP2K)
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--variant", default="tblite_default")
    parser.add_argument("--mixer", choices=("tblite", "cp2k"), default="tblite")
    parser.add_argument("--max-scf", type=int, default=1200)
    parser.add_argument("--damping", type=float, default=0.4)
    parser.add_argument("--memory", type=int, default=6)
    parser.add_argument("--alpha", type=float, default=0.2)
    parser.add_argument("--disable-diis", action="store_true")
    parser.add_argument("--omega0", type=float, default=0.01)
    parser.add_argument("--min-weight", type=float, default=1.0)
    parser.add_argument("--max-weight", type=float, default=100000.0)
    parser.add_argument("--weight-factor", type=float, default=0.01)
    parser.add_argument("--extrapolate-scc", action="store_true")
    parser.add_argument("--promote", action="store_true")
    args = parser.parse_args()

    refs = {ref.solid: ref for ref in base.REFERENCES}
    if args.solid not in refs:
        parser.error(f"Unknown solid {args.solid!r}")
    ref = refs[args.solid]

    start_scc_restart = args.start_scc_restart.resolve() if args.start_scc_restart else None
    if start_scc_restart is not None and not start_scc_restart.exists():
        parser.error(f"SCC restart file does not exist: {start_scc_restart}")

    if args.start_restart is None:
        _, _, restart, scc_restart = run_point(
            ref,
            args.method,
            args.mesh,
            args.start_scale,
            None,
            start_scc_restart,
            None,
            None,
            args.start_scale,
            False,
            args.cp2k,
            args.threads,
            args.variant,
            args.mixer,
            args.max_scf,
            args.damping,
            args.memory,
            args.alpha,
            args.disable_diis,
            args.omega0,
            args.min_weight,
            args.max_weight,
            args.weight_factor,
        )
    else:
        restart = args.start_restart.resolve()
        if not restart.exists():
            parser.error(f"Restart file does not exist: {restart}")
        scc_restart = start_scc_restart
    previous_scc_restart = args.previous_scc_restart.resolve() if args.previous_scc_restart else None
    if previous_scc_restart is not None and not previous_scc_restart.exists():
        parser.error(f"Previous SCC restart file does not exist: {previous_scc_restart}")
    if previous_scc_restart is not None and args.previous_scale is None:
        parser.error("--previous-scale is required with --previous-scc-restart")

    predecessor = args.start_scale
    previous_scale = args.previous_scale
    for scale in args.scale:
        current_scc_restart = scc_restart
        inp, out, restart, scc_restart = run_point(
            ref,
            args.method,
            args.mesh,
            scale,
            restart,
            scc_restart,
            previous_scc_restart,
            previous_scale,
            predecessor,
            args.extrapolate_scc,
            args.cp2k,
            args.threads,
            args.variant,
            args.mixer,
            args.max_scf,
            args.damping,
            args.memory,
            args.alpha,
            args.disable_diis,
            args.omega0,
            args.min_weight,
            args.max_weight,
            args.weight_factor,
        )
        if args.promote:
            promote(ref, args.method, args.mesh, scale, predecessor, inp, out)
        previous_scc_restart = current_scc_restart
        previous_scale = predecessor
        predecessor = scale
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
