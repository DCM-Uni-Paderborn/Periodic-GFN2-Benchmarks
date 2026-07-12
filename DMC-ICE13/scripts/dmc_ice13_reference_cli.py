#!/usr/bin/env python3
"""Compare CP2K-native Gamma results with the tblite command-line driver."""

from __future__ import annotations

import argparse
import concurrent.futures
import os
import sys
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
REPOSITORY = ROOT.parent
sys.path.insert(0, str(REPOSITORY / "scripts"))

from run_x23b_reference_cli_checks import (  # noqa: E402
    inject_reference_cli,
    parse_reference_cli,
    prepare_reference_cli_program,
    summarize,
    write_csv,
)


PHASES = [
    "Ih",
    "II",
    "III",
    "IV",
    "VI",
    "VII",
    "VIII",
    "IX",
    "XI",
    "XIII",
    "XIV",
    "XV",
    "XVII",
]
METHODS = ["GFN1", "GFN2"]


def run_case(method: str, phase: str, args: argparse.Namespace) -> dict[str, str]:
    source = ROOT / "inputs" / f"ice_{phase}_{method}.inp"
    run_dir = args.out / method / phase
    run_dir.mkdir(parents=True, exist_ok=True)
    prefix = f"refcli_dmc13_{method}_{phase}".replace("-", "_")
    inp = run_dir / f"ice_{phase}_{method}_reference_cli.inp"
    program = prepare_reference_cli_program(run_dir, args.tblite)
    inp.write_text(inject_reference_cli(source.read_text(), program, prefix, args.keep_files))
    out_file = run_dir / "cp2k.out"

    env = os.environ.copy()
    env.setdefault("OMP_NUM_THREADS", "1")
    env.setdefault("OPENBLAS_NUM_THREADS", "1")
    env.setdefault("MKL_NUM_THREADS", "1")
    with out_file.open("w") as handle:
        proc = subprocess.run(
            [str(args.cp2k), "-i", inp.name],
            cwd=run_dir,
            stdout=handle,
            stderr=subprocess.STDOUT,
            env=env,
            text=True,
        )

    row = {
        "source_kind": "dmc13_gamma",
        "method": method,
        "system": phase,
        "source": source.relative_to(REPOSITORY).as_posix(),
        "run_dir": run_dir.relative_to(args.out).as_posix(),
        "returncode": str(proc.returncode),
    }
    row.update(parse_reference_cli(out_file))
    return row


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cp2k", type=Path, required=True)
    parser.add_argument("--tblite", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--jobs", type=int, default=4)
    parser.add_argument("--keep-files", action="store_true")
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    rows: list[dict[str, str]] = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as pool:
        futures = {
            pool.submit(run_case, method, phase, args): (method, phase)
            for method in METHODS
            for phase in PHASES
        }
        for future in concurrent.futures.as_completed(futures):
            row = future.result()
            rows.append(row)
            print(
                row["method"],
                row["system"],
                "rc",
                row["returncode"],
                "gmax",
                row["gradient_diff_max"],
                "vmax",
                row["virial_diff_max"],
                flush=True,
            )

    rows.sort(key=lambda row: (row["method"], row["system"]))
    columns = [
        "source_kind",
        "method",
        "system",
        "returncode",
        "energy_cp2k_hartree",
        "energy_cli_hartree",
        "energy_absdiff_hartree",
        "gradient_diff_sum",
        "gradient_diff_max",
        "virial_diff_sum",
        "virial_diff_max",
        "exceeded_error_limit",
        "skipped",
        "source",
        "run_dir",
    ]
    write_csv(args.out / "reference_cli_rows.csv", rows, columns)
    summary = summarize(rows)
    write_csv(
        args.out / "reference_cli_summary.csv",
        summary,
        list(summary[0].keys()) if summary else ["source_kind", "method", "n"],
    )


if __name__ == "__main__":
    main()
