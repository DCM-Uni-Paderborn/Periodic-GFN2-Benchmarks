#!/usr/bin/env python3
from __future__ import annotations

import argparse
import concurrent.futures
import csv
import json
import math
import os
import re
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path


HARTREE_TO_KJMOL = 2625.499638
FLOAT = r"[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[Ee][-+]?\d+)?"
REPOSITORY = Path(__file__).resolve().parents[1]
LEGACY_METHODS = {"GFN1", "GFN2"}


KPOINTS = {
    "gamma": "",
    "k222": "\n".join(
        [
            "    &KPOINTS",
            "      SCHEME MACDONALD 2 2 2 0.25 0.25 0.25",
            "      SYMMETRY T",
            "      FULL_GRID F",
            "      SYMMETRY_BACKEND SPGLIB",
            "      SYMMETRY_REDUCTION_METHOD SPGLIB",
            "    &END KPOINTS",
        ]
    ),
    "k333": "\n".join(
        [
            "    &KPOINTS",
            "      SCHEME MACDONALD 3 3 3 0.0 0.0 0.0",
            "      SYMMETRY T",
            "      FULL_GRID F",
            "      SYMMETRY_BACKEND SPGLIB",
            "      SYMMETRY_REDUCTION_METHOD SPGLIB",
            "    &END KPOINTS",
        ]
    ),
}


@dataclass(frozen=True)
class Case:
    method: str
    system: str
    source_restart: Path


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, str]], columns: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)


def parse_energy(output: Path) -> float | None:
    energy = None
    if not output.exists():
        return energy
    for line in output.read_text(errors="ignore").splitlines():
        if "ENERGY| Total FORCE_EVAL" in line:
            energy = float(line.split()[-1])
    return energy


def parse_last_float(output: Path, label: str) -> float | None:
    if not output.exists():
        return None
    value = None
    pattern = re.compile(rf"{re.escape(label)}:?\s+({FLOAT})", re.I)
    for match in pattern.finditer(output.read_text(errors="ignore")):
        value = float(match.group(1))
    return value


def completed_cp2k(output: Path) -> bool:
    if not output.exists():
        return False
    text = output.read_text(errors="ignore")
    if "ABORT" in text or re.search(r"SCF.*NOT|NOT.*SCF|DID NOT CONVERGE|convergence failure", text, re.I):
        return False
    return "PROGRAM ENDED" in text


def metadata_by_system(root: Path) -> dict[str, dict[str, object]]:
    metadata = json.loads((root / "X23b" / "data" / "metadata.json").read_text())
    return {str(system["id"]): system for system in metadata["systems"]}


def gas_energy(root: Path, method: str, system: str) -> float | None:
    stem = f"{system}_{method}_mol_geoopt"
    output = root / "X23b" / "runs" / "molecule_geoopt" / method / stem / f"{stem}.out"
    return parse_energy(output)


def source_cases(source_csv: Path, variant: str, methods: set[str] | None, systems: set[str] | None) -> list[Case]:
    rows = read_csv(source_csv)
    cases: list[Case] = []
    for row in rows:
        if row.get("variant") != variant or row.get("opt_completed") != "True":
            continue
        method = row["method"]
        system = row["system"]
        if method not in LEGACY_METHODS:
            if methods is None or method in methods:
                raise ValueError(
                    f"legacy final-kpoint runner does not support {method}; "
                    "use X23b/scripts/x23b_final_kpoint_sp.py"
                )
            continue
        if methods is not None and method not in methods:
            continue
        if systems is not None and system not in systems:
            continue
        run_dir = Path(row["run_dir"])
        restarts = sorted(
            path
            for path in run_dir.glob("*-1.restart")
            if ".bak" not in path.name and path.is_file()
        )
        if not restarts:
            raise FileNotFoundError(f"No final restart found in {run_dir}")
        cases.append(Case(method, system, restarts[-1]))
    return sorted(cases, key=lambda case: (case.method, case.system))


def strip_motion(text: str) -> str:
    return re.sub(r"\n\s*&MOTION\b.*?\n\s*&END\s+MOTION\s*\n?", "\n", text, flags=re.I | re.S)


def replace_run_type(text: str) -> str:
    if re.search(r"^\s*RUN_TYPE\s+", text, re.M | re.I):
        return re.sub(r"(^\s*RUN_TYPE\s+).*$", r"\g<1>ENERGY", text, count=1, flags=re.M | re.I)
    return re.sub(r"(^\s*&GLOBAL\b.*?$)", r"\1\n   RUN_TYPE ENERGY", text, count=1, flags=re.M | re.I)


def replace_project(text: str, project: str) -> str:
    if re.search(r"^\s*PROJECT_NAME\s+", text, re.M | re.I):
        return re.sub(r"(^\s*PROJECT_NAME\s+).*$", rf'\g<1>"{project}"', text, count=1, flags=re.M | re.I)
    if re.search(r"^\s*PROJECT\s+", text, re.M | re.I):
        return re.sub(r"(^\s*PROJECT\s+).*$", rf"\g<1>{project}", text, count=1, flags=re.M | re.I)
    return re.sub(r"(^\s*&GLOBAL\b.*?$)", rf'\1\n   PROJECT_NAME "{project}"', text, count=1, flags=re.M | re.I)


def remove_kpoints(text: str) -> str:
    return re.sub(r"\n\s*&KPOINTS\b.*?\n\s*&END\s+KPOINTS\s*\n?", "\n", text, flags=re.I | re.S)


def insert_kpoints(text: str, mesh: str) -> str:
    block = KPOINTS[mesh]
    text = remove_kpoints(text)
    if not block:
        return text

    def repl(match: re.Match[str]) -> str:
        return f"{match.group(0)}\n{block}"

    return re.sub(r"(^\s*&END\s+QS\s*$)", repl, text, count=1, flags=re.M | re.I)


def make_input(case: Case, mesh: str, run_dir: Path) -> Path:
    project = f"{case.system}_{case.method}_{mesh}_cellopt_final_sp".replace("-", "_")
    text = case.source_restart.read_text()
    text = strip_motion(text)
    text = replace_run_type(text)
    text = replace_project(text, project)
    text = insert_kpoints(text, mesh)
    path = run_dir / f"{project}.inp"
    path.write_text(text)
    return path


def row_for_output(
    root: Path,
    meta: dict[str, dict[str, object]],
    case: Case,
    mesh: str,
    run_dir: Path,
    output: Path,
    returncode: int | None,
) -> dict[str, str]:
    system_meta = meta[case.system]
    n_mol = int(system_meta["molecules_per_cell"])
    ref_lattice = float(system_meta["ref_energy"])
    ref_volume = float(system_meta["x23b_same_cell_ref_volume"])
    energy = parse_energy(output)
    gas = gas_energy(root, case.method, case.system)
    lattice = None
    error = None
    if energy is not None and gas is not None:
        lattice = (gas - energy / n_mol) * HARTREE_TO_KJMOL
        error = lattice - ref_lattice
    volume = parse_last_float(output, "CELL| Volume [angstrom^3]")
    volume_error = None
    if volume is not None:
        volume_error = 100.0 * (volume - ref_volume) / ref_volume
    return {
        "method": case.method,
        "system": case.system,
        "mesh": mesh,
        "returncode": "" if returncode is None else str(returncode),
        "program_ended": str(completed_cp2k(output)),
        "energy_hartree": "" if energy is None else f"{energy:.12f}",
        "gas_energy_hartree": "" if gas is None else f"{gas:.12f}",
        "lattice_energy_kJmol": "" if lattice is None else f"{lattice:.6f}",
        "x23b_ref_lattice_energy_kJmol": f"{ref_lattice:.6f}",
        "error_kJmol": "" if error is None else f"{error:.6f}",
        "volume_A3": "" if volume is None else f"{volume:.6f}",
        "x23b_same_cell_ref_volume_A3": f"{ref_volume:.6f}",
        "volume_error_percent": "" if volume_error is None else f"{volume_error:.6f}",
        "source_restart": str(case.source_restart),
        "run_dir": str(run_dir),
        "output": str(output),
    }


def run_case(args: argparse.Namespace, meta: dict[str, dict[str, object]], case: Case, mesh: str) -> dict[str, str]:
    run_dir = args.out / "runs" / case.method / case.system / mesh
    out_file = run_dir / "cp2k.out"
    if args.resume and completed_cp2k(out_file):
        return row_for_output(args.benchmark_root, meta, case, mesh, run_dir, out_file, 0)
    if run_dir.exists():
        shutil.rmtree(run_dir)
    run_dir.mkdir(parents=True)
    inp = make_input(case, mesh, run_dir)
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = str(args.omp_threads)
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
    return row_for_output(args.benchmark_root, meta, case, mesh, run_dir, out_file, proc.returncode)


def stats(values: list[float]) -> tuple[float, float, float, float]:
    n = len(values)
    me = sum(values) / n
    mae = sum(abs(value) for value in values) / n
    rmse = math.sqrt(sum(value * value for value in values) / n)
    maxae = max(abs(value) for value in values)
    return me, mae, rmse, maxae


def summarize(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    summary: list[dict[str, str]] = []
    for method in sorted({row["method"] for row in rows}):
        for mesh in sorted({row["mesh"] for row in rows}, key=lambda m: ["gamma", "k222", "k333"].index(m)):
            subset = [
                row for row in rows
                if row["method"] == method and row["mesh"] == mesh and row["program_ended"] == "True" and row["error_kJmol"]
            ]
            if not subset:
                continue
            me, mae, rmse, maxae = stats([float(row["error_kJmol"]) for row in subset])
            summary.append(
                {
                    "quantity": "lattice_energy_kJmol",
                    "calculation": "cellopt_final_single_point",
                    "mesh": mesh,
                    "method": f"{method}-xTB",
                    "n": str(len(subset)),
                    "ME": f"{me:.6f}",
                    "MAE": f"{mae:.6f}",
                    "RMSE": f"{rmse:.6f}",
                    "MaxAE": f"{maxae:.6f}",
                }
            )
    return summary


def mesh_deltas(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    by_key = {(row["method"], row["system"], row["mesh"]): row for row in rows}
    out: list[dict[str, str]] = []
    for method, system, mesh in sorted(by_key):
        if mesh == "gamma":
            continue
        gamma = by_key.get((method, system, "gamma"))
        row = by_key[(method, system, mesh)]
        if gamma is None or not gamma["error_kJmol"] or not row["error_kJmol"]:
            continue
        gerr = float(gamma["error_kJmol"])
        merr = float(row["error_kJmol"])
        genergy = float(gamma["lattice_energy_kJmol"])
        menergy = float(row["lattice_energy_kJmol"])
        out.append(
            {
                "method": f"{method}-xTB",
                "system": system,
                "mesh": mesh,
                "gamma_lattice_energy_kJmol": f"{genergy:.6f}",
                "mesh_lattice_energy_kJmol": f"{menergy:.6f}",
                "delta_lattice_energy_kJmol": f"{menergy - genergy:.6f}",
                "gamma_error_kJmol": f"{gerr:.6f}",
                "mesh_error_kJmol": f"{merr:.6f}",
                "delta_error_kJmol": f"{merr - gerr:.6f}",
                "delta_abs_error_kJmol": f"{abs(merr) - abs(gerr):.6f}",
            }
        )
    return out


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cp2k", type=Path, required=True)
    parser.add_argument("--benchmark-root", type=Path, default=REPOSITORY)
    parser.add_argument("--source-csv", type=Path, required=True)
    parser.add_argument("--variant", default="cg_2pnt_keep_angles")
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--jobs", type=int, default=4)
    parser.add_argument("--omp-threads", type=int, default=1)
    parser.add_argument("--mesh", action="append", choices=sorted(KPOINTS))
    parser.add_argument("--method", action="append", choices=["GFN1", "GFN2"])
    parser.add_argument("--system", action="append")
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args()

    meshes = args.mesh if args.mesh else ["gamma", "k222", "k333"]
    methods = set(args.method) if args.method else None
    systems = set(args.system) if args.system else None
    args.out.mkdir(parents=True, exist_ok=True)
    meta = metadata_by_system(args.benchmark_root)
    cases = source_cases(args.source_csv, args.variant, methods, systems)

    rows: list[dict[str, str]] = []
    columns = [
        "method",
        "system",
        "mesh",
        "returncode",
        "program_ended",
        "energy_hartree",
        "gas_energy_hartree",
        "lattice_energy_kJmol",
        "x23b_ref_lattice_energy_kJmol",
        "error_kJmol",
        "volume_A3",
        "x23b_same_cell_ref_volume_A3",
        "volume_error_percent",
        "source_restart",
        "run_dir",
        "output",
    ]
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as pool:
        futures = [pool.submit(run_case, args, meta, case, mesh) for case in cases for mesh in meshes]
        for future in concurrent.futures.as_completed(futures):
            row = future.result()
            print(row["method"], row["system"], row["mesh"], "rc", row["returncode"], "err", row["error_kJmol"], flush=True)
            rows.append(row)
            write_csv(args.out / "x23b_cellopt_final_kpoint_sp_rows.partial.csv", rows, columns)

    rows.sort(key=lambda row: (row["method"], row["system"], row["mesh"]))
    write_csv(args.out / "x23b_cellopt_final_kpoint_sp_rows.csv", rows, columns)
    summary = summarize(rows)
    write_csv(args.out / "x23b_cellopt_final_kpoint_sp_summary.csv", summary, list(summary[0].keys()))
    deltas = mesh_deltas(rows)
    write_csv(args.out / "x23b_cellopt_final_kpoint_sp_mesh_deltas.csv", deltas, list(deltas[0].keys()))


if __name__ == "__main__":
    main()
