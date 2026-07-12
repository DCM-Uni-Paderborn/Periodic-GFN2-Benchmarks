#!/usr/bin/env python3
from __future__ import annotations

import argparse
import concurrent.futures
import csv
import json
import os
import re
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path


HARTREE_TO_KJMOL = 2625.499638
FLOAT = r"[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[Ee][-+]?\d+)?"
REPOSITORY = Path(__file__).resolve().parents[1]


@dataclass(frozen=True)
class Case:
    method: str
    system: str


@dataclass(frozen=True)
class Variant:
    name: str
    accuracy: float | None = None
    optimizer: str = "BFGS"
    trust_radius: float | None = None
    max_dr: float | None = None
    rms_dr: float | None = None
    max_force: float | None = None
    rms_force: float | None = None
    pressure_tol_bar: float | None = None
    keep_angles: bool = False
    keep_symmetry: bool = False
    lbfgs: bool = False
    cg_2pnt: bool = False


DEFAULT_CASES = [
    Case("GFN2", "adamantane"),
    Case("GFN2", "cytosine"),
    Case("GFN2", "oxalic_acid_beta"),
    Case("GFN1", "acetic_acid"),
]

VARIANTS = [
    Variant("acc005", accuracy=0.05),
    Variant("acc001", accuracy=0.01),
    Variant("bfgs_trust010", trust_radius=0.10),
    Variant("bfgs_trust005", trust_radius=0.05),
    Variant("tight_trust005", trust_radius=0.05, max_dr=1.0e-3, rms_dr=5.0e-4),
    Variant("acc001_trust005", accuracy=0.01, trust_radius=0.05),
    Variant("cg_2pnt", optimizer="CG", cg_2pnt=True),
    Variant("cg_2pnt_keep_angles", optimizer="CG", cg_2pnt=True, keep_angles=True),
    Variant("keep_angles", keep_angles=True),
    Variant("keep_sym", keep_symmetry=True),
    Variant("lbfgs_keep_sym", optimizer="LBFGS", lbfgs=True, keep_symmetry=True),
]

SYMMETRY = {
    ("GFN2", "adamantane"): "TETRAGONAL_AB",
    ("GFN2", "cytosine"): "ORTHORHOMBIC",
    ("GFN2", "oxalic_acid_beta"): "MONOCLINIC",
    ("GFN1", "acetic_acid"): "ORTHORHOMBIC",
}


def parse_energy(output: Path) -> float | None:
    if not output.exists():
        return None
    energy = None
    for line in output.read_text(errors="ignore").splitlines():
        if "ENERGY| Total FORCE_EVAL" in line:
            energy = float(line.split()[-1])
    return energy


def completed_optimization(output: Path) -> bool:
    if not output.exists():
        return False
    text = output.read_text(errors="ignore")
    if "MAXIMUM NUMBER OF OPTIMIZATION STEPS REACHED" in text:
        return False
    return "GEOMETRY OPTIMIZATION COMPLETED" in text or "CELL OPTIMIZATION COMPLETED" in text


def completed_cp2k_run(output: Path) -> bool:
    if not output.exists():
        return False
    text = output.read_text(errors="ignore")
    if "ABORT" in text or re.search(r"SCF.*NOT|NOT.*SCF|DID NOT CONVERGE|convergence failure", text, re.I):
        return False
    return "PROGRAM ENDED" in text


def parse_last_float(text: str, label: str) -> float | None:
    value = None
    pattern = re.compile(rf"{re.escape(label)}:?\s+({FLOAT})", re.I)
    for match in pattern.finditer(text):
        value = float(match.group(1))
    return value


def parse_output(output: Path) -> dict[str, object]:
    result: dict[str, object] = {
        "returncode": "",
        "program_ended": False,
        "opt_completed": False,
        "max_iter_reached": False,
        "energy_hartree": None,
        "volume_A3": None,
        "last_step": None,
        "last_pressure_bar": None,
        "last_max_step": None,
        "last_rms_step": None,
        "last_max_gradient": None,
        "last_rms_gradient": None,
    }
    if not output.exists():
        return result
    text = output.read_text(errors="ignore")
    result["program_ended"] = "PROGRAM ENDED" in text
    result["opt_completed"] = completed_optimization(output)
    result["max_iter_reached"] = "MAXIMUM NUMBER OF OPTIMIZATION STEPS REACHED" in text
    result["energy_hartree"] = parse_energy(output)
    result["volume_A3"] = parse_last_float(text, "CELL| Volume [angstrom^3]")
    if result["volume_A3"] is None:
        result["volume_A3"] = parse_last_float(text, "CELL| Volume")
    result["last_pressure_bar"] = parse_last_float(text, "OPT| Internal pressure [bar]")
    if result["last_pressure_bar"] is None:
        result["last_pressure_bar"] = parse_last_float(text, "OPT| Pressure deviation [bar]")
    result["last_max_step"] = parse_last_float(text, "OPT| Maximum step size")
    result["last_rms_step"] = parse_last_float(text, "OPT| RMS step size")
    result["last_max_gradient"] = parse_last_float(text, "OPT| Maximum gradient")
    result["last_rms_gradient"] = parse_last_float(text, "OPT| RMS gradient")
    step = None
    step_value = parse_last_float(text, "OPT| Step number")
    if step_value is not None:
        step = int(step_value)
    for match in re.finditer(r"^\s*Step number\s+(\d+)", text, re.M):
        step = int(match.group(1))
    if step is None:
        for match in re.finditer(r"^\s*Informations at step\s*=\s*(\d+)", text, re.M):
            step = int(match.group(1))
    result["last_step"] = step
    return result


def metadata_by_system(root: Path) -> dict[str, dict[str, object]]:
    metadata = json.loads((root / "X23b" / "data" / "metadata.json").read_text())
    return {str(system["id"]): system for system in metadata["systems"]}


def gas_energy(root: Path, method: str, system: str) -> float | None:
    stem = f"{system}_{method}_mol_geoopt"
    output = root / "X23b" / "runs" / "molecule_geoopt" / method / stem / f"{stem}.out"
    if not completed_optimization(output):
        return None
    return parse_energy(output)


def input_path(root: Path, method: str, system: str) -> Path:
    return root / "X23b" / "inputs" / "cellopt_gamma" / method / f"{system}_{method}_gamma_cellopt.inp"


def baseline_output_path(root: Path, method: str, system: str) -> Path:
    stem = f"{system}_{method}_gamma_cellopt"
    run_dir = root / "X23b" / "runs" / "cellopt_gamma" / method / stem
    continuation = run_dir / "continue_800.out"
    if completed_optimization(continuation):
        return continuation
    return run_dir / f"{stem}.out"


def strip_motion(text: str) -> str:
    return re.sub(r"\n&MOTION\b.*?\n&END\s+MOTION\s*\n?", "\n", text, flags=re.I | re.S)


def replace_project(text: str, project: str) -> str:
    return re.sub(r"(^\s*PROJECT\s+).*$", rf"\g<1>{project}", text, count=1, flags=re.M | re.I)


def replace_accuracy(text: str, accuracy: float | None) -> str:
    if accuracy is None:
        return text
    return re.sub(r"(^\s*ACCURACY\s+).*$", rf"\g<1>{accuracy:g}", text, count=1, flags=re.M | re.I)


def insert_cell_symmetry(text: str, symmetry: str) -> str:
    def repl(match: re.Match[str]) -> str:
        block = match.group(0)
        if re.search(r"^\s*SYMMETRY\s+", block, re.M | re.I):
            return re.sub(r"(^\s*SYMMETRY\s+).*$", rf"\g<1>{symmetry}", block, count=1, flags=re.M | re.I)
        return re.sub(r"(^\s*PERIODIC\s+XYZ\s*$)", rf"\1\n      SYMMETRY {symmetry}", block, count=1, flags=re.M | re.I)

    return re.sub(r"(^\s*&CELL\b.*?^\s*&END\s+CELL\s*$)", repl, text, count=1, flags=re.M | re.S | re.I)


def cell_opt_block(variant: Variant, max_iter: int) -> str:
    optimizer = "LBFGS" if variant.lbfgs else variant.optimizer
    lines = [
        "",
        "&MOTION",
        "  &CELL_OPT",
        f"    OPTIMIZER {optimizer}",
        f"    MAX_ITER {max_iter}",
        "    EXTERNAL_PRESSURE [bar] 0.0",
    ]
    if variant.keep_angles:
        lines.append("    KEEP_ANGLES T")
    if variant.keep_symmetry:
        lines.append("    KEEP_SYMMETRY T")
    if variant.max_dr is not None:
        lines.append(f"    MAX_DR {variant.max_dr:.8g}")
    if variant.rms_dr is not None:
        lines.append(f"    RMS_DR {variant.rms_dr:.8g}")
    if variant.max_force is not None:
        lines.append(f"    MAX_FORCE {variant.max_force:.8g}")
    if variant.rms_force is not None:
        lines.append(f"    RMS_FORCE {variant.rms_force:.8g}")
    if variant.pressure_tol_bar is not None:
        lines.append(f"    PRESSURE_TOLERANCE [bar] {variant.pressure_tol_bar:.8g}")
    if optimizer == "BFGS" and variant.trust_radius is not None:
        lines += [
            "    &BFGS",
            f"      TRUST_RADIUS [angstrom] {variant.trust_radius:.8g}",
            "    &END BFGS",
        ]
    if optimizer == "CG" and variant.cg_2pnt:
        lines += [
            "    &CG",
            "      &LINE_SEARCH",
            "        TYPE 2PNT",
            "      &END LINE_SEARCH",
            "    &END CG",
        ]
    lines += [
        "  &END CELL_OPT",
        "&END MOTION",
        "",
    ]
    return "\n".join(lines)


def make_input(root: Path, case: Case, variant: Variant, run_dir: Path, max_iter: int) -> Path:
    text = input_path(root, case.method, case.system).read_text()
    project = f"{case.system}_{case.method}_{variant.name}".replace("-", "_")
    text = replace_project(text, project)
    text = replace_accuracy(text, variant.accuracy)
    text = strip_motion(text).rstrip() + "\n"
    if variant.keep_symmetry:
        text = insert_cell_symmetry(text, SYMMETRY[(case.method, case.system)])
    text += cell_opt_block(variant, max_iter)
    path = run_dir / f"{project}.inp"
    path.write_text(text)
    return path


def format_float(value: object, digits: int = 6) -> str:
    if value is None or value == "":
        return ""
    return f"{float(value):.{digits}f}"


def row_for_output(
    root: Path,
    meta: dict[str, dict[str, object]],
    case: Case,
    variant: str,
    output: Path,
    run_dir: Path,
    source: str,
    returncode: int | None = None,
) -> dict[str, str]:
    info = parse_output(output)
    if returncode is not None:
        info["returncode"] = returncode
    system_meta = meta[case.system]
    n_mol = int(system_meta["molecules_per_cell"])
    ref_lattice = float(system_meta["ref_energy"])
    ref_volume = float(system_meta["x23b_same_cell_ref_volume"])
    gas = gas_energy(root, case.method, case.system)
    energy = info["energy_hartree"]
    volume = info["volume_A3"]
    lattice = None
    lattice_error = None
    volume_error = None
    if energy is not None and gas is not None:
        lattice = (float(gas) - float(energy) / n_mol) * HARTREE_TO_KJMOL
        lattice_error = lattice - ref_lattice
    if volume is not None:
        volume_error = 100.0 * (float(volume) - ref_volume) / ref_volume
    row = {
        "method": case.method,
        "system": case.system,
        "variant": variant,
        "source": source,
        "returncode": str(info["returncode"]),
        "program_ended": str(bool(info["program_ended"])),
        "opt_completed": str(bool(info["opt_completed"])),
        "max_iter_reached": str(bool(info["max_iter_reached"])),
        "last_step": "" if info["last_step"] is None else str(info["last_step"]),
        "energy_hartree": format_float(energy, 12),
        "gas_energy_hartree": format_float(gas, 12),
        "lattice_energy_kJmol": format_float(lattice, 6),
        "x23b_ref_lattice_energy_kJmol": format_float(ref_lattice, 6),
        "error_kJmol": format_float(lattice_error, 6),
        "volume_A3": format_float(volume, 6),
        "x23b_same_cell_ref_volume_A3": format_float(ref_volume, 6),
        "volume_error_percent": format_float(volume_error, 6),
        "last_pressure_bar": format_float(info["last_pressure_bar"], 6),
        "last_max_step": format_float(info["last_max_step"], 10),
        "last_rms_step": format_float(info["last_rms_step"], 10),
        "last_max_gradient": format_float(info["last_max_gradient"], 10),
        "last_rms_gradient": format_float(info["last_rms_gradient"], 10),
        "run_dir": str(run_dir),
        "output": str(output),
    }
    return row


def run_case(args: argparse.Namespace, meta: dict[str, dict[str, object]], case: Case, variant: Variant) -> dict[str, str]:
    run_dir = args.out / "runs" / case.method / case.system / variant.name
    out_file = run_dir / "cp2k.out"
    if args.resume and out_file.exists() and completed_cp2k_run(out_file):
        return row_for_output(args.benchmark_root, meta, case, variant.name, out_file, run_dir, "variant_resumed")
    if run_dir.exists():
        shutil.rmtree(run_dir)
    run_dir.mkdir(parents=True)
    inp = make_input(args.benchmark_root, case, variant, run_dir, args.max_iter)
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
    return row_for_output(args.benchmark_root, meta, case, variant.name, out_file, run_dir, "variant", proc.returncode)


def write_csv(path: Path, rows: list[dict[str, str]], columns: list[str]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)


def summarize(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    out: list[dict[str, str]] = []
    for variant in sorted({row["variant"] for row in rows}):
        subset = [row for row in rows if row["variant"] == variant and row["opt_completed"] == "True"]
        summary = {"variant": variant, "n_complete": str(len(subset))}
        for key in ["error_kJmol", "volume_error_percent"]:
            vals = [float(row[key]) for row in subset if row[key]]
            summary[f"{key}_ME"] = f"{sum(vals) / len(vals):.6f}" if vals else ""
            summary[f"{key}_MAE"] = f"{sum(abs(v) for v in vals) / len(vals):.6f}" if vals else ""
            summary[f"{key}_MaxAE"] = f"{max(abs(v) for v in vals):.6f}" if vals else ""
        out.append(summary)
    return out


def parse_cases(args: argparse.Namespace, meta: dict[str, dict[str, object]]) -> list[Case]:
    if args.all_systems:
        methods = [args.method] if args.method else ["GFN1", "GFN2"]
        cases = [Case(method, system) for method in methods for system in sorted(meta)]
    else:
        cases = DEFAULT_CASES
    if args.case:
        parsed = []
        for item in args.case:
            if ":" not in item:
                raise SystemExit("--case must be METHOD:SYSTEM, for example GFN2:adamantane")
            method, system = item.split(":", 1)
            parsed.append(Case(method, system))
        cases = parsed
    if args.method:
        cases = [case for case in cases if case.method == args.method]
    if args.system:
        allowed = set(args.system)
        cases = [case for case in cases if case.system in allowed]
    return cases


def parse_variants(args: argparse.Namespace) -> list[Variant]:
    variants = VARIANTS
    if args.variant:
        wanted = set(args.variant)
        variants = [variant for variant in VARIANTS if variant.name in wanted]
        missing = wanted.difference({variant.name for variant in variants})
        if missing:
            raise SystemExit(f"Unknown variant(s): {', '.join(sorted(missing))}")
    return variants


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cp2k", type=Path, required=True)
    parser.add_argument("--benchmark-root", type=Path, default=REPOSITORY)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--jobs", type=int, default=4)
    parser.add_argument("--max-iter", type=int, default=800)
    parser.add_argument("--omp-threads", type=int, default=1)
    parser.add_argument("--method", choices=["GFN1", "GFN2"])
    parser.add_argument("--system", action="append")
    parser.add_argument("--case", action="append")
    parser.add_argument("--all-systems", action="store_true")
    parser.add_argument("--variant", action="append")
    parser.add_argument("--skip-baseline", action="store_true")
    parser.add_argument("--skip-run", action="store_true")
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args()

    args.out.mkdir(parents=True, exist_ok=True)
    meta = metadata_by_system(args.benchmark_root)
    cases = parse_cases(args, meta)
    variants = parse_variants(args)
    rows: list[dict[str, str]] = []

    if not args.skip_baseline:
        for case in cases:
            output = baseline_output_path(args.benchmark_root, case.method, case.system)
            rows.append(row_for_output(args.benchmark_root, meta, case, "baseline_existing", output, output.parent, "existing"))

    if not args.skip_run:
        with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as pool:
            futures = [pool.submit(run_case, args, meta, case, variant) for case in cases for variant in variants]
            for future in concurrent.futures.as_completed(futures):
                rows.append(future.result())
                write_csv(args.out / "cellopt_variant_rows.partial.csv", rows, COLUMNS)

    rows.sort(key=lambda row: (row["method"], row["system"], row["variant"]))
    write_csv(args.out / "cellopt_variant_rows.csv", rows, COLUMNS)
    write_csv(args.out / "cellopt_variant_summary.csv", summarize(rows), SUMMARY_COLUMNS)


COLUMNS = [
    "method",
    "system",
    "variant",
    "source",
    "returncode",
    "program_ended",
    "opt_completed",
    "max_iter_reached",
    "last_step",
    "energy_hartree",
    "gas_energy_hartree",
    "lattice_energy_kJmol",
    "x23b_ref_lattice_energy_kJmol",
    "error_kJmol",
    "volume_A3",
    "x23b_same_cell_ref_volume_A3",
    "volume_error_percent",
    "last_pressure_bar",
    "last_max_step",
    "last_rms_step",
    "last_max_gradient",
    "last_rms_gradient",
    "run_dir",
    "output",
]

SUMMARY_COLUMNS = [
    "variant",
    "n_complete",
    "error_kJmol_ME",
    "error_kJmol_MAE",
    "error_kJmol_MaxAE",
    "volume_error_percent_ME",
    "volume_error_percent_MAE",
    "volume_error_percent_MaxAE",
]


if __name__ == "__main__":
    main()
