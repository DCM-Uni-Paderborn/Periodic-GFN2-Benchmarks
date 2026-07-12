#!/usr/bin/env python3
from __future__ import annotations

import argparse
import concurrent.futures
import csv
import os
import re
import shutil
import subprocess
from pathlib import Path


FLOAT = r"[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[Ee][-+]?\d+)?"
REPOSITORY = Path(__file__).resolve().parents[1]


def strip_motion(text: str) -> str:
    return re.sub(r"\n&MOTION\b.*?\n&END\s+MOTION\s*\n?", "\n", text, flags=re.I | re.S)


def ensure_analytical_stress(text: str) -> str:
    if re.search(r"^\s*STRESS_TENSOR\b", text, flags=re.I | re.M):
        return text

    method = re.compile(r"^(?P<indent>\s*)METHOD\s+(?:QS|QUICKSTEP)\s*$", flags=re.I | re.M)

    def add_stress(match: re.Match[str]) -> str:
        return f"{match.group(0)}\n{match.group('indent')}STRESS_TENSOR ANALYTICAL"

    updated, count = method.subn(add_stress, text, count=1)
    if count != 1:
        raise ValueError("Could not locate FORCE_EVAL METHOD line for stress-tensor injection")
    return updated


def prepare_reference_cli_program(run_dir: Path, executable: Path) -> str:
    target = executable.resolve(strict=True)
    link = run_dir / "tblite-reference-cli"
    if link.exists() or link.is_symlink():
        if link.resolve() != target:
            link.unlink()
    if not link.exists():
        link.symlink_to(target)
    return f"./{link.name}"


def inject_reference_cli(text: str, tblite: Path | str, prefix: str, keep_files: bool) -> str:
    text = re.sub(r"\bRUN_TYPE\s+\S+", "RUN_TYPE ENERGY_FORCE", text, count=1)
    text = strip_motion(text)
    text = ensure_analytical_stress(text)
    if "REFERENCE_CLI" in text:
        return text
    block = [
        "          &REFERENCE_CLI",
        f"            PROGRAM_NAME {tblite}",
        "            CHECK_ENERGY T",
        "            CHECK_FORCES T",
        "            CHECK_VIRIAL T",
        "            ERROR_LIMIT 1.0E-6",
        "            STOP_ON_ERROR F",
        f"            PREFIX {prefix}",
        "            WORK_DIRECTORY .",
        f"            KEEP_FILES {'T' if keep_files else 'F'}",
        "          &END REFERENCE_CLI",
    ]
    pattern = re.compile(r"(^\s*&END\s+TBLITE\s*$)", re.M | re.I)
    return pattern.sub("\n".join(block) + "\n\\1", text, count=1)


def parse_reference_cli(out_file: Path) -> dict[str, str]:
    text = out_file.read_text(errors="replace")
    result = {
        "energy_cp2k_hartree": "",
        "energy_cli_hartree": "",
        "energy_absdiff_hartree": "",
        "gradient_diff_sum": "",
        "gradient_diff_max": "",
        "virial_diff_sum": "",
        "virial_diff_max": "",
        "exceeded_error_limit": "False",
        "skipped": "False",
    }
    if "tblite reference CLI check skipped" in text:
        result["skipped"] = "True"
    m = re.search(r"Energy CP2K/CLI/absdiff:\s*(%s)\s+(%s)\s+(%s)" % (FLOAT, FLOAT, FLOAT), text)
    if m:
        result["energy_cp2k_hartree"], result["energy_cli_hartree"], result["energy_absdiff_hartree"] = m.groups()
    m = re.search(r"Gradient diff sum/max:\s*(%s)\s+(%s)" % (FLOAT, FLOAT), text)
    if m:
        result["gradient_diff_sum"], result["gradient_diff_max"] = m.groups()
    m = re.search(r"Virial diff sum/max:\s*(%s)\s+(%s)" % (FLOAT, FLOAT), text)
    if m:
        result["virial_diff_sum"], result["virial_diff_max"] = m.groups()
    result["exceeded_error_limit"] = str("tblite reference CLI deviation exceeded ERROR_LIMIT" in text)
    return result


def x23b_root(root: Path) -> Path:
    nested = root / "X23b"
    if nested.is_dir():
        return nested
    if (root / "runs").is_dir() and (root / "inputs").is_dir():
        return root
    raise ValueError(f"X23b benchmark tree not found below {root}")


def portable_source(path: Path, benchmark_root: Path) -> str:
    root = x23b_root(benchmark_root.resolve()).resolve()
    try:
        relative = path.resolve().relative_to(root)
    except ValueError:
        return path.name
    return (Path("X23b") / relative).as_posix()


def discover_initial(root: Path, method: str) -> list[tuple[str, str, str, Path]]:
    base = x23b_root(root) / "inputs" / "cellopt_gamma" / method
    suffix = f"_{method}_gamma_cellopt.inp"
    return [("initial", method, p.name.removesuffix(suffix), p) for p in sorted(base.glob(f"*{suffix}"))]


def discover_restarts(root: Path, method: str) -> list[tuple[str, str, str, Path]]:
    base = x23b_root(root) / "runs" / "cellopt_gamma" / method
    cases: list[tuple[str, str, str, Path]] = []
    suffix = f"_{method}_gamma_cellopt"
    project_method = method.replace("-", "_")
    for run_dir in sorted(base.glob(f"*{suffix}")):
        system = run_dir.name.removesuffix(suffix)
        restart = run_dir / f"{system}_{project_method}_gamma_cell_opt-1.restart"
        if not restart.exists():
            matches = sorted(run_dir.glob(f"*_{project_method}_gamma_cell_opt-1.restart"))
            restart = matches[0] if matches else restart
        if restart.exists():
            cases.append(("wsc_last_restart", method, system, restart))
    return cases


def latest_numbered_restart(run_dir: Path) -> Path | None:
    restarts: list[tuple[int, Path]] = []
    for path in run_dir.glob("*.restart"):
        match = re.search(r"_(\d+)\.restart$", path.name)
        if match:
            restarts.append((int(match.group(1)), path))
    if not restarts:
        return None
    return max(restarts)[1]


def discover_variant_csv(source_csv: Path, variant: str, source_kind: str) -> list[tuple[str, str, str, Path]]:
    cases: list[tuple[str, str, str, Path]] = []
    with source_csv.open(newline="") as handle:
        for row in csv.DictReader(handle):
            if row.get("variant") != variant or row.get("opt_completed") != "True":
                continue
            restart = latest_numbered_restart(Path(row["run_dir"]))
            if restart is None:
                raise RuntimeError(f"No numbered restart found in {row['run_dir']}")
            cases.append((source_kind, row["method"], row["system"], restart))
    return cases


def run_case(case: tuple[str, str, str, Path], args: argparse.Namespace) -> dict[str, str]:
    source_kind, method, system, source = case
    run_dir = args.out / source_kind / method / system
    run_dir.mkdir(parents=True, exist_ok=True)
    prefix = f"refcli_{source_kind}_{method}_{system}".replace("-", "_")
    inp = run_dir / f"{system}_{source_kind}_reference_cli.inp"
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
        "source_kind": source_kind,
        "method": method,
        "system": system,
        "source": portable_source(source, args.benchmark_root),
        "run_dir": run_dir.relative_to(args.out).as_posix(),
        "returncode": str(proc.returncode),
    }
    row.update(parse_reference_cli(out_file))
    return row


def summarize(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    out: list[dict[str, str]] = []
    for source_kind in sorted({row["source_kind"] for row in rows}):
        methods = sorted({row.get("method", "") for row in rows if row["source_kind"] == source_kind})
        for method in methods:
            subset = [
                row
                for row in rows
                if row["source_kind"] == source_kind
                and row.get("method", "") == method
                and row["returncode"] == "0"
                and row["skipped"] == "False"
            ]
            summary = {"source_kind": source_kind, "method": method, "n": str(len(subset))}
            for key in ["energy_absdiff_hartree", "gradient_diff_sum", "gradient_diff_max", "virial_diff_sum", "virial_diff_max"]:
                values = [float(row[key]) for row in subset if row.get(key)]
                summary[f"{key}_max"] = f"{max(values):.12e}" if values else ""
                summary[f"{key}_sum"] = f"{sum(values):.12e}" if values else ""
            out.append(summary)
    return out


def write_csv(path: Path, rows: list[dict[str, str]], columns: list[str]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cp2k", type=Path, required=True)
    parser.add_argument("--tblite", type=Path, required=True)
    parser.add_argument("--benchmark-root", type=Path, default=REPOSITORY)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--jobs", type=int, default=6)
    parser.add_argument("--method", choices=["GFN1", "GFN2"], default="GFN2")
    parser.add_argument("--only-initial", action="store_true")
    parser.add_argument("--only-restarts", action="store_true")
    parser.add_argument("--keep-files", action="store_true")
    parser.add_argument("--system", action="append")
    parser.add_argument("--source-csv", type=Path)
    parser.add_argument("--variant", default="cg_2pnt")
    parser.add_argument("--source-kind", default="variant_final")
    args = parser.parse_args()

    if args.out.exists():
        shutil.rmtree(args.out)
    args.out.mkdir(parents=True)

    cases: list[tuple[str, str, str, Path]] = []
    if args.source_csv:
        cases.extend(discover_variant_csv(args.source_csv, args.variant, args.source_kind))
    elif not args.only_restarts:
        cases.extend(discover_initial(args.benchmark_root, args.method))
    if not args.source_csv and not args.only_initial:
        cases.extend(discover_restarts(args.benchmark_root, args.method))
    if args.source_csv:
        cases = [case for case in cases if case[1] == args.method]
    if args.system:
        wanted = set(args.system)
        cases = [case for case in cases if case[2] in wanted]

    rows: list[dict[str, str]] = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as pool:
        future_map = {pool.submit(run_case, case, args): case for case in cases}
        for future in concurrent.futures.as_completed(future_map):
            row = future.result()
            rows.append(row)
            print(
                row["source_kind"],
                row["method"],
                row["system"],
                "rc",
                row["returncode"],
                "skipped",
                row["skipped"],
                "gmax",
                row["gradient_diff_max"],
                "vmax",
                row["virial_diff_max"],
                flush=True,
            )

    rows.sort(key=lambda row: (row["source_kind"], row["method"], row["system"]))
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
    write_csv(args.out / "reference_cli_summary.csv", summary, list(summary[0].keys()) if summary else ["source_kind", "method", "n"])
    print(args.out)
    if rows and all(row["skipped"] == "True" for row in rows):
        raise SystemExit("all selected REFERENCE_CLI checks were skipped")


if __name__ == "__main__":
    main()
