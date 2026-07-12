#!/usr/bin/env python3
"""Prepare, run, and collect native-Bloch X23b k-point cell optimizations."""

from __future__ import annotations

import argparse
import csv
import fcntl
import json
import os
import re
import shutil
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
METHODS = ("GFN1", "GFN2")
HARTREE_TO_KJMOL = 2625.499638
VARIANT = "k222_cellopt_keep_angles"
FIELDS = (
    "method",
    "system",
    "variant",
    "mesh",
    "source_variant",
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
    "source_restart",
    "run_dir",
    "output",
)

CELL_OPT_LIMITS = {
    "pressure": 100.0,
    "max_step": 0.003,
    "rms_step": 0.0015,
    "max_gradient": 0.00045,
    "rms_gradient": 0.0003,
}

CELL_OPT_PATTERNS = {
    "step": re.compile(r"^\s*OPT\| Step number\s+(\d+)\s*$"),
    "pressure": re.compile(
        r"^\s*OPT\| (?:Internal pressure|Pressure deviation) \[?bar\]?\s+([-+0-9.Ee]+)\s*$"
    ),
    "max_step": re.compile(r"^\s*OPT\| Maximum step size\s+([-+0-9.Ee]+)\s*$"),
    "rms_step": re.compile(r"^\s*OPT\| RMS step size\s+([-+0-9.Ee]+)\s*$"),
    "max_gradient": re.compile(r"^\s*OPT\| Maximum gradient\s+([-+0-9.Ee]+)\s*$"),
    "rms_gradient": re.compile(r"^\s*OPT\| RMS gradient\s+([-+0-9.Ee]+)\s*$"),
}


def systems() -> list[dict[str, object]]:
    metadata = json.loads((ROOT / "data" / "metadata.json").read_text())
    return list(metadata["systems"])


def section_bounds(lines: list[str], section: str) -> tuple[int, int]:
    target = f"&{section.upper()}"
    for start, line in enumerate(lines):
        if line.strip().upper().split(maxsplit=1)[0] != target:
            continue
        depth = 0
        for end in range(start, len(lines)):
            stripped = lines[end].strip().upper()
            if stripped.startswith("&END"):
                depth -= 1
            elif stripped.startswith("&"):
                depth += 1
            if depth == 0:
                return start, end
        break
    raise ValueError(f"section {section} not found or not closed")


def restart_to_k222_input(source: Path, project: str) -> str:
    lines = source.read_text().splitlines()
    try:
        global_start = next(index for index, line in enumerate(lines) if line.strip().upper() == "&GLOBAL")
    except StopIteration as exc:
        raise ValueError(f"GLOBAL section missing in {source}") from exc
    lines = lines[global_start:]

    motion_start, motion_end = section_bounds(lines, "MOTION")
    del lines[motion_start : motion_end + 1]

    project_pattern = re.compile(r"^(\s*PROJECT(?:_NAME)?\s+).*$", re.I)
    for index, line in enumerate(lines):
        match = project_pattern.match(line)
        if match:
            lines[index] = f'{match.group(1)}"{project}"'
            break
    else:
        raise ValueError(f"PROJECT keyword missing in {source}")

    dft_start, dft_end = section_bounds(lines, "DFT")
    if any("&KPOINTS" in line.upper() for line in lines[dft_start:dft_end]):
        raise ValueError(f"source restart already contains KPOINTS: {source}")
    indent = re.match(r"\s*", lines[dft_start]).group(0) + "  "
    kpoints = [
        f"{indent}&KPOINTS",
        f"{indent}  SCHEME MACDONALD 2 2 2 0.25 0.25 0.25",
        f"{indent}  EPS_SYMMETRY 1.0E-8",
        f"{indent}  SYMMETRY T",
        f"{indent}  FULL_GRID F",
        f"{indent}  SYMMETRY_BACKEND SPGLIB",
        f"{indent}  SYMMETRY_REDUCTION_METHOD SPGLIB",
        f"{indent}&END KPOINTS",
    ]
    lines[dft_end:dft_end] = kpoints

    cell_start, cell_end = section_bounds(lines, "CELL")
    if not any("CANONICALIZE" in line.upper() for line in lines[cell_start:cell_end]):
        indent = re.match(r"\s*", lines[cell_start]).group(0) + "  "
        lines.insert(cell_start + 1, f"{indent}CANONICALIZE TRUE")

    lines += [
        "",
        "&MOTION",
        "  &CELL_OPT",
        "    OPTIMIZER CG",
        "    MAX_ITER 500",
        "    EXTERNAL_PRESSURE [bar] 0.0",
        "    KEEP_ANGLES T",
        "    &CG",
        "      &LINE_SEARCH",
        "        TYPE 2PNT",
        "      &END LINE_SEARCH",
        "    &END CG",
        "  &END CELL_OPT",
        "&END MOTION",
    ]
    return "\n".join(lines) + "\n"


def parse_overrides(values: list[str]) -> dict[tuple[str, str], Path]:
    overrides: dict[tuple[str, str], Path] = {}
    for value in values:
        key, separator, path = value.partition("=")
        if not separator or "/" not in key:
            raise ValueError(f"invalid override {value!r}; expected METHOD/system=/path/restart")
        method, system = key.split("/", 1)
        overrides[(method, system)] = Path(path).resolve()
    return overrides


def find_gamma_restart(source_root: Path, method: str, system: str) -> Path:
    directories = (
        source_root / "runs" / method / system / "gamma_cellopt_keep_angles",
        source_root / "runs" / "cellopt_gamma" / method / f"{system}_{method}_gamma_cellopt",
    )
    for directory in directories:
        matches = sorted(directory.glob("*-1.restart"))
        if len(matches) == 1:
            return matches[0]
        if len(matches) > 1:
            raise ValueError(f"expected one final restart in {directory}, found {len(matches)}")
    raise ValueError("final Gamma restart not found in: " + ", ".join(str(path) for path in directories))


def prepare(args: argparse.Namespace) -> None:
    output_root = args.output_root.resolve()
    output_root.parent.mkdir(parents=True, exist_ok=True)
    lock_path = output_root.parent / ".x23b_k222_cellopt_manifest.lock"
    with lock_path.open("w") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        _prepare(args, output_root)


def _prepare(args: argparse.Namespace, output_root: Path) -> None:
    if args.clean and output_root.exists():
        shutil.rmtree(output_root)
    overrides = parse_overrides(args.override)
    selected_methods = (args.method,) if args.method else METHODS
    selected_systems = set(args.system) if args.system else None
    manifest = output_root.parent / "x23b_k222_cellopt_manifest.csv"
    manifest_rows: dict[tuple[str, str], dict[str, str]] = {}
    if manifest.is_file() and not args.clean:
        with manifest.open(newline="") as handle:
            for row in csv.DictReader(handle):
                manifest_rows[(row["method"], row["system"])] = row
    prepared = 0
    for system_data in systems():
        system = str(system_data["id"])
        if selected_systems is not None and system not in selected_systems:
            continue
        for method in selected_methods:
            source = overrides.get((method, system))
            if source is None:
                source = find_gamma_restart(args.gamma_root.resolve(), method, system)
            if not source.is_file():
                raise FileNotFoundError(source)
            project = f"{system}_{method}_{VARIANT}".replace("-", "_")
            run_dir = output_root / method / system / VARIANT
            run_dir.mkdir(parents=True, exist_ok=True)
            input_path = run_dir / f"{project}.inp"
            input_path.write_text(restart_to_k222_input(source, project))
            manifest_rows[(method, system)] = {
                "method": method,
                "system": system,
                "source_restart": str(source),
                "input": str(input_path),
                "run_dir": str(run_dir),
            }
            prepared += 1
    rows = [manifest_rows[key] for key in sorted(manifest_rows)]
    with manifest.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=("method", "system", "source_restart", "input", "run_dir"),
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)
    print(
        f"Prepared {prepared} native-Bloch k222 cell optimizations in {output_root}; "
        f"manifest contains {len(rows)} entries"
    )


def cp2k_completed(output: Path) -> bool:
    if not output.is_file():
        return False
    text = output.read_text(errors="ignore")
    return "PROGRAM ENDED" in text and "GEOMETRY OPTIMIZATION COMPLETED" in text


def cp2k_terminal(output: Path) -> bool:
    if not output.is_file():
        return False
    text = output.read_text(errors="ignore")
    return "PROGRAM ENDED" in text and (
        "GEOMETRY OPTIMIZATION COMPLETED" in text
        or "MAXIMUM NUMBER OF OPTIMIZATION STEPS REACHED" in text
    )


def run_one(input_path: Path, cp2k: Path, force: bool) -> tuple[Path, int, str]:
    run_dir = input_path.parent
    output = run_dir / "cp2k.out"
    with input_path.open() as input_lock:
        try:
            fcntl.flock(input_lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            return input_path, 0, "BUSY"
        if not force and cp2k_terminal(output):
            action = "SKIP_CONVERGED" if cp2k_completed(output) else "SKIP_MAX_ITER"
            return input_path, 0, action
        if force:
            for path in run_dir.iterdir():
                if path != input_path:
                    if path.is_dir():
                        shutil.rmtree(path)
                    else:
                        path.unlink()
        process = subprocess.run(
            [str(cp2k), "-i", input_path.name, "-o", output.name],
            cwd=run_dir,
            check=False,
        )
        (run_dir / "returncode.txt").write_text(f"{process.returncode}\n")
        if process.returncode != 0:
            action = "FAILED"
        elif cp2k_completed(output):
            action = "CONVERGED"
        elif cp2k_terminal(output):
            action = "MAX_ITER"
        else:
            action = "INCOMPLETE"
        return input_path, process.returncode, action


def run(args: argparse.Namespace) -> None:
    inputs = sorted(args.output_root.resolve().glob(f"*/**/{VARIANT}/*.inp"))
    if args.method:
        inputs = [path for path in inputs if args.method in path.parts]
    all_systems = sorted(str(row["id"]) for row in systems())
    if args.system and (args.start_system or args.end_system):
        raise ValueError("--system cannot be combined with --start-system or --end-system")
    if args.system:
        selected_systems = sorted(set(args.system))
    else:
        start = all_systems.index(args.start_system) if args.start_system else 0
        end = all_systems.index(args.end_system) + 1 if args.end_system else len(all_systems)
        if start >= end:
            raise ValueError("--start-system must not follow --end-system")
        selected_systems = all_systems[start:end]
    if args.system or args.start_system or args.end_system:
        inputs = [path for path in inputs if any(system in path.parts for system in selected_systems)]
    expected = len(selected_systems) * (1 if args.method else len(METHODS))
    if len(inputs) != expected:
        selection = f" for {args.method}" if args.method else ""
        raise ValueError(f"expected {expected} prepared inputs{selection}, found {len(inputs)}")
    environment_threads = str(args.threads_per_job)
    os.environ["OMP_NUM_THREADS"] = environment_threads
    os.environ["OPENBLAS_NUM_THREADS"] = environment_threads
    os.environ["MKL_NUM_THREADS"] = environment_threads
    failed = []
    max_iter = []
    with ThreadPoolExecutor(max_workers=args.jobs) as executor:
        futures = {executor.submit(run_one, path, args.cp2k.resolve(), args.force): path for path in inputs}
        for future in as_completed(futures):
            input_path, returncode, action = future.result()
            relative = input_path.relative_to(args.output_root.resolve())
            print(f"{action:14s} {relative} rc={returncode}", flush=True)
            if action in {"MAX_ITER", "SKIP_MAX_ITER"}:
                max_iter.append(relative)
            elif returncode != 0 or action == "INCOMPLETE":
                failed.append(relative)
    if max_iter:
        print(f"{len(max_iter)} optimization(s) require continue-maxiter", flush=True)
    if failed:
        raise SystemExit(f"{len(failed)} CP2K jobs failed")


def latest_numbered_restart(run_dir: Path) -> tuple[int, Path] | None:
    restarts: list[tuple[int, Path]] = []
    for path in run_dir.glob("*-1_*.restart"):
        match = re.search(r"-1_(\d+)\.restart$", path.name)
        if match:
            restarts.append((int(match.group(1)), path))
    return max(restarts) if restarts else None


def optimization_records(output: Path) -> list[dict[str, float | int]]:
    records: list[dict[str, float | int]] = []
    current: dict[str, float | int] = {}
    for line in output.read_text(errors="ignore").splitlines():
        step_match = CELL_OPT_PATTERNS["step"].match(line)
        if step_match:
            if "step" in current:
                records.append(current)
            current = {"step": int(step_match.group(1))}
            continue
        if "step" not in current:
            continue
        for name, pattern in CELL_OPT_PATTERNS.items():
            if name == "step":
                continue
            match = pattern.match(line)
            if match:
                current[name] = float(match.group(1))
                break
    if "step" in current:
        records.append(current)
    return records


def best_polish_restart(run_dir: Path) -> tuple[Path, dict[str, float | int]]:
    candidates: list[tuple[float, float, int, Path, dict[str, float | int]]] = []
    for output in sorted(run_dir.glob("cp2k*.out")):
        for record in optimization_records(output):
            if any(name not in record for name in CELL_OPT_LIMITS):
                continue
            step = int(record["step"])
            restarts = list(run_dir.glob(f"*-1_{step}.restart"))
            if not restarts:
                continue
            restart = max(restarts, key=lambda path: path.stat().st_mtime_ns)
            ratios = [
                abs(float(record[name])) / limit if name == "pressure" else float(record[name]) / limit
                for name, limit in CELL_OPT_LIMITS.items()
            ]
            score = max(ratios)
            norm = sum(value * value for value in ratios)
            candidates.append((score, norm, -step, restart, record))
    if not candidates:
        raise ValueError(f"no complete optimization record with a matching restart in {run_dir}")
    _, _, _, restart, record = min(candidates, key=lambda item: item[:3])
    return restart, record


def set_cell_opt_keyword(lines: list[str], keyword: str, value: str) -> None:
    start, end = section_bounds(lines, "CELL_OPT")
    pattern = re.compile(rf"^(\s*{re.escape(keyword)}\s+).*$", re.I)
    for index in range(start + 1, end):
        match = pattern.match(lines[index])
        if match:
            lines[index] = f"{match.group(1)}{value}"
            return
    indent = re.match(r"\s*", lines[start]).group(0) + "  "
    lines.insert(end, f"{indent}{keyword} {value}")


def bfgs_polish_input(source: Path, project: str, max_iter: int, trust_radius: float) -> str:
    lines = source.read_text().splitlines()
    global_start = next(index for index, line in enumerate(lines) if line.strip().upper() == "&GLOBAL")
    lines = lines[global_start:]

    project_pattern = re.compile(r"^(\s*PROJECT(?:_NAME)?\s+).*$", re.I)
    for index, line in enumerate(lines):
        match = project_pattern.match(line)
        if match:
            lines[index] = f'{match.group(1)}"{project}"'
            break
    else:
        raise ValueError(f"PROJECT keyword missing in {source}")

    cell_start, cell_end = section_bounds(lines, "CELL_OPT")
    for subsection in ("CG", "BFGS", "LBFGS"):
        local_lines = lines[cell_start : cell_end + 1]
        try:
            local_start, local_end = section_bounds(local_lines, subsection)
        except ValueError:
            continue
        del lines[cell_start + local_start : cell_start + local_end + 1]
        cell_start, cell_end = section_bounds(lines, "CELL_OPT")

    settings = {
        "OPTIMIZER": "BFGS",
        "MAX_ITER": str(max_iter),
        "STEP_START_VAL": "0",
        "MAX_DR": "0.003",
        "RMS_DR": "0.0015",
        "MAX_FORCE": "0.00045",
        "RMS_FORCE": "0.0003",
        "PRESSURE_TOLERANCE": "[bar] 100.0",
    }
    for keyword, value in settings.items():
        set_cell_opt_keyword(lines, keyword, value)

    cell_start, cell_end = section_bounds(lines, "CELL_OPT")
    indent = re.match(r"\s*", lines[cell_start]).group(0) + "  "
    lines[cell_end:cell_end] = [
        f"{indent}&BFGS",
        f"{indent}  TRUST_RADIUS [angstrom] {trust_radius:.12g}",
        f"{indent}&END BFGS",
    ]
    return "\n".join(lines) + "\n"


def continuation_input(source: Path, project: str, additional_steps: int) -> tuple[str, int]:
    text = source.read_text()
    step_match = re.search(r"^\s*STEP_START_VAL\s+(\d+)\s*$", text, flags=re.M)
    if step_match is None:
        raise ValueError(f"STEP_START_VAL missing in {source}")
    start_step = int(step_match.group(1))
    text, project_count = re.subn(
        r'^(\s*PROJECT(?:_NAME)?\s+).+$',
        rf'\1"{project}"',
        text,
        count=1,
        flags=re.I | re.M,
    )
    text, max_iter_count = re.subn(
        r"^(\s*MAX_ITER\s+)\d+\s*$",
        rf"\g<1>{start_step + additional_steps}",
        text,
        count=1,
        flags=re.I | re.M,
    )
    if project_count != 1 or max_iter_count != 1:
        raise ValueError(f"could not update continuation input from {source}")
    return text, start_step


def archive_path(path: Path, label: str) -> Path:
    candidate = path.with_name(f"{path.stem}.{label}{path.suffix}")
    index = 1
    while candidate.exists():
        candidate = path.with_name(f"{path.stem}.{label}.{index}{path.suffix}")
        index += 1
    return candidate


def continue_one(
    input_path: Path,
    cp2k: Path,
    additional_steps: int,
    rounds: int,
) -> tuple[Path, int, str]:
    run_dir = input_path.parent
    canonical_output = run_dir / "cp2k.out"
    with input_path.open() as input_lock:
        try:
            fcntl.flock(input_lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            return input_path, 0, "BUSY"
        if cp2k_completed(canonical_output):
            return input_path, 0, "SKIP"
        canonical_text = canonical_output.read_text(errors="ignore") if canonical_output.exists() else ""
        if "MAXIMUM NUMBER OF OPTIMIZATION STEPS REACHED" not in canonical_text:
            return input_path, 0, "WAIT"

        for round_index in range(1, rounds + 1):
            latest = latest_numbered_restart(run_dir)
            if latest is None:
                return input_path, 1, "NO_RESTART"
            step, restart = latest
            project = f"{input_path.stem}_continue_{step}"
            text, start_step = continuation_input(restart, project, additional_steps)
            if start_step != step:
                raise ValueError(f"restart step mismatch in {restart}: filename={step}, input={start_step}")
            continuation = run_dir / f"{project}.inp"
            continuation.write_text(text)
            output = run_dir / f"cp2k.continue_{step}.out"
            if output.exists():
                output.unlink()
            code = subprocess.run(
                [str(cp2k), "-i", continuation.name, "-o", output.name],
                cwd=run_dir,
                check=False,
            ).returncode
            if code != 0:
                return input_path, code, f"CONTINUE_{round_index}_FAILED"
            result = output.read_text(errors="ignore")
            if "GEOMETRY OPTIMIZATION COMPLETED" in result and "PROGRAM ENDED" in result:
                if canonical_output.exists():
                    canonical_output.replace(archive_path(canonical_output, f"precontinue_{step}"))
                output.replace(canonical_output)
                (run_dir / "returncode.txt").write_text("0\n")
                return input_path, 0, f"CONTINUE_{round_index}"
            if "MAXIMUM NUMBER OF OPTIMIZATION STEPS REACHED" not in result:
                return input_path, 1, f"CONTINUE_{round_index}_INCOMPLETE"
        return input_path, 1, "MAX_ITER"


def continue_maxiter(args: argparse.Namespace) -> None:
    inputs = sorted(args.output_root.resolve().glob(f"*/**/{VARIANT}/*.inp"))
    if args.method:
        inputs = [path for path in inputs if args.method in path.parts]
    if args.system:
        wanted = set(args.system)
        inputs = [path for path in inputs if any(system in path.parts for system in wanted)]
    environment_threads = str(args.threads_per_job)
    os.environ["OMP_NUM_THREADS"] = environment_threads
    os.environ["OPENBLAS_NUM_THREADS"] = environment_threads
    os.environ["MKL_NUM_THREADS"] = environment_threads
    failed = []
    with ThreadPoolExecutor(max_workers=args.jobs) as executor:
        futures = {
            executor.submit(
                continue_one,
                path,
                args.cp2k.resolve(),
                args.additional_steps,
                args.rounds,
            ): path
            for path in inputs
        }
        for future in as_completed(futures):
            input_path, returncode, action = future.result()
            relative = input_path.relative_to(args.output_root.resolve())
            print(f"{action:24s} {relative} rc={returncode}", flush=True)
            if returncode != 0:
                failed.append(relative)
    if failed:
        raise SystemExit(f"{len(failed)} continuation jobs failed")


def polish_one(
    input_path: Path,
    cp2k: Path,
    max_iter: int,
    trust_radius: float,
    force: bool,
) -> tuple[Path, int, str]:
    run_dir = input_path.parent
    canonical_output = run_dir / "cp2k.out"
    with input_path.open() as input_lock:
        try:
            fcntl.flock(input_lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            return input_path, 0, "BUSY"
        if cp2k_completed(canonical_output):
            return input_path, 0, "SKIP_CONVERGED"
        if not cp2k_terminal(canonical_output):
            return input_path, 0, "WAIT"

        restart, record = best_polish_restart(run_dir)
        step = int(record["step"])
        project = f"{input_path.stem}_bfgs_polish_{step}"
        polish_dir = run_dir / f"bfgs_polish_from_{step}"
        if force and polish_dir.exists():
            shutil.rmtree(polish_dir)
        polish_dir.mkdir(parents=True, exist_ok=True)
        polish_input = polish_dir / f"{project}.inp"
        polish_output = polish_dir / "cp2k.out"
        polish_input.write_text(bfgs_polish_input(restart, project, max_iter, trust_radius))

        if force or not cp2k_completed(polish_output):
            code = subprocess.run(
                [str(cp2k), "-i", polish_input.name, "-o", polish_output.name],
                cwd=polish_dir,
                check=False,
            ).returncode
        else:
            code = 0
        if code != 0 or not cp2k_completed(polish_output):
            return input_path, code or 1, "POLISH_FAILED"

        final_restarts = list(polish_dir.glob("*-1.restart"))
        if not final_restarts:
            return input_path, 1, "POLISH_NO_RESTART"
        final_restart = max(final_restarts, key=lambda path: path.stat().st_mtime_ns)
        archived_output = archive_path(canonical_output, f"prebfgs_{step}")
        canonical_output.replace(archived_output)
        shutil.copyfile(polish_output, canonical_output)
        promoted_restart = run_dir / f"{input_path.stem}_bfgs_polished-1.restart"
        shutil.copyfile(final_restart, promoted_restart)
        promoted_restart.touch()
        (run_dir / "returncode.txt").write_text("0\n")
        provenance = {
            "source_restart": str(restart),
            "source_step": step,
            "source_metrics": record,
            "optimizer": "BFGS",
            "trust_radius_angstrom": trust_radius,
            "max_iter": max_iter,
            "polish_input": str(polish_input),
            "polish_output": str(polish_output),
            "promoted_restart": str(promoted_restart),
            "archived_cg_output": str(archived_output),
        }
        (run_dir / "bfgs_polish_provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
        return input_path, 0, f"POLISHED_FROM_{step}"


def polish_bfgs(args: argparse.Namespace) -> None:
    inputs = sorted(args.output_root.resolve().glob(f"*/**/{VARIANT}/*.inp"))
    if args.method:
        inputs = [path for path in inputs if args.method in path.parts]
    if args.system:
        wanted = set(args.system)
        inputs = [path for path in inputs if any(system in path.parts for system in wanted)]
    if not inputs:
        raise ValueError("no prepared cell-optimization inputs selected")
    environment_threads = str(args.threads_per_job)
    os.environ["OMP_NUM_THREADS"] = environment_threads
    os.environ["OPENBLAS_NUM_THREADS"] = environment_threads
    os.environ["MKL_NUM_THREADS"] = environment_threads
    failed = []
    with ThreadPoolExecutor(max_workers=args.jobs) as executor:
        futures = {
            executor.submit(
                polish_one,
                path,
                args.cp2k.resolve(),
                args.max_iter,
                args.trust_radius,
                args.force,
            ): path
            for path in inputs
        }
        for future in as_completed(futures):
            input_path, returncode, action = future.result()
            relative = input_path.relative_to(args.output_root.resolve())
            print(f"{action:24s} {relative} rc={returncode}", flush=True)
            if returncode != 0:
                failed.append(relative)
    if failed:
        raise SystemExit(f"{len(failed)} BFGS polishing jobs failed")


def last_float(text: str, pattern: str) -> float | None:
    matches = re.findall(pattern, text, flags=re.M)
    return float(matches[-1]) if matches else None


def last_int(text: str, pattern: str) -> int | None:
    matches = re.findall(pattern, text, flags=re.M)
    return int(matches[-1]) if matches else None


def format_number(value: float | int | None, digits: int = 6) -> str:
    if value is None:
        return ""
    if isinstance(value, int):
        return str(value)
    return f"{value:.{digits}f}"


def load_gamma_rows(path: Path) -> dict[tuple[str, str], dict[str, str]]:
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    selected = {(row["method"], row["system"]): row for row in rows}
    expected = {(method, str(system["id"])) for system in systems() for method in METHODS}
    missing = expected - set(selected)
    if missing:
        raise ValueError(f"Gamma result table is missing {len(missing)} method/system rows")
    return selected


def load_molecule_rows(run_root: Path) -> dict[tuple[str, str], dict[str, str]]:
    rows: dict[tuple[str, str], dict[str, str]] = {}
    for system_data in systems():
        system = str(system_data["id"])
        for method in METHODS:
            stem = f"{system}_{method}_mol_geoopt"
            output = run_root / "runs" / "molecule_geoopt" / method / stem / f"{stem}.out"
            text = output.read_text(errors="ignore") if output.is_file() else ""
            energy = last_float(text, r"^\s*ENERGY\| Total FORCE_EVAL .*?([-+0-9.Ee]+)\s*$")
            if "GEOMETRY OPTIMIZATION COMPLETED" not in text or "PROGRAM ENDED" not in text or energy is None:
                raise ValueError(f"completed molecule optimization not found for {method}/{system}: {output}")
            rows[(method, system)] = {"gas_energy_hartree": f"{energy:.12f}"}
    return rows


def collect(args: argparse.Namespace) -> None:
    output_root = args.output_root.resolve()
    manifest_path = output_root.parent / "x23b_k222_cellopt_manifest.csv"
    with manifest_path.open(newline="") as handle:
        manifest = {(row["method"], row["system"]): row for row in csv.DictReader(handle)}
    if args.gamma_csv is not None:
        gamma_rows = load_gamma_rows(args.gamma_csv.resolve())
    else:
        gamma_rows = load_molecule_rows(args.molecule_run_root.resolve())
    metadata = {str(row["id"]): row for row in systems()}
    rows = []
    for method in METHODS:
        for system in sorted(metadata):
            source = gamma_rows[(method, system)]
            run_dir = output_root / method / system / VARIANT
            output = run_dir / "cp2k.out"
            text = output.read_text(errors="ignore") if output.is_file() else ""
            energy = last_float(text, r"^\s*ENERGY\| Total FORCE_EVAL .*?([-+0-9.Ee]+)\s*$")
            volume = last_float(text, r"^\s*CELL\| Volume.*?([-+0-9.Ee]+)\s*$")
            step = last_int(text, r"^\s*OPT\| Step number\s+(\d+)\s*$")
            pressure = last_float(text, r"^\s*OPT\| (?:Internal pressure|Pressure deviation) \[?bar\]?\s+([-+0-9.Ee]+)\s*$")
            max_step = last_float(text, r"^\s*OPT\| Maximum step size\s+([-+0-9.Ee]+)\s*$")
            rms_step = last_float(text, r"^\s*OPT\| RMS step size\s+([-+0-9.Ee]+)\s*$")
            max_gradient = last_float(text, r"^\s*OPT\| Maximum gradient\s+([-+0-9.Ee]+)\s*$")
            rms_gradient = last_float(text, r"^\s*OPT\| RMS gradient\s+([-+0-9.Ee]+)\s*$")
            gas_energy = float(source["gas_energy_hartree"])
            ref_energy = float(metadata[system]["ref_energy"])
            ref_volume = float(metadata[system]["x23b_same_cell_ref_volume"])
            n_molecules = int(metadata[system]["molecules_per_cell"])
            lattice = None if energy is None else (gas_energy - energy / n_molecules) * HARTREE_TO_KJMOL
            error = None if lattice is None else lattice - ref_energy
            volume_error = None if volume is None else 100.0 * (volume - ref_volume) / ref_volume
            returncode_file = run_dir / "returncode.txt"
            returncode = int(returncode_file.read_text()) if returncode_file.is_file() else None
            row = {
                "method": method,
                "system": system,
                "variant": VARIANT,
                "mesh": "k222",
                "source_variant": "gamma_cellopt_keep_angles",
                "source": "gamma_cellopt_restart",
                "returncode": "" if returncode is None else returncode,
                "program_ended": "PROGRAM ENDED" in text,
                "opt_completed": "GEOMETRY OPTIMIZATION COMPLETED" in text,
                "max_iter_reached": "MAXIMUM NUMBER OF OPTIMIZATION STEPS REACHED" in text,
                "last_step": format_number(step),
                "energy_hartree": format_number(energy, 12),
                "gas_energy_hartree": format_number(gas_energy, 12),
                "lattice_energy_kJmol": format_number(lattice),
                "x23b_ref_lattice_energy_kJmol": format_number(ref_energy),
                "error_kJmol": format_number(error),
                "volume_A3": format_number(volume),
                "x23b_same_cell_ref_volume_A3": format_number(ref_volume),
                "volume_error_percent": format_number(volume_error),
                "last_pressure_bar": format_number(pressure),
                "last_max_step": format_number(max_step, 10),
                "last_rms_step": format_number(rms_step, 10),
                "last_max_gradient": format_number(max_gradient, 10),
                "last_rms_gradient": format_number(rms_gradient, 10),
                "source_restart": manifest[(method, system)]["source_restart"],
                "run_dir": str(run_dir),
                "output": str(output),
            }
            rows.append(row)
    args.csv.parent.mkdir(parents=True, exist_ok=True)
    with args.csv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    for method in METHODS:
        selected = [row for row in rows if row["method"] == method]
        complete = sum(row["program_ended"] and row["opt_completed"] for row in selected)
        print(f"{method}: {complete}/{len(selected)} converged")


def main() -> None:
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="command", required=True)

    prepare_parser = subparsers.add_parser("prepare")
    prepare_parser.add_argument("--gamma-root", type=Path, required=True)
    prepare_parser.add_argument("--output-root", type=Path, required=True)
    prepare_parser.add_argument("--override", action="append", default=[])
    prepare_parser.add_argument("--method", choices=METHODS)
    prepare_parser.add_argument("--system", action="append", choices=sorted(str(row["id"]) for row in systems()))
    prepare_parser.add_argument("--clean", action="store_true")
    prepare_parser.set_defaults(function=prepare)

    run_parser = subparsers.add_parser("run")
    run_parser.add_argument("--output-root", type=Path, required=True)
    run_parser.add_argument("--cp2k", type=Path, required=True)
    run_parser.add_argument("--jobs", type=int, default=8)
    run_parser.add_argument("--threads-per-job", type=int, default=1)
    run_parser.add_argument("--method", choices=METHODS)
    run_parser.add_argument("--system", action="append", choices=sorted(str(row["id"]) for row in systems()))
    run_parser.add_argument("--start-system", choices=sorted(str(row["id"]) for row in systems()))
    run_parser.add_argument("--end-system", choices=sorted(str(row["id"]) for row in systems()))
    run_parser.add_argument("--force", action="store_true")
    run_parser.set_defaults(function=run)

    continue_parser = subparsers.add_parser("continue-maxiter")
    continue_parser.add_argument("--output-root", type=Path, required=True)
    continue_parser.add_argument("--cp2k", type=Path, required=True)
    continue_parser.add_argument("--jobs", type=int, default=2)
    continue_parser.add_argument("--threads-per-job", type=int, default=1)
    continue_parser.add_argument("--additional-steps", type=int, default=300)
    continue_parser.add_argument("--rounds", type=int, default=3)
    continue_parser.add_argument("--method", choices=METHODS)
    continue_parser.add_argument("--system", action="append")
    continue_parser.set_defaults(function=continue_maxiter)

    polish_parser = subparsers.add_parser("polish-bfgs")
    polish_parser.add_argument("--output-root", type=Path, required=True)
    polish_parser.add_argument("--cp2k", type=Path, required=True)
    polish_parser.add_argument("--jobs", type=int, default=2)
    polish_parser.add_argument("--threads-per-job", type=int, default=1)
    polish_parser.add_argument("--max-iter", type=int, default=300)
    polish_parser.add_argument("--trust-radius", type=float, default=0.002)
    polish_parser.add_argument("--method", choices=METHODS)
    polish_parser.add_argument("--system", action="append")
    polish_parser.add_argument("--force", action="store_true")
    polish_parser.set_defaults(function=polish_bfgs)

    collect_parser = subparsers.add_parser("collect")
    collect_parser.add_argument("--output-root", type=Path, required=True)
    energy_source = collect_parser.add_mutually_exclusive_group(required=True)
    energy_source.add_argument("--gamma-csv", type=Path)
    energy_source.add_argument("--molecule-run-root", type=Path)
    collect_parser.add_argument("--csv", type=Path, required=True)
    collect_parser.set_defaults(function=collect)

    args = parser.parse_args()
    args.function(args)


if __name__ == "__main__":
    main()
