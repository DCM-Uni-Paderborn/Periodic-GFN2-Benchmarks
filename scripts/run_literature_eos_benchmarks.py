#!/usr/bin/env python3
"""Run and analyse the Klimes Solids23 EOS benchmark.

The periodic calculations use CP2K's native Bloch k-point implementation and
the linked tblite library.  Every point is independent, restartable, and kept
with its generated CP2K input and compact output for provenance.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import re
import shutil
import subprocess
import sys
import time
from collections import Counter
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np


REPO = Path(__file__).resolve().parents[1]
CP2K_DEFAULT = Path("/Users/tkuehne/gxtb-local-build/install/cp2k/bin/cp2k.ssmp")
TBLITE_DEFAULT = Path("/Users/tkuehne/gxtb-local-build/install/tblite/bin/tblite")
HARTREE_TO_EV = 27.211386245988
HARTREE_PER_A3_TO_GPA = 29421.02648438959
DEFAULT_VOLUME_FACTORS = tuple(round(x, 2) for x in np.linspace(0.8, 1.2, 11))
METHODS = ("GFN1", "GFN2")


@dataclass(frozen=True)
class Structure:
    system: str
    label: str
    category: str
    cell: np.ndarray
    symbols: tuple[str, ...]
    fractional: np.ndarray
    kmesh: tuple[int, int, int]
    reference_volume: float

    @property
    def natoms(self) -> int:
        return len(self.symbols)

    @property
    def volume(self) -> float:
        return abs(float(np.linalg.det(self.cell)))


@dataclass(frozen=True)
class Job:
    dataset: str
    structure: Structure
    method: str
    volume_factor: float
    cp2k: Path
    threads: int
    full_grid: bool

    @property
    def root(self) -> Path:
        return dataset_root(self.dataset) / "results/raw" / self.method / self.structure.system / f"v_{self.volume_factor:.4f}"


def dataset_root(dataset: str) -> Path:
    if dataset == "klimes":
        return REPO / "Klimes-Solids23"
    raise ValueError(f"Unknown dataset: {dataset}")


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def fcc_sites() -> list[tuple[float, float, float]]:
    return [(0, 0, 0), (0, 0.5, 0.5), (0.5, 0, 0.5), (0.5, 0.5, 0)]


def klimes_structure(row: dict[str, str]) -> Structure:
    system = row["system"]
    a = float(row["experiment_ZPEC"])
    structure = row["structure"]
    cell = np.eye(3) * a
    species: dict[str, tuple[str, ...]] = {
        "LiF": ("Li", "F"),
        "LiCl": ("Li", "Cl"),
        "NaF": ("Na", "F"),
        "NaCl": ("Na", "Cl"),
        "MgO": ("Mg", "O"),
        "SiC": ("Si", "C"),
        "GaAs": ("Ga", "As"),
    }
    if structure == "fcc":
        fractional = fcc_sites()
        symbols = [system] * 4
    elif structure == "bcc":
        fractional = [(0, 0, 0), (0.5, 0.5, 0.5)]
        symbols = [system] * 2
    elif structure == "rocksalt":
        first = fcc_sites()
        second = [tuple((np.asarray(site) + (0.5, 0, 0)) % 1.0) for site in first]
        fractional = first + second
        symbols = [species[system][0]] * 4 + [species[system][1]] * 4
    elif structure in {"diamond", "zincblende"}:
        first = fcc_sites()
        second = [tuple((np.asarray(site) + 0.25) % 1.0) for site in first]
        fractional = first + second
        if structure == "diamond":
            symbols = [system] * 8
        else:
            symbols = [species[system][0]] * 4 + [species[system][1]] * 4
    else:
        raise ValueError(f"Unsupported Klimes structure: {structure}")
    metallic = row["category"] in {
        "transition_metal",
        "alkali_metal",
        "alkaline_earth_metal",
        "simple_metal",
    }
    k = 16 if metallic else 8
    return Structure(
        system=system,
        label=system,
        category=row["category"],
        cell=cell,
        symbols=tuple(symbols),
        fractional=np.asarray(fractional, dtype=float),
        kmesh=(k, k, k),
        reference_volume=a**3,
    )


def load_structures(dataset: str) -> list[Structure]:
    if dataset == "klimes":
        rows = read_csv(REPO / "Klimes-Solids23/references/lattice_constants.csv")
        return [klimes_structure(row) for row in rows]
    raise ValueError(f"Unknown dataset: {dataset}")


def scale_structure(structure: Structure, volume_factor: float) -> tuple[np.ndarray, np.ndarray]:
    linear_factor = volume_factor ** (1.0 / 3.0)
    return structure.cell * linear_factor, structure.fractional.copy()


def safe_project(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_-]", "_", value)[:72]


def cp2k_input(
    job: Job,
    damping: float,
    max_scf: int,
    mixer_memory: int,
    scc_mixer: str = "TBLITE",
) -> str:
    structure = job.structure
    cell, fractional = scale_structure(structure, job.volume_factor)
    project = safe_project(f"{job.dataset}_{structure.system}_{job.method}_v{job.volume_factor:.4f}")
    vectors = "\n".join(
        f"      {name} {vec[0]:.14f} {vec[1]:.14f} {vec[2]:.14f}"
        for name, vec in zip("ABC", cell)
    )
    coordinates = "\n".join(
        f"      {symbol:2s} {xyz[0]:.14f} {xyz[1]:.14f} {xyz[2]:.14f}"
        for symbol, xyz in zip(structure.symbols, fractional)
    )
    kx, ky, kz = structure.kmesh
    if job.full_grid:
        reduction = "      FULL_GRID T"
    else:
        reduction = "      FULL_GRID F\n      SYMMETRY F\n      INVERSION_SYMMETRY_ONLY T"
    scc_mixer = scc_mixer.upper()
    if scc_mixer == "TBLITE":
        xtb_mixer = f"""        &TBLITE_MIXER
          ITERATIONS {max_scf}
          MEMORY {mixer_memory}
          DAMPING {damping:.8f}
        &END TBLITE_MIXER"""
        scf_mixing = ""
    elif scc_mixer == "CP2K":
        xtb_mixer = ""
        scf_mixing = f"""      EPS_DIIS 1.0E-12
      &MIXING
        METHOD DIRECT_P_MIXING
        ALPHA {damping:.8f}
      &END MIXING
"""
    else:
        raise ValueError(f"unsupported SCC mixer: {scc_mixer}")
    return f"""&GLOBAL
  PROJECT {project}
  RUN_TYPE ENERGY
  PRINT_LEVEL LOW
&END GLOBAL

&FORCE_EVAL
  METHOD Quickstep
  &DFT
    &QS
      EPS_DEFAULT 1.0E-12
      METHOD xTB
      &XTB
        GFN_TYPE TBLITE
        SCC_MIXER {scc_mixer}
        &TBLITE
          METHOD {job.method}
          ACCURACY 0.05
        &END TBLITE
{xtb_mixer}
      &END XTB
    &END QS
    &KPOINTS
      SCHEME MACDONALD {kx} {ky} {kz} 0.0 0.0 0.0
{reduction}
    &END KPOINTS
    &SCF
      EPS_SCF 1.0E-9
      MAX_SCF {max_scf}
      SCF_GUESS MOPAC
      ADDED_MOS -1 -1
      &SMEAR ON
        METHOD FERMI_DIRAC
        ELECTRONIC_TEMPERATURE [K] 300.0
      &END SMEAR
{scf_mixing}      &PRINT
        &RESTART OFF
        &END RESTART
      &END PRINT
    &END SCF
  &END DFT
  &SUBSYS
    &CELL
      PERIODIC XYZ
{vectors}
    &END CELL
    &COORD
      SCALED
{coordinates}
    &END COORD
  &END SUBSYS
&END FORCE_EVAL
"""


def parse_cp2k_output(path: Path) -> tuple[bool, float | None, str]:
    if not path.exists():
        return False, None, "missing output"
    text = path.read_text(errors="replace")
    energy = None
    extrapolated_energy = None
    for line in text.splitlines():
        if "Total energy (extrapolated to T->0)" in line:
            try:
                extrapolated_energy = float(line.split()[-1])
            except ValueError:
                pass
        if "ENERGY| Total FORCE_EVAL" in line:
            try:
                energy = float(line.split()[-1])
            except ValueError:
                pass
    failures = ("ABORT", "SCF run NOT converged", "DID NOT CONVERGE")
    if any(token in text for token in failures):
        return False, extrapolated_energy if extrapolated_energy is not None else energy, "CP2K reported a convergence failure or abort"
    if "PROGRAM ENDED" not in text:
        return False, extrapolated_energy if extrapolated_energy is not None else energy, "CP2K did not reach PROGRAM ENDED"
    if energy is None:
        return False, None, "energy not found"
    return True, extrapolated_energy if extrapolated_energy is not None else energy, "ok"


def result_path(job: Job) -> Path:
    return job.root / "result.json"


def completed_job(job: Job) -> bool:
    path = result_path(job)
    if not path.exists():
        return False
    try:
        result = json.loads(path.read_text())
    except (json.JSONDecodeError, OSError):
        return False
    return bool(result.get("ok")) and result.get("energy_hartree") is not None


def cleanup_cp2k_scratch(root: Path) -> None:
    patterns = ("*.wfn*", "*.restart*", "*.ener", "*.BFGS", "mainLog.out")
    for pattern in patterns:
        for path in root.glob(pattern):
            if path.is_file():
                path.unlink()


def run_job(job: Job, start_attempt: int = 1, force: bool = False) -> dict[str, object]:
    root = job.root
    root.mkdir(parents=True, exist_ok=True)
    if completed_job(job) and not force:
        return json.loads(result_path(job).read_text())
    attempts = [
        ("TBLITE", 0.40, 500, 250),
        ("TBLITE", 0.20, 750, 100),
        ("TBLITE", 0.10, 1200, 50),
        ("TBLITE", 0.05, 1800, 25),
        ("CP2K", 0.10, 1500, 0),
        ("CP2K", 0.05, 2000, 0),
        ("CP2K", 0.02, 3000, 0),
    ]
    if not 1 <= start_attempt <= len(attempts):
        raise ValueError(f"start_attempt must be between 1 and {len(attempts)}")
    started = time.time()
    last_reason = "not run"
    last_energy: float | None = None
    return_code: int | None = None
    used_damping: float | None = None
    used_max_scf: int | None = None
    used_mixer_memory: int | None = None
    used_scc_mixer: str | None = None
    attempts_executed = 0
    for attempt_index in range(start_attempt - 1, len(attempts)):
        attempt = attempt_index + 1
        scc_mixer, damping, max_scf, mixer_memory = attempts[attempt_index]
        attempts_executed += 1
        used_damping = damping
        used_max_scf = max_scf
        used_mixer_memory = mixer_memory
        used_scc_mixer = scc_mixer
        input_path = root / "cp2k.inp"
        output_path = root / "cp2k.out"
        input_path.write_text(cp2k_input(job, damping, max_scf, mixer_memory, scc_mixer))
        if output_path.exists():
            output_path.unlink()
        env = os.environ.copy()
        env.update(
            {
                "OMP_NUM_THREADS": str(job.threads),
                "OMP_PROC_BIND": "false",
                "OPENBLAS_NUM_THREADS": "1",
                "MKL_NUM_THREADS": "1",
                "VECLIB_MAXIMUM_THREADS": "1",
            }
        )
        with (root / "launcher.log").open("w") as launcher:
            proc = subprocess.run(
                [str(job.cp2k), "-i", input_path.name, "-o", output_path.name],
                cwd=root,
                env=env,
                stdout=launcher,
                stderr=subprocess.STDOUT,
                check=False,
            )
        return_code = proc.returncode
        ok, last_energy, last_reason = parse_cp2k_output(output_path)
        cleanup_cp2k_scratch(root)
        if ok and return_code == 0:
            break
    else:
        ok = False
        attempt = len(attempts)
    cell, _ = scale_structure(job.structure, job.volume_factor)
    result: dict[str, object] = {
        "ok": ok,
        "reason": last_reason,
        "return_code": return_code,
        "attempts": attempts_executed,
        "mixer_profile_stage": attempt,
        "start_attempt": start_attempt,
        "damping": used_damping,
        "max_scf": used_max_scf,
        "tblite_mixer_damping": used_damping,
        "tblite_mixer_iterations": used_max_scf,
        "tblite_mixer_memory": used_mixer_memory if used_scc_mixer == "TBLITE" else None,
        "cp2k_mixer_alpha": used_damping if used_scc_mixer == "CP2K" else None,
        "scc_mixer": used_scc_mixer,
        "dataset": job.dataset,
        "system": job.structure.system,
        "label": job.structure.label,
        "category": job.structure.category,
        "method": job.method,
        "volume_factor": job.volume_factor,
        "volume_A3": abs(float(np.linalg.det(cell))),
        "energy_hartree": last_energy,
        "natoms": job.structure.natoms,
        "kmesh": list(job.structure.kmesh),
        "kpoint_path": "CP2K native Bloch",
        "full_grid": job.full_grid,
        "elapsed_seconds": time.time() - started,
    }
    result_path(job).write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    return result


def selected_structures(dataset: str, systems: set[str] | None) -> list[Structure]:
    structures = load_structures(dataset)
    if systems is None:
        return structures
    selected = [structure for structure in structures if structure.system in systems]
    missing = systems - {structure.system for structure in selected}
    if missing:
        raise ValueError(f"Unknown systems in {dataset}: {', '.join(sorted(missing))}")
    return selected


def make_jobs(
    datasets: Iterable[str],
    methods: Iterable[str],
    systems: set[str] | None,
    factors: Iterable[float],
    cp2k: Path,
    threads: int,
    full_grid: bool,
) -> list[Job]:
    jobs: list[Job] = []
    for dataset in datasets:
        for structure in selected_structures(dataset, systems):
            for method in methods:
                for factor in factors:
                    jobs.append(Job(dataset, structure, method, factor, cp2k, threads, full_grid))
    return jobs


def run_jobs(args: argparse.Namespace) -> int:
    if not args.cp2k.is_file():
        raise FileNotFoundError(args.cp2k)
    jobs = make_jobs(
        args.datasets,
        args.methods,
        set(args.systems) if args.systems else None,
        args.volume_factors,
        args.cp2k,
        args.threads,
        args.full_grid,
    )
    pending = jobs if args.force else [job for job in jobs if not completed_job(job)]
    print(f"EOS jobs: {len(jobs)} total, {len(jobs) - len(pending)} complete, {len(pending)} pending", flush=True)
    if not pending:
        return 0
    failures = 0
    completed = len(jobs) - len(pending)
    started = time.time()
    with ThreadPoolExecutor(max_workers=args.jobs) as executor:
        futures = {executor.submit(run_job, job, args.start_attempt, args.force): job for job in pending}
        for future in as_completed(futures):
            job = futures[future]
            try:
                result = future.result()
            except Exception as exc:  # keep the other independent jobs running
                failures += 1
                print(f"FAIL {job.dataset}/{job.method}/{job.structure.system}/{job.volume_factor:.2f}: {exc}", flush=True)
                continue
            completed += 1
            if not result["ok"]:
                failures += 1
            elapsed = time.time() - started
            rate = completed / max(elapsed, 1.0)
            remaining = (len(jobs) - completed) / rate if rate else math.nan
            state = "OK" if result["ok"] else "FAIL"
            print(
                f"[{completed:4d}/{len(jobs):4d}] {state:4s} {job.dataset:14s} {job.method} "
                f"{job.structure.system:18s} v={job.volume_factor:.2f} "
                f"({float(result['elapsed_seconds']):.1f}s; ETA {remaining / 60:.1f} min)",
                flush=True,
            )
    print(f"Finished with {failures} failed jobs", flush=True)
    return 1 if failures else 0


def load_result(job: Job) -> dict[str, object] | None:
    path = result_path(job)
    if not path.exists():
        return None
    try:
        return json.loads(path.read_text())
    except json.JSONDecodeError:
        return None


def status_jobs(args: argparse.Namespace) -> int:
    jobs = make_jobs(
        args.datasets,
        args.methods,
        set(args.systems) if args.systems else None,
        args.volume_factors,
        args.cp2k,
        args.threads,
        args.full_grid,
    )
    states = Counter()
    failures: list[str] = []
    seconds = 0.0
    for job in jobs:
        result = load_result(job)
        if result is None:
            states["pending"] += 1
        elif result.get("ok"):
            states["complete"] += 1
            seconds += float(result.get("elapsed_seconds", 0.0))
        else:
            states["failed"] += 1
            failures.append(f"{job.dataset}/{job.method}/{job.structure.system}/v={job.volume_factor:.2f}: {result.get('reason')}")
    print(json.dumps({"total": len(jobs), **states, "accumulated_job_hours": seconds / 3600.0}, indent=2))
    for failure in failures:
        print(f"FAIL {failure}")
    return 1 if failures else 0


def murnaghan_feature(volumes: np.ndarray, v0: float, b1: float) -> np.ndarray:
    """Return the coefficient multiplying B0 in the Murnaghan energy EOS."""
    volumes = np.asarray(volumes, dtype=float)
    log_ratio = np.log(v0 / volumes)
    if abs(b1) < 1.0e-8:
        return v0 - volumes - volumes * log_ratio
    if abs(b1 - 1.0) < 1.0e-8:
        return v0 * log_ratio + volumes - v0
    numerator = volumes * np.expm1(b1 * log_ratio) + b1 * (volumes - v0)
    return numerator / (b1 * (b1 - 1.0))


def murnaghan_feature_grid(volumes: np.ndarray, v0: float, b1: np.ndarray) -> np.ndarray:
    volumes = np.asarray(volumes, dtype=float)
    b1 = np.asarray(b1, dtype=float)
    values = np.empty((len(b1), len(volumes)), dtype=float)
    zero = np.abs(b1) < 1.0e-8
    one = np.abs(b1 - 1.0) < 1.0e-8
    regular = ~(zero | one)
    log_ratio = np.log(v0 / volumes)
    if np.any(regular):
        b = b1[regular, None]
        numerator = volumes[None, :] * np.expm1(b * log_ratio[None, :]) + b * (volumes[None, :] - v0)
        values[regular] = numerator / (b * (b - 1.0))
    if np.any(zero):
        values[zero] = v0 - volumes - volumes * log_ratio
    if np.any(one):
        values[one] = v0 * log_ratio + volumes - v0
    return values


def fit_murnaghan_at_grid(
    volumes: np.ndarray,
    energies: np.ndarray,
    v0: float,
    b1_values: np.ndarray,
) -> tuple[float, float, float, float]:
    features = murnaghan_feature_grid(volumes, v0, b1_values)
    centered_features = features - features.mean(axis=1, keepdims=True)
    centered_energies = energies - energies.mean()
    denominator = np.sum(centered_features**2, axis=1)
    with np.errstate(divide="ignore", invalid="ignore"):
        b0 = np.sum(centered_features * centered_energies[None, :], axis=1) / denominator
    e0 = energies.mean() - b0 * features.mean(axis=1)
    residual = energies[None, :] - (e0[:, None] + b0[:, None] * features)
    scores = np.sum(residual**2, axis=1)
    scores[(b0 <= 0.0) | ~np.isfinite(scores)] = math.inf
    index = int(np.argmin(scores))
    return float(scores[index]), float(e0[index]), float(b0[index]), float(b1_values[index])


def fit_murnaghan(volumes: np.ndarray, energies: np.ndarray) -> dict[str, float | bool]:
    order = np.argsort(volumes)
    volumes = np.asarray(volumes, dtype=float)[order]
    energies = np.asarray(energies, dtype=float)[order]
    full_v0_lower, full_v0_upper = 0.92 * volumes.min(), 1.08 * volumes.max()
    full_b1_lower, full_b1_upper = -20.0, 30.0
    v0_lower, v0_upper = full_v0_lower, full_v0_upper
    b1_lower, b1_upper = full_b1_lower, full_b1_upper
    best_v0 = float(volumes[np.argmin(energies)])
    best_b1 = 4.0
    best_e0 = float(energies.min())
    best_b0 = math.nan
    best_score = math.inf
    for _ in range(6):
        v0_grid = np.linspace(v0_lower, v0_upper, 41)
        b1_grid = np.linspace(b1_lower, b1_upper, 51)
        for candidate_v0 in v0_grid:
            score, e0, b0, b1 = fit_murnaghan_at_grid(volumes, energies, float(candidate_v0), b1_grid)
            if score < best_score:
                best_score = score
                best_v0 = float(candidate_v0)
                best_e0 = e0
                best_b0 = b0
                best_b1 = b1
        v0_step = float(v0_grid[1] - v0_grid[0])
        b1_step = float(b1_grid[1] - b1_grid[0])
        v0_lower = max(full_v0_lower, best_v0 - 2.0 * v0_step)
        v0_upper = min(full_v0_upper, best_v0 + 2.0 * v0_step)
        b1_lower = max(full_b1_lower, best_b1 - 2.0 * b1_step)
        b1_upper = min(full_b1_upper, best_b1 + 2.0 * b1_step)
    predicted = best_e0 + best_b0 * murnaghan_feature(volumes, best_v0, best_b1)
    residual = energies - predicted
    rmse = float(np.sqrt(np.mean(residual**2)))
    total = float(np.sum((energies - energies.mean()) ** 2))
    r_squared = 1.0 - float(np.sum(residual**2)) / total if total > 0 else 1.0
    return {
        "E0_hartree": best_e0,
        "V0_A3": best_v0,
        "B0_hartree_A3": best_b0,
        "B0_GPa": best_b0 * HARTREE_PER_A3_TO_GPA,
        "B1": best_b1,
        "rmse_hartree": rmse,
        "r_squared": r_squared,
        "minimum_inside_sample": bool(volumes.min() < best_v0 < volumes.max()),
        "sample_minimum_at_edge": bool(np.argmin(energies) in {0, len(energies) - 1}),
        "fit_score": best_score,
    }


def write_dict_csv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("")
        return
    fields: list[str] = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def available_result_factors(dataset: str, structure: Structure, method: str, requested: Iterable[float]) -> list[float]:
    factors = {float(factor) for factor in requested}
    root = dataset_root(dataset) / "results/raw" / method / structure.system
    for result_path in root.glob("v_*/result.json"):
        try:
            factors.add(float(result_path.parent.name.removeprefix("v_")))
        except ValueError:
            continue
    return sorted(factors)


def select_fit_window(points: list[dict[str, object]], max_points: int = 11) -> list[dict[str, object]]:
    ordered = sorted(points, key=lambda row: float(row["volume_A3"]))
    if len(ordered) <= max_points:
        return ordered
    energies = np.asarray([float(row["energy_hartree"]) for row in ordered])
    minimum = int(np.argmin(energies))
    start = max(0, minimum - max_points // 2)
    start = min(start, len(ordered) - max_points)
    return ordered[start : start + max_points]


def select_eos_fit(
    points: list[dict[str, object]],
    min_points: int = 7,
    max_points: int = 11,
) -> tuple[list[dict[str, object]], dict[str, float | bool]]:
    ordered = sorted(points, key=lambda row: float(row["volume_A3"]))
    energies = np.asarray([float(row["energy_hartree"]) for row in ordered])
    factors = np.asarray([float(row["volume_factor"]) for row in ordered])
    local_minima = [
        index
        for index in range(1, len(ordered) - 1)
        if energies[index] <= energies[index - 1] and energies[index] <= energies[index + 1]
    ]
    local_minima.sort(key=lambda index: abs(factors[index] - 1.0))
    for minimum in local_minima:
        for size in range(min(max_points, len(ordered)), min_points - 1, -1):
            start_min = max(0, minimum - size + 2)
            start_max = min(minimum - 1, len(ordered) - size)
            candidates: list[
                tuple[float, float, float, list[dict[str, object]], dict[str, float | bool]]
            ] = []
            for start in range(start_min, start_max + 1):
                window = ordered[start : start + size]
                volumes = np.asarray([float(row["volume_A3"]) for row in window])
                window_energies = np.asarray([float(row["energy_hartree"]) for row in window])
                fit = fit_murnaghan(volumes, window_energies)
                if eos_fit_reason(fit) != "ok":
                    continue
                midpoint = 0.5 * (float(window[0]["volume_factor"]) + float(window[-1]["volume_factor"]))
                span = float(window[-1]["volume_factor"]) - float(window[0]["volume_factor"])
                candidates.append(
                    (
                        span,
                        abs(midpoint - factors[minimum]),
                        float(fit["rmse_hartree"]),
                        window,
                        fit,
                    )
                )
            if candidates:
                _, _, _, window, fit = min(candidates, key=lambda item: (item[0], item[1], item[2]))
                return window, fit
    fallback = select_fit_window(ordered, max_points=max_points)
    volumes = np.asarray([float(row["volume_A3"]) for row in fallback])
    fallback_energies = np.asarray([float(row["energy_hartree"]) for row in fallback])
    return fallback, fit_murnaghan(volumes, fallback_energies)


def eos_fit_reason(fit: dict[str, float | bool]) -> str:
    if not fit["minimum_inside_sample"] or fit["sample_minimum_at_edge"]:
        return "EOS minimum not bracketed by sampled volumes"
    if not 0.0 < float(fit["B0_GPa"]) < 5000.0:
        return "nonphysical fitted EOS curvature"
    if not -20.0 < float(fit["B1"]) < 30.0:
        return "EOS pressure derivative hit fit boundary"
    if float(fit["r_squared"]) < 0.999:
        return "EOS fit R-squared below 0.999"
    return "ok"


def collect_and_fit(
    dataset: str,
    methods: list[str],
    factors: list[float],
    systems: set[str] | None = None,
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    points: list[dict[str, object]] = []
    fits: list[dict[str, object]] = []
    dummy_cp2k = CP2K_DEFAULT
    for structure in selected_structures(dataset, systems):
        for method in methods:
            method_points = []
            for factor in available_result_factors(dataset, structure, method, factors):
                job = Job(dataset, structure, method, factor, dummy_cp2k, 1, True)
                result = load_result(job)
                if result is None or not result.get("ok"):
                    continue
                row = {
                    "dataset": dataset,
                    "system": structure.system,
                    "label": structure.label,
                    "category": structure.category,
                    "method": method,
                    "volume_factor": factor,
                    "volume_A3": float(result["volume_A3"]),
                    "energy_hartree": float(result["energy_hartree"]),
                    "natoms": structure.natoms,
                    "kx": structure.kmesh[0],
                    "ky": structure.kmesh[1],
                    "kz": structure.kmesh[2],
                }
                points.append(row)
                method_points.append(row)
            if len(method_points) < 7:
                fits.append(
                    {
                        "dataset": dataset,
                        "system": structure.system,
                        "label": structure.label,
                        "category": structure.category,
                        "method": method,
                        "npoints": len(method_points),
                        "npoints_available": len(method_points),
                        "fit_ok": False,
                        "fit_reason": "fewer than 7 converged EOS points",
                    }
                )
                continue
            fit_points, fit = select_eos_fit(method_points)
            fit_reason = eos_fit_reason(fit)
            fit_ok = fit_reason == "ok"
            fits.append(
                {
                    "dataset": dataset,
                    "system": structure.system,
                    "label": structure.label,
                    "category": structure.category,
                    "method": method,
                    "npoints": len(fit_points),
                    "npoints_available": len(method_points),
                    "fit_volume_factor_min": min(float(row["volume_factor"]) for row in fit_points),
                    "fit_volume_factor_max": max(float(row["volume_factor"]) for row in fit_points),
                    "sample_minimum_volume_factor": min(
                        fit_points, key=lambda row: float(row["energy_hartree"])
                    )["volume_factor"],
                    "natoms": structure.natoms,
                    "reference_volume_A3": structure.reference_volume,
                    "eos_model": "Murnaghan",
                    "fit_ok": fit_ok,
                    "fit_reason": fit_reason,
                    "rmse_meV_atom": float(fit["rmse_hartree"]) * HARTREE_TO_EV * 1000.0 / structure.natoms,
                    **fit,
                }
            )
    return points, fits


def error_stats(values: list[float], references: list[float]) -> dict[str, float | int]:
    value = np.asarray(values, dtype=float)
    reference = np.asarray(references, dtype=float)
    error = value - reference
    relative = 100.0 * error / reference
    return {
        "n": len(value),
        "ME": float(error.mean()),
        "MAE": float(np.abs(error).mean()),
        "RMSE": float(np.sqrt(np.mean(error**2))),
        "MRE_percent": float(relative.mean()),
        "MARE_percent": float(np.abs(relative).mean()),
    }


def atomic_energies_path() -> Path:
    return REPO / "Klimes-Solids23/results/atomic_energies.json"


def analyse_klimes(fits: list[dict[str, object]]) -> None:
    lattice_rows = read_csv(REPO / "Klimes-Solids23/references/lattice_constants.csv")
    bulk_rows = {row["system"]: row for row in read_csv(REPO / "Klimes-Solids23/references/bulk_moduli.csv")}
    cohesive_rows = {row["system"]: row for row in read_csv(REPO / "Klimes-Solids23/references/cohesive_energies.csv")}
    fit_map = {(str(row["method"]), str(row["system"])): row for row in fits if row.get("fit_ok")}
    atom_data = json.loads(atomic_energies_path().read_text()) if atomic_energies_path().exists() else {}
    comparison: list[dict[str, object]] = []
    for reference in lattice_rows:
        system = reference["system"]
        for method in METHODS:
            fit = fit_map.get((method, system))
            if fit is None:
                continue
            a0 = float(fit["V0_A3"]) ** (1.0 / 3.0)
            cohesive = None
            if method in atom_data:
                structure = next(s for s in load_structures("klimes") if s.system == system)
                counts = Counter(structure.symbols)
                atomic_sum = sum(count * float(atom_data[method][symbol]["energy_hartree"]) for symbol, count in counts.items())
                cohesive = (atomic_sum - float(fit["E0_hartree"])) * HARTREE_TO_EV / structure.natoms
            row: dict[str, object] = {
                "system": system,
                "category": reference["category"],
                "structure": reference["structure"],
                "method": method,
                "a0_A": a0,
                "a_ref_A": float(reference["experiment_ZPEC"]),
                "a_error_A": a0 - float(reference["experiment_ZPEC"]),
                "B0_GPa": float(fit["B0_GPa"]),
                "B_ref_GPa": float(bulk_rows[system]["experiment"]),
                "B_error_GPa": float(fit["B0_GPa"]) - float(bulk_rows[system]["experiment"]),
                "cohesive_eV_atom": cohesive,
                "cohesive_ref_eV_atom": float(cohesive_rows[system]["experiment_ZPEC"]),
                "cohesive_error_eV_atom": None if cohesive is None else cohesive - float(cohesive_rows[system]["experiment_ZPEC"]),
                "fit_rmse_meV_atom": fit["rmse_meV_atom"],
            }
            comparison.append(row)
    root = REPO / "Klimes-Solids23/results"
    write_dict_csv(root / "comparison.csv", comparison)
    aggregate: list[dict[str, object]] = []
    property_specs = [
        ("lattice_constant_A", "a0_A", "a_ref_A"),
        ("bulk_modulus_GPa", "B0_GPa", "B_ref_GPa"),
        ("cohesive_energy_eV_atom", "cohesive_eV_atom", "cohesive_ref_eV_atom"),
    ]
    for property_name, value_key, ref_key in property_specs:
        for method in METHODS:
            rows = [row for row in comparison if row["method"] == method and row[value_key] is not None]
            if rows:
                aggregate.append(
                    {"property": property_name, "method": method, **error_stats([float(row[value_key]) for row in rows], [float(row[ref_key]) for row in rows])}
                )
    literature = ["revPBE_vdW", "rPW86_vdW2", "optPBE_vdW", "optB88_vdW", "optB86b_vdW", "LDA", "PBEsol", "PBE"]
    for method in literature:
        aggregate.append(
            {
                "property": "lattice_constant_A",
                "method": method,
                **error_stats([float(row[method]) for row in lattice_rows], [float(row["experiment_ZPEC"]) for row in lattice_rows]),
            }
        )
        bulk_list = list(bulk_rows.values())
        aggregate.append(
            {
                "property": "bulk_modulus_GPa",
                "method": method,
                **error_stats([float(row[method]) for row in bulk_list], [float(row["experiment"]) for row in bulk_list]),
            }
        )
        cohesive_list = list(cohesive_rows.values())
        aggregate.append(
            {
                "property": "cohesive_energy_eV_atom",
                "method": method,
                **error_stats([float(row[method]) for row in cohesive_list], [float(row["experiment_ZPEC"]) for row in cohesive_list]),
            }
        )
    for property_name, value_key, ref_key in property_specs:
        for method in METHODS:
            for category in sorted({str(row["category"]) for row in comparison}):
                rows = [
                    row
                    for row in comparison
                    if row["method"] == method and row["category"] == category and row[value_key] is not None
                ]
                if rows:
                    aggregate.append(
                        {
                            "property": property_name,
                            "method": f"{method}:{category}",
                            **error_stats([float(row[value_key]) for row in rows], [float(row[ref_key]) for row in rows]),
                        }
                    )
    write_dict_csv(root / "aggregate_statistics.csv", aggregate)


def analyse(args: argparse.Namespace) -> int:
    failed_fits = 0
    systems = set(args.systems) if args.systems else None
    for dataset in args.datasets:
        points, fits = collect_and_fit(dataset, args.methods, args.volume_factors, systems)
        root = dataset_root(dataset) / "results"
        write_dict_csv(root / "eos_points.csv", points)
        write_dict_csv(root / "eos_fits.csv", fits)
        failed = [row for row in fits if not row.get("fit_ok")]
        failed_fits += len(failed)
        print(f"{dataset}: {len(fits) - len(failed)}/{len(fits)} EOS fits accepted")
        for row in failed:
            print(f"  FIT FAIL {row['method']} {row['system']}: {row['fit_reason']}")
        analyse_klimes(fits)
    return 1 if failed_fits else 0


ATOMIC_SPINS = {
    "Ag": 1,
    "Al": 1,
    "As": 3,
    "Ba": 0,
    "C": 2,
    "Ca": 0,
    "Cl": 1,
    "Cs": 1,
    "Cu": 1,
    "F": 1,
    "Ga": 1,
    "Ge": 2,
    "K": 1,
    "Li": 1,
    "Mg": 0,
    "Na": 1,
    "O": 2,
    "Pd": 0,
    "Rb": 1,
    "Rh": 3,
    "Si": 2,
    "Sr": 0,
}


def run_atoms(args: argparse.Namespace) -> int:
    if not args.tblite.is_file():
        raise FileNotFoundError(args.tblite)
    structures = load_structures("klimes")
    elements = sorted({symbol for structure in structures for symbol in structure.symbols})
    root = REPO / "Klimes-Solids23/results/atoms"
    data: dict[str, dict[str, object]] = {}
    failures = 0
    for method in args.methods:
        data[method] = {}
        for element in elements:
            run_root = root / method / element
            run_root.mkdir(parents=True, exist_ok=True)
            xyz = run_root / f"{element}.xyz"
            xyz.write_text(f"1\n{element} atom\n{element} 0.0 0.0 0.0\n")
            output_json = run_root / "tblite.json"
            stdout = run_root / "tblite.out"
            command = [
                str(args.tblite),
                "run",
                "--method",
                method.lower(),
                "--spin",
                str(ATOMIC_SPINS[element]),
                "--acc",
                "0.05",
                "--iterations",
                "1000",
                "--no-restart",
                "--json",
                str(output_json),
                str(xyz),
            ]
            if ATOMIC_SPINS[element] > 0:
                command.insert(-3, "--spin-polarized")
            with stdout.open("w") as handle:
                proc = subprocess.run(command, cwd=run_root, stdout=handle, stderr=subprocess.STDOUT, check=False)
            if proc.returncode or not output_json.exists():
                failures += 1
                print(f"ATOM FAIL {method} {element}")
                continue
            result = json.loads(output_json.read_text())
            data[method][element] = {
                "energy_hartree": float(result["energy"]),
                "unpaired_electrons": ATOMIC_SPINS[element],
                "spin_polarized": ATOMIC_SPINS[element] > 0,
            }
            print(f"ATOM OK {method} {element}: {float(result['energy']):.12f} Eh")
    atomic_energies_path().write_text(json.dumps(data, indent=2, sort_keys=True) + "\n")
    return 1 if failures else 0


def parse_factors(value: str) -> list[float]:
    values = [float(field) for field in value.split(",")]
    if len(values) < 1 or any(value <= 0 for value in values):
        raise argparse.ArgumentTypeError("volume factors must be positive")
    return values


def add_selection_options(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--datasets",
        nargs="+",
        choices=("klimes",),
        default=["klimes"],
    )
    parser.add_argument("--methods", nargs="+", choices=METHODS, default=list(METHODS))
    parser.add_argument("--systems", nargs="+", help="optional exact system identifiers")
    parser.add_argument("--volume-factors", type=parse_factors, default=list(DEFAULT_VOLUME_FACTORS))
    parser.add_argument("--cp2k", type=Path, default=CP2K_DEFAULT)
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--inversion-reduced", dest="full_grid", action="store_false")
    parser.set_defaults(full_grid=True)


def parser() -> argparse.ArgumentParser:
    top = argparse.ArgumentParser(description=__doc__)
    subparsers = top.add_subparsers(dest="command", required=True)
    run = subparsers.add_parser("run", help="run missing CP2K EOS points")
    add_selection_options(run)
    run.add_argument("--jobs", type=int, default=max(1, min(8, os.cpu_count() or 1)))
    run.add_argument("--start-attempt", type=int, choices=(1, 2, 3, 4, 5, 6, 7), default=1)
    run.add_argument("--force", action="store_true", help="recompute selected points even when sidecars exist")
    run.set_defaults(function=run_jobs)
    status = subparsers.add_parser("status", help="show completion counts")
    add_selection_options(status)
    status.set_defaults(function=status_jobs)
    analysis = subparsers.add_parser("analyse", help="fit EOS curves and compare with literature")
    add_selection_options(analysis)
    analysis.set_defaults(function=analyse)
    atoms = subparsers.add_parser("atoms", help="run isolated atoms for Klimes cohesive energies")
    atoms.add_argument("--methods", nargs="+", choices=METHODS, default=list(METHODS))
    atoms.add_argument("--tblite", type=Path, default=TBLITE_DEFAULT)
    atoms.set_defaults(function=run_atoms)
    return top


def main() -> int:
    args = parser().parse_args()
    try:
        return int(args.function(args))
    except (FileNotFoundError, ValueError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
