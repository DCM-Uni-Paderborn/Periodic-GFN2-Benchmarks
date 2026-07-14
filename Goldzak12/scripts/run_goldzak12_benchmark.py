#!/usr/bin/env python3
from __future__ import annotations

import argparse
import concurrent.futures
import csv
import hashlib
import json
import math
import os
import re
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CP2K = Path(os.environ.get("CP2K", "cp2k.ssmp"))
DEFAULT_TBLITE = Path(os.environ.get("TBLITE", "tblite"))
DEFAULT_CP2K_SOURCE = Path(os.environ.get("CP2K_SOURCE", "../cp2k"))
DEFAULT_TBLITE_SOURCE = Path(os.environ.get("TBLITE_SOURCE", "../tblite"))
HARTREE_TO_EV = 27.211386245988
FLOAT = r"[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[Ee][-+]?\d+)?"


@dataclass(frozen=True)
class Reference:
    solid: str
    structure: str
    formula: tuple[tuple[str, int], ...]
    a_exp: float
    a_hf: float
    a_mp2: float
    a_scs_mp2: float
    a_sos_mp2: float
    ecoh_exp: float
    ecoh_hf: float
    ecoh_mp2: float
    ecoh_scs_mp2: float
    ecoh_sos_mp2: float


REFERENCES: tuple[Reference, ...] = (
    Reference("C", "diamond", (("C", 1),), 3.553, 3.547, 3.540, 3.550, 3.554, 7.55, 5.38, 7.98, 7.65, 7.50),
    Reference("Si", "diamond", (("Si", 1),), 5.421, 5.508, 5.399, 5.425, 5.437, 4.70, 3.03, 4.97, 4.69, 4.56),
    Reference("SiC", "zincblende", (("Si", 1), ("C", 1)), 4.347, 4.371, 4.350, 4.358, 4.362, 6.47, 4.53, 6.79, 6.49, 6.35),
    Reference("BN", "zincblende", (("B", 1), ("N", 1)), 3.593, 3.596, 3.596, 3.603, 3.606, 6.76, 4.78, 7.13, 6.92, 6.82),
    Reference("BP", "zincblende", (("B", 1), ("P", 1)), 4.525, 4.584, 4.495, 4.517, 4.528, 5.14, 3.42, 5.58, 5.30, 5.16),
    Reference("AlN", "zincblende", (("Al", 1), ("N", 1)), 4.368, 4.365, 4.388, 4.389, 4.389, 5.85, 3.86, 6.00, 5.85, 5.78),
    Reference("AlP", "zincblende", (("Al", 1), ("P", 1)), 5.448, 5.542, 5.444, 5.465, 5.475, 4.31, 2.71, 4.42, 4.23, 4.14),
    Reference("MgS", "rocksalt", (("Mg", 1), ("S", 1)), 5.188, 5.281, 5.171, 5.191, 5.201, 4.04, 2.78, 4.20, 3.97, 3.86),
    Reference("LiF", "rocksalt", (("Li", 1), ("F", 1)), 3.973, 3.964, 3.990, 3.992, 3.993, 4.46, 3.41, 4.58, 4.49, 4.44),
    Reference("LiCl", "rocksalt", (("Li", 1), ("Cl", 1)), 5.072, 5.253, 5.021, 5.059, 5.078, 3.58, 2.73, 3.69, 3.58, 3.52),
)

METHODS = ("GFN1", "GFN2")
ELEMENT_MULTIPLICITY = {
    "Li": 2,
    "B": 2,
    "C": 3,
    "N": 4,
    "F": 2,
    "Mg": 1,
    "Al": 2,
    "Si": 3,
    "P": 4,
    "S": 3,
    "Cl": 2,
}


def fcc_sites() -> list[tuple[float, float, float]]:
    return [
        (0.0, 0.0, 0.0),
        (0.0, 0.5, 0.5),
        (0.5, 0.0, 0.5),
        (0.5, 0.5, 0.0),
    ]


def frac_add(a: tuple[float, float, float], b: tuple[float, float, float]) -> tuple[float, float, float]:
    return tuple((a[i] + b[i]) % 1.0 for i in range(3))  # type: ignore[return-value]


def conventional_cell_atoms(ref: Reference) -> list[tuple[str, float, float, float]]:
    sites = fcc_sites()
    if ref.structure == "diamond":
        element = ref.formula[0][0]
        return [(element, *p) for p in sites] + [(element, *frac_add(p, (0.25, 0.25, 0.25))) for p in sites]
    if ref.structure == "zincblende":
        a, b = ref.formula[0][0], ref.formula[1][0]
        return [(a, *p) for p in sites] + [(b, *frac_add(p, (0.25, 0.25, 0.25))) for p in sites]
    if ref.structure == "rocksalt":
        a, b = ref.formula[0][0], ref.formula[1][0]
        return [(a, *p) for p in sites] + [(b, *frac_add(p, (0.5, 0.0, 0.0))) for p in sites]
    raise ValueError(ref.structure)


def atom_counts(ref: Reference) -> dict[str, int]:
    counts: dict[str, int] = {}
    for element, *_ in conventional_cell_atoms(ref):
        counts[element] = counts.get(element, 0) + 1
    return counts


def kpoint_block(mesh: str) -> list[str]:
    if not mesh.startswith("k") or not mesh[1:].isdigit():
        raise ValueError(f"Bad mesh {mesh!r}; expected k333, k444, ...")
    digits = mesh[1:]
    if len(digits) == 3 and len(set(digits)) == 1:
        n = int(digits[0])
    else:
        n = int(digits)
    shift = 0.0 if n % 2 else 1.0 / (2.0 * n)
    return [
        "    &KPOINTS",
        f"      SCHEME MACDONALD {n} {n} {n} {shift:.10g} {shift:.10g} {shift:.10g}",
        "      SYMMETRY T",
        "      FULL_GRID F",
        "      SYMMETRY_BACKEND SPGLIB",
        "      SYMMETRY_REDUCTION_METHOD SPGLIB",
        "    &END KPOINTS",
    ]


def quickstep_block(
    method: str,
    periodic: bool,
    mesh: str | None = None,
    multiplicity: int | None = None,
    added_mos: int | None = None,
    smear_off: bool = False,
) -> list[str]:
    lines = [
        "  METHOD Quickstep",
        "  STRESS_TENSOR ANALYTICAL",
        "  &DFT",
    ]
    if multiplicity is not None:
        lines += [
            "    UKS T",
            f"    MULTIPLICITY {multiplicity}",
        ]
    lines += [
        "    &QS",
        "      EPS_DEFAULT 1.0E-12",
        "      METHOD xTB",
        "      &XTB",
        "        GFN_TYPE TBLITE",
        "        &TBLITE",
        f"          METHOD {method}",
        "          ACCURACY 0.05",
        "        &END TBLITE",
        "      &END XTB",
        "    &END QS",
    ]
    if not periodic:
        lines += [
            "    &POISSON",
            "      PERIODIC NONE",
            "    &END POISSON",
        ]
    if mesh:
        lines += kpoint_block(mesh)
    lines += [
        "    &SCF",
        "      EPS_SCF 1.0E-9",
        "      MAX_SCF 300",
        "      SCF_GUESS MOPAC",
    ]
    if added_mos is not None:
        lines.append(f"      ADDED_MOS {added_mos}")
    if smear_off:
        lines += [
            "      &SMEAR OFF",
            "      &END SMEAR",
        ]
    lines += [
        "      &MIXING",
        "        METHOD DIRECT_P_MIXING",
        "        ALPHA 0.2",
        "      &END MIXING",
        "      &PRINT",
        "        &RESTART OFF",
        "        &END RESTART",
        "      &END PRINT",
        "    &END SCF",
        "  &END DFT",
    ]
    return lines


def solid_input(ref: Reference, method: str, run_type: str, mesh: str, lattice_a: float, project: str) -> str:
    atoms = conventional_cell_atoms(ref)
    lines = [
        "&GLOBAL",
        "  PRINT_LEVEL LOW",
        f"  PROJECT {project}",
        f"  RUN_TYPE {run_type}",
        "&END GLOBAL",
        "",
        "&FORCE_EVAL",
    ]
    lines += quickstep_block(method, periodic=True, mesh=mesh)
    lines += [
        "  &SUBSYS",
        "    &CELL",
        "      CANONICALIZE TRUE",
        f"      ABC {lattice_a:.12f} {lattice_a:.12f} {lattice_a:.12f}",
        "      PERIODIC XYZ",
        "      SYMMETRY CUBIC",
        "    &END CELL",
        "    &COORD",
        "      SCALED",
    ]
    for element, x, y, z in atoms:
        lines.append(f"      {element:<2} {x: .12f} {y: .12f} {z: .12f}")
    lines += [
        "    &END COORD",
        "  &END SUBSYS",
        "&END FORCE_EVAL",
    ]
    if run_type == "CELL_OPT":
        lines += [
            "",
            "&MOTION",
            "  &CELL_OPT",
            "    OPTIMIZER BFGS",
            "    MAX_ITER 160",
            "    EXTERNAL_PRESSURE [bar] 0.0",
            "    KEEP_ANGLES T",
            "    KEEP_SYMMETRY T",
            "    MAX_DR 2.0E-3",
            "    RMS_DR 1.0E-3",
            "    MAX_FORCE 6.0E-4",
            "    RMS_FORCE 3.0E-4",
            "    PRESSURE_TOLERANCE [bar] 150.0",
            "    &BFGS",
            "      TRUST_RADIUS [angstrom] 0.05",
            "    &END BFGS",
            "  &END CELL_OPT",
            "&END MOTION",
        ]
    return "\n".join(lines) + "\n"


def atom_input(element: str, method: str) -> str:
    multiplicity = ELEMENT_MULTIPLICITY[element]
    project = f"atom_{element}_{method}"
    lines = [
        "&GLOBAL",
        "  PRINT_LEVEL LOW",
        f"  PROJECT {project}",
        "  RUN_TYPE ENERGY",
        "&END GLOBAL",
        "",
        "&FORCE_EVAL",
    ]
    lines += quickstep_block(method, periodic=False, multiplicity=multiplicity, smear_off=True)
    lines += [
        "  &SUBSYS",
        "    &CELL",
        "      ABC 30.0 30.0 30.0",
        "      PERIODIC NONE",
        "    &END CELL",
        "    &COORD",
        f"      {element:<2} 0.0 0.0 0.0",
        "    &END COORD",
        "  &END SUBSYS",
        "&END FORCE_EVAL",
    ]
    return "\n".join(lines) + "\n"


def write_file(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)


def project_name(solid: str, method: str, kind: str, mesh: str) -> str:
    clean_solid = solid.replace("/", "_")
    return f"{clean_solid}_{method}_{kind}_{mesh}"


def run_cp2k(cp2k: Path, inp: Path, out: Path, threads: int) -> int:
    out.parent.mkdir(parents=True, exist_ok=True)
    if out.exists():
        out.unlink()
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = str(threads)
    env["OMP_PROC_BIND"] = "false"
    main_log = inp.parent / "mainLog.out"
    if main_log.exists():
        main_log.unlink()
    proc = subprocess.run(
        [str(cp2k), "-i", str(inp.name), "-o", str(out.name)],
        cwd=inp.parent,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.STDOUT,
        env=env,
    )
    if main_log.exists() and (not out.exists() or "PROGRAM ENDED" not in out.read_text(errors="ignore")):
        shutil.copyfile(main_log, out)
    return proc.returncode


def output_ok(output: Path, require_opt: bool = False) -> bool:
    if not output.exists():
        return False
    text = output.read_text(errors="ignore")
    if "PROGRAM ENDED" not in text:
        return False
    bad = ("ABORT", "DID NOT CONVERGE", "SCF run NOT converged")
    if any(token in text for token in bad):
        return False
    if require_opt and "CELL OPTIMIZATION COMPLETED" not in text and "GEOMETRY OPTIMIZATION COMPLETED" not in text:
        return False
    return True


def parse_energy(output: Path) -> float | None:
    extrapolated_energy = None
    force_eval_energy = None
    if not output.exists():
        return None
    for line in output.read_text(errors="ignore").splitlines():
        if "Total energy (extrapolated to T->0)" in line:
            extrapolated_energy = float(line.split()[-1])
        if "ENERGY| Total FORCE_EVAL" in line:
            force_eval_energy = float(line.split()[-1])
    if extrapolated_energy is not None:
        return extrapolated_energy
    return force_eval_energy


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def command_output(command: list[str], allow_empty: bool = False) -> str:
    try:
        proc = subprocess.run(command, text=True, capture_output=True, timeout=60, check=False)
    except (OSError, subprocess.TimeoutExpired) as exc:
        return f"unavailable: {exc}"
    text = (proc.stdout + proc.stderr).strip()
    if allow_empty and not text and proc.returncode == 0:
        return ""
    return text if text else f"exit status {proc.returncode}"


def shared_library_hashes(executable: Path) -> dict[str, str]:
    roots = (executable.resolve().parent, executable.resolve().parent.parent / "lib")
    patterns = ("libcp2k*.dylib", "libcp2k*.so*", "libtblite*.dylib", "libtblite*.so*")
    libraries: dict[str, str] = {}
    seen: set[Path] = set()
    for root in roots:
        if not root.is_dir():
            continue
        for pattern in patterns:
            for candidate in sorted(root.glob(pattern)):
                resolved = candidate.resolve()
                if not resolved.is_file() or resolved in seen:
                    continue
                seen.add(resolved)
                libraries[resolved.name] = sha256(resolved)
    return libraries


def version_summary(executable: Path) -> str:
    lines = command_output([str(executable), "--version"]).splitlines()
    summary = []
    for line in lines:
        if line.strip().lower().startswith("compiler options:"):
            break
        summary.append(line.rstrip())
    return "\n".join(summary)


def git_metadata(source: Path) -> dict[str, object]:
    if not source.exists():
        return {"available": False}
    revision = command_output(["git", "-C", str(source), "rev-parse", "HEAD"])
    branch = command_output(["git", "-C", str(source), "branch", "--show-current"])
    status = command_output(["git", "-C", str(source), "status", "--short"], allow_empty=True)
    diff = command_output(["git", "-C", str(source), "diff", "--binary"], allow_empty=True)
    return {
        "available": True,
        "revision": revision,
        "branch": branch,
        "dirty": bool(status),
        "working_tree_diff_sha256": hashlib.sha256(diff.encode()).hexdigest(),
    }


def repository_patch_metadata() -> dict[str, dict[str, str]]:
    patches = {
        "cp2k": ROOT.parent / "patches" / "cp2k_trunk_tblite_full_symmetry_scc.patch",
        "tblite": ROOT.parent / "patches" / "tblite_main_pr350_wsc_derivatives.patch",
    }
    return {
        name: {"path": f"../../patches/{path.name}", "sha256": sha256(path)}
        for name, path in patches.items()
        if path.is_file()
    }


def write_build_provenance(
    cp2k: Path,
    tblite: Path,
    cp2k_source: Path,
    tblite_source: Path,
    protocol: dict[str, object],
) -> None:
    payload = {
        "cp2k": {
            "executable": cp2k.name,
            "sha256": sha256(cp2k),
            "shared_library_sha256": shared_library_hashes(cp2k),
            "version": version_summary(cp2k),
            "source": git_metadata(cp2k_source),
        },
        "tblite": {
            "executable": tblite.name,
            "sha256": sha256(tblite),
            "shared_library_sha256": shared_library_hashes(tblite),
            "version": version_summary(tblite),
            "source": git_metadata(tblite_source),
        },
        "repository_patches": repository_patch_metadata(),
        "protocol": protocol,
    }
    path = ROOT / "data" / "build_provenance.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")


def parse_cell_from_restart(path: Path) -> tuple[list[float], list[float], list[float]] | None:
    if not path.exists():
        return None
    lines = path.read_text(errors="ignore").splitlines()
    in_cell = False
    vectors: dict[str, list[float]] = {}
    abc: list[float] | None = None
    for line in lines:
        stripped = line.strip()
        upper = stripped.upper()
        if upper.split()[0] == "&CELL":
            in_cell = True
            continue
        if in_cell and upper in {"&END CELL", "&END"}:
            break
        if not in_cell:
            continue
        parts = stripped.split()
        if len(parts) >= 4 and parts[0].upper() in {"A", "B", "C"}:
            vectors[parts[0].upper()] = [float(parts[1]), float(parts[2]), float(parts[3])]
        elif len(parts) >= 4 and parts[0].upper() == "ABC":
            abc = [float(parts[1]), float(parts[2]), float(parts[3])]
    if {"A", "B", "C"} <= vectors.keys():
        return vectors["A"], vectors["B"], vectors["C"]
    if abc is not None:
        return [abc[0], 0.0, 0.0], [0.0, abc[1], 0.0], [0.0, 0.0, abc[2]]
    return None


def latest_restart(run_dir: Path, project: str) -> Path | None:
    matches = list(run_dir.glob(f"{project}*.restart"))
    if not matches:
        return None
    return max(matches, key=lambda p: p.stat().st_mtime)


def norm(vec: list[float]) -> float:
    return math.sqrt(sum(x * x for x in vec))


def volume(a: list[float], b: list[float], c: list[float]) -> float:
    return abs(
        a[0] * (b[1] * c[2] - b[2] * c[1])
        - a[1] * (b[0] * c[2] - b[2] * c[0])
        + a[2] * (b[0] * c[1] - b[1] * c[0])
    )


def optimized_lattice(run_dir: Path, project: str) -> float | None:
    restart = latest_restart(run_dir, project)
    if restart is None:
        return None
    cell = parse_cell_from_restart(restart)
    if cell is None:
        return None
    return sum(norm(v) for v in cell) / 3.0


def reference_rows() -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for ref in REFERENCES:
        rows.append(
            {
                "solid": ref.solid,
                "structure": ref.structure,
                "formula": "".join(f"{el}{n if n != 1 else ''}" for el, n in ref.formula),
                "a_exp_A": ref.a_exp,
                "a_HF_A": ref.a_hf,
                "a_MP2_A": ref.a_mp2,
                "a_SCS_MP2_A": ref.a_scs_mp2,
                "a_SOS_MP2_A": ref.a_sos_mp2,
                "ecoh_exp_eV_per_atom": ref.ecoh_exp,
                "ecoh_HF_eV_per_atom": ref.ecoh_hf,
                "ecoh_MP2_eV_per_atom": ref.ecoh_mp2,
                "ecoh_SCS_MP2_eV_per_atom": ref.ecoh_scs_mp2,
                "ecoh_SOS_MP2_eV_per_atom": ref.ecoh_sos_mp2,
            }
        )
    return rows


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("")
        return
    fieldnames = list(rows[0].keys())
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def setup_inputs(cell_mesh: str, energy_meshes: list[str]) -> None:
    write_csv(ROOT / "data" / "reference_goldzak2022.csv", reference_rows())
    for ref in REFERENCES:
        for method in METHODS:
            project = project_name(ref.solid, method, "cellopt", cell_mesh)
            text = solid_input(ref, method, "CELL_OPT", cell_mesh, ref.a_exp, project)
            write_file(ROOT / "inputs" / "cellopt" / method / ref.solid / f"{project}.inp", text)
            for mesh in energy_meshes:
                sp_project = project_name(ref.solid, method, "sp", mesh)
                text = solid_input(ref, method, "ENERGY", mesh, ref.a_exp, sp_project)
                write_file(ROOT / "inputs" / "single_point_initial" / method / ref.solid / f"{sp_project}.inp", text)
    elements = sorted({el for ref in REFERENCES for el, _ in ref.formula})
    for method in METHODS:
        for element in elements:
            write_file(ROOT / "inputs" / "atoms" / method / f"atom_{element}_{method}.inp", atom_input(element, method))


def run_jobs(job_specs: list[tuple[str, Path, Path, bool]], cp2k: Path, jobs: int, threads: int, force: bool) -> None:
    pending: list[tuple[str, Path, Path, bool]] = []
    for label, inp, out, require_opt in job_specs:
        if not force and output_ok(out, require_opt=require_opt):
            continue
        pending.append((label, inp, out, require_opt))
    if not pending:
        print("No jobs pending.")
        return

    def worker(spec: tuple[str, Path, Path, bool]) -> tuple[str, int, bool]:
        label, inp, out, require_opt = spec
        code = run_cp2k(cp2k, inp, out, threads)
        return label, code, output_ok(out, require_opt=require_opt)

    print(f"Running {len(pending)} jobs with {jobs} worker(s), OMP_NUM_THREADS={threads}.")
    with concurrent.futures.ThreadPoolExecutor(max_workers=jobs) as pool:
        futures = {pool.submit(worker, spec): spec for spec in pending}
        done = 0
        for future in concurrent.futures.as_completed(futures):
            label, code, ok = future.result()
            done += 1
            status = "ok" if ok else f"failed rc={code}"
            print(f"[{done:3d}/{len(pending):3d}] {status:14s} {label}", flush=True)


def atom_job_specs() -> list[tuple[str, Path, Path, bool]]:
    specs: list[tuple[str, Path, Path, bool]] = []
    elements = sorted({el for ref in REFERENCES for el, _ in ref.formula})
    for method in METHODS:
        for element in elements:
            inp = ROOT / "inputs" / "atoms" / method / f"atom_{element}_{method}.inp"
            out = ROOT / "runs" / "atoms" / method / element / f"atom_{element}_{method}.out"
            run_inp = out.parent / inp.name
            if not run_inp.exists() or run_inp.read_text() != inp.read_text():
                write_file(run_inp, inp.read_text())
            specs.append((f"atom {method} {element}", run_inp, out, False))
    return specs


def run_tblite_atom_jobs(tblite: Path, jobs: int, force: bool) -> None:
    elements = sorted({el for ref in REFERENCES for el, _ in ref.formula})
    specs: list[tuple[str, str, Path]] = []
    for method in METHODS:
        for element in elements:
            run_dir = ROOT / "runs" / "atoms_cli" / method / element
            json_path = run_dir / f"atom_{element}_{method}.json"
            out_path = run_dir / f"atom_{element}_{method}.out"
            xyz_path = run_dir / f"atom_{element}.xyz"
            write_file(xyz_path, f"1\n{element} atom\n{element} 0.0 0.0 0.0\n")
            if not force and json_path.exists() and parse_tblite_json_energy(json_path) is not None:
                continue
            specs.append((method, element, run_dir))
    if not specs:
        print("No tblite atom jobs pending.")
        return

    def worker(spec: tuple[str, str, Path]) -> tuple[str, int, bool]:
        method, element, run_dir = spec
        spin = ELEMENT_MULTIPLICITY[element] - 1
        xyz_path = run_dir / f"atom_{element}.xyz"
        json_path = run_dir / f"atom_{element}_{method}.json"
        out_path = run_dir / f"atom_{element}_{method}.out"
        cmd = [
            str(tblite),
            "run",
            "--method",
            method.lower(),
            "--spin",
            str(spin),
            "--acc",
            "0.05",
            "--json",
            json_path.name,
            "--no-restart",
            xyz_path.name,
        ]
        with out_path.open("w") as handle:
            proc = subprocess.run(cmd, cwd=run_dir, stdout=handle, stderr=subprocess.STDOUT)
        ok = proc.returncode == 0 and parse_tblite_json_energy(json_path) is not None
        return f"atom-cli {method} {element}", proc.returncode, ok

    print(f"Running {len(specs)} tblite atom jobs with {jobs} worker(s).")
    with concurrent.futures.ThreadPoolExecutor(max_workers=jobs) as pool:
        futures = {pool.submit(worker, spec): spec for spec in specs}
        done = 0
        for future in concurrent.futures.as_completed(futures):
            label, code, ok = future.result()
            done += 1
            status = "ok" if ok else f"failed rc={code}"
            print(f"[{done:3d}/{len(specs):3d}] {status:14s} {label}", flush=True)


def parse_tblite_json_energy(path: Path) -> float | None:
    if not path.exists():
        return None
    try:
        data = json.loads(path.read_text())
    except json.JSONDecodeError:
        return None
    energy = data.get("energy")
    return float(energy) if energy is not None else None


def cellopt_job_specs(cell_mesh: str) -> list[tuple[str, Path, Path, bool]]:
    specs: list[tuple[str, Path, Path, bool]] = []
    for ref in REFERENCES:
        for method in METHODS:
            project = project_name(ref.solid, method, "cellopt", cell_mesh)
            inp = ROOT / "inputs" / "cellopt" / method / ref.solid / f"{project}.inp"
            out = ROOT / "runs" / "cellopt" / method / ref.solid / cell_mesh / f"{project}.out"
            # CP2K writes restart files into the input directory, so run from an isolated copy.
            run_inp = out.parent / inp.name
            if not run_inp.exists() or run_inp.read_text() != inp.read_text():
                write_file(run_inp, inp.read_text())
            specs.append((f"cellopt {method} {ref.solid} {cell_mesh}", run_inp, out, True))
    return specs


def generate_final_sp_inputs(cell_mesh: str, energy_meshes: list[str]) -> None:
    missing: list[str] = []
    for ref in REFERENCES:
        for method in METHODS:
            cell_project = project_name(ref.solid, method, "cellopt", cell_mesh)
            run_dir = ROOT / "runs" / "cellopt" / method / ref.solid / cell_mesh
            opt_output = run_dir / f"{cell_project}.out"
            if not output_ok(opt_output, require_opt=True):
                missing.append(f"{method} {ref.solid}")
                continue
            a_opt = optimized_lattice(run_dir, cell_project)
            if a_opt is None:
                missing.append(f"{method} {ref.solid}")
                continue
            for mesh in energy_meshes:
                sp_project = project_name(ref.solid, method, "sp", mesh)
                text = solid_input(ref, method, "ENERGY", mesh, a_opt, sp_project)
                write_file(ROOT / "runs" / "single_point" / method / ref.solid / mesh / f"{sp_project}.inp", text)
    if missing:
        print("Missing optimized cells for:", ", ".join(missing), file=sys.stderr)


def sp_job_specs(energy_meshes: list[str]) -> list[tuple[str, Path, Path, bool]]:
    specs: list[tuple[str, Path, Path, bool]] = []
    for ref in REFERENCES:
        for method in METHODS:
            for mesh in energy_meshes:
                project = project_name(ref.solid, method, "sp", mesh)
                inp = ROOT / "runs" / "single_point" / method / ref.solid / mesh / f"{project}.inp"
                if not inp.exists():
                    continue
                out = ROOT / "runs" / "single_point" / method / ref.solid / mesh / f"{project}.out"
                specs.append((f"sp {method} {ref.solid} {mesh}", inp, out, False))
    return specs


def atom_energies() -> dict[tuple[str, str], float]:
    energies: dict[tuple[str, str], float] = {}
    for method in METHODS:
        atom_root = ROOT / "runs" / "atoms_cli" / method
        for element_dir in atom_root.glob("*"):
            if not element_dir.is_dir():
                continue
            out = element_dir / f"atom_{element_dir.name}_{method}.json"
            energy = parse_tblite_json_energy(out)
            if energy is not None:
                energies[(method, element_dir.name)] = energy
    write_csv(
        ROOT / "data" / "atom_energies_tblite_cli.csv",
        [
            {"method": method, "element": element, "energy_hartree": f"{energy:.12f}", "source": "tblite_cli"}
            for (method, element), energy in sorted(energies.items())
        ],
    )
    return energies


def analyse(cell_mesh: str, energy_meshes: list[str], result_mesh: str) -> None:
    atom_e = atom_energies()
    rows: list[dict[str, object]] = []
    for ref in REFERENCES:
        n_atoms = len(conventional_cell_atoms(ref))
        counts = atom_counts(ref)
        for method in METHODS:
            cell_project = project_name(ref.solid, method, "cellopt", cell_mesh)
            cell_run = ROOT / "runs" / "cellopt" / method / ref.solid / cell_mesh
            opt_out = cell_run / f"{cell_project}.out"
            cellopt_completed = output_ok(opt_out, require_opt=True)
            a_opt = optimized_lattice(cell_run, cell_project) if cellopt_completed else None
            opt_energy = parse_energy(opt_out)
            atom_sum = None
            if all((method, el) in atom_e for el in counts):
                atom_sum = sum(atom_e[(method, el)] * count for el, count in counts.items())
            for mesh in energy_meshes:
                sp_project = project_name(ref.solid, method, "sp", mesh)
                sp_out = ROOT / "runs" / "single_point" / method / ref.solid / mesh / f"{sp_project}.out"
                sp_energy = parse_energy(sp_out)
                ecoh = None
                if atom_sum is not None and sp_energy is not None:
                    ecoh = (atom_sum - sp_energy) * HARTREE_TO_EV / n_atoms
                rows.append(
                    {
                        "solid": ref.solid,
                        "structure": ref.structure,
                        "method": method,
                        "cell_mesh": cell_mesh,
                        "energy_mesh": mesh,
                        "cellopt_completed": cellopt_completed,
                        "sp_completed": output_ok(sp_out, require_opt=False),
                        "a_calc_A": f"{a_opt:.8f}" if a_opt is not None else "",
                        "a_ref_exp_A": ref.a_exp,
                        "a_error_A": f"{(a_opt - ref.a_exp):.8f}" if a_opt is not None else "",
                        "a_abs_error_A": f"{abs(a_opt - ref.a_exp):.8f}" if a_opt is not None else "",
                        "ecoh_calc_eV_per_atom": f"{ecoh:.8f}" if ecoh is not None else "",
                        "ecoh_ref_exp_eV_per_atom": ref.ecoh_exp,
                        "ecoh_error_eV_per_atom": f"{(ecoh - ref.ecoh_exp):.8f}" if ecoh is not None else "",
                        "ecoh_abs_error_eV_per_atom": f"{abs(ecoh - ref.ecoh_exp):.8f}" if ecoh is not None else "",
                        "solid_energy_hartree": f"{sp_energy:.12f}" if sp_energy is not None else "",
                        "cellopt_last_energy_hartree": f"{opt_energy:.12f}" if opt_energy is not None else "",
                        "n_atoms_conventional_cell": n_atoms,
                        "atom_reference_source": "tblite_cli",
                    }
                )
    write_csv(ROOT / "data" / "results.csv", rows)

    final_rows = [r for r in rows if r["energy_mesh"] == result_mesh]
    summary: list[dict[str, object]] = []
    for method in METHODS:
        method_rows = [r for r in final_rows if r["method"] == method and r["cellopt_completed"] and r["sp_completed"]]
        a_err = [float(r["a_error_A"]) for r in method_rows if r["a_error_A"] != ""]
        e_err = [float(r["ecoh_error_eV_per_atom"]) for r in method_rows if r["ecoh_error_eV_per_atom"] != ""]
        summary.append(
            {
                "method": method,
                "n_complete": len(method_rows),
                "result_mesh": result_mesh,
                "cell_mesh": cell_mesh,
                "a_ME_A": f"{sum(a_err) / len(a_err):.8f}" if a_err else "",
                "a_MAE_A": f"{sum(abs(x) for x in a_err) / len(a_err):.8f}" if a_err else "",
                "a_RMSE_A": f"{math.sqrt(sum(x * x for x in a_err) / len(a_err)):.8f}" if a_err else "",
                "ecoh_ME_eV_per_atom": f"{sum(e_err) / len(e_err):.8f}" if e_err else "",
                "ecoh_MAE_eV_per_atom": f"{sum(abs(x) for x in e_err) / len(e_err):.8f}" if e_err else "",
                "ecoh_RMSE_eV_per_atom": f"{math.sqrt(sum(x * x for x in e_err) / len(e_err)):.8f}" if e_err else "",
            }
        )
    write_csv(ROOT / "data" / "summary.csv", summary)
    write_markdown(final_rows, summary, result_mesh)
    plot(final_rows, result_mesh)


def write_markdown(rows: list[dict[str, object]], summary: list[dict[str, object]], result_mesh: str) -> None:
    order = [ref.solid for ref in REFERENCES]
    by_key = {(r["solid"], r["method"]): r for r in rows}
    lines = [
        f"# LC10 CP2K/tblite results ({result_mesh} final energies)",
        "",
        "All GFN values use native Bloch k-points in CP2K. Cohesive energies are in eV per atom.",
        "",
        "## Summary",
        "",
        "| method | n | a ME (A) | a MAE (A) | Ecoh ME (eV/atom) | Ecoh MAE (eV/atom) |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for row in summary:
        lines.append(
            f"| {row['method']} | {row['n_complete']} | {row['a_ME_A']} | {row['a_MAE_A']} | "
            f"{row['ecoh_ME_eV_per_atom']} | {row['ecoh_MAE_eV_per_atom']} |"
        )
    lines += [
        "",
        "## Per-system comparison to experiment",
        "",
        "| solid | a exp | a GFN1 | da GFN1 | a GFN2 | da GFN2 | Ecoh exp | Ecoh GFN1 | dE GFN1 | Ecoh GFN2 | dE GFN2 |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    refs = {ref.solid: ref for ref in REFERENCES}
    for solid in order:
        ref = refs[solid]
        g1 = by_key.get((solid, "GFN1"), {})
        g2 = by_key.get((solid, "GFN2"), {})
        lines.append(
            "| {solid} | {aexp:.3f} | {a1} | {da1} | {a2} | {da2} | {eexp:.2f} | {e1} | {de1} | {e2} | {de2} |".format(
                solid=solid,
                aexp=ref.a_exp,
                a1=fmt(g1.get("a_calc_A"), 4),
                da1=fmt(g1.get("a_error_A"), 4),
                a2=fmt(g2.get("a_calc_A"), 4),
                da2=fmt(g2.get("a_error_A"), 4),
                eexp=ref.ecoh_exp,
                e1=fmt(g1.get("ecoh_calc_eV_per_atom"), 3),
                de1=fmt(g1.get("ecoh_error_eV_per_atom"), 3),
                e2=fmt(g2.get("ecoh_calc_eV_per_atom"), 3),
                de2=fmt(g2.get("ecoh_error_eV_per_atom"), 3),
            )
        )
    (ROOT / "data").mkdir(parents=True, exist_ok=True)
    (ROOT / "data" / "results.md").write_text("\n".join(lines) + "\n")


def fmt(value: object, digits: int) -> str:
    if value is None or value == "":
        return ""
    try:
        return f"{float(value):.{digits}f}"
    except (TypeError, ValueError):
        return str(value)


def plot(rows: list[dict[str, object]], result_mesh: str) -> None:
    import matplotlib.pyplot as plt
    import numpy as np

    complete = [r for r in rows if r["cellopt_completed"] and r["sp_completed"]]
    if not complete:
        return
    solids = [ref.solid for ref in REFERENCES]
    x = np.arange(len(solids))
    width = 0.36
    colors = {"GFN1": "#4C78A8", "GFN2": "#F58518"}
    for prop, ylabel, filename, key in [
        ("lattice", "lattice-constant error (A)", "goldzak12_lattice_errors", "a_error_A"),
        ("cohesive", "cohesive-energy error (eV/atom)", "goldzak12_cohesive_errors", "ecoh_error_eV_per_atom"),
    ]:
        fig, ax = plt.subplots(figsize=(10.5, 4.6))
        for idx, method in enumerate(METHODS):
            values = []
            for solid in solids:
                row = next((r for r in complete if r["solid"] == solid and r["method"] == method), None)
                values.append(float(row[key]) if row and row[key] != "" else np.nan)
            offset = (idx - 0.5) * width
            ax.bar(x + offset, values, width, label=method, color=colors[method])
        ax.axhline(0.0, color="black", linewidth=0.8)
        ax.set_ylabel(ylabel)
        ax.set_xticks(x)
        ax.set_xticklabels(solids, rotation=45, ha="right")
        ax.set_title(f"LC10 CP2K/tblite native-Bloch {result_mesh}: {prop} errors")
        ax.legend(frameon=False)
        ax.grid(axis="y", color="#d0d0d0", linewidth=0.6, alpha=0.7)
        fig.tight_layout()
        out_base = ROOT / "figures" / filename
        out_base.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out_base.with_suffix(".png"), dpi=220)
        fig.savefig(out_base.with_suffix(".pdf"))
        plt.close(fig)

    ref_by_solid = {ref.solid: ref for ref in REFERENCES}
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.4))
    for method in METHODS:
        a_ref = []
        a_calc = []
        e_ref = []
        e_calc = []
        for solid in solids:
            row = next((r for r in complete if r["solid"] == solid and r["method"] == method), None)
            if not row:
                continue
            a_ref.append(ref_by_solid[solid].a_exp)
            a_calc.append(float(row["a_calc_A"]))
            e_ref.append(ref_by_solid[solid].ecoh_exp)
            e_calc.append(float(row["ecoh_calc_eV_per_atom"]))
        axes[0].scatter(a_ref, a_calc, label=method, color=colors[method], s=44)
        axes[1].scatter(e_ref, e_calc, label=method, color=colors[method], s=44)
    for ax, label in zip(axes, ["lattice constant (A)", "cohesive energy (eV/atom)"]):
        lo, hi = ax.get_xlim()
        ylo, yhi = ax.get_ylim()
        mn, mx = min(lo, ylo), max(hi, yhi)
        ax.plot([mn, mx], [mn, mx], color="black", linewidth=0.8)
        ax.set_xlim(mn, mx)
        ax.set_ylim(mn, mx)
        ax.set_xlabel(f"experiment {label}")
        ax.set_ylabel(f"CP2K/tblite {label}")
        ax.grid(color="#d0d0d0", linewidth=0.6, alpha=0.7)
    axes[0].legend(frameon=False)
    fig.suptitle(f"LC10 CP2K/tblite native-Bloch {result_mesh} vs experiment")
    fig.tight_layout()
    out_base = ROOT / "figures" / "goldzak12_scatter"
    fig.savefig(out_base.with_suffix(".png"), dpi=220)
    fig.savefig(out_base.with_suffix(".pdf"))
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers(dest="command", required=True)

    def common(p: argparse.ArgumentParser) -> None:
        p.add_argument("--cell-mesh", default="k444")
        p.add_argument("--energy-mesh", action="append", default=[])

    p_setup = sub.add_parser("setup")
    common(p_setup)

    p_run = sub.add_parser("run")
    common(p_run)
    p_run.add_argument("--cp2k", type=Path, default=DEFAULT_CP2K)
    p_run.add_argument("--tblite", type=Path, default=DEFAULT_TBLITE)
    p_run.add_argument("--jobs", type=int, default=4)
    p_run.add_argument("--threads", type=int, default=1)
    p_run.add_argument("--force", action="store_true")

    p_sp = sub.add_parser("single-points")
    common(p_sp)
    p_sp.add_argument("--cp2k", type=Path, default=DEFAULT_CP2K)
    p_sp.add_argument("--jobs", type=int, default=4)
    p_sp.add_argument("--threads", type=int, default=1)
    p_sp.add_argument("--force", action="store_true")

    p_analyse = sub.add_parser("analyse")
    common(p_analyse)
    p_analyse.add_argument("--result-mesh", default="")

    args = parser.parse_args()
    energy_meshes = args.energy_mesh or ["k333", "k444", "k555"]

    if args.command == "setup":
        setup_inputs(args.cell_mesh, energy_meshes)
        return 0

    if args.command == "run":
        setup_inputs(args.cell_mesh, energy_meshes)
        run_tblite_atom_jobs(args.tblite, args.jobs, args.force)
        run_jobs(cellopt_job_specs(args.cell_mesh), args.cp2k, args.jobs, args.threads, args.force)
        generate_final_sp_inputs(args.cell_mesh, energy_meshes)
        run_jobs(sp_job_specs(energy_meshes), args.cp2k, args.jobs, args.threads, args.force)
        analyse(args.cell_mesh, energy_meshes, energy_meshes[-1])
        return 0

    if args.command == "single-points":
        generate_final_sp_inputs(args.cell_mesh, energy_meshes)
        run_jobs(sp_job_specs(energy_meshes), args.cp2k, args.jobs, args.threads, args.force)
        analyse(args.cell_mesh, energy_meshes, energy_meshes[-1])
        return 0

    if args.command == "analyse":
        result_mesh = args.result_mesh or energy_meshes[-1]
        analyse(args.cell_mesh, energy_meshes, result_mesh)
        return 0

    return 2


if __name__ == "__main__":
    raise SystemExit(main())
