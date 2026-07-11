#!/usr/bin/env python3
from __future__ import annotations

import argparse
import concurrent.futures
import csv
import json
import math
import os
import subprocess
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import run_goldzak12_benchmark as base  # noqa: E402


ROOT = base.ROOT
DEFAULT_SCALES = (0.82, 0.88, 0.94, 0.98, 1.00, 1.02, 1.06, 1.12, 1.20, 1.30, 1.45)
ADAPTIVE_SCALES = {("MgO", "GFN2"): (0.90, 0.92)}


def scale_tag(scale: float) -> str:
    return f"s{scale:.3f}".replace(".", "p")


def eos_project(solid: str, method: str, mesh: str, scale: float) -> str:
    return f"{solid}_{method}_eos_{mesh}_{scale_tag(scale)}"


def final_project(solid: str, method: str, mesh: str) -> str:
    return f"{solid}_{method}_eos_final_{mesh}"


def scales_for(solid: str, method: str, scales: tuple[float, ...]) -> tuple[float, ...]:
    return tuple(sorted(set(scales) | set(ADAPTIVE_SCALES.get((solid, method), ()))))


def eos_job_specs(mesh: str, scales: tuple[float, ...]) -> list[tuple[str, Path, Path, bool]]:
    specs: list[tuple[str, Path, Path, bool]] = []
    for ref in base.REFERENCES:
        for method in base.METHODS:
            for scale in scales_for(ref.solid, method, scales):
                project = eos_project(ref.solid, method, mesh, scale)
                a = ref.a_exp * scale
                run_dir = ROOT / "runs" / "eos" / method / ref.solid / mesh / scale_tag(scale)
                inp = run_dir / f"{project}.inp"
                out = run_dir / f"{project}.out"
                text = base.solid_input(ref, method, "ENERGY", mesh, a, project)
                base.write_file(inp, text)
                specs.append((f"eos {method} {ref.solid} {mesh} {scale_tag(scale)}", inp, out, False))
    return specs


def strategy_path(output: Path) -> Path:
    return output.with_suffix(".strategy.json")


def write_strategy(output: Path, strategy: str, completed: bool) -> None:
    strategy_path(output).write_text(json.dumps({"strategy": strategy, "completed": completed}, indent=2) + "\n")


def read_strategy(output: Path) -> str:
    path = strategy_path(output)
    if not path.exists():
        return "unknown"


def retries_exhausted(output: Path) -> bool:
    path = strategy_path(output)
    if not path.exists():
        return False
    try:
        data = json.loads(path.read_text())
    except json.JSONDecodeError:
        return False
    return data.get("strategy") == "retry_m1_d001" and not bool(data.get("completed"))
    try:
        return str(json.loads(path.read_text())["strategy"])
    except (json.JSONDecodeError, KeyError):
        return "unknown"


def retry_input(inp: Path, iterations: int, memory: int, damping: float, label: str) -> Path:
    text = inp.read_text()
    marker = "        &END TBLITE\n      &END XTB"
    replacement = (
        "        &END TBLITE\n"
        "        SCC_MIXER TBLITE\n"
        "        &TBLITE_MIXER\n"
        f"          ITERATIONS {iterations}\n"
        f"          MEMORY {memory}\n"
        f"          DAMPING {damping:.6f}\n"
        "        &END TBLITE_MIXER\n"
        "      &END XTB"
    )
    if marker not in text:
        raise ValueError(f"Cannot insert TBLITE_MIXER into {inp}")
    text = text.replace(marker, replacement, 1)
    text = text.replace("      MAX_SCF 300", f"      MAX_SCF {iterations}", 1)
    path = inp.with_name(f"{inp.stem}_{label}.inp")
    base.write_file(path, text)
    return path


def run_jobs(
    specs: list[tuple[str, Path, Path, bool]],
    cp2k: Path,
    jobs: int,
    threads: int,
    force: bool,
    retry_scf: bool = True,
) -> None:
    pending = [
        spec
        for spec in specs
        if force
        or (
            not base.output_ok(spec[2], require_opt=spec[3])
            and not retries_exhausted(spec[2])
        )
    ]
    if not pending:
        print("No EOS jobs pending.")
        return

    def worker(spec: tuple[str, Path, Path, bool]) -> tuple[str, int, bool, tuple[str, Path, Path, bool]]:
        label, inp, out, require_opt = spec
        code = base.run_cp2k(cp2k, inp, out, threads)
        ok = base.output_ok(out, require_opt=require_opt)
        write_strategy(out, "default_tblite_mixer", ok)
        return label, code, ok, spec

    print(f"Running {len(pending)} CP2K jobs with {jobs} worker(s), OMP_NUM_THREADS={threads}.")
    with concurrent.futures.ThreadPoolExecutor(max_workers=jobs) as pool:
        futures = {pool.submit(worker, spec): spec for spec in pending}
        done = 0
        failed: list[tuple[str, Path, Path, bool]] = []
        for future in concurrent.futures.as_completed(futures):
            label, code, ok, spec = future.result()
            done += 1
            status = "ok" if ok else f"failed rc={code}"
            print(f"[{done:3d}/{len(pending):3d}] {status:14s} {label}", flush=True)
            if not ok:
                failed.append(spec)

    if not retry_scf or not failed:
        return

    profiles = (
        ("retry_m1_d005", 1200, 1, 0.05),
        ("retry_m1_d001", 2400, 1, 0.01),
    )
    for profile, iterations, memory, damping in profiles:
        if not failed:
            break
        print(
            f"Retrying {len(failed)} failed job(s) with TBLITE_MIXER "
            f"MEMORY={memory}, DAMPING={damping}, ITERATIONS={iterations}."
        )

        def retry_worker(spec: tuple[str, Path, Path, bool]) -> tuple[str, int, bool, tuple[str, Path, Path, bool]]:
            label, inp, out, require_opt = spec
            robust_inp = retry_input(inp, iterations, memory, damping, profile)
            code = base.run_cp2k(cp2k, robust_inp, out, threads)
            ok = base.output_ok(out, require_opt=require_opt)
            write_strategy(out, profile, ok)
            return label, code, ok, spec

        next_failed: list[tuple[str, Path, Path, bool]] = []
        with concurrent.futures.ThreadPoolExecutor(max_workers=min(jobs, 4)) as pool:
            futures = {pool.submit(retry_worker, spec): spec for spec in failed}
            done = 0
            for future in concurrent.futures.as_completed(futures):
                label, code, ok, spec = future.result()
                done += 1
                status = "ok" if ok else f"failed rc={code}"
                print(f"[retry {done:2d}/{len(failed):2d}] {status:14s} {label}", flush=True)
                if not ok:
                    next_failed.append(spec)
        failed = next_failed

    if failed:
        print("SCF retries exhausted for: " + ", ".join(spec[0] for spec in failed), file=sys.stderr)


def load_eos_points(
    ref: base.Reference,
    method: str,
    mesh: str,
    scales: tuple[float, ...],
) -> list[tuple[float, float, float | None, bool]]:
    points: list[tuple[float, float, float | None, bool]] = []
    for scale in scales:
        project = eos_project(ref.solid, method, mesh, scale)
        out = ROOT / "runs" / "eos" / method / ref.solid / mesh / scale_tag(scale) / f"{project}.out"
        energy = base.parse_energy(out)
        ok = base.output_ok(out)
        points.append((ref.a_exp * scale, scale, energy, ok))
    return sorted(points)


def fit_eos(points: list[tuple[float, float, float, bool]]) -> dict[str, object]:
    if len(points) < 3:
        return {"a_eos_A": "", "energy_fit_hartree": "", "fit_status": "insufficient_points", "n_points": len(points)}
    energies = np.array([p[2] for p in points], dtype=float)
    local_minima = [
        i
        for i in range(1, len(points) - 1)
        if points[i][2] < points[i - 1][2] and points[i][2] < points[i + 1][2]
    ]
    if not local_minima:
        return {
            "a_eos_A": "",
            "energy_fit_hartree": "",
            "fit_status": "no_local_minimum",
            "n_points": len(points),
            "grid_min_a_A": f"{points[int(np.argmin(energies))][0]:.10f}",
            "grid_min_scale": f"{points[int(np.argmin(energies))][1]:.5f}",
            "grid_min_energy_hartree": f"{points[int(np.argmin(energies))][2]:.12f}",
        }
    preferred = [i for i in local_minima if 0.88 <= points[i][1] <= 1.12]
    if preferred:
        idx = min(preferred, key=lambda i: points[i][2])
    else:
        idx = min(local_minima, key=lambda i: abs(points[i][1] - 1.0))
    lo = max(0, idx - 2)
    hi = min(len(points), idx + 3)
    if hi - lo < 3:
        lo = max(0, min(lo, len(points) - 3))
        hi = min(len(points), lo + 3)
    fit_points = points[lo:hi]
    x = np.array([p[0] for p in fit_points], dtype=float)
    y = np.array([p[2] for p in fit_points], dtype=float)
    coeff = np.polyfit(x, y, 2)
    fit_rmse = float(np.sqrt(np.mean((np.polyval(coeff, x) - y) ** 2)))
    status = "quadratic"
    if coeff[0] <= 0:
        a_min = points[idx][0]
        e_min = points[idx][2]
        status = "grid_min_negative_curvature"
    else:
        a_min = float(-coeff[1] / (2.0 * coeff[0]))
        e_min = float(np.polyval(coeff, a_min))
        if a_min < points[0][0] or a_min > points[-1][0]:
            a_min = points[idx][0]
            e_min = points[idx][2]
            status = "grid_min_outside_fit"
        elif fit_rmse > 2.0e-2 or e_min > points[idx][2] + 2.0e-2:
            return {
                "a_eos_A": "",
                "energy_fit_hartree": "",
                "fit_status": "poor_quadratic_fit",
                "fit_rmse_hartree": f"{fit_rmse:.12f}",
                "n_points": len(points),
                "grid_min_a_A": f"{points[idx][0]:.10f}",
                "grid_min_scale": f"{points[idx][1]:.5f}",
                "grid_min_energy_hartree": f"{points[idx][2]:.12f}",
            }
    return {
        "a_eos_A": f"{a_min:.10f}",
        "energy_fit_hartree": f"{e_min:.12f}",
        "fit_status": status,
        "fit_rmse_hartree": f"{fit_rmse:.12f}",
        "n_points": len(points),
        "grid_min_a_A": f"{points[idx][0]:.10f}",
        "grid_min_scale": f"{points[idx][1]:.5f}",
        "grid_min_energy_hartree": f"{points[idx][2]:.12f}",
    }


def make_eos_table(mesh: str, scales: tuple[float, ...]) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    point_rows: list[dict[str, object]] = []
    for ref in base.REFERENCES:
        for method in base.METHODS:
            requested_scales = scales_for(ref.solid, method, scales)
            all_points = load_eos_points(ref, method, mesh, requested_scales)
            points = [(a, scale, energy, ok) for a, scale, energy, ok in all_points if energy is not None and ok]
            fit = fit_eos(points)
            rows.append(
                {
                    "solid": ref.solid,
                    "structure": ref.structure,
                    "method": method,
                    "eos_mesh": mesh,
                    "a_exp_A": ref.a_exp,
                    "n_requested": len(requested_scales),
                    "n_completed": len(points),
                    **fit,
                }
            )
            for a, scale, energy, ok in all_points:
                project = eos_project(ref.solid, method, mesh, scale)
                output = ROOT / "runs" / "eos" / method / ref.solid / mesh / scale_tag(scale) / f"{project}.out"
                point_rows.append(
                    {
                        "solid": ref.solid,
                        "method": method,
                        "mesh": mesh,
                        "scale": f"{scale:.5f}",
                        "a_A": f"{a:.10f}",
                        "energy_hartree": f"{energy:.12f}" if energy is not None else "",
                        "completed": ok,
                        "scf_strategy": read_strategy(output),
                    }
                )
    base.write_csv(ROOT / "data" / "eos_points.csv", point_rows)
    base.write_csv(ROOT / "data" / "eos_fits.csv", rows)
    return rows


def final_sp_specs(fits: list[dict[str, object]], meshes: list[str]) -> list[tuple[str, Path, Path, bool]]:
    refs = {ref.solid: ref for ref in base.REFERENCES}
    specs: list[tuple[str, Path, Path, bool]] = []
    for row in fits:
        a_text = row.get("a_eos_A", "")
        if a_text == "":
            continue
        ref = refs[str(row["solid"])]
        method = str(row["method"])
        a = float(a_text)
        for mesh in meshes:
            project = final_project(ref.solid, method, mesh)
            run_dir = ROOT / "runs" / "eos_final_sp" / method / ref.solid / mesh
            inp = run_dir / f"{project}.inp"
            out = run_dir / f"{project}.out"
            base.write_file(inp, base.solid_input(ref, method, "ENERGY", mesh, a, project))
            specs.append((f"eos-final {method} {ref.solid} {mesh}", inp, out, False))
    return specs


def collect_results(fits: list[dict[str, object]], result_meshes: list[str], result_mesh: str) -> tuple[list[dict[str, object]], list[dict[str, object]], list[dict[str, object]]]:
    refs = {ref.solid: ref for ref in base.REFERENCES}
    atom_e = base.atom_energies()
    rows: list[dict[str, object]] = []
    for fit in fits:
        if not fit.get("a_eos_A"):
            continue
        ref = refs[str(fit["solid"])]
        method = str(fit["method"])
        a_calc = float(fit["a_eos_A"])
        n_atoms = len(base.conventional_cell_atoms(ref))
        counts = base.atom_counts(ref)
        atom_sum = None
        if all((method, el) in atom_e for el in counts):
            atom_sum = sum(atom_e[(method, el)] * count for el, count in counts.items())
        for mesh in result_meshes:
            project = final_project(ref.solid, method, mesh)
            out = ROOT / "runs" / "eos_final_sp" / method / ref.solid / mesh / f"{project}.out"
            e_solid = base.parse_energy(out)
            ecoh = (atom_sum - e_solid) * base.HARTREE_TO_EV / n_atoms if atom_sum is not None and e_solid is not None else None
            rows.append(
                {
                    "solid": ref.solid,
                    "structure": ref.structure,
                    "method": method,
                    "eos_mesh": fit["eos_mesh"],
                    "energy_mesh": mesh,
                    "fit_status": fit["fit_status"],
                    "sp_completed": base.output_ok(out),
                    "sp_scf_strategy": read_strategy(out),
                    "a_calc_A": f"{a_calc:.8f}" if a_calc is not None else "",
                    "a_ref_exp_A": ref.a_exp,
                    "a_error_A": f"{(a_calc - ref.a_exp):.8f}" if a_calc is not None else "",
                    "a_abs_error_A": f"{abs(a_calc - ref.a_exp):.8f}" if a_calc is not None else "",
                    "ecoh_calc_eV_per_atom": f"{ecoh:.8f}" if ecoh is not None else "",
                    "ecoh_ref_exp_eV_per_atom": ref.ecoh_exp,
                    "ecoh_error_eV_per_atom": f"{(ecoh - ref.ecoh_exp):.8f}" if ecoh is not None else "",
                    "ecoh_abs_error_eV_per_atom": f"{abs(ecoh - ref.ecoh_exp):.8f}" if ecoh is not None else "",
                    "solid_energy_hartree": f"{e_solid:.12f}" if e_solid is not None else "",
                    "atom_reference_source": "tblite_cli",
                    "a_HF_A": ref.a_hf,
                    "a_MP2_A": ref.a_mp2,
                    "a_SCS_MP2_A": ref.a_scs_mp2,
                    "a_SOS_MP2_A": ref.a_sos_mp2,
                    "ecoh_HF_eV_per_atom": ref.ecoh_hf,
                    "ecoh_MP2_eV_per_atom": ref.ecoh_mp2,
                    "ecoh_SCS_MP2_eV_per_atom": ref.ecoh_scs_mp2,
                    "ecoh_SOS_MP2_eV_per_atom": ref.ecoh_sos_mp2,
                }
            )
    base.write_csv(ROOT / "data" / "eos_results.csv", rows)
    summary = summary_rows(rows, result_mesh)
    lit_summary = literature_summary()
    base.write_csv(ROOT / "data" / "eos_summary.csv", summary + lit_summary)
    convergence = kpoint_convergence(rows)
    base.write_csv(ROOT / "data" / "eos_kpoint_convergence.csv", convergence)
    write_markdown(rows, summary, lit_summary, convergence, result_mesh, fits)
    plot(rows, summary, lit_summary, result_mesh)
    plot_eos_diagnostics(fits)
    return rows, summary, convergence


def summary_rows(rows: list[dict[str, object]], result_mesh: str) -> list[dict[str, object]]:
    summary: list[dict[str, object]] = []
    for method in base.METHODS:
        selected = [r for r in rows if r["method"] == method and r["energy_mesh"] == result_mesh and r["sp_completed"]]
        a_err = [float(r["a_error_A"]) for r in selected if r["a_error_A"] != ""]
        e_err = [float(r["ecoh_error_eV_per_atom"]) for r in selected if r["ecoh_error_eV_per_atom"] != ""]
        summary.append(
            {
                "source": "CP2K/tblite EOS",
                "method": method,
                "n_complete": len(selected),
                "a_ME_A": mean(a_err),
                "a_MAE_A": mae(a_err),
                "a_RMSE_A": rmse(a_err),
                "ecoh_ME_eV_per_atom": mean(e_err),
                "ecoh_MAE_eV_per_atom": mae(e_err),
                "ecoh_RMSE_eV_per_atom": rmse(e_err),
            }
        )
    return summary


def literature_summary() -> list[dict[str, object]]:
    mapping = {
        "HF": ("a_hf", "ecoh_hf"),
        "MP2": ("a_mp2", "ecoh_mp2"),
        "SCS-MP2": ("a_scs_mp2", "ecoh_scs_mp2"),
        "SOS-MP2": ("a_sos_mp2", "ecoh_sos_mp2"),
    }
    rows: list[dict[str, object]] = []
    for name, (akey, ekey) in mapping.items():
        a_err = [getattr(ref, akey) - ref.a_exp for ref in base.REFERENCES]
        e_err = [getattr(ref, ekey) - ref.ecoh_exp for ref in base.REFERENCES]
        rows.append(
            {
                "source": "Goldzak2022",
                "method": name,
                "n_complete": len(base.REFERENCES),
                "a_ME_A": mean(a_err),
                "a_MAE_A": mae(a_err),
                "a_RMSE_A": rmse(a_err),
                "ecoh_ME_eV_per_atom": mean(e_err),
                "ecoh_MAE_eV_per_atom": mae(e_err),
                "ecoh_RMSE_eV_per_atom": rmse(e_err),
            }
        )
    return rows


def kpoint_convergence(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    by_key = {(r["solid"], r["method"], r["energy_mesh"]): r for r in rows}
    conv: list[dict[str, object]] = []
    for ref in base.REFERENCES:
        for method in base.METHODS:
            dense = by_key.get((ref.solid, method, "k555"))
            if not dense or dense["ecoh_calc_eV_per_atom"] == "":
                continue
            e_dense = float(dense["ecoh_calc_eV_per_atom"])
            for mesh in ("k333", "k444"):
                row = by_key.get((ref.solid, method, mesh))
                if not row or row["ecoh_calc_eV_per_atom"] == "":
                    continue
                e = float(row["ecoh_calc_eV_per_atom"])
                conv.append(
                    {
                        "solid": ref.solid,
                        "method": method,
                        "mesh": mesh,
                        "reference_mesh": "k555",
                        "delta_ecoh_eV_per_atom": f"{(e - e_dense):.8f}",
                    }
                )
    return conv


def mean(values: list[float]) -> str:
    return f"{(sum(values) / len(values)):.8f}" if values else ""


def mae(values: list[float]) -> str:
    return f"{(sum(abs(v) for v in values) / len(values)):.8f}" if values else ""


def rmse(values: list[float]) -> str:
    return f"{math.sqrt(sum(v * v for v in values) / len(values)):.8f}" if values else ""


def write_markdown(
    rows: list[dict[str, object]],
    summary: list[dict[str, object]],
    lit_summary: list[dict[str, object]],
    convergence: list[dict[str, object]],
    result_mesh: str,
    fits: list[dict[str, object]],
) -> None:
    selected = [r for r in rows if r["energy_mesh"] == result_mesh]
    by_key = {(r["solid"], r["method"]): r for r in selected}
    refs = {ref.solid: ref for ref in base.REFERENCES}
    lines = [
        f"# LC12 (Goldzak12) EOS results ({result_mesh} final energies)",
        "",
        "Solid energies use CP2K/tblite native Bloch k-points. Atomic references use the matching tblite CLI.",
        "",
        "## EOS fit coverage",
        "",
        "| method | valid fits | excluded EOS curves |",
        "|---|---:|---|",
    ]
    fit_labels = {
        "no_local_minimum": "no bracketed minimum",
        "poor_quadratic_fit": "discontinuous EOS",
    }
    for method in base.METHODS:
        method_fits = [fit for fit in fits if fit["method"] == method]
        excluded = [
            f"{fit['solid']} ({fit_labels.get(str(fit['fit_status']), str(fit['fit_status']))})"
            for fit in method_fits
            if fit.get("a_eos_A", "") == ""
        ]
        lines.append(f"| {method} | {len(method_fits) - len(excluded)}/{len(method_fits)} | {', '.join(excluded) or '-'} |")
    lines += [
        "",
        "## MAE comparison to experiment",
        "",
        "| source | method | n | a ME (A) | a MAE (A) | Ecoh ME (eV/atom) | Ecoh MAE (eV/atom) |",
        "|---|---|---:|---:|---:|---:|---:|",
    ]
    for row in summary + lit_summary:
        lines.append(
            f"| {row['source']} | {row['method']} | {row['n_complete']} | {row['a_ME_A']} | {row['a_MAE_A']} | "
            f"{row['ecoh_ME_eV_per_atom']} | {row['ecoh_MAE_eV_per_atom']} |"
        )
    lines += [
        "",
        "## Per-system GFN comparison",
        "",
        "| solid | a exp | a GFN1 | da GFN1 | a GFN2 | da GFN2 | Ecoh exp | Ecoh GFN1 | dE GFN1 | Ecoh GFN2 | dE GFN2 |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for ref in base.REFERENCES:
        g1 = by_key.get((ref.solid, "GFN1"), {})
        g2 = by_key.get((ref.solid, "GFN2"), {})
        lines.append(
            "| {solid} | {aexp:.3f} | {a1} | {da1} | {a2} | {da2} | {eexp:.2f} | {e1} | {de1} | {e2} | {de2} |".format(
                solid=ref.solid,
                aexp=ref.a_exp,
                a1=base.fmt(g1.get("a_calc_A"), 4),
                da1=base.fmt(g1.get("a_error_A"), 4),
                a2=base.fmt(g2.get("a_calc_A"), 4),
                da2=base.fmt(g2.get("a_error_A"), 4),
                eexp=ref.ecoh_exp,
                e1=base.fmt(g1.get("ecoh_calc_eV_per_atom"), 3),
                de1=base.fmt(g1.get("ecoh_error_eV_per_atom"), 3),
                e2=base.fmt(g2.get("ecoh_calc_eV_per_atom"), 3),
                de2=base.fmt(g2.get("ecoh_error_eV_per_atom"), 3),
            )
        )
    lines += [
        "",
        "## k-point convergence of cohesive energies",
        "",
        "| method | mesh vs k555 | mean abs delta (eV/atom) | max abs delta (eV/atom) |",
        "|---|---|---:|---:|",
    ]
    for method in base.METHODS:
        for mesh in ("k333", "k444"):
            vals = [abs(float(r["delta_ecoh_eV_per_atom"])) for r in convergence if r["method"] == method and r["mesh"] == mesh]
            lines.append(f"| {method} | {mesh} | {sum(vals) / len(vals):.6f} | {max(vals):.6f} |" if vals else f"| {method} | {mesh} |  |  |")
    (ROOT / "data" / "eos_results.md").write_text("\n".join(lines) + "\n")


def plot(rows: list[dict[str, object]], summary: list[dict[str, object]], lit_summary: list[dict[str, object]], result_mesh: str) -> None:
    selected = [r for r in rows if r["energy_mesh"] == result_mesh and r["sp_completed"]]
    solids = [ref.solid for ref in base.REFERENCES]
    x = np.arange(len(solids))
    width = 0.36
    colors = {"GFN1": "#4C78A8", "GFN2": "#F58518"}
    for key, ylabel, name in [
        ("a_error_A", "lattice-constant error (A)", "goldzak12_eos_lattice_errors"),
        ("ecoh_error_eV_per_atom", "cohesive-energy error (eV/atom)", "goldzak12_eos_cohesive_errors"),
    ]:
        fig, ax = plt.subplots(figsize=(10.5, 4.6))
        for i, method in enumerate(base.METHODS):
            vals = []
            missing_positions = []
            for solid in solids:
                row = next((r for r in selected if r["solid"] == solid and r["method"] == method), None)
                vals.append(float(row[key]) if row and row[key] != "" else np.nan)
            positions = x + (i - 0.5) * width
            missing_positions = [position for position, value in zip(positions, vals) if np.isnan(value)]
            ax.bar(positions, vals, width, label=method, color=colors[method])
            for position in missing_positions:
                ax.annotate(
                    "n/a",
                    (position, 0.0),
                    xytext=(0, 5),
                    textcoords="offset points",
                    ha="center",
                    va="bottom",
                    color="#666666",
                    fontsize=8,
                )
        ax.axhline(0, color="black", linewidth=0.8)
        ax.set_xticks(x)
        ax.set_xticklabels(solids, rotation=45, ha="right")
        ax.set_ylabel(ylabel)
        ax.set_title(f"LC12 (Goldzak12) EOS CP2K/tblite native-Bloch {result_mesh}")
        ax.legend(frameon=False)
        ax.grid(axis="y", color="#d0d0d0", linewidth=0.6, alpha=0.7)
        fig.tight_layout()
        out = ROOT / "figures" / name
        out.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out.with_suffix(".png"), dpi=220)
        fig.savefig(out.with_suffix(".pdf"))
        plt.close(fig)

    labels = [f"{r['method']}\n(n={r['n_complete']})" for r in lit_summary + summary]
    a_mae = [float(r["a_MAE_A"]) for r in lit_summary] + [float(r["a_MAE_A"]) for r in summary]
    e_mae = [float(r["ecoh_MAE_eV_per_atom"]) for r in lit_summary] + [float(r["ecoh_MAE_eV_per_atom"]) for r in summary]
    fig, axes = plt.subplots(1, 2, figsize=(10.8, 4.2))
    axes[0].bar(labels, a_mae, color=["#72B7B2"] * len(lit_summary) + ["#4C78A8", "#F58518"])
    axes[1].bar(labels, e_mae, color=["#72B7B2"] * len(lit_summary) + ["#4C78A8", "#F58518"])
    axes[0].set_ylabel("MAE a (A)")
    axes[1].set_ylabel("MAE Ecoh (eV/atom)")
    for ax in axes:
        ax.tick_params(axis="x", rotation=35)
        ax.grid(axis="y", color="#d0d0d0", linewidth=0.6, alpha=0.7)
    fig.suptitle("LC12 (Goldzak12) comparison to zero-point corrected experiment")
    fig.tight_layout()
    out = ROOT / "figures" / "goldzak12_eos_mae_comparison"
    fig.savefig(out.with_suffix(".png"), dpi=220)
    fig.savefig(out.with_suffix(".pdf"))
    plt.close(fig)


def plot_eos_diagnostics(fits: list[dict[str, object]]) -> None:
    with (ROOT / "data" / "eos_points.csv").open(newline="") as handle:
        points = list(csv.DictReader(handle))
    fit_by_key = {(str(row["solid"]), str(row["method"])): row for row in fits}
    systems = ("MgO", "LiH")
    fig, axes = plt.subplots(1, 2, figsize=(10.8, 4.4), sharey=True)
    for ax, solid in zip(axes, systems):
        selected = sorted(
            (row for row in points if row["solid"] == solid and row["method"] == "GFN2"),
            key=lambda row: float(row["scale"]),
        )
        completed = [row for row in selected if row["completed"] == "True" and row["energy_hartree"] != ""]
        failed = [row for row in selected if row["completed"] != "True"]
        energy_min = min(float(row["energy_hartree"]) for row in completed)
        scales = [float(row["scale"]) for row in completed]
        relative = [
            (float(row["energy_hartree"]) - energy_min) * base.HARTREE_TO_EV / 8.0 for row in completed
        ]
        ax.plot(scales, relative, color="#D97706", linewidth=1.2, alpha=0.75)
        ax.scatter(scales, relative, color="#D97706", s=38, label="converged")
        marker_height = max(relative) * 1.08 if max(relative) > 0 else 1.0
        for row in failed:
            scale = float(row["scale"])
            ax.axvline(scale, color="#888888", linewidth=0.8, linestyle="--", alpha=0.6)
            ax.scatter(scale, marker_height, color="#666666", marker="x", s=42, label="SCF failed")
        fit = fit_by_key[(solid, "GFN2")]
        label = "no bracketed minimum" if fit["fit_status"] == "no_local_minimum" else "discontinuous EOS"
        ax.set_title(f"{solid}: {label}")
        ax.set_xlabel("lattice scale (a / experimental a)")
        ax.grid(color="#d0d0d0", linewidth=0.6, alpha=0.7)
    axes[0].set_ylabel("relative energy (eV/atom)")
    handles, labels = axes[0].get_legend_handles_labels()
    unique = dict(zip(labels, handles))
    axes[0].legend(unique.values(), unique.keys(), frameon=False)
    fig.suptitle("LC12 GFN2 EOS diagnostics (native-Bloch k444)")
    fig.tight_layout()
    output = ROOT / "figures" / "goldzak12_gfn2_eos_diagnostics"
    fig.savefig(output.with_suffix(".png"), dpi=220)
    fig.savefig(output.with_suffix(".pdf"))
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cp2k", type=Path, default=base.DEFAULT_CP2K)
    parser.add_argument("--tblite", type=Path, default=base.DEFAULT_TBLITE)
    parser.add_argument("--cp2k-source", type=Path, default=base.DEFAULT_CP2K_SOURCE)
    parser.add_argument("--tblite-source", type=Path, default=base.DEFAULT_TBLITE_SOURCE)
    parser.add_argument("--jobs", type=int, default=6)
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--eos-mesh", default="k444")
    parser.add_argument("--energy-mesh", action="append", default=[])
    parser.add_argument("--result-mesh", default="k444")
    parser.add_argument("--scale", type=float, action="append", default=[])
    args = parser.parse_args()

    scales = tuple(args.scale) if args.scale else DEFAULT_SCALES
    energy_meshes = args.energy_mesh or ["k333", "k444", "k555"]
    if args.result_mesh not in energy_meshes:
        parser.error(f"--result-mesh {args.result_mesh} must also be supplied as --energy-mesh")
    base.write_build_provenance(
        args.cp2k,
        args.tblite,
        args.cp2k_source,
        args.tblite_source,
        {
            "benchmark": "LC12 (Goldzak12)",
            "cell_protocol": "cubic equation of state",
            "eos_mesh": args.eos_mesh,
            "energy_meshes": energy_meshes,
            "result_mesh": args.result_mesh,
            "scales": scales,
            "adaptive_scales": {f"{solid}/{method}": values for (solid, method), values in ADAPTIVE_SCALES.items()},
            "kpoint_scheme": "CP2K native Bloch MACDONALD FULL_GRID",
            "smearing_temperature_K": 300.0,
            "reported_energy": "Total energy extrapolated to T->0",
            "tblite_accuracy": 0.05,
            "default_scf_strategy": "tblite modified-Broyden defaults",
            "scf_retry_strategies": [
                "TBLITE_MIXER ITERATIONS 1200 MEMORY 1 DAMPING 0.05",
                "TBLITE_MIXER ITERATIONS 2400 MEMORY 1 DAMPING 0.01",
            ],
        },
    )
    base.setup_inputs(args.eos_mesh, energy_meshes)
    base.run_tblite_atom_jobs(args.tblite, args.jobs, args.force)
    run_jobs(eos_job_specs(args.eos_mesh, scales), args.cp2k, args.jobs, args.threads, args.force)
    fits = make_eos_table(args.eos_mesh, scales)
    run_jobs(final_sp_specs(fits, energy_meshes), args.cp2k, args.jobs, args.threads, args.force)
    collect_results(fits, energy_meshes, args.result_mesh)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
