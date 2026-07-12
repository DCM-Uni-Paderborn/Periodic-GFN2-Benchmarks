#!/usr/bin/env python3
"""Prepare and analyse k-point dependent DMC-ICE13 GFN benchmarks."""

from __future__ import annotations

import argparse
import csv
import json
import math
import shutil
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
FIGURES = ROOT / "figures"
HARTREE_TO_KJMOL = 2625.499638

PHASES = ["Ih", "II", "III", "IV", "VI", "VII", "VIII", "IX", "XI", "XIII", "XIV", "XV", "XVII"]
METHODS = ["GFN1", "GFN2"]

DMC_ABS_KJMOL = {
    "Ih": -59.45,
    "II": -59.14,
    "III": -58.20,
    "IV": -55.62,
    "VI": -57.67,
    "VII": -54.46,
    "VIII": -55.22,
    "IX": -58.85,
    "XI": -59.29,
    "XIII": -57.33,
    "XIV": -57.75,
    "XV": -57.71,
    "XVII": -57.70,
}

MESHES = [
    {
        "id": "gamma",
        "label": "Gamma",
        "nk_total": 1,
        "scheme": "GAMMA",
        "shift": "",
        "della_pia_role": "Gamma-only reference",
    },
    {
        "id": "k111",
        "label": "1x1x1",
        "nk_total": 1,
        "scheme": "MACDONALD 1 1 1 0.0 0.0 0.0",
        "shift": "0.0 0.0 0.0",
        "della_pia_role": "explicit one-point MacDonald mesh",
    },
    {
        "id": "k222",
        "label": "2x2x2",
        "nk_total": 8,
        "scheme": "MACDONALD 2 2 2 0.25 0.25 0.25",
        "shift": "0.25 0.25 0.25",
        "della_pia_role": "coarse Gamma-centered convergence check",
    },
    {
        "id": "k333",
        "label": "3x3x3",
        "nk_total": 27,
        "scheme": "MACDONALD 3 3 3 0.0 0.0 0.0",
        "shift": "0.0 0.0 0.0",
        "della_pia_role": "GGA, meta-GGA, and vdW single points",
    },
    {
        "id": "k444",
        "label": "4x4x4",
        "nk_total": 64,
        "scheme": "MACDONALD 4 4 4 0.375 0.375 0.375",
        "shift": "0.375 0.375 0.375",
        "della_pia_role": "hybrid-XC single points",
    },
    {
        "id": "k555",
        "label": "5x5x5",
        "nk_total": 125,
        "scheme": "MACDONALD 5 5 5 0.0 0.0 0.0",
        "shift": "0.0 0.0 0.0",
        "della_pia_role": "dense convergence check",
    },
]


def insert_kpoints(input_text: str, mesh: dict[str, object]) -> str:
    block = [
        "    &KPOINTS",
        f"      SCHEME {mesh['scheme']}",
    ]
    if mesh["id"] != "gamma":
        block += [
            "      SYMMETRY T",
            "      FULL_GRID F",
            "      SYMMETRY_BACKEND SPGLIB",
            "      SYMMETRY_REDUCTION_METHOD SPGLIB",
        ]
    block.append("    &END KPOINTS")
    return input_text.replace("    &END QS\n", "    &END QS\n" + "\n".join(block) + "\n", 1)


def enable_cell_canonicalization(input_text: str) -> str:
    return input_text.replace("    &CELL\n", "    &CELL\n      CANONICALIZE TRUE\n", 1)


def prepare_inputs() -> None:
    for mesh in MESHES:
        mesh_id = str(mesh["id"])
        out_dir = ROOT / "kpoint_inputs" / mesh_id
        out_dir.mkdir(parents=True, exist_ok=True)
        for method in METHODS:
            for phase in PHASES:
                base_path = ROOT / "inputs" / f"ice_{phase}_{method}.inp"
                text = base_path.read_text()
                text = enable_cell_canonicalization(text)
                project = f"ice_{phase}_{method}_{mesh_id}"
                text = text.replace(f"PROJECT ice_{phase}_{method}", f"PROJECT {project}")
                if mesh_id != "gamma":
                    text = insert_kpoints(text, mesh)
                else:
                    text = insert_kpoints(text, mesh)
                (out_dir / f"{project}.inp").write_text(text)


def parse_energy(output: Path) -> float | None:
    if not output.exists():
        return None
    energy = None
    for line in output.read_text(errors="ignore").splitlines():
        if "ENERGY| Total FORCE_EVAL" in line:
            energy = float(line.split()[-1])
    return energy


def output_path(mesh_id: str, method: str, phase: str) -> Path:
    if mesh_id == "gamma":
        return ROOT / "runs" / method / phase / f"ice_{phase}_{method}.out"
    return ROOT / "runs_kpoints" / mesh_id / method / phase / f"ice_{phase}_{method}_{mesh_id}.out"


def stats(errors: list[float]) -> dict[str, float]:
    return {
        "ME": sum(errors) / len(errors),
        "MAE": sum(abs(e) for e in errors) / len(errors),
        "RMSE": math.sqrt(sum(e * e for e in errors) / len(errors)),
        "MaxAE": max(abs(e) for e in errors),
    }


def analyse() -> dict[str, object]:
    geometries = json.loads((DATA / "geometries.json").read_text())
    dmc_rel = {phase: DMC_ABS_KJMOL[phase] - DMC_ABS_KJMOL["Ih"] for phase in PHASES}
    results: dict[str, object] = {"meshes": MESHES, "methods": METHODS, "results": {}}
    relative_rows: list[dict[str, object]] = []
    stats_rows: list[dict[str, object]] = []

    for mesh in MESHES:
        mesh_id = str(mesh["id"])
        mesh_results: dict[str, object] = {}
        for method in METHODS:
            energies = {phase: parse_energy(output_path(mesh_id, method, phase)) for phase in PHASES}
            complete = all(value is not None for value in energies.values())
            method_result: dict[str, object] = {"complete": complete, "energies_hartree": energies}
            if complete:
                per_h2o = {
                    phase: float(energies[phase]) / geometries[phase]["counts"]["O"]
                    for phase in PHASES
                }
                ih = per_h2o["Ih"]
                rel = {phase: (per_h2o[phase] - ih) * HARTREE_TO_KJMOL for phase in PHASES}
                err = {phase: rel[phase] - dmc_rel[phase] for phase in PHASES}
                err_nonref = [err[phase] for phase in PHASES if phase != "Ih"]
                s = stats(err_nonref)
                method_result.update(
                    {
                        "per_h2o_hartree": per_h2o,
                        "relative_kjmol": rel,
                        "relative_errors_kjmol": err,
                        "stats_nonreference": s,
                    }
                )
                stats_rows.append(
                    {
                        "mesh": mesh_id,
                        "mesh_label": mesh["label"],
                        "nk_total": mesh["nk_total"],
                        "method": f"{method}-xTB",
                        **{key: f"{value:.6f}" for key, value in s.items()},
                    }
                )
                for phase in PHASES:
                    relative_rows.append(
                        {
                            "mesh": mesh_id,
                            "mesh_label": mesh["label"],
                            "method": f"{method}-xTB",
                            "phase": phase,
                            "DMC_relative_kJmol": f"{dmc_rel[phase]:.6f}",
                            "relative_kJmol": f"{rel[phase]:.6f}",
                            "error_kJmol": f"{err[phase]:.6f}",
                        }
                    )
            mesh_results[method] = method_result
        results["results"][mesh_id] = mesh_results

    DATA.mkdir(exist_ok=True)
    (DATA / "kpoint_results.json").write_text(json.dumps(results, indent=2))
    write_csv(
        DATA / "dmc_ice13_kpoint_relative_energies.csv",
        relative_rows,
        ["mesh", "mesh_label", "method", "phase", "DMC_relative_kJmol", "relative_kJmol", "error_kJmol"],
    )
    write_csv(
        DATA / "dmc_ice13_kpoint_stats.csv",
        stats_rows,
        ["mesh", "mesh_label", "nk_total", "method", "ME", "MAE", "RMSE", "MaxAE"],
    )
    make_plots(stats_rows)
    return results


def write_csv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def make_plots(stats_rows: list[dict[str, object]]) -> None:
    if not stats_rows:
        return
    if shutil.which("gnuplot") is None or shutil.which("rsvg-convert") is None:
        return
    FIGURES.mkdir(exist_ok=True)
    dat = DATA / "dmc_ice13_kpoint_stats_for_plot.dat"
    stats_by_mesh_method = {
        (str(row["mesh"]), str(row["method"])): row for row in stats_rows
    }
    with dat.open("w") as handle:
        handle.write("# index mesh_label GFN1_MAE GFN2_MAE GFN1_RMSE GFN2_RMSE\n")
        for index, mesh in enumerate(MESHES, start=1):
            mesh_id = str(mesh["id"])
            gfn1 = stats_by_mesh_method.get((mesh_id, "GFN1-xTB"))
            gfn2 = stats_by_mesh_method.get((mesh_id, "GFN2-xTB"))
            if gfn1 is None or gfn2 is None:
                continue
            handle.write(
                f"{index} \"{mesh['label']}\" {gfn1['MAE']} {gfn2['MAE']} "
                f"{gfn1['RMSE']} {gfn2['RMSE']}\n"
            )

    script = f"""
set terminal svg enhanced font 'Helvetica,12' size 840,500
set object 1 rectangle from screen 0,0 to screen 1,1 fillcolor rgb 'white' behind
set output '{FIGURES / 'dmc_ice13_kpoint_mae.svg'}'
set border lw 1.2
set tics out nomirror
set grid ytics lc rgb '#d0d0d0' lw 0.6
set key top right spacing 1.2 samplen 2
set xlabel 'k-point mesh'
set ylabel 'Relative-energy error / kJ mol^{{-1}} per H_2O'
set xrange [0.75:6.25]
set yrange [0:*]
set xtics ('Gamma' 1, '1x1x1' 2, '2x2x2' 3, '3x3x3' 4, '4x4x4' 5, '5x5x5' 6)
set style line 1 lc rgb '#c44e52' lw 2.2 pt 7 ps 0.9
set style line 2 lc rgb '#4c72b0' lw 2.2 pt 9 ps 0.9
plot '{dat}' using 1:3 with linespoints ls 1 title 'GFN1-xTB MAE', \\
     '' using 1:4 with linespoints ls 2 title 'GFN2-xTB MAE'
"""
    subprocess.run(["gnuplot"], input=script.encode(), check=True)
    svg = FIGURES / "dmc_ice13_kpoint_mae.svg"
    svg.write_text(
        "\n".join(line.rstrip() for line in svg.read_text().splitlines()) + "\n"
    )
    subprocess.run(["rsvg-convert", str(svg), "-o", str(svg.with_suffix(".png"))], check=True)
    subprocess.run(["rsvg-convert", "-f", "pdf", str(svg), "-o", str(svg.with_suffix(".pdf"))], check=True)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("command", choices=["prepare", "analyse", "all"], nargs="?", default="all")
    args = parser.parse_args()
    if args.command in {"prepare", "all"}:
        prepare_inputs()
    if args.command in {"analyse", "all"}:
        analyse()


if __name__ == "__main__":
    main()
