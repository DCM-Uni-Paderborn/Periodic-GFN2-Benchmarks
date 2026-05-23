#!/usr/bin/env python3
"""Create DMC-ICE13 tables and plots for the periodic GFN benchmark."""

from __future__ import annotations

import csv
import json
import math
import re
import subprocess
from html import escape
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
FIGURES = ROOT / "figures"

PHASES = ["Ih", "II", "III", "IV", "VI", "VII", "VIII", "IX", "XI", "XIII", "XIV", "XV", "XVII"]

DMC_ABS = {
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

# Della Pia et al. report relative DMC lattice-energy uncertainties of
# 0.10 kJ mol-1 for the non-reference phases in Table I.
DMC_REL_ERROR = {phase: (0.0 if phase == "Ih" else 0.10) for phase in PHASES}

# Published single-point DFT values on the DMC-ICE13 structures, from
# Della Pia et al., J. Chem. Phys. 157, 134701 (2022), Table II.
PUBLISHED_TABLE = """
DMC                  -59.45 -59.14 -58.20 -55.62 -57.67 -54.46 -55.22 -58.85 -59.29 -57.33 -57.75 -57.71 -57.70
B3LYP-D4             -63.03 -61.11 -60.42 -58.10 -58.69 -53.78 -54.98 -61.33 -63.20 -60.18 -59.60 -58.65 -62.05
B3LYP-D3(BJ)atm      -63.63 -61.90 -61.11 -58.90 -59.60 -54.96 -56.21 -62.03 -63.80 -61.00 -60.43 -59.56 -62.68
B3LYP-D3(BJ)         -64.04 -62.52 -61.65 -59.57 -60.37 -55.86 -57.09 -62.62 -64.23 -61.68 -61.17 -60.31 -63.01
B3LYP-D3atm          -64.14 -64.02 -62.10 -60.89 -62.24 -59.89 -60.95 -63.39 -64.28 -62.99 -62.61 -62.16 -63.24
B3LYP-D3             -64.55 -64.65 -62.64 -61.56 -63.01 -60.78 -61.84 -63.98 -64.71 -63.67 -63.34 -62.91 -63.57
B3LYP                -52.71 -46.96 -48.22 -43.85 -42.88 -35.97 -37.36 -48.49 -52.73 -45.72 -44.46 -43.03 -52.09
revPBE0-D4           -56.31 -54.00 -53.65 -51.17 -51.58 -46.87 -48.09 -54.41 -56.42 -53.05 -52.44 -51.62 -55.27
revPBE0-D3(BJ)atm    -56.94 -55.02 -54.50 -52.22 -52.84 -48.38 -49.66 -55.32 -57.06 -54.15 -53.58 -52.86 -55.89
revPBE0-D3(BJ)       -57.35 -55.65 -55.04 -52.89 -53.60 -49.28 -50.55 -55.91 -57.48 -54.83 -54.32 -53.61 -56.22
revPBE0-D3atm        -57.03 -56.96 -54.81 -53.80 -55.43 -54.71 -55.78 -56.03 -57.09 -55.65 -55.34 -55.40 -56.11
revPBE0-D3           -57.45 -57.58 -55.35 -54.47 -56.20 -55.61 -56.67 -56.61 -57.51 -56.33 -56.08 -56.15 -56.44
revPBE0              -44.05 -37.58 -39.51 -34.72 -33.31 -25.88 -27.35 -39.63 -44.00 -36.41 -35.03 -33.58 -43.29
PBE0-MBD             -65.15 -63.08 -62.56 -60.11 -60.74 -55.92 -57.14 -63.40 -65.42 -62.18 -61.54 -60.64 -64.11
PBE0-TS              -64.67 -64.14 -62.66 -61.30 -62.52 -41.88 -58.46 -63.76 -64.95 -63.44 -63.13 -62.43 -63.62
PBE0-D4              -64.38 -61.43 -61.48 -58.49 -58.75 -53.32 -54.56 -62.12 -64.63 -60.48 -59.74 -58.70 -63.34
PBE0-D3(BJ)atm       -64.68 -61.80 -61.81 -58.87 -59.17 -53.88 -55.14 -62.45 -64.93 -60.86 -60.13 -59.12 -63.66
PBE0-D3(BJ)          -65.09 -62.42 -62.35 -59.54 -59.94 -54.78 -56.02 -63.03 -65.35 -61.54 -60.86 -59.87 -63.99
PBE0-D3atm           -65.80 -64.77 -63.44 -61.18 -61.90 -57.86 -58.99 -64.35 -66.05 -63.22 -62.62 -61.82 -64.81
PBE0-D3              -66.21 -64.77 -63.98 -61.85 -62.67 -58.75 -59.88 -64.94 -66.47 -63.91 -63.35 -62.57 -65.14
PBE0                 -57.47 -52.12 -53.27 -49.05 -48.39 -41.88 -43.23 -53.49 -57.63 -50.88 -49.74 -48.45 -56.77
SCAN+rVV10           -68.26 -67.83 -66.33 -65.17 -66.27 -64.07 -65.53 -67.38 -68.52 -67.14 -66.90 -66.61 -67.29
R2SCAN               -62.83 -61.34 -60.39 -58.60 -59.53 -56.39 -57.61 -61.00 -63.08 -60.40 -60.04 -59.43 -62.10
RSCAN                -61.40 -59.38 -58.67 -56.51 -57.35 -53.83 -55.04 -59.15 -61.68 -58.35 -57.90 -57.21 -60.75
SCAN                 -64.42 -62.84 -61.95 -60.05 -60.58 -57.58 -59.09 -62.77 -64.80 -61.95 -61.46 -61.00 -63.69
optB88-vdW           -67.90 -67.85 -66.91 -65.66 -66.99 -63.50 -64.61 -67.88 -68.32 -67.66 -67.52 -66.82 -66.40
optB86b-vdW          -68.69 -67.89 -67.47 -65.74 -66.83 -62.84 -63.93 -68.22 -69.18 -67.68 -67.43 -66.63 -67.21
optPBE-vdW           -65.55 -65.72 -64.79 -63.77 -64.84 -61.25 -62.46 -65.93 -65.82 -65.66 -65.50 -64.83 -63.87
rev-vdW-DF2          -66.37 -64.26 -64.27 -61.87 -62.52 -57.90 -59.03 -64.84 -66.84 -63.79 -63.32 -62.35 -65.26
vdW-DF2              -59.43 -60.16 -58.40 -57.93 -59.18 -56.49 -57.80 -59.87 -59.43 -59.92 -59.84 -59.30 -58.06
vdW-DF               -53.59 -53.78 -52.86 -52.03 -52.74 -49.13 -50.47 -54.19 -53.59 -53.81 -53.60 -52.96 -51.84
revPBE-D4            -59.96 -56.47 -57.12 -54.00 -53.97 -48.67 -49.73 -57.33 -60.29 -55.68 -54.98 -53.90 -58.96
revPBE-D3(BJ)atm     -58.86 -55.92 -56.39 -53.54 -53.66 -48.25 -49.43 -56.74 -59.17 -55.31 -54.64 -53.62 -57.80
revPBE-D3(BJ)        -59.27 -56.55 -56.92 -54.21 -54.42 -49.15 -50.31 -57.33 -59.59 -55.99 -55.37 -54.38 -58.13
revPBE-D3atm         -58.60 -57.12 -56.15 -54.37 -55.40 -53.93 -54.85 -56.82 -58.83 -56.02 -55.57 -55.31 -57.68
revPBE-D3            -59.01 -57.75 -56.69 -55.04 -56.16 -54.83 -55.74 -57.41 -59.25 -56.71 -56.30 -56.07 -58.00
revPBE               -43.86 -35.59 -38.96 -33.16 -30.88 -21.89 -23.30 -38.52 -43.97 -34.66 -33.02 -31.14 -43.11
PBE-MBD              -70.48 -67.29 -67.77 -64.69 -64.75 -58.87 -59.94 -68.13 -70.96 -66.64 -65.83 -64.57 -69.40
PBE-TS               -69.59 -67.75 -67.44 -65.29 -65.86 -58.84 -60.36 -68.04 -70.10 -67.30 -66.81 -65.71 -68.53
PBE-D4               -69.62 -65.55 -66.60 -62.97 -62.71 -56.21 -57.35 -66.75 -70.10 -64.82 -63.96 -62.58 -68.61
PBE-D3(BJ)atm        -69.95 -65.95 -66.96 -63.38 -63.19 -56.85 -58.02 -67.11 -70.43 -65.23 -64.38 -63.05 -68.96
PBE-D3(BJ)           -70.36 -66.58 -67.49 -64.05 -63.96 -57.75 -58.90 -67.69 -70.85 -65.92 -65.11 -63.80 -69.29
PBE-D3atm            -70.39 -67.75 -67.94 -65.15 -65.41 -59.76 -60.90 -68.38 -70.87 -67.07 -66.37 -65.25 -69.42
PBE-D3               -70.80 -68.37 -68.47 -65.82 -66.17 -60.66 -61.79 -68.97 -71.29 -67.76 -67.11 -66.00 -69.75
PBE                  -62.23 -55.54 -57.85 -52.85 -51.58 -43.79 -45.05 -57.56 -62.60 -54.54 -53.24 -51.56 -61.52
HF-D4                -42.89 -46.50 -41.70 -42.92 -45.69 -46.31 -47.99 -44.36 -42.18 -45.14 -45.47 -46.08 -41.82
HF-D3(BJ)atm         -50.53 -53.57 -49.36 -50.04 -53.03 -53.70 -55.15 -51.58 -49.97 -52.28 -52.77 -53.25 -49.53
HF-D3(BJ)            -50.94 -54.20 -49.90 -50.70 -53.79 -54.60 -56.03 -52.16 -50.39 -52.97 -53.50 -54.00 -49.86
HF-D3atm             -39.07 -44.94 -38.65 -41.32 -44.88 -47.30 -49.09 -41.88 -38.27 -43.58 -44.14 -45.31 -38.01
HF-D3                -39.49 -45.57 -39.19 -41.99 -45.64 -48.20 -49.97 -42.47 -38.69 -44.27 -44.88 -46.06 -38.34
HF                   -26.57 -25.53 -23.38 -22.01 -22.58 -19.54 -21.57 -25.42 -25.62 -24.06 -23.47 -23.29 -25.73
LDA                  -100.08 -94.37 -95.92 -90.94 -91.03 -83.78 -84.66 -95.37 -101.22 -93.04 -91.99 -90.32 -99.67
"""


def rel_from_abs(values: list[float]) -> dict[str, float]:
    ref = values[0]
    return {phase: value - ref for phase, value in zip(PHASES, values)}


def parse_published_abs() -> dict[str, list[float]]:
    data: dict[str, list[float]] = {}
    for line in PUBLISHED_TABLE.splitlines():
        stripped = line.strip()
        if not stripped:
            continue
        values = [float(x) for x in re.findall(r"-?\d+\.\d+", stripped)]
        if len(values) != len(PHASES):
            raise ValueError(f"Could not parse DMC-ICE13 table line: {line}")
        name = stripped[: stripped.index(f"{values[0]:.2f}")].strip()
        data[name] = values
    return data


def stats(errors: list[float]) -> dict[str, float]:
    return {
        "ME": sum(errors) / len(errors),
        "MAE": sum(abs(e) for e in errors) / len(errors),
        "RMSE": math.sqrt(sum(e * e for e in errors) / len(errors)),
        "MaxAE": max(abs(e) for e in errors),
    }


def load_primary_gfn_relative_energies(results: dict[str, object]) -> dict[str, dict[str, float]]:
    kpoint_path = DATA / "kpoint_results.json"
    if kpoint_path.exists():
        kpoint_results = json.loads(kpoint_path.read_text())
        mesh = "k333"
        mesh_results = kpoint_results["results"][mesh]
        return {
            "GFN1-xTB": mesh_results["GFN1"]["relative_kjmol"],
            "GFN2-xTB": mesh_results["GFN2"]["relative_kjmol"],
        }
    return {
        "GFN1-xTB": results["GFN1"]["relative_kjmol"],
        "GFN2-xTB": results["GFN2"]["relative_kjmol"],
    }


def load_gamma_gfn_relative_energies(results: dict[str, object]) -> dict[str, dict[str, float]]:
    kpoint_path = DATA / "kpoint_results.json"
    if kpoint_path.exists():
        kpoint_results = json.loads(kpoint_path.read_text())
        mesh_results = kpoint_results["results"]["gamma"]
        return {
            "GFN1-xTB (Gamma)": mesh_results["GFN1"]["relative_kjmol"],
            "GFN2-xTB (Gamma)": mesh_results["GFN2"]["relative_kjmol"],
        }
    return {
        "GFN1-xTB (Gamma)": results["GFN1"]["relative_kjmol"],
        "GFN2-xTB (Gamma)": results["GFN2"]["relative_kjmol"],
    }


def write_csv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def run_gnuplot(script: str) -> None:
    subprocess.run(["gnuplot"], input=script.encode(), check=True)


def svg_attr(attrs: dict[str, object]) -> str:
    converted = []
    for key, value in attrs.items():
        name = "class" if key == "class_" else key.replace("_", "-")
        converted.append(f'{name}="{escape(str(value))}"')
    return " ".join(converted)


def svg_text(x: float, y: float, content: str, **attrs: object) -> str:
    return f'<text x="{x:.1f}" y="{y:.1f}" {svg_attr(attrs)}>{escape(content)}</text>'


def svg_line(x1: float, y1: float, x2: float, y2: float, **attrs: object) -> str:
    return f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" {svg_attr(attrs)}/>'


def svg_circle(x: float, y: float, r: float, **attrs: object) -> str:
    return f'<circle cx="{x:.1f}" cy="{y:.1f}" r="{r:.1f}" {svg_attr(attrs)}/>'


def svg_rect(x: float, y: float, width: float, height: float, **attrs: object) -> str:
    return f'<rect x="{x:.1f}" y="{y:.1f}" width="{width:.1f}" height="{height:.1f}" {svg_attr(attrs)}/>'


def make_log_mae_plot(summary_rows: list[dict[str, object]]) -> None:
    rows = sorted(summary_rows, key=lambda row: float(row["MAE"]))
    width, height = 1280, 1600
    left, right, top, bottom = 285, 70, 112, 120
    plot_w = width - left - right
    row_gap = (height - top - bottom) / (len(rows) + 1)
    xmin, xmax = math.log10(0.3), math.log10(15.0)

    def xlog(value: float) -> float:
        return left + (math.log10(value) - xmin) * plot_w / (xmax - xmin)

    color_by_method = {
        "GFN2-xTB": "#4c72b0",
        "GFN2-xTB (Gamma)": "#8fb0df",
        "GFN1-xTB": "#c44e52",
        "GFN1-xTB (Gamma)": "#e6a0a3",
        "PBE-D4": "#55a868",
        "PBE0-D4": "#8172b3",
        "SCAN+rVV10": "#ccb974",
    }
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        """
<style>
  text { font-family: Helvetica, Arial, sans-serif; fill: #1f2933; }
  .title { font-size: 26px; font-weight: 700; }
  .subtitle { font-size: 15px; fill: #667085; }
  .axis { font-size: 16px; font-weight: 600; fill: #344054; }
  .tick { font-size: 13px; fill: #53616d; }
  .method { font-size: 12px; fill: #53616d; }
  .method-hi { font-size: 12px; font-weight: 700; fill: #1f2933; }
  .value { font-size: 11px; fill: #344054; }
  .value-hi { font-size: 11px; font-weight: 700; fill: #1f2933; }
</style>
""",
        svg_rect(0, 0, width, height, fill="#ffffff"),
        svg_text(56, 44, "DMC-ICE13 relative-energy MAE ranking", class_="title"),
        svg_text(
            56,
            70,
            "All values are relative to ice Ih; MAE is evaluated over the 12 non-reference polymorphs.",
            class_="subtitle",
        ),
        svg_rect(left, top, plot_w, height - top - bottom, fill="#fbfcfd", stroke="#d5d8dc", stroke_width="1.1"),
    ]

    for tick in [0.3, 0.5, 1, 2, 5, 10, 15]:
        tx = xlog(tick)
        parts.append(svg_line(tx, top, tx, height - bottom, stroke="#e5e7eb", stroke_width="1"))
        parts.append(svg_text(tx, height - bottom + 27, f"{tick:g}", text_anchor="middle", class_="tick"))
    for tick, label in [(1, "1"), (10, "10")]:
        tx = xlog(tick)
        parts.append(svg_line(tx, top, tx, height - bottom, stroke="#9aa5b1", stroke_width="1.4", stroke_dasharray="4 5"))
        parts.append(svg_text(tx + 5, top + 18, f"{label} kJ mol-1", class_="tick"))

    parts.append(svg_text(left + plot_w / 2, height - 52, "Relative-energy MAE vs DMC / kJ mol-1 (log scale)", text_anchor="middle", class_="axis"))
    parts.append(svg_text(56, top - 14, "Method", class_="axis"))
    parts.append(svg_text(width - right, top - 14, "MAE", text_anchor="end", class_="axis"))

    for index, row in enumerate(rows, start=1):
        method = str(row["method"])
        value = float(row["MAE"])
        y = top + row_gap * index
        is_primary_gfn = method in {"GFN1-xTB", "GFN2-xTB"}
        is_gamma_gfn = method in {"GFN1-xTB (Gamma)", "GFN2-xTB (Gamma)"}
        is_hi = is_primary_gfn or is_gamma_gfn
        is_selected = method in color_by_method
        method_label = method.replace("(Gamma)", "(\u0393-point)")
        color = color_by_method.get(method, "#9aa5b1")
        opacity = "0.95" if is_primary_gfn else "0.70" if is_gamma_gfn else "0.78" if is_selected else "0.46"
        width_line = "7.0" if is_primary_gfn else "4.6" if is_gamma_gfn else "5.0" if is_selected else "3.2"
        line_attrs = {
            "stroke": color,
            "stroke_width": width_line,
            "stroke_linecap": "round",
            "stroke_opacity": opacity,
        }
        if is_gamma_gfn:
            line_attrs["stroke_dasharray"] = "8 7"
        parts.append(svg_line(left, y, xlog(value), y, **line_attrs))
        parts.append(svg_circle(xlog(value), y, 5.7 if is_primary_gfn else 4.8 if is_gamma_gfn else 4.4 if is_selected else 3.1, fill="#ffffff", stroke=color, stroke_width=2.0 if is_hi else 1.4, stroke_opacity=opacity))
        parts.append(svg_text(left - 12, y + 4, method_label, text_anchor="end", class_="method-hi" if is_hi else "method"))
        parts.append(svg_text(xlog(value) + 9, y + 4, f"{value:.2f}", class_="value-hi" if is_hi else "value"))

    parts.append(svg_text(56, height - 20, "Published DFT data from Della Pia et al., J. Chem. Phys. 157, 134701 (2022); GFN entries from CP2K/tblite single points in this work.", class_="subtitle"))
    parts.append("</svg>\n")
    (FIGURES / "dmc_ice13_relative_mae_all_methods_log.svg").write_text("\n".join(parts))


def main() -> None:
    DATA.mkdir(exist_ok=True)
    FIGURES.mkdir(exist_ok=True)
    results = json.loads((DATA / "results.json").read_text())
    published_abs = parse_published_abs()

    dmc_rel = rel_from_abs(published_abs["DMC"])
    method_rel = load_primary_gfn_relative_energies(results)
    gamma_rel = load_gamma_gfn_relative_energies(results)
    method_rel.update({name: rel_from_abs(vals) for name, vals in published_abs.items() if name != "DMC"})

    published_rows = []
    for name, values in published_abs.items():
        row = {"method": name}
        row.update({phase: f"{value:.4f}" for phase, value in zip(PHASES, values)})
        published_rows.append(row)
    write_csv(DATA / "dmc_ice13_published_dft_absolute_energies.csv", published_rows, ["method", *PHASES])

    rows = []
    for index, phase in enumerate(PHASES, start=1):
        row = {
            "index": index,
            "phase": phase,
            "DMC_relative_kJmol": f"{dmc_rel[phase]:.4f}",
            "GFN1_relative_kJmol": f"{method_rel['GFN1-xTB'][phase]:.4f}",
            "GFN2_relative_kJmol": f"{method_rel['GFN2-xTB'][phase]:.4f}",
            "GFN1_error_kJmol": f"{method_rel['GFN1-xTB'][phase] - dmc_rel[phase]:.4f}",
            "GFN2_error_kJmol": f"{method_rel['GFN2-xTB'][phase] - dmc_rel[phase]:.4f}",
        }
        rows.append(row)
    write_csv(
        DATA / "dmc_ice13_relative_energies.csv",
        rows,
        ["index", "phase", "DMC_relative_kJmol", "GFN1_relative_kJmol", "GFN2_relative_kJmol", "GFN1_error_kJmol", "GFN2_error_kJmol"],
    )

    summary_rows = []
    summary_rel = dict(method_rel)
    summary_rel.update(gamma_rel)
    for method, rel in summary_rel.items():
        errors = [rel[phase] - dmc_rel[phase] for phase in PHASES if phase != "Ih"]
        row = {"method": method}
        row.update({key: f"{value:.4f}" for key, value in stats(errors).items()})
        summary_rows.append(row)
    summary_rows.sort(key=lambda row: float(row["MAE"]))
    write_csv(DATA / "dmc_ice13_relative_mae_comparison.csv", summary_rows, ["method", "ME", "MAE", "RMSE", "MaxAE"])
    make_log_mae_plot(summary_rows)

    rel_dat = DATA / "relative_energies_for_plot.dat"
    with rel_dat.open("w") as handle:
        handle.write("# index phase DMC DMC_error GFN1 GFN2 GFN1_Gamma GFN2_Gamma PBE-D4 PBE0-D4 SCAN+rVV10\n")
        for index, phase in enumerate(PHASES, start=1):
            handle.write(
                f"{index} {phase} {dmc_rel[phase]:.6f} {DMC_REL_ERROR[phase]:.6f} "
                f"{method_rel['GFN1-xTB'][phase]:.6f} {method_rel['GFN2-xTB'][phase]:.6f} "
                f"{gamma_rel['GFN1-xTB (Gamma)'][phase]:.6f} {gamma_rel['GFN2-xTB (Gamma)'][phase]:.6f} "
                f"{method_rel['PBE-D4'][phase]:.6f} {method_rel['PBE0-D4'][phase]:.6f} "
                f"{method_rel['SCAN+rVV10'][phase]:.6f}\n"
            )

    err_dat = DATA / "relative_errors_for_plot.dat"
    with err_dat.open("w") as handle:
        handle.write("# index phase GFN1_error GFN2_error\n")
        for index, phase in enumerate(PHASES, start=1):
            handle.write(
                f"{index} {phase} {method_rel['GFN1-xTB'][phase] - dmc_rel[phase]:.6f} "
                f"{method_rel['GFN2-xTB'][phase] - dmc_rel[phase]:.6f}\n"
            )

    mae_dat = DATA / "relative_mae_for_plot.dat"
    compact_methods = {
        "B3LYP-D3atm",
        "SCAN+rVV10",
        "RSCAN",
        "PBE-D4",
        "PBE0-D4",
        "GFN2-xTB",
        "GFN1-xTB",
    }
    compact_rows = [row for row in summary_rows if row["method"] in compact_methods]
    with mae_dat.open("w") as handle:
        handle.write("# index method MAE\n")
        for index, row in enumerate(compact_rows, start=1):
            handle.write(f"{index} \"{row['method']}\" {row['MAE']}\n")

    common = """
	set terminal svg enhanced dashed font 'Helvetica,14' size 1180,640
set object 1 rectangle from screen 0,0 to screen 1,1 fillcolor rgb 'white' behind
set border lw 1.3
set tics out nomirror
set grid ytics lc rgb '#d0d0d0' lw 0.6
set bars 0.35
set key outside right top spacing 1.15 samplen 1.8
set rmargin 25
set lmargin 9
set style line 1 lt 1 lc rgb '#111111' lw 2.6 pt 7 ps 0.85
set style line 2 lt 1 lc rgb '#c44e52' lw 2.3 pt 5 ps 0.75
set style line 3 lt 1 lc rgb '#4c72b0' lw 2.3 pt 9 ps 0.75
	set style line 4 lt 1 lc rgb '#55a868' lw 1.7 pt 13 ps 0.65
	set style line 5 lt 1 lc rgb '#8172b3' lw 1.7 pt 11 ps 0.65
	set style line 6 lt 1 lc rgb '#ccb974' lw 1.7 pt 15 ps 0.65
	set style line 7 lt 2 lc rgb '#e6a0a3' lw 1.6 pt 5 ps 0.55
	set style line 8 lt 2 lc rgb '#8fb0df' lw 1.6 pt 9 ps 0.55
	"""
    run_gnuplot(
        common
        + f"""
set output '{FIGURES / 'dmc_ice13_relative_energies.svg'}'
set ylabel 'Relative energy to ice Ih / kJ mol^{-1}'
set xlabel 'Ice polymorph'
set xrange [0.5:13.5]
set xtics rotate by -45
set yrange [-25:25]
set ytics 5
plot '{rel_dat}' using 1:3:4:xtic(2) w yerrorlines ls 1 title 'DMC', \\
	     '' using 1:5 w lp ls 2 title 'GFN1-xTB', \\
	     '' using 1:6 w lp ls 3 title 'GFN2-xTB', \\
	     '' using 1:7 w lp ls 7 title 'GFN1-xTB (\u0393-point)', \\
	     '' using 1:8 w lp ls 8 title 'GFN2-xTB (\u0393-point)', \\
	     '' using 1:9 w lp ls 4 title 'PBE-D4', \\
	     '' using 1:10 w lp ls 5 title 'PBE0-D4', \\
	     '' using 1:11 w lp ls 6 title 'SCAN+rVV10'
"""
    )

    run_gnuplot(
        common
        + f"""
set output '{FIGURES / 'dmc_ice13_relative_errors.svg'}'
set ylabel 'Relative-energy error vs DMC / kJ mol^{-1}'
set xlabel 'Ice polymorph'
set xrange [0.5:13.5]
set xtics rotate by -45
set yrange [-24:22]
set yzeroaxis lw 1.2 lc rgb '#333333'
plot '{err_dat}' using 1:3:xtic(2) w lp ls 2 title 'GFN1-xTB', \\
     '' using 1:4 w lp ls 3 title 'GFN2-xTB'
"""
    )

    run_gnuplot(
        """
set terminal svg enhanced font 'Helvetica,12' size 720,480
set object 1 rectangle from screen 0,0 to screen 1,1 fillcolor rgb 'white' behind
set output '"""
        + str(FIGURES / "dmc_ice13_relative_mae_comparison.svg")
        + """'
set border lw 1.2
set tics out nomirror
set grid ytics lc rgb '#d0d0d0' lw 0.6
set style fill solid 0.85 border -1
set boxwidth 0.72
set ylabel 'MAE of relative energies / kJ mol^{-1}'
set xlabel 'Method'
set xrange [0.3:10.7]
set yrange [0:12]
set xtics rotate by -45
unset key
plot '"""
        + str(mae_dat)
        + """' using 1:3:xtic(2) with boxes lc rgb '#4c72b0'
"""
    )

    for svg in FIGURES.glob("*.svg"):
        png = svg.with_suffix(".png")
        pdf = svg.with_suffix(".pdf")
        subprocess.run(["rsvg-convert", str(svg), "-o", str(png)], check=True)
        subprocess.run(["rsvg-convert", "-f", "pdf", str(svg), "-o", str(pdf)], check=True)


if __name__ == "__main__":
    main()
