#!/usr/bin/env python3
from __future__ import annotations

import csv
import importlib.util
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


REPO = Path(__file__).resolve().parents[1]
X23B = REPO / "X23b"
DATA = X23B / "data"
FIGURES = X23B / "figures"
LATTICE_CSV = DATA / "x23b_lattice_energies.csv"
VOLUME_CSV = DATA / "x23b_cell_volumes.csv"
SUMMARY_CSV = DATA / "x23b_summary.csv"
PIPELINE = X23B / "scripts" / "x23b_pipeline.py"
ENERGY_CALCULATION = "cell_opt_single_point"
ENERGY_MESH = "k333"
GEOMETRY_MESH = "k222"
ENERGY_LABEL = "k222 opt + k333 SP"


def load_pipeline_constants():
    spec = importlib.util.spec_from_file_location("x23b_pipeline_constants", PIPELINE)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def read_rows() -> list[dict[str, str]]:
    with LATTICE_CSV.open(newline="") as handle:
        lattice_rows = [
            row
            for row in csv.DictReader(handle)
            if row["calculation"] == ENERGY_CALCULATION and row["mesh"] == ENERGY_MESH
        ]
    with VOLUME_CSV.open(newline="") as handle:
        volumes = {
            (row["system"], row["method"]): row
            for row in csv.DictReader(handle)
            if row["calculation"] == "cell_opt" and row["mesh"] == GEOMETRY_MESH
        }
    if len(lattice_rows) != 46 or len(volumes) != 46:
        raise ValueError(
            "publication plots require 46 final-geometry energy rows and 46 k222 volume rows; "
            f"found {len(lattice_rows)} and {len(volumes)}"
        )
    lattice_keys = {(row["system"], row["method"]) for row in lattice_rows}
    if lattice_keys != set(volumes):
        raise ValueError("final-geometry energy and k222 volume rows do not cover identical systems")
    rows = []
    for row in lattice_rows:
        volume = volumes[(row["system"], row["method"])]
        combined = dict(row)
        combined["volume_error_percent"] = volume["volume_error_percent"]
        rows.append(combined)
    return rows


def read_summary() -> dict[str, dict[str, float]]:
    out: dict[str, dict[str, float]] = {}
    with SUMMARY_CSV.open(newline="") as handle:
        for row in csv.DictReader(handle):
            method = row["method"]
            if (
                row["quantity"] == "lattice_energy_kJmol"
                and row["calculation"] == ENERGY_CALCULATION
                and row["mesh"] == ENERGY_MESH
            ):
                prefix = "lattice_energy"
            elif (
                row["quantity"] == "volume_error_percent"
                and row["calculation"] == "cell_opt"
                and row["mesh"] == GEOMETRY_MESH
            ):
                prefix = "volume"
            else:
                continue
            out.setdefault(method, {})
            for key in ("ME", "MAE", "RMSE", "MaxAE"):
                out[method][f"{prefix}_{key}"] = float(row[key])
    required = {
        f"{prefix}_{key}"
        for prefix in ("lattice_energy", "volume")
        for key in ("ME", "MAE", "RMSE", "MaxAE")
    }
    for method in ("GFN1-xTB", "GFN2-xTB"):
        missing = sorted(required - set(out.get(method, {})))
        if missing:
            raise ValueError(f"incomplete publication summary for {method}: {', '.join(missing)}")
    return out


def save_all(fig: plt.Figure, stem: str) -> None:
    FIGURES.mkdir(parents=True, exist_ok=True)
    for suffix in ("pdf", "png", "svg"):
        fig.savefig(FIGURES / f"{stem}.{suffix}", bbox_inches="tight")


def method_lookup(rows: list[dict[str, str]], field: str) -> dict[tuple[str, str], float]:
    return {(row["system"], row["method"]): float(row[field]) for row in rows}


def make_lattice_profile(rows: list[dict[str, str]], constants) -> None:
    errors = method_lookup(rows, "error_kJmol")
    systems = sorted(constants.SYSTEMS, key=lambda item: float(item["ref_energy"]))
    x = np.arange(1, len(systems) + 1)
    labels = [constants.PLOT_LABELS[str(system["id"])] for system in systems]
    refs = np.array([float(system["ref_energy"]) for system in systems])
    dmc = np.array([constants.DMC_X23[str(system["id"])][0] for system in systems])
    dmc_err = np.array([constants.DMC_X23[str(system["id"])][1] for system in systems])
    mlcc = np.array([constants.MULTILEVEL_CC_X23[str(system["id"])] for system in systems])
    gfn1 = np.array([errors[(str(system["id"]), "GFN1-xTB")] for system in systems])
    gfn2 = np.array([errors[(str(system["id"]), "GFN2-xTB")] for system in systems])

    fig, axes = plt.subplots(2, 1, figsize=(11.2, 7.6), sharex=True, gridspec_kw={"height_ratios": [1, 1.25]})
    ax = axes[0]
    ax.plot(x, refs, "-o", color="#222222", linewidth=1.4, markersize=4.5, label="X23b reference")
    ax.errorbar(x, dmc, yerr=dmc_err, fmt="o", color="#4C78A8", markersize=4.5, linewidth=1.1, capsize=2, label="DMC X23")
    ax.plot(x, mlcc, "-s", color="#7E57C2", linewidth=1.2, markersize=4.2, label="ML-CCSD(T)/RPA+ph")
    ax.set_ylabel("Lattice-energy magnitude / kJ mol$^{-1}$")
    ax.grid(axis="y", color="#d8d8d8", linewidth=0.6)
    ax.legend(frameon=False, loc="upper left", ncol=3)

    ax = axes[1]
    ax.axhspan(-4.184, 4.184, color="#e6e6e6", zorder=0)
    ax.axhline(0.0, color="#555555", linewidth=1.0)
    ax.errorbar(x, dmc - refs, yerr=dmc_err, fmt="o", color="#4C78A8", markersize=4.5, linewidth=1.1, capsize=2, label="DMC X23 - X23b")
    ax.plot(x, mlcc - refs, "-s", color="#7E57C2", linewidth=1.2, markersize=4.2, label="ML-CCSD(T)/RPA+ph - X23b")
    ax.plot(x, gfn1, "-D", color="#E45756", linewidth=1.4, markersize=4.3, label=f"GFN1-xTB {ENERGY_LABEL} - X23b")
    ax.plot(x, gfn2, "-^", color="#54A24B", linewidth=1.4, markersize=4.5, label=f"GFN2-xTB {ENERGY_LABEL} - X23b")
    ax.set_ylabel("Deviation from X23b / kJ mol$^{-1}$")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=55, ha="right", fontsize=9)
    ax.grid(axis="y", color="#d8d8d8", linewidth=0.6)
    ax.legend(frameon=False, loc="upper left", ncol=2)
    fig.tight_layout()
    save_all(fig, "x23b_lattice_energy_prl_style")
    plt.close(fig)


def percentile(values: list[float], fraction: float) -> float:
    return float(np.percentile(np.array(values), 100.0 * fraction))


def make_error_ranges(rows: list[dict[str, str]], constants) -> None:
    gfn1 = [float(row["error_kJmol"]) for row in rows if row["method"] == "GFN1-xTB"]
    gfn2 = [float(row["error_kJmol"]) for row in rows if row["method"] == "GFN2-xTB"]
    mlcc = [constants.MULTILEVEL_CC_X23[str(system["id"])] - float(system["ref_energy"]) for system in constants.SYSTEMS]
    dmc = [constants.DMC_X23[str(system["id"])][0] - float(system["ref_energy"]) for system in constants.SYSTEMS]
    series = [
        (f"GFN1-xTB {ENERGY_LABEL}", gfn1),
        (f"GFN2-xTB {ENERGY_LABEL}", gfn2),
        ("ML-CCSD(T)/RPA+ph", mlcc),
        ("DMC-X23", dmc),
    ]

    fig, ax = plt.subplots(figsize=(9.6, 5.0))
    y = np.arange(len(series), 0, -1)
    ax.axvspan(-4.184, 4.184, color="#f0f0f0", zorder=0)
    ax.axvline(0.0, color="#555555", linewidth=1.1)
    for yi, (label, values) in zip(y, series):
        q1, med, q3 = percentile(values, 0.25), percentile(values, 0.50), percentile(values, 0.75)
        ax.hlines(yi, min(values), max(values), color="#9c9c9c", linewidth=1.6)
        ax.hlines(yi, q1, q3, color="#4C78A8", linewidth=5.0)
        ax.plot(med, yi, "o", color="#4C78A8", markersize=7, label="median" if yi == y[0] else None)
        ax.plot(sum(values) / len(values), yi, "D", color="#E45756", markersize=6, label="mean signed error" if yi == y[0] else None)
        ax.plot(sum(abs(value) for value in values) / len(values), yi, "^", color="#54A24B", markersize=7, label="MAE" if yi == y[0] else None)
    ax.set_yticks(y)
    ax.set_yticklabels([label for label, _ in series])
    ax.set_xlabel("Deviation from X23b / kJ mol$^{-1}$")
    ax.set_ylim(0.4, len(series) + 0.6)
    ax.grid(axis="x", color="#e2e2e2", linewidth=0.7)
    ax.legend(
        frameon=False,
        loc="lower left",
        bbox_to_anchor=(0.0, 1.02),
        ncol=3,
        borderaxespad=0.0,
        handletextpad=0.5,
        columnspacing=1.2,
    )
    fig.tight_layout()
    save_all(fig, "x23b_error_ranges_prl3_style")
    plt.close(fig)


def make_mae_summary(summary: dict[str, dict[str, float]]) -> None:
    methods = ["GFN1-xTB", "GFN2-xTB"]
    energy = [summary[method]["lattice_energy_MAE"] for method in methods]
    volume = [summary[method]["volume_MAE"] for method in methods]
    x = np.arange(len(methods))
    width = 0.36
    fig, ax1 = plt.subplots(figsize=(6.4, 4.3))
    ax2 = ax1.twinx()
    b1 = ax1.bar(x - width / 2, energy, width, color="#E45756", label="Lattice-energy MAE")
    b2 = ax2.bar(x + width / 2, volume, width, color="#54A24B", label="Volume MAE")
    ax1.set_xticks(x)
    ax1.set_xticklabels(methods)
    ax1.set_ylabel("Lattice-energy MAE / kJ mol$^{-1}$")
    ax2.set_ylabel("Volume MAE / %")
    fig.suptitle("X23b native Bloch k222 optimization; k333 energy", y=0.99)
    ax1.grid(axis="y", color="#e2e2e2", linewidth=0.7)
    for bars, axis in [(b1, ax1), (b2, ax2)]:
        for bar in bars:
            axis.text(bar.get_x() + bar.get_width() / 2, bar.get_height(), f"{bar.get_height():.2f}", ha="center", va="bottom", fontsize=9)
    fig.legend(
        [b1, b2],
        ["Lattice-energy MAE", "Volume MAE"],
        frameon=False,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.92),
        ncol=2,
    )
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.83))
    save_all(fig, "x23b_mae_summary")
    plt.close(fig)


def make_volume_comparison(rows: list[dict[str, str]], constants) -> None:
    data: list[tuple[str, list[float]]] = [
        ("GFN1-xTB k222 opt", [float(row["volume_error_percent"]) for row in rows if row["method"] == "GFN1-xTB"]),
        ("GFN2-xTB k222 opt", [float(row["volume_error_percent"]) for row in rows if row["method"] == "GFN2-xTB"]),
    ]
    ref_vol = {str(system["id"]): float(system["ref_volume"]) for system in constants.SYSTEMS}
    for method, volumes in constants.BOESE_DFT_D3_VOLUMES.items():
        data.append((method, [(float(volume) - ref_vol[system]) / ref_vol[system] * 100.0 for system, volume in volumes.items()]))

    fig, ax = plt.subplots(figsize=(9.0, 4.8))
    values = [values for _, values in data]
    labels = [label for label, _ in data]
    bp = ax.boxplot(values, patch_artist=True, vert=True, showmeans=True)
    colors = ["#E45756", "#54A24B", "#9C755F", "#4C78A8", "#F58518"]
    for patch, color in zip(bp["boxes"], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.65)
    ax.axhline(0.0, color="#555555", linewidth=1.0)
    ax.set_ylabel("Cell-volume error / %")
    ax.set_xticklabels(labels, rotation=25, ha="right")
    ax.grid(axis="y", color="#e2e2e2", linewidth=0.7)
    fig.tight_layout()
    save_all(fig, "x23b_volume_comparison_boese")
    plt.close(fig)


def main() -> None:
    constants = load_pipeline_constants()
    rows = read_rows()
    summary = read_summary()
    make_lattice_profile(rows, constants)
    make_error_ranges(rows, constants)
    make_mae_summary(summary)
    make_volume_comparison(rows, constants)


if __name__ == "__main__":
    main()
