#!/usr/bin/env python3
"""Build publication-ready Klimes-23 tables and figures from analysed results."""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.patches import Patch


REPO = Path(__file__).resolve().parents[1]
ROOT = REPO / "Klimes-Solids23"
RESULTS = ROOT / "results"
FIGURES = ROOT / "figures"
LITERATURE = (
    "revPBE_vdW",
    "rPW86_vdW2",
    "optPBE_vdW",
    "optB88_vdW",
    "optB86b_vdW",
    "LDA",
    "PBEsol",
    "PBE",
)
METHODS = ("GFN1", "GFN2")
COLORS = {"GFN1": "#C44E52", "GFN2": "#007C91", "DFT": "#777777"}
DISPLAY = {
    "GFN1": "GFN1-xTB",
    "GFN2": "GFN2-xTB",
    "revPBE_vdW": "revPBE-vdW",
    "rPW86_vdW2": "rPW86-vdW2",
    "optPBE_vdW": "optPBE-vdW",
    "optB88_vdW": "optB88-vdW",
    "optB86b_vdW": "optB86b-vdW",
    "LDA": "LDA",
    "PBEsol": "PBEsol",
    "PBE": "PBE",
}
PROPERTIES = {
    "lattice_constant_A": {
        "reference_file": "lattice_constants.csv",
        "reference_column": "experiment_ZPEC",
        "comparison_value": "a0_A",
        "comparison_reference": "a_ref_A",
        "label": "Lattice constant",
        "unit": r"$\AA$",
    },
    "bulk_modulus_GPa": {
        "reference_file": "bulk_moduli.csv",
        "reference_column": "experiment",
        "comparison_value": "B0_GPa",
        "comparison_reference": "B_ref_GPa",
        "label": "Bulk modulus",
        "unit": "GPa",
    },
    "cohesive_energy_eV_atom": {
        "reference_file": "cohesive_energies.csv",
        "reference_column": "experiment_ZPEC",
        "comparison_value": "cohesive_eV_atom",
        "comparison_reference": "cohesive_ref_eV_atom",
        "label": "Cohesive energy",
        "unit": r"eV atom$^{-1}$",
    },
}


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    fields: list[str] = []
    for row in rows:
        for field in row:
            if field not in fields:
                fields.append(field)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def number(value: str | float | None) -> float | None:
    if value is None or str(value).strip() == "":
        return None
    return float(value)


def reference_data() -> tuple[list[str], dict[str, str], dict[str, dict[str, dict[str, str]]]]:
    lattice = read_csv(ROOT / "references/lattice_constants.csv")
    systems = [row["system"] for row in lattice]
    categories = {row["system"]: row["category"] for row in lattice}
    data: dict[str, dict[str, dict[str, str]]] = {}
    for property_name, spec in PROPERTIES.items():
        rows = read_csv(ROOT / "references" / str(spec["reference_file"]))
        data[property_name] = {row["system"]: row for row in rows}
    return systems, categories, data


def build_long_table() -> list[dict[str, object]]:
    systems, categories, references = reference_data()
    comparison = read_csv(RESULTS / "comparison.csv")
    rows: list[dict[str, object]] = []
    for property_name, spec in PROPERTIES.items():
        for system in systems:
            reference_row = references[property_name][system]
            reference = float(reference_row[str(spec["reference_column"])])
            for method in LITERATURE:
                value = float(reference_row[method])
                rows.append(
                    {
                        "system": system,
                        "category": categories[system],
                        "property": property_name,
                        "method": method,
                        "value": value,
                        "reference": reference,
                        "error": value - reference,
                        "relative_error_percent": 100.0 * (value / reference - 1.0),
                        "source": "Klimes et al., Phys. Rev. B 83, 195131 (2011)",
                    }
                )
        for result in comparison:
            value = number(result[str(spec["comparison_value"])])
            reference = number(result[str(spec["comparison_reference"])])
            if value is None or reference is None:
                continue
            rows.append(
                {
                    "system": result["system"],
                    "category": result["category"],
                    "property": property_name,
                    "method": result["method"],
                    "value": value,
                    "reference": reference,
                    "error": value - reference,
                    "relative_error_percent": 100.0 * (value / reference - 1.0),
                    "source": "CP2K native Bloch + tblite",
                }
            )
    write_csv(RESULTS / "literature_comparison_long.csv", rows)
    return rows


def build_system_table() -> list[dict[str, object]]:
    systems, categories, references = reference_data()
    lattice = {row["system"]: row for row in read_csv(ROOT / "references/lattice_constants.csv")}
    comparison = {(row["system"], row["method"]): row for row in read_csv(RESULTS / "comparison.csv")}
    fits = {(row["system"], row["method"]): row for row in read_csv(RESULTS / "eos_fits.csv")}
    rows: list[dict[str, object]] = []
    for system in systems:
        row: dict[str, object] = {
            "system": system,
            "category": categories[system],
            "structure": lattice[system]["structure"],
            "a_exp_ZPEC_A": references["lattice_constant_A"][system]["experiment_ZPEC"],
            "B_exp_GPa": references["bulk_modulus_GPa"][system]["experiment"],
            "Ecoh_exp_ZPEC_eV_atom": references["cohesive_energy_eV_atom"][system]["experiment_ZPEC"],
        }
        for method in METHODS:
            result = comparison.get((system, method))
            fit = fits[(system, method)]
            row[f"{method}_EOS_status"] = "ok" if fit["fit_ok"] == "True" else fit["fit_reason"]
            row[f"{method}_a0_A"] = "" if result is None else result["a0_A"]
            row[f"{method}_B0_GPa"] = "" if result is None else result["B0_GPa"]
            row[f"{method}_Ecoh_eV_atom"] = "" if result is None else result["cohesive_eV_atom"]
        rows.append(row)
    write_csv(RESULTS / "system_summary.csv", rows)
    return rows


def error_statistics(values: np.ndarray, references: np.ndarray) -> dict[str, float]:
    errors = values - references
    relative = 100.0 * errors / references
    return {
        "ME": float(errors.mean()),
        "MAE": float(np.abs(errors).mean()),
        "RMSE": float(np.sqrt(np.mean(errors**2))),
        "MRE_percent": float(relative.mean()),
        "MARE_percent": float(np.abs(relative).mean()),
    }


def build_paired_summary() -> list[dict[str, object]]:
    systems, _, _ = reference_data()
    comparison = {(row["system"], row["method"]): row for row in read_csv(RESULTS / "comparison.csv")}
    paired_systems = [system for system in systems if all((system, method) in comparison for method in METHODS)]
    rows: list[dict[str, object]] = []
    for property_name, spec in PROPERTIES.items():
        references = np.asarray(
            [float(comparison[(system, "GFN1")][str(spec["comparison_reference"])]) for system in paired_systems]
        )
        values = {
            method: np.asarray(
                [float(comparison[(system, method)][str(spec["comparison_value"])]) for system in paired_systems]
            )
            for method in METHODS
        }
        first = error_statistics(values["GFN1"], references)
        second = error_statistics(values["GFN2"], references)
        first_absolute = np.abs(values["GFN1"] - references)
        second_absolute = np.abs(values["GFN2"] - references)
        rows.append(
            {
                "property": property_name,
                "n_paired": len(paired_systems),
                **{f"GFN1_{key}": value for key, value in first.items()},
                **{f"GFN2_{key}": value for key, value in second.items()},
                "GFN2_minus_GFN1_MAE": second["MAE"] - first["MAE"],
                "GFN2_minus_GFN1_MARE_percent_points": second["MARE_percent"] - first["MARE_percent"],
                "systems_GFN2_better": int(np.sum(second_absolute < first_absolute)),
                "systems_GFN2_worse": int(np.sum(second_absolute > first_absolute)),
                "systems_equal": int(np.sum(np.isclose(second_absolute, first_absolute))),
            }
        )
    write_csv(RESULTS / "paired_gfn_statistics.csv", rows)
    return rows


def markdown_tables(system_rows: list[dict[str, object]], paired_rows: list[dict[str, object]]) -> None:
    aggregate = [
        row
        for row in read_csv(RESULTS / "aggregate_statistics.csv")
        if row["method"] in METHODS + LITERATURE
    ]
    method_order = METHODS + LITERATURE
    property_order = tuple(PROPERTIES)
    aggregate_map = {(row["property"], row["method"]): row for row in aggregate}
    lines = [
        "# Klimes-23 benchmark tables",
        "",
        "## Paired GFN1-xTB/GFN2-xTB comparison",
        "",
        "| Property | n | GFN1 MAE | GFN2 MAE | GFN1 MARE (%) | GFN2 MARE (%) | GFN2 better/worse |",
        "|---|---:|---:|---:|---:|---:|---:|",
    ]
    for row in paired_rows:
        lines.append(
            f"| {row['property']} | {row['n_paired']} | {float(row['GFN1_MAE']):.4f} | "
            f"{float(row['GFN2_MAE']):.4f} | {float(row['GFN1_MARE_percent']):.2f} | "
            f"{float(row['GFN2_MARE_percent']):.2f} | {row['systems_GFN2_better']}/{row['systems_GFN2_worse']} |"
        )
    lines.extend(
        [
            "",
            "## Aggregate errors",
            "",
            "| Property | Method | n | ME | MAE | RMSE | MRE (%) | MARE (%) |",
            "|---|---:|---:|---:|---:|---:|---:|---:|",
        ]
    )
    for property_name in property_order:
        for method in method_order:
            row = aggregate_map[(property_name, method)]
            lines.append(
                f"| {property_name} | {method} | {int(float(row['n']))} | "
                f"{float(row['ME']):.4f} | {float(row['MAE']):.4f} | {float(row['RMSE']):.4f} | "
                f"{float(row['MRE_percent']):.2f} | {float(row['MARE_percent']):.2f} |"
            )
    lines.extend(
        [
            "",
            "## Per-system GFN results",
            "",
            "| Solid | a(exp) | a(GFN1) | a(GFN2) | B(exp) | B(GFN1) | B(GFN2) | Ecoh(exp) | Ecoh(GFN1) | Ecoh(GFN2) |",
            "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
        ]
    )
    for row in system_rows:
        def cell(key: str, digits: int) -> str:
            value = number(row[key])
            return "NA" if value is None else f"{value:.{digits}f}"

        lines.append(
            f"| {row['system']} | {cell('a_exp_ZPEC_A', 3)} | {cell('GFN1_a0_A', 3)} | {cell('GFN2_a0_A', 3)} | "
            f"{cell('B_exp_GPa', 1)} | {cell('GFN1_B0_GPa', 1)} | {cell('GFN2_B0_GPa', 1)} | "
            f"{cell('Ecoh_exp_ZPEC_eV_atom', 2)} | {cell('GFN1_Ecoh_eV_atom', 2)} | {cell('GFN2_Ecoh_eV_atom', 2)} |"
        )
    (RESULTS / "benchmark_tables.md").write_text("\n".join(lines) + "\n")


def style_axis(ax: plt.Axes) -> None:
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="y", color="#D0D0D0", linewidth=0.6, alpha=0.7)
    ax.set_axisbelow(True)


def save_figure(fig: plt.Figure, stem: str) -> None:
    FIGURES.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIGURES / f"{stem}.png", dpi=240, bbox_inches="tight", facecolor="white")
    fig.savefig(
        FIGURES / f"{stem}.pdf",
        bbox_inches="tight",
        facecolor="white",
        metadata={"Creator": "Periodic-GFN2-Benchmarks", "CreationDate": None, "ModDate": None},
    )
    plt.close(fig)


def plot_system_errors(long_rows: list[dict[str, object]]) -> None:
    systems, _, _ = reference_data()
    lookup = {
        (str(row["property"]), str(row["method"]), str(row["system"])): float(row["relative_error_percent"])
        for row in long_rows
    }
    x = np.arange(len(systems), dtype=float)
    fig, axes = plt.subplots(3, 1, figsize=(15.5, 10.5), sharex=True, constrained_layout=True)
    for ax, (property_name, spec) in zip(axes, PROPERTIES.items()):
        dft_errors = np.asarray(
            [[lookup[(property_name, method, system)] for method in LITERATURE] for system in systems]
        )
        ax.vlines(x, dft_errors.min(axis=1), dft_errors.max(axis=1), color=COLORS["DFT"], linewidth=3.0, alpha=0.38, label="DFT range")
        ax.scatter(x, np.median(dft_errors, axis=1), marker="_", s=55, color=COLORS["DFT"], linewidth=1.5, label="DFT median")
        for offset, method, marker in [(-0.14, "GFN1", "o"), (0.14, "GFN2", "s")]:
            present_x = []
            values = []
            for index, system in enumerate(systems):
                key = (property_name, method, system)
                if key in lookup:
                    present_x.append(index + offset)
                    values.append(lookup[key])
            ax.scatter(
                present_x,
                values,
                s=34,
                marker=marker,
                color=COLORS[method],
                edgecolor="white",
                linewidth=0.45,
                zorder=3,
                label=f"{DISPLAY[method]} (n={len(values)})",
            )
        ax.axhline(0.0, color="#333333", linewidth=0.8)
        ax.set_ylabel("Relative error (%)")
        ax.set_title(f"{spec['label']} ({spec['unit']})", loc="left", fontweight="normal")
        if property_name == "bulk_modulus_GPa":
            ax.set_yscale("symlog", linthresh=25.0)
        elif property_name == "cohesive_energy_eV_atom":
            ax.set_yscale("symlog", linthresh=20.0)
        style_axis(ax)
    axes[0].legend(ncol=4, frameon=False, loc="upper left")
    axes[-1].set_xticks(x, systems, rotation=55, ha="right")
    axes[-1].set_xlim(-0.7, len(systems) - 0.3)
    fig.suptitle("Klimes-23: system-wise errors against experiment", fontweight="normal")
    save_figure(fig, "klimes23-system-relative-errors")


def plot_aggregate_errors() -> None:
    rows = read_csv(RESULTS / "aggregate_statistics.csv")
    lookup = {
        (row["property"], row["method"]): row
        for row in rows
        if row["method"] in METHODS + LITERATURE
    }
    method_order = METHODS + LITERATURE
    y = np.arange(len(method_order))
    fig, axes = plt.subplots(1, 3, figsize=(15.5, 6.4), sharey=True, constrained_layout=True)
    for ax, (property_name, spec) in zip(axes, PROPERTIES.items()):
        values = [float(lookup[(property_name, method)]["MARE_percent"]) for method in method_order]
        colors = [COLORS.get(method, "#A8A8A8") for method in method_order]
        ax.barh(y, values, color=colors, height=0.68)
        for index, (method, value) in enumerate(zip(method_order, values)):
            n = int(float(lookup[(property_name, method)]["n"]))
            ax.text(value * 1.06, index, f"{value:.1f}% (n={n})", va="center", fontsize=8)
        ax.set_xscale("log")
        ax.set_xlim(max(0.1, min(values) * 0.65), max(values) * 2.8)
        ax.set_xlabel("Mean absolute relative error (%)")
        ax.set_title(str(spec["label"]), loc="left", fontweight="normal")
        ax.grid(axis="x", color="#D0D0D0", linewidth=0.6, alpha=0.7)
        ax.spines[["top", "right", "left"]].set_visible(False)
        ax.tick_params(axis="y", length=0)
    axes[0].set_yticks(y, [DISPLAY[method] for method in method_order])
    axes[0].invert_yaxis()
    fig.suptitle("Klimes-23: aggregate comparison with published DFT", fontweight="normal")
    save_figure(fig, "klimes23-aggregate-mare")


def plot_gfn2_change(long_rows: list[dict[str, object]]) -> None:
    systems, _, _ = reference_data()
    lookup = {
        (str(row["property"]), str(row["method"]), str(row["system"])): float(row["relative_error_percent"])
        for row in long_rows
        if row["method"] in METHODS
    }
    x = np.arange(len(systems))
    fig, axes = plt.subplots(3, 1, figsize=(15.5, 9.5), sharex=True, constrained_layout=True)
    for ax, (property_name, spec) in zip(axes, PROPERTIES.items()):
        changes = []
        colors = []
        missing = []
        for index, system in enumerate(systems):
            first = lookup.get((property_name, "GFN1", system))
            second = lookup.get((property_name, "GFN2", system))
            if first is None or second is None:
                changes.append(np.nan)
                colors.append("#999999")
                missing.append(index)
                continue
            change = abs(second) - abs(first)
            changes.append(change)
            colors.append("#158466" if change < 0 else "#C44E52")
        ax.bar(x, np.nan_to_num(changes), color=colors, width=0.72)
        if missing:
            ax.scatter(missing, np.zeros(len(missing)), marker="x", s=42, color="#555555", zorder=4, label="No paired EOS")
        ax.axhline(0.0, color="#333333", linewidth=0.8)
        ax.set_ylabel("Change in |relative error| (pp)")
        ax.set_title(str(spec["label"]), loc="left", fontweight="normal")
        style_axis(ax)
    axes[0].legend(
        handles=[
            Patch(facecolor="#158466", label="GFN2-xTB lower absolute error"),
            Patch(facecolor="#C44E52", label="GFN2-xTB higher absolute error"),
            Line2D([], [], color="#555555", marker="x", linestyle="None", label="No paired EOS"),
        ],
        frameon=False,
        ncol=3,
        loc="upper right",
    )
    axes[-1].set_xticks(x, systems, rotation=55, ha="right")
    axes[-1].set_xlim(-0.7, len(systems) - 0.3)
    fig.suptitle("GFN2-xTB versus GFN1-xTB: negative values indicate improvement", fontweight="normal")
    save_figure(fig, "klimes23-gfn2-change")


def main() -> int:
    FIGURES.mkdir(parents=True, exist_ok=True)
    long_rows = build_long_table()
    system_rows = build_system_table()
    paired_rows = build_paired_summary()
    markdown_tables(system_rows, paired_rows)
    plot_system_errors(long_rows)
    plot_aggregate_errors()
    plot_gfn2_change(long_rows)
    print(f"wrote {len(long_rows)} long-form rows, {len(system_rows)} system rows, and 3 figures")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
