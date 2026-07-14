#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import math
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
FIGURES = ROOT / "figures"
COMMON5 = ("C", "Si", "SiC", "LiF", "LiCl")

FAMILY_COLORS = {
    "GFN-xTB": "#E45756",
    "DFT": "#4C78A8",
    "post-HF": "#54A24B",
}


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def value(text: str | None) -> float | None:
    if text is None or text.strip() == "":
        return None
    return float(text)


def add_record(
    records: list[dict[str, object]],
    source: str,
    family: str,
    method: str,
    solid: str,
    a_calc: float | None,
    ecoh_calc: float | None,
) -> None:
    records.append(
        {
            "source": source,
            "family": family,
            "method": method,
            "solid": solid,
            "a_calc_A": a_calc,
            "ecoh_calc_eV_per_atom": ecoh_calc,
        }
    )


def provenance_result_mesh() -> str:
    provenance = json.loads((DATA / "build_provenance.json").read_text())
    return str(provenance["protocol"]["result_mesh"])


def collect_records(result_mesh: str) -> tuple[list[dict[str, object]], dict[str, dict[str, float]]]:
    goldzak = read_csv(DATA / "reference_goldzak2022.csv")
    refs = {
        row["solid"]: {
            "a": float(row["a_exp_A"]),
            "ecoh": float(row["ecoh_exp_eV_per_atom"]),
        }
        for row in goldzak
    }
    records: list[dict[str, object]] = []

    for row in read_csv(DATA / "eos_results.csv"):
        if row["energy_mesh"] != result_mesh or row["sp_completed"] != "True":
            continue
        add_record(
            records,
            f"This work: CP2K/tblite {result_mesh} EOS",
            "GFN-xTB",
            row["method"],
            row["solid"],
            value(row["a_calc_A"]),
            value(row["ecoh_calc_eV_per_atom"]),
        )

    post_hf = {
        "HF": ("a_HF_A", "ecoh_HF_eV_per_atom"),
        "MP2": ("a_MP2_A", "ecoh_MP2_eV_per_atom"),
        "SCS-MP2": ("a_SCS_MP2_A", "ecoh_SCS_MP2_eV_per_atom"),
        "SOS-MP2": ("a_SOS_MP2_A", "ecoh_SOS_MP2_eV_per_atom"),
    }
    for row in goldzak:
        for method, (a_key, e_key) in post_hf.items():
            add_record(
                records,
                "Goldzak et al. 2022",
                "post-HF",
                method,
                row["solid"],
                value(row[a_key]),
                value(row[e_key]),
            )

    mejia_methods = {
        "SCAN": ("a_SCAN_A", "ecoh_SCAN_eV_per_atom"),
        "SCAN-L": ("a_SCAN_L_A", "ecoh_SCAN_L_eV_per_atom"),
        "r2SCAN": ("a_r2SCAN_A", "ecoh_r2SCAN_eV_per_atom"),
        "r2SCAN-L": ("a_r2SCAN_L_A", "ecoh_r2SCAN_L_eV_per_atom"),
    }
    for row in read_csv(DATA / "reference_dft_mejia2020.csv"):
        for method, (a_key, e_key) in mejia_methods.items():
            add_record(
                records,
                "Mejia-Rodriguez and Trickey 2020",
                "DFT",
                method,
                row["solid"],
                value(row[a_key]),
                value(row[e_key]),
            )

    mo_methods = {
        "LSDA": ("a_LSDA_A", "ecoh_LSDA_eV_per_atom"),
        "PBE": ("a_PBE_A", "ecoh_PBE_eV_per_atom"),
        "PBEsol": ("a_PBEsol_A", "ecoh_PBEsol_eV_per_atom"),
        "TPSS": ("a_TPSS_A", "ecoh_TPSS_eV_per_atom"),
        "revTPSS": ("a_revTPSS_A", "ecoh_revTPSS_eV_per_atom"),
        "TM": ("a_TM_A", "ecoh_TM_eV_per_atom"),
        "HSE06": ("a_HSE06_A", "ecoh_HSE06_eV_per_atom"),
        "optB86b-vdW": ("a_optB86b_vdW_A", "ecoh_optB86b_vdW_eV_per_atom"),
    }
    for row in read_csv(DATA / "reference_dft_mo2017.csv"):
        for method, (a_key, e_key) in mo_methods.items():
            a_calc = value(row[a_key])
            ecoh_calc = value(row[e_key])
            if a_calc is None and ecoh_calc is None:
                continue
            add_record(
                records,
                "Mo et al. 2017",
                "DFT",
                method,
                row["solid"],
                a_calc,
                ecoh_calc,
            )

    for row in records:
        ref = refs[str(row["solid"])]
        a_calc = row["a_calc_A"]
        ecoh_calc = row["ecoh_calc_eV_per_atom"]
        row["a_error_A"] = None if a_calc is None else float(a_calc) - ref["a"]
        row["ecoh_error_eV_per_atom"] = None if ecoh_calc is None else float(ecoh_calc) - ref["ecoh"]
    return records, refs


def metric(values: list[float]) -> tuple[float, float, float]:
    array = np.asarray(values, dtype=float)
    return float(array.mean()), float(np.abs(array).mean()), float(np.sqrt(np.mean(array * array)))


def summarize(records: list[dict[str, object]], subset: tuple[str, ...] | None = None) -> list[dict[str, object]]:
    grouped: dict[tuple[str, str, str], list[dict[str, object]]] = defaultdict(list)
    for row in records:
        if subset is None or row["solid"] in subset:
            grouped[(str(row["source"]), str(row["family"]), str(row["method"]))].append(row)

    summary: list[dict[str, object]] = []
    for (source, family, method), rows in grouped.items():
        a_err = [float(row["a_error_A"]) for row in rows if row["a_error_A"] is not None]
        e_err = [float(row["ecoh_error_eV_per_atom"]) for row in rows if row["ecoh_error_eV_per_atom"] is not None]
        if subset is not None:
            if len(a_err) != len(subset) or len(e_err) != len(subset):
                continue
        a_me, a_mae, a_rmse = metric(a_err) if a_err else (math.nan, math.nan, math.nan)
        e_me, e_mae, e_rmse = metric(e_err) if e_err else (math.nan, math.nan, math.nan)
        summary.append(
            {
                "source": source,
                "family": family,
                "method": method,
                "n_a": len(a_err),
                "a_ME_A": a_me,
                "a_MAE_A": a_mae,
                "a_RMSE_A": a_rmse,
                "n_ecoh": len(e_err),
                "ecoh_ME_eV_per_atom": e_me,
                "ecoh_MAE_eV_per_atom": e_mae,
                "ecoh_RMSE_eV_per_atom": e_rmse,
            }
        )
    family_order = {"GFN-xTB": 0, "DFT": 1, "post-HF": 2}
    return sorted(summary, key=lambda row: (family_order[str(row["family"])], str(row["method"])))


def csv_value(item: object) -> object:
    if isinstance(item, float):
        return "" if math.isnan(item) else f"{item:.8f}"
    if item is None:
        return ""
    return item


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: csv_value(item) for key, item in row.items()})


def fmt(value_: object, digits: int = 4) -> str:
    if value_ is None:
        return ""
    number = float(value_)
    if math.isnan(number):
        return ""
    return f"{number:.{digits}f}"


def value_error(row: dict[str, object] | None, value_key: str, error_key: str, digits: int) -> str:
    if row is None or row[value_key] is None:
        return ""
    return f"{float(row[value_key]):.{digits}f} ({float(row[error_key]):+.{digits}f})"


def write_markdown(
    records: list[dict[str, object]],
    available: list[dict[str, object]],
    common5: list[dict[str, object]],
) -> None:
    lines = [
        "# LC10 literature comparison",
        "",
        "All errors below use the same zero-point-corrected experimental values from Goldzak et al. (2022).",
        "The source-specific experimental columns are retained in the raw reference CSV files but are not mixed into the MAEs.",
        "",
        "## Available literature coverage",
        "",
        "| family | method | n a | a ME (A) | a MAE (A) | n Ecoh | Ecoh ME (eV/atom) | Ecoh MAE (eV/atom) |",
        "|---|---|---:|---:|---:|---:|---:|---:|",
    ]
    for row in available:
        lines.append(
            f"| {row['family']} | {row['method']} | {row['n_a']} | {fmt(row['a_ME_A'])} | "
            f"{fmt(row['a_MAE_A'])} | {row['n_ecoh']} | {fmt(row['ecoh_ME_eV_per_atom'])} | "
            f"{fmt(row['ecoh_MAE_eV_per_atom'])} |"
        )
    lines += [
        "",
        "## Common five-system comparison",
        "",
        "Common subset: C, Si, SiC, LiF, and LiCl. Every listed method has both properties for all five systems.",
        "",
        "| family | method | a ME (A) | a MAE (A) | Ecoh ME (eV/atom) | Ecoh MAE (eV/atom) |",
        "|---|---|---:|---:|---:|---:|",
    ]
    for row in sorted(common5, key=lambda item: float(item["ecoh_MAE_eV_per_atom"])):
        lines.append(
            f"| {row['family']} | {row['method']} | {fmt(row['a_ME_A'])} | {fmt(row['a_MAE_A'])} | "
            f"{fmt(row['ecoh_ME_eV_per_atom'])} | {fmt(row['ecoh_MAE_eV_per_atom'])} |"
        )

    selected_methods = ("GFN1", "GFN2", "MP2", "SCS-MP2", "SCAN", "r2SCAN")
    by_key = {(str(row["method"]), str(row["solid"])): row for row in records}
    solids = ("C", "Si", "SiC", "BN", "BP", "AlN", "AlP", "MgS", "LiF", "LiCl")
    lines += [
        "",
        "## Per-system lattice constants",
        "",
        "Entries are value in A followed by the signed error in parentheses.",
        "",
        "| solid | " + " | ".join(selected_methods) + " |",
        "|---|" + "---:|" * len(selected_methods),
    ]
    for solid in solids:
        cells = [value_error(by_key.get((method, solid)), "a_calc_A", "a_error_A", 3) for method in selected_methods]
        lines.append(f"| {solid} | " + " | ".join(cells) + " |")

    lines += [
        "",
        "## Per-system cohesive energies",
        "",
        "Entries are value in eV/atom followed by the signed error in parentheses.",
        "",
        "| solid | " + " | ".join(selected_methods) + " |",
        "|---|" + "---:|" * len(selected_methods),
    ]
    for solid in solids:
        cells = [
            value_error(by_key.get((method, solid)), "ecoh_calc_eV_per_atom", "ecoh_error_eV_per_atom", 3)
            for method in selected_methods
        ]
        lines.append(f"| {solid} | " + " | ".join(cells) + " |")

    lines += [
        "",
        "## Sources",
        "",
        "- Goldzak et al., J. Chem. Phys. 157, 174112 (2022), https://doi.org/10.1063/5.0119633",
        "- Mejia-Rodriguez and Trickey, Phys. Rev. B 102, 121109(R) (2020), https://doi.org/10.1103/PhysRevB.102.121109",
        "- Mo et al., Phys. Rev. B 95, 035118 (2017), https://doi.org/10.1103/PhysRevB.95.035118",
    ]
    (DATA / "literature_comparison.md").write_text("\n".join(lines) + "\n")


def plot_available_mae(summary: list[dict[str, object]]) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(12.2, 8.4))
    specs = (
        ("a_MAE_A", "n_a", "Lattice-constant MAE (A)"),
        ("ecoh_MAE_eV_per_atom", "n_ecoh", "Cohesive-energy MAE (eV/atom)"),
    )
    for ax, (metric_key, n_key, xlabel) in zip(axes, specs):
        rows = [row for row in summary if not math.isnan(float(row[metric_key]))]
        rows.sort(key=lambda row: float(row[metric_key]), reverse=True)
        y = np.arange(len(rows))
        values = np.array([float(row[metric_key]) for row in rows])
        colors = [FAMILY_COLORS[str(row["family"])] for row in rows]
        ax.hlines(y, 0.0, values, color="#c7c7c7", linewidth=1.0)
        ax.scatter(values, y, c=colors, s=45, zorder=3)
        ax.set_yticks(y, [str(row["method"]) for row in rows])
        ax.set_xlabel(xlabel)
        ax.grid(axis="x", color="#d8d8d8", linewidth=0.6)
        ax.set_axisbelow(True)
        ax.set_xlim(left=0.0, right=max(values) * 1.19)
        for yi, (number, row) in enumerate(zip(values, rows)):
            ax.text(number + max(values) * 0.018, yi, f"n={row[n_key]}", va="center", fontsize=8)
    handles = [
        Line2D([0], [0], marker="o", linestyle="", color=color, label=family)
        for family, color in FAMILY_COLORS.items()
    ]
    fig.suptitle(
        "LC10: literature comparison against one common experimental reference",
        y=0.995,
    )
    fig.legend(
        handles=handles,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.965),
        ncol=3,
        frameon=False,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.925))
    FIGURES.mkdir(parents=True, exist_ok=True)
    out = FIGURES / "goldzak12_literature_mae_comparison"
    fig.savefig(out.with_suffix(".png"), dpi=220)
    fig.savefig(out.with_suffix(".pdf"))
    plt.close(fig)


def plot_common5_heatmap(records: list[dict[str, object]]) -> None:
    methods = (
        "GFN1",
        "GFN2",
        "HF",
        "MP2",
        "SCS-MP2",
        "PBE",
        "PBEsol",
        "HSE06",
        "optB86b-vdW",
        "SCAN",
        "SCAN-L",
        "r2SCAN",
        "r2SCAN-L",
    )
    by_key = {(str(row["method"]), str(row["solid"])): row for row in records}
    specs = (
        ("a_error_A", "Lattice-constant error (A)", 3),
        ("ecoh_error_eV_per_atom", "Cohesive-energy error (eV/atom)", 2),
    )
    fig, axes = plt.subplots(1, 2, figsize=(12.6, 8.1))
    for ax, (key, title, digits) in zip(axes, specs):
        matrix = np.array(
            [[float(by_key[(method, solid)][key]) for solid in COMMON5] for method in methods],
            dtype=float,
        )
        vmax = float(np.abs(matrix).max())
        image = ax.imshow(matrix, cmap="RdBu_r", vmin=-vmax, vmax=vmax, aspect="auto")
        ax.set_xticks(np.arange(len(COMMON5)), COMMON5)
        ax.set_yticks(np.arange(len(methods)), methods)
        ax.set_title(title)
        for i in range(matrix.shape[0]):
            for j in range(matrix.shape[1]):
                normalized = abs(matrix[i, j]) / vmax if vmax else 0.0
                color = "white" if normalized > 0.55 else "black"
                ax.text(j, i, f"{matrix[i, j]:+.{digits}f}", ha="center", va="center", color=color, fontsize=7.5)
        fig.colorbar(image, ax=ax, fraction=0.046, pad=0.04)
    fig.suptitle("LC10 common five-system subset: signed errors")
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    FIGURES.mkdir(parents=True, exist_ok=True)
    out = FIGURES / "goldzak12_common5_error_heatmap"
    fig.savefig(out.with_suffix(".png"), dpi=220)
    fig.savefig(out.with_suffix(".pdf"))
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mesh", help="final cohesive-energy mesh; defaults to build provenance")
    args = parser.parse_args()
    result_mesh = args.mesh or provenance_result_mesh()
    records, _ = collect_records(result_mesh)
    available = summarize(records)
    common5 = summarize(records, COMMON5)
    write_csv(DATA / "literature_comparison_records.csv", records)
    write_csv(DATA / "literature_comparison_summary.csv", available)
    write_csv(DATA / "literature_comparison_common5.csv", common5)
    write_markdown(records, available, common5)
    plot_available_mae(available)
    plot_common5_heatmap(records)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
