#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

import run_goldzak12_benchmark as base


ROOT = base.ROOT
DEFAULT_BASELINE = ROOT / "data" / "baseline_20260710" / "eos_results.csv"


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def number(row: dict[str, str] | None, key: str) -> float | None:
    if row is None or row.get(key, "") == "":
        return None
    return float(row[key])


def difference(new: float | None, old: float | None) -> float | None:
    if new is None or old is None:
        return None
    return new - old


def formatted(value: float | None, digits: int = 8) -> str:
    return "" if value is None else f"{value:.{digits}f}"


def selected(rows: list[dict[str, str]], mesh: str) -> dict[tuple[str, str], dict[str, str]]:
    return {(row["solid"], row["method"]): row for row in rows if row["energy_mesh"] == mesh}


def comparison_records(
    old_rows: list[dict[str, str]],
    new_rows: list[dict[str, str]],
    mesh: str,
) -> list[dict[str, object]]:
    old = selected(old_rows, mesh)
    new = selected(new_rows, mesh)
    records: list[dict[str, object]] = []
    for ref in base.REFERENCES:
        for method in base.METHODS:
            old_row = old.get((ref.solid, method))
            new_row = new.get((ref.solid, method))
            old_a = number(old_row, "a_calc_A")
            new_a = number(new_row, "a_calc_A")
            old_e = number(old_row, "ecoh_calc_eV_per_atom")
            new_e = number(new_row, "ecoh_calc_eV_per_atom")
            records.append(
                {
                    "solid": ref.solid,
                    "method": method,
                    "energy_mesh": mesh,
                    "a_ref_A": f"{ref.a_exp:.8f}",
                    "a_old_A": formatted(old_a),
                    "a_new_A": formatted(new_a),
                    "a_delta_new_minus_old_A": formatted(difference(new_a, old_a)),
                    "a_error_old_A": formatted(difference(old_a, ref.a_exp)),
                    "a_error_new_A": formatted(difference(new_a, ref.a_exp)),
                    "ecoh_ref_eV_per_atom": f"{ref.ecoh_exp:.8f}",
                    "ecoh_old_eV_per_atom": formatted(old_e),
                    "ecoh_new_eV_per_atom": formatted(new_e),
                    "ecoh_delta_new_minus_old_eV_per_atom": formatted(difference(new_e, old_e)),
                    "ecoh_error_old_eV_per_atom": formatted(difference(old_e, ref.ecoh_exp)),
                    "ecoh_error_new_eV_per_atom": formatted(difference(new_e, ref.ecoh_exp)),
                }
            )
    return records


def mae(values: list[float]) -> float:
    return sum(abs(value) for value in values) / len(values)


def me(values: list[float]) -> float:
    return sum(values) / len(values)


def rmse(values: list[float]) -> float:
    return math.sqrt(sum(value * value for value in values) / len(values))


def metric_summary(
    records: list[dict[str, object]],
    method: str,
    metric: str,
) -> dict[str, object]:
    if metric == "a":
        label = "lattice_constant_A"
        old_key, new_key = "a_old_A", "a_new_A"
        old_error_key, new_error_key = "a_error_old_A", "a_error_new_A"
    else:
        label = "cohesive_energy_eV_per_atom"
        old_key, new_key = "ecoh_old_eV_per_atom", "ecoh_new_eV_per_atom"
        old_error_key, new_error_key = "ecoh_error_old_eV_per_atom", "ecoh_error_new_eV_per_atom"

    rows = [row for row in records if row["method"] == method]
    old_errors = [float(row[old_error_key]) for row in rows if row[old_error_key] != ""]
    new_errors = [float(row[new_error_key]) for row in rows if row[new_error_key] != ""]
    common = [row for row in rows if row[old_key] != "" and row[new_key] != ""]
    old_common_errors = [float(row[old_error_key]) for row in common]
    new_common_errors = [float(row[new_error_key]) for row in common]
    changes = [float(row[new_key]) - float(row[old_key]) for row in common]

    return {
        "method": method,
        "metric": label,
        "n_old": len(old_errors),
        "n_new": len(new_errors),
        "n_common": len(common),
        "old_ME_reported": formatted(me(old_errors)),
        "new_ME_reported": formatted(me(new_errors)),
        "old_MAE_reported": formatted(mae(old_errors)),
        "new_MAE_reported": formatted(mae(new_errors)),
        "delta_MAE_reported": formatted(mae(new_errors) - mae(old_errors)),
        "old_RMSE_reported": formatted(rmse(old_errors)),
        "new_RMSE_reported": formatted(rmse(new_errors)),
        "old_MAE_common": formatted(mae(old_common_errors)),
        "new_MAE_common": formatted(mae(new_common_errors)),
        "delta_MAE_common": formatted(mae(new_common_errors) - mae(old_common_errors)),
        "mean_abs_value_change": formatted(mae(changes)),
        "max_abs_value_change": formatted(max(abs(value) for value in changes)),
    }


def make_summary(records: list[dict[str, object]]) -> list[dict[str, object]]:
    return [metric_summary(records, method, metric) for method in base.METHODS for metric in ("a", "ecoh")]


def write_markdown(records: list[dict[str, object]], summary: list[dict[str, object]], mesh: str) -> None:
    lines = [
        f"# LC12 old/new comparison ({mesh})",
        "",
        "Negative MAE changes indicate an improvement. The common-subset columns isolate numerical changes from coverage changes.",
        "",
        "## Summary",
        "",
        "| method | metric | n old/new/common | MAE old | MAE new | dMAE | dMAE common | mean abs value change | max abs value change |",
        "|---|---|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in summary:
        lines.append(
            f"| {row['method']} | {row['metric']} | {row['n_old']}/{row['n_new']}/{row['n_common']} | "
            f"{row['old_MAE_reported']} | {row['new_MAE_reported']} | {row['delta_MAE_reported']} | "
            f"{row['delta_MAE_common']} | {row['mean_abs_value_change']} | {row['max_abs_value_change']} |"
        )

    lines += [
        "",
        "## Per-system changes",
        "",
        "| solid | method | a old | a new | da new-old | Ecoh old | Ecoh new | dEcoh new-old |",
        "|---|---|---:|---:|---:|---:|---:|---:|",
    ]
    for row in records:
        lines.append(
            f"| {row['solid']} | {row['method']} | {base.fmt(row['a_old_A'], 4)} | {base.fmt(row['a_new_A'], 4)} | "
            f"{base.fmt(row['a_delta_new_minus_old_A'], 4)} | {base.fmt(row['ecoh_old_eV_per_atom'], 4)} | "
            f"{base.fmt(row['ecoh_new_eV_per_atom'], 4)} | {base.fmt(row['ecoh_delta_new_minus_old_eV_per_atom'], 4)} |"
        )
    (ROOT / "data" / "old_vs_new.md").write_text("\n".join(lines) + "\n")


def plot(records: list[dict[str, object]], mesh: str) -> None:
    solids = [ref.solid for ref in base.REFERENCES]
    x = np.arange(len(solids))
    width = 0.36
    colors = {"GFN1": "#3B6FB6", "GFN2": "#D97706"}
    fig, axes = plt.subplots(2, 1, figsize=(10.6, 7.2), sharex=True)
    panels = (
        ("a_delta_new_minus_old_A", "new - old lattice constant (A)"),
        ("ecoh_delta_new_minus_old_eV_per_atom", "new - old cohesive energy (eV/atom)"),
    )
    for ax, (key, ylabel) in zip(axes, panels):
        for index, method in enumerate(base.METHODS):
            values = []
            for solid in solids:
                row = next(item for item in records if item["solid"] == solid and item["method"] == method)
                values.append(float(row[key]) if row[key] != "" else np.nan)
            positions = x + (index - 0.5) * width
            ax.bar(positions, values, width, label=method, color=colors[method])
            for position, value in zip(positions, values):
                if np.isnan(value):
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
        ax.axhline(0.0, color="#222222", linewidth=0.8)
        ax.set_ylabel(ylabel)
        ax.grid(axis="y", color="#d0d0d0", linewidth=0.6, alpha=0.7)
    axes[0].legend(frameon=False, ncol=2)
    axes[-1].set_xticks(x)
    axes[-1].set_xticklabels(solids, rotation=45, ha="right")
    fig.suptitle(f"LC12 current stack versus previous results ({mesh} native Bloch)")
    fig.tight_layout()
    output = ROOT / "figures" / "goldzak12_old_vs_new_deltas"
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output.with_suffix(".png"), dpi=220)
    fig.savefig(output.with_suffix(".pdf"))
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--baseline", type=Path, default=DEFAULT_BASELINE)
    parser.add_argument("--current", type=Path, default=ROOT / "data" / "eos_results.csv")
    parser.add_argument("--mesh", default="k444")
    args = parser.parse_args()

    records = comparison_records(read_csv(args.baseline), read_csv(args.current), args.mesh)
    summary = make_summary(records)
    base.write_csv(ROOT / "data" / "old_vs_new_records.csv", records)
    base.write_csv(ROOT / "data" / "old_vs_new_summary.csv", summary)
    write_markdown(records, summary, args.mesh)
    plot(records, args.mesh)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
