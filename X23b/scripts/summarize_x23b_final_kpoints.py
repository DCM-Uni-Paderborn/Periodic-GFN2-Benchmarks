#!/usr/bin/env python3
"""Summarize k-point single points on final X23b k222 cell-opt geometries."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

import x23b_kpoint_cellopt as cellopt


ROW_FIELDS = (
    "method",
    "system",
    "k222_lattice_energy_kJmol",
    "k333_lattice_energy_kJmol",
    "k444_lattice_energy_kJmol",
    "k222_error_kJmol",
    "k333_error_kJmol",
    "k444_error_kJmol",
    "delta_k333_minus_k222_kJmol",
    "delta_k444_minus_k333_kJmol",
)

SUMMARY_FIELDS = (
    "method",
    "mesh",
    "N",
    "ME",
    "MAE",
    "RMSE",
    "MaxAE",
    "mean_abs_change_from_previous_kJmol",
    "max_abs_change_from_previous_kJmol",
)


def read_rows(path: Path) -> dict[tuple[str, str], dict[str, str]]:
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    selected = {(row["method"], row["system"]): row for row in rows}
    if len(selected) != len(rows):
        raise ValueError(f"duplicate method/system rows in {path}")
    return selected


def number(row: dict[str, str], field: str) -> float:
    value = float(row[field])
    if not math.isfinite(value):
        raise ValueError(f"non-finite {field} for {row['method']}/{row['system']}")
    return value


def fmt(value: float) -> str:
    return f"{value:.12f}"


def write_csv(path: Path, fields: tuple[str, ...], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--k333-csv", type=Path, required=True)
    parser.add_argument("--k444-csv", type=Path, required=True)
    parser.add_argument("--rows-csv", type=Path, required=True)
    parser.add_argument("--summary-csv", type=Path, required=True)
    parser.add_argument("--allow-incomplete", action="store_true")
    args = parser.parse_args()

    k333 = read_rows(args.k333_csv)
    k444 = read_rows(args.k444_csv)
    if set(k333) != set(k444):
        missing_333 = sorted(set(k444) - set(k333))
        missing_444 = sorted(set(k333) - set(k444))
        raise ValueError(f"mismatched result sets: missing k333={missing_333}, missing k444={missing_444}")

    expected = {
        (method, str(system["id"]))
        for method in cellopt.METHODS
        for system in cellopt.systems()
    }
    if not args.allow_incomplete and set(k333) != expected:
        missing = sorted(expected - set(k333))
        extra = sorted(set(k333) - expected)
        raise ValueError(f"complete 46-case result required: missing={missing}, extra={extra}")

    rows: list[dict[str, str]] = []
    values: dict[tuple[str, str], list[float]] = {}
    changes: dict[tuple[str, str], list[float]] = {}
    for key in sorted(k333):
        row3, row4 = k333[key], k444[key]
        if row3["program_ended"] != "True" or row4["program_ended"] != "True":
            raise ValueError(f"incomplete single point for {key[0]}/{key[1]}")
        source3 = number(row3, "source_lattice_energy_kJmol")
        source4 = number(row4, "source_lattice_energy_kJmol")
        if abs(source3 - source4) > 1.0e-8:
            raise ValueError(f"inconsistent k222 source energy for {key[0]}/{key[1]}")
        target3 = number(row3, "target_lattice_energy_kJmol")
        target4 = number(row4, "target_lattice_energy_kJmol")
        error3 = number(row3, "target_error_kJmol")
        error4 = number(row4, "target_error_kJmol")
        delta3 = number(row3, "delta_target_minus_source_kJmol")
        error2 = error3 - delta3
        delta4_from3 = target4 - target3
        rows.append(
            {
                "method": key[0],
                "system": key[1],
                "k222_lattice_energy_kJmol": fmt(source3),
                "k333_lattice_energy_kJmol": fmt(target3),
                "k444_lattice_energy_kJmol": fmt(target4),
                "k222_error_kJmol": fmt(error2),
                "k333_error_kJmol": fmt(error3),
                "k444_error_kJmol": fmt(error4),
                "delta_k333_minus_k222_kJmol": fmt(delta3),
                "delta_k444_minus_k333_kJmol": fmt(delta4_from3),
            }
        )
        values.setdefault((key[0], "k222"), []).append(error2)
        values.setdefault((key[0], "k333"), []).append(error3)
        values.setdefault((key[0], "k444"), []).append(error4)
        changes.setdefault((key[0], "k333"), []).append(delta3)
        changes.setdefault((key[0], "k444"), []).append(delta4_from3)

    summary: list[dict[str, str]] = []
    for method in cellopt.METHODS:
        for mesh in ("k222", "k333", "k444"):
            errors = values.get((method, mesh), [])
            if not errors:
                continue
            delta = changes.get((method, mesh), [])
            summary.append(
                {
                    "method": method,
                    "mesh": mesh,
                    "N": str(len(errors)),
                    "ME": fmt(sum(errors) / len(errors)),
                    "MAE": fmt(sum(abs(value) for value in errors) / len(errors)),
                    "RMSE": fmt(math.sqrt(sum(value * value for value in errors) / len(errors))),
                    "MaxAE": fmt(max(abs(value) for value in errors)),
                    "mean_abs_change_from_previous_kJmol": (
                        fmt(sum(abs(value) for value in delta) / len(delta)) if delta else ""
                    ),
                    "max_abs_change_from_previous_kJmol": (
                        fmt(max(abs(value) for value in delta)) if delta else ""
                    ),
                }
            )

    write_csv(args.rows_csv, ROW_FIELDS, rows)
    write_csv(args.summary_csv, SUMMARY_FIELDS, summary)
    for row in summary:
        print(
            f"{row['method']} {row['mesh']} N={row['N']} "
            f"MAE={float(row['MAE']):.6f} MaxAE={float(row['MaxAE']):.6f}"
        )


if __name__ == "__main__":
    main()
