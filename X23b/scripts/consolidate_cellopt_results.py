#!/usr/bin/env python3
"""Combine fully converged X23b cell-optimization result tables."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
METADATA = ROOT / "data" / "metadata.json"
METHODS = ("GFN1", "GFN2")
REQUIRED_NUMERIC_FIELDS = (
    "energy_hartree",
    "gas_energy_hartree",
    "lattice_energy_kJmol",
    "error_kJmol",
    "volume_A3",
    "volume_error_percent",
    "last_pressure_bar",
    "last_max_step",
    "last_rms_step",
    "last_max_gradient",
    "last_rms_gradient",
)


def read_rows(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            raise ValueError(f"missing CSV header: {path}")
        return list(reader.fieldnames), list(reader)


def is_finite_number(value: str) -> bool:
    try:
        number = float(value)
    except ValueError:
        return False
    return math.isfinite(number)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, action="append", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    systems = {str(row["id"]) for row in json.loads(METADATA.read_text())["systems"]}
    expected = {(method, system) for method in METHODS for system in systems}
    selected: dict[tuple[str, str], tuple[Path, dict[str, str]]] = {}
    fieldnames: list[str] | None = None

    for path in args.input:
        current_fields, rows = read_rows(path)
        if fieldnames is None:
            fieldnames = current_fields
        elif current_fields != fieldnames:
            raise ValueError(f"incompatible CSV columns: {path}")
        for row in rows:
            if row.get("opt_completed") != "True":
                continue
            if row.get("program_ended", "True") != "True":
                raise ValueError(f"optimization marked complete without a clean CP2K end in {path}")
            if row.get("mesh") != "k222":
                raise ValueError(f"unexpected cell-optimization mesh in {path}: {row.get('mesh')!r}")
            key = (row.get("method", ""), row.get("system", ""))
            if key not in expected:
                raise ValueError(f"unexpected method/system {key!r} in {path}")
            missing_values = [field for field in REQUIRED_NUMERIC_FIELDS if not is_finite_number(row.get(field, ""))]
            if missing_values:
                raise ValueError(f"non-numeric fields for {key!r} in {path}: {', '.join(missing_values)}")
            if key in selected:
                previous = selected[key][0]
                raise ValueError(f"duplicate converged result for {key!r}: {previous} and {path}")
            selected[key] = (path, row)

    missing = sorted(expected - set(selected))
    if missing:
        raise ValueError("missing converged results: " + ", ".join(f"{method}/{system}" for method, system in missing))
    assert fieldnames is not None

    rows = [selected[key][1] for key in sorted(selected)]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    for method in METHODS:
        count = sum(row["method"] == method for row in rows)
        print(f"{method}: {count}/{len(systems)} converged")


if __name__ == "__main__":
    main()
