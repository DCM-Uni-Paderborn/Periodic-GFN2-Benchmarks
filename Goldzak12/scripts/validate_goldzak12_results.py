#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import run_goldzak12_benchmark as base
import run_goldzak12_eos_benchmark as eos


ROOT = base.ROOT


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def truth(value: str) -> bool:
    return value.strip().lower() == "true"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--eos-mesh", default="k444")
    parser.add_argument("--energy-mesh", action="append", default=[])
    parser.add_argument("--result-mesh", default="k444")
    args = parser.parse_args()
    energy_meshes = args.energy_mesh or ["k333", "k444", "k555"]

    expected_pairs = {(ref.solid, method) for ref in base.REFERENCES for method in base.METHODS}
    problems: list[str] = []

    points = read_csv(ROOT / "data" / "eos_points.csv")
    expected_point_count = sum(
        len(eos.scales_for(solid, method, eos.DEFAULT_SCALES)) for solid, method in expected_pairs
    )
    if len(points) != expected_point_count:
        problems.append(f"EOS point rows: expected {expected_point_count}, found {len(points)}")
    failed_points = [row for row in points if not truth(row["completed"])]

    fits = read_csv(ROOT / "data" / "eos_fits.csv")
    fit_pairs = {(row["solid"], row["method"]) for row in fits if row["eos_mesh"] == args.eos_mesh}
    if fit_pairs != expected_pairs:
        problems.append(f"EOS fit coverage differs: missing {sorted(expected_pairs - fit_pairs)}")
    bad_fits = [row for row in fits if row["a_eos_A"] == "" or row["fit_status"] != "quadratic"]
    allowed_bad_fits = [
        row for row in bad_fits if row["fit_status"] in {"poor_quadratic_fit", "no_local_minimum"}
    ]
    unexpected_bad_fits = [row for row in bad_fits if row not in allowed_bad_fits]
    if unexpected_bad_fits:
        labels = ", ".join(f"{row['method']}/{row['solid']}={row['fit_status']}" for row in unexpected_bad_fits)
        problems.append(f"Unexpected invalid EOS fits ({len(unexpected_bad_fits)}): {labels}")

    results = read_csv(ROOT / "data" / "eos_results.csv")
    valid_pairs = {(row["solid"], row["method"]) for row in fits if row["a_eos_A"] != ""}
    expected_results = {(solid, method, mesh) for solid, method in valid_pairs for mesh in energy_meshes}
    result_keys = {(row["solid"], row["method"], row["energy_mesh"]) for row in results}
    if result_keys != expected_results:
        problems.append(f"Final result coverage differs: missing {sorted(expected_results - result_keys)}")
    failed_sp = [row for row in results if not truth(row["sp_completed"])]
    if failed_sp:
        labels = ", ".join(f"{row['method']}/{row['solid']}/{row['energy_mesh']}" for row in failed_sp)
        problems.append(f"Incomplete final single points ({len(failed_sp)}): {labels}")

    summary = read_csv(ROOT / "data" / "eos_summary.csv")
    gfn_summary = [row for row in summary if row["source"] == "CP2K/tblite EOS"]
    expected_method_counts = {
        method: sum(1 for solid, fit_method in valid_pairs if fit_method == method) for method in base.METHODS
    }
    for row in gfn_summary:
        if int(row["n_complete"]) != expected_method_counts[row["method"]]:
            problems.append(
                f"Summary coverage {row['method']}: {row['n_complete']}/{expected_method_counts[row['method']]}"
            )

    provenance_path = ROOT / "data" / "build_provenance.json"
    provenance = json.loads(provenance_path.read_text())
    if provenance["protocol"]["result_mesh"] != args.result_mesh:
        problems.append("Provenance result mesh does not match validation request")

    if problems:
        print("LC12 validation FAILED")
        for problem in problems:
            print(f"- {problem}")
        return 1

    print(
        f"LC12 validation passed: {len(points)} EOS points ({len(failed_points)} nonessential failures), "
        f"{len(valid_pairs)}/{len(fits)} valid fits, {len(results)} final single points."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
