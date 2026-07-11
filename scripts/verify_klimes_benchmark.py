#!/usr/bin/env python3
"""Verify the complete Klimes-23 result, provenance, table, and figure set."""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
ROOT = REPO / "Klimes-Solids23"
RESULTS = ROOT / "results"
RAW = RESULTS / "raw"
METHODS = ("GFN1", "GFN2")
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
METAL_CATEGORIES = {"transition_metal", "alkali_metal", "alkaline_earth_metal", "simple_metal"}


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> int:
    lattice = read_csv(ROOT / "references/lattice_constants.csv")
    bulk = read_csv(ROOT / "references/bulk_moduli.csv")
    cohesive = read_csv(ROOT / "references/cohesive_energies.csv")
    systems = [row["system"] for row in lattice]
    categories = {row["system"]: row["category"] for row in lattice}
    require(len(systems) == len(set(systems)) == 23, "reference set is not 23 unique systems")
    require([row["system"] for row in bulk] == systems, "bulk reference order differs")
    require([row["system"] for row in cohesive] == systems, "cohesive reference order differs")
    for rows in (lattice, bulk, cohesive):
        require(all(method in rows[0] for method in LITERATURE), "published DFT column missing")

    sidecars = sorted(RAW.glob("*/*/v_*/result.json"))
    inputs = sorted(RAW.glob("*/*/v_*/cp2k.inp"))
    require(len(sidecars) == len(inputs) == 867, "raw input/sidecar count mismatch")
    successful = 0
    failed: list[tuple[str, str, float]] = []
    for path in sidecars:
        result = json.loads(path.read_text())
        input_path = path.with_name("cp2k.inp")
        require(input_path.is_file(), f"missing input next to {path}")
        text = input_path.read_text()
        system = str(result["system"])
        expected = 16 if categories[system] in METAL_CATEGORIES else 8
        require(f"SCHEME MACDONALD {expected} {expected} {expected}" in text, f"wrong k-mesh: {path}")
        require("FULL_GRID T" in text, f"reduced grid in production input: {path}")
        require("MULTIPLE_UNIT_CELL" not in text and "SUPERCELL" not in text, f"BvK keyword found: {path}")
        if result.get("ok"):
            successful += 1
            require(result.get("energy_source") == "CP2K Total energy (extrapolated to T->0)", f"wrong energy source: {path}")
            require(result.get("kpoint_path") == "CP2K native Bloch", f"wrong k-point path: {path}")
            require(result.get("full_grid") is True, f"sidecar is not full-grid: {path}")
            require(result.get("provenance_schema") == 2, f"old provenance schema: {path}")
            require(float(result.get("smearing_temperature_K")) == 300.0, f"wrong smearing: {path}")
        else:
            failed.append((str(result["method"]), system, float(result["volume_factor"])))
    require(successful == 866, f"expected 866 successful sidecars, found {successful}")
    require(failed == [("GFN2", "GaAs", 0.8)], f"unexpected failed raw jobs: {failed}")

    points = read_csv(RESULTS / "eos_points.csv")
    fits = read_csv(RESULTS / "eos_fits.csv")
    comparison = read_csv(RESULTS / "comparison.csv")
    long_rows = read_csv(RESULTS / "literature_comparison_long.csv")
    system_summary = read_csv(RESULTS / "system_summary.csv")
    paired = read_csv(RESULTS / "paired_gfn_statistics.csv")
    require(len(points) == 866, "EOS point table does not contain all successful sidecars")
    require(len(fits) == 46, "EOS fit matrix is not 23 x 2")
    require(all(row.get("eos_model") == "Murnaghan" for row in fits), "non-Murnaghan fit found")
    accepted = [row for row in fits if row["fit_ok"] == "True"]
    rejected = {(row["method"], row["system"]) for row in fits if row["fit_ok"] != "True"}
    require(len(accepted) == 44, "expected 44 accepted EOS fits")
    require(rejected == {("GFN2", "Cu"), ("GFN2", "MgO")}, f"unexpected EOS failures: {rejected}")
    require(len(comparison) == 44, "comparison table is not aligned with accepted fits")
    require(len(long_rows) == 684, "long-form literature table has unexpected coverage")
    require(len(system_summary) == 23, "system summary does not contain 23 systems")
    require(len(paired) == 3 and all(int(row["n_paired"]) == 21 for row in paired), "paired statistics are not n=21")

    for stem in ("klimes23-system-relative-errors", "klimes23-aggregate-mare", "klimes23-gfn2-change"):
        for suffix in ("png", "pdf"):
            path = ROOT / "figures" / f"{stem}.{suffix}"
            require(path.is_file() and path.stat().st_size > 10_000, f"missing or empty figure: {path}")

    provenance = json.loads((ROOT / "provenance.json").read_text())
    require(provenance["counts"]["raw_sidecars"] == len(sidecars), "manifest sidecar count differs")
    require(provenance["counts"]["accepted_eos_fits"] == len(accepted), "manifest fit count differs")
    cp2k = Path(provenance["cp2k"]["binary_path_as_run"])
    tblite = Path(provenance["tblite"]["binary_path_as_run"])
    patch = (ROOT / provenance["cp2k"]["local_patch"]).resolve()
    if cp2k.is_file():
        require(sha256(cp2k) == provenance["cp2k"]["binary_sha256"], "CP2K binary hash differs")
    if tblite.is_file():
        require(sha256(tblite) == provenance["tblite"]["binary_sha256"], "tblite binary hash differs")
    require(sha256(patch) == provenance["cp2k"]["local_patch_sha256"], "CP2K patch hash differs")

    print(
        "Klimes-23 verification passed: "
        f"{successful}/{len(sidecars)} raw jobs successful, {len(accepted)}/{len(fits)} EOS fits accepted, "
        "23 systems, 8 published DFT comparators, native Bloch only"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
