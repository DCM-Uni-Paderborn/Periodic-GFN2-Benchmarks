#!/usr/bin/env python3
"""Normalize Klimes-23 result provenance from the saved CP2K inputs and outputs."""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
RAW = REPO / "Klimes-Solids23/results/raw"


def last_float(text: str, marker: str) -> float | None:
    value = None
    for line in text.splitlines():
        if marker not in line:
            continue
        try:
            value = float(line.split()[-1])
        except ValueError:
            continue
    return value


def keyword(text: str, name: str, cast, default=None):
    matches = re.findall(rf"^\s*{re.escape(name)}\s+([^\s]+)", text, flags=re.MULTILINE)
    if not matches:
        return default
    return cast(matches[-1])


def electronic_temperature(input_text: str) -> float | None:
    matches = re.findall(
        r"^\s*ELECTRONIC_TEMPERATURE(?:\s+\[[^]]+\])?\s+([^\s]+)",
        input_text,
        flags=re.MULTILINE,
    )
    return float(matches[-1]) if matches else None


def mixer_metadata(input_text: str) -> dict[str, object]:
    requested = keyword(input_text, "SCC_MIXER", str, "AUTO").upper()
    active = "TBLITE" if requested == "AUTO" else requested
    max_scf = keyword(input_text, "MAX_SCF", int)
    metadata: dict[str, object] = {
        "scc_mixer_requested": requested,
        "scc_mixer": active,
        "max_scf": max_scf,
    }
    if active == "TBLITE":
        explicit = "&TBLITE_MIXER" in input_text.upper()
        damping = keyword(input_text, "DAMPING", float, 0.4) if explicit else 0.4
        iterations = keyword(input_text, "ITERATIONS", int, 250) if explicit else 250
        memory = keyword(input_text, "MEMORY", int, iterations) if explicit else 250
        metadata.update(
            {
                "damping": damping,
                "tblite_mixer_damping": damping,
                "tblite_mixer_iterations": iterations,
                "tblite_mixer_memory": memory,
                "cp2k_mixer_alpha": None,
            }
        )
    elif active == "CP2K":
        alpha = keyword(input_text, "ALPHA", float)
        metadata.update(
            {
                "damping": alpha,
                "tblite_mixer_damping": None,
                "tblite_mixer_iterations": None,
                "tblite_mixer_memory": None,
                "cp2k_mixer_alpha": alpha,
            }
        )
    return metadata


def normalize(path: Path, dry_run: bool, require_output: bool = False) -> tuple[bool, float]:
    result = json.loads(path.read_text())
    if not result.get("ok"):
        return False, 0.0
    root = path.parent
    output_path = root / "cp2k.out"
    input_path = root / "cp2k.inp"
    if not input_path.is_file():
        raise RuntimeError(f"missing CP2K provenance next to {path}")
    if not output_path.is_file():
        normalized = (
            result.get("provenance_schema") == 2
            and result.get("energy_source") == "CP2K Total energy (extrapolated to T->0)"
            and result.get("free_energy_hartree") is not None
            and result.get("electronic_entropy_correction_hartree") is not None
        )
        if normalized and not require_output:
            return False, float(result["electronic_entropy_correction_hartree"])
        raise RuntimeError(f"missing CP2K output next to {path}")
    output = output_path.read_text(errors="replace")
    if "PROGRAM ENDED" not in output:
        raise RuntimeError(f"successful sidecar has incomplete CP2K output: {path}")
    free_energy = last_float(output, "ENERGY| Total FORCE_EVAL")
    zero_temperature = last_float(output, "Total energy (extrapolated to T->0)")
    if free_energy is None:
        raise RuntimeError(f"final CP2K energy not found: {output_path}")
    selected = zero_temperature if zero_temperature is not None else free_energy
    correction = selected - free_energy
    input_text = input_path.read_text(errors="replace")
    result.update(mixer_metadata(input_text))
    result.update(
        {
            "energy_hartree": selected,
            "free_energy_hartree": free_energy,
            "electronic_entropy_correction_hartree": correction,
            "energy_source": (
                "CP2K Total energy (extrapolated to T->0)"
                if zero_temperature is not None
                else "CP2K ENERGY| Total FORCE_EVAL"
            ),
            "smearing_temperature_K": electronic_temperature(input_text),
            "provenance_schema": 2,
        }
    )
    rendered = json.dumps(result, indent=2, sort_keys=True) + "\n"
    changed = rendered != path.read_text()
    if changed and not dry_run:
        path.write_text(rendered)
    return changed, correction


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--require-output", action="store_true")
    args = parser.parse_args()
    paths = sorted(RAW.glob("*/*/v_*/result.json"))
    changed = 0
    skipped = 0
    corrections: list[tuple[float, Path, float]] = []
    for path in paths:
        was_changed, correction = normalize(path, args.dry_run, args.require_output)
        if was_changed:
            changed += 1
        result = json.loads(path.read_text())
        if not result.get("ok"):
            skipped += 1
        corrections.append((abs(correction), path, correction))
    print(f"audited={len(paths)} changed={changed} failed_skipped={skipped} dry_run={args.dry_run}")
    for _, path, correction in sorted(corrections, reverse=True)[:10]:
        print(f"{path.relative_to(REPO)} dE(T->0-free)={correction:+.12e} Eh")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
