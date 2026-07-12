#!/usr/bin/env python3
"""Prepare and analyse CP2K/tblite DMC-ICE13 single-point benchmarks."""

from __future__ import annotations

import json
import math
import re
from pathlib import Path

HARTREE_TO_KJMOL = 2625.499638

PHASES = [
    "Ih",
    "II",
    "III",
    "IV",
    "VI",
    "VII",
    "VIII",
    "IX",
    "XI",
    "XIII",
    "XIV",
    "XV",
    "XVII",
]

DMC_ABS_KJMOL = {
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


def is_float_triplet(line: str) -> bool:
    parts = line.split()
    if len(parts) != 3:
        return False
    try:
        [float(x) for x in parts]
    except ValueError:
        return False
    return True


def read_triplet(line: str) -> list[float]:
    return [float(x) for x in line.split()]


def clean_lines(text: str) -> list[str]:
    return [line.strip().replace("\x0c", "").strip() for line in text.splitlines()]


def parse_poscars(si_text: Path) -> dict[str, dict]:
    lines = clean_lines(si_text.read_text(errors="ignore"))
    start = next(i for i, line in enumerate(lines) if line.startswith("9") and "Geometries" in line)
    lines = lines[start:]
    data: dict[str, dict] = {}
    for phase in PHASES:
        label = f"Ice {phase}"
        idx = next(i for i, line in enumerate(lines) if line == label)
        i = idx + 1
        while i < len(lines):
            try:
                scale = float(lines[i])
                break
            except ValueError:
                i += 1
        cell = []
        i += 1
        while len(cell) < 3:
            if is_float_triplet(lines[i]):
                cell.append(read_triplet(lines[i]))
            i += 1
        while i < len(lines) and lines[i].replace(" ", "") != "HO":
            i += 1
        i += 1
        while i < len(lines) and not re.match(r"^\d+\s+\d+$", lines[i]):
            i += 1
        n_h, n_o = [int(x) for x in lines[i].split()]
        i += 1
        while i < len(lines) and lines[i] not in {"Direct", "Cartesian"}:
            i += 1
        coord_mode = lines[i]
        i += 1
        coords: list[list[float]] = []
        while len(coords) < n_h + n_o:
            if is_float_triplet(lines[i]):
                coords.append(read_triplet(lines[i]))
            i += 1
        data[phase] = {
            "scale": scale,
            "cell": [[scale * x for x in vec] for vec in cell],
            "mode": coord_mode,
            "counts": {"H": n_h, "O": n_o},
            "coords": coords,
        }
    return data


def cp2k_input(phase: str, geom: dict, method: str) -> str:
    project = f"ice_{phase}_{method}".replace("-", "_")
    lines = [
        "&GLOBAL",
        "  PRINT_LEVEL LOW",
        f"  PROJECT {project}",
        "  RUN_TYPE ENERGY",
        "&END GLOBAL",
        "",
        "&FORCE_EVAL",
        "  METHOD Quickstep",
        "  &DFT",
        "    &QS",
        "      EPS_DEFAULT 1.0E-12",
        "      METHOD xTB",
        "      &XTB",
        "        GFN_TYPE TBLITE",
        "        &TBLITE",
        f"          METHOD {method}",
        "          ACCURACY 0.1",
        "        &END TBLITE",
        "      &END XTB",
        "    &END QS",
        "    &SCF",
        "      EPS_SCF 1.0E-9",
        "      MAX_SCF 300",
        "      SCF_GUESS MOPAC",
        "      &MIXING",
        "        METHOD DIRECT_P_MIXING",
        "        ALPHA 0.2",
        "      &END MIXING",
        "      &PRINT",
        "        &RESTART OFF",
        "        &END RESTART",
        "      &END PRINT",
        "    &END SCF",
        "  &END DFT",
        "  &SUBSYS",
        "    &CELL",
        "      PERIODIC XYZ",
    ]
    for key, vec in zip(("A", "B", "C"), geom["cell"]):
        lines.append(f"      {key} {vec[0]:.12f} {vec[1]:.12f} {vec[2]:.12f}")
    lines += [
        "    &END CELL",
        "    &COORD",
    ]
    if geom["mode"] == "Direct":
        lines.append("      SCALED")
    n_h = geom["counts"]["H"]
    coords = geom["coords"]
    for element, coord in [("H", xyz) for xyz in coords[:n_h]] + [("O", xyz) for xyz in coords[n_h:]]:
        lines.append(f"      {element:<2} {coord[0]: .12f} {coord[1]: .12f} {coord[2]: .12f}")
    lines += [
        "    &END COORD",
        "  &END SUBSYS",
        "&END FORCE_EVAL",
        "",
    ]
    return "\n".join(lines)


def write_inputs(data: dict[str, dict], root: Path) -> None:
    (root / "poscars").mkdir(parents=True, exist_ok=True)
    (root / "inputs").mkdir(parents=True, exist_ok=True)
    for phase, geom in data.items():
        poscar = [f"Ice {phase}", f"{geom['scale']:.16g}"]
        poscar += [" ".join(f"{x:.12f}" for x in vec) for vec in geom["cell"]]
        poscar += ["H O", f"{geom['counts']['H']} {geom['counts']['O']}", geom["mode"]]
        poscar += [" ".join(f"{x:.12f}" for x in coord) for coord in geom["coords"]]
        (root / "poscars" / f"ice_{phase}.vasp").write_text("\n".join(poscar) + "\n")
        for method in ("GFN1", "GFN2"):
            (root / "inputs" / f"ice_{phase}_{method}.inp").write_text(cp2k_input(phase, geom, method))


def parse_energy(output: Path) -> float | None:
    energy = None
    if not output.exists():
        return None
    for line in output.read_text(errors="ignore").splitlines():
        if "ENERGY| Total FORCE_EVAL" in line:
            energy = float(line.split()[-1])
    return energy


def analyse(root: Path) -> dict:
    geometry_path = root / "data" / "geometries.json"
    if not geometry_path.exists():
        geometry_path = root / "geometries.json"
    geometries = json.loads(geometry_path.read_text())
    results = {}
    for method in ("GFN1", "GFN2"):
        values = {}
        for phase in PHASES:
            out = root / "runs" / method / phase / f"ice_{phase}_{method}.out"
            energy = parse_energy(out)
            values[phase] = energy
        if any(v is None for v in values.values()):
            results[method] = {"energies_hartree": values, "complete": False}
            continue
        per_h2o = {}
        for phase in PHASES:
            geom = geometries[phase]
            n_h2o = geom["counts"]["O"]
            per_h2o[phase] = values[phase] / n_h2o
        ih = per_h2o["Ih"]
        rel = {phase: (per_h2o[phase] - ih) * HARTREE_TO_KJMOL for phase in PHASES}
        errors = {phase: rel[phase] - (DMC_ABS_KJMOL[phase] - DMC_ABS_KJMOL["Ih"]) for phase in PHASES}
        benchmark_phases = [phase for phase in PHASES if phase != "Ih"]
        abs_errors = [abs(errors[phase]) for phase in benchmark_phases]
        signed_errors = [errors[phase] for phase in benchmark_phases]
        results[method] = {
            "complete": True,
            "energies_hartree": values,
            "per_h2o_hartree": per_h2o,
            "relative_kjmol": rel,
            "relative_errors_kjmol": errors,
            "mae_relative_kjmol": sum(abs_errors) / len(abs_errors),
            "me_relative_kjmol": sum(signed_errors) / len(signed_errors),
            "rmse_relative_kjmol": math.sqrt(sum(e * e for e in signed_errors) / len(signed_errors)),
            "maxae_relative_kjmol": max(abs_errors),
        }
    return results


def main() -> None:
    root = Path(__file__).resolve().parents[1]
    data_dir = root / "data"
    data_dir.mkdir(exist_ok=True)
    si_candidates = [root / "SI_ucl.txt", Path(__file__).resolve().parent / "SI_ucl.txt"]
    si_text = next((path for path in si_candidates if path.exists()), None)
    if si_text is not None:
        data = parse_poscars(si_text)
        (data_dir / "geometries.json").write_text(json.dumps(data, indent=2))
        write_inputs(data, root)
    elif not (data_dir / "geometries.json").exists():
        raise FileNotFoundError("Need SI_ucl.txt or data/geometries.json to analyse DMC-ICE13.")
    results = analyse(root)
    (data_dir / "results.json").write_text(json.dumps(results, indent=2))
    print(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()
