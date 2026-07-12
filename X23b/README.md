# X23b Periodic GFN Benchmark

This directory contains CP2K/tblite calculations for the X23b molecular-crystal
benchmark of Dolgonos, Hoja, and Boese. The reference lattice energies are the
recommended experimental back-corrected values from Table 5, and the reference
cell volumes are the electronic reference volumes from Table 2 of that work.
The primary relaxed-cell benchmark in the current manuscript revision uses
native Bloch 2x2x2 CP2K `&KPOINTS` cell optimizations with full SPGLIB
symmetry reduction. Reported lattice energies are reevaluated on those final
geometries with a Gamma-centered 3x3x3 mesh. This is neither a Gamma-only cell
optimization nor a Born-von-Karman supercell calculation.

The crystal structures are taken from the open X23 `refdata` set. Hexamine is
the only special case: the open experimental CIF contains only heavy atoms, so
the complete X23 Quantum ESPRESSO crystal input is used for that system.

## Contents

- `structures/`: P1 CIF crystal structures and gas-phase molecular starting
  geometries.
- `inputs/`: CP2K input files for crystal single points, gas-phase molecular
  optimizations, and retained Gamma-point crystal cell optimizations.
- `runs/`: generated CP2K working directories, ignored by Git.
- `data/`: metadata, reference values, extracted energies, volume errors, and
  aggregate statistics, including the DMC-X23 comparison values used for the
  system-resolved lattice-energy figure. The final-geometry mesh convergence
  is retained in `x23b_final_geometry_kpoint_{rows,summary}.csv`. The
  exact source, patch, executable, and protocol record is
  `data/build_provenance.json`; the
  `x23b_reference_cli_{gfn1,gfn2}_{rows,summary}.csv` files record direct
  CP2K-native versus tblite CLI checks on all initial and final Gamma
  geometries.
- `figures/`: PDF, SVG, and PNG versions of the three X23b plots used in the
  revised manuscript and Supporting Information.
- `scripts/`: input generation, analysis, plotting, and run scripts.

## Run Defaults

The run script expects the CP2K executable through the `CP2K` environment
variable, or otherwise falls back to `cp2k.psmp`. The default execution mode is
many independent single-core jobs:

- `OMP_NUM_THREADS=1`
- `OPENBLAS_NUM_THREADS=1`
- `MKL_NUM_THREADS=1`
- `CP2K_PARALLEL_JOBS=20`

This was faster for the small DMC-ICE13 and X23b-style xTB jobs than hybrid
MPI/OpenMP execution.

## Current primary result

The final 23/23 converged X23b relaxed-cell data per method are stored in
`data/x23b_lattice_energies.csv`, `data/x23b_cell_volumes.csv`, and
`data/x23b_summary.csv`. The volume rows are `cell_opt,k222`; the manuscript
lattice-energy rows are `cell_opt_single_point,k333`. The raw k222 energies
and the k444 convergence checks remain in the same files.

| Quantity | Method | ME | MAE | RMSE | MaxAE |
|---|---|---:|---:|---:|---:|
| Lattice energy / kJ mol-1, k333 on k222 geometry | GFN1-xTB | 0.258871 | 11.345702 | 14.019344 | 30.935058 |
| Lattice energy / kJ mol-1, k333 on k222 geometry | GFN2-xTB | -12.018989 | 14.092104 | 21.341752 | 77.785392 |
| Cell volume / percent, k222 optimization | GFN1-xTB | -5.960071 | 7.514116 | 9.019708 | 19.236681 |
| Cell volume / percent, k222 optimization | GFN2-xTB | -1.657324 | 5.842296 | 7.530373 | 19.952589 |

The k333-to-k444 mean absolute energy changes on the final geometries are
0.079329 kJ mol-1 for GFN1-xTB and 0.084265 kJ mol-1 for GFN2-xTB. The
fixed-reference-geometry single-point rows remain as a separate diagnostic.

Recreate the curated tables and figures after collecting the cell
optimizations and final-geometry single points:

```bash
python3 X23b/scripts/x23b_pipeline.py analyse \
  --cellopt-csv /path/to/x23b_k222_cellopt_results.csv \
  --final-kpoint-csv X23b/data/x23b_final_geometry_kpoint_rows.csv
```
