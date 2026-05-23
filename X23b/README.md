# X23b Periodic GFN Benchmark

This directory contains CP2K/tblite calculations for the X23b molecular-crystal
benchmark of Dolgonos, Hoja, and Boese. The reference lattice energies are the
recommended experimental back-corrected values from Table 5, and the reference
cell volumes are the electronic reference volumes from Table 2 of that work.

The crystal structures are taken from the open X23 `refdata` set. Hexamine is
the only special case: the open experimental CIF contains only heavy atoms, so
the complete X23 Quantum ESPRESSO crystal input is used for that system.

## Contents

- `structures/`: P1 CIF crystal structures and gas-phase molecular starting
  geometries.
- `inputs/`: CP2K input files for crystal single points, gas-phase molecular
  optimizations, and Gamma-point crystal cell optimizations.
- `runs/`: generated CP2K working directories, ignored by Git.
- `data/`: metadata, reference values, extracted energies, volume errors, and
  aggregate statistics, including the DMC-X23 comparison values used for the
  system-resolved lattice-energy figure.
- `figures/`: plots generated from the benchmark data, including the
  system-resolved X23b lattice-energy profile.
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
