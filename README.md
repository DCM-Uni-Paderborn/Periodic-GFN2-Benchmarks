# Periodic GFN2 Benchmarks

This repository collects the paper-relevant benchmark inputs, curated output
data, analysis scripts, and the figures used in the revised manuscript and
Supporting Information for periodic GFN calculations in CP2K.

## Contents

- `DMC-ICE13/`: CP2K/tblite single-point benchmark for the DMC-ICE13 ice
  polymorph data set, comparing periodic GFN1-xTB and GFN2-xTB relative
  energies with the diffusion Monte Carlo reference energies. The current
  manuscript data use native Bloch k-point calculations through CP2K
  `&KPOINTS`.
- `X23b/`: CP2K/tblite molecular-crystal benchmark with gas-phase molecular
  optimizations, crystal single-point k-point tests, native Bloch 2x2x2
  crystal cell optimizations, converged 3x3x3/4x4x4 final-geometry energies,
  volume errors, and summary plots.
- `Goldzak12/`: LC10 equations of state, cohesive energies, and literature
  comparisons for ten cubic covalent and ionic solids.
- `patches/`: local CP2K and tblite patches used for the final benchmark
  revision.
- `scripts/`: helper scripts used for the final k-point, cell-optimization,
  and CP2K-native-vs-tblite-CLI checks.
- `FINAL_RESULTS.md`, `CODE_PATCHES.md`, and `paper_revision_numbers.csv`:
  compact provenance for the current paper revision.

Generated CP2K working directories, raw standard-output files, and optional
diagnostic plots are not tracked. They can be recreated from the versioned
inputs and scripts; the curated CSV, JSON, plotting data, and manuscript figure
files are the benchmark data used in the paper.

## Current revision snapshot

The final calculations use DCM-Uni-Paderborn CP2K development trunk revision
`faf9aae91266170dfee8a9f7171a5135bc5eb368` with tblite support. The tblite
build combines `main` revision `eb50bbfbe1c0869e2e18c9b7cc13144e5130b6df`
with PR 350 head `8c5e56255dc0f7001615489f24162ed770888d8b` in local merge
`8a9d09474b93d25c044d6f46ce920750c7fe4cf7`; PR 343 is already in the base.
The frozen CP2K and tblite executable SHA-256 hashes are
`f2b8e6e516b60d49af722997dd0bf06c10b54b2a2a221f786e5eaea38cccd8a5`
and `d50145af569a6ce4ea4e73e68d1cb004c3ca240105deb941c0244b7d431ed47f`.

Primary aggregate results:

| Benchmark | Setup | Method | MAE |
|---|---|---|---:|
| DMC-ICE13 relative energies | native Bloch 3x3x3 | GFN1-xTB | 8.005255 kJ mol-1 |
| DMC-ICE13 relative energies | native Bloch 3x3x3 | GFN2-xTB | 3.462919 kJ mol-1 |
| X23b lattice energies | k333 SP on native Bloch k222 cell opt | GFN1-xTB | 11.345702 kJ mol-1 |
| X23b lattice energies | k333 SP on native Bloch k222 cell opt | GFN2-xTB | 14.092104 kJ mol-1 |
| X23b cell volumes | native Bloch k222 cell opt | GFN1-xTB | 7.514116 percent |
| X23b cell volumes | native Bloch k222 cell opt | GFN2-xTB | 5.842296 percent |
| LC10 lattice constants | k444 EOS, 10/10 | GFN1-xTB | 0.145118 A |
| LC10 lattice constants | k444 EOS, 10/10 | GFN2-xTB | 0.062410 A |
| LC10 cohesive energies | k555 on k444 EOS minima, 10/10 | GFN1-xTB | 1.543851 eV atom-1 |
| LC10 cohesive energies | k555 on k444 EOS minima, 10/10 | GFN2-xTB | 1.299325 eV atom-1 |

All production k-point calculations use native Bloch sampling with full
SPGLIB symmetry reduction. The completed production counts are 156/156 for
DMC-ICE13, 46/46 for X23b k222 cell optimization, and 46/46 each for the
X23b k333 and k444 final-geometry single points.
