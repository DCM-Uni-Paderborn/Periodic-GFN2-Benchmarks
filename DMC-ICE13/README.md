# DMC-ICE13 Periodic GFN Benchmark

This directory contains CP2K/tblite single-point calculations for the
DMC-ICE13 ice polymorph benchmark. The calculations compare periodic
GFN1-xTB and GFN2-xTB relative energies against the diffusion Monte Carlo
reference values of Della Pia, Zen, Alfe, and Michaelides,
J. Chem. Phys. 157, 134701 (2022), DOI: 10.1063/5.0102645.

## Data included

- `poscars/`: POSCAR geometries for the 13 DMC-ICE13 polymorphs.
- `inputs/`: Gamma-only CP2K input files for GFN1-xTB and GFN2-xTB.
- `kpoint_inputs/`: explicit native Bloch 1x1x1, 2x2x2, 3x3x3, 4x4x4, and
  5x5x5 MacDonald k-point CP2K input files.
- `runs/`: generated Gamma-only CP2K working directories, ignored by Git.
- `runs_kpoints/`: generated k-point CP2K working directories, ignored by Git.
- `data/results.json`: raw CP2K total energies, per-water energies, relative
  energies with respect to ice Ih, and error statistics for the Gamma-only
  calculations.
- `data/kpoint_results.json`: raw and relative energies for the k-point
  dependent calculations.
- `data/dmc_ice13_relative_energies.csv`: 3x3x3 relative energies and GFN
  errors used as the primary manuscript values.
- `data/dmc_ice13_kpoint_stats.csv`: aggregate DMC-ICE13 error statistics as a
  function of k-point mesh.
- `data/dmc_ice13_kpoint_relative_energies.csv`: phase-resolved relative
  energies and errors as a function of k-point mesh.
- `data/previous_vs_full_pr350_mae.csv` and the companion Markdown file:
  explicit comparison with the earlier partial-PR350 manuscript stack.
- `data/dmc_ice13_relative_mae_comparison.csv`: comparison with the published
  DFT data from the DMC-ICE13 paper.
- `data/dmc_ice13_published_dft_absolute_energies.csv`: published DMC and DFT
  absolute lattice energies from the DMC-ICE13 paper, used to compute the
  relative-energy MAE ranking.
- `data/build_provenance.json`: source revisions, executable and shared-library
  hashes, patch hashes, build flags, and the completed-calculation count.
- `data/dmc_ice13_reference_cli_rows.csv` and
  `data/dmc_ice13_reference_cli_summary.csv`: direct CP2K-native versus tblite
  CLI energy, gradient, and virial checks for all 26 Gamma calculations.
- `figures/`: PDF, SVG, and PNG plots generated from the benchmark data.
- `scripts/`: input generation, extraction, analysis, plotting, and run scripts.

The original PDF and Supporting Information are not redistributed here. The
geometries and DMC reference values are documented through the paper DOI above.

## CP2K setup used

The calculations were run from CP2K development trunk, not from a numbered
release. The executable reports `2026.1 (Development Version)` and is
interfaced to tblite:

- CP2K source revision: `faf9aae91266170dfee8a9f7171a5135bc5eb368`
- CP2K flags reported by the executable: `omp no_statm_access spglib libdftd4
  dftd4_v4_2 s_dftd3 mctc-lib tblite`
- tblite: `8a9d09474b93d25c044d6f46ce920750c7fe4cf7` (`tblite` 0.6.0),
  merging current `main` with the complete tblite PR 350 series through
  `8c5e562`; the earlier PR 343 changes are also included
- CP2K working-tree additions: overlap-covariant native-Bloch full symmetry,
  analytical force/stress checks, separate Bloch-wavefunction and multipolar
  SCC restart handling, and Broyden-mixer history fixes
- `TBLITE/ACCURACY`: `0.1`
- `EPS_SCF`: `1.0E-9`
- run-script defaults: `OMP_NUM_THREADS=1`,
  `OPENBLAS_NUM_THREADS=1`, `MKL_NUM_THREADS=1`, and
  `CP2K_PARALLEL_JOBS=20`, i.e. independent single-core CP2K jobs are launched
  concurrently.

The primary comparison uses the Gamma-centered 3x3x3 k-point mesh, matching
the non-hybrid DFT single-point setup in the DMC-ICE13 reference. The explicit
1x1x1 mesh verifies equivalence to the Gamma-only calculation, the 2x2x2 mesh
documents the approach to convergence, and the 4x4x4 and 5x5x5 checks confirm
that the 3x3x3 aggregate statistics are converged. All energies in the CSV
summaries are relative to ice Ih and reported in kJ mol-1 per water molecule.
The independent reference-CLI checks use analytical CP2K stress tensors; all
26 calculations complete, with maximum force-component differences of
`3.57e-8` (GFN1) and `1.51e-7` atomic units and maximum virial-component
differences of `1.26e-6` and `2.29e-6` atomic units, respectively.

Current aggregate MAEs:

| Mesh | GFN1-xTB | GFN2-xTB |
|---|---:|---:|
| Gamma | 6.694624 | 5.578897 |
| 2x2x2 | 7.956838 | 3.510100 |
| 3x3x3 | 8.005255 | 3.462919 |
| 4x4x4 | 8.006494 | 3.461424 |
| 5x5x5 | 8.006485 | 3.461353 |
