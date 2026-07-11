# Code Patch Inventory

## Current Klimes-23 build (2026-07-11)

- CP2K trunk: `faf9aae91266170dfee8a9f7171a5135bc5eb368`.
- tblite: `8a9d09474b93d25c044d6f46ce920750c7fe4cf7`, which contains
  PR #343 (`a32675a`) and PR #350 (`8c5e562`).
- Applied local patch: `patches/cp2k_tblite_interface_local.patch`.
- `patches/tblite_wsc_multipole_ewald_local.patch` was **not** applied to the
  Klimes-23 production build. It is retained as the earlier experimental
  DMC13/X23b physics patch described below.

Exact Klimes-23 build and run provenance is recorded in
`Klimes-Solids23/provenance.json`.

## Earlier experimental tblite patch

Patch: `patches/tblite_wsc_multipole_ewald_local.patch`

Base/HEAD: `5b14b8430bb2ffb3c96808466ad670821f81f745`.

Files changed:

- `src/tblite/coulomb/multipole.f90`
- `src/tblite/coulomb/ewald.f90`
- `src/tblite/wignerseitz.f90`
- `src/tblite/cutoff.f90`

Purpose:

- Correct Wigner-Seitz image indexing and image weighting for multipolar electrostatics.
- Use multipole-aware Ewald real/reciprocal cutoff estimates.
- Use WSC images consistently in multipole matrix and gradient/virial paths.
- Respect directional periodicity masks for cutoff and central-cell wrapping.

## Current CP2K compatibility patch

Patch: `patches/cp2k_tblite_interface_local.patch`

Base/HEAD: `faf9aae91266170dfee8a9f7171a5135bc5eb368`.

File changed:

- `src/tblite_interface.F`

Purpose:

- Forward mixer settings into `tb%calc%mixer_input` for the current tblite API.

## Benchmark Scripts

Full helper scripts used for the final revision are in `scripts/`.
Some runner defaults preserve the local production paths used for the paper
revision; override `--benchmark-root`, `--out`, or `CP2K` when replaying them
in a different checkout. `scripts/update_x23b_k222_figures.py` is repo-relative
and regenerates the X23b figures from the versioned `X23b/data` files.

- `run_dmc13_kpoint_jobs.py`: DMC13 native Bloch k-point benchmark runner.
- `run_x23b_cellopt_variant_matrix.py`: X23b Gamma cellopt variant runner, including `cg_2pnt_keep_angles`.
- `run_x23b_cellopt_final_kpoint_sp.py`: X23b final-cellopt native Bloch k222/k333 single-point runner.
- `run_x23b_reference_cli_checks.py`: CP2K-native vs tblite CLI reference comparison.
- `run_literature_eos_benchmarks.py`: Klimes-23 native-Bloch EOS runner and Murnaghan analysis.
- `audit_klimes_results.py`: normalize saved energies to CP2K's extrapolated T=0 values.
- `build_klimes_report.py`: regenerate Klimes-23 tables and figures.
- `verify_klimes_benchmark.py`: validate all Klimes sidecars, fits, tables, figures, and hashes.
