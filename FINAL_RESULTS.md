# Periodic GFN2-xTB Paper Revision Package

> This file records the earlier DMC13/X23b production build. The later,
> independently versioned Klimes-23 build and results are documented in
> `Klimes-Solids23/README.md`; do not mix the two build manifests.

Generated: 2026-07-06

This package collects the final local results used to revise `main_jctc_revised.tex` for the paper "Advancing GFN2-xTB for Periodic Systems via multipolar Ewald Summation in CP2K".

## Build State

- CP2K binary: `/private/tmp/cp2k-pr-5533-build-cellhess4/bin/cp2k.psmp`
- CP2K revision: `518a50992f009b083c127372f294e6485306c05b` (`CP2K version 2026.1 (Development Version)`, trunk build)
- CP2K flags include: `spglib`, `mctc-lib`, `tblite`
- tblite source: `/private/tmp/gxtb-clean-pr350-20260702/tblite-src`
- tblite HEAD: `5b14b8430bb2ffb3c96808466ad670821f81f745`
- tblite includes PR343 commit `a32675a` and PR350 commit `c36ee59`
- tblite binary: `/private/tmp/gxtb-clean-pr350-20260702/tblite-install/bin/tblite` (`tblite version 0.6.0`)

## Code Patchsets

- `patches/tblite_wsc_multipole_ewald_local.patch`: Wigner-Seitz/multipole Ewald/cutoff fixes in tblite.
- `patches/cp2k_tblite_interface_local.patch`: CP2K tblite-interface mixer-setting forwarding for the current tblite API.
- `scripts/`: exact benchmark helper scripts used for DMC13, X23b cellopt variants, X23b final k-point single-points, and CP2K-native-vs-CLI checks.
- `tex/`: synchronized copies of the updated Overleaf manuscript and SI files.
- `figures/`: synchronized PDF figures regenerated from the final DMC13 and X23b data.

Interpretation for PR preparation:

- The tblite patch is the main physics fix candidate: WSC image indexing, multipolar Ewald cutoffs, direct/reciprocal multipole matrix and virial/gradient consistency.
- The CP2K patch is an interface-maintenance fix: mixer settings are forwarded into `tb%calc%mixer_input` for the current tblite API.
- The X23b `KEEP_ANGLES` improvement is a benchmark protocol choice, not a physics-code patch.

## Final DMC-ICE13 Result

Final DMC13 numbers are from `/Users/tkuehne/Documents/g-xTB/dmc13_nativebloch_final_20260705/DMC-ICE13/data/dmc_ice13_kpoint_stats.csv`.

| Mesh | Method | ME | MAE | RMSE | MaxAE |
|---|---|---:|---:|---:|---:|
| Gamma | GFN1-xTB | -5.288685 | 6.696681 | 7.155635 | 10.475948 |
| Gamma | GFN2-xTB | 1.372619 | 5.355715 | 10.382540 | 34.362751 |
| k333 | GFN1-xTB | -8.008187 | 8.008187 | 8.653532 | 13.666572 |
| k333 | GFN2-xTB | -2.383503 | 3.185301 | 3.648584 | 6.792092 |
| k555 | GFN1-xTB | -8.009417 | 8.009417 | 8.654820 | 13.666818 |
| k555 | GFN2-xTB | -2.385201 | 3.183706 | 3.646806 | 6.792092 |

Compared with the previous paper numbers, the k333 MAE improves from 13.112632 to 8.008187 kJ mol-1 for GFN1-xTB and from 8.877913 to 3.185301 kJ mol-1 for GFN2-xTB.

## Final X23b Cell Optimization Result

Recommended relaxed-cell protocol for this revision: CP2K-native Bloch \(2^3\)-mesh cell optimization with angle conservation. The final GFN2/cytosine point uses CP2K numerical stress and the relaxed maximum-gradient threshold discussed in the benchmark logs; CP2K writes the final CIF with `Optimization converged: TRUE`.

| Quantity | Method | N | ME | MAE | RMSE | MaxAE |
|---|---|---:|---:|---:|---:|---:|
| Lattice energy / kJ mol-1 | GFN1-xTB | 23 | -0.710350 | 11.129018 | 14.815022 | 38.461935 |
| Lattice energy / kJ mol-1 | GFN2-xTB | 23 | -11.986205 | 14.459836 | 23.427726 | 90.069662 |
| Cell volume / percent | GFN1-xTB | 23 | -5.815548 | 7.914787 | 9.478583 | 19.084442 |
| Cell volume / percent | GFN2-xTB | 23 | -1.791691 | 5.616637 | 7.363220 | 18.421974 |

Key point: full native-Bloch k-point cell optimization reduces the relaxed-cell lattice-energy MAEs strongly relative to the Gamma-only cellopt table. GFN2-xTB gives the better volume statistics, while GFN1-xTB remains better in aggregate lattice-energy MAE on this benchmark.

## Verification

- tblite CLI finite-difference checks: final CLI virial max errors of order 1e-8 to 1e-7; force max errors below 5e-9 in the small test set.
- CP2K-native finite-difference checks: Gamma and k-point stress/force checks pass at about 1e-6 to 4e-6 stress-component scale, with small force-component differences.
- CP2K-native vs tblite CLI on X23b keep-angle final geometries:
  - max energy difference: 1.3918e-08 Eh (GFN1), 1.1903e-08 Eh (GFN2)
  - max force-component difference: 2.8001e-07 (GFN1), 5.6153e-07 (GFN2)
  - max virial-component difference: 2.6810e-06 (GFN1), 6.6761e-06 (GFN2)

## Paper Guidance

Recommended wording changes:

- DMC13: emphasize a clear improvement from GFN1-xTB to GFN2-xTB after the current multipolar Ewald/k-point fixes. Use k333/k555 as converged native Bloch k-point data.
- X23b relaxed cells: use the native Bloch \(2^3\)-mesh cellopt table generated in `x23b_k222_cellopt_nativebloch_final_20260705_merged`.
- Code: present tblite multipolar Ewald/WSC changes as the main implementation correction; keep CP2K interface details and benchmark protocol separate.
