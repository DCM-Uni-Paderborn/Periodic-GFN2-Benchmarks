# Periodic GFN2-xTB Paper Revision Package

Generated: 2026-07-12

This package records the final data used for the revised manuscript
"Advancing GFN2-xTB for Periodic Systems via multipolar Ewald Summation in
CP2K".

## Build State

- CP2K development trunk: `faf9aae91266170dfee8a9f7171a5135bc5eb368`
- tblite `main`: `eb50bbfbe1c0869e2e18c9b7cc13144e5130b6df`
- tblite PR 350 head: `8c5e56255dc0f7001615489f24162ed770888d8b`
- tblite local merge: `8a9d09474b93d25c044d6f46ce920750c7fe4cf7`
- CP2K executable SHA-256: `f2b8e6e516b60d49af722997dd0bf06c10b54b2a2a221f786e5eaea38cccd8a5`
- tblite executable SHA-256: `d50145af569a6ce4ea4e73e68d1cb004c3ca240105deb941c0244b7d431ed47f`
- linked `libcp2k` SHA-256: `7656b6154614290d3e121f2f7a2d527b5e0e5d128eacf456e3c35e847894741d`
- linked `libtblite` SHA-256: `f7bd2a841543dcbb71da0954f2e2bf016b7a202dace377f26f72ac139035c3e6`

The two self-contained source diffs are
`patches/tblite_main_pr350_wsc_derivatives.patch` and
`patches/cp2k_trunk_tblite_full_symmetry_scc.patch`. All production k-point
calculations use native Bloch sampling with full SPGLIB symmetry reduction.

## DMC-ICE13

All 156/156 calculations completed. Statistics are relative to ice Ih over
the twelve non-reference phases, in kJ mol-1 per water molecule.

| Mesh | Method | ME | MAE | RMSE | MaxAE |
|---|---|---:|---:|---:|---:|
| Gamma | GFN1-xTB | -5.285749 | 6.694624 | 7.153219 | 10.470050 |
| Gamma | GFN2-xTB | 0.906983 | 5.578897 | 10.188726 | 33.462413 |
| k333 | GFN1-xTB | -8.005255 | 8.005255 | 8.650658 | 13.663936 |
| k333 | GFN2-xTB | -2.837920 | 3.462919 | 3.812981 | 6.346296 |
| k555 | GFN1-xTB | -8.006485 | 8.006485 | 8.651947 | 13.664182 |
| k555 | GFN2-xTB | -2.839590 | 3.461353 | 3.811621 | 6.346296 |

Relative to the previous manuscript values, the k333 GFN1 MAE changes from
8.008187 to 8.005255 kJ mol-1 and the GFN2 MAE from 3.185301 to 3.462919
kJ mol-1.

## X23b

All 46/46 native-Bloch k222 angle-conserving cell optimizations completed.
All 46/46 k333 and 46/46 k444 single points on the final k222 geometries also
completed. Volumes therefore come from k222 optimization, while the reported
lattice energies use the converged k333 reevaluation.

| Quantity | Method | N | ME | MAE | RMSE | MaxAE |
|---|---|---:|---:|---:|---:|---:|
| Lattice energy, k333 on k222 geometry / kJ mol-1 | GFN1-xTB | 23 | 0.258871 | 11.345702 | 14.019344 | 30.935058 |
| Lattice energy, k333 on k222 geometry / kJ mol-1 | GFN2-xTB | 23 | -12.018989 | 14.092104 | 21.341752 | 77.785392 |
| Cell volume, k222 optimization / percent | GFN1-xTB | 23 | -5.960071 | 7.514116 | 9.019708 | 19.236681 |
| Cell volume, k222 optimization / percent | GFN2-xTB | 23 | -1.657324 | 5.842296 | 7.530373 | 19.952589 |

The k333-to-k444 mean absolute changes are 0.079329 kJ mol-1 for GFN1-xTB
and 0.084265 kJ mol-1 for GFN2-xTB. The corresponding k444 MAEs are 11.366118
and 14.166994 kJ mol-1.

## LC10

The k444 equation-of-state fits are valid for all ten solids with both
methods. Cohesive energies use k555 single points on the k444 minima.

| Method | N | Lattice-constant MAE / A | Cohesive-energy MAE / eV atom-1 |
|---|---:|---:|---:|
| GFN1-xTB | 10 | 0.145118 | 1.543851 |
| GFN2-xTB | 10 | 0.062410 | 1.299325 |

## Verification

- tblite: 35/35 enabled no-ddX test groups pass.
- CP2K: 625/625 xTB regression matchers and 131/131 focused periodic
  matchers pass.
- Standalone tblite finite differences give maximum force/virial residuals of
  `2.19e-8`/`1.11e-7` atomic units for GFN1 and
  `6.36e-9`/`4.93e-8` atomic units for GFN2.
- Native-Bloch k-point analytical forces and stresses agree with finite
  differences to below 0.01 percent.
- On all final Gamma X23b geometries, CP2K-native and tblite CLI energies,
  gradients, and virials agree to the limits recorded in `CODE_PATCHES.md`.
