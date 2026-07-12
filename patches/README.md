# Code patches used for the final benchmark build

The benchmark reruns use CP2K development trunk and current tblite `main`,
including the complete changes from tblite PR 350. These two patches are the
exact source diffs used to build the executables recorded in the benchmark
provenance files. Obsolete intermediate backports are intentionally omitted.

## `tblite_main_pr350_wsc_derivatives.patch`

- Base revision: `eb50bbfbe1c0869e2e18c9b7cc13144e5130b6df` (tblite
  `main`; PR 343 is already included).
- PR 350 head: `8c5e56255dc0f7001615489f24162ed770888d8b`.
- Local merge revision: `8a9d09474b93d25c044d6f46ce920750c7fe4cf7`.
- Patch SHA-256:
  `1fbc37841a420cb766d3aeafac59fb1479c3bc2874152620356059518d7ff5c7`.

The patch is self-contained relative to the listed `main` revision: it includes
the complete PR 350 series and the additional Wigner-Seitz
cell corrections used by the final build: integer lattice-translation image
indices, smooth normalized image weights, and the corresponding analytical
force and stress derivatives. It also updates the affected unit and C-API
reference tests.

Apply from a tblite checkout at the base revision:

```sh
git apply /path/to/tblite_main_pr350_wsc_derivatives.patch
```

## `cp2k_trunk_tblite_full_symmetry_scc.patch`

- Base revision: `faf9aae91266170dfee8a9f7171a5135bc5eb368` (CP2K
  development trunk).
- Patch SHA-256:
  `a54705795ea3c5c3ffd6e5e2197e00329e49b5ced05d3945f7ce39e9d84e3a2a`.

The patch provides the native periodic CP2K/tblite Bloch path used in the
benchmarks, including analytical forces and stresses, full SPGLIB k-point
symmetry reduction, fractional translations of nonsymmorphic operations, and
robust overlap-phase handling for atoms on fractional-cell boundaries. It also
contains the corrected dynamic symmetry refresh for cell optimization,
independent tblite SCC restart controls, propagation of mixer settings, and the
CP2K-side Broyden history fix. Periodic atom gauges are derived from CP2K's
actual `pbc()` mapping and the two valid Bloch-phase directions are selected by
overlap covariance. This keeps full symmetry reduction stable when cell
optimization moves atoms across floating-point cell boundaries. Regression inputs cover Gamma and
k-point forces/stresses, skew cells, boundary atoms, a dynamic Urea cell
optimization, nonsymmorphic crystals, SCC restart, and a native DFT control
case.

Apply from a CP2K checkout at the base revision:

```sh
git apply /path/to/cp2k_trunk_tblite_full_symmetry_scc.patch
```

## Validation

Both patches pass `git apply --check` against the revisions listed above. All
35 non-C-API tblite test groups enabled in the no-ddX build pass. The final
CP2K validation includes 625/625 xTB matcher checks, including 16/16 focused
SPGLIB checks. Reduced and unreduced native-Bloch energies agree to about
`7e-15` hartree; the corresponding native-DFT control agrees to `1e-15`
hartree. Finite-difference checks give force/stress residuals below 0.01% for
native-Bloch k-point calculations; the skew-cell h-BN test gives summed
absolute residuals of `7e-8` and `5.21e-7` atomic units for forces and stress,
respectively. The monoclinic native-DFT control gives a summed stress residual
of `8.75e-10` atomic units. On the 46 final Gamma-optimized X23b structures,
CP2K-native and the current tblite CLI have maximum energy, gradient-component,
and virial-component differences of `1.34e-8` hartree, `2.55e-7` atomic units,
and `1.35e-6` atomic units for GFN1, and `1.01e-8` hartree, `8.95e-7` atomic
units, and `4.45e-6` atomic units for GFN2. This external CLI comparison is a
Gamma-path validation; active native-Bloch k points are validated directly
against finite differences and full-versus-reduced meshes.

The benchmark provenance records the executable and linked-library hashes in
addition to these source revisions and patch hashes.
