# Code Patch Inventory

The final benchmark build is reconstructed from the two self-contained patches
in `patches/`. Obsolete intermediate backports are not part of the branch.

## tblite

Patch: `patches/tblite_main_pr350_wsc_derivatives.patch`

Base: tblite `main` at `eb50bbfbe1c0869e2e18c9b7cc13144e5130b6df`.
The patch includes PR 350 at `8c5e56255dc0f7001615489f24162ed770888d8b`
and the final Wigner-Seitz corrections and analytical derivatives. PR 343 is
already present in the base revision.

Main changes:

- Preserve true lattice-translation indices for Wigner-Seitz images.
- Use smooth normalized weights for competing nearest images.
- Include weight derivatives in periodic charge and multipole forces and
  virials.
- Correct the damped dipole-dipole radial factor covered by PR 350.
- Add focused energy, force, virial, and boundary-continuity regression tests.

## CP2K

Patch: `patches/cp2k_trunk_tblite_full_symmetry_scc.patch`

Base: DCM-Uni-Paderborn CP2K development trunk at
`faf9aae91266170dfee8a9f7171a5135bc5eb368`.

Main changes:

- Native periodic CP2K/tblite Bloch energies, forces, and stress tensors.
- Full SPGLIB k-point symmetry reduction, including fractional translations
  for nonsymmorphic operations.
- PBC-consistent atomic gauges and overlap-covariant Bloch phases for atoms on
  fractional-cell boundaries.
- Dynamic symmetry refresh during cell optimization.
- Separate Bloch-wavefunction and tblite SCC restart handling.
- tblite and CP2K SCC mixer propagation and corrected Broyden history.
- Gamma, reduced/full-grid k-point, FD force/stress, boundary, restart,
  nonsymmorphic, cell-optimization, and native-DFT regression coverage.

## Validation

- Both patches pass `git apply --check` against their listed bases.
- tblite: 35/35 enabled non-C-API test groups pass in the no-ddX build.
- CP2K: 625/625 xTB regression matchers pass; the focused periodic block is
  131/131.
- Reduced and unreduced Urea native-Bloch energies agree within about
  `7e-15` hartree; the native-DFT control agrees within `1e-15` hartree.
- On all 23 final Gamma-optimized X23b structures per method, CP2K-native and
  the current tblite CLI agree to maximum energy differences of
  `1.34e-8` hartree (GFN1) and `1.01e-8` hartree (GFN2), maximum gradient
  component differences of `2.55e-7` and `8.95e-7` atomic units, and maximum
  virial component differences of `1.35e-6` and `4.45e-6` atomic units.
- Standalone tblite CLI analytical forces and virials agree with finite
  differences to maximum component residuals of `2.19e-8` and `1.11e-7`
  atomic units for GFN1 and `6.36e-9` and `4.93e-8` atomic units for GFN2 in
  the final periodic checks.
- Native-Bloch analytical forces and stresses agree with finite differences
  below 0.01% in the periodic regression matrix. CP2K's external
  `REFERENCE_CLI` diagnostic is a Gamma-path check and intentionally does not
  substitute a Born-von-Karman supercell for active k-point calculations.

See `patches/README.md` for exact SHA-256 hashes and application commands.

## Benchmark Scripts

The benchmark-specific runners are versioned under `DMC-ICE13/scripts/`,
`X23b/scripts/`, `Goldzak12/scripts/`, and `scripts/`. Runtime paths are command
line options or environment variables; curated tables and figures are generated
from the corresponding run trees.
