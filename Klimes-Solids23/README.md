# Klimes-23 native-Bloch benchmark

This directory contains only the 23-solids benchmark of Klimes, Bowler, and
Michaelides. No Buccheri benchmark data or analysis is included.

## Result at a glance

The fair GFN1-xTB/GFN2-xTB comparison uses the 21 systems for which both
methods have a valid Murnaghan EOS fit.

| Property | GFN1-xTB MAE | GFN2-xTB MAE | GFN1-xTB MARE | GFN2-xTB MARE | GFN2 better/worse |
|---|---:|---:|---:|---:|---:|
| Lattice constant | 0.4871 A | 0.5602 A | 9.98% | 10.65% | 11/10 |
| Bulk modulus | 498.00 GPa | 507.69 GPa | 544.08% | 1500.68% | 8/13 |
| Cohesive energy | 1.5815 eV/atom | 0.8034 eV/atom | 96.33% | 35.62% | 16/5 |

GFN2-xTB clearly improves cohesive energies on this set, approximately halving
the paired MAE. It does not improve the aggregate lattice constants, and the
bulk moduli remain much too large. The latter problem is especially severe for
the soft alkali and alkaline-earth metals. GFN2-xTB improves several individual
ionic and semiconductor results, including the LiF and NaF lattice constants
and the Ge and GaAs bulk moduli, but worsens K, Rb, and Cs strongly.

The eight published DFT methods remain substantially more accurate for all
three aggregate properties. This benchmark is therefore useful as a demanding
transferability diagnostic and gives a favorable GFN2-vs-GFN1 result only for
cohesive energies, not as an overall claim against DFT.

## Coverage

- GFN1-xTB: 23/23 accepted EOS fits.
- GFN2-xTB: 21/23 accepted EOS fits.
- GFN2-xTB Cu: electronically collapsed energy branch; no bracketed solid EOS.
- GFN2-xTB MgO: no bracketed minimum; additional compression runs exhibit SCC
  collapse and are retained only as compact diagnostics.
- Raw sidecars: 867 total, 866 successful and one failed compressed GaAs point.
- Final fit matrix: 44/46 accepted method/system combinations.

The GFN2-xTB GaAs fit is valid despite the failed `V/V0=0.80` point because its
continuous local minimum is bracketed by the independently converged
`V/V0=1.00-1.40` points.

## Protocol

- Conventional cells and reference structures from the original paper.
- CP2K native Bloch k-points, never a Born-von Karman supercell path.
- `16x16x16` full grids for metals and `8x8x8` full grids for ionic and
  covalent solids, matching the published protocol.
- Eleven standard volume factors from 0.80 to 1.20; adaptive extensions were
  added only where the local solid minimum was not bracketed.
- Independent CP2K energy calculations with 300 K Fermi-Dirac smearing;
  analysis uses CP2K's energy extrapolated to zero electronic temperature.
- Published Murnaghan equation of state, with at least seven points surrounding
  a local minimum, positive curvature, and `R^2 >= 0.999`.
- Spin-polarized isolated-atom tblite calculations for cohesive energies.

Some s-block systems develop disconnected or nearly degenerate electronic
branches at large volume. The analysis selects the continuous local solid
branch nearest the reference cell rather than treating a collapsed SCC branch
as an equilibrium phase. In particular, GFN1-xTB Ca has two minima separated by
less than 1 meV/atom. Those extended-cell results are numerically reproducible
but should be interpreted as branch-sensitive model behavior.

## Build provenance

- CP2K trunk revision: `faf9aae91266170dfee8a9f7171a5135bc5eb368`.
- CP2K banner: `2026.1 (Development Version)`; this is the trunk commit above,
  not the official 2026.1 release tarball.
- CP2K flags include `tblite`; binary SHA-256:
  `f5eb42ab68102490db4914b98601512b37bd2d283eb2e5fbf3a6ac3cd17c8785`.
- tblite revision: `8a9d09474b93d25c044d6f46ce920750c7fe4cf7`.
- This tblite revision contains PR #343 (`a32675a`) and PR #350 (`8c5e562`).
- tblite binary SHA-256:
  `1136763738cc7c0b095829344ccc2949d77bc76ad804365ec2cf5839437d109c`.
- Applied local CP2K API-compatibility patch:
  `../patches/cp2k_tblite_interface_local.patch`.
- The older experimental `tblite_wsc_multipole_ewald_local.patch` was not
  applied to this Klimes-23 build.

Machine-readable details are in [provenance.json](provenance.json).

## Outputs

- [Full per-system table](results/benchmark_tables.md)
- [Machine-readable system summary](results/system_summary.csv)
- [Paired GFN statistics](results/paired_gfn_statistics.csv)
- [All GFN and published DFT values in long form](results/literature_comparison_long.csv)
- [EOS fits](results/eos_fits.csv)
- [Aggregate statistics](results/aggregate_statistics.csv)
- [System-wise relative errors](figures/klimes23-system-relative-errors.pdf)
- [Aggregate comparison with DFT](figures/klimes23-aggregate-mare.pdf)
- [GFN2-xTB change relative to GFN1-xTB](figures/klimes23-gfn2-change.pdf)

## Reproduction

```bash
python3 scripts/audit_klimes_results.py --dry-run
python3 scripts/run_literature_eos_benchmarks.py analyse
python3 scripts/build_klimes_report.py
python3 scripts/verify_klimes_benchmark.py
python3 -m unittest discover -s tests -p 'test_*.py' -v
```

The runner is restartable and automatically skips successful sidecars. Its
default binary paths record the production machine; use `--cp2k` and `--tblite`
when replaying with another installation. To recompute selected points rather
than use the versioned sidecars, pass `--force`, for example:

```bash
python3 scripts/run_literature_eos_benchmarks.py run --force --systems LiF --methods GFN2 --jobs 8
python3 scripts/run_literature_eos_benchmarks.py atoms --methods GFN2
python3 scripts/audit_klimes_results.py --require-output
```
