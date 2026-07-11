# Periodic GFN2 Benchmarks

This repository collects the paper-relevant benchmark inputs, curated output
data, analysis scripts, and figures for periodic GFN calculations in CP2K.

## Contents

- `DMC-ICE13/`: CP2K/tblite single-point benchmark for the DMC-ICE13 ice
  polymorph data set, comparing periodic GFN1-xTB and GFN2-xTB relative
  energies with the diffusion Monte Carlo reference energies. The current
  manuscript data use native Bloch k-point calculations through CP2K
  `&KPOINTS`.
- `X23b/`: CP2K/tblite molecular-crystal benchmark with gas-phase molecular
  optimizations, crystal single-point k-point tests, native Bloch 2x2x2
  crystal cell optimizations, extracted X23b lattice energies, volume errors,
  and summary plots.
- `Klimes-Solids23/`: 23 inorganic solids with native-Bloch GFN1-xTB and
  GFN2-xTB lattice constants, bulk moduli, and cohesive energies, compared
  with experiment and eight published DFT methods. This directory contains
  only Klimes-23; no Buccheri benchmark is included.
- `patches/`: build-specific local CP2K and experimental tblite patches; see
  `CODE_PATCHES.md` before applying them.
- `scripts/`: helper scripts used for the final k-point, cell-optimization,
  and CP2K-native-vs-tblite-CLI checks.
- `FINAL_RESULTS.md`, `CODE_PATCHES.md`, and `paper_revision_numbers.csv`:
  compact provenance for the current paper revision.

Generated CP2K working directories and raw standard-output files are not
tracked. They can be recreated from the versioned inputs and scripts; the
curated CSV, JSON, and plotting data files are the benchmark data used in the
manuscript.

## DMC13/X23b revision snapshot

The final calculations use CP2K trunk revision
`518a50992f009b083c127372f294e6485306c05b` with tblite support and tblite
revision `5b14b8430bb2ffb3c96808466ad670821f81f745` (`tblite` 0.6.0),
including the changes corresponding to tblite PRs 343 and 350.

The later Klimes-23 calculations use CP2K trunk
`faf9aae91266170dfee8a9f7171a5135bc5eb368` and tblite
`8a9d09474b93d25c044d6f46ce920750c7fe4cf7`. Its independent manifest and
interpretation are in `Klimes-Solids23/README.md` and
`Klimes-Solids23/provenance.json`.

Primary aggregate results:

| Benchmark | Setup | Method | MAE |
|---|---|---|---:|
| DMC-ICE13 relative energies | native Bloch 3x3x3 | GFN1-xTB | 8.008187 kJ mol-1 |
| DMC-ICE13 relative energies | native Bloch 3x3x3 | GFN2-xTB | 3.185301 kJ mol-1 |
| X23b lattice energies | native Bloch 2x2x2 cell opt | GFN1-xTB | 11.129018 kJ mol-1 |
| X23b lattice energies | native Bloch 2x2x2 cell opt | GFN2-xTB | 14.459836 kJ mol-1 |
| X23b cell volumes | native Bloch 2x2x2 cell opt | GFN1-xTB | 7.914787 percent |
| X23b cell volumes | native Bloch 2x2x2 cell opt | GFN2-xTB | 5.616637 percent |
| Klimes-23 lattice constants, paired n=21 | native Bloch 8x8x8/16x16x16 EOS | GFN1-xTB | 0.487065 A |
| Klimes-23 lattice constants, paired n=21 | native Bloch 8x8x8/16x16x16 EOS | GFN2-xTB | 0.560231 A |
| Klimes-23 cohesive energies, paired n=21 | native Bloch 8x8x8/16x16x16 EOS | GFN1-xTB | 1.581526 eV atom-1 |
| Klimes-23 cohesive energies, paired n=21 | native Bloch 8x8x8/16x16x16 EOS | GFN2-xTB | 0.803420 eV atom-1 |
