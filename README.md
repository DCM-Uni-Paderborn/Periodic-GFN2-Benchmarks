# Periodic GFN2 Benchmarks

This repository collects the paper-relevant benchmark inputs, curated output
data, analysis scripts, and figures for periodic GFN calculations in CP2K.

## Contents

- `DMC-ICE13/`: CP2K/tblite single-point benchmark for the DMC-ICE13 ice
  polymorph data set, comparing periodic GFN1-xTB and GFN2-xTB relative
  energies with the diffusion Monte Carlo reference energies.
- `X23b/`: CP2K/tblite molecular-crystal benchmark with gas-phase molecular
  optimizations, crystal single-point k-point tests, Gamma-point cell
  optimizations, extracted X23b lattice energies, volume errors, and summary
  plots.

Generated CP2K working directories and raw standard-output files are not
tracked. They can be recreated from the versioned inputs and scripts; the
curated CSV, JSON, and plotting data files are the benchmark data used in the
manuscript.
