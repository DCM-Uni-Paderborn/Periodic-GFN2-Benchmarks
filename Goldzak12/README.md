LC10 native-Bloch CP2K/tblite benchmark
========================================

LC10 is a ten-solid subset of the cubic covalent and ionic systems reported
by Goldzak, Wang, Ye, and Berkelbach, J. Chem. Phys. 157, 174112 (2022). It
contains C, Si, SiC, BN, BP, AlN, AlP, MgS, LiF, and LiCl and compares
CP2K/tblite GFN1-xTB and GFN2-xTB with published post-HF, DFT, and
zero-point-corrected experimental lattice constants and cohesive energies.

Current production data
-----------------------

The underlying calculations use:

- CP2K trunk revision `faf9aae91266170dfee8a9f7171a5135bc5eb368` with the
  local CP2K/tblite interface patch recorded in `data/build_provenance.json`;
- tblite revision `8a9d09474b93d25c044d6f46ce920750c7fe4cf7`, combining current
  `main` with PR #350 and the previously merged PR #343;
- conventional cubic eight-atom cells and CP2K native Bloch sampling through
  `&KPOINTS`, with full SPGLIB symmetry reduction and no Born-von-Karman
  supercells;
- `k444` equations of state and `k555` cohesive-energy single points on the
  fitted minima;
- matching tblite CLI isolated-atom references.

All 20 equations of state and all 60 final `k333/k444/k555` single points are
complete. The reported LC10 statistics are:

| method | coverage | lattice MAE (A) | cohesive-energy MAE (eV/atom) |
|---|---:|---:|---:|
| GFN1-xTB | 10/10 | 0.145118 | 1.543851 |
| GFN2-xTB | 10/10 | 0.062410 | 1.299325 |

Analysis-only reproduction
--------------------------

The current tables and figures can be rebuilt from the existing outputs
without launching CP2K or tblite calculations:

```bash
python3 Goldzak12/scripts/run_goldzak12_eos_benchmark.py \
  --analysis-only --eos-mesh k444 \
  --energy-mesh k333 --energy-mesh k444 --energy-mesh k555 \
  --result-mesh k555

python3 Goldzak12/scripts/compare_goldzak12_results.py --mesh k555
python3 Goldzak12/scripts/plot_literature_comparison.py
python3 Goldzak12/scripts/validate_goldzak12_results.py \
  --eos-mesh k444 --energy-mesh k333 --energy-mesh k444 \
  --energy-mesh k555 --result-mesh k555
```

Raw calculations and generated inputs are kept below `Goldzak12/runs` and
`Goldzak12/inputs` and are ignored by Git. Curated tables, references, and
provenance are versioned in `Goldzak12/data`; manuscript figures are in
`Goldzak12/figures`.

Literature comparison
---------------------

`scripts/plot_literature_comparison.py` recomputes every error against the
same zero-point-corrected experimental reference. It includes the post-HF
values from Goldzak et al., SCAN-family values from Mejia-Rodriguez and
Trickey (2020), and additional solid-state DFT values from Mo et al. (2017).
