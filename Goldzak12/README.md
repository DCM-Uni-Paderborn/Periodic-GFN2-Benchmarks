LC12 (Goldzak12) native-Bloch CP2K/tblite benchmark
===================================================

This directory contains the 12 cubic covalent and ionic solids studied by
Goldzak, Wang, Ye, and Berkelbach, J. Chem. Phys. 157, 174112 (2022). It
compares CP2K/tblite GFN1-xTB and GFN2-xTB with the reported HF, MP2, SCS-MP2,
SOS-MP2, and zero-point-corrected experimental lattice constants and cohesive
energies.

Current production run
----------------------

The 2026-07-11 rerun used:

- CP2K trunk revision `faf9aae91266170dfee8a9f7171a5135bc5eb368` with the
  local CP2K/tblite interface patch recorded by hash in
  `data/build_provenance.json`;
- tblite revision `8a9d09474b93d25c044d6f46ce920750c7fe4cf7`, which combines current
  `main` with PR #350 and includes the previously merged PR #343;
- conventional cubic eight-atom cells and CP2K native Bloch sampling through
  `&KPOINTS` with `SCHEME MACDONALD` and `FULL_GRID T`; no Born-von-Karman
  supercells;
- a `k444` cubic equation of state, `k333/k444/k555` final single points, and
  `k444` as the reported cohesive-energy mesh;
- CP2K energies extrapolated to electronic temperature `T -> 0`;
- matching tblite CLI isolated-atom references with explicit atomic spins and
  `ACCURACY 0.05`.

Only completed SCF points enter an EOS fit. A quadratic fit is rejected when
its local RMSE exceeds 0.02 hartree or when its fitted minimum lies more than
0.02 hartree above the sampled local minimum. The current run gives 12/12
valid GFN1 fits and 10/12 valid GFN2 fits. GFN2/MgO has no bracketed stable
minimum on the compressed branch, while GFN2/LiH has a discontinuous EOS and
fails the general fit-quality criterion.

Current versus previous results
-------------------------------

| method | coverage | lattice MAE (A) | cohesive-energy MAE (eV/atom) |
|---|---:|---:|---:|
| GFN1 current | 12/12 | 0.137628 | 1.450563 |
| GFN1 previous | 12/12 | 0.164341 | 1.455859 |
| GFN2 current | 10/12 | 0.062599 | 1.293172 |
| GFN2 previous | 11/12 | 0.147638 | 1.731839 |

On the identical ten-system GFN2 subset, the lattice-constant MAE decreases
from 0.133264 to 0.062599 A and the cohesive-energy MAE decreases from 1.534124
to 1.293172 eV/atom. The frozen previous tables are in
`data/baseline_20260710`; `data/old_vs_new.md` and the associated CSV files
contain the complete per-system comparison.

Reproduction
------------

Run from the repository root:

```bash
python3 Goldzak12/scripts/run_goldzak12_eos_benchmark.py \
  --cp2k /Users/tkuehne/gxtb-local-build/install/cp2k/bin/cp2k.ssmp \
  --tblite /Users/tkuehne/gxtb-local-build/install/tblite/bin/tblite \
  --cp2k-source /Users/tkuehne/gxtb-local-build/cp2k \
  --tblite-source /Users/tkuehne/gxtb-local-build/tblite \
  --jobs 10 --threads 1 --eos-mesh k444 \
  --energy-mesh k333 --energy-mesh k444 --energy-mesh k555 \
  --result-mesh k444

python3 Goldzak12/scripts/compare_goldzak12_results.py --mesh k444
python3 Goldzak12/scripts/plot_literature_comparison.py
python3 Goldzak12/scripts/validate_goldzak12_results.py \
  --eos-mesh k444 --energy-mesh k333 --energy-mesh k444 \
  --energy-mesh k555 --result-mesh k444
```

Raw calculations and generated inputs are kept below `Goldzak12/runs` and
`Goldzak12/inputs` and are ignored by Git. Curated CSV/Markdown tables,
provenance, and publication-ready PNG/PDF figures are versioned in
`Goldzak12/data` and `Goldzak12/figures`.

Literature comparison
---------------------

`scripts/plot_literature_comparison.py` augments the current EOS results with
published DFT and post-HF values. All errors are recomputed against the same
zero-point-corrected experimental values from Goldzak et al.; experimental
columns from other sources are retained only for provenance.

The comparison includes SCAN, SCAN-L, r2SCAN, and r2SCAN-L for 11 common
solids from Mejia-Rodriguez and Trickey (2020), and LSDA, PBE, PBEsol, TPSS,
revTPSS, TM, HSE06, and optB86b-vdW data from Mo et al. (2017). Coverage is
reported explicitly because the literature subsets differ.
