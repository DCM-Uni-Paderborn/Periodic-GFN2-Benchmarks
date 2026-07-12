LC12 (Goldzak12) native-Bloch CP2K/tblite benchmark
===================================================

This directory contains the 12 cubic covalent and ionic solids studied by
Goldzak, Wang, Ye, and Berkelbach, J. Chem. Phys. 157, 174112 (2022). It
compares CP2K/tblite GFN1-xTB and GFN2-xTB with the reported HF, MP2, SCS-MP2,
SOS-MP2, and zero-point-corrected experimental lattice constants and cohesive
energies.

Current production run
----------------------

The 2026-07-12 rerun used:

- CP2K trunk revision `faf9aae91266170dfee8a9f7171a5135bc5eb368` with the
  local CP2K/tblite interface patch recorded by hash in
  `data/build_provenance.json`;
- tblite revision `8a9d09474b93d25c044d6f46ce920750c7fe4cf7`, which combines current
  `main` with PR #350 and includes the previously merged PR #343;
- conventional cubic eight-atom cells and CP2K native Bloch sampling through
  `&KPOINTS` with `SCHEME MACDONALD`, `SYMMETRY T`, and `FULL_GRID F`, using
  full SPGLIB symmetry reduction; no Born-von-Karman supercells;
- a `k444` cubic equation of state, `k333/k444/k555` final single points, and
  `k555` as the reported cohesive-energy mesh;
- CP2K energies extrapolated to electronic temperature `T -> 0`;
- matching tblite CLI isolated-atom references with explicit atomic spins and
  `ACCURACY 0.05`.

Only completed SCF points enter an EOS fit. A quadratic fit is rejected when
its local RMSE exceeds 0.02 hartree or when its fitted minimum lies more than
0.02 hartree above the sampled local minimum. The current run gives 12/12
valid GFN1 fits and 10/12 valid GFN2 fits. GFN2/MgO has no bracketed stable
minimum on the compressed branch, while GFN2/LiH has a discontinuous EOS and
fails the general fit-quality criterion.

Volume continuation with separate CP2K Bloch-wavefunction and native tblite
SCC restarts removes the earlier independent-start failures. LiH converges for
all 32 sampled points down to scale 0.71, but its energy continues to decrease
until the electronic branch collapses below that range. MgO can be followed in
fine steps to scale 0.926, where its energy is still decreasing; the additional
0.90 and 0.88 points enter the same charge-collapse branch even after damped
2400-step SCC retries. The 10/12 coverage therefore reflects missing physical
EOS minima, not unfinished production jobs.

Current versus previous results
-------------------------------

| method | coverage | lattice MAE (A) | cohesive-energy MAE (eV/atom) |
|---|---:|---:|---:|
| GFN1 current | 12/12 | 0.136650 | 1.457694 |
| GFN1 previous | 12/12 | 0.164341 | 1.457325 |
| GFN2 current | 10/12 | 0.062410 | 1.299325 |
| GFN2 previous | 11/12 | 0.147638 | 1.731839 |

On the identical ten-system GFN2 subset, the lattice-constant MAE decreases
from 0.133264 to 0.062410 A and the cohesive-energy MAE decreases from 1.534521
to 1.299325 eV/atom. The frozen previous tables are in
`data/baseline_20260710`; `data/old_vs_new.md` and the associated CSV files
contain the complete per-system comparison.

Reproduction
------------

Run from the repository root:

```bash
python3 Goldzak12/scripts/run_goldzak12_eos_benchmark.py \
  --cp2k /path/to/cp2k.ssmp \
  --tblite /path/to/tblite \
  --cp2k-source /path/to/cp2k \
  --tblite-source /path/to/tblite \
  --jobs 10 --threads 1 --eos-mesh k444 \
  --energy-mesh k333 --energy-mesh k444 --energy-mesh k555 \
  --result-mesh k555

python3 Goldzak12/scripts/compare_goldzak12_results.py --mesh k555
python3 Goldzak12/scripts/plot_literature_comparison.py
python3 Goldzak12/scripts/validate_goldzak12_results.py \
  --eos-mesh k444 --energy-mesh k333 --energy-mesh k444 \
  --energy-mesh k555 --result-mesh k555

python3 Goldzak12/scripts/continue_goldzak12_eos.py \
  --solid LiH --method GFN2 --mesh k444 --start-scale 0.94 \
  --scale 0.93 --scale 0.92 --variant lih_scc_continuation \
  --mixer tblite --memory 250 --damping 0.4 --promote
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
