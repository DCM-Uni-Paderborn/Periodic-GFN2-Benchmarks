# Reference data

Primary source:

J. Klimes, D. R. Bowler, and A. Michaelides, "Van der Waals density
functionals applied to solids," Phys. Rev. B 83, 195131 (2011),
https://doi.org/10.1103/PhysRevB.83.195131; preprint:
https://arxiv.org/abs/1102.1358.

The CSV files transcribe the complete 23-system tables from the paper:

- `lattice_constants.csv`: Table I, including ZPE-corrected experiment.
- `bulk_moduli.csv`: Table II; the experimental bulk moduli are not ZPE
  corrected, consistently with the source.
- `cohesive_energies.csv`: Table III atomization energies of the solids,
  reported here as positive cohesive energies, including ZPE-corrected
  experiment.

Each table includes all eight complete published DFT series:
revPBE-vdW, rPW86-vdW2, optPBE-vdW, optB88-vdW, optB86b-vdW, LDA, PBEsol,
and PBE. The partial RPA values shown in the paper's figures are not included in
aggregate statistics because a complete, consistently tabulated 23-system RPA
series is not provided there.

The transcribed values reproduce the paper's rounded aggregate statistics, for
example the PBEsol lattice-constant MAE of 0.033 A, the PBEsol bulk-modulus
MARE of 4.8%, and the optB88-vdW cohesive-energy MAE of 0.07 eV/atom.

The original computational protocol uses at least seven points around the
energy minimum, a Murnaghan EOS, conventional cells, `8x8x8` k-points for
semiconductors and ionic solids, and `16x16x16` k-points for metals.
