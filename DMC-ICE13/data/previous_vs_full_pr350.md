# DMC-ICE13: previous manuscript stack versus final corrected stack

The previous manuscript data used CP2K `518a509` and tblite `5b14b843`. The
final rerun uses CP2K trunk `faf9aae` and tblite merge `8a9d094`, containing
the complete PR 350 head `8c5e562` on top of current `main`, together with the
smooth Wigner--Seitz image weighting and analytic force/stress derivatives.
All k-point calculations use native Bloch wavefunctions and full SPGLIB
symmetry reduction.

All values below are relative-energy MAEs in kJ mol-1 per water molecule over
the twelve non-Ih phases.

| mesh | GFN1 previous | GFN1 current | change | GFN2 previous | GFN2 current | change |
|---|---:|---:|---:|---:|---:|---:|
| Gamma | 6.696681 | 6.694624 | -0.002057 | 5.355715 | 5.578897 | +0.223182 |
| 2x2x2 | 7.959770 | 7.956838 | -0.002932 | 3.233027 | 3.510100 | +0.277073 |
| 3x3x3 | 8.008187 | 8.005255 | -0.002932 | 3.185301 | 3.462919 | +0.277618 |
| 4x4x4 | 8.009427 | 8.006494 | -0.002933 | 3.183780 | 3.461424 | +0.277644 |
| 5x5x5 | 8.009417 | 8.006485 | -0.002932 | 3.183706 | 3.461353 | +0.277647 |

Negative changes are improvements. Relative to the manuscript values, GFN1 is
unchanged within 0.003 kJ mol-1, whereas GFN2 is higher by 0.278 kJ mol-1 on
the primary 3x3x3 mesh. GFN2 nevertheless retains a 4.54 kJ mol-1 lower MAE
than GFN1 at that mesh. All 156 calculations completed with zero return codes.
