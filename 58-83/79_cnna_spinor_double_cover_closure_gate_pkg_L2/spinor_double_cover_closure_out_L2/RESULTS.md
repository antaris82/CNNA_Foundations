# RESULTS — Spinor / double-cover closure gate

## Comparative table

| option | variant | beta | pairs | asm | Q harm | P harm | pair J-lock | motif lock | projective pass | strong pass | T²≈-αI | α | T⁴≈α²I | raw J²+I | QP mean | T²Q≈-Q | T²P≈-P | best Ω/g | used dβ? |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|---:|
| A_to_B_rescan_strong_reuseB | real_growth | (1,0,4,0) | 4 | 1 | 0.381795 | 0.248801 | 0.448553 | 0.543826 | 0 | 0 | 0.116759 | 0.704253 | 0.221012 | 0.339506 | 3.17298e-13 | 0.116759 | 0.116759 | union_pair_exchange_skew / union_C_square_metric | False |
| A_to_B_rescan_strong_reuseB | strict_symmetrized_control | (1,0,0,0) | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |  /  | False |
| A_to_B_rescan_strong_reuseB | no_backreaction | (1,0,8,0) | 4 | 2 | 0.405623 | 0.343367 | 0.339477 | 0.661108 | 26 | 0 | 0.0866198 | 0.538057 | 0.167829 | 0.471393 | 1.22678e-13 | 0.0866198 | 0.0866198 | union_pair_exchange_skew / union_C_square_metric | False |

## Gate definition

A primary projective spinor row passes only if:

```text
α > min_spinor_alpha,
T² is close to -α I,
T⁴ is close to α² I,
T² maps Q and P back to their own subspaces with negative sign,
T maps Q/P subspaces into each other sufficiently well.
```

This distinguishes three situations:

1. ordinary failed J: Q/P may exist, but T² is not scalar negative;
2. projective spinor-like closure: T²≈-αI and T⁴≈α²I, even if α≠1;
3. strong complex/spinor closure: T²≈-I and T⁴≈I.

No normalization or polar correction is used to force closure.
