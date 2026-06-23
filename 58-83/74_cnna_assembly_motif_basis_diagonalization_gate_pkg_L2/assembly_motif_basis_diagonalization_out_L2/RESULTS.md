# RESULTS — assembly motif basis diagonalization gate

## Comparative table

| option | variant | beta | pairs | asm | pair harm | Q harm | P harm | pair J-lock | signed | motif n | direct motif lock | union motif lock | direct pass | union pass | used dBeta? |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| A_to_B_rescan_strong_reuseB | real_growth | (1,0,4,0) | 4 | 1 | 0.381795 | 0.381795 | 0.248801 | 0.448553 | -0.624664 | 1 | 0.543826 | 0.543826 | 0 | 0 | False |
| A_to_B_rescan_strong_reuseB | strict_symmetrized_control | (1,0,0,0) | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | False |
| A_to_B_rescan_strong_reuseB | no_backreaction | (1,0,8,0) | 4 | 2 | 0.405623 | 0.405623 | 0.343367 | 0.339477 | -0.285265 | 2 | 0.661108 | 0.661108 | 0 | 0 | False |

## Gate criterion

The motif basis gate passes only if the combined A/B motif, not merely an individual pair, shows:

```text
J(span(A_Q,B_Q)) approximately subset span(A_P,B_P),
J(span(A_P,B_P)) approximately subset span(A_Q,B_Q),
low J-span leakage,
projected J^2 approximately -I,
strict_sym killed,
used_delta_beta = false.
```

This is a subspace-lock diagnostic.  It permits derived A/B mixing inside the motif planes, but it does not fit an arbitrary rotation and does not define P as JQ.
