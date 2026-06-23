# RESULTS — kappa-permuted candidate Pareto gate

## Comparative table

| variant | matched | A both | beta2 audit both | C-lock pass | kappa flip pass | all pass | all+beta2 | best C lock | best C flip | used dBeta? |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| real_growth | 206 | 206 | 88 | 19 | 23 | 3 | 0 | 0.128002 | 0.328323 | False |
| strict_symmetrized_control | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | False |
| no_backreaction | 236 | 236 | 108 | 15 | 28 | 2 | 0 | 0.0667835 | 1 | False |

## Gate logic

A successful candidate would satisfy:

```text
A_gate in identity and kappa-reflected model,
C-eigen J-lock residual below threshold in both,
Q/P support in both,
signed_birth flips under kappa reflection,
optional beta2-opening audit in both.
```

`delta_beta2` is audit-only and is never used to form the candidate space or decide a move.
