# RESULTS — signed J-lock role-coupling gate

## Comparative table

| option | variant | beta | pairs | asm | pair harm | Q harm | P harm | J-lock | signed birth | selected C | selected kappa | used dBeta? |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| A_to_B_rescan_connected_reuseB | real_growth | (1,0,4,0) | 4 | 1 | 0.381795 | 0.381795 | 0.248801 | 0.448553 | -0.624664 | 0.491433 | 1 | False |
| A_to_B_rescan_connected_reuseB | strict_symmetrized_control | (1,0,0,0) | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | False |
| A_to_B_rescan_connected_reuseB | no_backreaction | (1,0,8,0) | 4 | 2 | 0.405623 | 0.405623 | 0.343367 | 0.339477 | -0.285265 | 0.339477 | 0.506181 | False |
| A_to_B_rescan_connected_noReuseB | real_growth | (1,0,4,0) | 4 | 1 | 0.381795 | 0.381795 | 0.248801 | 0.448553 | -0.624664 | 0.491433 | 1 | False |
| A_to_B_rescan_connected_noReuseB | strict_symmetrized_control | (1,0,0,0) | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | False |
| A_to_B_rescan_connected_noReuseB | no_backreaction | (1,0,8,0) | 4 | 2 | 0.405623 | 0.405623 | 0.343367 | 0.339477 | -0.285265 | 0.339477 | 0.506181 | False |
| A_to_B_rescan_strong_reuseB | real_growth | (1,0,4,0) | 4 | 1 | 0.381795 | 0.381795 | 0.248801 | 0.448553 | -0.624664 | 0.491433 | 1 | False |
| A_to_B_rescan_strong_reuseB | strict_symmetrized_control | (1,0,0,0) | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | False |
| A_to_B_rescan_strong_reuseB | no_backreaction | (1,0,8,0) | 4 | 2 | 0.405623 | 0.405623 | 0.343367 | 0.339477 | -0.285265 | 0.339477 | 0.506181 | False |
| A_to_B_rescan_strong_noReuseB | real_growth | (1,0,4,0) | 4 | 1 | 0.381795 | 0.381795 | 0.248801 | 0.448553 | -0.624664 | 0.491433 | 1 | False |
| A_to_B_rescan_strong_noReuseB | strict_symmetrized_control | (1,0,0,0) | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | False |
| A_to_B_rescan_strong_noReuseB | no_backreaction | (1,0,8,0) | 4 | 2 | 0.405623 | 0.405623 | 0.343367 | 0.339477 | -0.285265 | 0.339477 | 0.506181 | False |

## Gate criterion

A constructive result would require the same non-strict row to show:

```text
beta2 open,
Q/P harmonic positive,
signed_birth strong,
J-lock low,
strict_sym killed,
decision_used_delta_beta_any = false.
```

The previous order/context ablation showed that signed orientation and J-lock separate.  This test asks whether explicitly coupling those two quantities in the B-role can make them coincide dynamically.
