# RESULTS — dual assembly order/context ablation gate

## Comparative table

| option | variant | beta | pairs | asm | pair harm | Q harm | P harm | J-lock | signed birth | used dBeta? |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| A_to_B_rescan_connected_reuseB | real_growth | (1,0,4,0) | 4 | 1 | 0.381795 | 0.381795 | 0.248801 | 0.448553 | -0.624664 | False |
| A_to_B_rescan_connected_reuseB | strict_symmetrized_control | (1,0,0,0) | 0 | 0 | 0 | 0 | 0 | 0 | 0 | False |
| A_to_B_rescan_connected_reuseB | no_backreaction | (1,0,6,0) | 4 | 2 | 0.302564 | 0.302564 | 0.304787 | 0.208886 | 0.0704225 | False |
| A_to_B_rescan_connected_noReuseB | real_growth | (1,0,4,0) | 4 | 1 | 0.381795 | 0.381795 | 0.248801 | 0.448553 | -0.624664 | False |
| A_to_B_rescan_connected_noReuseB | strict_symmetrized_control | (1,0,0,0) | 0 | 0 | 0 | 0 | 0 | 0 | 0 | False |
| A_to_B_rescan_connected_noReuseB | no_backreaction | (1,0,6,0) | 4 | 2 | 0.302564 | 0.302564 | 0.304787 | 0.208886 | 0.0704225 | False |
| B_to_A_rescan_connected_reuseB | real_growth | (1,0,4,0) | 4 | 1 | 0.339534 | 0.339534 | 0.284669 | 0.371324 | -0.618486 | False |
| B_to_A_rescan_connected_reuseB | strict_symmetrized_control | (1,0,0,0) | 0 | 0 | 0 | 0 | 0 | 0 | 0 | False |
| B_to_A_rescan_connected_reuseB | no_backreaction | (1,0,4,0) | 4 | 1 | 0.327035 | 0.327035 | 0.26254 | 0.334279 | -0.562168 | False |
| B_to_A_rescan_connected_noReuseB | real_growth | (1,0,4,0) | 4 | 1 | 0.339534 | 0.339534 | 0.284669 | 0.371324 | -0.618486 | False |
| B_to_A_rescan_connected_noReuseB | strict_symmetrized_control | (1,0,0,0) | 0 | 0 | 0 | 0 | 0 | 0 | 0 | False |
| B_to_A_rescan_connected_noReuseB | no_backreaction | (1,0,4,0) | 4 | 1 | 0.327035 | 0.327035 | 0.26254 | 0.334279 | -0.562168 | False |
| stale_same_scan_connected_reuseB | real_growth | (1,0,6,0) | 4 | 1 | 0.25221 | 0.25221 | 0.333058 | 0.325883 | -0.241618 | False |
| stale_same_scan_connected_reuseB | strict_symmetrized_control | (1,0,0,0) | 0 | 0 | 0 | 0 | 0 | 0 | 0 | False |
| stale_same_scan_connected_reuseB | no_backreaction | (1,0,4,0) | 4 | 0 | 0.345443 | 0.345443 | 0.293827 | 0.334831 | -0.410279 | False |
| stale_same_scan_connected_noReuseB | real_growth | (1,0,6,0) | 4 | 1 | 0.25221 | 0.25221 | 0.333058 | 0.325883 | -0.241618 | False |
| stale_same_scan_connected_noReuseB | strict_symmetrized_control | (1,0,0,0) | 0 | 0 | 0 | 0 | 0 | 0 | 0 | False |
| stale_same_scan_connected_noReuseB | no_backreaction | (1,0,4,0) | 4 | 0 | 0.345443 | 0.345443 | 0.293827 | 0.334831 | -0.410279 | False |

## Interpretation protocol

A robust constructive result would show one option with:

```text
strict_sym killed,
beta2 open,
Q/P harmonic positive,
signed_birth strong,
J-lock low,
decision_used_delta_beta_any false.
```

The ablation asks whether the signed/J-lock tension seen in the dynamic two-pair assembly is caused by order, stale legality, context restriction, or face reuse.
