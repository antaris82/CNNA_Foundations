# SUMMARY — dual assembly order/context ablation gate

Model label:
CNNA growing primal simplicial complex with deterministic sequential provenance growth, directed antisymmetric birth-transport operators, and dynamic/stale two-pair assembly ablations.

This package tests whether the two-pair assembly result depends on:

```text
A->B after rescan
B->A after rescan
stale same-scan A/B
connected versus strong context
allowing B to reuse A faces
```

No decision uses delta_beta, H2, complex scalars, Hodge, positivity, a physical adjoint, or final sym(M).

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

## Best-option summary

```json
{
  "best_signed": {
    "option": "A_to_B_rescan_connected_reuseB",
    "variant": "real_growth",
    "beta2": 4,
    "pairings": 4,
    "assemblies": 1,
    "pair_harm": 0.38179479579105513,
    "Q_harm": 0.38179479579105513,
    "P_harm": 0.24880104631557917,
    "J_lock": 0.4485530110320348,
    "signed_birth": -0.6246641039192367,
    "used_delta_beta": false
  },
  "best_J_lock": {
    "option": "A_to_B_rescan_connected_reuseB",
    "variant": "no_backreaction",
    "beta2": 6,
    "pairings": 4,
    "assemblies": 2,
    "pair_harm": 0.30256353746579273,
    "Q_harm": 0.30256353746579273,
    "P_harm": 0.30478665264400706,
    "J_lock": 0.20888636621708123,
    "signed_birth": 0.07042246953279575,
    "used_delta_beta": false
  },
  "best_beta2": {
    "option": "A_to_B_rescan_connected_reuseB",
    "variant": "no_backreaction",
    "beta2": 6,
    "pairings": 4,
    "assemblies": 2,
    "pair_harm": 0.30256353746579273,
    "Q_harm": 0.30256353746579273,
    "P_harm": 0.30478665264400706,
    "J_lock": 0.20888636621708123,
    "signed_birth": 0.07042246953279575,
    "used_delta_beta": false
  },
  "best_pair_harm": {
    "option": "A_to_B_rescan_connected_reuseB",
    "variant": "real_growth",
    "beta2": 4,
    "pairings": 4,
    "assemblies": 1,
    "pair_harm": 0.38179479579105513,
    "Q_harm": 0.38179479579105513,
    "P_harm": 0.24880104631557917,
    "J_lock": 0.4485530110320348,
    "signed_birth": -0.6246641039192367,
    "used_delta_beta": false
  }
}
```
