# SUMMARY — signed J-lock role-coupling gate

Model label:
CNNA growing primal simplicial complex with deterministic sequential provenance growth, directed antisymmetric birth-transport operators, and a dynamic two-pair assembly rule whose B-role couples C-lock with signed kappa/birth amplitude.

This test keeps the A-role as provenance/QP-carrier proxy and changes only the B-role scoring:

```text
B role = low C-eigen J-lock residual + signed birth/kappa amplitude + kappa flip,
         with shared face/edge context to A.
```

No decision uses delta_beta, H2, harmonic projections, complex scalars, Hodge, positivity, physical adjoint, or final sym(M).

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

## Best rows

```json
{
  "best_signed": {
    "option": "A_to_B_rescan_connected_reuseB",
    "variant": "real_growth",
    "beta": [
      1,
      0,
      4,
      0
    ],
    "pairings": 4,
    "assemblies": 1,
    "pair_harm": 0.38179479579105513,
    "Q_harm": 0.38179479579105513,
    "P_harm": 0.24880104631557917,
    "J_lock": 0.4485530110320348,
    "signed_birth": -0.6246641039192367,
    "selected_C_lock_avg": 0.49143315688516676,
    "selected_kappa_flip_avg": 0.9999999999980294,
    "used_delta_beta": false
  },
  "best_J_lock": {
    "option": "A_to_B_rescan_connected_reuseB",
    "variant": "no_backreaction",
    "beta": [
      1,
      0,
      8,
      0
    ],
    "pairings": 4,
    "assemblies": 2,
    "pair_harm": 0.40562305647622443,
    "Q_harm": 0.40562305647622443,
    "P_harm": 0.3433668689249218,
    "J_lock": 0.3394769439901298,
    "signed_birth": -0.2852651970622306,
    "selected_C_lock_avg": 0.3394769439901298,
    "selected_kappa_flip_avg": 0.5061807642415207,
    "used_delta_beta": false
  },
  "best_signed_over_J": {
    "option": "A_to_B_rescan_connected_reuseB",
    "variant": "real_growth",
    "beta": [
      1,
      0,
      4,
      0
    ],
    "pairings": 4,
    "assemblies": 1,
    "pair_harm": 0.38179479579105513,
    "Q_harm": 0.38179479579105513,
    "P_harm": 0.24880104631557917,
    "J_lock": 0.4485530110320348,
    "signed_birth": -0.6246641039192367,
    "selected_C_lock_avg": 0.49143315688516676,
    "selected_kappa_flip_avg": 0.9999999999980294,
    "used_delta_beta": false
  }
}
```
