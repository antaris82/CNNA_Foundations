# SUMMARY — shared edge link-cycle operator decision gate

Model label:
CNNA growing primal simplicial complex with deterministic sequential provenance growth, directed antisymmetric birth-transport operators, signed-Jlock two-pair assembly motifs, and a shared-edge link-cycle diagnostic.

This package is intentionally a **decision test**, not another open-ended move to a larger local environment.  It tests whether the shared-edge link around an A/B assembly edge provides the missing local orientation/alignment operator.

Positive decision gate requires all of the following:

```text
1. nonzero signed circulation on the shared-edge link,
2. circulation flips under sibling birth-order kappa mirror,
3. strict_sym remains null,
4. used_delta_beta remains false,
5. effective motif Q/P J-lock residual < 0.2,
6. projected J^2 + I residual < 0.25.
```

If the link-cycle still gives residuals around 0.4--0.6, the documented interpretation is not "try the next bigger local structure", but: the current local line supports Q/P carrier structure without deriving a good local J-orientation on the actual Q/P motif space.

| option | variant | beta | pairs | asm | pair harm | Q harm | P harm | pair J-lock | signed | union lock | edge-if lock | link-cycle best | active link lock | link J2 | active J2 | nonzero circ | kappa flips | decision pass | abort marker | used dBeta? |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| A_to_B_rescan_strong_reuseB | real_growth | (1,0,4,0) | 4 | 1 | 0.381795 | 0.381795 | 0.248801 | 0.448553 | -0.624664 | 0.543826 | 0.596246 | 0.543826 | 0.5534 | 0.339506 | 0.339946 | 72 | 0 | 0 | True | False |
| A_to_B_rescan_strong_reuseB | strict_symmetrized_control | (1,0,0,0) | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | True | False |
| A_to_B_rescan_strong_reuseB | no_backreaction | (1,0,8,0) | 4 | 2 | 0.405623 | 0.405623 | 0.343367 | 0.339477 | -0.285265 | 0.661108 | 0.689923 | 0.661108 | 0.666166 | 0.471393 | 0.568405 | 144 | 0 | 0 | True | False |

## Best non-strict row

```json
{
  "option": "A_to_B_rescan_strong_reuseB",
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
  "pair_local_J_lock": 0.4485530110320348,
  "signed_birth": -0.6246641039192367,
  "union_motif_lock": 0.5438258343030387,
  "edge_interface_lock": 0.596246290239973,
  "link_cycle_lock": 0.5438258343030387,
  "link_cycle_active_lock": 0.5533997685830272,
  "link_cycle_J2": 0.3395063127800071,
  "link_cycle_active_J2": 0.33994644433996885,
  "nonzero_circ": 72,
  "kappa_flip_count": 0,
  "decision_pass": 0,
  "local_negative_abort_marker": true,
  "best_order_mode": "geometric_angle",
  "best_block_mode": "identity",
  "best_carrier_mode": "motif_only",
  "used_delta_beta": false
}
```

## Local negative decision marker

```json
{"all_non_strict_variants_abort_marked": true}
```
