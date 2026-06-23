# SUMMARY — real symplectic before star gate

Model label:
CNNA growing primal simplicial complex with deterministic sequential provenance growth, directed antisymmetric birth-transport operators, signed-Jlock two-pair assembly motifs, and a real symplectic-first diagnostic.

Corrected ladder order:

```text
real Q/P carrier
→ real symplectic form Ω
→ possible derived #/*-structure
→ only then J / complex orientation
```

This package therefore does **not** use `J²≈-I` as the primary gate.  It asks first whether any already-derived skew operator candidate restricts to a nondegenerate real two-form on the actual Q/P or A/B motif space.

Primary Ω candidates:

```text
pair_exchange / union_pair_exchange
edge_interface_only
union_plus_edge_interface
link_cycle_only
union_plus_link_cycle
```

The `data_wedge_control` is logged only as a tautological upper-bound control because it uses Q/P data directly.

Positive primary symplectic gate requires:

```text
projected Ω skew residual <= 1e-08
projected Ω full rank on an even-dimensional Q/P span
nondegeneracy ratio >= 0.001
Q and P have equal positive rank
Ω(Q,P) is full rank with ratio >= 0.001
Q and P are approximately Ω-isotropic <= 0.35
strict_sym remains null
used_delta_beta remains false
```

| option | variant | beta | pairs | asm | pair harm | Q harm | P harm | pair J-lock | union J-lock | signed | primary symp pass | motif symp pass | best primary Ω | level | dim/rank | Ω nondeg ratio | Ω QP ratio | used dβ? |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|---|---:|---:|---:|---:|
| A_to_B_rescan_strong_reuseB | real_growth | (1,0,4,0) | 4 | 1 | 0.381795 | 0.381795 | 0.248801 | 0.448553 | 0.543826 | -0.624664 | 40 | 40 | union_plus_link_cycle_geometric_angle_identity_unit_(0, 2) | motif | 4/4 | 0.785587 | 0.785587 | False |
| A_to_B_rescan_strong_reuseB | strict_symmetrized_control | (1,0,0,0) | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |  |  | 0/0 | 0 | 0 | False |
| A_to_B_rescan_strong_reuseB | no_backreaction | (1,0,8,0) | 4 | 2 | 0.405623 | 0.405623 | 0.343367 | 0.339477 | 0.661108 | -0.285265 | 80 | 80 | union_plus_link_cycle_geometric_angle_identity_unit_(4, 7) | motif | 4/4 | 0.83832 | 0.83832 | False |

## High-level decision

```json
{
  "any_non_strict_primary_symplectic_pass": true,
  "any_non_strict_primary_motif_symplectic_pass": true
}
```
