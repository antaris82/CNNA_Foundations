# RESULTS — Kähler compatibility / star-from-symplectic gate

## Comparative table

| option | variant | beta | pairs | asm | pair harm | Q harm | P harm | symp pass | best Ω | Ω ratio | compat pass | valid Ω-g | best Ω/g | J2 resid | QP lock | metric orth | #J anti | star span | used dβ? |
|---|---|---:|---:|---:|---:|---:|---:|---:|---|---:|---:|---:|---|---:|---:|---:|---:|---:|---:|
| A_to_B_rescan_strong_reuseB | real_growth | (1,0,4,0) | 4 | 1 | 0.381795 | 0.381795 | 0.248801 | 40 | union_plus_link_cycle_geometric_angle_identity_unit_(0, 2) | 0.785587 | 0 | 156 | union_pair_exchange_skew / cochain_identity_metric | 0.339506 | 3.17297e-13 | 0.339506 | 2.54945e-16 | 6.4453e-16 | False |
| A_to_B_rescan_strong_reuseB | strict_symmetrized_control | (1,0,0,0) | 0 | 0 | 0 | 0 | 0 | 0 |  | 0 | 0 | 0 |  /  | 0 | 0 | 0 | 0 | 0 | False |
| A_to_B_rescan_strong_reuseB | no_backreaction | (1,0,8,0) | 4 | 2 | 0.405623 | 0.405623 | 0.343367 | 80 | union_plus_link_cycle_geometric_angle_identity_unit_(4, 7) | 0.83832 | 0 | 312 | union_pair_exchange_skew / cochain_identity_metric | 0.471393 | 1.22677e-13 | 0.471393 | 2.92083e-16 | 3.82069e-16 | False |

## Gate definition

A row only passes if all of the following hold:

```text
Ω is skew and nondegenerate on the actual Q/P motif span,
g is symmetric and nondegenerate on the same span,
J = g^-1 Ω satisfies J² ≈ -I,
J maps span(Q_motif) to span(P_motif),
J is g-orthogonal enough,
J# ≈ -J for the metric-adjoint #,
the tested operator family is # stable,
strict_sym remains null,
used_delta_beta remains false.
```

The test deliberately logs, but does not reward, the tautological facts that a metric adjoint is involutive and anti-multiplicative for invertible g.  The success condition is simultaneous compatibility of Ω, g, #, and the Q/P split.

## Interpretation rule

- If `compat pass > 0`, the ladder may proceed from real Ω to a candidate real #/*-algebra and then to J.
- If Ω passes but compatibility fails, the result is a Kähler-like compatibility obstruction: the ingredients exist separately but do not align.
- If strict_sym is nonzero, the result is invalid as a provenance asymmetry diagnostic.
