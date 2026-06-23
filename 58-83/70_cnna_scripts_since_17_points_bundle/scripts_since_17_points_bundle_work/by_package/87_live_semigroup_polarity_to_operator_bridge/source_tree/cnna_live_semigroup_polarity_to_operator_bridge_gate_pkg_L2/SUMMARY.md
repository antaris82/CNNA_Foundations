# CNNA live-semigroup polarity to operator bridge gate

Package: `cnna_live_semigroup_polarity_to_operator_bridge_gate_pkg_L2`

## Model status

- Growth model: established script-1/script-2 ternary sequential growth family.
- Boundary response: true real Laplace Schur/DtN matrices.
- Live layer: fixed-topology post-birth relaxation.
- Regime label: growing real complement network, true-Schur/DtN, Record/Live two-layer semantics.
- No `i`, no `J`, no spin target, no Hodge star, no physical `*`, no Hilbert positivity, no C*-norm, no Q/P target, no `delta_beta` decision.

This finite run uses `max_level=3`, `relax_steps=6`.  L4 in the same unoptimized true-Schur implementation is expensive because each birth/relaxation cut recomputes Schur complements on the growing full graph.

## Question

Can the directed Record→Live Schur/DtN semigroup polarity prepare a real operator bridge before any complex structure is asserted?

The gate is deliberately earlier than `J`:

```text
Birth/Record layer
→ Live Schur/DtN relaxation polarity
→ real boundary-operator bridge diagnostics
→ possible later # / passivity / positivity surrogate
→ only later J / i / Q/P compatibility
```

## Construction

For each birth event and each fixed boundary cut, the test computes a true DtN sequence on a fixed boundary port set:

```text
Lambda_0, Lambda_1, ..., Lambda_N
```

where `Lambda_0` is the birth-record DtN and `Lambda_N` is the live-relaxed DtN after fixed-topology relaxation.

The matrices are projected to the canonical zero-sum boundary-voltage subspace, because a Laplace DtN operator has the constant vector in its kernel.  Then:

```text
G     = projected Lambda_record
Delta = projected Lambda_live_final - projected Lambda_record
ell   = derived cut-vs-UV longitudinal boundary mode
A     = G^+ Delta
```

The test records:

- `G_adjoint_residual` for `A`, but does not count this alone as a success because such an adjunction is mostly linear algebra once `G` is fixed.
- One-sided signed spectrum of `axis_sign * Delta`.
- Longitudinal energy sign `ell^T Delta ell`.
- Decay of live increments `Delta Lambda_k`.
- Strict-sym and birth-only controls.

## Summary table

| variant | rows | bridge pass | passivity pass | bounded G | selfadjoint | decay | one-sided eig | incr sign bias | +long E | -long E |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| real_growth_linear_live_polarity_operator_bridge | 102 | 0.147 | 0.157 | 1.000 | 0.529 | 1.000 | 0.941 | 0.301 | 1.000 | 0.000 |
| log_growth_live_polarity_operator_bridge | 102 | 0.127 | 0.157 | 1.000 | 0.461 | 1.000 | 0.958 | 0.376 | 1.000 | 0.000 |
| saturating_growth_live_polarity_operator_bridge | 102 | 0.471 | 1.000 | 1.000 | 0.471 | 1.000 | 1.000 | 1.000 | 1.000 | 0.000 |
| kappa_reversed_birth_order_operator_bridge_control | 102 | 0.147 | 0.157 | 1.000 | 0.529 | 1.000 | 0.941 | 0.301 | 1.000 | 0.000 |
| longitudinal_axis_flip_operator_bridge_control | 102 | 0.147 | 0.157 | 1.000 | 0.529 | 1.000 | 0.941 | 0.301 | 0.000 | 1.000 |
| strict_symmetrized_response_operator_bridge_control | 102 | 0.000 | 0.000 | 1.000 | 1.000 | 1.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| birth_only_no_relax_operator_bridge_control | 0 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |

## Interpretation

The result is **partly positive**, but not a derived `*`-algebra and not a complex structure.

Positive findings:

```text
1. The Record→Live Schur/DtN drift is nonzero in the dynamic growth variants.
2. The projected record DtN metric G is bounded on all valid rows.
3. Live increments decay on all valid dynamic rows.
4. The signed Delta spectrum is strongly one-sided in the growth variants.
5. The longitudinal boundary-energy polarity is consistent.
6. Longitudinal axis flip reverses the interpreted sign.
7. strict_sym kills the nontrivial passivity/bridge signal.
8. birth_only_no_relax has no bridge rows, as expected.
```

The strongest case is `saturating_growth`, where every valid row passes the passivity surrogate and about 47% pass the full weak operator-bridge gate.  Linear/log growth already show strong one-sided spectral tendency, but weaker simultaneous bridge pass.

Negative / cautionary findings:

```text
1. Metric selfadjointness alone is not meaningful enough; it is partly expected from real symmetric Schur/DtN data.
2. The bridge is not robust across all rows in linear/log variants.
3. kappa-reversed sibling labels do not change the gate; the polarity is not primarily sibling-label chiral.
4. This does not yet define a physical star, positivity, C*-norm, J, i, or Q/P-compatible operator.
```

The cleanest reading is:

```text
The irreversible Live-Schur/DtN semigroup prepares a real boundary-operator polarity:
  bounded Record metric,
  decaying live increments,
  one-sided signed live drift,
  longitudinal boundary-energy orientation.

This is a plausible precursor to passivity/adjunction, but not yet an operator algebra.
```

## Methodological note

The `longitudinal_axis_flip` control does not change the dynamics. It flips the root/front interpretation of the longitudinal polarity. The result swaps `+long E` and `-long E`, which is exactly what one expects if the sign is relative to the chosen root/front axis rather than an absolute scalar sign.

## Next test

`test_live_semigroup_operator_family_star_candidate_gate.py`

Use only rows that pass the weak bridge/passivity gate and build a small real operator family from the Schur/DtN Record→Live drift operators.  Then test whether a candidate real `#` can be made stable **without** importing transpose as a physical star:

```text
## ≈ id
(AB)# ≈ B# A#
commutator compatibility
Record/Live polarity compatibility
strict_sym null
no J/i/QP target
```

This should be treated as an audit of a possible real adjunction layer, not as a C*-algebra claim.
