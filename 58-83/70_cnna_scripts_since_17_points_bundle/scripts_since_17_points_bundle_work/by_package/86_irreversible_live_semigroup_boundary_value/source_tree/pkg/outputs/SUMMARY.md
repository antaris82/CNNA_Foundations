# CNNA irreversible live semigroup / boundary-value gate

Package: `cnna_irreversible_live_semigroup_boundary_value_gate_pkg_L2`

## Question

Does the true Schur/DtN live evolution after a birth event behave as a directed, bounded, retarded semigroup layer?

This test no longer asks whether the newborn can be deleted.  It checks whether the post-birth live Schur/DtN states form a forward-oriented record-to-live evolution with bounded contraction, nonzero record/live gap, and boundary-value polarity.

## Method

The established script-1/script-2 ternary growth model is kept:

- sequential sibling births;
- incoming environment response to the newborn;
- newborn backreaction to parent line and older siblings;
- true real Laplace Schur/DtN maps on boundary ports;
- fixed-topology live relaxation between births.

For every birth and every fixed cut, the DtN matrix is converted to a real invariant vector at relax steps `0..N`.  A feature-space forward map `T_forward` and reverse map `T_reverse` are fit as diagnostics only.  The gate is based on bounded forward evolution, decaying DtN increments, nonzero record/live gap, and stable oriented boundary polarity.  The forward feature-map norm may be slightly above 1; the contraction evidence is therefore the decay of DtN increments, not a strict norm < 1 in every feature direction.  No J, i, Hodge, physical star, positivity, C*-norm, Q/P target, or delta-beta is used.

## Summary table

| variant | semigroup gate | fwd bounded/quasi-contracting | reverse singular | transitions | fwd norm | median last/first | record/live gap | halfplane bias |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| real_growth_linear_true_schur_live_semigroup | True | True | True | 3276 | 1.04 | 0.0432 | 0.032 | 1 |
| log_growth_true_schur_live_semigroup | True | True | True | 3276 | 1.05 | 0.0421 | 0.0175 | 1 |
| saturating_growth_true_schur_live_semigroup | True | True | True | 3276 | 1.09 | 0.101 | 0.278 | 1 |
| kappa_reversed_birth_order_true_schur_live_semigroup_control | True | True | True | 3276 | 1.04 | 0.0432 | 0.032 | 1 |
| longitudinal_axis_flip_live_semigroup_control | True | True | True | 3276 | 1.04 | 0.0432 | 0.032 | 1 |
| strict_symmetrized_response_control_true_schur_live_semigroup | False | True | True | 3276 | 1 | 0 | 0 | 0 |
| birth_only_no_relax_record_control_true_schur_semigroup | False | False | False | 0 | 0 | 0 | 0 | 0 |

## Interpretation

A positive semigroup gate means, conservatively:

```text
record DtN state  ->  live DtN state
```

is directed, bounded, and quasi-contracting in the live relaxation layer: the Schur/DtN increments decay, the record/live gap remains nonzero, and the reverse map is singular/ill-conditioned.  It is not a proof of a unique continuum semigroup generator.  It does not mean that a complex structure, J, spin, a *-operation, or Q/P compatibility has been derived.

A negative strict-sym control means that topology alone is not enough: the response/backreaction layer is required.

The longitudinal-axis-flip control uses the same live evolution but reverses the interpretation of the root-front polarity.  This is not a new physical dynamics; it audits whether the boundary-value sign is relative to the chosen longitudinal orientation.

## Next test

`test_live_semigroup_polarity_to_operator_bridge_gate.py`

Use the directed live-semigroup polarity as an input to an operator bridge, but keep it separate from J.  The next gate should ask whether this record-to-live semigroup polarity can induce a candidate real adjunction or passivity/positivity surrogate before any complex structure is asserted.

## Additional audit note

The kappa-reversed birth-label control stays positive for this gate.  That means the live semigroup polarity is not primarily a sibling-label chirality.  It is carried by the record-to-live response direction.

The longitudinal-axis-flip control also has absolute halfplane bias `1`, but the sign fractions swap: the same live process is being read with the opposite root/front orientation.  Therefore this control confirms relative boundary polarity rather than a new dynamics.

The strict-sym control can have a formally bounded feature map because its feature rows are static/degenerate, but it has zero record/live gap and zero halfplane bias, so the semigroup boundary gate is false.
