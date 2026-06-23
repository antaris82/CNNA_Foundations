# CNNA live-semigroup polarity to operator bridge gate

Package: `cnna_live_semigroup_polarity_to_operator_bridge_gate_pkg_L2`

## Question

Can the directed Record->Live Schur/DtN live-semigroup polarity prepare a real operator bridge before any J, i, spin, Hodge, physical star, Q/P target, or C*-positivity is asserted?

This is a deliberately weak but simultaneous gate.  It does not count a linear-algebraic adjoint by itself.  It asks whether the following appear together:

```text
nontrivial Record->Live DtN drift
bounded record DtN metric on the zero-sum boundary-voltage subspace
metric self-adjointness of the drift operator A = G^+ Delta
one-sided dissipative / passivity-surrogate spectrum of signed Delta
stable longitudinal polarity of the relaxation increments
bounded decay of live Schur/DtN increments
```

## Method

For each birth event and fixed boundary cut:

```text
G     = zero-sum projection of Lambda_record
Delta = zero-sum projection of Lambda_live_final - Lambda_record
ell   = derived cut-vs-UV longitudinal boundary mode
A     = G^+ Delta
```

Then the test records:

- `G_adjoint_residual` for A;
- one-sided signed spectrum of `axis_sign * Delta`;
- longitudinal energy sign `ell^T Delta ell`;
- decay of live increments `Delta Lambda_k`;
- strict-sym and birth-only controls.

This is not a physical positivity theorem.  It is a bridge audit from live-boundary polarity toward possible later real adjunction/passivity structure.

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

## Interpretation guide

- `bridge pass` is the simultaneous weak operator-bridge gate.
- `passivity pass` means signed live drift has one-sided spectral tendency and stable longitudinal increment polarity.
- `selfadjoint` alone is expected because the construction uses real Schur/DtN symmetric matrices; it is logged but not a success criterion by itself.
- A positive `strict_sym` result would be suspicious.  The expected result is zero rows/pass there.

## Next test

If the bridge gate is positive, the next test should examine whether the bridge induces a stable real `#` candidate on an operator family.  If only passivity/one-sided polarity is positive but bridge pass is weak, the next step is to strengthen the boundary-value/semigroup layer before claiming any operator structure.
