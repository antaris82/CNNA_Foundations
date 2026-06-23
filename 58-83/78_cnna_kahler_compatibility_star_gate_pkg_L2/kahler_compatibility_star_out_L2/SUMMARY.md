# SUMMARY — Kähler compatibility / star-from-symplectic gate

Model tag: `CQNM/s=-1 saturated geometry reference, provenance-growth L2 diagnostic`.

This package corrects the previous symplectic test: a nondegenerate skew form alone is not counted as a breakthrough.  The gate now asks whether a primary Ω candidate and a primary g candidate are simultaneously compatible enough that

```text
J = g^-1 Ω
```

has small `J²+I`, maps the derived Q-subspace to the derived P-subspace, is metric-compatible, and yields a stable metric-adjoint `#` on the tested operator family.

| option | variant | beta | pairs | asm | pair harm | Q harm | P harm | symp pass | best Ω | Ω ratio | compat pass | valid Ω-g | best Ω/g | J2 resid | QP lock | metric orth | #J anti | star span | used dβ? |
|---|---|---:|---:|---:|---:|---:|---:|---:|---|---:|---:|---:|---|---:|---:|---:|---:|---:|---:|
| A_to_B_rescan_strong_reuseB | real_growth | (1,0,4,0) | 4 | 1 | 0.381795 | 0.381795 | 0.248801 | 40 | union_plus_link_cycle_geometric_angle_identity_unit_(0, 2) | 0.785587 | 0 | 156 | union_pair_exchange_skew / cochain_identity_metric | 0.339506 | 3.17297e-13 | 0.339506 | 2.54945e-16 | 6.4453e-16 | False |
| A_to_B_rescan_strong_reuseB | strict_symmetrized_control | (1,0,0,0) | 0 | 0 | 0 | 0 | 0 | 0 |  | 0 | 0 | 0 |  /  | 0 | 0 | 0 | 0 | 0 | False |
| A_to_B_rescan_strong_reuseB | no_backreaction | (1,0,8,0) | 4 | 2 | 0.405623 | 0.405623 | 0.343367 | 80 | union_plus_link_cycle_geometric_angle_identity_unit_(4, 7) | 0.83832 | 0 | 312 | union_pair_exchange_skew / cochain_identity_metric | 0.471393 | 1.22677e-13 | 0.471393 | 2.92083e-16 | 3.82069e-16 | False |

Decision:

```json
{
  "any_non_strict_primary_compatibility_pass": false
}
```
