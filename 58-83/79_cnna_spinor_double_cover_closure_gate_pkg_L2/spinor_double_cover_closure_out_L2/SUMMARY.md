# SUMMARY — Spinor / double-cover closure gate

Model tag: `CQNM/s=-1 saturated geometry reference, provenance-growth L2 diagnostic`.

This package tests the user's spin-1/2 suspicion in a strict derived-only form.  It does not set spin.  For each existing Ω/g-derived operator `T = g^-1 Ω` on the actual Q/P motif span it audits whether the operator is projectively order-four:

```text
T(Q) -> P,
T² ≈ -α I,
T⁴ ≈ α² I,
```

with α > 0.  This is the real double-cover/spinor-like signature: a sign after the second application and closure after the fourth application.  The strong unscaled gate additionally requires `T²≈-I` and `T⁴≈I`.

| option | variant | beta | pairs | asm | Q harm | P harm | pair J-lock | motif lock | projective pass | strong pass | T²≈-αI | α | T⁴≈α²I | raw J²+I | QP mean | T²Q≈-Q | T²P≈-P | best Ω/g | used dβ? |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|---:|
| A_to_B_rescan_strong_reuseB | real_growth | (1,0,4,0) | 4 | 1 | 0.381795 | 0.248801 | 0.448553 | 0.543826 | 0 | 0 | 0.116759 | 0.704253 | 0.221012 | 0.339506 | 3.17298e-13 | 0.116759 | 0.116759 | union_pair_exchange_skew / union_C_square_metric | False |
| A_to_B_rescan_strong_reuseB | strict_symmetrized_control | (1,0,0,0) | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |  /  | False |
| A_to_B_rescan_strong_reuseB | no_backreaction | (1,0,8,0) | 4 | 2 | 0.405623 | 0.343367 | 0.339477 | 0.661108 | 26 | 0 | 0.0866198 | 0.538057 | 0.167829 | 0.471393 | 1.22678e-13 | 0.0866198 | 0.0866198 | union_pair_exchange_skew / union_C_square_metric | False |

Decision:

```json
{
  "any_non_strict_primary_projective_spinor_pass": true,
  "any_non_strict_primary_strong_unscaled_spinor_pass": false
}
```
