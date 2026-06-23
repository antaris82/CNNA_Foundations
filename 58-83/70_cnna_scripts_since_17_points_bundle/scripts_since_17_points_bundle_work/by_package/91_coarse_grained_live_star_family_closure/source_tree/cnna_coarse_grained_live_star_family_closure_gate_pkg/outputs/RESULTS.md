# CNNA coarse-grained live #/*-family closure gate

## Purpose

This package tests whether the missing closure of the live semigroup operator family is a coarse-graining issue.

It compares:

```text
raw_per_cut
depth_group_average
principal_stable_modes
mean_plus_principal
```

on actual Record→Live true Schur/DtN operators, not just stored scalar residuals.

No `J`, `i`, Hodge, physical Hilbert adjoint, C*-norm, Q/P target, or delta-beta decision is used.

## Deepest-level summary

| variant | coarse level | rows | weak # family pass | strong stable pass | product resid | commutator resid | commutator norm | rank |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| real_growth_linear_star_candidate | depth_group_average | 30 | 0.233 | 0.233 | 0.014 | 0.287 | 0.002 | 3.83 |
| real_growth_linear_star_candidate | mean_plus_principal | 30 | 0.233 | 0.233 | 0.009 | 0.162 | 0.002 | 5.43 |
| real_growth_linear_star_candidate | principal_stable_modes | 30 | 0.233 | 0.233 | 0.180 | 0.244 | 0.250 | 3.30 |
| real_growth_linear_star_candidate | raw_per_cut | 198 | 0.000 | 0.000 | 0.049 | 0.546 | 0.003 | 4.36 |
| saturating_growth_star_candidate | depth_group_average | 59 | 0.119 | 0.119 | 0.060 | 0.427 | 0.003 | 4.41 |
| saturating_growth_star_candidate | mean_plus_principal | 59 | 0.119 | 0.119 | 0.034 | 0.242 | 0.003 | 6.69 |
| saturating_growth_star_candidate | principal_stable_modes | 59 | 0.119 | 0.119 | 0.353 | 0.373 | 0.344 | 3.64 |
| saturating_growth_star_candidate | raw_per_cut | 198 | 0.015 | 0.000 | 0.053 | 0.523 | 0.003 | 4.36 |
| strict_symmetrized_response_star_control | raw_per_cut | 198 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 1.00 |

## Interpretation

A successful coarse-graining result would require lower product and commutator residuals than the raw per-cut family, while preserving `#` closure and keeping strict-sym null.

`principal_stable_modes` and `mean_plus_principal` are the important tests: if commutator closure improves there, the raw failure was likely a non-coarse-grained noise/mode-mixing problem. If not, the missing algebraic closure is structural or needs a different operator family.

## Next step

If coarse-graining improves only the product residual but not commutator residual, the next test should isolate the missing commutator generator rather than claim a real `*`-algebra.  A natural next package would be:

```text
test_commutator_generator_completion_gate.py
```

where one adds only derived commutator modes and checks whether closure stabilizes without importing `J` or complex structure.
