# CNNA completed operator-family depth stability gate

## Purpose

This package tests whether the commutator-completed, coarse-grained Record→Live true-Schur/DtN operator family remains stable across finite approximant depth and parent-generation.

It does **not** introduce `J`, `i`, Hodge, physical Hilbert adjoint, C*-norm, Q/P targets, positivity axioms, or delta-beta decisions.

## Model and limitation

The test reuses the established script-1/script-2 true-Schur/DtN growth line.  Full true-Schur/DtN growth is run up to `full_until_level = 3`.  Above that level, frontier parents are selected by deterministic sampling with cap `8`.  Therefore L4 here is a sampled/frontier approximant if `full_until_level < 4`; this package does **not** claim a full infinite-tree theorem.

## Best deepest-level rows

| variant | coarse family | completion | level | strong pass | product resid | commutator resid | rank by level | stable? |
|---|---|---|---:|---:|---:|---:|---|---|
| real_growth_linear_star_candidate | depth_group_average | commutator_completed_2 | 4 | 0.567 | 0.001285 | 2.19e-08 | 2.5000 6.2500 9.0333 | True |
| saturating_growth_star_candidate | depth_group_average | commutator_completed_2 | 4 | 0.322 | 0.01599 | 2.304e-07 | 4.6667 12.2121 13.6102 | False |
| real_growth_linear_star_candidate | mean_plus_principal | commutator_completed_2 | 4 | 0.567 | 0.0006001 | 1.963e-05 | 2.5000 6.2500 10.1333 | True |
| real_growth_linear_star_candidate | mean_plus_principal | commutator_completed_1 | 4 | 0.567 | 0.002246 | 8.892e-05 | 2.5000 6.0625 8.4667 | True |
| saturating_growth_star_candidate | mean_plus_principal | commutator_completed_2 | 4 | 0.322 | 0.0082 | 0.0002221 | 4.6667 14.3030 16.3390 | False |
| saturating_growth_star_candidate | mean_plus_principal | commutator_completed_1 | 4 | 0.322 | 0.01991 | 0.000418 | 4.6667 10.7879 11.9322 | False |
| real_growth_linear_star_candidate | principal_stable_modes | commutator_completed_2 | 4 | 0.400 | 0.01757 | 0.006854 | 2.5000 6.0625 8.1667 | False |
| saturating_growth_star_candidate | principal_stable_modes | commutator_completed_2 | 4 | 0.237 | 0.05208 | 0.007055 | 4.6667 10.3030 11.3220 | False |
| saturating_growth_star_candidate | depth_group_average | commutator_completed_1 | 4 | 0.237 | 0.05867 | 0.0386 | 4.6667 8.1212 8.7458 | False |
| real_growth_linear_star_candidate | depth_group_average | commutator_completed_1 | 4 | 0.400 | 0.04299 | 0.07019 | 2.5000 5.3125 6.7000 | False |

## Level-stability passes

| variant | coarse family | completion | last strong | last product | last commutator |
|---|---|---|---:|---:|---:|
| real_growth_linear_star_candidate | depth_group_average | commutator_completed_2 | 0.567 | 0.001285 | 2.19e-08 |
| real_growth_linear_star_candidate | mean_plus_principal | commutator_completed_1 | 0.567 | 0.002246 | 8.892e-05 |
| real_growth_linear_star_candidate | mean_plus_principal | commutator_completed_2 | 0.567 | 0.0006001 | 1.963e-05 |

## Generation-stability passes at deepest level

| variant | coarse family | completion | median parent strong | max parent product | max parent comm | parent levels |
|---|---|---|---:|---:|---:|---|
| real_growth_linear_star_candidate | depth_group_average | commutator_completed_2 | 0.500 | 0.002142 | 3.642e-08 | 1 2 3 |
| real_growth_linear_star_candidate | mean_plus_principal | commutator_completed_1 | 0.500 | 0.003743 | 0.0001482 | 1 2 3 |
| real_growth_linear_star_candidate | mean_plus_principal | commutator_completed_2 | 0.500 | 0.001 | 3.271e-05 | 1 2 3 |
| real_growth_linear_star_candidate | principal_stable_modes | commutator_completed_2 | 0.500 | 0.02797 | 0.01047 | 1 2 3 |

## Strict-sym audit

```json
{
  "strict_sym_rows": 0,
  "strict_sym_any_strong_pass": false,
  "strict_sym_any_weak_pass": false,
  "strict_sym_max_product_resid": 0.0,
  "strict_sym_max_comm_resid": 0.0,
  "strict_sym_null_gate": true,
  "used_delta_beta_any": false
}
```

## Interpretation

A positive row means only this:

```text
The completed coarse-grained real operator family is stable under #_G, products,
and commutators across the finite/sampled depth audit.
```

It does **not** mean:

```text
C*-algebra, Hilbert positivity, J, i, or a physical quantum algebra.
```

If the best pass rows are present but rely on sampled L4, the correct next step is to harden the depth evidence with either a more efficient Schur/DtN implementation or a narrower full-L4 subset chosen by a predeclared rule.
