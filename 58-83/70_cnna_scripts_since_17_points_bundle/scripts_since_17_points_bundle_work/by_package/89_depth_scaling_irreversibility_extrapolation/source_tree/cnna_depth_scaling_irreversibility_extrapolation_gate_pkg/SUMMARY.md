# CNNA depth-scaling irreversibility extrapolation gate

## Purpose

This package tests the depth-scaling claim for the established script-1/script-2 true Schur/DtN growth + live-relaxation model:

```text
Irreversibility is not merely a property of one birth event.
Birth effects distribute into old conductances and then mix with live relaxation;
therefore the Record/Live gap should grow with depth/generation.
```

No `J`, `i`, Hodge, `*`, C*-claim, positivity axiom, Q/P target or `delta_beta` decision is used.

## Compute regime

- `full_L2_L4`: complete ternary finite approximants through L4.
- `sampled_L5_L6`: deterministic sampled frontier approximants for L5/L6, because full true-Schur/DtN L5/L6 is too expensive here. Sampled rows are marked and must not be read as full-tree averages.

## Summary at L4 and L6

| run | level | variant | events | sampled? | mean Record/Live gap | mean birth ΔDtN | live/birth ratio | reverse nonreconstructability proxy |
|---|---:|---|---:|---:|---:|---:|---:|---:|
| full_L2_L4 | L4 | real_growth_linear_depth_scaling | 120 | 0 | 0.03047 | 0.2929 | 0.1051 | 1.000 |
| full_L2_L4 | L4 | log_growth_depth_scaling | 120 | 0 | 0.01706 | 0.1728 | 0.102 | 1.000 |
| full_L2_L4 | L4 | saturating_growth_depth_scaling | 120 | 0 | 0.2619 | 0.6544 | 0.4152 | 1.000 |
| full_L2_L4 | L4 | strict_sym_depth_control | 120 | 0 | 0 | 0 | 0 | 0.000 |
| sampled_L5_L6 | L6 | real_growth_linear_depth_scaling | 111 | 1 | 0.03081 | 0.3033 | 0.1026 | 1.000 |
| sampled_L5_L6 | L6 | log_growth_depth_scaling | 111 | 1 | 0.01654 | 0.1718 | 0.09991 | 1.000 |
| sampled_L5_L6 | L6 | saturating_growth_depth_scaling | 111 | 1 | 0.2293 | 0.6592 | 0.3633 | 1.000 |
| sampled_L5_L6 | L6 | strict_sym_depth_control | 111 | 1 | 0 | 0 | 0 | 0.000 |

## Depth-slope audit

| run | scope | variant | depth groups | log gap slope | deep/shallow gap ratio | log birth-shock slope |
|---|---|---|---:|---:|---:|---:|
| full_L2_L4 | deepest_run_only | log_growth_depth_scaling | 4 | 0.109 | 1.414 | -0.078 |
| full_L2_L4 | all_runs_aggregate | log_growth_depth_scaling | 9 | 0.125 | 1.414 | -0.083 |
| full_L2_L4 | deepest_run_only | real_growth_linear_depth_scaling | 4 | 0.237 | 2.077 | 0.068 |
| full_L2_L4 | all_runs_aggregate | real_growth_linear_depth_scaling | 9 | 0.258 | 2.077 | 0.070 |
| full_L2_L4 | deepest_run_only | saturating_growth_depth_scaling | 4 | 0.149 | 1.609 | -0.111 |
| full_L2_L4 | all_runs_aggregate | saturating_growth_depth_scaling | 9 | 0.162 | 1.609 | -0.105 |
| full_L2_L4 | deepest_run_only | strict_sym_depth_control | 4 | 0.000 | 0.000 | 0.000 |
| full_L2_L4 | all_runs_aggregate | strict_sym_depth_control | 9 | 0.000 | 0.000 | 0.000 |
| sampled_L5_L6 | deepest_run_only | log_growth_depth_scaling | 6 | 0.048 | 1.355 | -0.050 |
| sampled_L5_L6 | all_runs_aggregate | log_growth_depth_scaling | 11 | 0.053 | 1.355 | -0.054 |
| sampled_L5_L6 | deepest_run_only | real_growth_linear_depth_scaling | 6 | 0.149 | 2.264 | 0.062 |
| sampled_L5_L6 | all_runs_aggregate | real_growth_linear_depth_scaling | 11 | 0.159 | 2.264 | 0.063 |
| sampled_L5_L6 | deepest_run_only | saturating_growth_depth_scaling | 6 | 0.030 | 1.265 | -0.075 |
| sampled_L5_L6 | all_runs_aggregate | saturating_growth_depth_scaling | 11 | 0.039 | 1.265 | -0.080 |
| sampled_L5_L6 | deepest_run_only | strict_sym_depth_control | 6 | 0.000 | 0.000 | 0.000 |
| sampled_L5_L6 | all_runs_aggregate | strict_sym_depth_control | 11 | 0.000 | 0.000 | 0.000 |

## Interpretation

The full L2-L4 runs already show positive depth slopes in all nontrivial growth variants and a null strict-sym control. The sampled L5/L6 runs preserve the same qualitative trend but are weaker evidence because the frontier is capped.

The key result is not an operator or `J` claim. It is the finite/sampled scaling statement:

```text
Record/Live Schur-DtN irreversibility grows with depth in the response growth variants,
while strict_sym stays null.
```

This supports the user's point that the earlier L2-local orientation searches may have been too close to the nearly reversible root regime.

## Next package

`cnna_stable_small_real_star_algebra_gate_pkg`: after securing the depth trend, test whether deep/passivity-positive Record→Live rows form a stable small real #/*-algebra-like operator family.
