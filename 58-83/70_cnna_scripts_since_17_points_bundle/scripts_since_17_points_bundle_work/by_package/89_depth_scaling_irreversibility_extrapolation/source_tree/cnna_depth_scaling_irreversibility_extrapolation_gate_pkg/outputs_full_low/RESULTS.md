# CNNA depth-scaling irreversibility extrapolation gate

## Purpose

Test whether Record/Live irreversibility in the established script-1/script-2 true Schur/DtN growth model scales with depth rather than being a constant single-birth feature.

Finite true-Schur/DtN computation is complete through L4.  For higher levels this package uses deterministic sampled frontier approximants with cap `18`; all such rows are marked.  This avoids pretending that the high-L data are full ternary trees.

No `J`, `i`, Hodge, `*`, C*-claim, positivity axiom, Q/P target or `delta_beta` decision is used.

## Deepest-run summary

| variant | events | relax rows | sampled? | mean record/live gap | mean birth ΔDtN | live/birth ratio | reverse nonreconstructability proxy |
|---|---:|---:|---:|---:|---:|---:|---:|
| real_growth_linear_depth_scaling | 120 | 120 | 0 | 0.03047 | 0.2929 | 0.1051 | 1.000 |
| log_growth_depth_scaling | 120 | 120 | 0 | 0.01706 | 0.1728 | 0.102 | 1.000 |
| saturating_growth_depth_scaling | 120 | 120 | 0 | 0.2619 | 0.6544 | 0.4152 | 1.000 |
| strict_sym_depth_control | 120 | 120 | 0 | 0 | 0 | 0 | 0.000 |

## Slope audit

| scope | variant | depth groups | log gap slope vs parent level | deep/shallow gap ratio | log birth shock slope |
|---|---|---:|---:|---:|---:|
| deepest_run_only | log_growth_depth_scaling | 4 | 0.109 | 1.414 | -0.078 |
| all_runs_aggregate | log_growth_depth_scaling | 9 | 0.125 | 1.414 | -0.083 |
| deepest_run_only | real_growth_linear_depth_scaling | 4 | 0.237 | 2.077 | 0.068 |
| all_runs_aggregate | real_growth_linear_depth_scaling | 9 | 0.258 | 2.077 | 0.070 |
| deepest_run_only | saturating_growth_depth_scaling | 4 | 0.149 | 1.609 | -0.111 |
| all_runs_aggregate | saturating_growth_depth_scaling | 9 | 0.162 | 1.609 | -0.105 |
| deepest_run_only | strict_sym_depth_control | 4 | 0.000 | 0.000 | 0.000 |
| all_runs_aggregate | strict_sym_depth_control | 9 | 0.000 | 0.000 | 0.000 |

## Interpretation

The key diagnostic is not a local `J` or orientation lock.  It is the depth behavior of the true Schur/DtN Record→Live gap after birth-plus-relaxation.

A positive log-gap slope and a deep/shallow ratio above one indicate that irreversibility grows with birth depth in these finite/sampled approximants.  The strict-sym control must stay null.  High-L sampled rows are evidence for a trend, not a full infinite-tree theorem.

## Next test

`test_stable_small_real_star_algebra_gate.py` should use the depth trend to filter for sufficiently deep, bridge-/passivity-positive Record→Live operator rows and then test whether the generated small real operator family becomes stable under the Record-DtN metric adjoint `#`, products and commutators.  This should remain weaker than a physical `*`-algebra or C*-claim.
