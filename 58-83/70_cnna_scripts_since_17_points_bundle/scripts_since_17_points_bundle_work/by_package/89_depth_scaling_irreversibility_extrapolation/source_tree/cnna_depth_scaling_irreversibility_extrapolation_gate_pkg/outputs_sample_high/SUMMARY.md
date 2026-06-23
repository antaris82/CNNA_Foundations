# CNNA depth-scaling irreversibility extrapolation gate

## Purpose

Test whether Record/Live irreversibility in the established script-1/script-2 true Schur/DtN growth model scales with depth rather than being a constant single-birth feature.

Finite true-Schur/DtN computation is complete through L3.  For higher levels this package uses deterministic sampled frontier approximants with cap `8`; all such rows are marked.  This avoids pretending that the high-L data are full ternary trees.

No `J`, `i`, Hodge, `*`, C*-claim, positivity axiom, Q/P target or `delta_beta` decision is used.

## Deepest-run summary

| variant | events | relax rows | sampled? | mean record/live gap | mean birth ΔDtN | live/birth ratio | reverse nonreconstructability proxy |
|---|---:|---:|---:|---:|---:|---:|---:|
| real_growth_linear_depth_scaling | 111 | 111 | 1 | 0.03081 | 0.3033 | 0.1026 | 1.000 |
| log_growth_depth_scaling | 111 | 111 | 1 | 0.01654 | 0.1718 | 0.09991 | 1.000 |
| saturating_growth_depth_scaling | 111 | 111 | 1 | 0.2293 | 0.6592 | 0.3633 | 1.000 |
| strict_sym_depth_control | 111 | 111 | 1 | 0 | 0 | 0 | 0.000 |

## Slope audit

| scope | variant | depth groups | log gap slope vs parent level | deep/shallow gap ratio | log birth shock slope |
|---|---|---:|---:|---:|---:|
| deepest_run_only | log_growth_depth_scaling | 6 | 0.048 | 1.355 | -0.050 |
| all_runs_aggregate | log_growth_depth_scaling | 11 | 0.053 | 1.355 | -0.054 |
| deepest_run_only | real_growth_linear_depth_scaling | 6 | 0.149 | 2.264 | 0.062 |
| all_runs_aggregate | real_growth_linear_depth_scaling | 11 | 0.159 | 2.264 | 0.063 |
| deepest_run_only | saturating_growth_depth_scaling | 6 | 0.030 | 1.265 | -0.075 |
| all_runs_aggregate | saturating_growth_depth_scaling | 11 | 0.039 | 1.265 | -0.080 |
| deepest_run_only | strict_sym_depth_control | 6 | 0.000 | 0.000 | 0.000 |
| all_runs_aggregate | strict_sym_depth_control | 11 | 0.000 | 0.000 | 0.000 |

## Interpretation

The key diagnostic is not a local `J` or orientation lock.  It is the depth behavior of the true Schur/DtN Record→Live gap after birth-plus-relaxation.

A positive log-gap slope and a deep/shallow ratio above one indicate that irreversibility grows with birth depth in these finite/sampled approximants.  The strict-sym control must stay null.  High-L sampled rows are evidence for a trend, not a full infinite-tree theorem.

## Next test

`test_stable_small_real_star_algebra_gate.py` should use the depth trend to filter for sufficiently deep, bridge-/passivity-positive Record→Live operator rows and then test whether the generated small real operator family becomes stable under the Record-DtN metric adjoint `#`, products and commutators.  This should remain weaker than a physical `*`-algebra or C*-claim.
