# Results: double-history parent-fan plaquette J-holonomy

This test combines two previously separate structures:

1. parent-fan plaquette J-candidates from full DtN face commutators;
2. double-history suffix identification across root sectors.

It deliberately avoids minimal SO(3) axis transport.  The primary gluing law is
provenance-directed: root sector 1 -> 2 -> 3 -> 1, with face labels shifted by
the directed root-sector difference.

## Current run

```text
CNNA DOUBLE-HISTORY PARENT-FAN PLAQUETTE J-HOLONOMY DIAGNOSTIC
max_level=5, mode=linear, min_omega=1e-10

MODEL SUMMARIES
  real_growth: nodes=364, cells=121, jnets=39, double_classes=13, three_sector_cycles=13
  symmetrized_birth: nodes=364, cells=121, jnets=39, double_classes=13, three_sector_cycles=13
  no_backreaction: nodes=364, cells=121, jnets=39, double_classes=13, three_sector_cycles=13

SELECTED SUMMARIES
  double_history_parent_fan_cycle | real_growth | live | full | directed_cyclic_root_shift: count=13, valid=1.000, K=0.00156178, Jfixed=0.6684, J±=0.005235, axis_abs=0.9998, z2obs=0, z2prod=1, track±=0.005235
  double_history_parent_fan_cycle | real_growth | live | full | same_label: count=13, valid=1.000, K=0.00156178, Jfixed=0.004999, J±=0.004999, axis_abs=0.9999, z2obs=0, z2prod=1, track±=0.004999
  double_history_parent_fan_cycle | real_growth | live | full | reverse_cyclic_root_shift: count=13, valid=1.000, K=0.00156178, Jfixed=0.6684, J±=0.005104, axis_abs=0.9998, z2obs=0, z2prod=1, track±=0.005104
  double_history_parent_fan_cycle | symmetrized_birth | live | full | directed_cyclic_root_shift: count=13, valid=1.000, K=6.4518e-06, Jfixed=0.6667, J±=5.595e-06, axis_abs=1.0000, z2obs=0, z2prod=1, track±=5.595e-06
  double_history_parent_fan_cycle | real_growth | handoff | full | directed_cyclic_root_shift: count=13, valid=1.000, K=0.0001036, Jfixed=0.6676, J±=0.002076, axis_abs=1.0000, z2obs=0, z2prod=1, track±=0.002076
  identical_history_control | real_growth | live | full | directed_cyclic_root_shift: count=39, valid=1.000, K=0.00156178, Jfixed=0, J±=0, axis_abs=1.0000, z2obs=0, z2prod=1, track±=0
  random_same_level_cycle_baseline | real_growth | live | full | directed_cyclic_root_shift: count=40, valid=1.000, K=0.000990745, Jfixed=0.4313, J±=0.007472, axis_abs=0.9998, z2obs=0, z2prod=1, track±=0.007472
  double_history_parent_fan_cycle | real_growth | live | diagonal | directed_cyclic_root_shift: count=13, valid=0.000, K=0, Jfixed=nan, J±=nan, axis_abs=nan, z2obs=nan, z2prod=nan, track±=nan
  double_history_parent_fan_cycle | real_growth | live | trace_scalar | directed_cyclic_root_shift: count=13, valid=0.000, K=0, Jfixed=nan, J±=nan, axis_abs=nan, z2obs=nan, z2prod=nan, track±=nan

READING RULE
  directed_cyclic_root_shift uses root-sector birth order 1->2->3->1 and shifts
  parent-fan face labels by the directed sector difference. This is the primary
  provenance-identification test; it is not a minimal SO(3) best-fit transport.
  identical_history_control should have zero/near-zero residuals and no obstruction.
  diagonal/trace/common-commuting reductions should invalidate/kill K.
  symmetrized_birth should suppress the K amplitude and any robust holonomy signal.
  A nonzero z2obs would indicate a sign-lift obstruction around the root-sector
  double-history cycle; a small J± with high axis_abs but large Jfixed indicates
  signless plane gluing without fixed-sign locking.
```

## Output files

- `double_history_J_cycle_rows.csv`
- `double_history_J_edge_rows.csv`
- `double_history_J_summary_main.csv`
- `double_history_J_summary_by_level.csv`
- `model_summaries.csv`
- `SUMMARY.txt`

## Interpretation status

This is still not a derivation of a global complex structure.  It tests whether
the local parent-fan Z3/J sectors remain nontrivial when glued by the same
suffix-forgetting double-history quotient that previously exposed DtN handoff
defects.

A positive Stufe-4 candidate would require more than local J-validity: it would
need a double-history-specific sign/holonomy obstruction that vanishes for
identical history, collapses under symmetrized birth and commuting reductions,
and is not reproduced by random same-level cycles.
