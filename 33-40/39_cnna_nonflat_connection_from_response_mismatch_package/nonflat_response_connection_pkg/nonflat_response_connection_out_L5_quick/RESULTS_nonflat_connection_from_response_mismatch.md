# Results: nonflat connection from same-label response mismatch

This test accepts identity/same-label face gluing and asks whether the remaining
DtN/response mismatch defines a nonflat connection around double-history cycles.

```text
CNNA NONFLAT CONNECTION FROM SAME-LABEL RESPONSE MISMATCH
max_level=5, mode=linear, min_omega=1e-10, ridge=1e-09

MODEL SUMMARIES
  real_growth: nodes=364, cells=121, jnets=39, double_classes=13, three_sector_cycles=13
  symmetrized_birth: nodes=364, cells=121, jnets=39, double_classes=13, three_sector_cycles=13
  no_backreaction: nodes=364, cells=121, jnets=39, double_classes=13, three_sector_cycles=13

SELECTED SUMMARIES
  double_history_suffix_cycle | real_growth | live | full | mean_face: count=39, valid=1.000, K=0.00156178, rawJ±=0.004999, plane=0.004998, amp=0.0426, metric=0.01353, Gcomm=0.0001657, T=0.01517, metricWilson=5.925e-11, KWilson=1.968e-15, GcommWilson=2.962e-09, mixedWilson=0, hol=1.723e-15, hol_skew=3.104e-16, hol_complex=0.000, loopJ±=9.098e-16, transJ±=0.00641, cond=6.048
  double_history_suffix_cycle | real_growth | live | full | edge_energy: count=39, valid=1.000, K=0.00156178, rawJ±=0.004999, plane=0.004998, amp=0.0426, metric=0.04197, Gcomm=0.0003453, T=0.04736, metricWilson=5.655e-11, KWilson=1.968e-15, GcommWilson=2.827e-09, mixedWilson=0, hol=1.561e-15, hol_skew=2.402e-16, hol_complex=0.000, loopJ±=9.444e-16, transJ±=0.02217, cond=6.256
  double_history_suffix_cycle | symmetrized_birth | live | full | mean_face: count=39, valid=1.000, K=6.4518e-06, rawJ±=4.7e-06, plane=4.7e-06, amp=0.002368, metric=0.0007551, Gcomm=8.692e-07, T=0.0009119, metricWilson=1.581e-15, KWilson=4.862e-25, GcommWilson=7.786e-14, mixedWilson=0, hol=1.353e-15, hol_skew=2.241e-16, hol_complex=0.000, loopJ±=8.039e-16, transJ±=0.0002107, cond=5.727
  double_history_suffix_cycle | no_backreaction | live | full | mean_face: count=39, valid=1.000, K=0.00321731, rawJ±=0.001406, plane=0.001406, amp=0.04037, metric=0.01496, Gcomm=0.0002811, T=0.01717, metricWilson=2.712e-10, KWilson=6.594e-15, GcommWilson=1.356e-08, mixedWilson=0, hol=1.684e-15, hol_skew=3.36e-16, hol_complex=0.000, loopJ±=9.411e-16, transJ±=0.004014, cond=6.004
  identical_history_clone_control | real_growth | live | full | mean_face: count=117, valid=1.000, K=0.00156178, rawJ±=0, plane=0, amp=0, metric=0, Gcomm=0, T=8.872e-16, metricWilson=0, KWilson=0, GcommWilson=0, mixedWilson=0, hol=2.661e-15, hol_skew=4.694e-16, hol_complex=0.000, loopJ±=1.431e-15, transJ±=4.772e-16, cond=6.048
  random_same_level_cycle_baseline | real_growth | live | full | mean_face: count=120, valid=1.000, K=0.00121732, rawJ±=0.006644, plane=0.006643, amp=0.06046, metric=0.02108, Gcomm=0.0002454, T=0.02328, metricWilson=2.876e-09, KWilson=3.619e-14, GcommWilson=1.438e-07, mixedWilson=0, hol=1.535e-15, hol_skew=2.641e-16, hol_complex=0.000, loopJ±=8.184e-16, transJ±=0.008484, cond=5.632
  double_history_suffix_cycle | real_growth | live | diagonal | mean_face: count=39, valid=0.000, K=0, rawJ±=nan, plane=nan, amp=0, metric=0.01381, Gcomm=0, T=0.01383, metricWilson=1.034e-16, KWilson=0, GcommWilson=0, mixedWilson=0, hol=1.425e-16, hol_skew=0, hol_complex=0.000, loopJ±=nan, transJ±=nan, cond=1.171
  double_history_suffix_cycle | real_growth | live | trace_scalar | mean_face: count=39, valid=0.000, K=0, rawJ±=nan, plane=nan, amp=0, metric=0.0138, Gcomm=0, T=0.01384, metricWilson=8.825e-17, KWilson=0, GcommWilson=0, mixedWilson=0, hol=9.964e-17, hol_skew=0, hol_complex=0.000, loopJ±=nan, transJ±=nan, cond=1
  missing
  missing

READING RULE
  Identity/same-label face gluing is taken as the boundary-resolved result
  from the previous test.  This script does not search over face labels and
  does not impose a cyclic Z3 shift.
  The response connection uses only symmetric DtN face metrics G and the
  canonical metric transport T_ij = G_j^{-1/2} G_i^{1/2}.  Its loop holonomy
  H = T_31 T_23 T_12 is nontrivial only if the response metrics fail to
  telescope/commute around the double-history cycle.
  A positive nonflat-connection gate would need real_growth/full/live to have
  holonomy or transported-J loop residual clearly above identical-history,
  diagonal/trace and random controls, and preferably stronger than
  symmetrized_birth.  Otherwise the local parent-fan J sector remains a local
  Stufe-2/3 phenomenon without Stufe-4 locking.
```

## Output files

- `response_connection_loop_rows.csv`
- `response_connection_edge_rows.csv`
- `response_connection_summary_main.csv`
- `response_connection_summary_by_label.csv`
- `response_connection_summary_by_level.csv`
- `model_summaries.csv`
- `SUMMARY.txt`

## Status

The primary row is `double_history_suffix_cycle | real_growth | live | full |
mean_face`.  The strongest anti-smuggling controls are identical-history,
diagonal/trace reduction and random same-level cycles.
