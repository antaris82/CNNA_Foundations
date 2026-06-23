# Results: nonflat connection from same-label response mismatch

This test accepts identity/same-label face gluing and asks whether the remaining
DtN/response mismatch defines a nonflat connection around double-history cycles.

```text
CNNA NONFLAT CONNECTION FROM SAME-LABEL RESPONSE MISMATCH
max_level=6, mode=linear, min_omega=1e-10, ridge=1e-09

MODEL SUMMARIES
  real_growth: nodes=1093, cells=364, jnets=120, double_classes=40, three_sector_cycles=40
  symmetrized_birth: nodes=1093, cells=364, jnets=120, double_classes=40, three_sector_cycles=40
  no_backreaction: nodes=1093, cells=364, jnets=120, double_classes=40, three_sector_cycles=40

SELECTED SUMMARIES
  double_history_suffix_cycle | real_growth | live | full | mean_face: count=120, valid=1.000, K=0.00186993, rawJ±=0.003562, plane=0.003561, amp=0.03211, metric=0.01092, Gcomm=0.0001223, T=0.01205, metricWilson=3.501e-11, KWilson=1.003e-15, GcommWilson=1.75e-09, mixedWilson=0, hol=1.558e-15, hol_skew=2.767e-16, hol_complex=0.000, loopJ±=9.269e-16, transJ±=0.00457, cond=5.833
  missing
  double_history_suffix_cycle | symmetrized_birth | live | full | mean_face: count=120, valid=1.000, K=7.08903e-06, rawJ±=2.65e-05, plane=2.65e-05, amp=0.003168, metric=0.0005777, Gcomm=6.549e-07, T=0.0006899, metricWilson=9.187e-16, KWilson=2.298e-25, GcommWilson=4.511e-14, mixedWilson=0, hol=1.494e-15, hol_skew=2.281e-16, hol_complex=0.000, loopJ±=8.476e-16, transJ±=0.0001639, cond=5.511
  double_history_suffix_cycle | no_backreaction | live | full | mean_face: count=120, valid=1.000, K=0.00363105, rawJ±=0.00128, plane=0.001279, amp=0.03158, metric=0.01179, Gcomm=0.000204, T=0.01328, metricWilson=1.337e-10, KWilson=2.209e-15, GcommWilson=6.684e-09, mixedWilson=0, hol=1.554e-15, hol_skew=2.863e-16, hol_complex=0.000, loopJ±=9.21e-16, transJ±=0.003127, cond=5.79
  identical_history_clone_control | real_growth | live | full | mean_face: count=120, valid=1.000, K=0.000249136, rawJ±=0, plane=0, amp=0, metric=0, Gcomm=0, T=8.819e-16, metricWilson=0, KWilson=0, GcommWilson=0, mixedWilson=0, hol=2.646e-15, hol_skew=4.022e-16, hol_complex=0.000, loopJ±=1.578e-15, transJ±=5.3e-16, cond=4.156
  random_same_level_cycle_baseline | real_growth | live | full | mean_face: count=120, valid=1.000, K=0.000569692, rawJ±=0.01726, plane=0.01725, amp=0.07379, metric=0.0226, Gcomm=0.0002272, T=0.02392, metricWilson=3.75e-09, KWilson=2.939e-14, GcommWilson=1.875e-07, mixedWilson=0, hol=1.508e-15, hol_skew=2.273e-16, hol_complex=0.000, loopJ±=8.904e-16, transJ±=0.018, cond=4.476
  double_history_suffix_cycle | real_growth | live | diagonal | mean_face: count=120, valid=0.000, K=0, rawJ±=nan, plane=nan, amp=0, metric=0.01112, Gcomm=0, T=0.01114, metricWilson=1.047e-16, KWilson=0, GcommWilson=0, mixedWilson=0, hol=1.468e-16, hol_skew=0, hol_complex=0.000, loopJ±=nan, transJ±=nan, cond=1.173
  double_history_suffix_cycle | real_growth | live | trace_scalar | mean_face: count=120, valid=0.000, K=0, rawJ±=nan, plane=nan, amp=0, metric=0.01112, Gcomm=0, T=0.01114, metricWilson=8.234e-17, KWilson=0, GcommWilson=0, mixedWilson=0, hol=9.807e-17, hol_skew=0, hol_complex=0.000, loopJ±=nan, transJ±=nan, cond=1
  double_history_suffix_cycle | real_growth | handoff | full | mean_face: count=120, valid=1.000, K=0.000115026, rawJ±=0.0003323, plane=0.0003323, amp=0.01832, metric=0.009805, Gcomm=1.523e-05, T=0.009832, metricWilson=8.87e-12, KWilson=5.814e-17, GcommWilson=4.435e-10, mixedWilson=0, hol=1.483e-15, hol_skew=6.79e-17, hol_complex=0.000, loopJ±=8.565e-16, transJ±=0.0003452, cond=1.617
  double_history_suffix_cycle | real_growth | aging | full | mean_face: count=120, valid=1.000, K=0.000111542, rawJ±=0.0003678, plane=0.0003678, amp=0.01832, metric=0.009788, Gcomm=1.575e-05, T=0.009818, metricWilson=8.577e-12, KWilson=5.74e-17, GcommWilson=4.289e-10, mixedWilson=0, hol=1.517e-15, hol_skew=6.587e-17, hol_complex=0.000, loopJ±=8.383e-16, transJ±=0.0003761, cond=1.613

READING RULE
  Identity/same-label face gluing is taken as the boundary-resolved result
  from the previous test.  This script does not search over face labels and
  does not impose a cyclic Z3 shift.
  The response connection uses only symmetric DtN face metrics G and response
  mismatch generators on same-label faces.  The exact metric transport
  T_ij = G_j^{-1/2} G_i^{1/2} is included only as a gauge audit, because it
  can telescope by construction.  The Wilson columns metricWilson, KWilson,
  GcommWilson and mixedWilson test the response-mismatch generators themselves.
  A positive nonflat-connection gate would need real_growth/full/live to have
  Wilson curvature, transported-J loop residual or holonomy clearly above identical-history,
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
