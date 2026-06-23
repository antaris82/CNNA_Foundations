# Results: parent-fan plaquette skew-sector J-candidate diagnostic

This is the first J-near test after the face-level plaquette-frustration and
orientation/coboundary audits.

It still does **not** introduce a complex Hilbert space, a C*-norm, GNS, AQFT,
or a global complex algebra.  The only candidate is extracted from the measured
real skew commutator

```text
K_abc = [A_ab, A_bc]
```

on oriented provenance-generated 2-simplices.

For each nonzero 3x3 skew `K`, the test extracts its rank-2 image plane and
normalizes `J = K / omega` on that plane.  The nontrivial issue is not the
linear algebra fact itself, but whether the robust parent-fan/live/full-DtN
carrier survives controls and has the right orientation behavior.

## Current run

```text
CNNA PARENT-FAN PLAQUETTE SKEW-SECTOR J-CANDIDATE DIAGNOSTIC
max_level=6, mode=linear, min_omega=1e-10

MODEL SUMMARIES
  real_growth: nodes=1093, completed=364, sibling=121, parent_fan=363, random=80
  symmetrized_birth: nodes=1093, completed=364, sibling=121, parent_fan=363, random=80
  no_backreaction: nodes=1093, completed=364, sibling=121, parent_fan=363, random=80

SELECTED FULL-DTN J-SECTOR SUMMARIES
  no_backreaction | parent_fan_triangle | aging | full: count=363, K=0, omega=0, J_valid=0.000, J2=nan, leak=nan, rank2=0, kgap=0, Jsign=nan, Ksign=0, axis_cos=0.000, plane_rev=nan, first=0.00e+00
  no_backreaction | parent_fan_triangle | handoff | full: count=363, K=9.17998e-20, omega=6.49123e-20, J_valid=0.000, J2=nan, leak=nan, rank2=6.49119e-08, kgap=0, Jsign=nan, Ksign=0, axis_cos=0.000, plane_rev=nan, first=0.00e+00
  no_backreaction | parent_fan_triangle | live | full: count=363, K=0.0036013, omega=0.0025465, J_valid=1.000, J2=4.05e-09, leak=1.8e-16, rank2=1, kgap=4.45e-09, Jsign=2.02e-15, Ksign=2.48e-15, axis_cos=1.000, plane_rev=2.03e-15, first=0.00e+00
  no_backreaction | random_same_level_triangle | aging | full: count=80, K=0, omega=0, J_valid=0.000, J2=nan, leak=nan, rank2=0, kgap=0, Jsign=nan, Ksign=0, axis_cos=0.000, plane_rev=nan, first=0.00e+00
  no_backreaction | random_same_level_triangle | handoff | full: count=80, K=0, omega=0, J_valid=0.000, J2=nan, leak=nan, rank2=0, kgap=0, Jsign=nan, Ksign=0, axis_cos=0.000, plane_rev=nan, first=0.00e+00
  no_backreaction | random_same_level_triangle | live | full: count=80, K=6.36467e-06, omega=4.5005e-06, J_valid=1.000, J2=6e-06, leak=2.44e-16, rank2=0.999997, kgap=5.54e-09, Jsign=5.63e-13, Ksign=7.21e-13, axis_cos=0.708, plane_rev=5.63e-13, first=0.00e+00
  no_backreaction | sibling_triangle | aging | full: count=121, K=0, omega=0, J_valid=0.000, J2=nan, leak=nan, rank2=0, kgap=0, Jsign=nan, Ksign=0, axis_cos=0.000, plane_rev=nan, first=0.00e+00
  no_backreaction | sibling_triangle | handoff | full: count=121, K=3.39936e-21, omega=2.40371e-21, J_valid=0.000, J2=nan, leak=nan, rank2=2.40371e-09, kgap=0, Jsign=nan, Ksign=0, axis_cos=0.000, plane_rev=nan, first=0.00e+00
  no_backreaction | sibling_triangle | live | full: count=121, K=4.16893e-06, omega=2.94788e-06, J_valid=1.000, J2=6.8e-07, leak=2.14e-16, rank2=1, kgap=3.93e-09, Jsign=2.56e-13, Ksign=3.21e-13, axis_cos=0.896, plane_rev=2.56e-13, first=0.00e+00
  real_growth | parent_fan_triangle | aging | full: count=363, K=0.000110915, omega=7.84285e-05, J_valid=1.000, J2=4.14e-08, leak=1.85e-16, rank2=1, kgap=4.12e-09, Jsign=1.71e-15, Ksign=2.16e-15, axis_cos=0.999, plane_rev=1.73e-15, first=0.00e+00
  real_growth | parent_fan_triangle | handoff | full: count=363, K=0.000114368, omega=8.08704e-05, J_valid=1.000, J2=4.14e-08, leak=1.9e-16, rank2=1, kgap=3.81e-09, Jsign=1.67e-15, Ksign=2.05e-15, axis_cos=0.999, plane_rev=1.69e-15, first=0.00e+00
  real_growth | parent_fan_triangle | live | full: count=363, K=0.00185628, omega=0.00131259, J_valid=1.000, J2=8.03e-09, leak=1.84e-16, rank2=1, kgap=4.4e-09, Jsign=4.64e-15, Ksign=6.05e-15, axis_cos=1.000, plane_rev=4.65e-15, first=0.00e+00
  real_growth | random_same_level_triangle | aging | full: count=80, K=2.13374e-06, omega=1.50878e-06, J_valid=0.975, J2=0.000582, leak=1.72e-16, rank2=0.999152, kgap=4.58e-09, Jsign=3.89e-14, Ksign=5.18e-14, axis_cos=0.427, plane_rev=3.89e-14, first=0.00e+00
  real_growth | random_same_level_triangle | handoff | full: count=80, K=2.15391e-06, omega=1.52304e-06, J_valid=1.000, J2=5.36e-05, leak=1.7e-16, rank2=0.999973, kgap=2.88e-09, Jsign=3.81e-14, Ksign=4.75e-14, axis_cos=0.430, plane_rev=3.81e-14, first=0.00e+00
  real_growth | random_same_level_triangle | live | full: count=80, K=1.71109e-05, omega=1.20992e-05, J_valid=1.000, J2=1.25e-06, leak=1.71e-16, rank2=0.999999, kgap=5.18e-09, Jsign=3.03e-13, Ksign=3.56e-13, axis_cos=0.805, plane_rev=3.03e-13, first=0.00e+00
  real_growth | sibling_triangle | aging | full: count=121, K=2.61871e-07, omega=1.85171e-07, J_valid=1.000, J2=0.00118, leak=2.24e-16, rank2=0.999411, kgap=4.35e-09, Jsign=4.18e-14, Ksign=5.36e-14, axis_cos=0.078, plane_rev=4.18e-14, first=0.00e+00
  real_growth | sibling_triangle | handoff | full: count=121, K=2.66732e-07, omega=1.88608e-07, J_valid=1.000, J2=0.00055, leak=1.98e-16, rank2=0.999725, kgap=3.74e-09, Jsign=3.24e-14, Ksign=4.08e-14, axis_cos=0.080, plane_rev=3.24e-14, first=0.00e+00
  real_growth | sibling_triangle | live | full: count=121, K=2.95477e-06, omega=2.08934e-06, J_valid=1.000, J2=9.87e-07, leak=1.68e-16, rank2=1, kgap=4.77e-09, Jsign=4.37e-13, Ksign=5.53e-13, axis_cos=0.804, plane_rev=4.37e-13, first=0.00e+00
  symmetrized_birth | parent_fan_triangle | aging | full: count=363, K=4.00871e-06, omega=2.83458e-06, J_valid=1.000, J2=4.94e-06, leak=1.85e-16, rank2=0.999998, kgap=3.85e-09, Jsign=1.66e-14, Ksign=1.99e-14, axis_cos=0.660, plane_rev=1.66e-14, first=0.00e+00
  symmetrized_birth | parent_fan_triangle | handoff | full: count=363, K=4.03911e-06, omega=2.85608e-06, J_valid=1.000, J2=4.85e-06, leak=1.94e-16, rank2=0.999998, kgap=3.76e-09, Jsign=1.65e-14, Ksign=2.06e-14, axis_cos=0.661, plane_rev=1.65e-14, first=0.00e+00
  symmetrized_birth | parent_fan_triangle | live | full: count=363, K=7.06818e-06, omega=4.99796e-06, J_valid=1.000, J2=1.79e-06, leak=1.86e-16, rank2=0.999999, kgap=3.33e-09, Jsign=2.67e-14, Ksign=3.44e-14, axis_cos=0.747, plane_rev=2.67e-14, first=0.00e+00
  symmetrized_birth | random_same_level_triangle | aging | full: count=80, K=4.61437e-09, omega=3.26285e-09, J_valid=0.537, J2=0.00136, leak=1.67e-16, rank2=0.962953, kgap=4.51e-09, Jsign=6.72e-14, Ksign=2.92e-12, axis_cos=0.000, plane_rev=6.72e-14, first=0.00e+00
  symmetrized_birth | random_same_level_triangle | handoff | full: count=80, K=4.70753e-09, omega=3.32872e-09, J_valid=0.575, J2=0.00185, leak=1.95e-16, rank2=0.97851, kgap=3.35e-09, Jsign=7.6e-14, Ksign=2.93e-12, axis_cos=0.000, plane_rev=7.6e-14, first=0.00e+00
  symmetrized_birth | random_same_level_triangle | live | full: count=80, K=3.23043e-08, omega=2.28426e-08, J_valid=0.750, J2=0.00109, leak=2.04e-16, rank2=0.987949, kgap=3.32e-09, Jsign=1.87e-12, Ksign=1.76e-11, axis_cos=0.001, plane_rev=1.87e-12, first=0.00e+00
  symmetrized_birth | sibling_triangle | aging | full: count=121, K=4.22534e-11, omega=2.98776e-11, J_valid=0.000, J2=nan, leak=nan, rank2=0.966505, kgap=3.48e-09, Jsign=nan, Ksign=2.99e-12, axis_cos=0.000, plane_rev=nan, first=0.00e+00
  symmetrized_birth | sibling_triangle | handoff | full: count=121, K=4.39066e-11, omega=3.10466e-11, J_valid=0.000, J2=nan, leak=nan, rank2=0.967976, kgap=3.28e-09, Jsign=nan, Ksign=2.87e-12, axis_cos=0.000, plane_rev=nan, first=0.00e+00
  symmetrized_birth | sibling_triangle | live | full: count=121, K=2.33912e-10, omega=1.65401e-10, J_valid=0.992, J2=0.012, leak=1.9e-16, rank2=0.993936, kgap=3.52e-09, Jsign=1.47e-11, Ksign=1.86e-11, axis_cos=0.000, plane_rev=1.47e-11, first=0.00e+00

REAL-GROWTH LIVE FULL VS COMMUTING ABLATIONS
  parent_fan_triangle | diagonal: K_remain=0, Jvalid_full=1.000, Jvalid_abl=0.000, full_K=0.00185628, abl_K=0
  parent_fan_triangle | trace_scalar: K_remain=0, Jvalid_full=1.000, Jvalid_abl=0.000, full_K=0.00185628, abl_K=0
  parent_fan_triangle | common_mean_diagonal: K_remain=3.82e-14, Jvalid_full=1.000, Jvalid_abl=0.000, full_K=0.00185628, abl_K=7.09332e-17
  random_same_level_triangle | diagonal: K_remain=0, Jvalid_full=1.000, Jvalid_abl=0.000, full_K=1.71109e-05, abl_K=0
  random_same_level_triangle | trace_scalar: K_remain=0, Jvalid_full=1.000, Jvalid_abl=0.000, full_K=1.71109e-05, abl_K=0
  random_same_level_triangle | common_mean_diagonal: K_remain=3.05e-12, Jvalid_full=1.000, Jvalid_abl=0.000, full_K=1.71109e-05, abl_K=5.21494e-17
  sibling_triangle | diagonal: K_remain=0, Jvalid_full=1.000, Jvalid_abl=0.000, full_K=2.95477e-06, abl_K=0
  sibling_triangle | trace_scalar: K_remain=0, Jvalid_full=1.000, Jvalid_abl=0.000, full_K=2.95477e-06, abl_K=0
  sibling_triangle | common_mean_diagonal: K_remain=1.59e-11, Jvalid_full=1.000, Jvalid_abl=0.000, full_K=2.95477e-06, abl_K=4.71005e-17

INTERPRETATION RULE
  A useful local J-candidate carrier requires K_norm and omega above threshold,
  rank2_balance near 1, kernel_gap near 0, J2_plane_residual near 0,
  orientation_J_sign_residual near 0 under reversed orientation, and collapse under
  diagonal/trace/common-commuting reductions.  These are necessary diagnostic gates,
  not yet a derived complex algebra or global J field.
```

## Output files

- `J_candidate_face_rows.csv`
- `J_candidate_summary_main.csv`
- `J_candidate_summary_by_face_level.csv`
- `J_candidate_summary_by_source.csv`
- `J_candidate_ablation_summary.csv`
- `model_summaries.csv`
- `SUMMARY.txt`

## Reading rule

Positive diagnostic signs are:

```text
real_growth / parent_fan_triangle / live / full:
  K_norm, omega nonzero
  frac_J_valid near 1
  rank2_balance near 1
  kernel_gap near 0
  J2_plane_residual near 0
  orientation_J_sign_residual near 0
  diagonal/trace/common-commuting ablations collapse the signal
```

Even if these hold, the output is still a local skew-sector J-candidate, not a
reelle *-Algebra, not a C*-completion, and not a GNS/AQFT net.
