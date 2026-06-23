# Results: oriented plaquette sign and coboundary controls

This test audits the first simplicial DtN plaquette-frustration signal.

It still does **not** set `i`, `J`, a Hilbert space, a C*-norm, or a *-algebra.

## Gate questions

1. Orientation reversal: `(a,b,c)` versus `(a,c,b)` should flip the commutator/skew sign.
2. Coboundary controls: first-order edge circulation is exact and telescopes; diagonal, trace-scalar, and common-mean-diagonal ablations should remove any scalar/commuting explanation.
3. Face separation: sibling triangles and parent-fan triangles are kept as distinct geometric carriers.

## Current run

```text
CNNA ORIENTED PLAQUETTE SIGN AND COBOUNDARY CONTROLS
max_level=6, mode=linear, step=0.04

MODEL SUMMARIES
  real_growth: nodes=1093, completed=364, sibling=121, parent_fan=363, random=80
  symmetrized_birth: nodes=1093, completed=364, sibling=121, parent_fan=363, random=80
  no_backreaction: nodes=1093, completed=364, sibling=121, parent_fan=363, random=80

SELECTED FULL-DTN FACE SUMMARIES
  no_backreaction | parent_fan_triangle | aging | full: count=363, comm=0, comm_sign=0, comm_cos=0.000, plaq=1.33538e-16, skew=0, skew_sign=0, inv=0.000305, complex_frac=0.000, first=0.00e+00
  no_backreaction | parent_fan_triangle | handoff | full: count=363, comm=9.17998e-20, comm_sign=0, comm_cos=0.000, plaq=1.36639e-16, skew=6.22229e-22, skew_sign=2.42e-20, inv=0.000311, complex_frac=0.000, first=0.00e+00
  no_backreaction | parent_fan_triangle | live | full: count=363, comm=0.0036013, comm_sign=2.48e-15, comm_cos=1.000, plaq=2.88107e-06, skew=2.88106e-06, skew_sign=2.62e-10, inv=2.7e-09, complex_frac=0.997, first=0.00e+00
  no_backreaction | random_same_level_triangle | aging | full: count=80, comm=0, comm_sign=0, comm_cos=0.000, plaq=1.22754e-16, skew=0, skew_sign=0, inv=0.000273, complex_frac=0.000, first=0.00e+00
  no_backreaction | random_same_level_triangle | handoff | full: count=80, comm=0, comm_sign=0, comm_cos=0.000, plaq=1.22754e-16, skew=0, skew_sign=0, inv=0.000273, complex_frac=0.000, first=0.00e+00
  no_backreaction | random_same_level_triangle | live | full: count=80, comm=6.36467e-06, comm_sign=7.21e-13, comm_cos=0.708, plaq=5.09192e-09, skew=5.09184e-09, skew_sign=2.27e-07, inv=6.35e-06, complex_frac=0.100, first=0.00e+00
  no_backreaction | sibling_triangle | aging | full: count=121, comm=0, comm_sign=0, comm_cos=0.000, plaq=5.4072e-17, skew=0, skew_sign=0, inv=0.000102, complex_frac=0.000, first=0.00e+00
  no_backreaction | sibling_triangle | handoff | full: count=121, comm=3.39936e-21, comm_sign=0, comm_cos=0.000, plaq=6.44528e-17, skew=2.71949e-24, skew_sign=9.77e-22, inv=0.000123, complex_frac=0.000, first=0.00e+00
  no_backreaction | sibling_triangle | live | full: count=121, comm=4.16893e-06, comm_sign=3.21e-13, comm_cos=0.896, plaq=3.3352e-09, skew=3.33518e-09, skew_sign=1.91e-08, inv=4.18e-07, complex_frac=0.000, first=0.00e+00
  real_growth | parent_fan_triangle | aging | full: count=363, comm=0.000110915, comm_sign=2.16e-15, comm_cos=0.999, plaq=8.87419e-08, skew=8.87371e-08, skew_sign=1.1e-09, inv=2.5e-08, complex_frac=1.000, first=0.00e+00
  real_growth | parent_fan_triangle | handoff | full: count=363, comm=0.000114368, comm_sign=2.05e-15, comm_cos=0.999, plaq=9.15049e-08, skew=9.14999e-08, skew_sign=1.21e-09, inv=2.64e-08, complex_frac=1.000, first=0.00e+00
  real_growth | parent_fan_triangle | live | full: count=363, comm=0.00185628, comm_sign=6.05e-15, comm_cos=1.000, plaq=1.48509e-06, skew=1.48506e-06, skew_sign=3.79e-10, inv=4.93e-09, complex_frac=1.000, first=0.00e+00
  real_growth | random_same_level_triangle | aging | full: count=80, comm=2.13374e-06, comm_sign=5.18e-14, comm_cos=0.427, plaq=1.70699e-09, skew=1.70699e-09, skew_sign=1.29e-05, inv=0.000215, complex_frac=0.000, first=0.00e+00
  real_growth | random_same_level_triangle | handoff | full: count=80, comm=2.15391e-06, comm_sign=4.75e-14, comm_cos=0.430, plaq=1.72313e-09, skew=1.72313e-09, skew_sign=1.21e-06, inv=2.44e-05, complex_frac=0.000, first=0.00e+00
  real_growth | random_same_level_triangle | live | full: count=80, comm=1.71109e-05, comm_sign=3.56e-13, comm_cos=0.805, plaq=1.36891e-08, skew=1.36889e-08, skew_sign=3.52e-08, inv=8.91e-07, complex_frac=0.388, first=0.00e+00
  real_growth | sibling_triangle | aging | full: count=121, comm=2.61871e-07, comm_sign=5.36e-14, comm_cos=0.078, plaq=2.09497e-10, skew=2.09497e-10, skew_sign=2.28e-05, inv=0.000474, complex_frac=0.000, first=0.00e+00
  real_growth | sibling_triangle | handoff | full: count=121, comm=2.66732e-07, comm_sign=4.08e-14, comm_cos=0.080, plaq=2.13386e-10, skew=2.13386e-10, skew_sign=1.18e-05, inv=0.000259, complex_frac=0.000, first=0.00e+00
  real_growth | sibling_triangle | live | full: count=121, comm=2.95477e-06, comm_sign=5.53e-13, comm_cos=0.804, plaq=2.36387e-09, skew=2.36385e-09, skew_sign=2.99e-08, inv=5.63e-07, complex_frac=0.000, first=0.00e+00
  symmetrized_birth | parent_fan_triangle | aging | full: count=363, comm=4.00871e-06, comm_sign=1.99e-14, comm_cos=0.660, plaq=3.20718e-09, skew=3.20708e-09, skew_sign=1.44e-07, inv=4.25e-06, complex_frac=0.000, first=0.00e+00
  symmetrized_birth | parent_fan_triangle | handoff | full: count=363, comm=4.03911e-06, comm_sign=2.06e-14, comm_cos=0.661, plaq=3.2315e-09, skew=3.2314e-09, skew_sign=1.35e-07, inv=4.28e-06, complex_frac=0.000, first=0.00e+00
  symmetrized_birth | parent_fan_triangle | live | full: count=363, comm=7.06818e-06, comm_sign=3.44e-14, comm_cos=0.747, plaq=5.65478e-09, skew=5.65467e-09, skew_sign=5.42e-08, inv=1.16e-06, complex_frac=0.000, first=0.00e+00
  symmetrized_birth | random_same_level_triangle | aging | full: count=80, comm=4.61437e-09, comm_sign=2.92e-12, comm_cos=0.000, plaq=3.69153e-12, skew=3.69149e-12, skew_sign=6.99e-05, inv=0.00168, complex_frac=0.000, first=0.00e+00
  symmetrized_birth | random_same_level_triangle | handoff | full: count=80, comm=4.70753e-09, comm_sign=2.93e-12, comm_cos=0.000, plaq=3.76603e-12, skew=3.76602e-12, skew_sign=6.42e-05, inv=0.000865, complex_frac=0.000, first=0.00e+00
  symmetrized_birth | random_same_level_triangle | live | full: count=80, comm=3.23043e-08, comm_sign=1.76e-11, comm_cos=0.001, plaq=2.58434e-11, skew=2.58434e-11, skew_sign=3.25e-05, inv=0.0007, complex_frac=0.000, first=0.00e+00
  symmetrized_birth | sibling_triangle | aging | full: count=121, comm=4.22534e-11, comm_sign=2.99e-12, comm_cos=0.000, plaq=3.38338e-14, skew=3.38084e-14, skew_sign=0.000111, inv=0.00208, complex_frac=0.000, first=0.00e+00
  symmetrized_birth | sibling_triangle | handoff | full: count=121, comm=4.39066e-11, comm_sign=2.87e-12, comm_cos=0.000, plaq=3.51566e-14, skew=3.51326e-14, skew_sign=0.000105, inv=0.00203, complex_frac=0.000, first=0.00e+00
  symmetrized_birth | sibling_triangle | live | full: count=121, comm=2.33912e-10, comm_sign=1.86e-11, comm_cos=0.000, plaq=1.87134e-13, skew=1.87129e-13, skew_sign=9.14e-05, inv=0.00178, complex_frac=0.000, first=0.00e+00

REAL-GROWTH LIVE COBOUNDARY/COMMUTING ABLATIONS
  parent_fan_triangle | diagonal: comm_remain=0, skew_remain=0, full_comm=0.00185628, abl_comm=0
  parent_fan_triangle | trace_scalar: comm_remain=0, skew_remain=0, full_comm=0.00185628, abl_comm=0
  parent_fan_triangle | common_mean_diagonal: comm_remain=4.67e-14, skew_remain=5.85e-11, full_comm=0.00185628, abl_comm=8.67648e-17
  random_same_level_triangle | diagonal: comm_remain=0, skew_remain=0, full_comm=1.71109e-05, abl_comm=0
  random_same_level_triangle | trace_scalar: comm_remain=0, skew_remain=0, full_comm=1.71109e-05, abl_comm=0
  random_same_level_triangle | common_mean_diagonal: comm_remain=3.74e-12, skew_remain=6e-09, full_comm=1.71109e-05, abl_comm=6.39421e-17
  sibling_triangle | diagonal: comm_remain=0, skew_remain=0, full_comm=2.95477e-06, abl_comm=0
  sibling_triangle | trace_scalar: comm_remain=0, skew_remain=0, full_comm=2.95477e-06, abl_comm=0
  sibling_triangle | common_mean_diagonal: comm_remain=2.09e-11, skew_remain=3.93e-08, full_comm=2.95477e-06, abl_comm=6.16584e-17

INTERPRETATION
  Good orientation behavior means commutator_sign_residual and skew_sign_residual are near 0,
  and axial_cos_reversal is near +1 when comparing forward with negative reversed orientation.
  Diagonal/trace/common-basis ablations test whether the signal is merely scalar/commuting coboundary data.
  first_order_closure is forced to telescope and should remain numerical zero.
  This is still not a J-test; it audits the face-level noncommutative carrier before any J extraction.
```

## Output files

- `oriented_pair_rows.csv`
- `orientation_summary_main.csv`
- `orientation_summary_by_face_level.csv`
- `orientation_summary_by_source.csv`
- `coboundary_ablation_summary.csv`
- `model_summaries.csv`
- `SUMMARY.txt`

## Reading rule

A useful positive sign audit is:

```text
commutator_sign_residual ~ 0
plaquette_skew_sign_residual ~ 0
commutator_axial_cos_reversal ~ 1
first_order_closure ~ 0
full DtN signal >> diagonal/trace/common-commuting ablations
```

If these fail, the prior plaquette signal is not yet reliable enough for any J-sector projection.
