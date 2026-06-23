# CNNA tetrahedral boundary orientation closure

max_level: 6
rows: 7599

## Primary reading
- count: 121
- K_fan_mean: 0.0018562815899274571
- K_base: 2.9547692634046493e-06
- base_to_fan_ratio: 0.006617949231204465
- fan_sum_rel: 0.0022059830901373263
- boundary_sum_rel: 6.479906666190294e-15
- boundary_gain_vs_fan: 0.0022059830901308467
- axis_side_abs_cos: 0.9998896851140475
- axis_side_chirality_abs: 1.7206588166308e-06
- axis_base_side_chirality_abs: 0.00012079981968124566
- Jfixed: 0.5398885249443797
- Jpm: 0.041672973206009004
- fan_wilson: 6.823749724783187e-08
- boundary_wilson: 2.987484448132693e-13
- fan_wilson_normed: 0.023092257600960984
- boundary_wilson_normed: 0.007410406375660452

## Gate interpretation
- `K_base/base_to_fan_ratio` tests whether the sibling base participates in the tetrahedral boundary.
- `boundary_sum_rel < fan_sum_rel` would indicate base-assisted oriented boundary closure.
- nonzero `axis_*_chirality_abs` would indicate genuine 3D/chiral orientation rather than parallel flat axes.
- diagonal/trace/common-mean reductions should kill the noncommutative K carrier.
- symmetrized birth should suppress the dynamical amplitude if the carrier is growth-asymmetry dependent.
