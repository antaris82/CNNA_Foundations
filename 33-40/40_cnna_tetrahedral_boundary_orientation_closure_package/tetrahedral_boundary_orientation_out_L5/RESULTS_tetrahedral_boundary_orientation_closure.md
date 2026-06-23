# CNNA tetrahedral boundary orientation closure

max_level: 5
rows: 3000

## Primary reading
- count: 40
- K_fan_mean: 0.0015280509987882038
- K_base: 2.322140923881806e-06
- base_to_fan_ratio: 0.005377896161562551
- fan_sum_rel: 0.001792632061256006
- boundary_sum_rel: 6.813743567063538e-15
- boundary_gain_vs_fan: 0.001792632061249192
- axis_side_abs_cos: 0.9999354623829285
- axis_side_chirality_abs: 1.4150350746522124e-06
- axis_base_side_chirality_abs: 0.00016914742286436394
- Jfixed: 0.544876930747652
- Jpm: 0.05276698115903231
- fan_wilson: 5.362756867760475e-08
- boundary_wilson: 2.3635049646042634e-13
- fan_wilson_normed: 0.023092734659638486
- boundary_wilson_normed: 0.007904940183174392

## Gate interpretation
- `K_base/base_to_fan_ratio` tests whether the sibling base participates in the tetrahedral boundary.
- `boundary_sum_rel < fan_sum_rel` would indicate base-assisted oriented boundary closure.
- nonzero `axis_*_chirality_abs` would indicate genuine 3D/chiral orientation rather than parallel flat axes.
- diagonal/trace/common-mean reductions should kill the noncommutative K carrier.
- symmetrized birth should suppress the dynamical amplitude if the carrier is growth-asymmetry dependent.
