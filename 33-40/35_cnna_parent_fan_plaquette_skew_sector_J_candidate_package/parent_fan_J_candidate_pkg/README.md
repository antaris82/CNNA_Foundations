# CNNA parent-fan plaquette skew-sector J-candidate diagnostic

This package is the first deliberately J-near test after the simplicial DtN
plaquette-frustration and oriented-sign/coboundary audits.

It still does **not** derive a complex Hilbert space, C*-norm, GNS
representation, AQFT net, or global complex algebra.  It asks the narrower gate
question:

> Does the robust parent-fan plaquette commutator carry a stable rank-2 real
> skew sector that behaves locally like a rotation plane?

The candidate is not freely set.  It is extracted from

```text
K_abc = [A_ab, A_bc]
```

where `A_ab = S_b - S_a` and `S_a,S_b,S_c` are full symmetric DtN vertex
operators on a provenance-generated oriented 2-simplex.

For nonzero 3x3 skew `K`, the test extracts its rank-2 image plane and
normalizes `J = K / omega` on that plane.  The nontrivial gate is therefore not
the linear algebra fact that a nonzero 3x3 skew matrix defines a plane rotation,
but whether the robust carrier survives controls and shows correct orientation
behavior.

## Run

```bash
python3 test_parent_fan_plaquette_skew_sector_J_candidate.py --max-level 6 --outdir parent_fan_plaquette_J_candidate_out_L6
```

## Outputs

- `J_candidate_face_rows.csv`
- `J_candidate_summary_main.csv`
- `J_candidate_summary_by_face_level.csv`
- `J_candidate_summary_by_source.csv`
- `J_candidate_ablation_summary.csv`
- `model_summaries.csv`
- `SUMMARY.txt`
- `RESULTS_parent_fan_plaquette_skew_sector_J_candidate.md`

## Positive diagnostic signs

A useful local J-candidate carrier requires:

- nonzero `K_norm` and `omega`,
- `frac_J_valid` near 1,
- `rank2_balance` near 1,
- `kernel_gap` near 0,
- `J2_plane_residual` near 0,
- `orientation_J_sign_residual` near 0 under reversed orientation,
- collapse under diagonal/trace/common-commuting ablations.

Even if this holds, it is still a local skew-sector J-candidate only.
