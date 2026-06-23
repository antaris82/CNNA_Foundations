# RESULTS — retarded boundary-value selection gate

| variant | beta | assemblies | pair_harm | Q_harm | P_harm | good rows | long-ret rows | causal signal | boundary signal | retarded gate | advanced boundary | best long-ret QP | max halfplane bias | fwd/rev pairs | fwd better | used Δβ? |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| real_growth | (1,0,4,0) | 1 | 0.278528 | 0.278528 | 0.348841 | 175 | 50 | 125 | 9 | 9 | 0 | 0.89726 | 1 | 25 | 15 | False |
| strict_symmetrized_control | (1,0,0,0) | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | False |
| no_backreaction | (1,0,6,0) | 2 | 0.298469 | 0.298469 | 0.279644 | 230 | 76 | 190 | 14 | 14 | 0 | 0.821893 | 1 | 20 | 14 | False |

## Interpretation rule

The test does **not** search for spin or `J`.  It asks whether the growth semi-order itself creates a robust boundary-value/retarded selection analogous to choosing an upper/lower half-plane prescription.

A positive row must be:

- well-conditioned, bounded, non-artifact, and active;
- forward/retarded rather than reverse/advanced or symmetrized;
- longitudinal/radial rather than purely transverse;
- directionally biased by sibling response circulation;
- half-plane biased on the Q/P projection;
- not dependent on `delta_beta` as a decision input.

## Compact JSON

```json
{
  "real_growth": {
    "transfer_row_count": 210,
    "ok_transfer_row_count": 210,
    "good_bounded_nonartifact_active_count": 175,
    "good_retarded_count": 125,
    "good_advanced_count": 25,
    "good_symmetric_count": 25,
    "good_longitudinal_retarded_count": 50,
    "good_transverse_retarded_count": 75,
    "causal_order_signal_count": 125,
    "longitudinal_retarded_signal_count": 45,
    "boundary_value_signal_count": 9,
    "retarded_boundary_value_gate_pass_count": 9,
    "advanced_boundary_signal_count": 0,
    "symmetric_boundary_signal_count": 0,
    "bounded_contraction_good_count": 165,
    "stable_halfplane_good_count": 20,
    "best_longitudinal_retarded_QP_resid": 0.8972598442619204,
    "best_longitudinal_retarded_J2_diag_resid": 0.9904358575715427,
    "max_longitudinal_retarded_halfplane_abs_bias": 0.99999999999975,
    "mean_longitudinal_retarded_circ_bias": 0.9999999999995688,
    "max_longitudinal_retarded_circ_bias": 0.9999999999995689,
    "forward_reverse_pair_count": 25,
    "forward_qp_better_than_reverse_count": 15,
    "forward_reverse_halfplane_opposed_count": 0,
    "best_retarded_row": {
      "bridge_mode": "radial",
      "response_op": "signed_C3_skew",
      "lift_mode": "raw_unbounded_audit__transfer_post_opnorm_bounded",
      "causal_class": "retarded",
      "bridge_axis_class": "longitudinal",
      "QP_transfer_mean_resid": 0.8972598442619204,
      "projected_J2_plus_I_resid": 1.0233879898761369,
      "projected_QP_opnorm": 0.3758279694961683,
      "projected_QP_normality_resid": 0.4386906291167282,
      "halfplane_bias_pos_minus_neg": -0.99999999999975,
      "response_circulation_directional_bias": 0.9999999999995689,
      "signed_mean_response_circulation": 1.1599222839686487,
      "retarded_boundary_value_gate_pass": true,
      "boundary_value_signal_pass": true,
      "longitudinal_retarded_signal_pass": true,
      "context": "shared_edge"
    },
    "first_boundary_pass_row": {
      "bridge_mode": "radial",
      "response_op": "cycle_difference",
      "lift_mode": "raw_unbounded_audit__transfer_post_opnorm_bounded",
      "causal_class": "retarded",
      "bridge_axis_class": "longitudinal",
      "QP_transfer_mean_resid": 0.8994894188245021,
      "projected_J2_plus_I_resid": 1.0334164766682559,
      "projected_QP_opnorm": 0.4824135101623705,
      "projected_QP_normality_resid": 0.5139770159986751,
      "halfplane_bias_pos_minus_neg": -0.99999999999975,
      "response_circulation_directional_bias": 0.9999999999995689,
      "signed_mean_response_circulation": 1.1599222839686487,
      "retarded_boundary_value_gate_pass": true,
      "boundary_value_signal_pass": true,
      "longitudinal_retarded_signal_pass": true,
      "context": "shared_edge"
    }
  },
  "strict_symmetrized_control": {
    "transfer_row_count": 0,
    "ok_transfer_row_count": 0,
    "good_bounded_nonartifact_active_count": 0,
    "good_retarded_count": 0,
    "good_advanced_count": 0,
    "good_symmetric_count": 0,
    "good_longitudinal_retarded_count": 0,
    "good_transverse_retarded_count": 0,
    "causal_order_signal_count": 0,
    "longitudinal_retarded_signal_count": 0,
    "boundary_value_signal_count": 0,
    "retarded_boundary_value_gate_pass_count": 0,
    "advanced_boundary_signal_count": 0,
    "symmetric_boundary_signal_count": 0,
    "bounded_contraction_good_count": 0,
    "stable_halfplane_good_count": 0,
    "best_longitudinal_retarded_QP_resid": 0.0,
    "best_longitudinal_retarded_J2_diag_resid": 0.0,
    "max_longitudinal_retarded_halfplane_abs_bias": 0.0,
    "mean_longitudinal_retarded_circ_bias": 0.0,
    "max_longitudinal_retarded_circ_bias": 0.0,
    "forward_reverse_pair_count": 0,
    "forward_qp_better_than_reverse_count": 0,
    "forward_reverse_halfplane_opposed_count": 0,
    "best_retarded_row": {},
    "first_boundary_pass_row": {}
  },
  "no_backreaction": {
    "transfer_row_count": 420,
    "ok_transfer_row_count": 420,
    "good_bounded_nonartifact_active_count": 230,
    "good_retarded_count": 190,
    "good_advanced_count": 20,
    "good_symmetric_count": 20,
    "good_longitudinal_retarded_count": 76,
    "good_transverse_retarded_count": 114,
    "causal_order_signal_count": 190,
    "longitudinal_retarded_signal_count": 61,
    "boundary_value_signal_count": 14,
    "retarded_boundary_value_gate_pass_count": 14,
    "advanced_boundary_signal_count": 0,
    "symmetric_boundary_signal_count": 2,
    "bounded_contraction_good_count": 200,
    "stable_halfplane_good_count": 51,
    "best_longitudinal_retarded_QP_resid": 0.8218925826458671,
    "best_longitudinal_retarded_J2_diag_resid": 0.9912972001341251,
    "max_longitudinal_retarded_halfplane_abs_bias": 0.99999999999975,
    "mean_longitudinal_retarded_circ_bias": 0.9999999999999809,
    "max_longitudinal_retarded_circ_bias": 0.9999999999999806,
    "forward_reverse_pair_count": 20,
    "forward_qp_better_than_reverse_count": 14,
    "forward_reverse_halfplane_opposed_count": 0,
    "best_retarded_row": {
      "bridge_mode": "radial",
      "response_op": "signed_C3_skew",
      "lift_mode": "raw_unbounded_audit__transfer_post_opnorm_bounded",
      "causal_class": "retarded",
      "bridge_axis_class": "longitudinal",
      "QP_transfer_mean_resid": 0.8632158485929573,
      "projected_J2_plus_I_resid": 1.000594549174874,
      "projected_QP_opnorm": 0.47323288023379556,
      "projected_QP_normality_resid": 0.6342847832024037,
      "halfplane_bias_pos_minus_neg": -0.99999999999975,
      "response_circulation_directional_bias": 0.9999999999999806,
      "signed_mean_response_circulation": 25.819395609712473,
      "retarded_boundary_value_gate_pass": true,
      "boundary_value_signal_pass": true,
      "longitudinal_retarded_signal_pass": true,
      "context": "shared_edge"
    },
    "first_boundary_pass_row": {
      "bridge_mode": "radial",
      "response_op": "signed_C3_skew",
      "lift_mode": "raw_unbounded_audit__transfer_post_opnorm_bounded",
      "causal_class": "retarded",
      "bridge_axis_class": "longitudinal",
      "QP_transfer_mean_resid": 0.9032302235933073,
      "projected_J2_plus_I_resid": 1.0282242336272993,
      "projected_QP_opnorm": 0.4104956348315417,
      "projected_QP_normality_resid": 0.4343005403895866,
      "halfplane_bias_pos_minus_neg": -0.99999999999975,
      "response_circulation_directional_bias": 0.9999999999999806,
      "signed_mean_response_circulation": 25.819395609712473,
      "retarded_boundary_value_gate_pass": true,
      "boundary_value_signal_pass": true,
      "longitudinal_retarded_signal_pass": true,
      "context": "shared_edge"
    }
  }
}
```

## Conservative reading

If the retarded/boundary signal is nonzero but does not order Q/P, then the result is not `i`; it is only a candidate pre-causal boundary-value layer.  This directly tracks the distinction raised by the logarithm/branch-cut discussion: the hidden convention may live in a retarded/advanced prescription rather than in an explicit geometric orientation.
