# RESULTS: local patch operator-system closure

## Purpose

This diagnostic enlarges the carrier from a single Double-History pair to a same-suffix local patch of Record/Live boundary cells. It asks whether DtN-adjoint and finite product closure become visible only at patch level.

## Parameters

```json
{
  "levels": [
    4
  ],
  "cases": [
    "real_growth"
  ],
  "tiers": [
    "base_patch"
  ],
  "patch_size": 3,
  "degree": 2,
  "word_cap": 250,
  "max_patches_per_control": 2,
  "operator_mode": "triangular_handoff",
  "longitudinal_mode": "triangular_record_live",
  "metric_source": "record_live_block"
}
```

## Primary summaries

### real_growth:L4:same_suffix_multi_history_patch:base_patch:patch3:degree2

```text
count                                      2
patch_dim_mean                             18.0
seed_count_mean                            28.0
word_space_dim_mean                        79.0
word_space_fraction_full_mean              0.24382716049382716
handoff_T_norm_mean_mean                   0.11152776454161809
patch_comm_C_norm_mean_mean                0.1657010400649197
patch_comm_C_rel_mean_mean                 0.2745094347504924
star_basis_rel_mean_mean                   0.46568316649475233
star_basis_rel_p95_mean                    0.5221963110742732
star_basis_rel_lt_0p25_fraction_mean       0.5
star_seed_rel_mean_mean                    0.2038877921222176
star_seed_rel_lt_0p25_fraction_mean        0.5714285714285714
mult_left_rel_mean_mean                    0.39781668569878603
mult_left_rel_lt_0p25_fraction_mean        0.3466666666666667
patch_metric_global_cond_mean              15.402090862666228
```

## Interpretation

A meaningful positive pre-* gate would require star residuals and multiplication residuals to decrease at patch level without saturating the full matrix algebra. If residuals remain large, the current patch is still too small, the generated operators are still incomplete, or *-closure is a limiting/local-net phenomenon rather than a finite-patch property.