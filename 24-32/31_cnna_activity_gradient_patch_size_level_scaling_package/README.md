# CNNA activity-gradient + patch-size star-closure scaling

This package continues the real local operator-system diagnostic chain after the local patch closure test.

It tests two linked hypotheses:

1. **Growth activity is local in the instantaneous/marginal sense.**  The largest last-shell conductance updates occur near the active parent frontier; the root receives the smallest shell-normalized marginal update.
2. **A real local *-operator system, if present, may be a local-net/limit phenomenon.**  Therefore closure should be tested on growth-defined local patches of increasing size/level, not on one isolated Double-History pair.

No `i`, no `J`, no `J² = -I`, no complex phase, no C*-norm, and no GNS representation are inserted.

Run:

```bash
python test_activity_gradient_and_patch_scaling_star_closure.py \
  --levels 4 5 6 \
  --activity-level 6 \
  --patch-sizes 2 3 \
  --degree 2 \
  --cases real_growth \
  --max-patches-per-control 1 \
  --word-cap 100 \
  --mult-sample-cap 40 \
  --out activity_gradient_patch_scaling_out
```

Main output:

```text
activity_gradient_patch_scaling_out/RESULTS_activity_gradient_patch_scaling.md
```
