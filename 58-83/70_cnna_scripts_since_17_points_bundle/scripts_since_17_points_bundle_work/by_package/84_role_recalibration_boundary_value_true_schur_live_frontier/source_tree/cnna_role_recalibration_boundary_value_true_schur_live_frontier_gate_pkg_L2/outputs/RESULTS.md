# RESULTS — CNNA true Schur/DtN live frontier-distance boundary-value audit

## Purpose

This package continues the established script-1/script-2 growth model and the true Schur/DtN birth+relaxation calculation.  It adds the realistic infinite-growth constraint raised in the discussion:

```text
In the realistic model growth is unbounded.
Very old/interior cuts are far from the active growth front.
For such cuts, direct discrete node-birth shocks should become small,
while fixed-topology live Schur/DtN relaxation remains the relevant measurement layer.
The newest children are arbitrarily far from the root in the infinite-depth limit.
```

The finite L2/L3/L4 runs are not the infinite limit.  They are a controlled approximant audit of the separation:

```text
birth_record shock:
  DtN jump caused by a new boundary role / UV port.

live measurement gap:
  DtN drift caused by relaxation at fixed topology and fixed boundary ports.

front distance:
  graph distance from the cut to the newborn/front event.
```

No J, i, Hodge, star, positivity, Q/P target, or delta-beta is used.

## L4 headline table

| variant | level | nodes | events | front attenuation gate | old interior live gate | fixed-topology relax gate | log birth slope vs distance | log child-coupling slope | far/near birth shock | old interior live/birth | advanced leakage | child own UV-tail |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| real_growth_linear_true_schur_live_frontier | L4 | 121 | 120 | 1 | 1 | 1 | -0.305 | -0.479 | 0.558 | 0.244 | 0.000 | 0.000 |
| log_growth_true_schur_live_frontier | L4 | 121 | 120 | 1 | 1 | 1 | -0.268 | -0.392 | 0.552 | 0.229 | 0.000 | 0.000 |
| saturating_growth_true_schur_live_frontier | L4 | 121 | 120 | 1 | 1 | 1 | -0.147 | -0.245 | 0.613 | 0.962 | 0.000 | 0.000 |
| strict_sym_true_schur_live_frontier_control | L4 | 121 | 120 | 0 | 0 | 0 | 0.000 | 0.000 | 0 | 0 | 0.000 | 0.000 |

## Interpretation

The test separates two effects that should not be conflated:

```text
Birth role update:
  A newborn changes topology and boundary-port roles.
  Parent/ancestor cuts gain a new UV-tail port.

Live relaxation:
  No new node is born and boundary ports are fixed.
  Conductances and true Schur/DtN responses continue to relax.
```

The finite-frontier audit supports the qualitative infinite-growth reading: old/interior cuts are not well described as repeatedly receiving large local birth shocks.  They retain nonzero live Schur/DtN gaps under fixed-topology relaxation.  In finite L4 this is only an approximant signal, but it correctly separates record and live layers.

The strict-sym control is expected to kill the response/live layer: topology can still grow, but with `alpha_env=br_ancestor=br_sibling=0` there is no directed conductance response, no live gap, and no meaningful boundary-value polarity.

## Consequence for the CNNA interpretation

The dynamic boundary-role layer should now be read as two-layer and scale-dependent:

```text
Near active frontier:
  discrete birth events dominate;
  new UV ports and boundary-role changes are visible.

Old/interior / far from frontier:
  discrete births are increasingly remote;
  the effective measurement is mostly live Schur/DtN relaxation.

Infinite-growth idealization:
  the newest children are infinitely far from the root;
  root/interior physics should be live-response dominated,
  not a sequence of local birth shocks.
```

This strengthens the record/live split:

```text
Record layer:
  immutable birth/event/boundary-role history.

Live layer:
  continually relaxing Schur/DtN response field on the already-grown network.
```

## Limits

- L2/L3/L4 are finite approximants, not a proof of the infinite-depth limit.
- The direct birth-shock attenuation is measured by finite front distance and is sensitive to the chosen deterministic growth/relaxation rule.
- DtN matrices are true Schur complements of the real conductance Laplacian; the relaxation law is still a deterministic model rule within script-1/script-2 quantities.
- The package does not claim J, i, spin, star, positivity, Q/P compatibility, or modular structure.

## Next test

`test_true_schur_live_boundary_polarity_flip_gate.py`

Use the separated record/live rows and test whether the live Schur/DtN boundary polarity is stable or flips under:

```text
1. κ birth-order mirror,
2. reverse/advanced response control,
3. longitudinal/root-front axis flip,
4. strict_sym control.
```

The primary gate should remain boundary-value polarity, not J/QP-lock.
