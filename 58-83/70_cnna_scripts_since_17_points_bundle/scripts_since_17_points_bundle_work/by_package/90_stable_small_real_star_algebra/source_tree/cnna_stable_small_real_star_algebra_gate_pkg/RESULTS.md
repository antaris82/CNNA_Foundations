# CNNA stable small real #/*-algebra-like operator-family gate

## Purpose

This package tests the next requested step after the depth-scaling check: whether bridge-/passivity-positive Record→Live Schur/DtN relaxation rows generate a small real operator family stable under:

```text
#_G  = record-DtN metric adjoint,
products,
commutators.
```

This is **not** a physical `*`-algebra and not a C*-algebra claim. The metric adjoint exists by linear algebra and is not counted as success unless the generated finite operator family also closes approximately under products and commutators.

No `J`, `i`, Hodge, physical Hilbert adjoint, C*-norm, Q/P target or `delta_beta` gate is used.

## Compute regime

- L3/L4: low/full-style runs with complete early frontier before sampling.
- L5/L6: deterministic sampled frontier approximants with smaller cap for feasibility. These are trend probes, not full-tree averages.

## Summary at L4 and L6

| run | level | variant | valid rows | candidate+ | weak # family pass | strong stable pass | product residual | commutator residual | commutator norm |
|---|---:|---|---:|---:|---:|---:|---:|---:|---:|
| full_or_low_L3_L4 | L4 | real_growth_linear_star_candidate | 198 | 0.354 | 0.106 | 0.000 | 0.010 | 0.402 | 0.002 |
| full_or_low_L3_L4 | L4 | saturating_growth_star_candidate | 198 | 0.995 | 0.116 | 0.000 | 0.053 | 0.526 | 0.003 |
| full_or_low_L3_L4 | L4 | strict_symmetrized_response_star_control | 198 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| sampled_L5_L6 | L6 | real_growth_linear_star_candidate | 237 | 0.418 | 0.084 | 0.000 | 0.019 | 0.466 | 0.002 |
| sampled_L5_L6 | L6 | saturating_growth_star_candidate | 237 | 0.996 | 0.084 | 0.000 | 0.047 | 0.535 | 0.003 |
| sampled_L5_L6 | L6 | strict_symmetrized_response_star_control | 237 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |

## Interpretation

The live semigroup supplies a real adjunction/passivity precursor, but a robust small real `*`-algebra-like family is **not yet uniformly established**.

The strong gate remains hard because product/commutator closure is the nontrivial part. A `#_G` adjoint by itself is not enough.

The result therefore remains ladder-correct:

```text
Record/Live semigroup polarity: yes.
Passivity-/adjunction-like precursor: partial yes.
Stable small real *-algebra-like family: not robust yet.
```

## Next test

The next step should not claim a C*-structure. It should identify whether the missing closure term is a coarse-graining issue:

```text
test_coarse_grained_live_star_family_closure_gate.py
```

Use generation/depth grouped operator averages or principal stable modes instead of raw per-relax-step operators, then retest #, product and commutator closure.
