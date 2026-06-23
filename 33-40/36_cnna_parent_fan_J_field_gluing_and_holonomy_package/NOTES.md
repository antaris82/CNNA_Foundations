# Notes

## Diagnostic status

This is a gluing/holonomy diagnostic after the local parent-fan J-candidate test.

The previous test established a robust local face-level carrier:

```text
real_growth / parent_fan_triangle / live / full DtN
K_norm nonzero
J^2_plane residual near zero
orientation reversal sends J -> -J
commuting reductions kill K
```

This package asks whether the three parent-fan faces around one parent

```text
(p,c1,c2), (p,c2,c3), (p,c3,c1)
```

carry mutually coherent J-planes.

## Reading rule

Positive gluing signs:

- valid faces near 1,
- small per-face J2 residual,
- small plane projector residual across the three faces,
- small signless J residual,
- high absolute axis cosine,
- small J-pair commutator,
- collapse under diagonal/trace/common-commuting reductions.

The signed fixed-order J residual is expected to expose whether a further sign-lift/orientation bundle is still needed.  A small signless residual but large fixed-sign residual means: the local planes glue, but the ±J sign is not globally locked by this test alone.

## Main limitation

This is not a proof of a global J field.  It tests local parent-fan coherence.  The next gate is a sign-lift/cocycle test: can the ±J choices be made consistently over the growth-defined face net, or does a nontrivial Z2 obstruction remain?
