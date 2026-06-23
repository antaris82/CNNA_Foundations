# Notes

## What this test joins

- Growth tests 1/2: sequential birth and backreaction.
- Geometry tests 3/9: 2-simplices and the failure of minimal-axis geometry alone.
- DtN tests 22/23: dynamic record/live DtN operators.
- Operator-closure tests 24+: node/edge provenance alone does not close a *-structure.
- Plaquette tests: full DtN matrices on oriented 2-simplices produce a noncommuting second-order face effect.

## Why the test is not tautological

Every nonzero real 3x3 skew matrix has a rank-2 rotation plane.  Therefore the
interesting question is not whether `K` can be normalized to something with
`J^2=-I` on its image plane.  The interesting question is whether a nonzero,
controlled, provenance/geometric `K` exists in the right carrier.

The audit gates are:

1. parent-fan/full/live carrier much stronger than sibling/random/symmetrized controls;
2. orientation reversal sends `K` and `J` to their negatives;
3. diagonal/trace/common-commuting reductions kill the signal;
4. handoff/aging signals vanish under no-backreaction as expected;
5. no global complex algebra or *-closure is inferred from the local plane alone.

## Main current finding

The real-growth parent-fan live/full carrier gives a stable local rank-2 skew
sector with `J2_plane_residual` near numerical zero and clean orientation flip.
Diagonal/scalar/common-commuting controls collapse.

This is evidence for a local face-level J-candidate carrier, not yet for a
canonical global J field or a complex local algebra.
