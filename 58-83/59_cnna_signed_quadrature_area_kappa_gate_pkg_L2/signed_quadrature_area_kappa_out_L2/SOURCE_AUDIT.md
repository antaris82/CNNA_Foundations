# SOURCE AUDIT

Package 58 found a robust Q/P split, but its `symplectic_area_candidate` was
`np.linalg.norm(Q x P)`.  This package tests Claude's objection directly by replacing
that magnitude-only quantity with signed projections against derived local references:
face birth-order normal, outward normal, and pair-axis.

A positive abs-area with near-zero or non-flipping signed area is not a J/i/symplectic
orientation derivation.  It is only a nonzero two-quadrature magnitude.
