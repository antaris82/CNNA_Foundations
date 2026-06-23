# SOURCE AUDIT

The previous package corrected the magnitude-only quadrature area test by measuring
signed area.  The signed spatial area did not flip coherently.  This package therefore
moves from spatial area to real pair algebra:

- pair exchange/conjugation C;
- oriented cochain exchange J_pair;
- their algebraic residuals;
- whether the actual derived Q/P channels are locked by that local J_pair;
- signed real commutator diagnostics.

This tests the user's hypothesis that the missing structure is closer to an
erzeuger/vernichter or quadrature-adjunction split than to a spatial face orientation.
