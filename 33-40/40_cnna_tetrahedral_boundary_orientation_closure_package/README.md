# CNNA tetrahedral boundary orientation closure

This package tests whether the full local 3-simplex `parent + three children`
contains an oriented boundary / chirality / closure obstruction that is invisible
in the parent-fan side faces alone.

Run for the default quick scan:

```bash
python3 test_tetrahedral_boundary_orientation_closure.py --max-levels 5 6
```

Primary output directories are `tetrahedral_boundary_orientation_out_L5/` and
`tetrahedral_boundary_orientation_out_L6/`.
