# SOURCE AUDIT

Previous gate found a local pair algebra:

```text
C_pair^2 = +I
J_pair^2 = -I
C_pair J_pair C_pair = -J_pair
```

but the actual raw Q/P channels were not locked by this J.  This package searches for
alignment only among derived channels already present in the pair data, not among arbitrary
rotations or user-chosen complex structures.
