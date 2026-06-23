# Notes

The test intentionally defines edge generators by exact vertex-potential differences:

```text
A_ab = S_b - S_a
```

Therefore

```text
A_ab + A_bc + A_ca = 0
```

is forced. Any surviving plaquette signal is second-order/noncommutative:

```text
[A_ab,A_bc]
```

The `common_mean_diagonal` ablation projects all three vertex operators into a shared eigenbasis of their mean. It is not a theorem-level best commuting approximation, but it is a strong diagnostic control: if the full signal is absent after common-basis diagonalization, the full signal is carried by noncommuting offdiagonal DtN structure rather than scalar/height-like coboundary data.

`commutator_sign_residual ~ 0` and `plaquette_skew_sign_residual ~ 0` mean forward and reversed orientations are negatives of each other. Axial cosine values are meaningful only when the underlying signal is not numerically tiny.
