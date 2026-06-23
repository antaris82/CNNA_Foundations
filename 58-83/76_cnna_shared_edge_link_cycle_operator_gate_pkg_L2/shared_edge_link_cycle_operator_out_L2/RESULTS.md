# RESULTS — shared edge link-cycle operator decision gate

## Comparative table

| option | variant | beta | pairs | asm | pair harm | Q harm | P harm | pair J-lock | signed | union lock | edge-if lock | link-cycle best | active link lock | link J2 | active J2 | nonzero circ | kappa flips | decision pass | abort marker | used dBeta? |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| A_to_B_rescan_strong_reuseB | real_growth | (1,0,4,0) | 4 | 1 | 0.381795 | 0.381795 | 0.248801 | 0.448553 | -0.624664 | 0.543826 | 0.596246 | 0.543826 | 0.5534 | 0.339506 | 0.339946 | 72 | 0 | 0 | True | False |
| A_to_B_rescan_strong_reuseB | strict_symmetrized_control | (1,0,0,0) | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | True | False |
| A_to_B_rescan_strong_reuseB | no_backreaction | (1,0,8,0) | 4 | 2 | 0.405623 | 0.405623 | 0.343367 | 0.339477 | -0.285265 | 0.661108 | 0.689923 | 0.661108 | 0.666166 | 0.471393 | 0.568405 | 144 | 0 | 0 | True | False |

## Interpretation

`link-cycle lock` is a residual.  Smaller is better.  Values around 0.4--0.6 mean the tested operator does not act as a good derived complex structure on the actual Q/P motif space.

This gate is stricter than the previous edge-interface package.  It requires signed circulation and kappa flip in addition to a low Q/P-J residual and projected J^2 behavior.  Merely finding nonzero magnitude or beta2 is not counted as success.

## Stop/continue rule

- If `decision pass > 0`, the shared-edge link is a serious candidate for the missing local alignment operator.
- If `decision pass = 0` and the best non-strict residual remains >= 0.35, this package marks the local link-cycle path as negative for J derivation on the actual Q/P motif space.  The next scientific step should be interpretation/formal obstruction documentation, not automatic escalation to a larger fan.
