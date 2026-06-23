# RESULTS — dual pairing two-edge assembly gate

## Comparative table

| variant | candidates | role A beta/QP | role B C/kappa | single-pair all | two-edge same-scan | connected | strong | used Δβ? |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| real_growth | 206 | 80 | 3 | 0 | 38 | 38 | 26 | False |
| strict_symmetrized_control | 0 | 0 | 0 | 0 | 0 | 0 | 0 | False |
| no_backreaction | 236 | 103 | 2 | 0 | 22 | 22 | 16 | False |

## Interpretation

This audit tests whether the roles that did not coincide on a single pair can be distributed across two same-scan pair candidates:

```text
Pair A: beta2/QP carrier
Pair B: C-lock/kappa-flip carrier
Assembly: same scan plus geometric context
```

The script does not introduce a new growth rule and does not use beta2 as a move decision.  `delta_beta2` is an audit label inherited from the candidate evaluation.

A positive two-edge signal requires nonzero same-scan assemblies; a stronger signal requires connected or strong face/edge context.  If connected assemblies exist, the single-pair obstruction is not the end of the line; the relevant object is at least a two-edge assembly.
