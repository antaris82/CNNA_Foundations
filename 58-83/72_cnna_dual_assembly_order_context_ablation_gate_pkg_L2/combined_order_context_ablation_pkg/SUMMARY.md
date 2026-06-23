# SUMMARY — dual assembly order/context ablation gate

Model label:
CNNA growing primal simplicial complex with deterministic sequential provenance growth, directed antisymmetric birth-transport operators, and two-pair assembly ablations.

This package compares:

```text
A_to_B_rescan
B_to_A_rescan
stale_same_scan
connected context
strong context
allow_B_reuse_A_faces on/off
```

Anti-smuggling constraints:

```text
no i
no global J
no Hodge
no *
no positivity
no final sym(M)
no delta_beta as a selection criterion
```

The strong-context run was executed for `real_growth` and `no_backreaction`; the connected run includes `strict_symmetrized_control`, which remains fully killed.

## Main pattern

```text
A_to_B_rescan:
  strongest beta2/QP/pair-harmonic carrier.
  no_backreaction gives the best J-lock (~0.209), but signed_birth is small (~0.070).
  real_growth gives the strongest signed_birth (~-0.625), but J-lock is poor (~0.449).

B_to_A_rescan:
  lowers beta2 in no_backreaction and worsens the best J-lock.
  signed_birth remains large and negative in both real/no-back variants.

stale_same_scan:
  can raise beta2 to 6 in real_growth,
  but is not a reliable dynamic legality rule and weakens signed/J-lock balance.
```

## Current conclusion

The two-pair assembly motif is real, but the order/context ablation does not yet align:

```text
strong H2/QP carrier
+ strong signed orientation
+ low J-lock
```

in the same dynamic variant. The tension is now order-sensitive, not simply a missing-pair problem.
