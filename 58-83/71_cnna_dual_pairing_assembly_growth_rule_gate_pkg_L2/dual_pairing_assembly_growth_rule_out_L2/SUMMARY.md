# SUMMARY — dual-pairing assembly growth rule gate

Model label:
CNNA growing primal simplicial complex with deterministic sequential provenance growth, directed antisymmetric birth-transport operators, and a dynamic two-step pair assembly rule.

This package tests the motif found by the previous audit:

```text
Pair A: provenance/QP carrier proxy, selected without delta-beta.
Pair B: C-lock/kappa context proxy, selected after Pair A by a rescan.
```

No decision uses delta_beta, H2, complex scalars, Hodge, positivity, or a physical adjoint.  The vertex operator uses the antisymmetric birth-transport term and no final sym(M).

| variant | beta | pairs | assemblies | pair harm | Q harm | P harm | best J-lock | signed birth | selected C-lock | selected kappa-flip | used dBeta? |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| real_growth | (1,0,4,0) | 4 | 1 | 0.381795 | 0.381795 | 0.248801 | 0.448553 | -0.624664 | 0.491433 | 1 | False |
| strict_symmetrized_control | (1,0,0,0) | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | False |
| no_backreaction | (1,0,6,0) | 4 | 2 | 0.302564 | 0.302564 | 0.304787 | 0.208886 | 0.0704225 | 0.331361 | 0.80838 | False |
