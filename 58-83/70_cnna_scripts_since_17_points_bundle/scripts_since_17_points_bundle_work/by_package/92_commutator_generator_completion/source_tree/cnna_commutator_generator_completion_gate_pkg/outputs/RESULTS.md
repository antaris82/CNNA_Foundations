# CNNA commutator-generator completion gate

## Model status

This package continues the established script-1/script-2 true-Schur/DtN Record→Live growth line.

Implemented model class: **growing real complement network with true Schur/DtN live relaxation**.

It is not SG-like geometry, not a naive static NGF attachment, and not a CQNM/s=-1 saturated geometry claim.  The geometry here is the existing deterministic ternary growth model with sequential births, incoming environment response, newborn backreaction, fixed-topology live relaxation, and true Schur complement/DtN boundary matrices.

No `J`, no `i`, no Hodge star, no physical Hilbert adjoint, no C*-norm, no Q/P target, no positivity axiom, and no `delta_beta` decision is used.

## Purpose

The previous coarse-graining package showed that generation/principal-mode coarse-graining improves the small real `#/*`-like operator-family closure, but the commutator residual remained the main obstruction.

This test starts from the same coarse-grained Record→Live true-Schur/DtN operator families and adds **only derived commutator modes**:

```text
[A_i,A_j] = A_i A_j - A_j A_i
```

No arbitrary generator is inserted.  The test compares:

```text
base
commutator_completed_1
commutator_completed_2
```

The strongest source family remains `mean_plus_principal`.

## Runtime note

The package runs L3 fully and L4 as a deterministic sampled-frontier approximant after `full_until_level=3`.  A full true-Schur/DtN L4 rerun was attempted but was too expensive in this runtime.  Therefore the L4 numbers should be read as a deterministic sampled scaling diagnostic, not as a full-tree theorem.

## Deepest sampled level: L4

### real_growth_linear

```text
mean_plus_principal / base:
  strong pass fraction        ≈ 0.233
  product residual            ≈ 0.0092
  commutator residual         ≈ 0.1955
  commutator max residual     ≈ 0.2725
  rank                        ≈ 5.43

mean_plus_principal / commutator_completed_1:
  strong pass fraction        ≈ 0.567
  product residual            ≈ 0.00225
  commutator residual         ≈ 0.000089
  commutator max residual     ≈ 0.000719
  commutator modes added      ≈ 3.27
  rank                        ≈ 8.47

mean_plus_principal / commutator_completed_2:
  strong pass fraction        ≈ 0.567
  product residual            ≈ 0.00060
  commutator residual         ≈ 0.000020
  commutator max residual     ≈ 0.000159
  rank                        ≈ 10.13
```

### saturating_growth

```text
mean_plus_principal / base:
  strong pass fraction        ≈ 0.119
  product residual            ≈ 0.0339
  commutator residual         ≈ 0.2907
  commutator max residual     ≈ 0.4328
  rank                        ≈ 6.69

mean_plus_principal / commutator_completed_1:
  strong pass fraction        ≈ 0.322
  product residual            ≈ 0.0199
  commutator residual         ≈ 0.000418
  commutator max residual     ≈ 0.00251
  commutator modes added      ≈ 5.39
  rank                        ≈ 11.93

mean_plus_principal / commutator_completed_2:
  strong pass fraction        ≈ 0.322
  product residual            ≈ 0.00820
  commutator residual         ≈ 0.000222
  commutator max residual     ≈ 0.00132
  rank                        ≈ 16.34
```

### strict_sym control

The strict-symmetrized control produces no nontrivial commutator-completion rows.  Nontrivial closure does not arise from topology alone in this test.

## Interpretation

This is a strong partial positive result.

The earlier obstruction was specifically that product closure improved under coarse-graining while commutator closure remained too weak.  Here, adding only derived commutator modes almost eliminates the commutator residual in the best coarse family:

```text
real_growth mean_plus_principal:
  commutator residual 0.1955 → 0.000089 → 0.000020

saturating mean_plus_principal:
  commutator residual 0.2907 → 0.000418 → 0.000222
```

The product residual also improves rather than collapses:

```text
real_growth mean_plus_principal:
  product residual 0.0092 → 0.00225 → 0.00060

saturating mean_plus_principal:
  product residual 0.0339 → 0.0199 → 0.00820
```

So the missing commutator closure was at least partly a **missing derived-generator** issue, not a hard obstruction of the Record→Live operator family.

## Conservative status

This still does **not** prove a physical `*`-algebra or C*-algebra.

What is now supported:

```text
Record→Live true-Schur/DtN operators exist.
Coarse-graining stabilizes product closure.
Derived commutator completion stabilizes commutator closure.
The G-metric adjoint # is compatible at the diagnostic level.
strict_sym remains null.
```

What is not shown:

```text
No J.
No i.
No physical Hilbert adjoint.
No positivity theorem.
No C*-norm.
No GNS/Tomita construction.
No full L4/L5/L6 theorem.
```

## Methodological conclusion

The previous phrase "stable small real *-algebra-like operator family: no" must now be updated:

```text
raw family: no
coarse-grained family: partial
coarse-grained + derived commutator completion: strong candidate, but sampled and pre-C* only
```

## Next test

The next test should not jump to complex structure.  It should audit the new completed family:

```text
test_completed_operator_family_depth_stability_gate.py
```

Goal:

```text
1. Track whether the completed family is stable across depth/generation.
2. Compare L3 full vs L4 sampled vs higher sampled levels.
3. Measure whether the same generator roles recur or drift.
4. Check # closure, product closure, commutator closure, and strict_sym null.
5. Do not introduce J, i, positivity, or C*-norm.
```
