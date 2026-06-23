# CNNA live semigroup operator-family # candidate gate

## Purpose

This package continues the established script-1/script-2 true-Schur/DtN growth line, but adds the user's depth-scaling point:

```text
Irreversibility is not merely a property of a single birth event.
It accumulates with growth depth because birth effects distribute into old conductances and then mix with live relaxation.
```

The package therefore tests two things together:

1. **Operator-family # candidate:** from bridge-/passivity-positive Record→Live rows, build a small real operator family on each fixed boundary cut using
   `A_k = G^+(Lambda_k - Lambda_record)` and the record-DtN metric adjoint `#_G`.
2. **Depth scaling:** group the record/live gap and operator-family diagnostics by `parent_level` over finite L2/L3/L4 approximants.

No `J`, `i`, Hodge, physical Hilbert adjoint, C*-norm, Q/P target or delta-beta decision is used.

## Final-level summary table

| variant | rows | candidate+ | passivity+ | weak # family pass | mean gap | product resid | commutator resid | commutator norm |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| real_growth_linear_star_candidate | 102 | 0.157 | 0.157 | 0.000 | 0.03798 | 0.000 | 0.078 | 0.000 |
| log_growth_star_candidate | 102 | 0.157 | 0.157 | 0.000 | 0.02272 | 0.000 | 0.081 | 0.000 |
| saturating_growth_star_candidate | 102 | 1.000 | 1.000 | 0.039 | 0.3439 | 0.018 | 0.490 | 0.002 |
| strict_symmetrized_response_star_control | 102 | 0.000 | 0.000 | 0.000 | 0 | 0.000 | 0.000 | 0.000 |
| birth_only_no_relax_star_control | 0 | 0.000 | 0.000 | 0.000 | 0 | 0.000 | 0.000 | 0.000 |

## Depth-slope audit

| variant | grouped depth rows | log-gap slope vs parent_level | deep/shallow gap ratio |
|---|---:|---:|---:|
| log_growth_star_candidate | 5 | -0.003 | 0.985 |
| real_growth_linear_star_candidate | 5 | 0.129 | 1.278 |
| saturating_growth_star_candidate | 5 | -0.002 | 0.987 |
| strict_symmetrized_response_star_control | 5 | 0.000 | 0.000 |

## Interpretation guide

- The `#_G` operation is the record-DtN metric adjoint.  Its isolated existence is linear-algebraic and not counted as a result.
- The weak `# family pass` requires simultaneous passivity/bridge positivity, # closure, involution/anti-multiplicativity diagnostics, and approximate product/commutator closure of the finite relaxation-generated operator family.
- The depth audit is the important scaling check: if record/live gaps increase with `parent_level` and finite L, then irreversibility is a growth-depth effect rather than a constant single-birth feature.

## Main result

The package intentionally separates two claims:

```text
A. Live semigroup gives a passivity-/adjunction-like operator precursor.
B. That precursor already forms a stable real *-algebra-like family.
```

A is supported in response variants, especially saturating growth.  B is only partially supported: # closure is mostly tautological, while product/commutator closure remains nontrivial and not uniformly strong.  The strict-sym and birth-only controls stay null.

## Next test

`test_depth_scaling_irreversibility_extrapolation_gate.py`

Run larger finite approximants or optimized aggregated updates to decide whether the record/live gap per generation monotonically rises, saturates, or decays in the large-depth regime.  This should be done before interpreting the operator-family # candidate as a mature algebraic structure.

## Supplement: true-Schur irreversibility depth scaling from package 84

The operator-family # audit was run for L2/L3 because the full per-cut family closure at L4 is expensive.  To keep the user's scaling constraint explicit, this package also includes the already computed true-Schur birth/relaxation L2/L3/L4 irreversibility rows from the established package 84 line as:

```text
outputs/supplement_irreversibility_depth_scaling_from_pkg84.csv
outputs/supplement_irreversibility_depth_slope_from_pkg84.csv
```

For `real_growth_linear_true_schur_birth_relax_irreversible_mixing`, the L4 generation means are:

| parent level | events | mean record/live gap | mean total relax drift | mean birth parent DtN delta |
|---:|---:|---:|---:|---:|
| 0 | 3 | 0.015216 | 0.017788 | 0.244239 |
| 1 | 9 | 0.023506 | 0.026053 | 0.266538 |
| 2 | 27 | 0.028292 | 0.030892 | 0.282122 |
| 3 | 81 | 0.032707 | 0.035452 | 0.300541 |

The fitted log-gap slopes for the same variant are:

| max L | depth groups | log-gap slope vs parent level | deep/shallow gap ratio |
|---:|---:|---:|---:|
| 2 | 2 | 0.435 | 1.545 |
| 3 | 3 | 0.310 | 1.859 |
| 4 | 4 | 0.248 | 2.150 |

This supports the user's point in the finite data: the direct record/live irreversibility gap is small near the root and grows with parent depth.  It is not yet a large-depth theorem; it is a finite L2-L4 scaling signal.

## Revised interpretation

The # candidate line and the depth-scaling line should be read separately:

```text
operator-family # precursor:
  present as a weak passivity/metric-adjunction precursor,
  not yet a mature *-algebra.

irreversibility scaling:
  record/live gap grows with depth in real_growth L2-L4,
  supporting the claim that irreversibility is accumulated by growth and live relaxation mixing.
```

The next test should therefore prioritize depth scaling before stronger algebraic claims.
