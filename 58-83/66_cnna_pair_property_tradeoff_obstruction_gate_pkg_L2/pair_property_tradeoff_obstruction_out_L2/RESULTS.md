# RESULTS — pair property tradeoff obstruction gate

## Comparative table

| variant | candidates | beta2 | C-lock | kappa-flip | Q/P | beta+C+flip | all four | corr beta↔flip | corr C↔flip | corr QP↔flip |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| real_growth | 206 | 88 | 19 | 23 | 195 | 0 | 0 | -0.0773 | 0.0154 | -0.0843 |
| no_backreaction | 236 | 108 | 15 | 28 | 225 | 0 | 0 | -0.12 | -0.00329 | -0.132 |
| strict_symmetrized_control | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |

## Interpretation

The test quantifies the pair-candidate tradeoffs:

```text
beta2-opening vs C-lock
beta2-opening vs κ-flip
C-lock vs κ-flip
Q/P-support vs κ-flip
directed_imbalance vs κ-flip
transverse_complementarity vs C-lock
```

The hard all-four gate remains empty in both nontrivial variants.  This is not a mere ranking failure: good C-lock candidates and good signed-flip candidates occur in different parts of the candidate space.

The most important qualitative split is:

```text
C-lock pass candidates: signed_flip_abs is typically near 1.
κ-flip pass candidates: C_lock_max_worst is typically much larger.
beta2-opening candidates: usually have poor signed flip.
Q/P support is common and therefore not the limiting constraint.
```

Thus the obstruction is a compatibility split among properties, not absence of all ingredients.

## Conservative status

This package does not prove a fundamental no-go.  It only shows that, in the current L2 A-gated matched candidate space, the required properties do not co-localize on one candidate.

## Next test

`test_dual_pairing_two_edge_assembly_gate.py`

Rationale: if no single pair can carry all roles, test whether CNNA naturally requires a two-edge assembly: one pair opens beta2 / QP carrier, another supplies κ-flip / C-lock, and the compatibility must be checked at the assembled cycle level rather than on a single pair.
