# CNNA test: double-history parent-fan plaquette J-holonomy

This package is the first bridge test that combines:

1. the `parent_fan_triangle` plaquette J-candidate from full DtN face commutators;
2. the double-history suffix quotient from the dynamic DtN gluing tests.

It deliberately avoids minimal SO(3) axis transport. The primary gluing law is provenance-directed:

```text
root sector 1 -> 2 -> 3 -> 1
same suffix after forgetting the first root-sector label
face labels shifted by the directed sector difference
```

The test asks whether a local Z3/parent-fan J sector remains nontrivial under this provenance identification.

It does not derive a global J-field, complex Hilbert space, C*-algebra, GNS representation, or AQFT net.

## Run

```bash
python3 test_double_history_parent_fan_plaquette_J_holonomy.py --max-level 6
```

Outputs are written to `double_history_parent_fan_J_holonomy_out_L6/`.

## Main controls

- `identical_history_control`: must vanish.
- `diagonal`, `trace_scalar`, `common_mean_diagonal`: must kill the noncommutative DtN carrier.
- `symmetrized_birth`: should suppress K/J amplitude.
- `random_same_level_cycle_baseline`: checks whether the quotient signal is structured.

## Reading rule

A genuine Stufe-4 sign/holonomy candidate would need more than local J-validity. It should produce a double-history-specific sign/holonomy obstruction that vanishes for identical history, collapses under symmetrized birth and commuting reductions, and is not reproduced by random same-level cycles.
