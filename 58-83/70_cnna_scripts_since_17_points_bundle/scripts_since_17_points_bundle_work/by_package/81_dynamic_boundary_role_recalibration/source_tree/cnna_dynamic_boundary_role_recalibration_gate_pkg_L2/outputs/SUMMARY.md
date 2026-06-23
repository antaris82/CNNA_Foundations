# RESULTS — CNNA dynamic boundary role recalibration gate, L2

## Model provenance

This package reuses the established script-1/script-2 dynamic birth model:

- ternary sequential births;
- newborn conductance derived from parent-line + older-sibling environment;
- existing environment edges into the newborn;
- newborn backreaction to the parent line and older siblings;
- completed sibling-triple monodromy diagnostics.

The new diagnostic is event-resolved boundary-role recalibration.  It logs, per birth, how the same newborn has no own UV-tail from its own perspective while immediately becoming a UV-tail/backreaction source for the parent line.

## Gate table

| variant | events | triples | topological role | response role | retarded frac | advanced frac | conductance frac | mean sibling increment | mean |log circ| | full Markov complex frac |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| real_growth_linear_script1_2 | 12 | 4 | 1 | 1 | 1.000 | 0.000 | 1.000 | 0.375692 | 1.141519 | 1.000 |
| log_growth_script1_2 | 12 | 4 | 1 | 1 | 1.000 | 0.000 | 1.000 | 0.076662 | 0.981801 | 1.000 |
| saturating_growth_script1_2 | 12 | 4 | 1 | 1 | 1.000 | 0.000 | 1.000 | 0.086730 | 2.456621 | 1.000 |
| strict_symmetrized_response_control | 12 | 4 | 1 | 0 | 0.000 | 0.000 | 0.000 | 0.000000 | 0.000000 | 0.000 |
| reverse_birth_label_order_control | 12 | 4 | 1 | 1 | 1.000 | 0.000 | 1.000 | 0.375692 | 1.141519 | 1.000 |

## Interpretation

The gate separates two claims that should not be conflated:

1. **Topological boundary polarity**: every newborn has no own UV-tail at birth and simultaneously becomes a new UV-tail for the parent cut.  This remains true even in the strict-sym response control because it is part of the birth event itself.
2. **Conductance/response role recalibration**: parent-line and older-sibling conductances are updated directionally by the newborn.  This is positive in the real/log/saturating/reverse-label growth variants and collapses in the strict-sym response control.

The result therefore supports a pre-causal boundary-role layer, not a J/i/spin claim.  Growth creates an irreversible role update:

```text
child self-perspective: no own UV-tail
parent-line perspective: child is new UV-tail / backreaction source
older sibling perspective: already-born cells are updated by later births
ancestor/root perspective: descendant birth changes effective response
```

This is the CNNA analogue of a retarded boundary-value selection candidate.  It is not yet a Q/P, *, J, or complex structure.

## Important limitations

- The `cut_response_delta` / `Schur_delta` columns are scalar conductance-response surrogates, not full matrix Schur complements.
- The strict-sym control kills conductance/response asymmetry, but not the topological fact that a birth adds a child to the parent cut.
- The test uses L2 only; deeper levels should check whether old-interior role deltas decay while frontier updates dominate.
- No delta-beta, H², Q/P, J, Hodge, physical adjoint, or positivity criterion is used as a decision gate.

## Next test

`test_role_recalibration_to_boundary_value_bridge_gate.py`

Purpose: use the event-resolved role rows as the source layer and test whether the retarded boundary-role polarity transfers to the later Q/P/operator layer more robustly than the previous response-monodromy bridge.
