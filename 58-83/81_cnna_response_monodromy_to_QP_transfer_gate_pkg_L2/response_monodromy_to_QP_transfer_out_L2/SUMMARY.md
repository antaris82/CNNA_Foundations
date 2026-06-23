# SUMMARY — response monodromy to Q/P transfer gate

This package is deliberately result-open.  It does not ask only for a Spin-1/2/double-cover signature.  It asks whether the lower sibling-response monodromy/circulation layer can be transferred to the later Q/P assembly motif layer in any clear form:

- Q/P-lock transfer,
- complex-eigen/rotational transfer,
- C3-like monodromy transfer,
- double-cover-like sign transfer,
- identity/contraction/leakage diagnostics.

| variant | beta | assemblies | pair harm | Q harm | P harm | base J-lock | response transfer active | QP-lock pass | complex-like pass | C3-like pass | double-cover pass | best QP transfer | max imag | best C3 resid | used Δβ? |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| real_growth | (1,0,4,0) | 1 | 0.381795 | 0.381795 | 0.248801 | 0.448553 | 35 | 0 | 16 | 0 | 0 | 0.866456 | 0.486059 | 0.70974 | False |
| strict_symmetrized_control | (1,0,0,0) | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | False |
| no_backreaction | (1,0,8,0) | 2 | 0.405623 | 0.405623 | 0.343367 | 0.339477 | 70 | 0 | 28 | 0 | 0 | 0.816354 | 4.77306e+06 | 0.79597 | False |

The test uses actual response data from the directed sibling fans (`directed_edges`) and bridges it to the Q/P motif carrier through face/node incidence and birth-frame projections.  It does not use matrix powers of a prebuilt Q/P operator as a substitute for monodromy.
