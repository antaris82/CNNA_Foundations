# RESULTS — dual-pairing assembly growth rule gate

## Comparative table

| variant | beta | pairs | assemblies | pair harm | Q harm | P harm | best J-lock | signed birth | selected C-lock | selected kappa-flip | used dBeta? |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| real_growth | (1,0,4,0) | 4 | 1 | 0.381795 | 0.381795 | 0.248801 | 0.448553 | -0.624664 | 0.491433 | 1 | False |
| strict_symmetrized_control | (1,0,0,0) | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | False |
| no_backreaction | (1,0,6,0) | 4 | 2 | 0.302564 | 0.302564 | 0.304787 | 0.208886 | 0.0704225 | 0.331361 | 0.80838 | False |

## Interpretation protocol

A constructive success would require:

```text
1. strict_sym remains zero/killed;
2. real_growth opens beta2;
3. Q/P and pair-transport harmonic channels remain positive;
4. J-lock improves or remains competitive;
5. signed orientation improves without phase/i/Hodge/* input;
6. decision_used_delta_beta_any remains false.
```

The rule is deliberately two-step: Pair B is chosen after applying Pair A and rescanning the updated complex.  This avoids pretending that a stale same-scan shared-face candidate is automatically still legal after Pair A.
