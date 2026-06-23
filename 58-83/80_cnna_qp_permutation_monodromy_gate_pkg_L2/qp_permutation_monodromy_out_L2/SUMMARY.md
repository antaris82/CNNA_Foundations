# CNNA Q/P permutation monodromy gate — SUMMARY

| option | variant | beta | pairs | asm | pair harm | Q harm | P harm | pair J-lock | signed | valid mono | mono pass | polar -I | polar +I | polar double | raw -I | raw +I | raw double | used dBeta? |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| A_to_B_rescan_strong_reuseB | real_growth | (1,0,4,0) | 4 | 1 | 0.381795 | 0.381795 | 0.248801 | 0.448553 | -0.624664 | 1 | 0 | 1 | 5.03074e-14 | 1.00615e-13 | 0.966246 | 0.520857 | 0.726345 | False |
| A_to_B_rescan_strong_reuseB | strict_symmetrized_control | (1,0,0,0) | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | False |
| A_to_B_rescan_strong_reuseB | no_backreaction | (1,0,8,0) | 4 | 2 | 0.405623 | 0.405623 | 0.343367 | 0.339477 | -0.285265 | 2 | 0 | 0.866025 | 8.63562e-14 | 6.98922e-14 | 0.909288 | 0.476767 | 0.691778 | False |

## Decision

This package tests the spinor/double-cover suspicion properly: it does not square a single operator.  It builds the closed sibling-label path

```text
id -> rho -> rho^2 -> id, rho: 1 -> 2 -> 3 -> 1
```

and recomputes the same A/B assembly motif Q/P carrier at each label configuration.  The monodromy is the composed basis/subspace transport around this closed path.

A positive spinor-like monodromy would require one loop close near `-I` and two loops close near `+I`, with good step overlaps and strict_sym null.  The observed result has zero monodromy pass count in all variants.
