# SUMMARY — pair property tradeoff obstruction gate

This package audits the matched identity/κ candidate space from the prior κ-permutation test.

It does not introduce a new pairing rule and does not use beta, H2, kappa, positivity, Hodge, `i`, `J`, or `*` as a move decision.

## Comparative result

| variant | candidates | beta2 | C-lock | kappa-flip | Q/P | beta+C+flip | all four | corr beta↔flip | corr C↔flip | corr QP↔flip |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| real_growth | 206 | 88 | 19 | 23 | 195 | 0 | 0 | -0.0773 | 0.0154 | -0.0843 |
| no_backreaction | 236 | 108 | 15 | 28 | 225 | 0 | 0 | -0.12 | -0.00329 | -0.132 |
| strict_symmetrized_control | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |

Main readout: the four properties exist separately, but the current L2 candidate space contains no candidate that combines beta2-opening, C-eigen J-lock, Q/P support, and signed κ-flip.
