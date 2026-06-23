# SUMMARY — dual pairing two-edge assembly gate

This package follows the tradeoff obstruction gate.  The previous result showed that beta2-opening, Q/P support, C-eigen lock and kappa-signed flip exist in the candidate space but do not coincide on one candidate.

This package checks whether they can be split across two coupled pairings in the same scan.

| variant | candidates | role A beta/QP | role B C/kappa | single-pair all | two-edge same-scan | connected | strong | used Δβ? |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| real_growth | 206 | 80 | 3 | 0 | 38 | 38 | 26 | False |
| strict_symmetrized_control | 0 | 0 | 0 | 0 | 0 | 0 | 0 | False |
| no_backreaction | 236 | 103 | 2 | 0 | 22 | 22 | 16 | False |

Conservative reading: a nonzero connected two-edge assembly count is not a derivation of `i`, `J`, `*`, positivity or a C*-structure.  It only means the single-pair-local obstruction was too strict and the next object is a multi-edge assembly.
