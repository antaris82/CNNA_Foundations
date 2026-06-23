# RESULTS — Q/P permutation monodromy gate

## Main table

| option | variant | beta | pairs | asm | pair harm | Q harm | P harm | pair J-lock | signed | valid mono | mono pass | polar -I | polar +I | polar double | raw -I | raw +I | raw double | used dBeta? |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| A_to_B_rescan_strong_reuseB | real_growth | (1,0,4,0) | 4 | 1 | 0.381795 | 0.381795 | 0.248801 | 0.448553 | -0.624664 | 1 | 0 | 1 | 5.03074e-14 | 1.00615e-13 | 0.966246 | 0.520857 | 0.726345 | False |
| A_to_B_rescan_strong_reuseB | strict_symmetrized_control | (1,0,0,0) | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | False |
| A_to_B_rescan_strong_reuseB | no_backreaction | (1,0,8,0) | 4 | 2 | 0.405623 | 0.405623 | 0.343367 | 0.339477 | -0.285265 | 2 | 0 | 0.866025 | 8.63562e-14 | 6.98922e-14 | 0.909288 | 0.476767 | 0.691778 | False |

## Interpretation

The old script 2 already showed a positive response-layer result: sequential birth plus backreaction produces nonzero log-circulation and complex local directed Markov sectors, while symmetrized and path-only controls are real/degenerate.  It also showed that kappa flips the selected forward-cycle J but does not preserve birth order.

This package asks a different question on the later Q/P assembly layer: does the closed sibling-label cycle produce a holonomy/monodromy `-I` on the Q/P motif carrier?

Result: no.  In the tested L2 assembly path the closed rho-cycle is either near identity/weakly nontrivial or far from a clean `-I`; it never passes the double-cover gate.  Thus the previous alpha-power test was not a valid spinor diagnostic, and the correct monodromy test does not currently support a Spin-1/2/double-cover claim on the Q/P motif layer.
