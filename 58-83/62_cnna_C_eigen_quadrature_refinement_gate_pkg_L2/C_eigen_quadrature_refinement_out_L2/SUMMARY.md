# SUMMARY — C-eigen quadrature refinement gate

This package audits whether the current nonlinear asymmetry-gated pairing rule misses beta2-opening candidates that have a better native C-pair even/odd J-lock.

It evaluates legal pair candidates at each scan state before the current selection rule applies.  The gate does not fit a new rotation and does not introduce i, Hodge, positivity, a physical adjoint, or a final symmetrization.

| variant phase | beta auto | pairings | ok candidates | beta2 candidates | beta2 + C-lock | best beta2 lock | selected lock | best beta2 selected? | best beta2 cos |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| real_growth_phaseplus1 | (1,0,2,0) | 2 | 13 | 9 | 1 | 0.138006 / 0.138006 | 0.288786 | False | -0.0602768 |
| real_growth_phaseminus1 | (1,0,2,0) | 2 | 13 | 9 | 1 | 0.173092 / 0.173092 | 0.212941 | False | 0.192677 |
| strict_symmetrized_control_phaseplus1 | (1,0,0,0) | 0 | 0 | 0 | 0 | 0 / 0 | 0 | False | 0 |
| strict_symmetrized_control_phaseminus1 | (1,0,0,0) | 0 | 0 | 0 | 0 | 0 / 0 | 0 | False | 0 |
| no_backreaction_phaseplus1 | (1,0,2,0) | 2 | 13 | 9 | 2 | 0.103391 / 0.103391 | 0.168755 | False | -0.0697727 |
| no_backreaction_phaseminus1 | (1,0,2,0) | 2 | 13 | 9 | 4 | 0.154455 / 0.154455 | 0.181179 | False | 0.306041 |

Main question:

```text
Do there exist candidates with delta_beta2 > 0 and C-eigen J-lock residual below the threshold?
```
