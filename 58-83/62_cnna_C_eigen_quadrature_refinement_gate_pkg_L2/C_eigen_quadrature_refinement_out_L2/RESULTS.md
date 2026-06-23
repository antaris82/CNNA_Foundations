# RESULTS — C-eigen quadrature refinement gate

## Comparative table

| variant phase | beta auto | pairings | ok candidates | beta2 candidates | beta2 + C-lock | best beta2 lock | selected lock | best beta2 selected? | best beta2 cos |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| real_growth_phaseplus1 | (1,0,2,0) | 2 | 13 | 9 | 1 | 0.138006 / 0.138006 | 0.288786 | False | -0.0602768 |
| real_growth_phaseminus1 | (1,0,2,0) | 2 | 13 | 9 | 1 | 0.173092 / 0.173092 | 0.212941 | False | 0.192677 |
| strict_symmetrized_control_phaseplus1 | (1,0,0,0) | 0 | 0 | 0 | 0 | 0 / 0 | 0 | False | 0 |
| strict_symmetrized_control_phaseminus1 | (1,0,0,0) | 0 | 0 | 0 | 0 | 0 / 0 | 0 | False | 0 |
| no_backreaction_phaseplus1 | (1,0,2,0) | 2 | 13 | 9 | 2 | 0.103391 / 0.103391 | 0.168755 | False | -0.0697727 |
| no_backreaction_phaseminus1 | (1,0,2,0) | 2 | 13 | 9 | 4 | 0.154455 / 0.154455 | 0.181179 | False | 0.306041 |

## Conservative interpretation

The test separates two possibilities:

```text
1. The current pairing rule selects the wrong pairs.
2. beta2-opening and C-eigen J-lock are structurally in tension in the available L2 candidate space.
```

A positive hit would require at least one beta2-opening candidate with both mean and max C-eigen J-lock residual below the configured thresholds.
