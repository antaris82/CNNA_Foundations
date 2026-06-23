# SUMMARY — pair J-alignment search gate

Model label:
CNNA growing primal simplicial complex with deterministic sequential provenance growth,
nonlinear asymmetry-gated complement pairing, directed antisymmetric birth-transport
operators, and local pair-exchange algebra.

This package does **not** fit an arbitrary rotation and does **not** define P as JQ.  It tests
only already-derived candidate pairs:

```text
raw Q/P, unit Q/P, norm-matched Q/P,
harmonic Q/P, harmonic unit/norm-matched Q/P,
Q vs commutator, P vs commutator,
C-pair eigen even/odd projections of raw pair data.
```

Gate:
`J_pair(Q') ≈ P'` and `J_pair(P') ≈ -Q'` must hold for an allowed candidate, while
strict_sym remains killed.

| variant/phase | beta | pairs | Q harm | P harm | best candidate | best mean resid | best max resid | per-pair best mean | comm signed birth | used dBeta? |
|---|---:|---:|---:|---:|---|---:|---:|---:|---:|---:|
| real_growth_phaseplus1 | (1,0,2,0) | 2 | 0.229731 | 0.222362 | unit_C_eigen_transported_even_odd | 0.375414 | 0.375414 | 0.375414 | 0.0305323 | False |
| real_growth_phaseminus1 | (1,0,2,0) | 2 | 0.20605 | 0.207661 | unit_C_eigen_transported_even_odd | 0.229976 | 0.229976 | 0.229976 | 0.151928 | False |
| strict_symmetrized_control_phaseplus1 | (1,0,0,0) | 0 | 0 | 0 |  | 0 | 0 | 0 | 0 | False |
| strict_symmetrized_control_phaseminus1 | (1,0,0,0) | 0 | 0 | 0 |  | 0 | 0 | 0 | 0 | False |
| no_backreaction_phaseplus1 | (1,0,2,0) | 2 | 0.223795 | 0.220329 | unit_C_eigen_transported_even_odd | 0.203483 | 0.203483 | 0.203483 | 0.0362835 | False |
| no_backreaction_phaseminus1 | (1,0,2,0) | 2 | 0.203952 | 0.205505 | unit_C_eigen_raw_even_odd | 0.190115 | 0.190115 | 0.190115 | 0.115044 | False |

## Phase-sign comparison

| variant | best + | best - | mean + | mean - | pair mean + | pair mean - | comm + | comm - | flip-score |
|---|---|---|---:|---:|---:|---:|---:|---:|---:|
| no_backreaction | unit_C_eigen_transported_even_odd | unit_C_eigen_raw_even_odd | 0.203483 | 0.190115 | 0.203483 | 0.190115 | 0.0362835 | 0.115044 | 1 |
| real_growth | unit_C_eigen_transported_even_odd | unit_C_eigen_transported_even_odd | 0.375414 | 0.229976 | 0.375414 | 0.229976 | 0.0305323 | 0.151928 | 1 |
| strict_symmetrized_control |  |  | 0 | 0 | 0 | 0 | 0 | 0 | 0 |

Conservative reading: if all allowed candidates keep large lock residuals, the local pair
algebra exists but does not dynamically select the Q/P or a/a† split.
