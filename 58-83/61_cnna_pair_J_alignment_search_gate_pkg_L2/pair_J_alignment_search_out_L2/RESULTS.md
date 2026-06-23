# RESULTS — pair J-alignment search gate

## Comparative table

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

## Interpretation protocol

A positive result would require an allowed, derived candidate channel with small residuals:

```text
J_pair(Q') -> P'
J_pair(P') -> -Q'
```

This package deliberately rejects the circular construction `P' := J_pair(Q')`.  The question
is whether Q/P-like channels already present in the data align with the local pair J.

## Anti-smuggling conditions

- no `i`, no complex scalars;
- no imported Hodge star, positivity, physical adjoint, or norm axiom;
- no final `sym(M)` in the directed birth-transport operator;
- no arbitrary fitted rotation;
- no topology/H2/kappa used in move decisions;
- `decision_used_delta_beta_any` must remain false.
