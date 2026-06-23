# RESULTS — signed quadrature area kappa gate

## Comparative table

| variant/phase | beta | pairs | Q harm | P harm | abs area | signed birth/abs | signed outward/abs | Q kappa | P kappa | used dBeta? |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| real_growth_phaseplus1 | (1,0,2,0) | 2 | 0.229731 | 0.222362 | 0.771655 | 0.0305323 | 0.0305323 | 0.0318445 | 0.0639773 | False |
| real_growth_phaseminus1 | (1,0,2,0) | 2 | 0.20605 | 0.207661 | 0.495372 | 0.151928 | 0.151928 | 0.0101165 | 0.0215065 | False |
| strict_symmetrized_control_phaseplus1 | (1,0,0,0) | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | False |
| strict_symmetrized_control_phaseminus1 | (1,0,0,0) | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | False |
| no_backreaction_phaseplus1 | (1,0,2,0) | 2 | 0.223795 | 0.220329 | 0.639563 | 0.0362835 | 0.0362835 | 0.0293466 | 0.0673154 | False |
| no_backreaction_phaseminus1 | (1,0,2,0) | 2 | 0.203952 | 0.205505 | 0.373249 | 0.115044 | 0.115044 | 0.00808436 | 0.0223135 | False |

## Phase-sign flip comparison

| variant | birth plus | birth minus | birth flip-score | outward plus | outward minus | outward flip-score |
|---|---:|---:|---:|---:|---:|---:|
| no_backreaction | 0.0362835 | 0.115044 | 1 | 0.0362835 | 0.115044 | 1 |
| real_growth | 0.0305323 | 0.151928 | 1 | 0.0305323 | 0.151928 | 1 |
| strict_symmetrized_control | 0 | 0 | 0 | 0 | 0 | 0 |

## Interpretation protocol

This package explicitly distinguishes:

```text
abs area     = ||Q x P||
signed area  = <Q x P, n_ref>
```

The abs area is not allowed to count as oriented symplectic structure.  The gate is
only positive if the signed ratios are nontrivial and transform coherently under the
phase/kappa proxy while strict_sym remains killed.

## Anti-smuggling conditions

- no `i`, no `J`, no Hodge star, no imported adjoint, no positivity;
- no final `sym(M)` in the antisymmetric birth-transport operator;
- `decision_used_delta_beta_any` must remain false;
- abs-area is treated as magnitude only.
