# RESULTS — C-eigen guided pairing rule gate

## Comparative table

| variant/rule/phase | beta | pairs | selected C-lock | pair harm | Q harm | P harm | best J-lock | signed birth | used dBeta? |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| real_growth_A_rank_rule_phaseplus1 | (1,0,2,0) | 2 | 0.375414 | 0.229731 | 0.229731 | 0.222362 | 0.375414 | 0.0305323 | False |
| real_growth_C_eigen_guided_rule_phaseplus1 | (1,0,3,0) | 2 | 0.251636 | 0.286611 | 0.286611 | 0.169823 | 0.251636 | -0.37836 | False |
| real_growth_A_rank_rule_phaseminus1 | (1,0,2,0) | 2 | 0.229976 | 0.20605 | 0.20605 | 0.207661 | 0.229976 | 0.151928 | False |
| real_growth_C_eigen_guided_rule_phaseminus1 | (1,0,3,0) | 2 | 0.41447 | 0.259184 | 0.259184 | 0.235788 | 0.41447 | 0.073591 | False |
| strict_symmetrized_control_A_rank_rule_phaseplus1 | (1,0,0,0) | 0 | 0 | 0 | 0 | 0 | 0 | 0 | False |
| strict_symmetrized_control_C_eigen_guided_rule_phaseplus1 | (1,0,0,0) | 0 | 0 | 0 | 0 | 0 | 0 | 0 | False |
| strict_symmetrized_control_A_rank_rule_phaseminus1 | (1,0,0,0) | 0 | 0 | 0 | 0 | 0 | 0 | 0 | False |
| strict_symmetrized_control_C_eigen_guided_rule_phaseminus1 | (1,0,0,0) | 0 | 0 | 0 | 0 | 0 | 0 | 0 | False |
| no_backreaction_A_rank_rule_phaseplus1 | (1,0,2,0) | 2 | 0.203483 | 0.223795 | 0.223795 | 0.220329 | 0.203483 | 0.0362835 | False |
| no_backreaction_C_eigen_guided_rule_phaseplus1 | (1,0,3,0) | 2 | 0.314669 | 0.278848 | 0.278848 | 0.168394 | 0.314669 | -0.447414 | False |
| no_backreaction_A_rank_rule_phaseminus1 | (1,0,2,0) | 2 | 0.190115 | 0.203952 | 0.203952 | 0.205505 | 0.190115 | 0.115044 | False |
| no_backreaction_C_eigen_guided_rule_phaseminus1 | (1,0,3,0) | 2 | 0.200924 | 0.266294 | 0.266294 | 0.230453 | 0.200924 | 0.0729985 | False |

## Phase flip comparison

| variant/rule | signed + | signed - | flip score | pair harm + | pair harm - | J-lock + | J-lock - |
|---|---:|---:|---:|---:|---:|---:|---:|
| no_backreaction / A_rank_rule | 0.0362835 | 0.115044 | 1 | 0.223795 | 0.203952 | 0.203483 | 0.190115 |
| no_backreaction / C_eigen_guided_rule | -0.447414 | 0.0729985 | -0.719459 | 0.278848 | 0.266294 | 0.314669 | 0.200924 |
| real_growth / A_rank_rule | 0.0305323 | 0.151928 | 1 | 0.229731 | 0.20605 | 0.375414 | 0.229976 |
| real_growth / C_eigen_guided_rule | -0.37836 | 0.073591 | -0.674341 | 0.286611 | 0.259184 | 0.251636 | 0.41447 |
| strict_symmetrized_control / A_rank_rule | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| strict_symmetrized_control / C_eigen_guided_rule | 0 | 0 | 0 | 0 | 0 | 0 | 0 |

## Interpretation protocol

A constructive success would require:

```text
1. beta2 still opens,
2. pair_transport_harmonic_ratio remains positive,
3. Q and P harmonic channels remain positive,
4. C-eigen J-lock residual improves relative to A_rank_rule,
5. strict_sym stays killed,
6. decision_used_delta_beta_any remains false.
```

A strong J/i claim is still forbidden.  This test only checks whether the previous obstruction came from selecting the wrong pairings.
