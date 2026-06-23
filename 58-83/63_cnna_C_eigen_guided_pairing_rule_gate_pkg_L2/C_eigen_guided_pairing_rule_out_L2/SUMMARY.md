# SUMMARY — C-eigen guided pairing rule gate

Model label:
CNNA growing primal simplicial complex with deterministic sequential provenance growth, nonlinear asymmetry-gated complement-pair candidate space, directed antisymmetric birth-transport operators, and local C/J pair algebra.

This package compares two selection rules:

```text
A_rank_rule:
  inherited nonlinear cascade ranking.

C_eigen_guided_rule:
  select only from legal A-gated candidates,
  but rank by native C-eigen J-lock residual.
```

The guided rule does not use delta_beta/H2 as a decision input.  Delta-beta columns are audit-only after scan enumeration.

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
