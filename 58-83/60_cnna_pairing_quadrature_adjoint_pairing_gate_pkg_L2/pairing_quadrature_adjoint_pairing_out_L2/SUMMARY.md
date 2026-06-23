# SUMMARY — pairing quadrature adjoint-pairing gate

Model label:
CNNA growing primal simplicial complex with deterministic sequential provenance growth,
nonlinear asymmetry-gated complement pairing, and directed antisymmetric birth-transport
operators.  This package tests a real pair-conjugation/quadrature structure.  It does not
claim complex scalars, a physical adjoint, positivity, norm, or C*-structure.

Core real pair operators for each actual glued face-pair:

```text
C_pair = [[0, R], [ R^T, 0]]       C_pair^2 = +I
J_pair = [[0, R], [-R^T, 0]]       J_pair^2 = -I
C_pair J_pair C_pair = -J_pair
```

Here R is the orientation-reversing transport determined by the actual face-pair gluing.
This is a local cochain-pair algebra, not a global i.

| variant/phase | beta | pairs | Q harm | P harm | J^2+I resid | C^2-I resid | CJC+J resid | JQ->P resid | comm signed birth | block comm grade | used dBeta? |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| real_growth_phaseplus1 | (1,0,2,0) | 2 | 0.229731 | 0.222362 | 4.2e-13 | 4.2e-13 | 4.2e-13 | 0.881133 | 0.0305323 | 1 | False |
| real_growth_phaseminus1 | (1,0,2,0) | 2 | 0.20605 | 0.207661 | 4.2e-13 | 4.2e-13 | 4.2e-13 | 0.966969 | 0.151928 | 1 | False |
| strict_symmetrized_control_phaseplus1 | (1,0,0,0) | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | False |
| strict_symmetrized_control_phaseminus1 | (1,0,0,0) | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | False |
| no_backreaction_phaseplus1 | (1,0,2,0) | 2 | 0.223795 | 0.220329 | 4.2e-13 | 4.2e-13 | 4.2e-13 | 0.932541 | 0.0362835 | 1 | False |
| no_backreaction_phaseminus1 | (1,0,2,0) | 2 | 0.203952 | 0.205505 | 4.2e-13 | 4.2e-13 | 4.2e-13 | 0.981961 | 0.115044 | 1 | False |

## Phase-sign flip comparison

| variant | comm birth + | comm birth - | comm flip-score | block grade + | block grade - | grade flip-score |
|---|---:|---:|---:|---:|---:|---:|
| no_backreaction | 0.0362835 | 0.115044 | 1 | 1 | 1 | 1 |
| real_growth | 0.0305323 | 0.151928 | 1 | 1 | 1 | 1 |
| strict_symmetrized_control | 0 | 0 | 0 | 0 | 0 | 0 |

Conservative reading:
If the local pair algebra is exact but `JQ->P` remains large and signed commutators do not
flip, then the test found a local real conjugation/J scaffold but not a derived Q/P
complex lock.
