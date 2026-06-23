# RESULTS — pairing quadrature adjoint-pairing gate

## Comparative table

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

## Interpretation protocol

This package separates four claims:

```text
1. Pair conjugation exists:       C^2 = +I.
2. Oriented cochain exchange:     J^2 = -I locally on each glued face-pair.
3. Conjugation compatibility:     C J C = -J.
4. Actual Q/P lock:               J(Q) = P and J(P) = -Q.
```

Only (1)-(3) are structural local pair-algebra tests.  Claim (4) is the nontrivial
quadrature-lock gate.  A nonzero commutator magnitude is not enough; signed/kappa-like
behavior is logged separately.

## Anti-smuggling conditions

- no `i`, no imported `J`, no Hodge star, no physical adjoint, no positivity;
- no final `sym(M)` in the directed birth-transport operator;
- `J_pair` is the oriented cochain-exchange map induced by the actual paired faces;
- local `J_pair^2=-I` is not interpreted as a global complex structure unless Q/P lock
  and coherence gates also pass;
- `decision_used_delta_beta_any` must remain false.
