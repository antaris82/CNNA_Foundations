# RESULTS — pairing quadrature split symplectic defect gate

## Comparative table

| variant | beta | pairs | Q harm | P harm | P/Q | P_H/Q_H | sympl area | Q kappa | P kappa | used dBeta? |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| real_growth | (1,0,2,0) | 2 | 0.229731 | 0.222362 | 1.13589 | 1.09946 | 0.771655 | 0.0318445 | 0.0639773 | False |
| strict_symmetrized_control | (1,0,0,0) | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | False |
| no_backreaction | (1,0,2,0) | 2 | 0.223795 | 0.220329 | 1.02504 | 1.00916 | 0.639563 | 0.0293466 | 0.0673154 | False |

## Interpretation protocol

```text
Q_even = ka + R(kb)
P_odd  = ka - R(kb)
```

where R transports face b into the orientation-reversed gluing chart of face a.

The decisive question is whether P survives when Q is antikohärent / cancellation-prone.
A positive P channel is not a proof of creation-annihilation structure.  It is only a
derived real pre-structure which could later be tested for operator composition,
anti-automorphism, CCR-like or symplectic closure.

## Anti-smuggling checks

- `decision_used_delta_beta_any` must remain false.
- `strict_symmetrized_control` must not produce the same Q/P signal.
- The antisymmetric transport term is derived from birth_order via q and h.
- No final symmetrization is applied to the directed vertex operator.
