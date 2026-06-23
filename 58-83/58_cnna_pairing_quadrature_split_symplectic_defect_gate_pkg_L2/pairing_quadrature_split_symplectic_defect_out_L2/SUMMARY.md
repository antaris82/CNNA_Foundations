# SUMMARY — pairing quadrature split symplectic defect gate

Model label:
CNNA growing primal simplicial complex with deterministic sequential provenance growth,
nonlinear asymmetry-gated complement pairing, and a directed antisymmetric birth-transport
operator.  This is not SG/ST as a global geometry, not a finished NGF/CQNM model,
and not a complex/J/*/positivity derivation.

Purpose:

```text
Test whether the apparent cancellation in the paired channel is only the even
real observable-like projection Q = ka + transport(kb), while the odd channel
P = ka - transport(kb) carries a second real quadrature.
```

No `i`, no `J`, no Hodge star, no positive-frequency split, no imported adjoint,
no positivity, and no final sym(M) is used.  Norms are finite-dimensional diagnostics
only, not axioms.

| variant | beta | pairs | Q harm | P harm | P/Q | P_H/Q_H | sympl area | Q kappa | P kappa | used dBeta? |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| real_growth | (1,0,2,0) | 2 | 0.229731 | 0.222362 | 1.13589 | 1.09946 | 0.771655 | 0.0318445 | 0.0639773 | False |
| strict_symmetrized_control | (1,0,0,0) | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | False |
| no_backreaction | (1,0,2,0) | 2 | 0.223795 | 0.220329 | 1.02504 | 1.00916 | 0.639563 | 0.0293466 | 0.0673154 | False |

Conservative reading:
If P is nonzero and harmonic while strict_sym remains zero, the pair cancellation is
not just a failure; it exposes a candidate real two-quadrature split.  This does not
derive a complex structure.  It only moves the question from spatial orientation-lock
to a possible real operator/quadrature split.
