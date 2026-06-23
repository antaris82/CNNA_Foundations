# SUMMARY — bounded response bridge sanity gate

This package follows the response-monodromy-to-Q/P transfer package, but it does **not** treat every complex-like projected eigenvalue as evidence.  It first asks whether the response→Q/P bridge is numerically and structurally sane:

- bridge singular values and condition number are audited;
- ill-conditioned or pseudo-inverse-sensitive rows are separated;
- response operators are bounded by operator norm / spectral radius without fitting a target structure;
- Q/P-lock, complex-like, C3-like, and minus-identity-like signatures are only counted on well-conditioned bounded non-artifact rows;
- contraction/leakage/half-plane bias are logged as possible semigroup/boundary-value traces, not as complex structure.

| variant | beta | assemblies | pair_harm | Q_harm | P_harm | base J-lock | good bounded rows | structural pass | QP pass | complex pass | C3 pass | -I pass | best good QP resid | max good ImEig | best C3 resid | max halfplane bias | used Δβ? |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| real_growth | (1,0,4,0) | 1 | 0.278528 | 0.278528 | 0.348841 | 0.449108 | 175 | 0 | 0 | 0 | 0 | 0 | 0.865711 | 0.521013 | 0.728974 | 1 | False |
| strict_symmetrized_control | (1,0,0,0) | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | False |
| no_backreaction | (1,0,6,0) | 2 | 0.298469 | 0.298469 | 0.279644 | 0.337379 | 230 | 0 | 0 | 0 | 0 | 0 | 0.783058 | 0.484363 | 0.799401 | 1 | False |
