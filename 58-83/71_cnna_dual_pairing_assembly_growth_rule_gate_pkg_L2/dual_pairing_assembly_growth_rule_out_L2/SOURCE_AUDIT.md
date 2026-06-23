# SOURCE AUDIT

Carried forward:

- Pair-property tradeoff gate showed single-pair locality is obstructed: beta2/QP, C-lock, and kappa-flip split across different candidates.
- Dual-pair assembly audit showed same-scan connected two-pair motifs exist.

This package turns that audit into a dynamic growth-rule test.

Anti-smuggling constraints:

- Pair A is selected using provenance/QP proxies, not delta_beta.
- Pair B is selected using C-lock/kappa/context proxies after a rescan, not delta_beta.
- delta_beta, H2 and harmonic diagnostics are audit-only after application.
- no i, no global J, no Hodge star, no positivity, no C*-norm, no final sym(M).

Limitations:

- This is L2 only.
- The kappa diagnostic is label-kappa on a fixed grown model, not a fully regrown reverse-sibling universe.
- If Pair B cannot be found after Pair A, that is an actual dynamic-legality obstruction, not a code failure.
