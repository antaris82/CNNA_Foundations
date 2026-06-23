# SOURCE AUDIT

Carried forward:

- Single-pair property tradeoff gate: beta2/QP, C-lock, and kappa-flip split across candidates.
- Dual-pair assembly audit: same-scan two-pair motifs exist.
- Dual-pair dynamic growth gate: dynamic A->B assemblies open larger beta2 and Q/P channels but signed orientation and J-lock do not stabilize together.

This package does not introduce a new final rule. It ablates the dynamic assembly mechanics:

- A_to_B_rescan
- B_to_A_rescan
- stale_same_scan
- connected/strong/any context
- allow_B_reuse_A_faces on/off

Anti-smuggling constraints:

- delta_beta/H2/harmonic data remain audit-only.
- no i, no global J, no Hodge star, no positivity, no C*-norm, no physical adjoint, no final sym(M).
- stale same-scan is explicitly marked as a legality diagnostic, not as a claimed growth law.
