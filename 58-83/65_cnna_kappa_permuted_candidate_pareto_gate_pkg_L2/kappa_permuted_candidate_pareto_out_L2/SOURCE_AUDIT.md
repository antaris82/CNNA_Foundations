# SOURCE AUDIT

This test is a derived-only candidate audit.

Anti-smuggling constraints:

- no complex scalar and no `i`;
- no Hodge star, no positivity, no physical adjoint;
- no final `sym(M)` in the directed birth-transport vertex operator;
- no arbitrary fitted rotation;
- no delta_beta/H2/kappa used as a move decision input;
- kappa is implemented as the concrete sibling birth-order reflection 1<->3 on `Node.birth_order`, not as a mere phase-sign flag.

Limitation:
The kappa audit keeps geometry, conductances, directed_edges and node IDs fixed to allow exact candidate matching.  It is stronger than a phase-sign flip but weaker than a full re-growth under reversed birth sequence.  A full address-permuted re-growth would require a separate address-level candidate matcher.
