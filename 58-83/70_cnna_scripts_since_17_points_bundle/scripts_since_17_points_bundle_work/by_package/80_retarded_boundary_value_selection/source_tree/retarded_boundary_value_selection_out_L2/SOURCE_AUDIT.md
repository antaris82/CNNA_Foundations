# SOURCE AUDIT

This package is derived from `test_bounded_response_bridge_sanity_gate.py` and reuses its bounded bridge rows, but changes the primary question.

It does not:

- set `i`, global `J`, Hodge star, physical adjoint, positivity, C*-norm, or a target spin signature;
- use `delta_beta`/H2 as a selection input;
- count Q/P-lock or `J²=-I` as primary success.

It does:

- distinguish retarded/forward response operators from advanced/reverse and symmetrized cycle operators;
- distinguish longitudinal/radial bridge modes from transverse/birth-q bridge modes;
- audit bounded contraction, bridge conditioning, half-plane bias, and response-circulation directional bias;
- keep strict-symmetrized control as a null test.

Model label: CQNM/s=-1 inspired growing real complement-network diagnostic.  The growth order / longitudinal axis is treated as a pre-causal provenance structure, not as physical time.
