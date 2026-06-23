# CNNA diagnostic: port-resolved rank-one adjoint closure

This package contains a numerical diagnostic for the growing real complement-network path.

It follows the real local operator-calculus and 2D role-plane adjoint tests.  The purpose is to keep the full 3-port DtN boundary structure longer, rather than compressing each local cell immediately to an Env/UV role plane.

Main script:

```bash
python test_port_resolved_rank_one_adjoint_closure.py --max-level 6 --out port_resolved_rank_one_adjoint_out_L6
```

Primary output:

```text
port_resolved_rank_one_adjoint_out_L6/RESULTS_port_resolved_rank_one_adjoint_closure.md
```

The diagnostic does not derive `J`, does not test `J²=-I`, and does not claim a C*-algebra or GNS representation.

It tests whether rank-one transverse handoff maps between full 3-port local DtN boundary spaces remain visible under the candidate adjoint defined by the full DtN-energy pairing.
