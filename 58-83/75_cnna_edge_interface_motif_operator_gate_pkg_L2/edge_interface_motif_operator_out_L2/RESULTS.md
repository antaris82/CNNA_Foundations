# RESULTS — edge interface motif operator gate

## Comparative table

| option | variant | beta | pairs | asm | pair harm | Q harm | P harm | pair J-lock | signed | base union | edge-if lock | edge-if J2 | pass | best mode | scale | used dBeta? |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|---|---:|
| A_to_B_rescan_strong_reuseB | real_growth | (1,0,4,0) | 4 | 1 | 0.381795 | 0.381795 | 0.248801 | 0.448553 | -0.624664 | 0.543826 | 0.596246 | 0.345392 | 0 | edge_projector | unit | False |
| A_to_B_rescan_strong_reuseB | strict_symmetrized_control | (1,0,0,0) | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |  |  | False |
| A_to_B_rescan_strong_reuseB | no_backreaction | (1,0,8,0) | 4 | 2 | 0.405623 | 0.405623 | 0.343367 | 0.339477 | -0.285265 | 0.661108 | 0.689923 | 0.475744 | 0 | edge_projector | unit | False |

## Gate criterion

A positive edge-interface result requires the effective operator

```text
J_eff = J_pair_union + lambda * J_edge_interface
```

to improve the motif Q/P subspace lock and give acceptable projected J^2 behavior, while strict_sym stays null and used_delta_beta remains false.

The tested interface signs come only from oriented boundary incidence of faces on their shared edge plus birth/provenance signatures.  The primary `incidence_identity` mode avoids Hodge and uses no metric edge rotation.  Edge-direction projector modes are included only as secondary geometric diagnostics.

## Interpretation rules

- If `edge_interface_lock` drops strongly below `base_union_lock` and projected J2 improves, the obstruction was likely a missing derived edge-interface handoff.
- If Q/P and beta2 remain positive but `edge_interface_lock` stays high, then merely adding incidence-level edge coupling is not enough.
- If strict_sym is zero, the path remains tied to nonsymmetric provenance growth.
