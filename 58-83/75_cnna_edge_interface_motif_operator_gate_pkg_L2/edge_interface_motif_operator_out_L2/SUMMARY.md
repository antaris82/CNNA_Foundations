# SUMMARY — edge interface motif operator gate

Model label:
CNNA growing primal simplicial complex with deterministic sequential provenance growth, directed antisymmetric birth-transport operators, signed-Jlock two-pair assembly motifs, and a derived shared-edge interface diagnostic.

This test checks whether the previous motif-basis obstruction is caused by missing A/B interface coupling.  The new operator is built on the actual union of faces of each complete A/B assembly:

```text
Face -> shared edge -> Face
```

Primary interface mode:

```text
incidence_identity:
  skew handoff from face to face through the common edge,
  sign = boundary-incidence sign * birth/provenance signature sign.
```

Secondary diagnostic modes:

```text
edge_projector
edge_complement_projector
```

No i, global J, Hodge star, *, positivity, C*-norm, final sym(M), or delta-beta/H2 decision is introduced.  Delta-beta/H2/harmonic quantities are measured after the fact only.

| option | variant | beta | pairs | asm | pair harm | Q harm | P harm | pair J-lock | signed | base union | edge-if lock | edge-if J2 | pass | best mode | scale | used dBeta? |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|---|---:|
| A_to_B_rescan_strong_reuseB | real_growth | (1,0,4,0) | 4 | 1 | 0.381795 | 0.381795 | 0.248801 | 0.448553 | -0.624664 | 0.543826 | 0.596246 | 0.345392 | 0 | edge_projector | unit | False |
| A_to_B_rescan_strong_reuseB | strict_symmetrized_control | (1,0,0,0) | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |  |  | False |
| A_to_B_rescan_strong_reuseB | no_backreaction | (1,0,8,0) | 4 | 2 | 0.405623 | 0.405623 | 0.343367 | 0.339477 | -0.285265 | 0.661108 | 0.689923 | 0.475744 | 0 | edge_projector | unit | False |

## Best non-strict row

```json
{
  "option": "A_to_B_rescan_strong_reuseB",
  "variant": "real_growth",
  "beta": [
    1,
    0,
    4,
    0
  ],
  "pairings": 4,
  "assemblies": 1,
  "pair_harm": 0.38179479579105513,
  "Q_harm": 0.38179479579105513,
  "P_harm": 0.24880104631557917,
  "pair_local_J_lock": 0.4485530110320348,
  "signed_birth": -0.6246641039192367,
  "motif_count": 1,
  "union_motif_lock": 0.5438258343030387,
  "base_union_lock": 0.5438258343030387,
  "edge_interface_lock": 0.596246290239973,
  "edge_interface_J2": 0.3453922919084411,
  "edge_interface_pass": 0,
  "edge_interface_improve_count": 0,
  "best_mode": "edge_projector",
  "best_scale": "unit",
  "used_delta_beta": false
}
```
