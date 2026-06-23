# SUMMARY — assembly motif basis diagonalization gate

Model label:
CNNA growing primal simplicial complex with deterministic sequential provenance growth, directed antisymmetric birth-transport operators, and signed-Jlock two-pair assembly motifs.

This test asks whether the previous J-lock obstruction is a measurement-basis artifact.  Instead of testing each pair separately, each complete A/B assembly is represented by the native motif basis

```text
A_Q, A_P, B_Q, B_P
```

Two derived motif operators are tested:

```text
direct_sum: pair-local block-diagonal J/C on A and B separately
union_sum: shared-face motif operator obtained by summing the two pair-exchange maps on the actual union of faces
```

A positive result would require J to map the motif Q-plane span(A_Q,B_Q) into the motif P-plane span(A_P,B_P), with low span leakage and reasonable projected J^2 = -I behavior.  No i, global J, Hodge, *, positivity, C*-norm, final sym(M), or delta-beta decision is introduced.

| option | variant | beta | pairs | asm | pair harm | Q harm | P harm | pair J-lock | signed | motif n | direct motif lock | union motif lock | direct pass | union pass | used dBeta? |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| A_to_B_rescan_strong_reuseB | real_growth | (1,0,4,0) | 4 | 1 | 0.381795 | 0.381795 | 0.248801 | 0.448553 | -0.624664 | 1 | 0.543826 | 0.543826 | 0 | 0 | False |
| A_to_B_rescan_strong_reuseB | strict_symmetrized_control | (1,0,0,0) | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | False |
| A_to_B_rescan_strong_reuseB | no_backreaction | (1,0,8,0) | 4 | 2 | 0.405623 | 0.405623 | 0.343367 | 0.339477 | -0.285265 | 2 | 0.661108 | 0.661108 | 0 | 0 | False |

## Best rows

```json
{
  "best_direct_sum": {
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
    "direct_sum_motif_lock": 0.5438258343030384,
    "union_sum_motif_lock": 0.5438258343030387,
    "direct_sum_gate_pass": 0,
    "union_sum_gate_pass": 0,
    "used_delta_beta": false
  },
  "best_union_sum": {
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
    "direct_sum_motif_lock": 0.5438258343030384,
    "union_sum_motif_lock": 0.5438258343030387,
    "direct_sum_gate_pass": 0,
    "union_sum_gate_pass": 0,
    "used_delta_beta": false
  }
}
```
