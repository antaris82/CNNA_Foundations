# SUMMARY — kappa-permuted candidate Pareto gate

Model label:
CNNA growing primal simplicial complex with deterministic sequential provenance growth, A-gated complement-pair candidate space, directed antisymmetric birth-transport operators, local C/J pair algebra, and explicit sibling birth-order reflection `1 <-> 3` as a model-level kappa audit.

This package replaces the earlier internal `phase_sign +/-1` audit with a concrete model-label transformation:

```text
identity model:       original birth_order labels
kappa-reflected model: Node.birth_order is reflected inside every sibling fan, 1<->3, 2->2
```

The reflection is applied to the model's birth_order fields, not to a standalone phase-sign parameter.  Geometry, conductances, directed_edges and vertex IDs are kept fixed so the same face-pair candidates can be matched exactly.

| variant | matched | A both | beta2 audit both | C-lock pass | kappa flip pass | all pass | all+beta2 | best C lock | best C flip | used dBeta? |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| real_growth | 206 | 206 | 88 | 19 | 23 | 3 | 0 | 0.128002 | 0.328323 | False |
| strict_symmetrized_control | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | False |
| no_backreaction | 236 | 236 | 108 | 15 | 28 | 2 | 0 | 0.0667835 | 1 | False |

Read conservatively: this is a label-kappa audit on the same grown model, not a full re-growth with reversed birth sequence.
