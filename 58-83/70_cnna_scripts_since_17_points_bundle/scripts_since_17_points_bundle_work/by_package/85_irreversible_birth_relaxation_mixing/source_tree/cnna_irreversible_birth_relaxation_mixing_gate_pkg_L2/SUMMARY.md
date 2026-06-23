# RESULTS: CNNA irreversible birth-relaxation mixing gate
## Model class
Established script-1/script-2 ternary sequential growth with true Schur/DtN matrices and fixed-topology live relaxation.  This is still a finite approximant; it is not an infinite-growth theorem.
## Primary question
Does the newborn's effect become irreversible because it is distributed into old nodes/old edges/old-cut DtN maps and then mixed with live relaxation?
## Headline L4
- `real_growth_irreversible_nonlocal_residue_fraction`: `1.0`
- `real_growth_live_relaxation_mixing_fraction`: `1.0`
- `real_growth_child_delete_recovers_old_state_fraction`: `0.0`
- `real_growth_mean_record_vs_live_gap_fro`: `0.030586228979614175`
- `real_growth_mean_distributed_old_diff_fraction`: `0.016708549511684538`
- `birth_only_irreversible_nonlocal_residue_fraction`: `1.0`
- `birth_only_live_relaxation_mixing_fraction`: `0.0`
- `pure_topology_child_delete_recovers_old_state_fraction`: `1.0`
- `strict_sym_child_delete_recovers_old_state_fraction`: `1.0`
- `saturating_mean_record_vs_live_gap_fro`: `0.2811681354847869`

## Interpretation
- Pure topology without response is reversible after deleting the newborn: this control separates mere node insertion from response dynamics.
- Birth-only response is already nonrecoverable: newborn backreaction changes old conductances/edges and old-cut Schur/DtN maps.
- Birth plus live relaxation adds a nonzero record/live gap: the Birth-Record is no longer the Live-State after fixed-topology relaxation.
- Therefore the growth process is not reversible by deleting the newest node. Its effect has been distributed into the old network and mixed with relaxation.

## Methodological status
No J, i, Hodge, physical star, positivity/C*-norm, Q/P target, or delta-beta gate is used.  DtN maps are true real Laplace Schur complements.

## Next test
`test_irreversible_live_semigroup_boundary_value_gate.py`: use the nonrecoverable live Schur/DtN evolution as a semigroup candidate and test retarded/advanced polarity under reverse, kappa and longitudinal controls.
