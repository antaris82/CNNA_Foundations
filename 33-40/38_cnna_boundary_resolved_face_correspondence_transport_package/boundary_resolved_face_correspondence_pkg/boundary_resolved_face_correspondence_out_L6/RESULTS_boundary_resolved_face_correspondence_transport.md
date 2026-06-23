# Results: boundary-resolved face-correspondence transport

This test infers the face correspondence between double-history parent-fan nets
from boundary/port fingerprints rather than from face labels or a prescribed Z3
shift.

```text
CNNA BOUNDARY-RESOLVED FACE-CORRESPONDENCE TRANSPORT
max_level=6, mode=linear, min_omega=1e-10, scale=unit

MODEL SUMMARIES
  real_growth: nodes=1093, cells=364, jnets=120, double_classes=40, three_sector_cycles=40
  symmetrized_birth: nodes=1093, cells=364, jnets=120, double_classes=40, three_sector_cycles=40
  no_backreaction: nodes=1093, cells=364, jnets=120, double_classes=40, three_sector_cycles=40

SELECTED SUMMARIES
  double_history_suffix_cycle | real_growth | live | full | boundary_ports: count=40, valid=1.000, K=0.00186993, id_match=1.000, dir_match=0.000, cyclic=0.000, refl=0.000, best=0.003, id_cost=0.003, dir_cost=0.2453, gap=0.1551, adv_dir_vs_id=-0.2423, cycle_id=1.000, cycle_nontriv=0.000
  missing
  missing
  missing
  double_history_suffix_cycle | real_growth | live | full | K_signed: count=40, valid=1.000, K=0.00186993, id_match=1.000, dir_match=0.000, cyclic=0.000, refl=0.000, best=0.003562, id_cost=0.003562, dir_cost=0.6682, gap=0.0005684, adv_dir_vs_id=-0.6646, cycle_id=1.000, cycle_nontriv=0.000
  missing
  double_history_suffix_cycle | symmetrized_birth | live | full | boundary_ports: count=40, valid=1.000, K=7.08903e-06, id_match=1.000, dir_match=0.000, cyclic=0.000, refl=0.000, best=0.0002156, id_cost=0.0002156, dir_cost=0.2412, gap=0.1602, adv_dir_vs_id=-0.2409, cycle_id=1.000, cycle_nontriv=0.000
  double_history_suffix_cycle | no_backreaction | live | full | boundary_ports: count=40, valid=1.000, K=0.00363105, id_match=1.000, dir_match=0.000, cyclic=0.000, refl=0.000, best=0.003857, id_cost=0.003857, dir_cost=0.2316, gap=0.1436, adv_dir_vs_id=-0.2278, cycle_id=1.000, cycle_nontriv=0.000
  identical_history_clone_control | real_growth | live | full | boundary_ports: count=5, valid=1.000, K=0.000111653, id_match=1.000, dir_match=0.000, cyclic=0.000, refl=0.000, best=0, id_cost=0, dir_cost=0.2052, gap=0.1323, adv_dir_vs_id=-0.2052, cycle_id=1.000, cycle_nontriv=0.000
  random_same_level_cycle_baseline | real_growth | live | full | boundary_ports: count=5, valid=1.000, K=0.00168523, id_match=1.000, dir_match=0.000, cyclic=0.000, refl=0.000, best=0.005364, id_cost=0.005364, dir_cost=0.2417, gap=0.1514, adv_dir_vs_id=-0.2364, cycle_id=1.000, cycle_nontriv=0.000
  double_history_suffix_cycle | real_growth | live | diagonal | boundary_ports: count=40, valid=0.000, K=0, id_match=1.000, dir_match=0.000, cyclic=0.000, refl=0.000, best=0.003007, id_cost=0.003007, dir_cost=0.2468, gap=0.1561, adv_dir_vs_id=-0.2438, cycle_id=1.000, cycle_nontriv=0.000
  double_history_suffix_cycle | real_growth | live | trace_scalar | boundary_ports: count=40, valid=0.000, K=0, id_match=1.000, dir_match=0.000, cyclic=0.000, refl=0.000, best=0.003008, id_cost=0.003008, dir_cost=0.2469, gap=0.1562, adv_dir_vs_id=-0.2439, cycle_id=1.000, cycle_nontriv=0.000
  missing
  missing

READING RULE
  The primary anti-smuggling fingerprint is boundary_ports: it uses DtN
  port/vertex boundary profiles and edge-response magnitudes, but not J/K
  as the matching feature. K/J fingerprints are audits only.
  The script does not restrict correspondences to Z3 label shifts. It tests
  all six S3 permutations and classifies the best one as identity, cyclic
  plus/minus or reflection.
  A positive boundary-resolved nontrivial transport would need real_growth
  double-history cycles to prefer cyclic/reflection permutations over identity,
  with identical-history and random controls not doing the same. Diagonal/trace
  reductions should kill the noncommutative K-sector even if boundary labels
  still have trivial identity structure.
```

## Output files

- `boundary_face_cycle_rows.csv`
- `boundary_face_edge_rows.csv`
- `boundary_face_summary_main.csv`
- `boundary_face_summary_by_level.csv`
- `model_summaries.csv`
- `SUMMARY.txt`

## Status

The primary row is `double_history_suffix_cycle | real_growth | live | full |
boundary_ports`.  It is the cleanest anti-smuggling gate because it excludes
J/K from the correspondence feature.
