# Results: boundary-resolved face-correspondence transport

This test infers the face correspondence between double-history parent-fan nets
from boundary/port fingerprints rather than from face labels or a prescribed Z3
shift.

```text
CNNA BOUNDARY-RESOLVED FACE-CORRESPONDENCE TRANSPORT
max_level=5, mode=linear, min_omega=1e-10, scale=unit

MODEL SUMMARIES
  real_growth: nodes=364, cells=121, jnets=39, double_classes=13, three_sector_cycles=13
  symmetrized_birth: nodes=364, cells=121, jnets=39, double_classes=13, three_sector_cycles=13
  no_backreaction: nodes=364, cells=121, jnets=39, double_classes=13, three_sector_cycles=13

SELECTED SUMMARIES
  double_history_suffix_cycle | real_growth | live | full | boundary_ports: count=13, valid=1.000, K=0.00156178, id_match=1.000, dir_match=0.000, cyclic=0.000, refl=0.000, best=0.003439, id_cost=0.003439, dir_cost=0.2466, gap=0.1557, adv_dir_vs_id=-0.2431, cycle_id=1.000, cycle_nontriv=0.000
  double_history_suffix_cycle | real_growth | live | full | boundary_ports_no_edges: count=13, valid=1.000, K=0.00156178, id_match=1.000, dir_match=0.000, cyclic=0.000, refl=0.000, best=0.003249, id_cost=0.003249, dir_cost=0.2475, gap=0.1564, adv_dir_vs_id=-0.2442, cycle_id=1.000, cycle_nontriv=0.000
  double_history_suffix_cycle | real_growth | live | full | record_live_boundary: count=13, valid=1.000, K=0.00156178, id_match=1.000, dir_match=0.000, cyclic=0.000, refl=0.000, best=0.003602, id_cost=0.003602, dir_cost=0.1924, gap=0.1211, adv_dir_vs_id=-0.1888, cycle_id=1.000, cycle_nontriv=0.000
  double_history_suffix_cycle | real_growth | live | full | face_dtn_full: count=13, valid=1.000, K=0.00156178, id_match=1.000, dir_match=0.000, cyclic=0.000, refl=0.000, best=0.003122, id_cost=0.003122, dir_cost=0.01948, gap=0.00531, adv_dir_vs_id=-0.01636, cycle_id=1.000, cycle_nontriv=0.000
  double_history_suffix_cycle | real_growth | live | full | K_signed: count=13, valid=1.000, K=0.00156178, id_match=1.000, dir_match=0.000, cyclic=0.000, refl=0.000, best=0.004999, id_cost=0.004999, dir_cost=0.6684, gap=0.0002631, adv_dir_vs_id=-0.6634, cycle_id=1.000, cycle_nontriv=0.000
  double_history_suffix_cycle | real_growth | live | full | K_signless: count=13, valid=1.000, K=0.00156178, id_match=0.795, dir_match=0.000, cyclic=0.000, refl=0.205, best=0.004999, id_cost=0.004999, dir_cost=0.005235, gap=1.136e-05, adv_dir_vs_id=-0.0002362, cycle_id=0.846, cycle_nontriv=0.154
  double_history_suffix_cycle | symmetrized_birth | live | full | boundary_ports: count=13, valid=1.000, K=6.4518e-06, id_match=1.000, dir_match=0.000, cyclic=0.000, refl=0.000, best=0.0002609, id_cost=0.0002609, dir_cost=0.2465, gap=0.1637, adv_dir_vs_id=-0.2462, cycle_id=1.000, cycle_nontriv=0.000
  double_history_suffix_cycle | no_backreaction | live | full | boundary_ports: count=13, valid=1.000, K=0.00321731, id_match=1.000, dir_match=0.000, cyclic=0.000, refl=0.000, best=0.004473, id_cost=0.004473, dir_cost=0.234, gap=0.1448, adv_dir_vs_id=-0.2295, cycle_id=1.000, cycle_nontriv=0.000
  identical_history_clone_control | real_growth | live | full | boundary_ports: count=10, valid=1.000, K=0.000137152, id_match=1.000, dir_match=0.000, cyclic=0.000, refl=0.000, best=0, id_cost=0, dir_cost=0.2144, gap=0.138, adv_dir_vs_id=-0.2144, cycle_id=1.000, cycle_nontriv=0.000
  random_same_level_cycle_baseline | real_growth | live | full | boundary_ports: count=10, valid=1.000, K=0.000971025, id_match=1.000, dir_match=0.000, cyclic=0.000, refl=0.000, best=0.005912, id_cost=0.005912, dir_cost=0.2349, gap=0.1468, adv_dir_vs_id=-0.229, cycle_id=1.000, cycle_nontriv=0.000
  double_history_suffix_cycle | real_growth | live | diagonal | boundary_ports: count=13, valid=0.000, K=0, id_match=1.000, dir_match=0.000, cyclic=0.000, refl=0.000, best=0.003445, id_cost=0.003445, dir_cost=0.2482, gap=0.1568, adv_dir_vs_id=-0.2447, cycle_id=1.000, cycle_nontriv=0.000
  double_history_suffix_cycle | real_growth | live | trace_scalar | boundary_ports: count=13, valid=0.000, K=0, id_match=1.000, dir_match=0.000, cyclic=0.000, refl=0.000, best=0.003445, id_cost=0.003445, dir_cost=0.2482, gap=0.1568, adv_dir_vs_id=-0.2448, cycle_id=1.000, cycle_nontriv=0.000
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
