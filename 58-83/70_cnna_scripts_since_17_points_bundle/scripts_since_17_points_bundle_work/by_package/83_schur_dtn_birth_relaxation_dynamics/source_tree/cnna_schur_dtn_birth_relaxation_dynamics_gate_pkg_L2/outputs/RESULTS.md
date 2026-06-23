# RESULTS — CNNA Schur/DtN birth + live relaxation dynamics gate, L2

## Model provenance

This package keeps the established script-1/script-2 dynamic birth conductance model and the true Schur/DtN boundary-role calculation.  It adds fixed-topology live relaxation after every birth event.

The split is explicit:

```text
birth_delta_DtN:
  topology and boundary roles change;
  newborn becomes a new UV-tail port for parent/ancestor cuts.

relax_delta_DtN:
  no new node is born;
  boundary ports stay fixed;
  live conductances and live edge weights relax;
  true Schur/DtN matrices can drift further.

record_vs_live_gap:
  difference between the immutable after-birth DtN record and the current live DtN state.
```

The live relaxation rule is deterministic and local: it uses only existing `birth_g`, current `g`, and incoming/outgoing edge response loads inherited from the script-1/script-2 growth model.  It does not use J, i, Hodge, physical star, positivity, a C*-norm, or delta-beta.

## Gate table

| variant | events | relax events | birth Schur/DtN gate | live relax gate | retarded | advanced | mean birth ΔDtN | first relax ΔDtN | last relax ΔDtN | total relax drift | record/live gap | last/first | relax active | quasi-equilibrium | full Markov complex |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| real_growth_linear_script1_2_true_schur_birth_plus_relax | 12 | 12 | 1 | 1 | 1.000 | 0.000 | 0.178083 | 0.0133824 | 0.00066569 | 0.023987 | 0.0214332 | 0.048 | 0.694 | 1.000 | 1.000 |
| log_growth_script1_2_true_schur_birth_plus_relax | 12 | 12 | 1 | 1 | 1.000 | 0.000 | 0.134361 | 0.00910151 | 0.000443109 | 0.0165773 | 0.0150021 | 0.049 | 0.694 | 1.000 | 1.000 |
| saturating_growth_script1_2_true_schur_birth_plus_relax | 12 | 12 | 1 | 1 | 1.000 | 0.000 | 0.557052 | 0.0940137 | 0.0103495 | 0.233961 | 0.233433 | 0.111 | 0.694 | 1.000 | 1.000 |
| strict_symmetrized_response_control_true_schur_birth_plus_relax | 12 | 12 | 0 | 0 | 0.000 | 0.000 | 0 | 0 | 0 | 0 | 0 | 0.000 | 0.000 | 0.000 | 0.000 |
| birth_only_no_relax_record_control_true_schur | 12 | 12 | 1 | 0 | 1.000 | 0.000 | 0.185046 | 0 | 0 | 0 | 0 | 0.000 | 0.000 | 0.000 | 1.000 |

## Interpretation

The birth-event Schur/DtN role gate remains positive in the nontrivial growth variants: the newborn still has no own UV-tail from its self-perspective, while parent/ancestor cuts receive a new UV-tail boundary port and a true Schur/DtN response change.

The new result is the fixed-topology live layer.  In the nontrivial growth variants, the after-birth DtN record does not remain equal to the later live DtN state.  Relaxation at fixed boundary produces nonzero `relax_delta_DtN` and nonzero `record_vs_live_gap`.  The strict symmetrized response control keeps the topology birth role but kills the response/live relaxation layer.

So the corrected interpretation is:

```text
Birth creates boundary-role/port changes.
Relaxation then changes live Schur/DtN responses without new births.
The network is not static during growth; birth records and live responses separate.
```

This supports the two-layer reading:

```text
Record layer:
  immutable birth/role/DtN record at the moment of boundary-role creation.

Live layer:
  continuing conductance/response relaxation under existing directed response structure.
```

## Limits

- This is still an L2 finite approximant.
- The DtN computation is a true matrix Schur complement of the real conductance Laplacian.
- The relaxation rule is deterministic and local, but it is a model choice within the existing script-1/script-2 quantities.  It is not yet a theorem forcing one unique relaxation law.
- The package does not claim J, i, spin, a *-algebra, positivity, or Q/P compatibility.

## Next test

`test_role_recalibration_to_boundary_value_bridge_true_schur_live_gate.py`

Use the event-resolved true Schur/DtN birth records and live-relaxed Schur/DtN states as separate inputs.  Test whether the difference `record_vs_live_gap` carries a robust retarded/advanced boundary-value polarity that can later bridge to Q/P or an operator adjunction.
