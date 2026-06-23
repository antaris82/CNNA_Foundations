# CNNA Growth Test: parent-fan J-field gluing and holonomy

This package follows `test_parent_fan_plaquette_skew_sector_J_candidate.py`.

It tests whether the local parent-fan plaquette skew-sector candidates

```text
K_abc = [A_ab, A_bc]
J_abc = K_abc / omega
```

glue coherently over the same-parent parent-fan face net.

The test uses:

- real sequential growth and backreaction from the dynamic birth tests,
- provenance-generated parent-fan 2-simplices,
- full DtN Record/Live vertex operators,
- orientation-aware face extraction,
- diagonal/trace/common-commuting reductions as ablations,
- symmetrized-birth and no-backreaction controls.

No global complex structure, C*-norm, GNS representation, AQFT net, or physical `i` is introduced.

## Run

```bash
python3 test_parent_fan_J_field_gluing_and_holonomy.py --max-level 6 --random-nets 30 --outdir parent_fan_J_field_gluing_out_L6
```

## Main files

- `test_parent_fan_J_field_gluing_and_holonomy.py`
- `record_live_block_base.py`
- `parent_fan_J_field_gluing_out_L6/RESULTS_parent_fan_J_field_gluing_and_holonomy.md`
- `parent_fan_J_field_gluing_out_L6/SUMMARY.txt`
- CSV detail and summary tables in the output directory.
