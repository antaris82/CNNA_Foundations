# CNNA history-intertwiner operator-enrichment diagnostic

This package continues the pre-* chain after `generated_real_operator_system_closure`.

It checks whether provenance-derived history intertwiners enrich the local L/T Record/Live block operator system enough to improve DtN-adjoint closure, without prematurely testing or setting `J`.

Run:

```bash
python3 test_history_intertwiner_operator_enrichment.py --out history_intertwiner_operator_enrichment_out --levels 4 5 6 --max-pairs-per-control 8
```

Main output:

- `RESULTS_history_intertwiner_operator_enrichment.md`
- `summary_table_all.csv`
- `pair_rows_all.csv`

Interpretation:

- positive only if adjoint/product residuals improve without near full-space saturation;
- saturation of the full 12x12 matrix space is an overlarge envelope, not meaningful closure;
- improvement with level at comparable dimension would support the hypothesis that a real *-structure is a local limiting object.
