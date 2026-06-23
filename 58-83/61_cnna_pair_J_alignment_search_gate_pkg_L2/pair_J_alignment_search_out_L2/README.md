# Pair J-alignment search gate

Run:

```bash
python3 test_pair_J_alignment_search_gate.py
```

Outputs:

- comparative_summary.json
- comparative_alignment_search_summary.csv
- comparative_candidate_summary.csv
- phase_flip_comparison.csv
- RESULTS.md
- SUMMARY.md
- SOURCE_AUDIT.md
- per-variant candidate logs

Positive result requires a derived candidate pair aligned by local J_pair and killed by
strict_sym.  A local algebra alone is not sufficient.
