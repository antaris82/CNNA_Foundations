# Signed quadrature area kappa gate

Run:

```bash
python3 test_signed_quadrature_area_kappa_gate.py
```

Outputs:

- comparative_summary.json
- comparative_signed_quadrature_area_summary.csv
- phase_flip_comparison.csv
- RESULTS.md
- SUMMARY.md
- SOURCE_AUDIT.md
- per-variant signed pair/face logs

Next gate depends on this result:

- if signed area is small/non-flipping: treat Q/P as magnitude-only and search for a
  native operator involution/quadrature pairing instead of spatial orientation;
- if signed area is nontrivial and flips: build a strict derived real symplectic-form
  closure test.
