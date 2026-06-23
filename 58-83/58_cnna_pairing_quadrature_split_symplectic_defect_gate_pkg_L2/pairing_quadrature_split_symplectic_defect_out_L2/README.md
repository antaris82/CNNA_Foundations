# Pairing quadrature split symplectic defect gate

Run:

```bash
python3 test_pairing_quadrature_split_symplectic_defect_gate.py
```

Outputs:

- comparative_summary.json
- comparative_pairing_quadrature_split_summary.csv
- RESULTS.md
- SUMMARY.md
- SOURCE_AUDIT.md
- per-variant pair and face support logs

Default L2 run uses variants:

```text
real_growth
strict_symmetrized_control
no_backreaction
```

Next test suggested by this package:
`test_real_quadrature_operator_closure_gate.py`, but only if P survives with a
nontrivial harmonic component and strict_sym remains killed.
