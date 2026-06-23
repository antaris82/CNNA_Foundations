# Pairing quadrature adjoint-pairing gate

Run:

```bash
python3 test_pairing_quadrature_adjoint_pairing_gate.py
```

Outputs:

- comparative_summary.json
- comparative_adjoint_pairing_summary.csv
- phase_flip_comparison.csv
- RESULTS.md
- SUMMARY.md
- SOURCE_AUDIT.md
- per-variant pair/face logs

Positive result would require more than local `J_pair^2=-I`: the actual Q/P channels
must be locked by the pair algebra and the signed commutator should transform coherently
under the phase/kappa proxy while strict_sym remains killed.
