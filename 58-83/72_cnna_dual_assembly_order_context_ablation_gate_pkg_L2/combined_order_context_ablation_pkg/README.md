# CNNA dual assembly order/context ablation gate

Run connected ablation:

```bash
python3 test_dual_assembly_order_context_ablation_gate.py \
  --orders A_to_B_rescan B_to_A_rescan stale_same_scan \
  --context-modes connected \
  --reuse-modes reuseB noReuseB \
  --max-eval-candidates 80
```

Run strong-context ablation:

```bash
python3 test_dual_assembly_order_context_ablation_gate.py \
  --orders A_to_B_rescan B_to_A_rescan stale_same_scan \
  --context-modes strong \
  --reuse-modes reuseB noReuseB \
  --variants real_growth no_backreaction \
  --max-eval-candidates 50 \
  --out dual_assembly_order_context_ablation_strong_L2 \
  --zip cnna_dual_assembly_order_context_ablation_strong_gate_pkg_L2.zip
```

Outputs include scripts, per-option logs, comparative JSON/CSV, RESULTS.md, SUMMARY.md, and SOURCE_AUDIT.md.
