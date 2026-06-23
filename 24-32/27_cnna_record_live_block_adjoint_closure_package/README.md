# CNNA record/live block adjoint closure diagnostic

This package contains `test_record_live_block_adjoint_closure.py`.

It follows the port-resolved rank-one adjoint-closure diagnostic and tests whether
the missing adjoint directions become visible when the local boundary space is
typed as

```text
W_p = B_p^record ⊕ B_p^live,   B_p = R^3.
```

The test uses true event-resolved growth inherited from the dynamic birth tests:

- ternary children are born sequentially;
- a newborn senses parent-line plus older siblings;
- the newborn backreacts on parent-line and older siblings;
- descendant shell-loads age older cells.

No complex phase, no `J`, no `J²=-I`, no C*-norm, no GNS representation, and no
AQFT net are assumed or derived.

## Run

```bash
python3 test_record_live_block_adjoint_closure.py --max-level 5 \
  --out record_live_block_adjoint_out_L5 \
  --operator-modes triangular_handoff \
  --longitudinal-modes triangular_record_live \
  --metric-sources record_live_block
```

The script supports larger levels, but the block-envelope least-squares closure
checks are significantly heavier than the previous live-only diagnostics.

## Primary output

```text
record_live_block_adjoint_out_L5/RESULTS_record_live_block_adjoint_closure.md
```
