# CNNA nonflat connection from same-label response mismatch

This package tests the next gate after the boundary-resolved face-correspondence transport test.

The previous package inferred identity/same-label gluing between double-history parent-fan nets from boundary/port/DtN data.  This script therefore stops searching over face labels and asks whether the remaining same-label response mismatch itself defines a nonflat connection around a double-history root-sector cycle.

The primary loop is, for each suffix class and each face label `F12`, `F23`, `F31`:

```text
root sector 1 -> root sector 2 -> root sector 3 -> root sector 1
```

The script measures:

- raw same-label mismatch of local parent-fan J sectors;
- amplitude, plane, metric and DtN mismatch;
- exact metric-identifying transport as a gauge audit;
- Wilson-loop diagnostics built from response-mismatch generators themselves.

No complex scalar, no physical `i`, no global `J`, no `*`-algebra, no C*-norm, no GNS representation and no AQFT net are introduced.

## Run

```bash
python3 test_nonflat_connection_from_response_mismatch.py --max-level 6 --metric-modes mean_face --random-cycles 40 --identical-cycles 40 --outdir nonflat_response_connection_out_L6
```

A faster exploratory run:

```bash
python3 test_nonflat_connection_from_response_mismatch.py --max-level 6 --quick-suite --random-cycles 40 --identical-cycles 40 --outdir nonflat_response_connection_out_L6_quick
```

## Main outputs

- `response_connection_loop_rows.csv`
- `response_connection_edge_rows.csv`
- `response_connection_summary_main.csv`
- `response_connection_summary_by_label.csv`
- `response_connection_summary_by_level.csv`
- `model_summaries.csv`
- `SUMMARY.txt`
