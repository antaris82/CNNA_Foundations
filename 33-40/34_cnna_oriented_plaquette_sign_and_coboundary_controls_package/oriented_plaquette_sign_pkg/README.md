# CNNA oriented plaquette sign and coboundary controls

This package audits the first simplicial DtN plaquette-frustration signal.

It keeps the same event-resolved ternary growth and dynamic Record/Live DtN data, but now tests whether the face signal behaves as an oriented noncommutative plaquette effect:

1. reversing `(a,b,c)` to `(a,c,b)` must flip the commutator/skew sign;
2. first-order edge circulation must telescope;
3. diagonal, trace-scalar, and common-basis commuting ablations must kill the signal;
4. sibling and parent-fan faces are kept separate.

No `i`, no `J`, no Hilbert space, no C*-norm, and no *-algebra are introduced.

## Run

```bash
python3 test_oriented_plaquette_sign_and_coboundary_controls.py --max-level 6 --outdir oriented_plaquette_sign_out_L6
```

## Main output

- `oriented_pair_rows.csv`
- `orientation_summary_main.csv`
- `orientation_summary_by_face_level.csv`
- `orientation_summary_by_source.csv`
- `coboundary_ablation_summary.csv`
- `model_summaries.csv`
- `SUMMARY.txt`
- `RESULTS_oriented_plaquette_sign_and_coboundary_controls.md`
