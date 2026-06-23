# CNNA simplicial DtN plaquette frustration test

First bridge test after the pure-provenance operator-closure obstruction.

It combines:

- event-resolved deterministic growth from tests 1/2,
- provenance-generated 2-simplices from test 3,
- the test-9 obstruction that minimal axis geometry alone is flat,
- dynamic Record/Live DtN germs from tests 22/23,
- the post-23 lesson that node/edge provenance alone remains gradient-like.

The script does not set `i`, `J`, a Hilbert space, a `*`-algebra, or a `C*` norm.
It tests only whether response/DtN operators on provenance-generated faces carry non-telescoping plaquette frustration.

Run:

```bash
python test_simplicial_dtn_plaquette_frustration.py --max-level 6 --mode linear --step 0.04 --random-faces 80 --outdir simplicial_dtn_plaquette_out_L6
```

Primary output:

- `simplicial_dtn_plaquette_out_L6/RESULTS_simplicial_dtn_plaquette_frustration.md`
- `simplicial_dtn_plaquette_out_L6/SUMMARY.txt`
- CSV tables for all face rows and grouped summaries.
