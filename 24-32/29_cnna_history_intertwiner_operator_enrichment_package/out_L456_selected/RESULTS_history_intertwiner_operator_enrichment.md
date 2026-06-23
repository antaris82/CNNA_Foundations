# RESULTS: history-intertwiner operator-enrichment closure

## Status

This is still not a J-test.  It checks whether provenance-derived history intertwiners enrich the visible L/T operator system enough to approach DtN-adjoint and finite product closure.

## Parameters

- levels: `[4, 5, 6]`
- tiers: `['base', 'completion_live', 'all_enriched']`
- degree: `2`
- max_pairs_per_control: `1`
- word_cap: `60`

## Primary double-history scale table

| L | tier | dim/full | star basis mean | star seed mean | mult mean | note |
|---:|---|---:|---:|---:|---:|---|
| 4 | `base` | 0.194 | 0.0702 | 0.00191 | 0.474 | partial adjoint visibility |
| 4 | `completion_live` | 0.257 | 0.677 | 0.307 | 0.593 | not closed |
| 4 | `all_enriched` | 0.299 | 0.651 | 0.327 | 0.463 | not closed |
| 5 | `base` | 0.194 | 0.868 | 0.436 | 0.368 | not closed |
| 5 | `completion_live` | 0.257 | 0.922 | 0.371 | 0.324 | not closed |
| 5 | `all_enriched` | 0.299 | 0.934 | 0.406 | 0.512 | not closed |
| 6 | `base` | 0.188 | 0.0759 | 0.00159 | 0.381 | partial adjoint visibility |
| 6 | `completion_live` | 0.257 | 0.696 | 0.3 | 0.722 | not closed |
| 6 | `all_enriched` | 0.299 | 0.644 | 0.317 | 0.569 | not closed |

## Interpretation

- If enrichment only works by saturating most of the 12x12 matrix space, that is not meaningful *-closure; it is an overlarge envelope.
- If residuals improve with L at comparable span dimension, that supports the idea that the *-structure is a local limiting phenomenon rather than present in the small finite window.
- If residuals do not improve with L, the current provenance-derived intertwiners are still missing necessary directions.
