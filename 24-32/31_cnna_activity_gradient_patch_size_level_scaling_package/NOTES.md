# Notes

## Key distinction

The diagnostic separates two notions of change:

- **Last-shell marginal update:** compares the same provenance node before and after the latest growth shell. This is the clean test of locality of current growth. It confirms that frontier/active-parent nodes change most and the root least.
- **Cumulative aging since completion:** compares a completed local cell's current live DtN data to its completion snapshot. This can be large in older interior cells because memory/backreaction accumulates over many later births. This is not the same as instantaneous local growth activity.

## Interpretation

The results support a two-scale reading:

- current update is local/frontier-dominated;
- old interior cells store accumulated memory;
- a *-operator candidate should not be sought on a root-only or isolated old-pair subsystem, but on growth-defined local patches and eventually a local-net/limit construction.

## Current status

Positive:

- last-shell marginal update gradient is clean and monotone toward the active parent frontier;
- same-suffix growth patches remain structured and nonrandom.

Negative/open:

- finite degree-2 patch closure is still not a stable *-algebra;
- patch-size/level trend is not yet sufficient as a limit signal;
- no J-candidate is extracted here.

## Next test

`test_active_shell_patch_net_refinement.py`

Goal: restrict patch construction to active and near-active shells, then compare closure residuals across local refinements while treating old interior cells as accumulated memory/background rather than as the primary carrier.
