# RESULTS — dual assembly order/context ablation gate

## Comparative table

| option | variant | β2 | pairs | asm | pair_harm | Q | P | J-lock | signed | used Δβ |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| A_to_B_rescan_connected_reuseB | real_growth | 4 | 4 | 1 | 0.382 | 0.382 | 0.249 | 0.449 | -0.625 | False |
| A_to_B_rescan_connected_reuseB | strict_symmetrized_control | 0 | 0 | 0 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | False |
| A_to_B_rescan_connected_reuseB | no_backreaction | 6 | 4 | 2 | 0.303 | 0.303 | 0.305 | 0.209 | 0.070 | False |
| A_to_B_rescan_connected_noReuseB | real_growth | 4 | 4 | 1 | 0.382 | 0.382 | 0.249 | 0.449 | -0.625 | False |
| A_to_B_rescan_connected_noReuseB | strict_symmetrized_control | 0 | 0 | 0 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | False |
| A_to_B_rescan_connected_noReuseB | no_backreaction | 6 | 4 | 2 | 0.303 | 0.303 | 0.305 | 0.209 | 0.070 | False |
| B_to_A_rescan_connected_reuseB | real_growth | 4 | 4 | 1 | 0.340 | 0.340 | 0.285 | 0.371 | -0.618 | False |
| B_to_A_rescan_connected_reuseB | strict_symmetrized_control | 0 | 0 | 0 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | False |
| B_to_A_rescan_connected_reuseB | no_backreaction | 4 | 4 | 1 | 0.327 | 0.327 | 0.263 | 0.334 | -0.562 | False |
| B_to_A_rescan_connected_noReuseB | real_growth | 4 | 4 | 1 | 0.340 | 0.340 | 0.285 | 0.371 | -0.618 | False |
| B_to_A_rescan_connected_noReuseB | strict_symmetrized_control | 0 | 0 | 0 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | False |
| B_to_A_rescan_connected_noReuseB | no_backreaction | 4 | 4 | 1 | 0.327 | 0.327 | 0.263 | 0.334 | -0.562 | False |
| stale_same_scan_connected_reuseB | real_growth | 6 | 4 | 1 | 0.252 | 0.252 | 0.333 | 0.326 | -0.242 | False |
| stale_same_scan_connected_reuseB | strict_symmetrized_control | 0 | 0 | 0 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | False |
| stale_same_scan_connected_reuseB | no_backreaction | 4 | 4 | 0 | 0.345 | 0.345 | 0.294 | 0.335 | -0.410 | False |
| stale_same_scan_connected_noReuseB | real_growth | 6 | 4 | 1 | 0.252 | 0.252 | 0.333 | 0.326 | -0.242 | False |
| stale_same_scan_connected_noReuseB | strict_symmetrized_control | 0 | 0 | 0 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | False |
| stale_same_scan_connected_noReuseB | no_backreaction | 4 | 4 | 0 | 0.345 | 0.345 | 0.294 | 0.335 | -0.410 | False |
| A_to_B_rescan_strong_reuseB | real_growth | 4 | 4 | 1 | 0.382 | 0.382 | 0.249 | 0.449 | -0.625 | False |
| A_to_B_rescan_strong_reuseB | no_backreaction | 6 | 4 | 2 | 0.303 | 0.303 | 0.305 | 0.209 | 0.070 | False |
| A_to_B_rescan_strong_noReuseB | real_growth | 4 | 4 | 1 | 0.382 | 0.382 | 0.249 | 0.449 | -0.625 | False |
| A_to_B_rescan_strong_noReuseB | no_backreaction | 6 | 4 | 2 | 0.303 | 0.303 | 0.305 | 0.209 | 0.070 | False |
| B_to_A_rescan_strong_reuseB | real_growth | 4 | 4 | 1 | 0.340 | 0.340 | 0.285 | 0.371 | -0.618 | False |
| B_to_A_rescan_strong_reuseB | no_backreaction | 4 | 4 | 1 | 0.327 | 0.327 | 0.263 | 0.334 | -0.562 | False |
| B_to_A_rescan_strong_noReuseB | real_growth | 4 | 4 | 1 | 0.340 | 0.340 | 0.285 | 0.371 | -0.618 | False |
| B_to_A_rescan_strong_noReuseB | no_backreaction | 4 | 4 | 1 | 0.327 | 0.327 | 0.263 | 0.334 | -0.562 | False |
| stale_same_scan_strong_reuseB | real_growth | 6 | 4 | 1 | 0.252 | 0.252 | 0.333 | 0.326 | -0.242 | False |
| stale_same_scan_strong_reuseB | no_backreaction | 4 | 4 | 0 | 0.345 | 0.345 | 0.294 | 0.335 | -0.410 | False |
| stale_same_scan_strong_noReuseB | real_growth | 6 | 4 | 1 | 0.252 | 0.252 | 0.333 | 0.326 | -0.242 | False |
| stale_same_scan_strong_noReuseB | no_backreaction | 4 | 4 | 0 | 0.345 | 0.345 | 0.294 | 0.335 | -0.410 | False |

## Interpretation

The ablation confirms that the assembly result is not a mere stale same-scan artifact: dynamic A->B and B->A rescan paths both create nontrivial Q/P and H2 carriers, and strict_sym remains killed in the connected run.

However, the key quantities still split:

```text
best signed_birth:
  A_to_B_rescan real_growth
  signed_birth ≈ -0.625
  J-lock ≈ 0.449

best J-lock:
  A_to_B_rescan no_backreaction
  J-lock ≈ 0.209
  signed_birth ≈ 0.070

B_to_A:
  signed_birth remains large,
  but J-lock does not improve enough.

stale_same_scan:
  improves beta2 in real_growth,
  but is only a legality diagnostic, not a preferred dynamic rule.
```

The strong-context results match the connected-context main pattern for the nontrivial variants: the same A->B options dominate the carrier metrics, and face reuse does not materially change L2 results.

## Next test

```text
test_signed_Jlock_role_coupling_gate.py
```

Purpose:

```text
Do not add more pairings.
Do not change topology target.
Instead explicitly couple the selected B-role to signed orientation and J-lock jointly,
and test whether signed_birth and J-lock can be co-optimized without delta_beta/H2 input.
```

Candidate ranking should remain derived-only:

```text
A role:
  Q/P carrier proxies + provenance asymmetry

B role:
  C-lock residual + signed κ/ birth signed amplitude jointly
  with shared-face/edge context to A

Audit only:
  beta2, H2, harmonic, delta_beta
```

If this fails, the obstruction is no longer pair availability or order legality; it is a conflict between signed orientation and J-lock in the current two-pair assembly motif.
