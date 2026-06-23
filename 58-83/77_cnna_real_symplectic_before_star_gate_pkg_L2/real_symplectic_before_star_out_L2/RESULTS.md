# RESULTS — real symplectic before star gate

## Comparative table

| option | variant | beta | pairs | asm | pair harm | Q harm | P harm | pair J-lock | union J-lock | signed | primary symp pass | motif symp pass | best primary Ω | level | dim/rank | Ω nondeg ratio | Ω QP ratio | used dβ? |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|---|---:|---:|---:|---:|
| A_to_B_rescan_strong_reuseB | real_growth | (1,0,4,0) | 4 | 1 | 0.381795 | 0.381795 | 0.248801 | 0.448553 | 0.543826 | -0.624664 | 40 | 40 | union_plus_link_cycle_geometric_angle_identity_unit_(0, 2) | motif | 4/4 | 0.785587 | 0.785587 | False |
| A_to_B_rescan_strong_reuseB | strict_symmetrized_control | (1,0,0,0) | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |  |  | 0/0 | 0 | 0 | False |
| A_to_B_rescan_strong_reuseB | no_backreaction | (1,0,8,0) | 4 | 2 | 0.405623 | 0.405623 | 0.343367 | 0.339477 | 0.661108 | -0.285265 | 80 | 80 | union_plus_link_cycle_geometric_angle_identity_unit_(4, 7) | motif | 4/4 | 0.83832 | 0.83832 | False |

## Interpretation rule

This package separates three things that earlier tests tended to entangle:

```text
1. Q/P carrier existence,
2. real symplectic nondegeneracy of Ω on Q/P,
3. J-lock / complex orientation.
```

A positive Ω result does not mean `i` or `J` has been derived.  It means the ladder can advance to a derived real #/* search.  Important audit detail: `union_pair_exchange` already passes on the motif in the non-strict variants; the link-cycle variants are not needed to make Ω nondegenerate, even when they tie or duplicate the best score.  A negative Ω result means the obstruction is even earlier than # or J: the Q/P carrier exists, but no tested derived real symplectic form has been found on it.

## Stop/continue rule

- If a primary Ω candidate passes in non-strict variants and strict_sym is null, the next test should be `test_real_star_from_symplectic_gate.py`.
- If only `data_wedge_control` passes, the Q/P plane is independent but the symplectic form is not derived; do not count it as success.
- If no primary Ω passes, document a symplectic-level obstruction before returning to J-search.
