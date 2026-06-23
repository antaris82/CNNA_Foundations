# Notes

## Why this test exists

Tests 33-36 showed that parent-fan 2-simplices carry a local noncommutative DtN plaquette sector. They also showed that the local J-planes glue signlessly and almost flatly under minimal SO(3) axis transport. The untested combination was:

```text
parent-fan plaquette J-sector + double-history suffix identification
```

This package tests that combination.

## What counts as success

A strong positive result would be:

```text
real_growth / directed double-history quotient:
  nonzero K amplitude
  nonzero sign/holonomy obstruction
  identity control vanishes
  diagonal/trace controls vanish
  symmetrized_birth collapses
  random baseline is weaker/incoherent
```

## Current caveat

The primary L5/L6 runs show a persistent directed-cyclic `Jfixed ≈ 2/3` signature, but no nontrivial Z2 sign-lift obstruction (`z2obs = 0`). The same-label quotient glues much more fixed-sign coherently than the cyclic root-shift. Therefore this package supports the presence of a local Z3 label/gluing signature, but does not yet prove Stufe-4 locking.
