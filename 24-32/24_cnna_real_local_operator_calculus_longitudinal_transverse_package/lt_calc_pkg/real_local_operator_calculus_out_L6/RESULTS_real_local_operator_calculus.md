# RESULTS: real local operator calculus, longitudinal/transverse

## Status

This is a real operator-calculus diagnostic. It does not derive J and does not test J²=-I.
It checks whether true event-resolved growth produces composable longitudinal and transverse operators.

## Parameters

- max_level: `6`
- conductance mode: `linear`
- source modes: `handoff, aging, live, record`
- random seed: `20260621`

## Primary handoff gate

- double-history handoff pairs: `1089`
- mean transverse T norm: `0.00416198550219`
- mean composition defect ||L_q T - T L_p||: `0.0253828655862`
- mean relative composition defect: `0.571854029657`
- nonzero T fraction: `1`
- nonzero composition-defect fraction: `1`

## Controls

### identical-history
- rows: `1089`
- mean T norm: `0`
- mean composition defect: `0`
- mean relative composition defect: `0`

### real-growth random same-level/port baseline
- rows: `1089`
- mean T norm: `0.0203288313541`
- mean composition defect: `0.123439475087`
- mean relative composition defect: `0.593605825498`

### no-backreaction double-history
- rows: `1089`
- mean T norm: `0.00384199061001`
- mean composition defect: `0.0183675071613`
- mean relative composition defect: `0.200979272915`

### symmetrized-birth double-history
- rows: `1089`
- mean T norm: `0.000150080459731`
- mean composition defect: `0.000600018208061`
- mean relative composition defect: `0.556812130492`

## Interpretation

Positive result criterion for this gate: double-history handoff rows have nonzero transverse T and nonzero mixed composition defect, while identical-history controls vanish. This would mean that Script 23's DtN defect can be promoted to a real longitudinal/transverse operator-calculus candidate.

Negative result criterion: if T or the composition defect vanishes under real growth, then Script 23 remains a scalar/vector handoff obstruction, not yet an operator-calculus object.

Next gate after a positive result: define a canonical pairing/energy form and test whether the generated real operator system admits a stable adjoint operation A -> A*.
