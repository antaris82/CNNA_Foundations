# RESULTS: record/live block adjoint closure

## Status

This is a Record/Live block adjoint-closure diagnostic. It does not derive J and does not test J²=-I.
It asks whether the adjoint directions missing in the live-only 3-port test become visible after typing each local boundary space as record ⊕ live.

## Parameters

- max_level: `5`
- conductance mode: `linear`
- operator modes: `triangular_handoff`
- longitudinal modes: `triangular_record_live`
- metric sources: `record_live_block`
- ridge: `1e-10`
- random seed: `20260621`

## Primary gate

- primary label: `real_growth:double_history_suffix_quotient:triangular_handoff:triangular_record_live:record_live_block`
- rows: `360`
- mean T norm: `0.0902866087177`
- mean mixed composition norm: `0.144199426456`
- mean mixed composition relative norm: `0.282858268684`
- mean T* same-port reverse residual: `0.974577406765`
- mean T* all-ports reverse residual: `0.974134084785`
- mean T* block-envelope residual: `8.23961797422e-15`
- mean C* same-port reverse residual: `0.984303208337`
- mean C* all-ports reverse residual: `0.965574549992`
- mean C* block-envelope residual: `2.99532201515e-15`
- fraction T* all-ports residual < 0.25: `0`
- fraction T* block-envelope residual < 0.25: `1`
- fraction C* all-ports residual < 0.25: `0`
- fraction C* block-envelope residual < 0.25: `1`
- mean metric condition number: `13.0340716765`
- mean metric min eigenvalue: `0.26756308083`

## Controls

### identical-history
- rows: `360`
- mean T norm: `0`
- mean C norm: `0`
- mean T* all-ports reverse residual: `0`
- mean T* block-envelope residual: `0`
- mean C* all-ports reverse residual: `0`
- mean C* block-envelope residual: `0`
- mean metric condition number: `12.8569839937`

### symmetrized-birth
- rows: `360`
- mean T norm: `0.00350710683491`
- mean C norm: `0.00386684545201`
- mean T* all-ports reverse residual: `0.966141703173`
- mean T* block-envelope residual: `3.48275037596e-15`
- mean C* all-ports reverse residual: `0.934011857225`
- mean C* block-envelope residual: `3.54758395239e-15`
- mean metric condition number: `11.818815799`

### no-backreaction
- rows: `360`
- mean T norm: `0.0816834774379`
- mean C norm: `0.107971040023`
- mean T* all-ports reverse residual: `0.974676115398`
- mean T* block-envelope residual: `6.40630279681e-15`
- mean C* all-ports reverse residual: `0.94836157599`
- mean C* block-envelope residual: `3.85605871219e-15`
- mean metric condition number: `13.9082143817`

### random same-level/port
- rows: `360`
- mean T norm: `0.203017302094`
- mean C norm: `0.322600600632`
- mean T* all-ports reverse residual: `0.972689286312`
- mean T* block-envelope residual: `7.76777364731e-15`
- mean C* all-ports reverse residual: `0.973129736463`
- mean C* block-envelope residual: `3.07641895367e-15`
- mean metric condition number: `14.3347548998`

## Interpretation

The test is positive only if the block-typed reverse operator system closes significantly better than the live-only reverse system. A small block-envelope residual but a large all-ports generated residual means that the adjoint lives in the typed block rank-one universe, but not yet in the actually growth-generated reverse handoff family.

If closure improves, the next gate is a finite generated real operator-system closure under addition, multiplication, and the block adjoint. If closure remains weak, the boundary history space must be enriched beyond record/live blocks before attempting any J extraction.
