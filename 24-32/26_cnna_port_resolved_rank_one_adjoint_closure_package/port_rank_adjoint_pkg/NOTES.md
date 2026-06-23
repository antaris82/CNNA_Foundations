# Notes

Interpretation of the primary L6 run:

- Keeping the full 3-port boundary space greatly improves the metric conditioning compared with the earlier compressed Env/UV 2D role-plane test.
- The metric condition number drops from an effectively degenerate role-plane scale to roughly order 10 in the primary run.
- However, the adjoint of the real-growth double-history transverse map does not close in the simple reverse all-port handoff span.
- The broad rank-one envelope closes almost trivially, but this envelope is too large to count as a generated operator system.
- The no-backreaction control closes much better in the all-port reverse span, which suggests that the live aging/backreaction layer introduces missing operator channels not captured by the simple rank-one reverse handoffs.

Methodological consequence:

The next test should keep a larger record/live block operator space rather than compressing record and live into one port map.  The likely missing structure is an adjoint pairing between record and live sectors, not merely within a single live 3-port boundary space.
