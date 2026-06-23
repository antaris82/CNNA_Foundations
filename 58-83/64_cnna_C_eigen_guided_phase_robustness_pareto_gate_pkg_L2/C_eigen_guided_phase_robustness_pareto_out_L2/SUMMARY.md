# SUMMARY — C-eigen guided phase robustness Pareto gate

Model label:
CNNA growing primal simplicial complex with deterministic sequential provenance growth, A-gated nonlinear complement-pair candidate space, directed antisymmetric birth-transport operators, and local C/J pair algebra.

This package does **not** apply a new growth rule. It audits the already enumerated A-gated candidate space from the C-eigen-guided pairing package and matches the same candidates across phase_sign +1/-1.

Anti-smuggling rule:
`delta_beta2` is included only as an audit column. It is not used in the Pareto score and not used as a decision input.

## Comparative summary

| variant/rule | matched | A-gated | beta2 audit | C-pass | flip-pass | all-pass | best score | best C-lock | best flip |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| no_backreaction / A_rank_rule | 9 | 9 | 9 | 2 | 0 | 0 | 1.30904e-12 | 0.154455 | 1 |
| no_backreaction / C_eigen_guided_rule | 9 | 9 | 9 | 2 | 0 | 0 | 1.30904e-12 | 0.154455 | 1 |
| real_growth / A_rank_rule | 9 | 9 | 9 | 1 | 0 | 0 | 8.45042e-13 | 0.173092 | 1 |
| real_growth / C_eigen_guided_rule | 9 | 9 | 9 | 1 | 0 | 0 | 8.45042e-13 | 0.173092 | 1 |

## Top candidates by Pareto audit score

| variant/rule | cand | C-lock worst | flip abs | signed amp min | Q/P proxy | A-rank avg | directed avg | transv avg | beta2 audit | selected? |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| no_backreaction / A_rank_rule | 32 | 0.154455 | 1 | 0.260742 | 0.905946 | 23.951 | 0.780578 | 0.336121 | True | False |
| no_backreaction / A_rank_rule | 36 | 0.44351 | 1 | 0.547273 | 0.833675 | 23.4294 | 0.590267 | 0.34173 | True | False |
| no_backreaction / A_rank_rule | 46 | 0.222113 | 1 | 0.275504 | 0.738668 | 22.5835 | 0.566622 | 0.328772 | True | False |
| no_backreaction / A_rank_rule | 44 | 0.238211 | 1 | 0.0560382 | 0.919146 | 22.8267 | 0.604231 | 0.315936 | True | True |
| no_backreaction / A_rank_rule | 34 | 0.181179 | 1 | 0.0139046 | 0.71899 | 23.549 | 0.590267 | 0.419414 | True | True |
| no_backreaction / A_rank_rule | 38 | 0.808002 | 1 | 0.52183 | 0.443592 | 18.8569 | 0.122062 | 0.265502 | True | False |
| no_backreaction / A_rank_rule | 52 | 0.808002 | 1 | 0.52183 | 0.443592 | 18.8569 | 0.122062 | 0.265502 | True | False |
| no_backreaction / A_rank_rule | 42 | 0.710251 | 1 | 0.333564 | 0.372949 | 20.1077 | 0.19557 | 0.281228 | True | False |
| no_backreaction / A_rank_rule | 40 | 0.908275 | 1 | 0.6577 | 0.212454 | 18.5595 | 0.0935159 | 0.296181 | True | False |
| no_backreaction / C_eigen_guided_rule | 32 | 0.154455 | 1 | 0.260742 | 0.905946 | 23.951 | 0.780578 | 0.336121 | True | True |
| no_backreaction / C_eigen_guided_rule | 36 | 0.44351 | 1 | 0.547273 | 0.833675 | 23.4294 | 0.590267 | 0.34173 | True | False |
| no_backreaction / C_eigen_guided_rule | 34 | 0.181179 | 1 | 0.0139046 | 0.71899 | 23.549 | 0.590267 | 0.419414 | True | False |
| no_backreaction / C_eigen_guided_rule | 38 | 0.808002 | 1 | 0.52183 | 0.443592 | 18.8569 | 0.122062 | 0.265502 | True | False |
| no_backreaction / C_eigen_guided_rule | 62 | 0.808002 | 1 | 0.52183 | 0.443592 | 18.8569 | 0.122062 | 0.265502 | True | True |
| no_backreaction / C_eigen_guided_rule | 42 | 0.710251 | 1 | 0.333564 | 0.372949 | 20.1077 | 0.19557 | 0.281228 | True | False |
| no_backreaction / C_eigen_guided_rule | 66 | 0.710251 | 1 | 0.333564 | 0.372949 | 20.1077 | 0.19557 | 0.281228 | True | True |
| no_backreaction / C_eigen_guided_rule | 40 | 0.908275 | 1 | 0.6577 | 0.212454 | 18.5595 | 0.0935159 | 0.296181 | True | False |
| no_backreaction / C_eigen_guided_rule | 64 | 0.908275 | 1 | 0.6577 | 0.212454 | 18.5595 | 0.0935159 | 0.296181 | True | False |
| real_growth / A_rank_rule | 32 | 0.173092 | 1 | 0.184157 | 0.929409 | 19.5051 | 0.479784 | 0.336121 | True | False |
| real_growth / A_rank_rule | 36 | 0.567192 | 1 | 0.46862 | 0.848112 | 19.3645 | 0.38115 | 0.34173 | True | False |
