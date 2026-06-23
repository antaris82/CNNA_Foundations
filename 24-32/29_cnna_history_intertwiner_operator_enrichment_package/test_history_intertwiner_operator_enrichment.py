"""
CNNA / growing real complement network
History-intertwiner operator-enrichment diagnostic.

Purpose
-------
The previous generated real operator-system closure test showed that the
Record/Live L/T generator set is still too small: the finite word-space is not
stable under the DtN block adjoint.  This diagnostic checks whether the missing
adjoint directions are already present in provenance-derived history
intertwiners, or whether the visible growth window is still too small and the
*-structure should be regarded as a local limiting object.

No J is set or tested.  This is a pre-* gate.

Provenance-derived enrichment tiers
-----------------------------------
All generators live on X_ab = W_a \oplus W_b, with W_p = B_p^record \oplus
B_p^live, B_p = R^3.

base:
  I, diag(L_a,L_b), cell projectors, all visible transverse handoffs T_ab^k and
  T_ba^k.

aging:
  base plus record->live aging injections built from A_p = D_live - D_live_at_completion
  and H_p = D_live - D_record.

completion_live:
  aging plus least-squares DtN transports record<->live derived from the pair
  (D_record, D_live), with only numerical ridge used for stability.

suffix:
  completion_live plus suffix-gluing identity intertwiners W_a <-> W_b.

birth_order:
  completion_live plus directed birth-order successor transports 1->2->3 inside
  record/live blocks.  This is still a directed growth object, not a cyclic J.

all_enriched:
  all of the above.

Interpretation warning
----------------------
A residual improvement caused only by high span saturation is not evidence for
*-closure; the diagnostic reports the generated dimension as a fraction of the
full 12x12 matrix space.  A meaningful positive gate would require lower star
and multiplication residuals without trivial full-matrix saturation.
"""
from __future__ import annotations

import argparse
import csv
import json
import random
import time
import zipfile
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

import numpy as np

from record_live_block_base import (
    EPS,
    RealGrowth,
    LocalCell,
    build_cells,
    double_history_pairs,
    identical_pairs,
    random_pairs,
    block_handoff,
    block_longitudinal,
    block_metric,
    adjoint,
    fro,
    mean,
    std,
    perc,
    write_csv,
)


def orthonormal_span(mats: List[np.ndarray], tol: float) -> Tuple[np.ndarray, List[np.ndarray], List[float]]:
    good = [np.asarray(M, dtype=float).reshape(-1) for M in mats if np.linalg.norm(M) > tol]
    if not good:
        return np.zeros((0, 0), dtype=float), [], []
    A = np.stack(good, axis=1)
    U, s, _ = np.linalg.svd(A, full_matrices=False)
    if len(s) == 0:
        return np.zeros((A.shape[0], 0), dtype=float), [], []
    keep = s > tol * max(A.shape) * max(float(s[0]), 1.0)
    Q = U[:, keep]
    dim = int(round(np.sqrt(A.shape[0])))
    basis = [Q[:, i].reshape(dim, dim) for i in range(Q.shape[1])]
    return Q, basis, [float(x) for x in s]


def residual_to_Q(M: np.ndarray, Q: np.ndarray) -> Tuple[float, float]:
    v = np.asarray(M, dtype=float).reshape(-1)
    n = float(np.linalg.norm(v))
    if n < EPS:
        return 0.0, 0.0
    if Q.size == 0 or Q.shape[1] == 0:
        return n, 1.0
    proj = Q @ (Q.T @ v)
    r = v - proj
    rn = float(np.linalg.norm(r))
    return rn, rn / (n + EPS)


def embed_diag(A: np.ndarray, B: np.ndarray) -> np.ndarray:
    Z = np.zeros_like(A)
    return np.block([[A, Z], [Z, B]])


def embed_ab(T: np.ndarray) -> np.ndarray:
    Z = np.zeros_like(T)
    return np.block([[Z, Z], [T, Z]])


def embed_ba(T: np.ndarray) -> np.ndarray:
    Z = np.zeros_like(T)
    return np.block([[Z, T], [Z, Z]])


def block_pair_metric(a: LocalCell, b: LocalCell, metric_source: str, ridge: float) -> Tuple[np.ndarray, Dict[str, float]]:
    Ga, sa = block_metric(a, metric_source, ridge)
    Gb, sb = block_metric(b, metric_source, ridge)
    Z = np.zeros_like(Ga)
    return np.block([[Ga, Z], [Z, Gb]]), {
        "metric_a_cond": sa["metric_cond"],
        "metric_b_cond": sb["metric_cond"],
        "metric_a_min_eig": sa["metric_min_eig"],
        "metric_b_min_eig": sb["metric_min_eig"],
    }


def safe_pinv_metric(M: np.ndarray, ridge: float) -> np.ndarray:
    A = 0.5 * (M + M.T)
    w, V = np.linalg.eigh(A)
    scale = max(float(np.max(np.abs(w))), 1.0)
    w2 = np.maximum(w, ridge * scale)
    return (V / w2) @ V.T


def cell_aging_injection(c: LocalCell, mode: str) -> np.ndarray:
    Z = np.zeros((3, 3), dtype=float)
    if mode == "aging":
        A = c.aging.copy()
    elif mode == "handoff":
        A = c.handoff.copy()
    else:
        raise ValueError(mode)
    return np.block([[Z, Z], [A, Z]])


def cell_completion_live_transport(c: LocalCell, direction: str, ridge: float) -> np.ndarray:
    R = c.record["D_total"].copy()
    L = c.live1["D_total"].copy()
    Z = np.zeros((3, 3), dtype=float)
    if direction == "record_to_live":
        Q = L @ safe_pinv_metric(R, ridge)
        return np.block([[Z, Z], [Q, Z]])
    if direction == "live_to_record":
        Q = R @ safe_pinv_metric(L, ridge)
        return np.block([[Z, Q], [Z, Z]])
    raise ValueError(direction)


def birth_order_shift(forward: bool = True) -> np.ndarray:
    S = np.zeros((3, 3), dtype=float)
    if forward:
        S[1, 0] = 1.0
        S[2, 1] = 1.0
    else:
        S[0, 1] = 1.0
        S[1, 2] = 1.0
    return S


def cell_birth_order_transport(c: LocalCell, include_reverse: bool) -> List[np.ndarray]:
    Z = np.zeros((3, 3), dtype=float)
    S = birth_order_shift(True)
    out = [np.block([[S, Z], [Z, S]])]
    if include_reverse:
        R = birth_order_shift(False)
        out.append(np.block([[R, Z], [Z, R]]))
    return out


def suffix_identity_maps() -> Tuple[np.ndarray, np.ndarray]:
    I6 = np.eye(6)
    return embed_ab(I6), embed_ba(I6)


def build_enriched_generators(
    a: LocalCell,
    b: LocalCell,
    tier: str,
    operator_mode: str,
    longitudinal_mode: str,
    ridge: float,
    include_projectors: bool,
    include_longitudinal_powers: bool,
) -> Tuple[List[np.ndarray], Dict[str, float]]:
    La = block_longitudinal(a, longitudinal_mode)
    Lb = block_longitudinal(b, longitudinal_mode)
    I6 = np.eye(6)
    Z6 = np.zeros((6, 6), dtype=float)
    I12 = np.eye(12)
    Pa = np.block([[I6, Z6], [Z6, Z6]])
    Pb = np.block([[Z6, Z6], [Z6, I6]])
    D = embed_diag(La, Lb)
    gens: List[np.ndarray] = [I12, D]
    if include_projectors:
        gens.extend([Pa, Pb])
    if include_longitudinal_powers:
        gens.extend([D @ D, D @ D @ D])

    for r in (1, 2, 3):
        gens.append(embed_ab(block_handoff(a, b, r, operator_mode)[0]))
        gens.append(embed_ba(block_handoff(b, a, r, operator_mode)[0]))

    def add_cell_diag(ops_a: List[np.ndarray], ops_b: List[np.ndarray]) -> None:
        for A in ops_a:
            gens.append(embed_diag(A, np.zeros_like(A)))
        for B in ops_b:
            gens.append(embed_diag(np.zeros_like(B), B))
        for A, B in zip(ops_a, ops_b):
            gens.append(embed_diag(A, B))

    if tier in {"aging", "completion_live", "suffix", "birth_order", "all_enriched"}:
        add_cell_diag(
            [cell_aging_injection(a, "aging"), cell_aging_injection(a, "handoff")],
            [cell_aging_injection(b, "aging"), cell_aging_injection(b, "handoff")],
        )

    if tier in {"completion_live", "suffix", "birth_order", "all_enriched"}:
        add_cell_diag(
            [
                cell_completion_live_transport(a, "record_to_live", ridge),
                cell_completion_live_transport(a, "live_to_record", ridge),
            ],
            [
                cell_completion_live_transport(b, "record_to_live", ridge),
                cell_completion_live_transport(b, "live_to_record", ridge),
            ],
        )

    if tier in {"suffix", "all_enriched"}:
        Sab, Sba = suffix_identity_maps()
        gens.extend([Sab, Sba])

    if tier in {"birth_order", "all_enriched"}:
        # forward birth-order is growth-derived. reverse is not included by default;
        # the adjoint test will indicate whether it is forced.
        add_cell_diag(cell_birth_order_transport(a, include_reverse=False), cell_birth_order_transport(b, include_reverse=False))

    stats = {
        "seed_count": len(gens),
        "La_norm": fro(La),
        "Lb_norm": fro(Lb),
        "tier": tier,
    }
    return gens, stats


def word_space(gens: List[np.ndarray], max_degree: int, tol: float, cap: int) -> Tuple[np.ndarray, List[np.ndarray], Dict[str, float]]:
    all_mats: List[np.ndarray] = list(gens)
    prev_words: List[np.ndarray] = list(gens)
    stats: Dict[str, float] = {"degree_1_raw_words": len(gens)}
    for d in range(2, max_degree + 1):
        new_words: List[np.ndarray] = []
        for A in prev_words:
            for G in gens:
                new_words.append(A @ G)
                if len(new_words) >= cap:
                    break
            if len(new_words) >= cap:
                break
        all_mats.extend(new_words)
        Qd, basis_d, _ = orthonormal_span(all_mats, tol)
        stats[f"degree_{d}_raw_words"] = len(new_words)
        stats[f"degree_{d}_span_dim"] = len(basis_d)
        prev_words = basis_d[:cap]
        all_mats = basis_d[:cap]
    Q, basis, svals = orthonormal_span(all_mats, tol)
    stats["word_space_dim"] = len(basis)
    stats["word_space_fraction_full_144"] = len(basis) / 144.0
    stats["word_space_sval_max"] = max(svals, default=0.0)
    stats["word_space_sval_min"] = min(svals, default=0.0) if svals else 0.0
    return Q, basis, stats


def star_stats(basis: List[np.ndarray], Q: np.ndarray, G: np.ndarray) -> Dict[str, float]:
    rels = []
    for A in basis:
        _, rel = residual_to_Q(adjoint(A, G, G), Q)
        rels.append(rel)
    return {
        "star_basis_rel_mean": mean(rels),
        "star_basis_rel_p50": perc(rels, 50),
        "star_basis_rel_p95": perc(rels, 95),
        "star_basis_rel_max": max(rels, default=float("nan")),
        "star_basis_rel_lt_0p05_fraction": mean([float(x < 0.05) for x in rels]),
        "star_basis_rel_lt_0p25_fraction": mean([float(x < 0.25) for x in rels]),
    }


def seed_star_stats(gens: List[np.ndarray], Q: np.ndarray, G: np.ndarray) -> Dict[str, float]:
    rels = []
    for A in gens:
        _, rel = residual_to_Q(adjoint(A, G, G), Q)
        rels.append(rel)
    return {
        "star_seed_rel_mean": mean(rels),
        "star_seed_rel_p95": perc(rels, 95),
        "star_seed_rel_lt_0p05_fraction": mean([float(x < 0.05) for x in rels]),
        "star_seed_rel_lt_0p25_fraction": mean([float(x < 0.25) for x in rels]),
    }


def mult_stats(gens: List[np.ndarray], basis: List[np.ndarray], Q: np.ndarray, sample_cap: int, seed: int) -> Dict[str, float]:
    rng = random.Random(seed)
    pairs = [(i, j) for i in range(len(gens)) for j in range(len(basis))]
    if len(pairs) > sample_cap:
        pairs = rng.sample(pairs, sample_cap)
    rels = []
    for i, j in pairs:
        _, rel = residual_to_Q(gens[i] @ basis[j], Q)
        rels.append(rel)
    return {
        "mult_left_samples": len(rels),
        "mult_left_rel_mean": mean(rels),
        "mult_left_rel_p50": perc(rels, 50),
        "mult_left_rel_p95": perc(rels, 95),
        "mult_left_rel_lt_0p05_fraction": mean([float(x < 0.05) for x in rels]),
        "mult_left_rel_lt_0p25_fraction": mean([float(x < 0.25) for x in rels]),
    }


def row_for_pair(a: LocalCell, b: LocalCell, port: int, control: str, pair_index: int, tier: str, args: argparse.Namespace, level: int) -> dict:
    gens, gst = build_enriched_generators(
        a, b, tier, args.operator_mode, args.longitudinal_mode, args.ridge,
        args.include_projectors, args.include_longitudinal_powers,
    )
    Q, basis, wst = word_space(gens, args.degree, args.tol, args.word_cap)
    G, mst = block_pair_metric(a, b, args.metric_source, args.ridge)
    Tab = block_handoff(a, b, port, args.operator_mode)[0]
    Tba = block_handoff(b, a, port, args.operator_mode)[0]
    La = block_longitudinal(a, args.longitudinal_mode)
    Lb = block_longitudinal(b, args.longitudinal_mode)
    Cab = Lb @ Tab - Tab @ La
    return {
        "global_level": level,
        "control": control,
        "tier": tier,
        "pair_index": pair_index,
        "parent_a": a.parent,
        "parent_b": b.parent,
        "port": port,
        "cell_level": a.level,
        "age_a": a.age,
        "age_b": b.age,
        "root_sector_a": a.root_sector,
        "root_sector_b": b.root_sector,
        "suffix": ".".join(map(str, a.suffix)),
        "T_norm": fro(Tab),
        "T_reverse_norm": fro(Tba),
        "C_norm": fro(Cab),
        **gst,
        **wst,
        **mst,
        **star_stats(basis, Q, G),
        **seed_star_stats(gens, Q, G),
        **mult_stats(gens, basis, Q, args.mult_sample_cap, args.seed + 997 * level + pair_index),
    }


def truncate_pairs(pairs: List[Tuple[LocalCell, LocalCell, int]], max_pairs: int, seed: int) -> List[Tuple[LocalCell, LocalCell, int]]:
    if max_pairs <= 0 or len(pairs) <= max_pairs:
        return pairs
    rng = random.Random(seed)
    return rng.sample(pairs, max_pairs)


def summarize(rows: List[dict], label: str) -> dict:
    out: Dict[str, float | str | int] = {"label": label, "count": len(rows)}
    keys = [
        "seed_count", "word_space_dim", "word_space_fraction_full_144",
        "T_norm", "C_norm",
        "star_basis_rel_mean", "star_basis_rel_p50", "star_basis_rel_p95", "star_basis_rel_max",
        "star_basis_rel_lt_0p05_fraction", "star_basis_rel_lt_0p25_fraction",
        "star_seed_rel_mean", "star_seed_rel_p95", "star_seed_rel_lt_0p05_fraction", "star_seed_rel_lt_0p25_fraction",
        "mult_left_rel_mean", "mult_left_rel_p50", "mult_left_rel_p95", "mult_left_rel_lt_0p05_fraction", "mult_left_rel_lt_0p25_fraction",
        "metric_a_cond", "metric_b_cond", "metric_a_min_eig", "metric_b_min_eig",
    ]
    for k in keys:
        vals = [r[k] for r in rows if k in r]
        out[f"{k}_mean"] = mean(vals)
        out[f"{k}_p50"] = perc(vals, 50)
        out[f"{k}_p95"] = perc(vals, 95)
        finite = [float(v) for v in vals if np.isfinite(float(v))]
        out[f"{k}_max"] = max(finite, default=float("nan"))
    return out


def build_model(max_level: int, case: str) -> RealGrowth:
    if case == "real_growth":
        m = RealGrowth(growth_rule="sequential", br_ancestor=0.045, br_sibling=0.035, shell_normalized=True)
    elif case == "symmetrized_birth":
        m = RealGrowth(growth_rule="symmetrized_birth", br_ancestor=0.045, br_sibling=0.035, shell_normalized=True)
    elif case == "no_backreaction":
        m = RealGrowth(growth_rule="sequential", br_ancestor=0.0, br_sibling=0.0, shell_normalized=True)
    else:
        raise ValueError(case)
    m.grow(max_level)
    return m


def run(args: argparse.Namespace) -> Tuple[List[dict], List[dict]]:
    all_rows: List[dict] = []
    summaries: List[dict] = []
    for L in args.levels:
        for case in args.cases:
            model = build_model(L, case)
            cells = build_cells(model, L)
            dh = truncate_pairs(double_history_pairs(cells), args.max_pairs_per_control, args.seed + 10 * L)
            ident = truncate_pairs(identical_pairs(cells), args.max_pairs_per_control, args.seed + 20 * L)
            rand = random_pairs(cells, min(args.max_pairs_per_control, max(1, len(dh))), args.seed + 30 * L)
            pair_sets = [("double_history_suffix_quotient", dh), ("identical_history_control", ident), ("random_same_level_port_baseline", rand)]
            summaries.append({"label": f"{case}:L{L}:meta", "count": len(cells), "nodes": len(model.nodes), "completed_cells": len(cells), "double_history_pairs_used": len(dh), "identical_pairs_used": len(ident), "random_pairs_used": len(rand)})
            for tier in args.tiers:
                for control, pairs in pair_sets:
                    rows = [row_for_pair(a, b, p, f"{case}:{control}", i, tier, args, L) for i, (a, b, p) in enumerate(pairs)]
                    all_rows.extend(rows)
                    summaries.append(summarize(rows, f"{case}:L{L}:{control}:{tier}:degree{args.degree}"))
    return all_rows, summaries


def get_summary(summaries: List[dict], label: str) -> dict | None:
    for s in summaries:
        if s.get("label") == label:
            return s
    return None


def write_results(out: Path, args: argparse.Namespace, summaries: List[dict]) -> None:
    lines: List[str] = []
    lines.append("# RESULTS: history-intertwiner operator-enrichment closure")
    lines.append("")
    lines.append("## Status")
    lines.append("")
    lines.append("This is still not a J-test.  It checks whether provenance-derived history intertwiners enrich the visible L/T operator system enough to approach DtN-adjoint and finite product closure.")
    lines.append("")
    lines.append("## Parameters")
    lines.append("")
    lines.append(f"- levels: `{args.levels}`")
    lines.append(f"- tiers: `{args.tiers}`")
    lines.append(f"- degree: `{args.degree}`")
    lines.append(f"- max_pairs_per_control: `{args.max_pairs_per_control}`")
    lines.append(f"- word_cap: `{args.word_cap}`")
    lines.append("")
    lines.append("## Primary double-history scale table")
    lines.append("")
    lines.append("| L | tier | dim/full | star basis mean | star seed mean | mult mean | note |")
    lines.append("|---:|---|---:|---:|---:|---:|---|")
    for L in args.levels:
        for tier in args.tiers:
            s = get_summary(summaries, f"real_growth:L{L}:double_history_suffix_quotient:{tier}:degree{args.degree}")
            if not s:
                continue
            frac = float(s.get("word_space_fraction_full_144_mean", float("nan")))
            star = float(s.get("star_basis_rel_mean_mean", float("nan")))
            seed = float(s.get("star_seed_rel_mean_mean", float("nan")))
            mult = float(s.get("mult_left_rel_mean_mean", float("nan")))
            note = ""
            if frac > 0.90:
                note = "near full-space saturation"
            elif star < 0.10 and mult < 0.10:
                note = "strong local closure candidate"
            elif star < 0.25:
                note = "partial adjoint visibility"
            else:
                note = "not closed"
            lines.append(f"| {L} | `{tier}` | {frac:.3f} | {star:.3g} | {seed:.3g} | {mult:.3g} | {note} |")
    lines.append("")
    lines.append("## Interpretation")
    lines.append("")
    lines.append("- If enrichment only works by saturating most of the 12x12 matrix space, that is not meaningful *-closure; it is an overlarge envelope.")
    lines.append("- If residuals improve with L at comparable span dimension, that supports the idea that the *-structure is a local limiting phenomenon rather than present in the small finite window.")
    lines.append("- If residuals do not improve with L, the current provenance-derived intertwiners are still missing necessary directions.")
    lines.append("")
    out.write_text("\n".join(lines), encoding="utf-8")


def write_csv_local(path: Path, rows: List[dict]) -> None:
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    keys = sorted({k for r in rows for k in r.keys()})
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", type=Path, default=Path("history_intertwiner_operator_enrichment_out"))
    ap.add_argument("--levels", type=int, nargs="+", default=[4, 5, 6])
    ap.add_argument("--cases", nargs="+", default=["real_growth", "symmetrized_birth", "no_backreaction"])
    ap.add_argument("--tiers", nargs="+", default=["base", "aging", "completion_live", "suffix", "birth_order", "all_enriched"])
    ap.add_argument("--degree", type=int, default=3)
    ap.add_argument("--max-pairs-per-control", type=int, default=8)
    ap.add_argument("--word-cap", type=int, default=500)
    ap.add_argument("--mult-sample-cap", type=int, default=200)
    ap.add_argument("--operator-mode", default="triangular_handoff")
    ap.add_argument("--longitudinal-mode", default="triangular_record_live")
    ap.add_argument("--metric-source", default="record_live_block")
    ap.add_argument("--ridge", type=float, default=1e-9)
    ap.add_argument("--tol", type=float, default=1e-9)
    ap.add_argument("--seed", type=int, default=20260621)
    ap.add_argument("--include-projectors", action="store_true", default=True)
    ap.add_argument("--include-longitudinal-powers", action="store_true", default=True)
    args = ap.parse_args()

    t0 = time.time()
    args.out.mkdir(parents=True, exist_ok=True)
    rows, summaries = run(args)
    write_csv_local(args.out / "pair_rows_all.csv", rows)
    write_csv_local(args.out / "summary_table_all.csv", summaries)
    write_results(args.out / "RESULTS_history_intertwiner_operator_enrichment.md", args, summaries)
    summary = {
        "elapsed_seconds": time.time() - t0,
        "rows": len(rows),
        "summaries": len(summaries),
        "levels": args.levels,
        "tiers": args.tiers,
    }
    (args.out / "SUMMARY.txt").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    zip_path = args.out.with_suffix(".zip")
    with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as z:
        for p in args.out.rglob("*"):
            z.write(p, p.relative_to(args.out.parent))
    print(json.dumps({"out": str(args.out), "zip": str(zip_path), **summary}, indent=2))


if __name__ == "__main__":
    main()
