"""
CNNA / growing real complement network
Local patch operator-system closure diagnostic.

Purpose
-------
The pair-level Record/Live L/T operator tests showed that a single
Double-History handoff pair is too small to carry a stable DtN-adjoint closure.
This diagnostic enlarges the carrier from a pair

    X_ab = W_a \oplus W_b

to a growth-defined local patch

    X_P = W_{p1} \oplus ... \oplus W_{pn},
    W_p = B_p^record \oplus B_p^live, B_p = R^3.

Patch definition
----------------
The primary patch is a same-suffix / multi-root-sector Double-History patch:
all completed local cells with the same (level, suffix) but different first
root-sector.  In the ternary growth runs this is typically a 3-cell patch.
This is the minimal patch that contains several different birth histories that
become identified after forgetting the first root-sector.

Generators on X_P
-----------------
No J, no i, and no complex phase are inserted.

The generated real operator system contains:
- identity on the patch;
- block-diagonal longitudinal Record/Live DtN operators L_p;
- cell projectors;
- all growth-visible transverse Record/Live handoffs T_pq^k between all ordered
  cells in the patch and all ports k=1,2,3;
- optional suffix identity intertwiners W_p -> W_q;
- optional intra-cell Record/Live transports from completion/live DtN data.

Gate
----
The test forms finite word-spaces up to a chosen degree and asks:
1. Is the generated patch word-space stable under the DtN block adjoint?
2. Is it stable under one more left multiplication by generators?
3. Does closure improve when moving from pair-space to patch-space?

Interpretation warning
----------------------
A positive result would still be pre-J and pre-C*-algebra.  It would only say
that a local real operator-system candidate is emerging on patches.  A negative
result means that either the patch is still too small, the generator family is
wrong/too small, or *-algebra is an infinite/local limiting object.
"""
from __future__ import annotations

import argparse
import csv
import json
import random
import time
import zipfile
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple

import numpy as np

from record_live_block_base import (
    EPS,
    RealGrowth,
    LocalCell,
    build_cells,
    block_handoff,
    block_longitudinal,
    block_metric,
    adjoint,
    fro,
    mean,
    perc,
    write_csv,
)


def std(xs: Iterable[float]) -> float:
    vals = [float(x) for x in xs if np.isfinite(float(x))]
    return float(np.std(vals)) if vals else float("nan")


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


def write_csv_local(path: Path, rows: List[dict]) -> None:
    write_csv(path, rows)


def embed_block(op: np.ndarray, i: int, j: int, n: int, wdim: int = 6) -> np.ndarray:
    M = np.zeros((n * wdim, n * wdim), dtype=float)
    M[i * wdim:(i + 1) * wdim, j * wdim:(j + 1) * wdim] = op
    return M


def block_diag_metric(cells: Sequence[LocalCell], metric_source: str, ridge: float) -> Tuple[np.ndarray, Dict[str, float]]:
    blocks = []
    conds = []
    mins = []
    for c in cells:
        G, st = block_metric(c, metric_source, ridge)
        blocks.append(G)
        conds.append(st["metric_cond"])
        mins.append(st["metric_min_eig"])
    Z = np.zeros_like(blocks[0])
    rows = []
    for i, A in enumerate(blocks):
        rows.append([A if i == j else Z.copy() for j in range(len(blocks))])
    Gp = np.block(rows)
    eig = np.linalg.eigvalsh(0.5 * (Gp + Gp.T))
    return Gp, {
        "patch_metric_cond_mean": mean(conds),
        "patch_metric_cond_max": max(conds, default=float("nan")),
        "patch_metric_min_eig_mean": mean(mins),
        "patch_metric_global_min_eig": float(eig[0]),
        "patch_metric_global_cond": float(eig[-1] / max(eig[0], EPS)),
    }


def safe_pinv_metric(M: np.ndarray, ridge: float) -> np.ndarray:
    A = 0.5 * (M + M.T)
    w, V = np.linalg.eigh(A)
    scale = max(float(np.max(np.abs(w))), 1.0)
    w2 = np.maximum(w, ridge * scale)
    return (V / w2) @ V.T


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


def cell_aging_injection(c: LocalCell, mode: str) -> np.ndarray:
    Z = np.zeros((3, 3), dtype=float)
    A = c.aging.copy() if mode == "aging" else c.handoff.copy()
    return np.block([[Z, Z], [A, Z]])


def patch_generators(
    cells: Sequence[LocalCell],
    tier: str,
    operator_mode: str,
    longitudinal_mode: str,
    metric_source: str,
    ridge: float,
    include_projectors: bool,
    include_diag_powers: bool,
) -> Tuple[List[np.ndarray], Dict[str, float]]:
    n = len(cells)
    dim = 6 * n
    I = np.eye(dim)
    gens: List[np.ndarray] = [I]

    Ls = [block_longitudinal(c, longitudinal_mode) for c in cells]
    D = sum(embed_block(Ls[i], i, i, n) for i in range(n))
    gens.append(D)
    if include_diag_powers:
        gens.extend([D @ D, D @ D @ D])

    if include_projectors:
        for i in range(n):
            gens.append(embed_block(np.eye(6), i, i, n))

    # Local longitudinal generators separately, not only as the total diagonal D.
    for i, L in enumerate(Ls):
        gens.append(embed_block(L, i, i, n))

    # Transverse handoffs between all ordered distinct cells and all ports.
    T_norms = []
    for i, a in enumerate(cells):
        for j, b in enumerate(cells):
            if i == j:
                continue
            for port in (1, 2, 3):
                T, _ = block_handoff(a, b, port, operator_mode)
                T_norms.append(fro(T))
                gens.append(embed_block(T, j, i, n))

    if tier in {"suffix_identity", "completion_live", "full_patch"}:
        I6 = np.eye(6)
        for i in range(n):
            for j in range(n):
                if i != j:
                    gens.append(embed_block(I6, j, i, n))

    if tier in {"completion_live", "full_patch"}:
        for i, c in enumerate(cells):
            for A in [
                cell_aging_injection(c, "aging"),
                cell_aging_injection(c, "handoff"),
                cell_completion_live_transport(c, "record_to_live", ridge),
                cell_completion_live_transport(c, "live_to_record", ridge),
            ]:
                gens.append(embed_block(A, i, i, n))

    if tier == "full_patch":
        # Derived directed birth-order shift inside each local record/live block.
        S = np.zeros((3, 3), dtype=float)
        S[1, 0] = 1.0
        S[2, 1] = 1.0
        Z = np.zeros((3, 3), dtype=float)
        B = np.block([[S, Z], [Z, S]])
        for i in range(n):
            gens.append(embed_block(B, i, i, n))

    return gens, {
        "seed_count": len(gens),
        "patch_cells": n,
        "patch_dim": dim,
        "handoff_T_norm_mean": mean(T_norms),
        "handoff_T_norm_p95": perc(T_norms, 95),
        "diag_L_norm_mean": mean([fro(L) for L in Ls]),
    }


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
    dim = gens[0].shape[0]
    stats["word_space_dim"] = len(basis)
    stats["word_space_fraction_full"] = len(basis) / float(dim * dim)
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


def commutator_stats(cells: Sequence[LocalCell], operator_mode: str, longitudinal_mode: str) -> Dict[str, float]:
    n = len(cells)
    Ls = [block_longitudinal(c, longitudinal_mode) for c in cells]
    vals = []
    rels = []
    for i, a in enumerate(cells):
        for j, b in enumerate(cells):
            if i == j:
                continue
            for port in (1, 2, 3):
                T, _ = block_handoff(a, b, port, operator_mode)
                C = Ls[j] @ T - T @ Ls[i]
                vals.append(fro(C))
                rels.append(fro(C) / (fro(Ls[j] @ T) + fro(T @ Ls[i]) + EPS))
    return {
        "patch_comm_C_norm_mean": mean(vals),
        "patch_comm_C_norm_p95": perc(vals, 95),
        "patch_comm_C_rel_mean": mean(rels),
        "patch_comm_C_rel_p95": perc(rels, 95),
    }


def patch_key(cells: Sequence[LocalCell]) -> str:
    return ";".join(".".join(map(str, c.address)) for c in cells)


def same_suffix_patches(cells: Dict[int, LocalCell], patch_size: int) -> List[List[LocalCell]]:
    groups: Dict[Tuple[int, Tuple[int, ...]], List[LocalCell]] = defaultdict(list)
    for c in cells.values():
        if c.level < 1:
            continue
        groups[(c.level, c.suffix)].append(c)
    patches = []
    for _, gs in sorted(groups.items(), key=lambda kv: (kv[0][0], kv[0][1])):
        roots = {g.root_sector for g in gs}
        if len(gs) >= patch_size and len(roots) >= min(2, patch_size):
            gs2 = sorted(gs, key=lambda c: (c.root_sector, c.parent))[:patch_size]
            patches.append(gs2)
    return patches


def identical_patches(cells: Dict[int, LocalCell], patch_size: int) -> List[List[LocalCell]]:
    out = []
    for c in sorted(cells.values(), key=lambda x: (x.level, x.parent)):
        if c.level >= 1:
            out.append([c for _ in range(patch_size)])
    return out


def random_patches(cells: Dict[int, LocalCell], patch_size: int, n: int, seed: int) -> List[List[LocalCell]]:
    rng = random.Random(seed)
    buckets: Dict[int, List[LocalCell]] = defaultdict(list)
    for c in cells.values():
        if c.level >= 1:
            buckets[c.level].append(c)
    levels = [k for k, v in buckets.items() if len(v) >= patch_size]
    out: List[List[LocalCell]] = []
    attempts = 0
    while len(out) < n and attempts < max(1000, 50 * n):
        attempts += 1
        L = rng.choice(levels)
        gs = rng.sample(buckets[L], patch_size)
        # Prefer not same suffix: a genuine nonquotient/random control.
        if len({g.suffix for g in gs}) == 1:
            continue
        out.append(sorted(gs, key=lambda c: c.parent))
    return out


def truncate_patches(patches: List[List[LocalCell]], max_patches: int, seed: int) -> List[List[LocalCell]]:
    if max_patches <= 0 or len(patches) <= max_patches:
        return patches
    rng = random.Random(seed)
    return rng.sample(patches, max_patches)


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


def row_for_patch(patch: Sequence[LocalCell], control: str, patch_index: int, tier: str, args: argparse.Namespace, level: int) -> dict:
    gens, gst = patch_generators(
        patch,
        tier,
        args.operator_mode,
        args.longitudinal_mode,
        args.metric_source,
        args.ridge,
        args.include_projectors,
        args.include_diag_powers,
    )
    Q, basis, wst = word_space(gens, args.degree, args.tol, args.word_cap)
    G, mst = block_diag_metric(patch, args.metric_source, args.ridge)
    return {
        "global_level": level,
        "control": control,
        "tier": tier,
        "patch_index": patch_index,
        "patch_key": patch_key(patch),
        "patch_cell_level_min": min(c.level for c in patch),
        "patch_cell_level_max": max(c.level for c in patch),
        "patch_age_mean": mean([c.age for c in patch]),
        "patch_root_sectors": ".".join(map(str, sorted({c.root_sector for c in patch}))),
        "patch_suffixes": "|".join(sorted({".".join(map(str, c.suffix)) for c in patch})),
        **gst,
        **wst,
        **mst,
        **commutator_stats(patch, args.operator_mode, args.longitudinal_mode),
        **star_stats(basis, Q, G),
        **seed_star_stats(gens, Q, G),
        **mult_stats(gens, basis, Q, args.mult_sample_cap, args.seed + 997 * level + patch_index),
    }


def summarize(rows: List[dict], label: str) -> dict:
    out: Dict[str, float | str | int] = {"label": label, "count": len(rows)}
    keys = [
        "patch_cells", "patch_dim", "seed_count", "word_space_dim", "word_space_fraction_full",
        "handoff_T_norm_mean", "patch_comm_C_norm_mean", "patch_comm_C_rel_mean",
        "star_basis_rel_mean", "star_basis_rel_p50", "star_basis_rel_p95", "star_basis_rel_max",
        "star_basis_rel_lt_0p05_fraction", "star_basis_rel_lt_0p25_fraction",
        "star_seed_rel_mean", "star_seed_rel_p95", "star_seed_rel_lt_0p05_fraction", "star_seed_rel_lt_0p25_fraction",
        "mult_left_rel_mean", "mult_left_rel_p50", "mult_left_rel_p95", "mult_left_rel_lt_0p05_fraction", "mult_left_rel_lt_0p25_fraction",
        "patch_metric_global_cond", "patch_metric_global_min_eig",
    ]
    for k in keys:
        vals = [r[k] for r in rows if k in r]
        out[f"{k}_mean"] = mean(vals)
        out[f"{k}_p50"] = perc(vals, 50)
        out[f"{k}_p95"] = perc(vals, 95)
        finite = [float(v) for v in vals if np.isfinite(float(v))]
        out[f"{k}_max"] = max(finite, default=float("nan"))
    return out


def run(args: argparse.Namespace) -> Tuple[List[dict], List[dict]]:
    all_rows: List[dict] = []
    summaries: List[dict] = []
    for L in args.levels:
        for case in args.cases:
            model = build_model(L, case)
            cells = build_cells(model, L)
            same = truncate_patches(same_suffix_patches(cells, args.patch_size), args.max_patches_per_control, args.seed + 10 * L)
            ident = truncate_patches(identical_patches(cells, args.patch_size), args.max_patches_per_control, args.seed + 20 * L)
            rand = random_patches(cells, args.patch_size, len(same), args.seed + 30 * L)
            patch_sets = [
                ("same_suffix_multi_history_patch", same),
                ("identical_history_patch_control", ident),
                ("random_same_level_patch_baseline", rand),
            ]
            summaries.append({
                "label": f"{case}:L{L}:meta",
                "count": len(cells),
                "nodes": len(model.nodes),
                "completed_cells": len(cells),
                "same_suffix_patches_used": len(same),
                "identical_patches_used": len(ident),
                "random_patches_used": len(rand),
            })
            for tier in args.tiers:
                for control, patches in patch_sets:
                    rows = [row_for_patch(p, f"{case}:{control}", i, tier, args, L) for i, p in enumerate(patches)]
                    all_rows.extend(rows)
                    summaries.append(summarize(rows, f"{case}:L{L}:{control}:{tier}:patch{args.patch_size}:degree{args.degree}"))
    return all_rows, summaries


def get_summary(summaries: List[dict], label: str) -> dict | None:
    for s in summaries:
        if s.get("label") == label:
            return s
    return None


def write_results(path: Path, args: argparse.Namespace, summaries: List[dict]) -> None:
    lines: List[str] = []
    lines.append("# RESULTS: local patch operator-system closure")
    lines.append("")
    lines.append("## Purpose")
    lines.append("")
    lines.append("This diagnostic enlarges the carrier from a single Double-History pair to a same-suffix local patch of Record/Live boundary cells. It asks whether DtN-adjoint and finite product closure become visible only at patch level.")
    lines.append("")
    lines.append("## Parameters")
    lines.append("")
    lines.append("```json")
    lines.append(json.dumps({
        "levels": args.levels,
        "cases": args.cases,
        "tiers": args.tiers,
        "patch_size": args.patch_size,
        "degree": args.degree,
        "word_cap": args.word_cap,
        "max_patches_per_control": args.max_patches_per_control,
        "operator_mode": args.operator_mode,
        "longitudinal_mode": args.longitudinal_mode,
        "metric_source": args.metric_source,
    }, indent=2))
    lines.append("```")
    lines.append("")
    lines.append("## Primary summaries")
    lines.append("")
    for case in args.cases:
        for L in args.levels:
            for tier in args.tiers:
                label = f"{case}:L{L}:same_suffix_multi_history_patch:{tier}:patch{args.patch_size}:degree{args.degree}"
                s = get_summary(summaries, label)
                if not s:
                    continue
                lines.append(f"### {label}")
                lines.append("")
                keys = [
                    "count",
                    "patch_dim_mean",
                    "seed_count_mean",
                    "word_space_dim_mean",
                    "word_space_fraction_full_mean",
                    "handoff_T_norm_mean_mean",
                    "patch_comm_C_norm_mean_mean",
                    "patch_comm_C_rel_mean_mean",
                    "star_basis_rel_mean_mean",
                    "star_basis_rel_p95_mean",
                    "star_basis_rel_lt_0p25_fraction_mean",
                    "star_seed_rel_mean_mean",
                    "star_seed_rel_lt_0p25_fraction_mean",
                    "mult_left_rel_mean_mean",
                    "mult_left_rel_lt_0p25_fraction_mean",
                    "patch_metric_global_cond_mean",
                ]
                lines.append("```text")
                for k in keys:
                    lines.append(f"{k:42s} {s.get(k)}")
                lines.append("```")
                lines.append("")
    lines.append("## Interpretation")
    lines.append("")
    lines.append("A meaningful positive pre-* gate would require star residuals and multiplication residuals to decrease at patch level without saturating the full matrix algebra. If residuals remain large, the current patch is still too small, the generated operators are still incomplete, or *-closure is a limiting/local-net phenomenon rather than a finite-patch property.")
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--levels", nargs="+", type=int, default=[4])
    ap.add_argument("--cases", nargs="+", default=["real_growth"])
    ap.add_argument("--tiers", nargs="+", default=["base_patch"])
    ap.add_argument("--patch-size", type=int, default=3)
    ap.add_argument("--max-patches-per-control", type=int, default=2)
    ap.add_argument("--degree", type=int, default=2)
    ap.add_argument("--word-cap", type=int, default=180)
    ap.add_argument("--mult-sample-cap", type=int, default=100)
    ap.add_argument("--operator-mode", default="triangular_handoff")
    ap.add_argument("--longitudinal-mode", default="triangular_record_live")
    ap.add_argument("--metric-source", default="record_live_block")
    ap.add_argument("--ridge", type=float, default=1e-9)
    ap.add_argument("--tol", type=float, default=1e-9)
    ap.add_argument("--seed", type=int, default=20260621)
    ap.add_argument("--out", type=Path, default=Path("local_patch_operator_system_out_L4"))
    ap.add_argument("--include-projectors", action="store_true", default=True)
    ap.add_argument("--include-diag-powers", action="store_true", default=True)
    args = ap.parse_args()

    t0 = time.time()
    args.out.mkdir(parents=True, exist_ok=True)
    rows, summaries = run(args)
    write_csv_local(args.out / "patch_rows_all.csv", rows)
    write_csv_local(args.out / "summary_table_all.csv", summaries)
    write_results(args.out / "RESULTS_local_patch_operator_system_closure.md", args, summaries)
    summary = {
        "elapsed_seconds": time.time() - t0,
        "rows": len(rows),
        "summaries": len(summaries),
        "levels": args.levels,
        "tiers": args.tiers,
        "patch_size": args.patch_size,
    }
    (args.out / "SUMMARY.txt").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    zip_path = args.out.with_suffix(".zip")
    with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as z:
        for p in args.out.rglob("*"):
            z.write(p, p.relative_to(args.out.parent))
    print(json.dumps({"out": str(args.out), "zip": str(zip_path), **summary}, indent=2))


if __name__ == "__main__":
    main()
