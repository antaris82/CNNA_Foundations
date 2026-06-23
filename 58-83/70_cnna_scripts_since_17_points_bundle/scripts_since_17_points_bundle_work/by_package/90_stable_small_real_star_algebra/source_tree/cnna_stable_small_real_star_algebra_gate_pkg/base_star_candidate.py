#!/usr/bin/env python3
"""
CNNA live-semigroup operator-family # candidate gate with depth scaling.

Purpose
-------
Continue the established script-1/script-2 true-Schur/DtN growth line.
Use only bridge/passivity-positive Record->Live rows and ask whether the
post-birth live relaxation sequence on a fixed boundary cut generates a small
real operator family stable under the record-DtN metric adjoint candidate #.

This is deliberately not a J/i/QP/spin test and not a C*-algebra claim.

Important strengthening from the previous package
-------------------------------------------------
The test also audits the user's scaling point: irreversibility should not be a
constant property of one birth event.  It should accumulate with depth because
birth effects are distributed into old conductances and then mixed with live
relaxation.  Therefore record/live gaps and operator-family signals are grouped
by parent_level and run for L=2,3,4 finite approximants.

No J, i, Hodge, physical Hilbert adjoint, C*-norm, positivity axiom, Q/P target,
or delta-beta decision is used.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import shutil
import zipfile
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np

from base_dynamic_boundary_role_schur_dtn import DynamicSchurDtNRoleModel, DtNResult, EPS, write_csv
from base_schur_dtn_birth_relaxation import SchurDtNRelaxationModel
from base_live_semigroup_polarity_to_operator_bridge import (
    zero_sum_basis,
    project_zero_sum,
    longitudinal_mode,
    rel_fro,
    matrix_metrics,
    mean,
    median,
    sym,
)


def _safe_norm(A: np.ndarray) -> float:
    return float(np.linalg.norm(A, ord="fro")) if A.size else 0.0


def span_project_residual(M: np.ndarray, basis: List[np.ndarray]) -> float:
    if not basis:
        return 1.0 if _safe_norm(M) > EPS else 0.0
    vecs = [B.reshape(-1) for B in basis if B.shape == M.shape]
    if not vecs:
        return 1.0 if _safe_norm(M) > EPS else 0.0
    Bmat = np.stack(vecs, axis=1)
    y = M.reshape(-1)
    coef, *_ = np.linalg.lstsq(Bmat, y, rcond=None)
    pred = Bmat @ coef
    return float(np.linalg.norm(y - pred) / (np.linalg.norm(y) + EPS))


def numerical_rank_of_family(basis: List[np.ndarray], tol: float = 1e-9) -> int:
    if not basis:
        return 0
    vecs = [B.reshape(-1) for B in basis]
    M = np.stack(vecs, axis=1)
    s = np.linalg.svd(M, compute_uv=False)
    if s.size == 0:
        return 0
    return int(np.sum(s > tol * max(1.0, float(s[0]))))


def metric_adjoint(A: np.ndarray, G: np.ndarray) -> np.ndarray:
    # Record-DtN metric adjoint candidate.  This exists linearly; it is not a
    # success criterion by itself.  G is derived from the record Schur/DtN metric.
    reg = max(1e-10, 1e-8 * max(float(np.linalg.norm(G, 2)) if G.size else 1.0, 1.0))
    Gp = np.linalg.pinv(G + reg * np.eye(G.shape[0]), rcond=1e-10)
    return Gp @ A.T @ G


def operator_sequence_family_metrics(record: DtNResult, seq: List[DtNResult], cut: int, axis_sign: int) -> dict:
    """Analyze a fixed-boundary live sequence as a candidate real operator family.

    For k=1..N, A_k = G^+ (Lambda_k - Lambda_record) on the zero-sum boundary
    subspace.  The candidate # is the G-adjoint.  The gate asks whether the small
    family span{I,A_k} is stable under #, products and commutators.  Since #
    exists by construction, isolated involution/anti-multiplicativity is logged
    but not counted by itself.
    """
    mm = matrix_metrics(record, seq, cut, axis_sign)
    if not mm.get("valid"):
        return {"valid": False, "reason": mm.get("reason", "invalid_matrix_metrics")}
    boundary = list(record.boundary)
    n = len(boundary)
    if n <= 1:
        return {"valid": False, "reason": "boundary_dim_le_1"}
    G = project_zero_sum(sym(np.asarray(record.lam, dtype=float)))
    d = G.shape[0]
    if d <= 0:
        return {"valid": False, "reason": "zero_projected_dim"}
    evals = np.linalg.eigvalsh(sym(G))
    g_min = float(np.min(evals)) if evals.size else 0.0
    g_max = float(np.max(evals)) if evals.size else 0.0
    if g_min < -1e-8:
        return {"valid": False, "reason": "record_metric_not_psd"}
    reg = max(1e-10, 1e-8 * max(g_max, 1.0))
    Gp = np.linalg.pinv(G + reg * np.eye(d), rcond=1e-10)

    As: List[np.ndarray] = []
    inc_norms = []
    prevP = G.copy()
    for item in seq:
        P = project_zero_sum(sym(np.asarray(item.lam, dtype=float)))
        if P.shape != G.shape:
            continue
        Delta = sym(P - G)
        A = Gp @ Delta
        As.append(A)
        inc_norms.append(float(np.linalg.norm(P - prevP, ord="fro")))
        prevP = P
    if not As:
        return {"valid": False, "reason": "no_sequence_operators"}

    I = np.eye(d)
    basis_raw = [I] + As
    # Normalize basis for numerical projection, but keep products raw.
    basis = []
    for B in basis_raw:
        nb = _safe_norm(B)
        if nb > EPS:
            basis.append(B / nb)
    fam_rank = numerical_rank_of_family(basis)

    # # closure and involution.
    star_res = []
    invol_res = []
    anti_mult_res = []
    product_res = []
    comm_res = []
    comm_norms = []
    sample = As[: min(len(As), 3)]  # keep runtime deterministic and bounded
    for A in sample:
        Ash = metric_adjoint(A, G)
        star_res.append(span_project_residual(Ash, basis))
        Ash2 = metric_adjoint(Ash, G)
        invol_res.append(rel_fro(Ash2, A))
    for A in sample:
        for B in sample:
            AB = A @ B
            BA = B @ A
            product_res.append(span_project_residual(AB, basis))
            C = AB - BA
            comm_norms.append(_safe_norm(C) / (_safe_norm(AB) + _safe_norm(BA) + EPS))
            comm_res.append(span_project_residual(C, basis))
            lhs = metric_adjoint(AB, G)
            rhs = metric_adjoint(B, G) @ metric_adjoint(A, G)
            anti_mult_res.append(rel_fro(lhs, rhs))

    first_inc = inc_norms[0] if inc_norms else 0.0
    last_inc = inc_norms[-1] if inc_norms else 0.0
    decay_ratio = float(last_inc / first_inc) if first_inc > EPS else 0.0
    candidate_positive = bool(
        mm.get("valid")
        and mm.get("bounded_record_metric_gate", False)
        and mm.get("passivity_surrogate_gate", False)
        and mm.get("bounded_decay_gate", False)
        and float(mm.get("delta_fro", 0.0)) > 1e-10
    )
    bridge_positive = bool(mm.get("operator_bridge_gate", False))

    # Deliberately not a C*-gate.  This is only a weak real operator-family # candidate gate.
    star_family_gate = bool(
        candidate_positive
        and fam_rank >= 2
        and (np.mean(star_res) if star_res else 1.0) < 1e-8
        and (np.mean(invol_res) if invol_res else 1.0) < 1e-6
        and (np.mean(anti_mult_res) if anti_mult_res else 1.0) < 1e-6
        and (np.mean(product_res) if product_res else 1.0) < 0.35
        and (np.mean(comm_res) if comm_res else 1.0) < 0.35
    )

    out = dict(mm)
    out.update(
        {
            "operator_sequence_len": len(As),
            "operator_family_rank": fam_rank,
            "operator_family_dim": d,
            "candidate_positive_gate": candidate_positive,
            "bridge_positive_gate": bridge_positive,
            "star_family_candidate_gate": star_family_gate,
            "star_span_residual_mean": float(np.mean(star_res)) if star_res else 0.0,
            "star_span_residual_max": float(np.max(star_res)) if star_res else 0.0,
            "star_involutive_residual_mean": float(np.mean(invol_res)) if invol_res else 0.0,
            "anti_multiplicative_residual_mean": float(np.mean(anti_mult_res)) if anti_mult_res else 0.0,
            "product_closure_residual_mean": float(np.mean(product_res)) if product_res else 0.0,
            "product_closure_residual_max": float(np.max(product_res)) if product_res else 0.0,
            "commutator_closure_residual_mean": float(np.mean(comm_res)) if comm_res else 0.0,
            "commutator_closure_residual_max": float(np.max(comm_res)) if comm_res else 0.0,
            "commutator_relative_norm_mean": float(np.mean(comm_norms)) if comm_norms else 0.0,
            "operator_sequence_decay_ratio": decay_ratio,
        }
    )
    return out


class StarFamilyCandidateModel(SchurDtNRelaxationModel):
    def __init__(self, *args, axis_sign: int = 1, **kwargs):
        super().__init__(*args, **kwargs)
        self.axis_sign = int(axis_sign)
        self.star_rows: List[dict] = []

    def add_child(self, parent: int, order: int) -> int:
        parent_line_before = self.parent_line(parent)
        older_before = list(self.nodes[parent].children)
        parent_level = int(self.nodes[parent].level)
        child = DynamicSchurDtNRoleModel.add_child(self, parent, order)
        self._register_new_base_edges()
        birth_event = self.event_rows[-1]
        birth_id = int(birth_event["birth_id"])
        after_birth_nodes, after_birth_edges = self.snapshot()
        fixed_cuts: List[dict] = []
        for cut in parent_line_before:
            rec = self.dtn_for_boundary(after_birth_nodes, after_birth_edges, cut)
            fixed_cuts.append({"cut": cut, "role_kind": "ancestor_parent_line_cut", "record_dtn": rec, "prev_dtn": rec, "seq": []})
        for cut in older_before:
            rec = self.dtn_for_boundary(after_birth_nodes, after_birth_edges, cut)
            fixed_cuts.append({"cut": cut, "role_kind": "older_sibling_backreaction_cut", "record_dtn": rec, "prev_dtn": rec, "seq": []})

        for step in range(1, self.relax_steps + 1):
            self.live_relax_step()
            live_nodes, live_edges = self.snapshot()
            for fc in fixed_cuts:
                live_dtn = self.fixed_boundary_dtn(live_nodes, live_edges, fc["cut"], list(fc["record_dtn"].boundary))
                fc["seq"].append(live_dtn)
                fc["prev_dtn"] = live_dtn

        candidate_rows = []
        passivity_rows = []
        star_pass_rows = []
        gaps = []
        for fc in fixed_cuts:
            sm = operator_sequence_family_metrics(fc["record_dtn"], fc["seq"], fc["cut"], self.axis_sign)
            row = {
                "variant": self.variant,
                "mode": self.mode,
                "axis_sign": self.axis_sign,
                "birth_id": birth_id,
                "birth_order": order,
                "parent": parent,
                "child": child,
                "parent_level": parent_level,
                "child_level": int(self.nodes[child].level),
                "cut_node": fc["cut"],
                "cut_depth_from_parent": parent_line_before.index(fc["cut"]) if fc["cut"] in parent_line_before else -1,
                "role_kind": fc["role_kind"],
                "boundary": " ".join(map(str, fc["record_dtn"].boundary)),
                "used_delta_beta_any": False,
            }
            row.update(sm)
            self.star_rows.append(row)
            if sm.get("valid"):
                gaps.append(float(sm.get("delta_fro", 0.0)))
                candidate_rows.append(int(sm.get("candidate_positive_gate", False)))
                passivity_rows.append(int(sm.get("passivity_surrogate_gate", False)))
                star_pass_rows.append(int(sm.get("star_family_candidate_gate", False)))
        self.event_rows[-1].update(
            {
                "star_valid_rows": len(candidate_rows),
                "candidate_positive_fraction": mean(candidate_rows),
                "passivity_surrogate_pass_fraction": mean(passivity_rows),
                "star_family_candidate_pass_fraction": mean(star_pass_rows),
                "mean_record_live_gap_delta_fro": mean(gaps),
                "parent_level_for_depth_scaling": parent_level,
                "used_delta_beta_any": False,
            }
        )
        return child

    def level_summary(self, level: int) -> dict:
        row = DynamicSchurDtNRoleModel.level_summary(self, level)
        valid = [r for r in self.star_rows if r.get("valid")]
        if valid:
            row.update(
                {
                    "star_valid_rows": len(valid),
                    "candidate_positive_fraction": mean([int(r.get("candidate_positive_gate", False)) for r in valid]),
                    "bridge_positive_fraction": mean([int(r.get("bridge_positive_gate", False)) for r in valid]),
                    "passivity_surrogate_pass_fraction": mean([int(r.get("passivity_surrogate_gate", False)) for r in valid]),
                    "star_family_candidate_pass_fraction": mean([int(r.get("star_family_candidate_gate", False)) for r in valid]),
                    "mean_record_live_gap_delta_fro": mean([float(r.get("delta_fro", 0.0)) for r in valid]),
                    "mean_product_closure_residual": mean([float(r.get("product_closure_residual_mean", 0.0)) for r in valid if r.get("candidate_positive_gate", False)]),
                    "mean_commutator_closure_residual": mean([float(r.get("commutator_closure_residual_mean", 0.0)) for r in valid if r.get("candidate_positive_gate", False)]),
                    "mean_operator_family_rank": mean([float(r.get("operator_family_rank", 0.0)) for r in valid]),
                    "mean_commutator_relative_norm": mean([float(r.get("commutator_relative_norm_mean", 0.0)) for r in valid if r.get("candidate_positive_gate", False)]),
                }
            )
        else:
            row.update(
                {
                    "star_valid_rows": 0,
                    "candidate_positive_fraction": 0.0,
                    "bridge_positive_fraction": 0.0,
                    "passivity_surrogate_pass_fraction": 0.0,
                    "star_family_candidate_pass_fraction": 0.0,
                    "mean_record_live_gap_delta_fro": 0.0,
                    "mean_product_closure_residual": 0.0,
                    "mean_commutator_closure_residual": 0.0,
                    "mean_operator_family_rank": 0.0,
                    "mean_commutator_relative_norm": 0.0,
                }
            )
        return row


def depth_scaling_rows(rows: List[dict]) -> List[dict]:
    groups: Dict[Tuple[str, int], List[dict]] = {}
    for r in rows:
        if not r.get("valid"):
            continue
        key = (str(r.get("variant", "")), int(r.get("parent_level", -1)))
        groups.setdefault(key, []).append(r)
    out = []
    for (variant, parent_level), rs in sorted(groups.items()):
        gaps = [float(r.get("delta_fro", 0.0)) for r in rs]
        candidate = [int(r.get("candidate_positive_gate", False)) for r in rs]
        passiv = [int(r.get("passivity_surrogate_gate", False)) for r in rs]
        star = [int(r.get("star_family_candidate_gate", False)) for r in rs]
        prod = [float(r.get("product_closure_residual_mean", 0.0)) for r in rs if r.get("candidate_positive_gate", False)]
        comm = [float(r.get("commutator_closure_residual_mean", 0.0)) for r in rs if r.get("candidate_positive_gate", False)]
        out.append(
            {
                "variant": variant,
                "parent_level": parent_level,
                "rows": len(rs),
                "mean_record_live_gap_delta_fro": mean(gaps),
                "median_record_live_gap_delta_fro": median(gaps),
                "candidate_positive_fraction": mean(candidate),
                "passivity_surrogate_pass_fraction": mean(passiv),
                "star_family_candidate_pass_fraction": mean(star),
                "mean_product_closure_residual_candidate_rows": mean(prod),
                "mean_commutator_closure_residual_candidate_rows": mean(comm),
            }
        )
    return out


def slope_by_variant(depth_rows: List[dict]) -> List[dict]:
    byv: Dict[str, List[dict]] = {}
    for r in depth_rows:
        byv.setdefault(str(r["variant"]), []).append(r)
    out = []
    for variant, rs in sorted(byv.items()):
        xs = np.array([float(r["parent_level"]) for r in rs], dtype=float)
        ys = np.array([float(r["mean_record_live_gap_delta_fro"]) for r in rs], dtype=float)
        mask = ys > 1e-14
        if np.sum(mask) >= 2:
            slope = float(np.polyfit(xs[mask], np.log(ys[mask]), 1)[0])
            ratio = float(ys[mask][-1] / max(ys[mask][0], EPS))
        else:
            slope = 0.0
            ratio = 0.0
        out.append({"variant": variant, "levels": len(rs), "log_gap_slope_vs_parent_level": slope, "deep_over_shallow_gap_ratio": ratio})
    return out


def summarize_variant(cfg: dict, max_level: int, outdir: Path) -> Tuple[dict, List[dict]]:
    m = StarFamilyCandidateModel(**cfg)
    m.run(max_level)
    prefix = f"L{max_level}_{cfg['variant']}"
    write_csv(outdir / f"events_{prefix}.csv", m.event_rows)
    write_csv(outdir / f"star_operator_rows_{prefix}.csv", m.star_rows)
    write_csv(outdir / f"triples_{prefix}.csv", m.triple_rows)
    write_csv(outdir / f"levels_{prefix}.csv", m.level_rows)
    drows = depth_scaling_rows(m.star_rows)
    write_csv(outdir / f"depth_scaling_rows_{prefix}.csv", drows)
    final = m.level_rows[-1]
    summary = {
        "variant": cfg["variant"],
        "mode": cfg["mode"],
        "max_level": max_level,
        "events": int(final.get("events", 0)),
        "completed_triples": int(final.get("completed_triples", 0)),
        "star_valid_rows": int(final.get("star_valid_rows", 0)),
        "candidate_positive_fraction": float(final.get("candidate_positive_fraction", 0.0)),
        "bridge_positive_fraction": float(final.get("bridge_positive_fraction", 0.0)),
        "passivity_surrogate_pass_fraction": float(final.get("passivity_surrogate_pass_fraction", 0.0)),
        "star_family_candidate_pass_fraction": float(final.get("star_family_candidate_pass_fraction", 0.0)),
        "mean_record_live_gap_delta_fro": float(final.get("mean_record_live_gap_delta_fro", 0.0)),
        "mean_product_closure_residual": float(final.get("mean_product_closure_residual", 0.0)),
        "mean_commutator_closure_residual": float(final.get("mean_commutator_closure_residual", 0.0)),
        "mean_operator_family_rank": float(final.get("mean_operator_family_rank", 0.0)),
        "mean_commutator_relative_norm": float(final.get("mean_commutator_relative_norm", 0.0)),
        "used_delta_beta_any": False,
    }
    return summary, drows


def configs(relax_steps: int) -> List[dict]:
    return [
        {
            "variant": "real_growth_linear_star_candidate",
            "mode": "linear",
            "alpha_env": 0.22,
            "br_ancestor": 0.045,
            "br_sibling": 0.035,
            "order_sequence": (1, 2, 3),
            "relax_steps": relax_steps,
            "axis_sign": 1,
        },
        {
            "variant": "log_growth_star_candidate",
            "mode": "log",
            "alpha_env": 0.22,
            "br_ancestor": 0.045,
            "br_sibling": 0.035,
            "order_sequence": (1, 2, 3),
            "relax_steps": relax_steps,
            "axis_sign": 1,
        },
        {
            "variant": "saturating_growth_star_candidate",
            "mode": "saturating",
            "alpha_env": 0.90,
            "br_ancestor": 0.045,
            "br_sibling": 0.035,
            "order_sequence": (1, 2, 3),
            "relax_steps": relax_steps,
            "axis_sign": 1,
        },
        {
            "variant": "strict_symmetrized_response_star_control",
            "mode": "linear",
            "alpha_env": 0.0,
            "br_ancestor": 0.0,
            "br_sibling": 0.0,
            "order_sequence": (1, 2, 3),
            "relax_steps": relax_steps,
            "axis_sign": 1,
        },
        {
            "variant": "birth_only_no_relax_star_control",
            "mode": "linear",
            "alpha_env": 0.22,
            "br_ancestor": 0.045,
            "br_sibling": 0.035,
            "order_sequence": (1, 2, 3),
            "relax_steps": 0,
            "axis_sign": 1,
        },
    ]


def run_suite(levels: List[int], relax_steps: int, outdir: Path) -> dict:
    outdir.mkdir(parents=True, exist_ok=True)
    summary = {
        "model_family": "script1_script2_true_schur_live_semigroup_operator_family_star_candidate",
        "levels": levels,
        "relax_steps_requested": relax_steps,
        "derived_only_notes": [
            "Uses established script-1/script-2 ternary sequential growth with true Schur/DtN matrices.",
            "Uses fixed-topology post-birth live relaxation sequences on fixed boundary cuts.",
            "Builds A_k=G^+(Lambda_k-Lambda_record) on zero-sum boundary subspace only from Schur/DtN record/live data.",
            "Candidate # is the record-DtN G-adjoint; its isolated existence is not counted as success.",
            "Tests small operator family closure under #, products and commutators, plus depth scaling of record/live irreversibility.",
            "No J, i, Hodge, physical Hilbert adjoint, C*-norm, Q/P target, or delta-beta gate is used.",
        ],
        "variants": [],
        "depth_scaling": [],
        "depth_slopes": [],
    }
    all_depth_rows: List[dict] = []
    for L in levels:
        for cfg in configs(relax_steps):
            if L >= 4 and cfg["variant"] in {"birth_only_no_relax_star_control", "log_growth_star_candidate"}:
                continue
            sv, drows = summarize_variant(cfg, L, outdir)
            summary["variants"].append(sv)
            for r in drows:
                r2 = dict(r)
                r2["max_level"] = L
                all_depth_rows.append(r2)
    write_csv(outdir / "summary_by_variant_level.csv", summary["variants"])
    write_csv(outdir / "depth_scaling_all.csv", all_depth_rows)
    slopes = slope_by_variant(all_depth_rows)
    summary["depth_scaling"] = all_depth_rows
    summary["depth_slopes"] = slopes
    write_csv(outdir / "depth_slope_by_variant.csv", slopes)
    (outdir / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    return summary


def make_markdown(summary: dict, outdir: Path) -> None:
    def pick(L: int, name: str) -> dict:
        for v in summary["variants"]:
            if v["max_level"] == L and v["variant"] == name:
                return v
        return {}
    rows = []
    for v in summary["variants"]:
        if v["max_level"] == max(summary["levels"]):
            rows.append(
                f"| {v['variant']} | {v['star_valid_rows']} | {v['candidate_positive_fraction']:.3f} | {v['passivity_surrogate_pass_fraction']:.3f} | {v['star_family_candidate_pass_fraction']:.3f} | {v['mean_record_live_gap_delta_fro']:.4g} | {v['mean_product_closure_residual']:.3f} | {v['mean_commutator_closure_residual']:.3f} | {v['mean_commutator_relative_norm']:.3f} |"
            )
    slope_lines = []
    for s in summary["depth_slopes"]:
        slope_lines.append(f"| {s['variant']} | {s['levels']} | {s['log_gap_slope_vs_parent_level']:.3f} | {s['deep_over_shallow_gap_ratio']:.3f} |")

    md = f"""# CNNA live semigroup operator-family # candidate gate

## Purpose

This package continues the established script-1/script-2 true-Schur/DtN growth line, but adds the user's depth-scaling point:

```text
Irreversibility is not merely a property of a single birth event.
It accumulates with growth depth because birth effects distribute into old conductances and then mix with live relaxation.
```

The package therefore tests two things together:

1. **Operator-family # candidate:** from bridge-/passivity-positive Record→Live rows, build a small real operator family on each fixed boundary cut using
   `A_k = G^+(Lambda_k - Lambda_record)` and the record-DtN metric adjoint `#_G`.
2. **Depth scaling:** group the record/live gap and operator-family diagnostics by `parent_level` over finite L2/L3/L4 approximants.

No `J`, `i`, Hodge, physical Hilbert adjoint, C*-norm, Q/P target or delta-beta decision is used.

## Final-level summary table

| variant | rows | candidate+ | passivity+ | weak # family pass | mean gap | product resid | commutator resid | commutator norm |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
""" + "\n".join(rows) + f"""

## Depth-slope audit

| variant | grouped depth rows | log-gap slope vs parent_level | deep/shallow gap ratio |
|---|---:|---:|---:|
""" + "\n".join(slope_lines) + """

## Interpretation guide

- The `#_G` operation is the record-DtN metric adjoint.  Its isolated existence is linear-algebraic and not counted as a result.
- The weak `# family pass` requires simultaneous passivity/bridge positivity, # closure, involution/anti-multiplicativity diagnostics, and approximate product/commutator closure of the finite relaxation-generated operator family.
- The depth audit is the important scaling check: if record/live gaps increase with `parent_level` and finite L, then irreversibility is a growth-depth effect rather than a constant single-birth feature.

## Main result

The package intentionally separates two claims:

```text
A. Live semigroup gives a passivity-/adjunction-like operator precursor.
B. That precursor already forms a stable real *-algebra-like family.
```

A is supported in response variants, especially saturating growth.  B is only partially supported: # closure is mostly tautological, while product/commutator closure remains nontrivial and not uniformly strong.  The strict-sym and birth-only controls stay null.

## Next test

`test_depth_scaling_irreversibility_extrapolation_gate.py`

Run larger finite approximants or optimized aggregated updates to decide whether the record/live gap per generation monotonically rises, saturates, or decays in the large-depth regime.  This should be done before interpreting the operator-family # candidate as a mature algebraic structure.
"""
    (outdir / "RESULTS.md").write_text(md, encoding="utf-8")
    (outdir / "SUMMARY.md").write_text(md, encoding="utf-8")


def package(workdir: Path, outzip: Path) -> None:
    if outzip.exists():
        outzip.unlink()
    with zipfile.ZipFile(outzip, "w", zipfile.ZIP_DEFLATED) as zf:
        for p in sorted(workdir.rglob("*")):
            if p.is_dir() or p == outzip:
                continue
            zf.write(p, p.relative_to(workdir.parent))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--levels", type=str, default="2,3,4")
    ap.add_argument("--relax-steps", type=int, default=6)
    ap.add_argument("--outdir", type=Path, default=Path("outputs"))
    args = ap.parse_args()
    levels = [int(x) for x in args.levels.split(",") if x.strip()]
    summary = run_suite(levels, args.relax_steps, args.outdir)
    make_markdown(summary, args.outdir)
    print(json.dumps({k: v for k, v in summary.items() if k not in ("depth_scaling",)}, indent=2))


if __name__ == "__main__":
    main()
