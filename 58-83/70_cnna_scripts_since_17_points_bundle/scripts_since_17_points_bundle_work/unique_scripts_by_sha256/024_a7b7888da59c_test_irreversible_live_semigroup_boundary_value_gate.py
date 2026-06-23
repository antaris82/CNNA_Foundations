#!/usr/bin/env python3
"""
CNNA irreversible live semigroup / boundary-value gate.

Builds on the established script-1/script-2 growth model with true Schur/DtN
matrices and fixed-topology live relaxation.  The gate does not ask whether a
child can be deleted.  It asks whether the live Schur/DtN evolution after birth
behaves like a directed, bounded, retarded semigroup layer rather than a stable
invertible/reversible dynamics.

No J, i, Hodge, physical star, positivity/C*-norm, Q/P target, or delta-beta is
used as a gate.  All boundary response maps are real Laplace Schur complements.
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

from base_dynamic_boundary_role_schur_dtn import DynamicSchurDtNRoleModel, Node, DtNResult, EPS, write_csv
from base_schur_dtn_birth_relaxation import SchurDtNRelaxationModel


def mean(xs):
    return float(np.mean(xs)) if xs else 0.0


def median(xs):
    return float(np.median(xs)) if xs else 0.0


def safe_norm(v):
    a = np.asarray(v, dtype=float)
    return float(np.linalg.norm(a)) if a.size else 0.0


def dtn_feature(r: DtNResult) -> np.ndarray:
    lam = np.asarray(r.lam, dtype=float)
    if lam.size:
        fro = float(np.linalg.norm(lam, ord="fro"))
        op = float(np.linalg.norm(lam, ord=2))
        trace = float(np.trace(lam))
        offdiag = lam.copy()
        np.fill_diagonal(offdiag, 0.0)
        offdiag_norm = float(np.linalg.norm(offdiag, ord="fro"))
        neg_offdiag_mass = float(np.sum(np.maximum(0.0, -offdiag)))
        pos_diag_mass = float(np.sum(np.maximum(0.0, np.diag(lam))))
    else:
        fro = op = trace = offdiag_norm = neg_offdiag_mass = pos_diag_mass = 0.0
    return np.array(
        [
            float(len(r.boundary)),
            float(len(r.uv_ports)),
            float(r.eff_tail_conductance),
            float(r.total_cut_to_uv_coupling),
            fro,
            op,
            trace,
            offdiag_norm,
            neg_offdiag_mass,
            pos_diag_mass,
            float(r.matrix_rank),
            float(r.min_eig),
            float(r.max_eig),
            float(r.cond_like),
            float(r.singular_flag),
        ],
        dtype=float,
    )


FEATURE_NAMES = [
    "boundary_dim",
    "uv_port_count",
    "eff_tail_conductance",
    "total_cut_to_uv_coupling",
    "fro_norm",
    "op_norm",
    "trace",
    "offdiag_norm",
    "neg_offdiag_mass",
    "pos_diag_mass",
    "matrix_rank",
    "min_eig",
    "max_eig",
    "cond_like",
    "singular_flag",
]


class SemigroupCaptureModel(SchurDtNRelaxationModel):
    def __init__(self, *args, axis_sign: int = 1, **kwargs):
        super().__init__(*args, **kwargs)
        self.axis_sign = int(axis_sign)
        self.semigroup_rows: List[dict] = []

    def add_child(self, parent: int, order: int) -> int:
        parent_line_before = self.parent_line(parent)
        older_before = list(self.nodes[parent].children)
        child = DynamicSchurDtNRoleModel.add_child(self, parent, order)
        self._register_new_base_edges()

        birth_event = self.event_rows[-1]
        birth_id = int(birth_event["birth_id"])
        after_birth_nodes, after_birth_edges = self.snapshot()
        fixed_cuts: List[dict] = []
        for cut in parent_line_before:
            rec = self.dtn_for_boundary(after_birth_nodes, after_birth_edges, cut)
            fixed_cuts.append(
                {
                    "cut": cut,
                    "role_kind": "ancestor_parent_line_cut",
                    "fixed_boundary": list(rec.boundary),
                    "record_dtn": rec,
                    "prev_dtn": rec,
                    "first_step_delta_fro": None,
                    "last_step_delta_fro": 0.0,
                    "total_relax_drift_fro": 0.0,
                }
            )
        for cut in older_before:
            rec = self.dtn_for_boundary(after_birth_nodes, after_birth_edges, cut)
            fixed_cuts.append(
                {
                    "cut": cut,
                    "role_kind": "older_sibling_backreaction_cut",
                    "fixed_boundary": list(rec.boundary),
                    "record_dtn": rec,
                    "prev_dtn": rec,
                    "first_step_delta_fro": None,
                    "last_step_delta_fro": 0.0,
                    "total_relax_drift_fro": 0.0,
                }
            )

        # Step-0 feature rows are the birth records.  These are the record layer.
        for fc in fixed_cuts:
            feat = dtn_feature(fc["record_dtn"])
            row = {
                "variant": self.variant,
                "mode": self.mode,
                "axis_sign": self.axis_sign,
                "birth_id": birth_id,
                "relax_step": 0,
                "cut_node": fc["cut"],
                "role_kind": fc["role_kind"],
                "parent": parent,
                "child": child,
                "fixed_boundary_dim": len(fc["fixed_boundary"]),
                "fixed_boundary": " ".join(map(str, fc["fixed_boundary"])),
                "record_eff_tail_conductance": fc["record_dtn"].eff_tail_conductance,
                "live_eff_tail_conductance": fc["record_dtn"].eff_tail_conductance,
                "record_to_live_eff_gap": 0.0,
                "oriented_eff_gap": 0.0,
                "relax_delta_dtn_fro": 0.0,
                "record_vs_live_gap_fro": 0.0,
                "node_delta_norm": 0.0,
                "edge_delta_norm": 0.0,
            }
            for name, val in zip(FEATURE_NAMES, feat):
                row[f"feat_{name}"] = float(val)
            self.semigroup_rows.append(row)

        for step in range(1, self.relax_steps + 1):
            relax_state = self.live_relax_step()
            live_nodes, live_edges = self.snapshot()
            for fc in fixed_cuts:
                live_dtn = self.fixed_boundary_dtn(live_nodes, live_edges, fc["cut"], fc["fixed_boundary"])
                step_fro, step_op, ok_step = self.dtn_delta(fc["prev_dtn"], live_dtn)
                gap_fro, gap_op, ok_gap = self.dtn_delta(fc["record_dtn"], live_dtn)
                if fc["first_step_delta_fro"] is None:
                    fc["first_step_delta_fro"] = step_fro
                fc["last_step_delta_fro"] = step_fro
                fc["total_relax_drift_fro"] += step_fro
                eff_gap = live_dtn.eff_tail_conductance - fc["record_dtn"].eff_tail_conductance
                row = {
                    "variant": self.variant,
                    "mode": self.mode,
                    "axis_sign": self.axis_sign,
                    "birth_id": birth_id,
                    "relax_step": step,
                    "cut_node": fc["cut"],
                    "role_kind": fc["role_kind"],
                    "parent": parent,
                    "child": child,
                    "fixed_boundary_dim": len(fc["fixed_boundary"]),
                    "fixed_boundary": " ".join(map(str, fc["fixed_boundary"])),
                    "record_eff_tail_conductance": fc["record_dtn"].eff_tail_conductance,
                    "live_eff_tail_conductance": live_dtn.eff_tail_conductance,
                    "record_to_live_eff_gap": eff_gap,
                    "oriented_eff_gap": self.axis_sign * eff_gap,
                    "relax_delta_dtn_fro": step_fro,
                    "relax_delta_dtn_op": step_op,
                    "record_vs_live_gap_fro": gap_fro,
                    "record_vs_live_gap_op": gap_op,
                    "dtn_delta_shape_ok": int(ok_step and ok_gap),
                    "node_delta_norm": relax_state["node_delta_norm"],
                    "edge_delta_norm": relax_state["edge_delta_norm"],
                    "node_state_norm": relax_state["node_state_norm"],
                    "edge_state_norm": relax_state["edge_state_norm"],
                    "topology_changed_during_relaxation": 0,
                    "boundary_ports_changed_during_relaxation": 0,
                }
                feat = dtn_feature(live_dtn)
                for name, val in zip(FEATURE_NAMES, feat):
                    row[f"feat_{name}"] = float(val)
                self.semigroup_rows.append(row)
                self.relax_rows.append(row.copy())
                fc["prev_dtn"] = live_dtn

        final_nodes, final_edges = self.snapshot()
        firsts = [fc["first_step_delta_fro"] or 0.0 for fc in fixed_cuts]
        lasts = [fc["last_step_delta_fro"] for fc in fixed_cuts]
        totals = [fc["total_relax_drift_fro"] for fc in fixed_cuts]
        gaps = []
        eff_gaps = []
        for fc in fixed_cuts:
            live_dtn = self.fixed_boundary_dtn(final_nodes, final_edges, fc["cut"], fc["fixed_boundary"])
            gap_fro, _gap_op, ok = self.dtn_delta(fc["record_dtn"], live_dtn)
            if ok:
                gaps.append(gap_fro)
            eff_gaps.append(live_dtn.eff_tail_conductance - fc["record_dtn"].eff_tail_conductance)
        convergence_ratios = [l / f for f, l in zip(firsts, lasts) if f > EPS]
        signs = [math.copysign(1.0, self.axis_sign * g) for g in eff_gaps if abs(g) > 1e-12]
        halfplane_bias = abs(sum(signs)) / len(signs) if signs else 0.0
        self.relax_event_rows.append(
            {
                "variant": self.variant,
                "mode": self.mode,
                "axis_sign": self.axis_sign,
                "birth_id": birth_id,
                "parent": parent,
                "child": child,
                "relax_steps": self.relax_steps,
                "fixed_cut_count": len(fixed_cuts),
                "birth_delta_parent_dtn_eff": birth_event.get("parent_schur_dtn_eff_delta", 0.0),
                "birth_new_child_port_coupling": birth_event.get("parent_new_child_port_coupling", 0.0),
                "mean_first_relax_delta_dtn_fro": mean(firsts),
                "mean_last_relax_delta_dtn_fro": mean(lasts),
                "mean_total_relax_drift_dtn_fro": mean(totals),
                "mean_record_vs_live_gap_fro": mean(gaps),
                "max_record_vs_live_gap_fro": max(gaps) if gaps else 0.0,
                "mean_relax_convergence_ratio_last_over_first": mean(convergence_ratios),
                "relaxation_active_fraction": mean([int(g > 1e-10) for g in gaps]),
                "quasi_equilibrium_fraction": mean([int((l / f) < 0.35) for f, l in zip(firsts, lasts) if f > EPS]),
                "halfplane_bias_oriented_eff_gap": halfplane_bias,
                "positive_oriented_eff_gap_fraction": mean([int(self.axis_sign * g > 1e-12) for g in eff_gaps]),
                "negative_oriented_eff_gap_fraction": mean([int(self.axis_sign * g < -1e-12) for g in eff_gaps]),
                "topology_changed_during_relaxation": 0,
                "boundary_ports_changed_during_relaxation": 0,
                "used_delta_beta_any": False,
            }
        )
        self.event_rows[-1].update(self.relax_event_rows[-1])
        return child

    def level_summary(self, level: int) -> dict:
        row = DynamicSchurDtNRoleModel.level_summary(self, level)
        relevant = [r for r in self.relax_event_rows]
        if relevant:
            row.update(
                {
                    "relax_events": len(relevant),
                    "mean_first_relax_delta_dtn_fro": mean([float(r["mean_first_relax_delta_dtn_fro"]) for r in relevant]),
                    "mean_last_relax_delta_dtn_fro": mean([float(r["mean_last_relax_delta_dtn_fro"]) for r in relevant]),
                    "mean_total_relax_drift_dtn_fro": mean([float(r["mean_total_relax_drift_dtn_fro"]) for r in relevant]),
                    "mean_record_vs_live_gap_fro": mean([float(r["mean_record_vs_live_gap_fro"]) for r in relevant]),
                    "max_record_vs_live_gap_fro": max(float(r["max_record_vs_live_gap_fro"]) for r in relevant),
                    "mean_relax_convergence_ratio_last_over_first": mean([float(r["mean_relax_convergence_ratio_last_over_first"]) for r in relevant]),
                    "relaxation_active_fraction": mean([float(r["relaxation_active_fraction"]) for r in relevant]),
                    "quasi_equilibrium_fraction": mean([float(r["quasi_equilibrium_fraction"]) for r in relevant]),
                    "mean_halfplane_bias_oriented_eff_gap": mean([float(r["halfplane_bias_oriented_eff_gap"]) for r in relevant]),
                    "mean_positive_oriented_eff_gap_fraction": mean([float(r["positive_oriented_eff_gap_fraction"]) for r in relevant]),
                    "mean_negative_oriented_eff_gap_fraction": mean([float(r["negative_oriented_eff_gap_fraction"]) for r in relevant]),
                }
            )
        else:
            row.update(
                {
                    "relax_events": 0,
                    "mean_first_relax_delta_dtn_fro": 0.0,
                    "mean_last_relax_delta_dtn_fro": 0.0,
                    "mean_total_relax_drift_dtn_fro": 0.0,
                    "mean_record_vs_live_gap_fro": 0.0,
                    "max_record_vs_live_gap_fro": 0.0,
                    "mean_relax_convergence_ratio_last_over_first": 0.0,
                    "relaxation_active_fraction": 0.0,
                    "quasi_equilibrium_fraction": 0.0,
                    "mean_halfplane_bias_oriented_eff_gap": 0.0,
                    "mean_positive_oriented_eff_gap_fraction": 0.0,
                    "mean_negative_oriented_eff_gap_fraction": 0.0,
                }
            )
        return row


def feature_matrix(rows: List[dict], step: int = None) -> np.ndarray:
    if step is not None:
        rows = [r for r in rows if int(r["relax_step"]) == step]
    X = []
    for r in rows:
        X.append([float(r.get(f"feat_{name}", 0.0)) for name in FEATURE_NAMES])
    return np.asarray(X, dtype=float)


def fit_map(X: np.ndarray, Y: np.ndarray) -> dict:
    if X.shape[0] < 4 or X.shape != Y.shape:
        return {
            "ok": False,
            "rows": int(X.shape[0]),
            "forward_rel_err": 0.0,
            "spectral_norm": 0.0,
            "min_singular": 0.0,
            "condition": 0.0,
            "rank": 0,
        }
    # Standardize feature columns by X scale only: this is an audit, not a target fit.
    scale = np.sqrt(np.mean(X * X, axis=0)) + 1e-12
    Xs = X / scale
    Ys = Y / scale
    F = np.linalg.pinv(Xs, rcond=1e-10) @ Ys
    pred = Xs @ F
    err = float(np.linalg.norm(pred - Ys, ord="fro") / (np.linalg.norm(Ys, ord="fro") + EPS))
    s = np.linalg.svd(F, compute_uv=False)
    return {
        "ok": True,
        "rows": int(X.shape[0]),
        "forward_rel_err": err,
        "spectral_norm": float(np.max(s)) if s.size else 0.0,
        "min_singular": float(np.min(s)) if s.size else 0.0,
        "condition": float((np.max(s) / max(np.min(s), 1e-12))) if s.size else 0.0,
        "rank": int(np.sum(s > 1e-8)),
        "singular_values": s.tolist(),
    }


def analyze_semigroup(rows: List[dict], relax_steps: int) -> dict:
    if not rows or relax_steps <= 0:
        return {
            "transition_rows": 0,
            "forward_fit_rel_err": 0.0,
            "forward_spectral_norm": 0.0,
            "forward_min_singular": 0.0,
            "forward_condition": 0.0,
            "forward_rank": 0,
            "reverse_fit_rel_err": 0.0,
            "reverse_spectral_norm": 0.0,
            "reverse_min_singular": 0.0,
            "reverse_condition": 0.0,
            "reverse_rank": 0,
            "median_last_over_first_relax_delta": 0.0,
            "mean_last_over_first_relax_delta": 0.0,
            "mean_final_record_vs_live_gap_fro": 0.0,
            "max_final_record_vs_live_gap_fro": 0.0,
            "halfplane_bias_final_oriented_eff_gap": 0.0,
            "positive_final_oriented_eff_gap_fraction": 0.0,
            "negative_final_oriented_eff_gap_fraction": 0.0,
            "median_feature_norm_ratio_step_to_next": 0.0,
            "forward_bounded_contracting_gate": False,
            "reverse_unstable_or_nonreconstructable_gate": False,
            "directed_live_semigroup_boundary_gate": False,
            "semigroup_interpretation": "no_semigroup_gate",
        }
    pairX = []
    pairY = []
    # Match by event/cut/role between consecutive relax steps.
    by_key = {}
    for r in rows:
        key = (r["birth_id"], r["cut_node"], r["role_kind"], int(r["relax_step"]))
        by_key[key] = r
    step_delta_ratios = []
    record_gap_final = []
    eff_gap_signs = []
    norm_ratios = []
    for r in rows:
        k = int(r["relax_step"])
        if k >= relax_steps:
            continue
        key2 = (r["birth_id"], r["cut_node"], r["role_kind"], k + 1)
        if key2 not in by_key:
            continue
        r2 = by_key[key2]
        x = np.array([float(r.get(f"feat_{name}", 0.0)) for name in FEATURE_NAMES], dtype=float)
        y = np.array([float(r2.get(f"feat_{name}", 0.0)) for name in FEATURE_NAMES], dtype=float)
        pairX.append(x)
        pairY.append(y)
        nx = float(np.linalg.norm(x))
        ny = float(np.linalg.norm(y))
        if nx > EPS:
            norm_ratios.append(ny / nx)
    X = np.asarray(pairX, dtype=float)
    Y = np.asarray(pairY, dtype=float)
    fwd = fit_map(X, Y)
    rev = fit_map(Y, X)
    # Event-level decay: compare first and last nonzero relaxation deltas per event/cut.
    event_cut = {}
    for r in rows:
        key = (r["birth_id"], r["cut_node"], r["role_kind"])
        event_cut.setdefault(key, []).append(r)
    for key, seq in event_cut.items():
        seq = sorted(seq, key=lambda z: int(z["relax_step"]))
        nonzero = [float(s.get("relax_delta_dtn_fro", 0.0)) for s in seq if int(s["relax_step"]) > 0]
        if nonzero and nonzero[0] > EPS:
            step_delta_ratios.append(nonzero[-1] / nonzero[0])
        final = seq[-1]
        record_gap_final.append(float(final.get("record_vs_live_gap_fro", 0.0)))
        eg = float(final.get("oriented_eff_gap", 0.0))
        if abs(eg) > 1e-12:
            eff_gap_signs.append(math.copysign(1.0, eg))
    halfplane_bias = abs(sum(eff_gap_signs)) / len(eff_gap_signs) if eff_gap_signs else 0.0
    positive_fraction = mean([int(x > 0) for x in eff_gap_signs])
    negative_fraction = mean([int(x < 0) for x in eff_gap_signs])
    median_decay = median(step_delta_ratios)
    mean_decay = mean(step_delta_ratios)
    fwd_contract = bool(fwd.get("ok") and fwd["spectral_norm"] <= 1.15 and median_decay < 0.35)
    # Reverse is unstable/nonreconstructable if reverse condition is large, reverse fit is poor,
    # or forward has collapsed singular directions while record/live gap is nonzero.
    reverse_bad = bool(
        rev.get("ok")
        and (
            rev["condition"] > 25.0
            or rev["forward_rel_err"] > max(0.20, 1.75 * fwd.get("forward_rel_err", 0.0))
            or fwd.get("min_singular", 1.0) < 0.08
        )
    )
    directed_gate = bool(
        fwd_contract
        and mean(record_gap_final) > 1e-10
        and halfplane_bias > 0.80
        and len(eff_gap_signs) > 0
    )
    return {
        "transition_rows": int(X.shape[0]),
        "forward_fit_rel_err": float(fwd.get("forward_rel_err", 0.0)),
        "forward_spectral_norm": float(fwd.get("spectral_norm", 0.0)),
        "forward_min_singular": float(fwd.get("min_singular", 0.0)),
        "forward_condition": float(fwd.get("condition", 0.0)),
        "forward_rank": int(fwd.get("rank", 0)),
        "reverse_fit_rel_err": float(rev.get("forward_rel_err", 0.0)),
        "reverse_spectral_norm": float(rev.get("spectral_norm", 0.0)),
        "reverse_min_singular": float(rev.get("min_singular", 0.0)),
        "reverse_condition": float(rev.get("condition", 0.0)),
        "reverse_rank": int(rev.get("rank", 0)),
        "median_last_over_first_relax_delta": float(median_decay),
        "mean_last_over_first_relax_delta": float(mean_decay),
        "mean_final_record_vs_live_gap_fro": float(mean(record_gap_final)),
        "max_final_record_vs_live_gap_fro": float(max(record_gap_final) if record_gap_final else 0.0),
        "halfplane_bias_final_oriented_eff_gap": float(halfplane_bias),
        "positive_final_oriented_eff_gap_fraction": float(positive_fraction),
        "negative_final_oriented_eff_gap_fraction": float(negative_fraction),
        "median_feature_norm_ratio_step_to_next": float(median(norm_ratios)),
        "forward_bounded_contracting_gate": bool(fwd_contract),
        "reverse_unstable_or_nonreconstructable_gate": bool(reverse_bad),
        "directed_live_semigroup_boundary_gate": bool(directed_gate),
        "semigroup_interpretation": "bounded_retarded_live_semigroup_candidate" if directed_gate else "no_semigroup_gate",
    }


def run_variant(cfg: dict, max_level: int, outdir: Path) -> Tuple[dict, SemigroupCaptureModel]:
    m = SemigroupCaptureModel(**cfg)
    m.run(max_level)
    prefix = f"L{max_level}_{cfg['variant']}"
    write_csv(outdir / f"events_{prefix}.csv", m.event_rows)
    write_csv(outdir / f"boundary_role_rows_{prefix}.csv", m.role_rows)
    write_csv(outdir / f"dtn_schur_rows_{prefix}.csv", m.dtn_rows)
    write_csv(outdir / f"semigroup_rows_{prefix}.csv", m.semigroup_rows)
    write_csv(outdir / f"relax_event_rows_{prefix}.csv", m.relax_event_rows)
    write_csv(outdir / f"triples_{prefix}.csv", m.triple_rows)
    write_csv(outdir / f"levels_{prefix}.csv", m.level_rows)
    sg = analyze_semigroup(m.semigroup_rows, int(cfg.get("relax_steps", 0)))
    final = m.level_rows[-1]
    summary = {
        "variant": cfg["variant"],
        "mode": cfg["mode"],
        "max_level": max_level,
        "events": int(final.get("events", 0)),
        "relax_events": int(final.get("relax_events", 0)),
        "completed_triples": int(final.get("completed_triples", 0)),
        "retarded_event_fraction": float(final.get("retarded_event_fraction", 0.0)),
        "advanced_leakage_fraction": float(final.get("advanced_leakage_fraction", 0.0)),
        "mean_abs_log_circulation": float(final.get("mean_abs_log_circulation", 0.0)),
        "frac_full_markov_complex": float(final.get("frac_full_markov_complex", 0.0)),
        "mean_record_vs_live_gap_fro": float(final.get("mean_record_vs_live_gap_fro", 0.0)),
        "mean_total_relax_drift_dtn_fro": float(final.get("mean_total_relax_drift_dtn_fro", 0.0)),
        "mean_relax_convergence_ratio_last_over_first": float(final.get("mean_relax_convergence_ratio_last_over_first", 0.0)),
        "mean_halfplane_bias_oriented_eff_gap": float(final.get("mean_halfplane_bias_oriented_eff_gap", 0.0)),
        "used_delta_beta_any": False,
    }
    summary.update(sg)
    return summary, m


def run_suite(max_level: int, relax_steps: int, outdir: Path) -> dict:
    outdir.mkdir(parents=True, exist_ok=True)
    configs = [
        {
            "variant": "real_growth_linear_true_schur_live_semigroup",
            "mode": "linear",
            "alpha_env": 0.22,
            "br_ancestor": 0.045,
            "br_sibling": 0.035,
            "order_sequence": (1, 2, 3),
            "relax_steps": relax_steps,
            "axis_sign": 1,
        },
        {
            "variant": "log_growth_true_schur_live_semigroup",
            "mode": "log",
            "alpha_env": 0.22,
            "br_ancestor": 0.045,
            "br_sibling": 0.035,
            "order_sequence": (1, 2, 3),
            "relax_steps": relax_steps,
            "axis_sign": 1,
        },
        {
            "variant": "saturating_growth_true_schur_live_semigroup",
            "mode": "saturating",
            "alpha_env": 0.90,
            "br_ancestor": 0.045,
            "br_sibling": 0.035,
            "order_sequence": (1, 2, 3),
            "relax_steps": relax_steps,
            "axis_sign": 1,
        },
        {
            "variant": "kappa_reversed_birth_order_true_schur_live_semigroup_control",
            "mode": "linear",
            "alpha_env": 0.22,
            "br_ancestor": 0.045,
            "br_sibling": 0.035,
            "order_sequence": (3, 2, 1),
            "relax_steps": relax_steps,
            "axis_sign": 1,
        },
        {
            "variant": "longitudinal_axis_flip_live_semigroup_control",
            "mode": "linear",
            "alpha_env": 0.22,
            "br_ancestor": 0.045,
            "br_sibling": 0.035,
            "order_sequence": (1, 2, 3),
            "relax_steps": relax_steps,
            "axis_sign": -1,
        },
        {
            "variant": "strict_symmetrized_response_control_true_schur_live_semigroup",
            "mode": "linear",
            "alpha_env": 0.0,
            "br_ancestor": 0.0,
            "br_sibling": 0.0,
            "order_sequence": (1, 2, 3),
            "relax_steps": relax_steps,
            "axis_sign": 1,
        },
        {
            "variant": "birth_only_no_relax_record_control_true_schur_semigroup",
            "mode": "linear",
            "alpha_env": 0.22,
            "br_ancestor": 0.045,
            "br_sibling": 0.035,
            "order_sequence": (1, 2, 3),
            "relax_steps": 0,
            "axis_sign": 1,
        },
    ]
    summary = {
        "model_family": "script1_script2_true_schur_dtn_irreversible_live_semigroup_boundary_value",
        "max_level": max_level,
        "relax_steps_requested": relax_steps,
        "derived_only_notes": [
            "Uses established ternary sequential script-1/script-2 growth model with true Schur/DtN matrices.",
            "Treats fixed-topology post-birth relaxation as the candidate live semigroup layer.",
            "Builds real DtN feature transition maps from record/live Schur-complement responses; no J/i/Hodge/star/QP/delta-beta gate is used.",
            "Boundary-value polarity is audited through oriented record-to-live DtN effective-tail gaps under normal and longitudinal-axis-flipped interpretations.",
        ],
        "variants": [],
    }
    all_sg_rows = []
    for cfg in configs:
        s, m = run_variant(cfg, max_level, outdir)
        summary["variants"].append(s)
        all_sg_rows.extend(m.semigroup_rows)
    write_csv(outdir / "all_semigroup_rows.csv", all_sg_rows)
    write_csv(outdir / "summary_by_variant.csv", summary["variants"])
    (outdir / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    return summary


def make_markdown(summary: dict, outdir: Path) -> None:
    rows = []
    for v in summary["variants"]:
        rows.append(
            f"| {v['variant']} | {v['directed_live_semigroup_boundary_gate']} | {v['forward_bounded_contracting_gate']} | {v['reverse_unstable_or_nonreconstructable_gate']} | {v['transition_rows']} | {v['forward_spectral_norm']:.3g} | {v['median_last_over_first_relax_delta']:.3g} | {v['mean_final_record_vs_live_gap_fro']:.3g} | {v['halfplane_bias_final_oriented_eff_gap']:.3g} |"
        )
    md = f"""# CNNA irreversible live semigroup / boundary-value gate

Package: `cnna_irreversible_live_semigroup_boundary_value_gate_pkg_L2`

## Question

Does the true Schur/DtN live evolution after a birth event behave as a directed, bounded, retarded semigroup layer?

This test no longer asks whether the newborn can be deleted.  It checks whether the post-birth live Schur/DtN states form a forward-oriented record-to-live evolution with bounded contraction, nonzero record/live gap, and boundary-value polarity.

## Method

The established script-1/script-2 ternary growth model is kept:

- sequential sibling births;
- incoming environment response to the newborn;
- newborn backreaction to parent line and older siblings;
- true real Laplace Schur/DtN maps on boundary ports;
- fixed-topology live relaxation between births.

For every birth and every fixed cut, the DtN matrix is converted to a real invariant vector at relax steps `0..N`.  A feature-space forward map `T_forward` and reverse map `T_reverse` are fit as diagnostics only.  The gate is based on bounded forward evolution, decaying DtN increments, nonzero record/live gap, and stable oriented boundary polarity.  No J, i, Hodge, physical star, positivity, C*-norm, Q/P target, or delta-beta is used.

## Summary table

| variant | semigroup gate | fwd bounded | reverse bad | transitions | fwd norm | median last/first | record/live gap | halfplane bias |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
""" + "\n".join(rows) + f"""

## Interpretation

A positive semigroup gate means:

```text
record DtN state  ->  live DtN state
```

is directed and bounded/contractive in the live relaxation layer.  It does not mean that a complex structure, J, spin, a *-operation, or Q/P compatibility has been derived.

A negative strict-sym control means that topology alone is not enough: the response/backreaction layer is required.

The longitudinal-axis-flip control uses the same live evolution but reverses the interpretation of the root-front polarity.  This is not a new physical dynamics; it audits whether the boundary-value sign is relative to the chosen longitudinal orientation.

## Next test

`test_live_semigroup_polarity_to_operator_bridge_gate.py`

Use the directed live-semigroup polarity as an input to an operator bridge, but keep it separate from J.  The next gate should ask whether this record-to-live semigroup polarity can induce a candidate real adjunction or passivity/positivity surrogate before any complex structure is asserted.
"""
    (outdir / "RESULTS.md").write_text(md, encoding="utf-8")
    (outdir / "SUMMARY.md").write_text(md, encoding="utf-8")


def package(workdir: Path, outzip: Path) -> None:
    if outzip.exists():
        outzip.unlink()
    with zipfile.ZipFile(outzip, "w", zipfile.ZIP_DEFLATED) as zf:
        for p in sorted(workdir.rglob("*")):
            if p == outzip or p.is_dir():
                continue
            zf.write(p, p.relative_to(workdir.parent))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--max-level", type=int, default=4)
    ap.add_argument("--relax-steps", type=int, default=6)
    ap.add_argument("--outdir", type=Path, default=Path("outputs"))
    args = ap.parse_args()
    summary = run_suite(args.max_level, args.relax_steps, args.outdir)
    make_markdown(summary, args.outdir)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
