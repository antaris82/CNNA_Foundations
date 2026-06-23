#!/usr/bin/env python3
"""
CNNA live-semigroup polarity -> real operator bridge gate, L2/L4 default.

Purpose
-------
Use the directed Record->Live Schur/DtN semigroup polarity as input for a
real boundary-operator bridge, but do not search for J, i, spin, Hodge, a
physical star, or a C*-algebra.

The gate is intentionally weaker and earlier in the ladder:

  dynamic boundary roles + true Schur/DtN + live semigroup polarity
      -> candidate real adjunction/passivity/positivity-surrogate diagnostics
      -> only later: Q/P, #, J, complex structure.

What is tested
--------------
For each birth event and each fixed boundary cut, compute true DtN matrices
Lambda_0, Lambda_1, ..., Lambda_N on the same boundary port set during post-birth
live relaxation.  Project DtN matrices to the canonical zero-sum boundary voltage
subspace, because Laplace DtN maps have the constant vector in the kernel.

From these matrices derive, without target fitting:

- record metric surrogate G = projected Lambda_0;
- live drift Delta = projected Lambda_N - projected Lambda_0;
- longitudinal cut-vs-UV mode ell;
- polarity signed by the root/front axis orientation;
- G-adjoint residual for A = G^+ Delta;
- passivity surrogate: eigenvalue sign of axis_sign * Delta;
- longitudinal boundary-value energy sign: ell^T Delta ell;
- contraction/decay of live increments dLambda_k.

Important: a G-adjoint exists by linear algebra.  It is not counted as success by
itself.  The gate requires simultaneous nontriviality, bounded semigroup decay,
consistent longitudinal polarity, and a predominantly one-sided dissipative
spectral surrogate.  This is still not a physical positivity theorem.

No J, i, Hodge, physical star, Hilbert positivity, C*-norm, Q/P target, or
delta-beta is used as a gate.
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


def mean(xs):
    return float(np.mean(xs)) if xs else 0.0


def median(xs):
    return float(np.median(xs)) if xs else 0.0


def sym(A: np.ndarray) -> np.ndarray:
    return 0.5 * (A + A.T)


def zero_sum_basis(n: int) -> np.ndarray:
    """Return an n x (n-1) basis for {x: sum x = 0}."""
    if n <= 1:
        return np.zeros((n, 0), dtype=float)
    # Simple deterministic basis e_i - e_n, QR-orthonormalized.
    M = np.zeros((n, n - 1), dtype=float)
    for i in range(n - 1):
        M[i, i] = 1.0
        M[n - 1, i] = -1.0
    Q, _ = np.linalg.qr(M)
    return Q[:, : n - 1]


def project_zero_sum(L: np.ndarray) -> np.ndarray:
    Z = zero_sum_basis(L.shape[0])
    if Z.shape[1] == 0:
        return np.zeros((0, 0), dtype=float)
    return sym(Z.T @ sym(L) @ Z)


def longitudinal_mode(boundary: List[int], cut: int) -> np.ndarray:
    n = len(boundary)
    v = np.zeros(n, dtype=float)
    if n <= 1:
        return v
    try:
        ci = boundary.index(cut)
    except ValueError:
        ci = 0
    v[:] = -1.0 / (n - 1)
    v[ci] = 1.0
    # Already zero-sum. Normalize for scale-invariant energy signs.
    norm = float(np.linalg.norm(v))
    if norm > EPS:
        v /= norm
    return v


def rel_fro(A: np.ndarray, B: np.ndarray | None = None) -> float:
    if B is None:
        return float(np.linalg.norm(A, ord="fro"))
    return float(np.linalg.norm(A - B, ord="fro") / (np.linalg.norm(B, ord="fro") + EPS))


def matrix_metrics(record: DtNResult, seq: List[DtNResult], cut: int, axis_sign: int) -> dict:
    """Analyze one fixed-boundary Record->Live sequence."""
    if not seq:
        return {"valid": False}
    boundary = list(record.boundary)
    n = len(boundary)
    d = max(0, n - 1)
    if n <= 1 or record.lam.shape[0] != n:
        return {"valid": False, "reason": "boundary_dim_le_1"}
    L0 = sym(np.asarray(record.lam, dtype=float))
    LN = sym(np.asarray(seq[-1].lam, dtype=float))
    G = project_zero_sum(L0)
    H = project_zero_sum(LN)
    if G.size == 0:
        return {"valid": False, "reason": "zero_projected_dim"}
    Delta = sym(H - G)
    # Derived longitudinal mode projected into zero-sum orthonormal basis.
    Z = zero_sum_basis(n)
    ell_full = longitudinal_mode(boundary, cut)
    ell = Z.T @ ell_full
    if np.linalg.norm(ell) > EPS:
        ell = ell / np.linalg.norm(ell)

    # Regularize only for diagnostic inverse in the zero-sum subspace.
    evals_G = np.linalg.eigvalsh(sym(G)) if G.size else np.array([])
    g_min = float(np.min(evals_G)) if evals_G.size else 0.0
    g_max = float(np.max(evals_G)) if evals_G.size else 0.0
    g_cond = float(g_max / max(g_min, 1e-12)) if evals_G.size else 0.0
    reg = max(1e-10, 1e-8 * max(g_max, 1.0))
    Gp = np.linalg.pinv(G + reg * np.eye(G.shape[0]), rcond=1e-10)
    A = Gp @ Delta
    GA = G @ A
    adj_res = float(np.linalg.norm(GA - A.T @ G, ord="fro") / (np.linalg.norm(GA, ord="fro") + EPS)) if GA.size else 0.0
    anti_adj_res = float(np.linalg.norm(GA + A.T @ G, ord="fro") / (np.linalg.norm(GA, ord="fro") + EPS)) if GA.size else 0.0

    signed_Delta = float(axis_sign) * Delta
    ev = np.linalg.eigvalsh(sym(signed_Delta)) if signed_Delta.size else np.array([])
    pos = int(np.sum(ev > 1e-10))
    neg = int(np.sum(ev < -1e-10))
    nonzero = pos + neg
    pos_frac = float(pos / nonzero) if nonzero else 0.0
    neg_frac = float(neg / nonzero) if nonzero else 0.0
    eig_balance = float(abs(pos - neg) / nonzero) if nonzero else 0.0
    # Passivity polarity is whichever one-sided sign dominates.  The axis sign tells
    # which direction is retarded; the result may be mostly positive or mostly negative.
    one_sided_frac = max(pos_frac, neg_frac)
    dominant_sign = 1 if pos_frac >= neg_frac else -1

    ell_energy = float(axis_sign * (ell.T @ Delta @ ell)) if ell.size else 0.0
    ell_energy_raw = float(ell.T @ Delta @ ell) if ell.size else 0.0
    # The dissipative direction can be positive or negative by convention; gate records bias not a fixed sign.

    # Increment decay and sign consistency along relaxation sequence.
    inc_norms = []
    inc_ell = []
    prev = G
    for item in seq:
        P = project_zero_sum(sym(np.asarray(item.lam, dtype=float)))
        if P.shape == prev.shape:
            D = sym(P - prev)
            inc_norms.append(float(np.linalg.norm(D, ord="fro")))
            inc_ell.append(float(axis_sign * (ell.T @ D @ ell)) if ell.size else 0.0)
            prev = P
    first_inc = inc_norms[0] if inc_norms else 0.0
    last_inc = inc_norms[-1] if inc_norms else 0.0
    decay_ratio = float(last_inc / first_inc) if first_inc > EPS else 0.0
    inc_signs = [math.copysign(1.0, x) for x in inc_ell if abs(x) > 1e-12]
    inc_sign_bias = float(abs(sum(inc_signs)) / len(inc_signs)) if inc_signs else 0.0

    # Operator bridge pass is intentionally stringent but not a claim of positivity theorem.
    nontrivial = float(np.linalg.norm(Delta, ord="fro")) > 1e-10
    bounded_record_metric = bool(g_min > -1e-8 and g_cond < 1e10)
    selfadjoint_like = adj_res < 1e-8  # expected for A=G^+Delta with symmetric Delta; logged but not enough.
    passivity_surrogate = bool(one_sided_frac >= 0.80 and abs(ell_energy) > 1e-12 and inc_sign_bias >= 0.80)
    bounded_decay = bool(decay_ratio < 0.35 or first_inc < 1e-12)
    bridge_pass = bool(nontrivial and bounded_record_metric and selfadjoint_like and passivity_surrogate and bounded_decay)

    return {
        "valid": True,
        "boundary_dim": n,
        "projected_dim": d,
        "g_min_eig": g_min,
        "g_max_eig": g_max,
        "g_cond": g_cond,
        "delta_fro": float(np.linalg.norm(Delta, ord="fro")),
        "delta_op": float(np.linalg.norm(Delta, ord=2)) if Delta.size else 0.0,
        "G_adjoint_residual": adj_res,
        "G_anti_adjoint_residual": anti_adj_res,
        "signed_delta_pos_frac": pos_frac,
        "signed_delta_neg_frac": neg_frac,
        "signed_delta_one_sided_frac": one_sided_frac,
        "signed_delta_dominant_sign": dominant_sign,
        "signed_delta_eig_balance": eig_balance,
        "longitudinal_energy_raw": ell_energy_raw,
        "oriented_longitudinal_energy": ell_energy,
        "increment_first_fro": first_inc,
        "increment_last_fro": last_inc,
        "increment_decay_ratio": decay_ratio,
        "increment_longitudinal_sign_bias": inc_sign_bias,
        "bounded_record_metric_gate": bounded_record_metric,
        "selfadjoint_metric_gate": selfadjoint_like,
        "passivity_surrogate_gate": passivity_surrogate,
        "bounded_decay_gate": bounded_decay,
        "operator_bridge_gate": bridge_pass,
    }


class OperatorBridgeModel(SchurDtNRelaxationModel):
    def __init__(self, *args, axis_sign: int = 1, **kwargs):
        super().__init__(*args, **kwargs)
        self.axis_sign = int(axis_sign)
        self.bridge_rows: List[dict] = []

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
            fixed_cuts.append({"cut": cut, "role_kind": "ancestor_parent_line_cut", "record_dtn": rec, "prev_dtn": rec, "seq": []})
        for cut in older_before:
            rec = self.dtn_for_boundary(after_birth_nodes, after_birth_edges, cut)
            fixed_cuts.append({"cut": cut, "role_kind": "older_sibling_backreaction_cut", "record_dtn": rec, "prev_dtn": rec, "seq": []})

        for step in range(1, self.relax_steps + 1):
            relax_state = self.live_relax_step()
            live_nodes, live_edges = self.snapshot()
            for fc in fixed_cuts:
                live_dtn = self.fixed_boundary_dtn(live_nodes, live_edges, fc["cut"], list(fc["record_dtn"].boundary))
                fc["seq"].append(live_dtn)
                step_fro, step_op, ok_step = self.dtn_delta(fc["prev_dtn"], live_dtn)
                gap_fro, gap_op, ok_gap = self.dtn_delta(fc["record_dtn"], live_dtn)
                self.relax_rows.append(
                    {
                        "variant": self.variant,
                        "mode": self.mode,
                        "axis_sign": self.axis_sign,
                        "birth_id": birth_id,
                        "relax_step": step,
                        "cut_node": fc["cut"],
                        "role_kind": fc["role_kind"],
                        "parent": parent,
                        "child": child,
                        "fixed_boundary_dim": len(fc["record_dtn"].boundary),
                        "fixed_boundary": " ".join(map(str, fc["record_dtn"].boundary)),
                        "record_eff_tail_conductance": fc["record_dtn"].eff_tail_conductance,
                        "live_eff_tail_conductance": live_dtn.eff_tail_conductance,
                        "record_to_live_eff_gap": live_dtn.eff_tail_conductance - fc["record_dtn"].eff_tail_conductance,
                        "oriented_eff_gap": self.axis_sign * (live_dtn.eff_tail_conductance - fc["record_dtn"].eff_tail_conductance),
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
                )
                fc["prev_dtn"] = live_dtn

        final_nodes, final_edges = self.snapshot()
        bridge_passes = []
        passivity_passes = []
        one_sideds = []
        ell_energies = []
        decay_ratios = []
        for fc in fixed_cuts:
            mm = matrix_metrics(fc["record_dtn"], fc["seq"], fc["cut"], self.axis_sign)
            row = {
                "variant": self.variant,
                "mode": self.mode,
                "axis_sign": self.axis_sign,
                "birth_id": birth_id,
                "parent": parent,
                "child": child,
                "cut_node": fc["cut"],
                "role_kind": fc["role_kind"],
                "boundary": " ".join(map(str, fc["record_dtn"].boundary)),
                "used_delta_beta_any": False,
            }
            row.update(mm)
            self.bridge_rows.append(row)
            if mm.get("valid"):
                bridge_passes.append(int(mm.get("operator_bridge_gate", False)))
                passivity_passes.append(int(mm.get("passivity_surrogate_gate", False)))
                one_sideds.append(float(mm.get("signed_delta_one_sided_frac", 0.0)))
                ell_energies.append(float(mm.get("oriented_longitudinal_energy", 0.0)))
                decay_ratios.append(float(mm.get("increment_decay_ratio", 0.0)))

        # Attach event-level bridge summary.
        signs = [math.copysign(1.0, e) for e in ell_energies if abs(e) > 1e-12]
        pol_bias = abs(sum(signs)) / len(signs) if signs else 0.0
        self.event_rows[-1].update(
            {
                "operator_bridge_valid_rows": len(bridge_passes),
                "operator_bridge_pass_fraction": mean(bridge_passes),
                "passivity_surrogate_pass_fraction": mean(passivity_passes),
                "mean_signed_delta_one_sided_frac": mean(one_sideds),
                "mean_increment_decay_ratio": mean(decay_ratios),
                "longitudinal_energy_polarity_bias": pol_bias,
                "used_delta_beta_any": False,
            }
        )
        return child

    def level_summary(self, level: int) -> dict:
        row = DynamicSchurDtNRoleModel.level_summary(self, level)
        valid = [r for r in self.bridge_rows if r.get("valid")]
        if valid:
            row.update(
                {
                    "operator_bridge_rows": len(valid),
                    "operator_bridge_pass_fraction": mean([int(r.get("operator_bridge_gate", False)) for r in valid]),
                    "passivity_surrogate_pass_fraction": mean([int(r.get("passivity_surrogate_gate", False)) for r in valid]),
                    "bounded_record_metric_fraction": mean([int(r.get("bounded_record_metric_gate", False)) for r in valid]),
                    "selfadjoint_metric_fraction": mean([int(r.get("selfadjoint_metric_gate", False)) for r in valid]),
                    "bounded_decay_fraction": mean([int(r.get("bounded_decay_gate", False)) for r in valid]),
                    "mean_signed_delta_one_sided_frac": mean([float(r.get("signed_delta_one_sided_frac", 0.0)) for r in valid]),
                    "mean_increment_longitudinal_sign_bias": mean([float(r.get("increment_longitudinal_sign_bias", 0.0)) for r in valid]),
                    "mean_increment_decay_ratio": mean([float(r.get("increment_decay_ratio", 0.0)) for r in valid]),
                    "mean_abs_oriented_longitudinal_energy": mean([abs(float(r.get("oriented_longitudinal_energy", 0.0))) for r in valid]),
                    "positive_longitudinal_energy_fraction": mean([int(float(r.get("oriented_longitudinal_energy", 0.0)) > 1e-12) for r in valid]),
                    "negative_longitudinal_energy_fraction": mean([int(float(r.get("oriented_longitudinal_energy", 0.0)) < -1e-12) for r in valid]),
                }
            )
        else:
            row.update(
                {
                    "operator_bridge_rows": 0,
                    "operator_bridge_pass_fraction": 0.0,
                    "passivity_surrogate_pass_fraction": 0.0,
                    "bounded_record_metric_fraction": 0.0,
                    "selfadjoint_metric_fraction": 0.0,
                    "bounded_decay_fraction": 0.0,
                    "mean_signed_delta_one_sided_frac": 0.0,
                    "mean_increment_longitudinal_sign_bias": 0.0,
                    "mean_increment_decay_ratio": 0.0,
                    "mean_abs_oriented_longitudinal_energy": 0.0,
                    "positive_longitudinal_energy_fraction": 0.0,
                    "negative_longitudinal_energy_fraction": 0.0,
                }
            )
        return row


def summarize_variant(cfg: dict, max_level: int, outdir: Path) -> dict:
    m = OperatorBridgeModel(**cfg)
    m.run(max_level)
    prefix = f"L{max_level}_{cfg['variant']}"
    write_csv(outdir / f"events_{prefix}.csv", m.event_rows)
    write_csv(outdir / f"boundary_role_rows_{prefix}.csv", m.role_rows)
    write_csv(outdir / f"dtn_schur_rows_{prefix}.csv", m.dtn_rows)
    write_csv(outdir / f"relax_rows_{prefix}.csv", m.relax_rows)
    write_csv(outdir / f"operator_bridge_rows_{prefix}.csv", m.bridge_rows)
    write_csv(outdir / f"triples_{prefix}.csv", m.triple_rows)
    write_csv(outdir / f"levels_{prefix}.csv", m.level_rows)
    final = m.level_rows[-1]
    summary = {
        "variant": cfg["variant"],
        "mode": cfg["mode"],
        "max_level": max_level,
        "events": int(final.get("events", 0)),
        "completed_triples": int(final.get("completed_triples", 0)),
        "retarded_event_fraction": float(final.get("retarded_event_fraction", 0.0)),
        "advanced_leakage_fraction": float(final.get("advanced_leakage_fraction", 0.0)),
        "mean_abs_log_circulation": float(final.get("mean_abs_log_circulation", 0.0)),
        "frac_full_markov_complex": float(final.get("frac_full_markov_complex", 0.0)),
        "operator_bridge_rows": int(final.get("operator_bridge_rows", 0)),
        "operator_bridge_pass_fraction": float(final.get("operator_bridge_pass_fraction", 0.0)),
        "passivity_surrogate_pass_fraction": float(final.get("passivity_surrogate_pass_fraction", 0.0)),
        "bounded_record_metric_fraction": float(final.get("bounded_record_metric_fraction", 0.0)),
        "selfadjoint_metric_fraction": float(final.get("selfadjoint_metric_fraction", 0.0)),
        "bounded_decay_fraction": float(final.get("bounded_decay_fraction", 0.0)),
        "mean_signed_delta_one_sided_frac": float(final.get("mean_signed_delta_one_sided_frac", 0.0)),
        "mean_increment_longitudinal_sign_bias": float(final.get("mean_increment_longitudinal_sign_bias", 0.0)),
        "mean_increment_decay_ratio": float(final.get("mean_increment_decay_ratio", 0.0)),
        "mean_abs_oriented_longitudinal_energy": float(final.get("mean_abs_oriented_longitudinal_energy", 0.0)),
        "positive_longitudinal_energy_fraction": float(final.get("positive_longitudinal_energy_fraction", 0.0)),
        "negative_longitudinal_energy_fraction": float(final.get("negative_longitudinal_energy_fraction", 0.0)),
        "used_delta_beta_any": False,
    }
    return summary


def run_suite(max_level: int, relax_steps: int, outdir: Path) -> dict:
    outdir.mkdir(parents=True, exist_ok=True)
    configs = [
        {
            "variant": "real_growth_linear_live_polarity_operator_bridge",
            "mode": "linear",
            "alpha_env": 0.22,
            "br_ancestor": 0.045,
            "br_sibling": 0.035,
            "order_sequence": (1, 2, 3),
            "relax_steps": relax_steps,
            "axis_sign": 1,
        },
        {
            "variant": "log_growth_live_polarity_operator_bridge",
            "mode": "log",
            "alpha_env": 0.22,
            "br_ancestor": 0.045,
            "br_sibling": 0.035,
            "order_sequence": (1, 2, 3),
            "relax_steps": relax_steps,
            "axis_sign": 1,
        },
        {
            "variant": "saturating_growth_live_polarity_operator_bridge",
            "mode": "saturating",
            "alpha_env": 0.90,
            "br_ancestor": 0.045,
            "br_sibling": 0.035,
            "order_sequence": (1, 2, 3),
            "relax_steps": relax_steps,
            "axis_sign": 1,
        },
        {
            "variant": "kappa_reversed_birth_order_operator_bridge_control",
            "mode": "linear",
            "alpha_env": 0.22,
            "br_ancestor": 0.045,
            "br_sibling": 0.035,
            "order_sequence": (3, 2, 1),
            "relax_steps": relax_steps,
            "axis_sign": 1,
        },
        {
            "variant": "longitudinal_axis_flip_operator_bridge_control",
            "mode": "linear",
            "alpha_env": 0.22,
            "br_ancestor": 0.045,
            "br_sibling": 0.035,
            "order_sequence": (1, 2, 3),
            "relax_steps": relax_steps,
            "axis_sign": -1,
        },
        {
            "variant": "strict_symmetrized_response_operator_bridge_control",
            "mode": "linear",
            "alpha_env": 0.0,
            "br_ancestor": 0.0,
            "br_sibling": 0.0,
            "order_sequence": (1, 2, 3),
            "relax_steps": relax_steps,
            "axis_sign": 1,
        },
        {
            "variant": "birth_only_no_relax_operator_bridge_control",
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
        "model_family": "script1_script2_true_schur_dtn_live_semigroup_polarity_to_operator_bridge",
        "max_level": max_level,
        "relax_steps_requested": relax_steps,
        "derived_only_notes": [
            "Uses established ternary sequential script-1/script-2 growth model with true Schur/DtN matrices.",
            "Uses Record->Live live-semigroup polarity as source, not J/i/spin/QP target.",
            "Projects DtN matrices to zero-sum boundary voltage subspace; derives G=Lambda_record and Delta=Lambda_live-Lambda_record.",
            "Tests only candidate real metric-adjunction, passivity/dissipation polarity, and positivity-surrogate diagnostics.",
            "No Hodge, physical star, Hilbert positivity, C*-norm, Fourier convention, i, J, Q/P target, or delta-beta gate is used.",
        ],
        "variants": [],
    }
    for cfg in configs:
        summary["variants"].append(summarize_variant(cfg, max_level, outdir))
    write_csv(outdir / "summary_by_variant.csv", summary["variants"])
    (outdir / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    return summary


def make_markdown(summary: dict, outdir: Path) -> None:
    lines = []
    for v in summary["variants"]:
        lines.append(
            "| {variant} | {rows} | {bridge:.3f} | {passiv:.3f} | {metric:.3f} | {selfadj:.3f} | {decay:.3f} | {onesided:.3f} | {signbias:.3f} | {pos:.3f} | {neg:.3f} |".format(
                variant=v["variant"],
                rows=v["operator_bridge_rows"],
                bridge=v["operator_bridge_pass_fraction"],
                passiv=v["passivity_surrogate_pass_fraction"],
                metric=v["bounded_record_metric_fraction"],
                selfadj=v["selfadjoint_metric_fraction"],
                decay=v["bounded_decay_fraction"],
                onesided=v["mean_signed_delta_one_sided_frac"],
                signbias=v["mean_increment_longitudinal_sign_bias"],
                pos=v["positive_longitudinal_energy_fraction"],
                neg=v["negative_longitudinal_energy_fraction"],
            )
        )
    md = f"""# CNNA live-semigroup polarity to operator bridge gate

Package: `cnna_live_semigroup_polarity_to_operator_bridge_gate_pkg_L2`

## Question

Can the directed Record->Live Schur/DtN live-semigroup polarity prepare a real operator bridge before any J, i, spin, Hodge, physical star, Q/P target, or C*-positivity is asserted?

This is a deliberately weak but simultaneous gate.  It does not count a linear-algebraic adjoint by itself.  It asks whether the following appear together:

```text
nontrivial Record->Live DtN drift
bounded record DtN metric on the zero-sum boundary-voltage subspace
metric self-adjointness of the drift operator A = G^+ Delta
one-sided dissipative / passivity-surrogate spectrum of signed Delta
stable longitudinal polarity of the relaxation increments
bounded decay of live Schur/DtN increments
```

## Method

For each birth event and fixed boundary cut:

```text
G     = zero-sum projection of Lambda_record
Delta = zero-sum projection of Lambda_live_final - Lambda_record
ell   = derived cut-vs-UV longitudinal boundary mode
A     = G^+ Delta
```

Then the test records:

- `G_adjoint_residual` for A;
- one-sided signed spectrum of `axis_sign * Delta`;
- longitudinal energy sign `ell^T Delta ell`;
- decay of live increments `Delta Lambda_k`;
- strict-sym and birth-only controls.

This is not a physical positivity theorem.  It is a bridge audit from live-boundary polarity toward possible later real adjunction/passivity structure.

## Summary table

| variant | rows | bridge pass | passivity pass | bounded G | selfadjoint | decay | one-sided eig | incr sign bias | +long E | -long E |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
""" + "\n".join(lines) + """

## Interpretation guide

- `bridge pass` is the simultaneous weak operator-bridge gate.
- `passivity pass` means signed live drift has one-sided spectral tendency and stable longitudinal increment polarity.
- `selfadjoint` alone is expected because the construction uses real Schur/DtN symmetric matrices; it is logged but not a success criterion by itself.
- A positive `strict_sym` result would be suspicious.  The expected result is zero rows/pass there.

## Next test

If the bridge gate is positive, the next test should examine whether the bridge induces a stable real `#` candidate on an operator family.  If only passivity/one-sided polarity is positive but bridge pass is weak, the next step is to strengthen the boundary-value/semigroup layer before claiming any operator structure.
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
    ap.add_argument("--max-level", type=int, default=4)
    ap.add_argument("--relax-steps", type=int, default=6)
    ap.add_argument("--outdir", type=Path, default=Path("outputs"))
    args = ap.parse_args()
    summary = run_suite(args.max_level, args.relax_steps, args.outdir)
    make_markdown(summary, args.outdir)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
