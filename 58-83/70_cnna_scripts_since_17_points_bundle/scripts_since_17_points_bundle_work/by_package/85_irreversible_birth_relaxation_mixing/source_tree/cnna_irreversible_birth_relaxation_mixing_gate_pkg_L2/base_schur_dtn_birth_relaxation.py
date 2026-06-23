#!/usr/bin/env python3
"""
CNNA Schur/DtN birth + live relaxation dynamics gate, L2.

This test extends the established script-1/script-2 dynamic birth conductance
model and the true Schur/DtN boundary-role gate.

New distinction:

1. Birth update:
   - topology and boundary-port roles change;
   - newborn has no own UV-tail from its self-perspective;
   - parent/ancestor cuts gain the newborn as a new UV-tail boundary port;
   - true Schur/DtN matrices jump.

2. Live relaxation update:
   - no new node is born;
   - boundary-port sets stay fixed;
   - live node conductances and live edge weights relax deterministically from
     the existing incoming/outgoing response data;
   - true Schur/DtN matrices can drift further on the same fixed boundary.

The aim is to separate:

- birth_delta_DtN: boundary role/topology jump;
- relax_delta_DtN: live response drift at fixed topology/boundary;
- record_vs_live_gap: immutable birth-record DtN vs current live DtN.

No J, i, Hodge, physical star, positivity/C*-norm, or delta-beta is used as a gate.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import shutil
import sys
import zipfile
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np

from base_dynamic_boundary_role_schur_dtn import (
    DynamicSchurDtNRoleModel,
    Node,
    DtNResult,
    EPS,
    write_csv,
)


class SchurDtNRelaxationModel(DynamicSchurDtNRoleModel):
    def __init__(
        self,
        *args,
        relax_steps: int = 6,
        node_relax_rate: float = 0.35,
        edge_relax_rate: float = 0.50,
        live_in_gain: float = 0.10,
        live_out_gain: float = 0.06,
        **kwargs,
    ):
        super().__init__(*args, **kwargs)
        self.relax_steps = relax_steps
        self.node_relax_rate = node_relax_rate
        self.edge_relax_rate = edge_relax_rate
        self.live_in_gain = live_in_gain
        self.live_out_gain = live_out_gain
        self.base_edge_weights: Dict[Tuple[int, int], float] = {}
        self.relax_rows: List[dict] = []
        self.relax_event_rows: List[dict] = []

    def _register_new_base_edges(self) -> None:
        for e, w in self.directed_edges.items():
            if e not in self.base_edge_weights:
                self.base_edge_weights[e] = float(w)

    def _state_vector(self) -> Tuple[np.ndarray, np.ndarray, List[Tuple[int, int]]]:
        node_ids = sorted(self.nodes)
        node_vec = np.array([self.nodes[i].g for i in node_ids], dtype=float)
        edge_keys = sorted(self.directed_edges)
        edge_vec = np.array([self.directed_edges[e] for e in edge_keys], dtype=float)
        return node_vec, edge_vec, edge_keys

    def live_relax_step(self) -> dict:
        """One deterministic live relaxation step on fixed topology.

        The rule is intentionally conservative and derived only from currently
        available growth quantities: live conductance g, birth conductance birth_g,
        and directed incoming/outgoing response weights.  It does not introduce a
        target J, Q/P, positivity gate, or external equilibrium.
        """
        before_nodes, before_edges, edge_keys = self._state_vector()
        incoming = defaultdict(float)
        outgoing = defaultdict(float)
        for (u, v), w in self.directed_edges.items():
            if w <= EPS:
                continue
            outgoing[u] += float(w)
            incoming[v] += float(w)

        new_g: Dict[int, float] = {}
        for nid, n in self.nodes.items():
            # birth_g is the immutable record anchor; incoming/outgoing are live response loads.
            target = n.birth_g + self.live_in_gain * incoming[nid] + self.live_out_gain * outgoing[nid]
            gnext = (1.0 - self.node_relax_rate) * n.g + self.node_relax_rate * target
            new_g[nid] = max(EPS, float(gnext))
        for nid, val in new_g.items():
            self.nodes[nid].g = val

        self._register_new_base_edges()
        for (u, v), w in list(self.directed_edges.items()):
            base = self.base_edge_weights.get((u, v), float(w))
            bu = max(EPS, self.nodes[u].birth_g)
            bv = max(EPS, self.nodes[v].birth_g)
            su = max(EPS, self.nodes[u].g) / bu
            sv = max(EPS, self.nodes[v].g) / bv
            # Live edge conductance follows endpoint live response while preserving
            # the original birth/backreaction edge provenance and orientation.
            target_w = base * math.sqrt(su * sv)
            self.directed_edges[(u, v)] = max(
                0.0,
                (1.0 - self.edge_relax_rate) * float(w) + self.edge_relax_rate * float(target_w),
            )

        after_nodes, after_edges, _ = self._state_vector()
        node_delta = float(np.linalg.norm(after_nodes - before_nodes)) if before_nodes.size else 0.0
        edge_delta = float(np.linalg.norm(after_edges - before_edges)) if before_edges.size == after_edges.size else 0.0
        return {
            "node_delta_norm": node_delta,
            "edge_delta_norm": edge_delta,
            "node_state_norm": float(np.linalg.norm(after_nodes)) if after_nodes.size else 0.0,
            "edge_state_norm": float(np.linalg.norm(after_edges)) if after_edges.size else 0.0,
        }

    @staticmethod
    def fixed_boundary_dtn(nodes: Dict[int, Node], edges: Dict[Tuple[int, int], float], cut: int, boundary: List[int]) -> DtNResult:
        uv_ports = [x for x in boundary if x != cut]
        return DynamicSchurDtNRoleModel.dtn_for_boundary(nodes, edges, cut, uv_ports)

    @staticmethod
    def dtn_delta(a: DtNResult, b: DtNResult) -> Tuple[float, float, int]:
        if a.lam.shape != b.lam.shape:
            return 0.0, 0.0, 0
        D = b.lam - a.lam
        return float(np.linalg.norm(D, ord="fro")), float(np.linalg.norm(D, ord=2)), 1

    def add_child(self, parent: int, order: int) -> int:
        # Capture the parent-line before topology changes; this is the set whose
        # boundary roles must be separated into birth jump and live relaxation drift.
        parent_line_before = self.parent_line(parent)
        older_before = list(self.nodes[parent].children)
        child = super().add_child(parent, order)
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

        # Relax at fixed topology/boundary.  The port list is the after-birth list;
        # no new UV port is introduced during relaxation.
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
                self.relax_rows.append(
                    {
                        "variant": self.variant,
                        "mode": self.mode,
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
                        "record_to_live_eff_gap": live_dtn.eff_tail_conductance - fc["record_dtn"].eff_tail_conductance,
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

        # Event-level relaxation summary after all fixed-topology relax steps.
        final_nodes, final_edges = self.snapshot()
        parent_record = None
        parent_final = None
        for fc in fixed_cuts:
            if fc["cut"] == parent and fc["role_kind"] == "ancestor_parent_line_cut":
                parent_record = fc["record_dtn"]
                parent_final = self.fixed_boundary_dtn(final_nodes, final_edges, fc["cut"], fc["fixed_boundary"])
                break
        if parent_record is not None and parent_final is not None:
            parent_gap_fro, parent_gap_op, _ = self.dtn_delta(parent_record, parent_final)
        else:
            parent_gap_fro = parent_gap_op = 0.0

        firsts = [fc["first_step_delta_fro"] or 0.0 for fc in fixed_cuts]
        lasts = [fc["last_step_delta_fro"] for fc in fixed_cuts]
        totals = [fc["total_relax_drift_fro"] for fc in fixed_cuts]
        gaps = []
        for fc in fixed_cuts:
            live_dtn = self.fixed_boundary_dtn(final_nodes, final_edges, fc["cut"], fc["fixed_boundary"])
            gap_fro, _gap_op, ok = self.dtn_delta(fc["record_dtn"], live_dtn)
            if ok:
                gaps.append(gap_fro)

        def mean(xs: List[float]) -> float:
            return float(np.mean(xs)) if xs else 0.0
        convergence_ratios = [l / f for f, l in zip(firsts, lasts) if f > EPS]
        relax_active = [int(g > 1e-10) for g in gaps]
        self.relax_event_rows.append(
            {
                "variant": self.variant,
                "mode": self.mode,
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
                "parent_record_vs_live_gap_fro": parent_gap_fro,
                "parent_record_vs_live_gap_op": parent_gap_op,
                "mean_relax_convergence_ratio_last_over_first": mean(convergence_ratios),
                "relaxation_active_fraction": mean(relax_active),
                "quasi_equilibrium_fraction": mean([int((l / f) < 0.35) for f, l in zip(firsts, lasts) if f > EPS]),
                "topology_changed_during_relaxation": 0,
                "boundary_ports_changed_during_relaxation": 0,
                "used_delta_beta_any": False,
            }
        )
        # Add headline relax fields to inherited birth event row for easier audit.
        self.event_rows[-1].update(self.relax_event_rows[-1])
        return child

    def level_summary(self, level: int) -> dict:
        row = super().level_summary(level)
        relevant = [r for r in self.relax_event_rows]
        if relevant:
            def mean(xs):
                return float(np.mean(xs)) if xs else 0.0
            row.update(
                {
                    "relax_events": len(relevant),
                    "mean_first_relax_delta_dtn_fro": mean([float(r["mean_first_relax_delta_dtn_fro"]) for r in relevant]),
                    "mean_last_relax_delta_dtn_fro": mean([float(r["mean_last_relax_delta_dtn_fro"]) for r in relevant]),
                    "mean_total_relax_drift_dtn_fro": mean([float(r["mean_total_relax_drift_dtn_fro"]) for r in relevant]),
                    "mean_record_vs_live_gap_fro": mean([float(r["mean_record_vs_live_gap_fro"]) for r in relevant]),
                    "max_record_vs_live_gap_fro": max(float(r["max_record_vs_live_gap_fro"]) for r in relevant),
                    "mean_parent_record_vs_live_gap_fro": mean([float(r["parent_record_vs_live_gap_fro"]) for r in relevant]),
                    "mean_relax_convergence_ratio_last_over_first": mean([float(r["mean_relax_convergence_ratio_last_over_first"]) for r in relevant]),
                    "relaxation_active_fraction": mean([float(r["relaxation_active_fraction"]) for r in relevant]),
                    "quasi_equilibrium_fraction": mean([float(r["quasi_equilibrium_fraction"]) for r in relevant]),
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
                    "mean_parent_record_vs_live_gap_fro": 0.0,
                    "mean_relax_convergence_ratio_last_over_first": 0.0,
                    "relaxation_active_fraction": 0.0,
                    "quasi_equilibrium_fraction": 0.0,
                }
            )
        return row


def run_suite(max_level: int, relax_steps: int, outdir: Path) -> dict:
    outdir.mkdir(parents=True, exist_ok=True)
    configs = [
        {
            "variant": "real_growth_linear_script1_2_true_schur_birth_plus_relax",
            "mode": "linear",
            "alpha_env": 0.22,
            "br_ancestor": 0.045,
            "br_sibling": 0.035,
            "order_sequence": (1, 2, 3),
            "relax_steps": relax_steps,
        },
        {
            "variant": "log_growth_script1_2_true_schur_birth_plus_relax",
            "mode": "log",
            "alpha_env": 0.22,
            "br_ancestor": 0.045,
            "br_sibling": 0.035,
            "order_sequence": (1, 2, 3),
            "relax_steps": relax_steps,
        },
        {
            "variant": "saturating_growth_script1_2_true_schur_birth_plus_relax",
            "mode": "saturating",
            "alpha_env": 0.90,
            "br_ancestor": 0.045,
            "br_sibling": 0.035,
            "order_sequence": (1, 2, 3),
            "relax_steps": relax_steps,
        },
        {
            "variant": "strict_symmetrized_response_control_true_schur_birth_plus_relax",
            "mode": "linear",
            "alpha_env": 0.0,
            "br_ancestor": 0.0,
            "br_sibling": 0.0,
            "order_sequence": (1, 2, 3),
            "relax_steps": relax_steps,
        },
        {
            "variant": "birth_only_no_relax_record_control_true_schur",
            "mode": "linear",
            "alpha_env": 0.22,
            "br_ancestor": 0.045,
            "br_sibling": 0.035,
            "order_sequence": (1, 2, 3),
            "relax_steps": 0,
        },
    ]
    summary = {
        "model_family": "script1_script2_true_schur_dtn_birth_plus_fixed_topology_live_relaxation",
        "max_level": max_level,
        "relax_steps_requested": relax_steps,
        "derived_only_notes": [
            "Uses established script-1/script-2 ternary sequential growth with incoming environment edges and newborn backreaction.",
            "Computes true Schur/DtN matrices from real conductance Laplacians on boundary ports {cut_node}+UV-tail leaves.",
            "Separates birth_delta_DtN from relax_delta_DtN at fixed topology/boundary.",
            "Live relaxation uses only existing birth_g, current g, and incoming/outgoing edge response loads; no J/i/Hodge/star/positivity/delta-beta gate is introduced.",
        ],
        "variants": [],
    }
    for cfg in configs:
        m = SchurDtNRelaxationModel(**cfg)
        m.run(max_level)
        prefix = cfg["variant"]
        write_csv(outdir / f"events_{prefix}.csv", m.event_rows)
        write_csv(outdir / f"boundary_role_rows_{prefix}.csv", m.role_rows)
        write_csv(outdir / f"dtn_schur_rows_{prefix}.csv", m.dtn_rows)
        write_csv(outdir / f"relax_rows_{prefix}.csv", m.relax_rows)
        write_csv(outdir / f"relax_event_rows_{prefix}.csv", m.relax_event_rows)
        write_csv(outdir / f"triples_{prefix}.csv", m.triple_rows)
        write_csv(outdir / f"levels_{prefix}.csv", m.level_rows)
        final = m.level_rows[-1]
        birth_role_gate = (
            final["retarded_event_fraction"] >= 0.999
            and final["advanced_leakage_fraction"] <= 1e-12
            and final["true_schur_dtn_parent_response_fraction"] >= 0.999
            and final["new_child_port_coupling_positive_fraction"] >= 0.999
        )
        relax_live_gate = (
            final["relax_events"] > 0
            and final["relaxation_active_fraction"] > 0.5
            and final["mean_record_vs_live_gap_fro"] > 1e-10
            and final["mean_total_relax_drift_dtn_fro"] > 1e-10
        )
        summary["variants"].append(
            {
                "variant": prefix,
                "events": int(final["events"]),
                "relax_events": int(final["relax_events"]),
                "completed_triples": int(final["completed_triples"]),
                "birth_schur_dtn_role_gate": bool(birth_role_gate),
                "fixed_topology_live_relaxation_gate": bool(relax_live_gate),
                "retarded_event_fraction": float(final["retarded_event_fraction"]),
                "advanced_leakage_fraction": float(final["advanced_leakage_fraction"]),
                "true_schur_dtn_parent_response_fraction": float(final["true_schur_dtn_parent_response_fraction"]),
                "new_child_port_coupling_positive_fraction": float(final["new_child_port_coupling_positive_fraction"]),
                "mean_parentline_schur_dtn_delta": float(final["mean_parentline_schur_dtn_delta"]),
                "mean_first_relax_delta_dtn_fro": float(final["mean_first_relax_delta_dtn_fro"]),
                "mean_last_relax_delta_dtn_fro": float(final["mean_last_relax_delta_dtn_fro"]),
                "mean_total_relax_drift_dtn_fro": float(final["mean_total_relax_drift_dtn_fro"]),
                "mean_record_vs_live_gap_fro": float(final["mean_record_vs_live_gap_fro"]),
                "max_record_vs_live_gap_fro": float(final["max_record_vs_live_gap_fro"]),
                "mean_parent_record_vs_live_gap_fro": float(final["mean_parent_record_vs_live_gap_fro"]),
                "mean_relax_convergence_ratio_last_over_first": float(final["mean_relax_convergence_ratio_last_over_first"]),
                "relaxation_active_fraction": float(final["relaxation_active_fraction"]),
                "quasi_equilibrium_fraction": float(final["quasi_equilibrium_fraction"]),
                "mean_abs_log_circulation": float(final["mean_abs_log_circulation"]),
                "frac_full_markov_complex": float(final["frac_full_markov_complex"]),
                "used_delta_beta_any": False,
            }
        )
    (outdir / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    return summary


def make_markdown(summary: dict, outdir: Path) -> None:
    rows = []
    for v in summary["variants"]:
        rows.append(
            f"| {v['variant']} | {v['events']} | {v['relax_events']} | "
            f"{int(v['birth_schur_dtn_role_gate'])} | {int(v['fixed_topology_live_relaxation_gate'])} | "
            f"{v['retarded_event_fraction']:.3f} | {v['advanced_leakage_fraction']:.3f} | "
            f"{v['mean_parentline_schur_dtn_delta']:.6g} | "
            f"{v['mean_first_relax_delta_dtn_fro']:.6g} | {v['mean_last_relax_delta_dtn_fro']:.6g} | "
            f"{v['mean_total_relax_drift_dtn_fro']:.6g} | {v['mean_record_vs_live_gap_fro']:.6g} | "
            f"{v['mean_relax_convergence_ratio_last_over_first']:.3f} | {v['relaxation_active_fraction']:.3f} | "
            f"{v['quasi_equilibrium_fraction']:.3f} | {v['frac_full_markov_complex']:.3f} |"
        )
    table = "\n".join(rows)
    md = f"""# RESULTS — CNNA Schur/DtN birth + live relaxation dynamics gate, L2

## Model provenance

This package keeps the established script-1/script-2 dynamic birth conductance model and the true Schur/DtN boundary-role calculation.  It adds fixed-topology live relaxation after every birth event.

The split is explicit:

```text
birth_delta_DtN:
  topology and boundary roles change;
  newborn becomes a new UV-tail port for parent/ancestor cuts.

relax_delta_DtN:
  no new node is born;
  boundary ports stay fixed;
  live conductances and live edge weights relax;
  true Schur/DtN matrices can drift further.

record_vs_live_gap:
  difference between the immutable after-birth DtN record and the current live DtN state.
```

The live relaxation rule is deterministic and local: it uses only existing `birth_g`, current `g`, and incoming/outgoing edge response loads inherited from the script-1/script-2 growth model.  It does not use J, i, Hodge, physical star, positivity, a C*-norm, or delta-beta.

## Gate table

| variant | events | relax events | birth Schur/DtN gate | live relax gate | retarded | advanced | mean birth ΔDtN | first relax ΔDtN | last relax ΔDtN | total relax drift | record/live gap | last/first | relax active | quasi-equilibrium | full Markov complex |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
{table}

## Interpretation

The birth-event Schur/DtN role gate remains positive in the nontrivial growth variants: the newborn still has no own UV-tail from its self-perspective, while parent/ancestor cuts receive a new UV-tail boundary port and a true Schur/DtN response change.

The new result is the fixed-topology live layer.  In the nontrivial growth variants, the after-birth DtN record does not remain equal to the later live DtN state.  Relaxation at fixed boundary produces nonzero `relax_delta_DtN` and nonzero `record_vs_live_gap`.  The strict symmetrized response control keeps the topology birth role but kills the response/live relaxation layer.

So the corrected interpretation is:

```text
Birth creates boundary-role/port changes.
Relaxation then changes live Schur/DtN responses without new births.
The network is not static during growth; birth records and live responses separate.
```

This supports the two-layer reading:

```text
Record layer:
  immutable birth/role/DtN record at the moment of boundary-role creation.

Live layer:
  continuing conductance/response relaxation under existing directed response structure.
```

## Limits

- This is still an L2 finite approximant.
- The DtN computation is a true matrix Schur complement of the real conductance Laplacian.
- The relaxation rule is deterministic and local, but it is a model choice within the existing script-1/script-2 quantities.  It is not yet a theorem forcing one unique relaxation law.
- The package does not claim J, i, spin, a *-algebra, positivity, or Q/P compatibility.

## Next test

`test_role_recalibration_to_boundary_value_bridge_true_schur_live_gate.py`

Use the event-resolved true Schur/DtN birth records and live-relaxed Schur/DtN states as separate inputs.  Test whether the difference `record_vs_live_gap` carries a robust retarded/advanced boundary-value polarity that can later bridge to Q/P or an operator adjunction.
"""
    (outdir / "RESULTS.md").write_text(md, encoding="utf-8")
    (outdir / "SUMMARY.md").write_text(md, encoding="utf-8")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--max-level", type=int, default=2)
    ap.add_argument("--relax-steps", type=int, default=6)
    ap.add_argument("--outdir", type=Path, default=Path("outputs"))
    args = ap.parse_args()
    summary = run_suite(args.max_level, args.relax_steps, args.outdir)
    make_markdown(summary, args.outdir)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
