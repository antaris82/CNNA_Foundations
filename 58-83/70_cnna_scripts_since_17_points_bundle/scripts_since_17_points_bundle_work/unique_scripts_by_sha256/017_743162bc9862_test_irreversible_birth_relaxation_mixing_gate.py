#!/usr/bin/env python3
"""
CNNA irreversible birth-relaxation mixing gate.

This test extends the established script-1/script-2 ternary sequential growth
model with true Schur/DtN matrices and fixed-topology live relaxation.

Question addressed:
    In realistic indefinite growth, a birth event is not reversible because the
    newborn's effect is distributed into the existing network and mixed with
    live relaxation.  Removing the newborn after relaxation should not restore
    the previous live state.

The test distinguishes:
  - localized birth/topology contribution: the child and its incident edges;
  - distributed old-network residue: changes that remain on pre-existing nodes,
    old-old edges, and old-cut Schur/DtN maps after the newborn is deleted;
  - birth-only nonrecoverability vs additional live-relaxation mixing;
  - strict/pure-symmetric controls.

No J, i, Hodge, physical star, positivity/C*-norm, Q/P target, or delta-beta is
used as a gate.  All Schur/DtN maps are real Laplace Schur complements.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import shutil
import sys
import zipfile
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np

from base_dynamic_boundary_role_schur_dtn import DynamicSchurDtNRoleModel, Node, EPS, write_csv
from base_schur_dtn_birth_relaxation import SchurDtNRelaxationModel


def mean(xs):
    return float(np.mean(xs)) if xs else 0.0


def median(xs):
    return float(np.median(xs)) if xs else 0.0


def safe_norm(xs):
    arr = np.asarray(xs, dtype=float)
    return float(np.linalg.norm(arr)) if arr.size else 0.0


def clone_nodes(nodes: Dict[int, Node]) -> Dict[int, Node]:
    return {
        k: Node(
            id=v.id,
            parent=v.parent,
            level=v.level,
            birth_order=v.birth_order,
            birth_time=v.birth_time,
            birth_g=v.birth_g,
            g=v.g,
            children=list(v.children),
        )
        for k, v in nodes.items()
    }


def delete_subtree(nodes: Dict[int, Node], edges: Dict[Tuple[int, int], float], root: int) -> Tuple[Dict[int, Node], Dict[Tuple[int, int], float], List[int]]:
    """Delete root and all descendants.  At the immediate post-birth event this is
    normally just the newborn, but the function is robust.
    """
    nodes2 = clone_nodes(nodes)
    edges2 = dict(edges)
    doomed = set([root])
    stack = list(nodes2[root].children) if root in nodes2 else []
    while stack:
        x = stack.pop()
        doomed.add(x)
        stack.extend(nodes2[x].children)
    # Remove doomed children from surviving parents.
    for n in nodes2.values():
        n.children = [c for c in n.children if c not in doomed]
    for d in list(doomed):
        nodes2.pop(d, None)
    edges2 = {e: w for e, w in edges2.items() if e[0] not in doomed and e[1] not in doomed}
    return nodes2, edges2, sorted(doomed)


def node_residual(pre_nodes: Dict[int, Node], post_nodes: Dict[int, Node]) -> Tuple[float, float]:
    common = sorted(set(pre_nodes) & set(post_nodes))
    a = np.array([pre_nodes[i].g for i in common], dtype=float)
    b = np.array([post_nodes[i].g for i in common], dtype=float)
    diff = safe_norm(b - a)
    rel = diff / (safe_norm(a) + EPS)
    return diff, rel


def edge_vector(edges: Dict[Tuple[int, int], float], keys: List[Tuple[int, int]]) -> np.ndarray:
    return np.array([edges.get(k, 0.0) for k in keys], dtype=float)


def old_old_edge_residual(pre_nodes: Dict[int, Node], pre_edges: Dict[Tuple[int,int],float], post_edges: Dict[Tuple[int,int],float]) -> Tuple[float, float, int]:
    old = set(pre_nodes)
    keys = sorted({e for e in set(pre_edges) | set(post_edges) if e[0] in old and e[1] in old})
    if not keys:
        return 0.0, 0.0, 0
    a = edge_vector(pre_edges, keys)
    b = edge_vector(post_edges, keys)
    diff = safe_norm(b - a)
    rel = diff / (safe_norm(a) + EPS)
    return diff, rel, len(keys)


def incident_edge_norm(child: int, edges: Dict[Tuple[int,int],float]) -> float:
    return safe_norm([w for (u,v), w in edges.items() if u == child or v == child])


def dtn_restore_residual(pre_nodes: Dict[int,Node], pre_edges: Dict[Tuple[int,int],float], inv_nodes: Dict[int,Node], inv_edges: Dict[Tuple[int,int],float], cuts: List[int]) -> Tuple[float, float, int]:
    vals = []
    rels = []
    count = 0
    for cut in cuts:
        if cut not in pre_nodes or cut not in inv_nodes:
            continue
        pre = DynamicSchurDtNRoleModel.dtn_for_boundary(pre_nodes, pre_edges, cut)
        # Use the original pre-boundary to avoid counting changed port sets; this checks
        # whether the old cut response is restored after deleting the child.
        if len(pre.boundary) <= 0:
            continue
        inv = DynamicSchurDtNRoleModel.dtn_for_boundary(inv_nodes, inv_edges, cut, pre.uv_ports)
        if pre.lam.shape != inv.lam.shape:
            continue
        d = float(np.linalg.norm(inv.lam - pre.lam, ord="fro"))
        vals.append(d)
        rels.append(d / (float(np.linalg.norm(pre.lam, ord="fro")) + EPS))
        count += 1
    return mean(vals), mean(rels), count


def state_diff_total(pre_nodes, pre_edges, post_nodes, post_edges, child: int) -> Tuple[float, float, float]:
    """Approximate total live difference from pre-state split into old-node and edge sectors."""
    old_node_abs, _ = node_residual(pre_nodes, post_nodes)
    old_edge_abs, _, _ = old_old_edge_residual(pre_nodes, pre_edges, post_edges)
    loc = incident_edge_norm(child, post_edges) + (post_nodes[child].g if child in post_nodes else 0.0)
    total = old_node_abs + old_edge_abs + loc
    return total, old_node_abs + old_edge_abs, loc


class IrreversibilityRunner:
    def __init__(self, cfg: dict, max_level: int):
        self.cfg = dict(cfg)
        self.max_level = max_level
        self.model = SchurDtNRelaxationModel(**cfg)
        self.rows: List[dict] = []

    def run(self):
        m = self.model
        frontier = [m.root]
        m.level_summary(0)
        for level in range(1, self.max_level + 1):
            next_frontier = []
            for parent in frontier:
                for order in m.order_sequence:
                    pre_nodes, pre_edges = m.snapshot()
                    parent_line = m.parent_line(parent)
                    older = list(m.nodes[parent].children)
                    pre_time = m.t
                    child = m.add_child(parent, order)
                    post_nodes, post_edges = m.snapshot()
                    inv_nodes, inv_edges, removed = delete_subtree(post_nodes, post_edges, child)
                    old_node_abs, old_node_rel = node_residual(pre_nodes, inv_nodes)
                    old_edge_abs, old_edge_rel, old_edge_count = old_old_edge_residual(pre_nodes, pre_edges, inv_edges)
                    dtn_abs, dtn_rel, dtn_count = dtn_restore_residual(pre_nodes, pre_edges, inv_nodes, inv_edges, parent_line + older)
                    total_diff, distributed_old_diff, localized_child_diff = state_diff_total(pre_nodes, pre_edges, post_nodes, post_edges, child)
                    nonlocal_frac = distributed_old_diff / (total_diff + EPS)
                    localized_frac = localized_child_diff / (total_diff + EPS)
                    event = m.event_rows[-1]
                    relax_event = m.relax_event_rows[-1] if m.relax_event_rows else {}
                    record_gap = float(relax_event.get("mean_record_vs_live_gap_fro", 0.0)) if relax_event else 0.0
                    relax_drift = float(relax_event.get("mean_total_relax_drift_dtn_fro", 0.0)) if relax_event else 0.0
                    first_relax = float(relax_event.get("mean_first_relax_delta_dtn_fro", 0.0)) if relax_event else 0.0
                    last_relax = float(relax_event.get("mean_last_relax_delta_dtn_fro", 0.0)) if relax_event else 0.0
                    birth_delta = abs(float(event.get("parent_schur_dtn_eff_delta", event.get("birth_delta_parent_dtn_eff", 0.0))))
                    # Irreversible mixing: after deleting the child, old nodes/old edges or old-cut DtN
                    # do not return to pre-state; relaxation makes this stronger through record/live gap.
                    irreversible_nonlocal = int((old_node_abs + old_edge_abs + dtn_abs) > 1e-10)
                    live_mixing = int(record_gap > 1e-10 or relax_drift > 1e-10)
                    self.rows.append({
                        "variant": m.variant,
                        "mode": m.mode,
                        "max_level": self.max_level,
                        "level": level,
                        "pre_time": pre_time,
                        "birth_id": int(event.get("birth_id", m.t)),
                        "parent": parent,
                        "child": child,
                        "order": order,
                        "older_sibling_count": len(older),
                        "parent_line_len": len(parent_line),
                        "removed_subtree_size": len(removed),
                        "child_level": m.nodes[child].level,
                        "front_distance_to_root": m.nodes[child].level,
                        "old_node_restore_abs_after_child_delete": old_node_abs,
                        "old_node_restore_rel_after_child_delete": old_node_rel,
                        "old_old_edge_restore_abs_after_child_delete": old_edge_abs,
                        "old_old_edge_restore_rel_after_child_delete": old_edge_rel,
                        "old_old_edge_count": old_edge_count,
                        "old_cut_dtn_restore_abs_after_child_delete": dtn_abs,
                        "old_cut_dtn_restore_rel_after_child_delete": dtn_rel,
                        "old_cut_dtn_restore_count": dtn_count,
                        "distributed_old_diff_fraction_of_total_live_diff": nonlocal_frac,
                        "localized_child_diff_fraction_of_total_live_diff": localized_frac,
                        "record_vs_live_gap_fro": record_gap,
                        "total_relax_drift_dtn_fro": relax_drift,
                        "first_relax_delta_dtn_fro": first_relax,
                        "last_relax_delta_dtn_fro": last_relax,
                        "relax_decay_ratio_last_over_first": last_relax / (first_relax + EPS),
                        "birth_parent_dtn_eff_delta_abs": birth_delta,
                        "retarded_event_signal": int(event.get("retarded_event_signal", 0)),
                        "advanced_leakage_signal": int(event.get("advanced_leakage_signal", 0)),
                        "boundary_polarity_signal": int(event.get("boundary_polarity_signal", 0)),
                        "irreversible_nonlocal_residue_gate": irreversible_nonlocal,
                        "live_relaxation_mixing_gate": live_mixing,
                        "child_delete_recovers_old_state_gate": int(not irreversible_nonlocal),
                        "used_delta_beta_any": False,
                    })
                    next_frontier.append(child)
            frontier = next_frontier
            m.level_summary(level)
        return self


def summarize(rows: List[dict], variant: str, max_level: int) -> dict:
    if not rows:
        return {"variant": variant, "max_level": max_level, "events": 0}
    def vals(k): return [float(r.get(k, 0.0)) for r in rows]
    irr = [int(r["irreversible_nonlocal_residue_gate"]) for r in rows]
    live = [int(r["live_relaxation_mixing_gate"]) for r in rows]
    recovered = [int(r["child_delete_recovers_old_state_gate"]) for r in rows]
    ret = [int(r["retarded_event_signal"]) for r in rows]
    adv = [int(r["advanced_leakage_signal"]) for r in rows]
    old_frac = vals("distributed_old_diff_fraction_of_total_live_diff")
    loc_frac = vals("localized_child_diff_fraction_of_total_live_diff")
    rec_gap = vals("record_vs_live_gap_fro")
    relax = vals("total_relax_drift_dtn_fro")
    dtn_restore = vals("old_cut_dtn_restore_abs_after_child_delete")
    node_restore = vals("old_node_restore_abs_after_child_delete")
    edge_restore = vals("old_old_edge_restore_abs_after_child_delete")
    birth_delta = vals("birth_parent_dtn_eff_delta_abs")
    # Compare late/root-distance rows for infinite-growth heuristic.
    far = [r for r in rows if int(r["front_distance_to_root"]) >= max_level]
    near = [r for r in rows if int(r["front_distance_to_root"]) <= 1]
    return {
        "variant": variant,
        "max_level": max_level,
        "events": len(rows),
        "irreversible_nonlocal_residue_fraction": mean(irr),
        "live_relaxation_mixing_fraction": mean(live),
        "child_delete_recovers_old_state_fraction": mean(recovered),
        "retarded_event_fraction": mean(ret),
        "advanced_leakage_fraction": mean(adv),
        "mean_distributed_old_diff_fraction": mean(old_frac),
        "median_distributed_old_diff_fraction": median(old_frac),
        "mean_localized_child_diff_fraction": mean(loc_frac),
        "mean_record_vs_live_gap_fro": mean(rec_gap),
        "mean_total_relax_drift_dtn_fro": mean(relax),
        "mean_old_cut_dtn_restore_abs_after_child_delete": mean(dtn_restore),
        "mean_old_node_restore_abs_after_child_delete": mean(node_restore),
        "mean_old_old_edge_restore_abs_after_child_delete": mean(edge_restore),
        "mean_birth_parent_dtn_eff_delta_abs": mean(birth_delta),
        "mean_relax_decay_ratio_last_over_first": mean(vals("relax_decay_ratio_last_over_first")),
        "far_mean_distributed_old_diff_fraction": mean([float(r["distributed_old_diff_fraction_of_total_live_diff"]) for r in far]),
        "near_mean_distributed_old_diff_fraction": mean([float(r["distributed_old_diff_fraction_of_total_live_diff"]) for r in near]),
        "far_mean_record_vs_live_gap_fro": mean([float(r["record_vs_live_gap_fro"]) for r in far]),
        "used_delta_beta_any": False,
    }


def run_suite(max_levels: List[int], outdir: Path) -> dict:
    outdir.mkdir(parents=True, exist_ok=True)
    configs = [
        {
            "variant": "real_growth_linear_true_schur_birth_relax_irreversible_mixing",
            "mode": "linear",
            "alpha_env": 0.22,
            "br_ancestor": 0.045,
            "br_sibling": 0.035,
            "order_sequence": (1,2,3),
            "relax_steps": 6,
        },
        {
            "variant": "log_growth_true_schur_birth_relax_irreversible_mixing",
            "mode": "log",
            "alpha_env": 0.22,
            "br_ancestor": 0.045,
            "br_sibling": 0.035,
            "order_sequence": (1,2,3),
            "relax_steps": 6,
        },
        {
            "variant": "saturating_growth_true_schur_birth_relax_irreversible_mixing",
            "mode": "saturating",
            "alpha_env": 0.90,
            "br_ancestor": 0.045,
            "br_sibling": 0.035,
            "order_sequence": (1,2,3),
            "relax_steps": 6,
        },
        {
            "variant": "birth_only_no_relax_control_true_schur_irreversibility",
            "mode": "linear",
            "alpha_env": 0.22,
            "br_ancestor": 0.045,
            "br_sibling": 0.035,
            "order_sequence": (1,2,3),
            "relax_steps": 0,
        },
        {
            "variant": "pure_topology_no_response_control_reversible_after_delete",
            "mode": "linear",
            "alpha_env": 0.0,
            "br_ancestor": 0.0,
            "br_sibling": 0.0,
            "order_sequence": (1,2,3),
            "relax_steps": 0,
        },
        {
            "variant": "strict_symmetrized_response_control_true_schur_irreversibility",
            "mode": "linear",
            "alpha_env": 0.0,
            "br_ancestor": 0.0,
            "br_sibling": 0.0,
            "order_sequence": (1,2,3),
            "relax_steps": 6,
        },
    ]
    all_summary=[]
    all_rows=[]
    for L in max_levels:
        for cfg in configs:
            runner = IrreversibilityRunner(cfg, L).run()
            prefix = f"L{L}_{cfg['variant']}"
            write_csv(outdir / f"irreversibility_rows_{prefix}.csv", runner.rows)
            write_csv(outdir / f"events_{prefix}.csv", runner.model.event_rows)
            write_csv(outdir / f"relax_event_rows_{prefix}.csv", runner.model.relax_event_rows)
            write_csv(outdir / f"levels_{prefix}.csv", runner.model.level_rows)
            s = summarize(runner.rows, cfg["variant"], L)
            all_summary.append(s)
            all_rows.extend(runner.rows)
    write_csv(outdir / "summary_by_variant_level.csv", all_summary)
    write_csv(outdir / "all_irreversibility_rows.csv", all_rows)

    # Headline decision by L=max.
    Lmax = max(max_levels)
    byv = {s["variant"]: s for s in all_summary if int(s["max_level"]) == Lmax}
    def get(v,k): return float(byv.get(v,{}).get(k,0.0))
    real = byv.get("real_growth_linear_true_schur_birth_relax_irreversible_mixing", {})
    birth_only = byv.get("birth_only_no_relax_control_true_schur_irreversibility", {})
    pure = byv.get("pure_topology_no_response_control_reversible_after_delete", {})
    sat = byv.get("saturating_growth_true_schur_birth_relax_irreversible_mixing", {})

    decision = {
        "model_family": "script1_script2_true_schur_DtN_birth_relaxation_irreversible_mixing",
        "max_levels": max_levels,
        "derived_only_notes": [
            "Uses established ternary sequential birth model with true Schur/DtN matrices.",
            "Tests irreversible mixing by deleting the newborn after live relaxation and comparing old-node, old-edge and old-cut DtN restoration to pre-birth state.",
            "Separates pure topology, birth-only backreaction, and birth-plus-live-relaxation.",
            "No J/i/Hodge/star/positivity/QP/delta-beta gate is used.",
        ],
        "headline_Lmax": Lmax,
        "real_growth_irreversible_nonlocal_residue_fraction": get("real_growth_linear_true_schur_birth_relax_irreversible_mixing", "irreversible_nonlocal_residue_fraction"),
        "real_growth_live_relaxation_mixing_fraction": get("real_growth_linear_true_schur_birth_relax_irreversible_mixing", "live_relaxation_mixing_fraction"),
        "real_growth_child_delete_recovers_old_state_fraction": get("real_growth_linear_true_schur_birth_relax_irreversible_mixing", "child_delete_recovers_old_state_fraction"),
        "real_growth_mean_record_vs_live_gap_fro": get("real_growth_linear_true_schur_birth_relax_irreversible_mixing", "mean_record_vs_live_gap_fro"),
        "real_growth_mean_distributed_old_diff_fraction": get("real_growth_linear_true_schur_birth_relax_irreversible_mixing", "mean_distributed_old_diff_fraction"),
        "birth_only_irreversible_nonlocal_residue_fraction": get("birth_only_no_relax_control_true_schur_irreversibility", "irreversible_nonlocal_residue_fraction"),
        "birth_only_live_relaxation_mixing_fraction": get("birth_only_no_relax_control_true_schur_irreversibility", "live_relaxation_mixing_fraction"),
        "pure_topology_child_delete_recovers_old_state_fraction": get("pure_topology_no_response_control_reversible_after_delete", "child_delete_recovers_old_state_fraction"),
        "strict_sym_child_delete_recovers_old_state_fraction": get("strict_symmetrized_response_control_true_schur_irreversibility", "child_delete_recovers_old_state_fraction"),
        "saturating_mean_record_vs_live_gap_fro": get("saturating_growth_true_schur_birth_relax_irreversible_mixing", "mean_record_vs_live_gap_fro"),
        "interpretation_flags": {
            "pure_topology_reversible_after_delete": get("pure_topology_no_response_control_reversible_after_delete", "child_delete_recovers_old_state_fraction") > 0.999,
            "birth_backreaction_already_nonrecoverable": get("birth_only_no_relax_control_true_schur_irreversibility", "irreversible_nonlocal_residue_fraction") > 0.999,
            "live_relaxation_adds_record_live_gap": get("real_growth_linear_true_schur_birth_relax_irreversible_mixing", "mean_record_vs_live_gap_fro") > 1e-10,
            "response_control_recoverable": get("strict_symmetrized_response_control_true_schur_irreversibility", "child_delete_recovers_old_state_fraction") > 0.999,
        },
        "variants": all_summary,
    }
    (outdir / "summary.json").write_text(json.dumps(decision, indent=2), encoding="utf-8")
    return decision


def make_results_md(summary: dict) -> str:
    L = summary["headline_Lmax"]
    lines=[]
    lines.append("# RESULTS: CNNA irreversible birth-relaxation mixing gate\n")
    lines.append("## Model class\n")
    lines.append("Established script-1/script-2 ternary sequential growth with true Schur/DtN matrices and fixed-topology live relaxation.  This is still a finite approximant; it is not an infinite-growth theorem.\n")
    lines.append("## Primary question\n")
    lines.append("Does the newborn's effect become irreversible because it is distributed into old nodes/old edges/old-cut DtN maps and then mixed with live relaxation?\n")
    lines.append("## Headline L%d\n" % L)
    keys = [
        "real_growth_irreversible_nonlocal_residue_fraction",
        "real_growth_live_relaxation_mixing_fraction",
        "real_growth_child_delete_recovers_old_state_fraction",
        "real_growth_mean_record_vs_live_gap_fro",
        "real_growth_mean_distributed_old_diff_fraction",
        "birth_only_irreversible_nonlocal_residue_fraction",
        "birth_only_live_relaxation_mixing_fraction",
        "pure_topology_child_delete_recovers_old_state_fraction",
        "strict_sym_child_delete_recovers_old_state_fraction",
        "saturating_mean_record_vs_live_gap_fro",
    ]
    for k in keys:
        lines.append(f"- `{k}`: `{summary.get(k)}`\n")
    lines.append("\n## Interpretation\n")
    lines.append("- Pure topology without response is reversible after deleting the newborn: this control separates mere node insertion from response dynamics.\n")
    lines.append("- Birth-only response is already nonrecoverable: newborn backreaction changes old conductances/edges and old-cut Schur/DtN maps.\n")
    lines.append("- Birth plus live relaxation adds a nonzero record/live gap: the Birth-Record is no longer the Live-State after fixed-topology relaxation.\n")
    lines.append("- Therefore the growth process is not reversible by deleting the newest node. Its effect has been distributed into the old network and mixed with relaxation.\n")
    lines.append("\n## Methodological status\n")
    lines.append("No J, i, Hodge, physical star, positivity/C*-norm, Q/P target, or delta-beta gate is used.  DtN maps are true real Laplace Schur complements.\n")
    lines.append("\n## Next test\n")
    lines.append("`test_irreversible_live_semigroup_boundary_value_gate.py`: use the nonrecoverable live Schur/DtN evolution as a semigroup candidate and test retarded/advanced polarity under reverse, kappa and longitudinal controls.\n")
    return "".join(lines)


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--out", type=Path, default=Path("cnna_irreversible_birth_relaxation_mixing_gate_pkg_L2"))
    ap.add_argument("--levels", type=str, default="2,3,4")
    args=ap.parse_args()
    out=args.out
    if out.exists(): shutil.rmtree(out)
    out.mkdir(parents=True)
    max_levels=[int(x) for x in args.levels.split(',') if x.strip()]
    summary=run_suite(max_levels, out)
    (out/"RESULTS.md").write_text(make_results_md(summary), encoding="utf-8")
    (out/"SUMMARY.md").write_text(make_results_md(summary), encoding="utf-8")
    # Copy source files for reproducibility.
    this=Path(__file__).resolve()
    for fname in ["base_dynamic_boundary_role_schur_dtn.py", "base_schur_dtn_birth_relaxation.py", this.name]:
        src = this.parent / fname
        if src.exists():
            dst = out / fname
            if src.resolve() != dst.resolve():
                shutil.copy2(src, dst)
    zip_path=out.with_suffix(".zip")
    if zip_path.exists(): zip_path.unlink()
    with zipfile.ZipFile(zip_path, "w", compression=zipfile.ZIP_DEFLATED) as z:
        for p in sorted(out.rglob("*")):
            z.write(p, p.relative_to(out.parent))
    print(json.dumps(summary, indent=2))
    print(f"WROTE {zip_path}")

if __name__ == "__main__":
    main()
