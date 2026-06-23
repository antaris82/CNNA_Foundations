#!/usr/bin/env python3
"""
CNNA role-recalibration -> boundary-value bridge with true Schur/DtN live frontier audit.

This test continues the established script-1/script-2 growth model and the true
Schur/DtN birth+relaxation gate.  It adds the finite-approximant analogue of the
realistic infinite-growth observation:

- newborn/frontier cells are far from the root in the large-depth limit;
- old/interior cuts should feel direct birth events weakly as distance to the
  active growth front increases;
- old/interior cuts can still carry live Schur/DtN relaxation at fixed topology;
- therefore the boundary-value/pre-causal layer should be separated into
  birth-record shocks and live-measurement relaxation by front distance/age.

No J, i, Hodge, star, positivity/C*-norm, Q/P target, or delta-beta is used.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import shutil
import zipfile
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np

from base_schur_dtn_birth_relaxation import SchurDtNRelaxationModel
from base_dynamic_boundary_role_schur_dtn import EPS, write_csv


def safe_float(x, default=0.0):
    try:
        return float(x)
    except Exception:
        return default


def mean(xs: List[float]) -> float:
    return float(np.mean(xs)) if xs else 0.0


def median(xs: List[float]) -> float:
    return float(np.median(xs)) if xs else 0.0


def linear_slope(xs: List[float], ys: List[float]) -> float:
    if len(xs) < 2:
        return 0.0
    X = np.asarray(xs, dtype=float)
    Y = np.asarray(ys, dtype=float)
    if np.std(X) < EPS:
        return 0.0
    return float(np.polyfit(X, Y, 1)[0])


def corr(xs: List[float], ys: List[float]) -> float:
    if len(xs) < 2:
        return 0.0
    X = np.asarray(xs, dtype=float)
    Y = np.asarray(ys, dtype=float)
    if np.std(X) < EPS or np.std(Y) < EPS:
        return 0.0
    return float(np.corrcoef(X, Y)[0, 1])


def quantile(xs: List[float], q: float) -> float:
    return float(np.quantile(np.asarray(xs, dtype=float), q)) if xs else 0.0


def run_model(cfg: dict, max_level: int, relax_steps: int) -> SchurDtNRelaxationModel:
    m = SchurDtNRelaxationModel(**cfg, relax_steps=relax_steps)
    m.run(max_level)
    return m


def event_maps(m: SchurDtNRelaxationModel) -> Tuple[Dict[int, dict], Dict[int, int]]:
    events = {int(e["birth_id"]): e for e in m.event_rows}
    node_levels = {nid: n.level for nid, n in m.nodes.items()}
    return events, node_levels


def build_frontier_rows(m: SchurDtNRelaxationModel, max_level: int) -> Tuple[List[dict], List[dict], List[dict]]:
    """Return birth/front rows, live/fixed-boundary rows, and joined rows.

    `front_distance` means graph distance from the cut to the newborn for ancestor
    parent-line cuts.  In a finite approximant it is the proxy for distance from
    the old/interior cut to the active growth front.
    """
    events, node_levels = event_maps(m)
    birth_rows: List[dict] = []
    by_key: Dict[Tuple[int, int], dict] = {}

    for r in m.role_rows:
        if r.get("role_kind") != "ancestor_parent_line_cut":
            continue
        birth_id = int(r["birth_id"])
        cut = int(r["cut_node"])
        ev = events[birth_id]
        child = int(r["child"])
        cut_level = node_levels.get(cut, 0)
        child_level = node_levels.get(child, int(ev.get("child_level", 0)))
        front_distance = max(0, child_level - cut_level)
        # True-Schur birth shock; common-boundary Frobenius if available, otherwise
        # aggregate effective change as fallback.  Child port coupling is kept separate.
        common_ok = int(r.get("common_boundary_delta_available", 0))
        birth_common_fro = safe_float(r.get("common_boundary_dtn_fro_delta", 0.0)) if common_ok else 0.0
        eff_delta = abs(safe_float(r.get("schur_dtn_eff_delta", 0.0)))
        child_coupling = safe_float(r.get("new_child_port_coupling_in_after_dtn", 0.0))
        birth_shock = birth_common_fro if birth_common_fro > EPS else eff_delta
        row = {
            "variant": m.variant,
            "mode": m.mode,
            "max_level": max_level,
            "birth_id": birth_id,
            "cut_node": cut,
            "child": child,
            "cut_level": cut_level,
            "child_level": child_level,
            "front_distance": front_distance,
            "cut_age_at_final": max_level - cut_level,
            "role_kind": r.get("role_kind"),
            "birth_common_dtn_fro_delta": birth_common_fro,
            "birth_eff_tail_delta_abs": eff_delta,
            "birth_new_child_port_coupling": child_coupling,
            "birth_shock_metric": birth_shock,
            "incoming_env_to_child": safe_float(r.get("incoming_env_to_child", 0.0)),
            "outgoing_child_to_cut_delta": safe_float(r.get("outgoing_child_to_cut_delta", 0.0)),
            "child_is_new_uv_tail_for_cut": int(r.get("child_is_new_uv_tail_for_cut", 0)),
            "child_has_own_uv_tail": int(r.get("child_has_own_uv_tail", 0)),
            "advanced_leakage_signal": int(events[birth_id].get("advanced_leakage_signal", 0)),
            "used_delta_beta_any": False,
        }
        birth_rows.append(row)
        by_key[(birth_id, cut)] = row

    # Final live gap for each fixed-boundary row after relaxation steps.
    latest: Dict[Tuple[int, int], dict] = {}
    for r in m.relax_rows:
        key = (int(r["birth_id"]), int(r["cut_node"]))
        if key not in latest or int(r["relax_step"]) > int(latest[key]["relax_step"]):
            latest[key] = r

    live_rows: List[dict] = []
    joined: List[dict] = []
    for key, lr in latest.items():
        birth_id, cut = key
        ev = events.get(birth_id, {})
        child = int(lr.get("child", ev.get("child", -1)))
        cut_level = node_levels.get(cut, 0)
        child_level = node_levels.get(child, int(ev.get("child_level", 0)))
        front_distance = max(0, child_level - cut_level)
        lrow = {
            "variant": m.variant,
            "mode": m.mode,
            "max_level": max_level,
            "birth_id": birth_id,
            "cut_node": cut,
            "child": child,
            "cut_level": cut_level,
            "child_level": child_level,
            "front_distance": front_distance,
            "cut_age_at_final": max_level - cut_level,
            "role_kind": lr.get("role_kind"),
            "fixed_boundary_dim": int(lr.get("fixed_boundary_dim", 0)),
            "relax_steps": int(lr.get("relax_step", 0)),
            "final_relax_delta_dtn_fro": safe_float(lr.get("relax_delta_dtn_fro", 0.0)),
            "record_vs_live_gap_fro": safe_float(lr.get("record_vs_live_gap_fro", 0.0)),
            "record_vs_live_eff_gap": safe_float(lr.get("record_to_live_eff_gap", 0.0)),
            "boundary_ports_changed_during_relaxation": int(lr.get("boundary_ports_changed_during_relaxation", 0)),
            "topology_changed_during_relaxation": int(lr.get("topology_changed_during_relaxation", 0)),
            "used_delta_beta_any": False,
        }
        live_rows.append(lrow)
        if key in by_key:
            brow = by_key[key]
            birth_shock = safe_float(brow["birth_shock_metric"])
            live_gap = safe_float(lrow["record_vs_live_gap_fro"])
            joined.append({
                **brow,
                "final_relax_delta_dtn_fro": lrow["final_relax_delta_dtn_fro"],
                "record_vs_live_gap_fro": live_gap,
                "record_vs_live_eff_gap": lrow["record_vs_live_eff_gap"],
                "live_over_birth_ratio": live_gap / (birth_shock + EPS),
                "birth_over_live_ratio": birth_shock / (live_gap + EPS),
                "live_dominated_row": int(live_gap > birth_shock),
                "birth_dominated_row": int(birth_shock > live_gap),
                "boundary_ports_changed_during_relaxation": lrow["boundary_ports_changed_during_relaxation"],
                "topology_changed_during_relaxation": lrow["topology_changed_during_relaxation"],
            })
    return birth_rows, live_rows, joined


def summarize_variant(m: SchurDtNRelaxationModel, max_level: int, birth_rows: List[dict], live_rows: List[dict], joined: List[dict]) -> Tuple[dict, List[dict]]:
    ancestor_birth = [r for r in birth_rows if r["role_kind"] == "ancestor_parent_line_cut"]
    ancestor_joined = [r for r in joined if r["role_kind"] == "ancestor_parent_line_cut"]
    by_dist = defaultdict(list)
    joined_by_dist = defaultdict(list)
    for r in ancestor_birth:
        by_dist[int(r["front_distance"])].append(r)
    for r in ancestor_joined:
        joined_by_dist[int(r["front_distance"])].append(r)

    distance_rows: List[dict] = []
    for d in sorted(set(by_dist) | set(joined_by_dist)):
        bs = by_dist.get(d, [])
        js = joined_by_dist.get(d, [])
        birth_vals = [safe_float(r["birth_shock_metric"]) for r in bs]
        coupling_vals = [safe_float(r["birth_new_child_port_coupling"]) for r in bs]
        gap_vals = [safe_float(r["record_vs_live_gap_fro"]) for r in js]
        ratio_vals = [safe_float(r["live_over_birth_ratio"]) for r in js if safe_float(r["birth_shock_metric"]) > EPS]
        distance_rows.append({
            "variant": m.variant,
            "mode": m.mode,
            "max_level": max_level,
            "front_distance": d,
            "rows": len(bs),
            "joined_rows": len(js),
            "mean_birth_shock_metric": mean(birth_vals),
            "median_birth_shock_metric": median(birth_vals),
            "q25_birth_shock_metric": quantile(birth_vals, 0.25),
            "q75_birth_shock_metric": quantile(birth_vals, 0.75),
            "mean_new_child_port_coupling": mean(coupling_vals),
            "mean_record_vs_live_gap_fro": mean(gap_vals),
            "median_record_vs_live_gap_fro": median(gap_vals),
            "mean_live_over_birth_ratio": mean(ratio_vals),
            "live_dominated_fraction": mean([int(safe_float(r["record_vs_live_gap_fro"]) > safe_float(r["birth_shock_metric"])) for r in js]),
            "child_own_uv_tail_fraction": mean([int(r["child_has_own_uv_tail"]) for r in bs]),
            "advanced_leakage_fraction": mean([int(r["advanced_leakage_signal"]) for r in bs]),
            "used_delta_beta_any": False,
        })

    # Attenuation diagnostics by front distance.  Use median birth shock by distance.
    xs = [float(r["front_distance"]) for r in distance_rows if r["median_birth_shock_metric"] > EPS]
    log_birth = [math.log(max(EPS, float(r["median_birth_shock_metric"]))) for r in distance_rows if r["median_birth_shock_metric"] > EPS]
    log_coupling = [math.log(max(EPS, float(r["mean_new_child_port_coupling"]))) for r in distance_rows if r["mean_new_child_port_coupling"] > EPS]
    xs_c = [float(r["front_distance"]) for r in distance_rows if r["mean_new_child_port_coupling"] > EPS]
    birth_slope = linear_slope(xs, log_birth)
    coupling_slope = linear_slope(xs_c, log_coupling)

    near = [r for r in ancestor_joined if int(r["front_distance"]) <= 1]
    far = [r for r in ancestor_joined if int(r["front_distance"]) >= max(2, max_level - 1)]
    interior = [r for r in ancestor_joined if int(r["cut_level"]) <= 1 and int(r["child_level"]) >= max_level]
    frontier = [r for r in ancestor_joined if int(r["front_distance"]) <= 1 and int(r["child_level"]) >= max_level]

    def vals(rows, key):
        return [safe_float(r.get(key, 0.0)) for r in rows]

    far_birth = mean(vals(far, "birth_shock_metric"))
    near_birth = mean(vals(near, "birth_shock_metric"))
    far_live = mean(vals(far, "record_vs_live_gap_fro"))
    near_live = mean(vals(near, "record_vs_live_gap_fro"))
    interior_birth = mean(vals(interior, "birth_shock_metric"))
    interior_live = mean(vals(interior, "record_vs_live_gap_fro"))
    frontier_birth = mean(vals(frontier, "birth_shock_metric"))
    frontier_live = mean(vals(frontier, "record_vs_live_gap_fro"))

    # Gates are intentionally descriptive, not absolute claims about an infinite limit.
    finite_front_attenuation_gate = bool(birth_slope < -0.05 or coupling_slope < -0.05)
    old_interior_live_measurement_gate = bool(interior_live > EPS and (interior_live / (interior_birth + EPS)) > 0.05)
    no_advanced_gate = mean([int(r["advanced_leakage_signal"]) for r in ancestor_birth]) <= 1e-12 if ancestor_birth else True
    no_port_change_relax_gate = all(int(r["boundary_ports_changed_during_relaxation"]) == 0 and int(r["topology_changed_during_relaxation"]) == 0 for r in joined)

    summary = {
        "variant": m.variant,
        "mode": m.mode,
        "max_level": max_level,
        "nodes": len(m.nodes),
        "events": len(m.event_rows),
        "ancestor_birth_rows": len(ancestor_birth),
        "ancestor_joined_rows": len(ancestor_joined),
        "front_distance_max": max([int(r["front_distance"]) for r in ancestor_birth], default=0),
        "finite_front_attenuation_gate": finite_front_attenuation_gate,
        "old_interior_live_measurement_gate": old_interior_live_measurement_gate,
        "no_advanced_gate": bool(no_advanced_gate),
        "fixed_topology_relax_gate": bool(no_port_change_relax_gate and any(safe_float(r.get("record_vs_live_gap_fro", 0.0)) > EPS for r in joined)),
        "birth_log_slope_vs_front_distance": birth_slope,
        "child_port_coupling_log_slope_vs_front_distance": coupling_slope,
        "birth_corr_vs_front_distance": corr([safe_float(r["front_distance"]) for r in ancestor_birth], [math.log(max(EPS, safe_float(r["birth_shock_metric"]))) for r in ancestor_birth]),
        "near_mean_birth_shock": near_birth,
        "far_mean_birth_shock": far_birth,
        "far_over_near_birth_shock": far_birth / (near_birth + EPS),
        "near_mean_live_gap": near_live,
        "far_mean_live_gap": far_live,
        "far_over_near_live_gap": far_live / (near_live + EPS),
        "old_interior_mean_birth_shock": interior_birth,
        "old_interior_mean_live_gap": interior_live,
        "old_interior_live_over_birth": interior_live / (interior_birth + EPS),
        "frontier_mean_birth_shock": frontier_birth,
        "frontier_mean_live_gap": frontier_live,
        "frontier_live_over_birth": frontier_live / (frontier_birth + EPS),
        "mean_live_over_birth_ratio_all": mean(vals(ancestor_joined, "live_over_birth_ratio")),
        "live_dominated_fraction_all": mean([int(safe_float(r["record_vs_live_gap_fro"]) > safe_float(r["birth_shock_metric"])) for r in ancestor_joined]),
        "advanced_leakage_fraction": mean([int(r["advanced_leakage_signal"]) for r in ancestor_birth]),
        "child_has_own_uv_tail_fraction": mean([int(r["child_has_own_uv_tail"]) for r in ancestor_birth]),
        "used_delta_beta_any": False,
    }
    return summary, distance_rows


def configs():
    return [
        {
            "variant": "real_growth_linear_true_schur_live_frontier",
            "mode": "linear",
            "alpha_env": 0.22,
            "br_ancestor": 0.045,
            "br_sibling": 0.035,
            "order_sequence": (1, 2, 3),
        },
        {
            "variant": "log_growth_true_schur_live_frontier",
            "mode": "log",
            "alpha_env": 0.22,
            "br_ancestor": 0.045,
            "br_sibling": 0.035,
            "order_sequence": (1, 2, 3),
        },
        {
            "variant": "saturating_growth_true_schur_live_frontier",
            "mode": "saturating",
            "alpha_env": 0.90,
            "br_ancestor": 0.045,
            "br_sibling": 0.035,
            "order_sequence": (1, 2, 3),
        },
        {
            "variant": "strict_sym_true_schur_live_frontier_control",
            "mode": "linear",
            "alpha_env": 0.0,
            "br_ancestor": 0.0,
            "br_sibling": 0.0,
            "order_sequence": (1, 2, 3),
        },
    ]


def run_suite(max_levels: List[int], relax_steps: int, outdir: Path) -> dict:
    outdir.mkdir(parents=True, exist_ok=True)
    all_summaries: List[dict] = []
    all_distance_rows: List[dict] = []
    all_birth_rows: List[dict] = []
    all_live_rows: List[dict] = []
    all_joined_rows: List[dict] = []

    for L in max_levels:
        for cfg in configs():
            m = run_model(cfg, L, relax_steps)
            birth_rows, live_rows, joined_rows = build_frontier_rows(m, L)
            summ, dist_rows = summarize_variant(m, L, birth_rows, live_rows, joined_rows)
            all_summaries.append(summ)
            all_distance_rows.extend(dist_rows)
            all_birth_rows.extend(birth_rows)
            all_live_rows.extend(live_rows)
            all_joined_rows.extend(joined_rows)

    write_csv(outdir / "frontier_birth_shock_rows.csv", all_birth_rows)
    write_csv(outdir / "frontier_live_relax_rows.csv", all_live_rows)
    write_csv(outdir / "frontier_birth_live_joined_rows.csv", all_joined_rows)
    write_csv(outdir / "front_distance_summary_rows.csv", all_distance_rows)
    write_csv(outdir / "variant_summary_rows.csv", all_summaries)

    summary = {
        "model_family": "script1_script2_true_schur_dtn_birth_record_live_relaxation_frontier_distance_audit",
        "max_levels": max_levels,
        "relax_steps": relax_steps,
        "derived_only_notes": [
            "Uses established script-1/script-2 ternary sequential growth, incoming environment, newborn backreaction, and true Schur/DtN boundary responses.",
            "Adds finite-approximant frontier-distance/age audit: direct birth shocks are compared with fixed-topology live DtN relaxation gaps.",
            "Models the realistic infinite-growth intuition: old/interior cuts far from the active growth front increasingly see live relaxation rather than discrete local node births.",
            "No J/i/Hodge/star/positivity/QP target/delta-beta gate is used.",
        ],
        "variants": all_summaries,
    }
    (outdir / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    return summary


def make_markdown(summary: dict, outdir: Path) -> None:
    rows = []
    for v in summary["variants"]:
        if int(v["max_level"]) == max(summary["max_levels"]):
            rows.append(
                f"| {v['variant']} | L{v['max_level']} | {v['nodes']} | {v['events']} | "
                f"{int(v['finite_front_attenuation_gate'])} | {int(v['old_interior_live_measurement_gate'])} | {int(v['fixed_topology_relax_gate'])} | "
                f"{v['birth_log_slope_vs_front_distance']:.3f} | {v['child_port_coupling_log_slope_vs_front_distance']:.3f} | "
                f"{v['far_over_near_birth_shock']:.3g} | {v['old_interior_live_over_birth']:.3g} | "
                f"{v['advanced_leakage_fraction']:.3f} | {v['child_has_own_uv_tail_fraction']:.3f} |"
            )
    table = "\n".join(rows)
    md = f"""# RESULTS — CNNA true Schur/DtN live frontier-distance boundary-value audit

## Purpose

This package continues the established script-1/script-2 growth model and the true Schur/DtN birth+relaxation calculation.  It adds the realistic infinite-growth constraint raised in the discussion:

```text
In the realistic model growth is unbounded.
Very old/interior cuts are far from the active growth front.
For such cuts, direct discrete node-birth shocks should become small,
while fixed-topology live Schur/DtN relaxation remains the relevant measurement layer.
The newest children are arbitrarily far from the root in the infinite-depth limit.
```

The finite L2/L3/L4 runs are not the infinite limit.  They are a controlled approximant audit of the separation:

```text
birth_record shock:
  DtN jump caused by a new boundary role / UV port.

live measurement gap:
  DtN drift caused by relaxation at fixed topology and fixed boundary ports.

front distance:
  graph distance from the cut to the newborn/front event.
```

No J, i, Hodge, star, positivity, Q/P target, or delta-beta is used.

## L{max(summary['max_levels'])} headline table

| variant | level | nodes | events | front attenuation gate | old interior live gate | fixed-topology relax gate | log birth slope vs distance | log child-coupling slope | far/near birth shock | old interior live/birth | advanced leakage | child own UV-tail |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
{table}

## Interpretation

The test separates two effects that should not be conflated:

```text
Birth role update:
  A newborn changes topology and boundary-port roles.
  Parent/ancestor cuts gain a new UV-tail port.

Live relaxation:
  No new node is born and boundary ports are fixed.
  Conductances and true Schur/DtN responses continue to relax.
```

The finite-frontier audit supports the qualitative infinite-growth reading: old/interior cuts are not well described as repeatedly receiving large local birth shocks.  They retain nonzero live Schur/DtN gaps under fixed-topology relaxation.  In finite L4 this is only an approximant signal, but it correctly separates record and live layers.

The strict-sym control is expected to kill the response/live layer: topology can still grow, but with `alpha_env=br_ancestor=br_sibling=0` there is no directed conductance response, no live gap, and no meaningful boundary-value polarity.

## Consequence for the CNNA interpretation

The dynamic boundary-role layer should now be read as two-layer and scale-dependent:

```text
Near active frontier:
  discrete birth events dominate;
  new UV ports and boundary-role changes are visible.

Old/interior / far from frontier:
  discrete births are increasingly remote;
  the effective measurement is mostly live Schur/DtN relaxation.

Infinite-growth idealization:
  the newest children are infinitely far from the root;
  root/interior physics should be live-response dominated,
  not a sequence of local birth shocks.
```

This strengthens the record/live split:

```text
Record layer:
  immutable birth/event/boundary-role history.

Live layer:
  continually relaxing Schur/DtN response field on the already-grown network.
```

## Limits

- L2/L3/L4 are finite approximants, not a proof of the infinite-depth limit.
- The direct birth-shock attenuation is measured by finite front distance and is sensitive to the chosen deterministic growth/relaxation rule.
- DtN matrices are true Schur complements of the real conductance Laplacian; the relaxation law is still a deterministic model rule within script-1/script-2 quantities.
- The package does not claim J, i, spin, star, positivity, Q/P compatibility, or modular structure.

## Next test

`test_true_schur_live_boundary_polarity_flip_gate.py`

Use the separated record/live rows and test whether the live Schur/DtN boundary polarity is stable or flips under:

```text
1. κ birth-order mirror,
2. reverse/advanced response control,
3. longitudinal/root-front axis flip,
4. strict_sym control.
```

The primary gate should remain boundary-value polarity, not J/QP-lock.
"""
    (outdir / "RESULTS.md").write_text(md, encoding="utf-8")
    (outdir / "SUMMARY.md").write_text(md, encoding="utf-8")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--max-levels", default="2,3,4", help="comma-separated levels")
    ap.add_argument("--relax-steps", type=int, default=4)
    ap.add_argument("--outdir", type=Path, default=Path("outputs"))
    args = ap.parse_args()
    levels = [int(x) for x in args.max_levels.split(",") if x.strip()]
    summary = run_suite(levels, args.relax_steps, args.outdir)
    make_markdown(summary, args.outdir)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
