"""
CNNA / growing real complement network
Activity-gradient + patch-size scaling diagnostic.

Purpose
-------
This test adds an explicit growth-activity audit to the patch-size / level
star-closure gate.

Hypothesis under test
---------------------
Most live response change should occur where nodes are actively growing or have
just completed; the root/old interior should change least.  If a local *-like
operator system is a limiting/local-net phenomenon, closure should be sought in
the correct active/local patches, not in a single old/root pair.

No J, no i, no complex phase, and no C*-norm are inserted.
"""
from __future__ import annotations

import argparse
import csv
import json
import time
import zipfile
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Any

import numpy as np

from record_live_block_base import RealGrowth, build_cells, fro, mean, perc, address
import local_patch_core as core


def write_csv(path: Path, rows: List[dict]) -> None:
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    keys: List[str] = []
    for r in rows:
        for k in r.keys():
            if k not in keys:
                keys.append(k)
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        for r in rows:
            w.writerow(r)


def build_model(level: int, case: str) -> RealGrowth:
    if case == "real_growth":
        m = RealGrowth(growth_rule="sequential", br_ancestor=0.045, br_sibling=0.035, shell_normalized=True)
    elif case == "symmetrized_birth":
        m = RealGrowth(growth_rule="symmetrized_birth", br_ancestor=0.045, br_sibling=0.035, shell_normalized=True)
    elif case == "no_backreaction":
        m = RealGrowth(growth_rule="sequential", br_ancestor=0.0, br_sibling=0.0, shell_normalized=True)
    else:
        raise ValueError(case)
    m.grow(level)
    return m


def node_level_counts(model: RealGrowth) -> Dict[int, int]:
    out: Dict[int, int] = defaultdict(int)
    for n in model.nodes.values():
        out[int(n.level)] += 1
    return dict(out)



def node_address_map(model: RealGrowth) -> Dict[tuple, Any]:
    return {address(model, nid): n for nid, n in model.nodes.items()}


def last_growth_delta_rows(level: int, case: str) -> List[dict]:
    """Compare the same provenance nodes before and after the last growth shell."""
    if level < 1:
        return []
    prev = build_model(level - 1, case)
    curr = build_model(level, case)
    prev_map = node_address_map(prev)
    curr_map = node_address_map(curr)
    rows: List[dict] = []
    for addr, n0 in prev_map.items():
        n1 = curr_map.get(addr)
        if n1 is None:
            continue
        rows.append({
            "global_level": level,
            "case": case,
            "address": ".".join(map(str, addr)) if addr else "root",
            "node_level": int(n0.level),
            "distance_to_active_parent_frontier": int((level - 1) - n0.level),
            "g_prev": float(n0.g),
            "g_curr": float(n1.g),
            "g_delta": float(n1.g - n0.g),
            "g_delta_abs": abs(float(n1.g - n0.g)),
            "relative_delta": abs(float(n1.g - n0.g)) / (abs(float(n0.g)) + 1e-12),
            "is_root": int(n0.level == 0),
            "is_active_parent_frontier": int(n0.level == level - 1),
        })
    return rows


def cell_activity_rows(level: int, case: str) -> List[dict]:
    model = build_model(level, case)
    cells = build_cells(model, level)
    counts = node_level_counts(model)
    rows: List[dict] = []
    for c in cells.values():
        completion_to_live = c.live1["D_total"] - c.live0["D_total"]
        child_drift = model.child_g_vector(c.parent) - model.completion_child_g[c.parent]
        row = {
            "global_level": level,
            "case": case,
            "parent": c.parent,
            "address": ".".join(map(str, c.address)) if c.address else "root",
            "cell_level": c.level,
            "completion_level": c.completion_level,
            "age": c.age,
            "frontier_distance": c.age,
            "nodes_at_cell_level": counts.get(c.level, 0),
            "aging_norm": fro(c.aging),
            "handoff_norm": fro(c.handoff),
            "completion_to_live_norm": fro(completion_to_live),
            "child_g_drift_norm": float(np.linalg.norm(child_drift)),
            "desc_load_norm": float(np.linalg.norm(model.child_descendant_loads(c.parent))),
            "ancestor_env_completion": float(model.completion_ancestor_env[c.parent]),
            "ancestor_env_live": float(model.ancestor_env_load(c.parent)),
            "ancestor_env_drift": float(model.ancestor_env_load(c.parent) - model.completion_ancestor_env[c.parent]),
            "is_root": int(c.level == 0),
            "is_active_completion_shell": int(c.age == 0),
        }
        rows.append(row)
    return rows


def summarize_by(rows: List[dict], key: str) -> List[dict]:
    groups: Dict[Any, List[dict]] = defaultdict(list)
    for r in rows:
        groups[r[key]].append(r)
    out: List[dict] = []
    metrics = [
        "aging_norm", "handoff_norm", "completion_to_live_norm",
        "child_g_drift_norm", "desc_load_norm", "ancestor_env_drift",
    ]
    for kval, gs in sorted(groups.items(), key=lambda kv: kv[0]):
        row: dict = {key: kval, "count": len(gs)}
        for m in metrics:
            vals = [float(g[m]) for g in gs]
            row[f"{m}_mean"] = mean(vals)
            row[f"{m}_p50"] = perc(vals, 50)
            row[f"{m}_p95"] = perc(vals, 95)
            row[f"{m}_max"] = max(vals) if vals else float("nan")
        out.append(row)
    return out



def summarize_metrics_by(rows: List[dict], key: str, metrics: List[str]) -> List[dict]:
    groups: Dict[Any, List[dict]] = defaultdict(list)
    for r in rows:
        groups[r[key]].append(r)
    out: List[dict] = []
    for kval, gs in sorted(groups.items(), key=lambda kv: kv[0]):
        row: dict = {key: kval, "count": len(gs)}
        for m in metrics:
            vals = [float(g[m]) for g in gs if m in g and np.isfinite(float(g[m]))]
            row[f"{m}_mean"] = mean(vals)
            row[f"{m}_p50"] = perc(vals, 50)
            row[f"{m}_p95"] = perc(vals, 95)
            row[f"{m}_max"] = max(vals) if vals else float("nan")
        out.append(row)
    return out


def patch_scaling_rows(levels: List[int], patch_sizes: List[int], cases: List[str], args: argparse.Namespace) -> tuple[List[dict], List[dict], List[dict]]:
    all_rows: List[dict] = []
    all_summaries: List[dict] = []
    primary_rows: List[dict] = []
    for level in levels:
        for patch_size in patch_sizes:
            ns = argparse.Namespace(
                levels=[level],
                cases=cases,
                tiers=[args.tier],
                patch_size=patch_size,
                max_patches_per_control=args.max_patches_per_control,
                degree=args.degree,
                word_cap=args.word_cap,
                mult_sample_cap=args.mult_sample_cap,
                operator_mode="triangular_handoff",
                longitudinal_mode="triangular_record_live",
                metric_source="record_live_block",
                ridge=1e-9,
                tol=1e-9,
                seed=args.seed + 100 * level + patch_size,
                out=args.out / f"tmp_patch_L{level}_k{patch_size}",
                include_projectors=True,
                include_diag_powers=True,
            )
            rows, summaries = core.run(ns)
            for r in rows:
                r["scan_level"] = level
                r["scan_patch_size"] = patch_size
                r["scan_degree"] = args.degree
                r["scan_tier"] = args.tier
                # active shell measure for selected patches
                r["patch_active_fraction"] = float(r.get("patch_age_mean", 0.0) == 0.0)
                r["patch_mean_age"] = r.get("patch_age_mean", float("nan"))
            for s in summaries:
                s["scan_level"] = level
                s["scan_patch_size"] = patch_size
                s["scan_degree"] = args.degree
                s["scan_tier"] = args.tier
            all_rows.extend(rows)
            all_summaries.extend(summaries)
            for s in summaries:
                label = str(s.get("label", ""))
                if ":same_suffix_multi_history_patch:" in label:
                    primary_rows.append({
                        "label": label,
                        "level": level,
                        "patch_size": patch_size,
                        "degree": args.degree,
                        "tier": args.tier,
                        "count": s.get("count"),
                        "word_dim_mean": s.get("word_space_dim_mean"),
                        "word_frac_mean": s.get("word_space_fraction_full_mean"),
                        "star_basis_rel_mean": s.get("star_basis_rel_mean_mean"),
                        "star_seed_rel_mean": s.get("star_seed_rel_mean_mean"),
                        "mult_left_rel_mean": s.get("mult_left_rel_mean_mean"),
                        "comm_rel_mean": s.get("patch_comm_C_rel_mean_mean"),
                        "metric_cond_mean": s.get("patch_metric_global_cond_mean"),
                    })
    return all_rows, all_summaries, primary_rows


def write_results(path: Path, args: argparse.Namespace, activity_by_age: List[dict], activity_by_level: List[dict], delta_by_node_level: List[dict], delta_by_distance: List[dict], primary_patch: List[dict]) -> None:
    lines: List[str] = []
    lines.append("# RESULTS: activity-gradient + patch-size star-closure scaling")
    lines.append("")
    lines.append("## Gate")
    lines.append("")
    lines.append("This diagnostic checks whether the strongest live DtN/Record drift is localized near the active growth/completion shell, and whether local patch operator-system closure improves with patch size/level. It remains pre-J and pre-C*.")
    lines.append("")
    lines.append("## Last-shell incremental node-change gradient")
    lines.append("")
    lines.append("This compares the same provenance nodes before and after the last growth shell. It is the cleaner test of your point: newest parent/frontier nodes receive the largest current update; the root receives the smallest shell-normalized marginal update.")
    lines.append("")
    lines.append("```text")
    lines.append("node_level count mean_abs_delta mean_rel_delta max_abs_delta")
    for r in delta_by_node_level:
        lines.append(
            f"{int(r['node_level']):<10} {int(r['count']):<5} "
            f"{float(r['g_delta_abs_mean']):<14.5g} "
            f"{float(r['relative_delta_mean']):<14.5g} "
            f"{float(r['g_delta_abs_max']):<14.5g}"
        )
    lines.append("```")
    lines.append("")
    lines.append("```text")
    lines.append("distance_to_frontier count mean_abs_delta mean_rel_delta max_abs_delta")
    for r in delta_by_distance:
        lines.append(
            f"{int(r['distance_to_active_parent_frontier']):<20} {int(r['count']):<5} "
            f"{float(r['g_delta_abs_mean']):<14.5g} "
            f"{float(r['relative_delta_mean']):<14.5g} "
            f"{float(r['g_delta_abs_max']):<14.5g}"
        )
    lines.append("```")
    lines.append("")

    lines.append("## Activity gradient by age")
    lines.append("")
    lines.append("Age is measured from completion to the final growth level. Age 0 is the active/latest completed shell; root/interior cells have large age.")
    lines.append("")
    lines.append("```text")
    lines.append("age count aging_mean handoff_mean child_g_drift desc_load ancestor_env_drift")
    for r in activity_by_age:
        lines.append(
            f"{int(r['age']):<3} {int(r['count']):<5} "
            f"{float(r['aging_norm_mean']):<11.5g} "
            f"{float(r['handoff_norm_mean']):<12.5g} "
            f"{float(r['child_g_drift_norm_mean']):<13.5g} "
            f"{float(r['desc_load_norm_mean']):<9.5g} "
            f"{float(r['ancestor_env_drift_mean']):<12.5g}"
        )
    lines.append("```")
    lines.append("")
    lines.append("## Activity gradient by cell level")
    lines.append("")
    lines.append("```text")
    lines.append("cell_level count aging_mean handoff_mean child_g_drift desc_load")
    for r in activity_by_level:
        lines.append(
            f"{int(r['cell_level']):<10} {int(r['count']):<5} "
            f"{float(r['aging_norm_mean']):<11.5g} "
            f"{float(r['handoff_norm_mean']):<12.5g} "
            f"{float(r['child_g_drift_norm_mean']):<13.5g} "
            f"{float(r['desc_load_norm_mean']):<9.5g}"
        )
    lines.append("```")
    lines.append("")
    lines.append("## Primary same-suffix patch closure scan")
    lines.append("")
    lines.append("```text")
    lines.append("L k count word_dim star_basis star_seed mult comm_rel")
    for r in primary_patch:
        lines.append(
            f"{int(r['level']):<2} {int(r['patch_size']):<2} {int(r['count']):<5} "
            f"{float(r['word_dim_mean']):<8.1f} "
            f"{float(r['star_basis_rel_mean']):<10.4f} "
            f"{float(r['star_seed_rel_mean']):<9.4f} "
            f"{float(r['mult_left_rel_mean']):<7.4f} "
            f"{float(r['comm_rel_mean']):<8.4f}"
        )
    lines.append("```")
    lines.append("")
    lines.append("## Interpretation")
    lines.append("")
    lines.append("The activity table should be read before the closure table: if the live changes concentrate near the active/latest shell while the root changes least, then a local operator-system candidate must be tested on active growth-defined patches, not on an isolated old/root subsystem. A positive local-net/*-limit signal would require star and multiplication residuals to decrease consistently as patch size and level grow. At this finite sampled scale, the scan is only diagnostic.")
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--levels", nargs="+", type=int, default=[4, 5, 6])
    ap.add_argument("--activity-level", type=int, default=6)
    ap.add_argument("--patch-sizes", nargs="+", type=int, default=[2, 3])
    ap.add_argument("--degree", type=int, default=2)
    ap.add_argument("--tier", default="base_patch")
    ap.add_argument("--cases", nargs="+", default=["real_growth"])
    ap.add_argument("--max-patches-per-control", type=int, default=1)
    ap.add_argument("--word-cap", type=int, default=100)
    ap.add_argument("--mult-sample-cap", type=int, default=40)
    ap.add_argument("--seed", type=int, default=20260621)
    ap.add_argument("--out", type=Path, default=Path("activity_gradient_patch_scaling_out"))
    args = ap.parse_args()

    t0 = time.time()
    args.out.mkdir(parents=True, exist_ok=True)

    activity_rows: List[dict] = []
    for case in args.cases:
        activity_rows.extend(cell_activity_rows(args.activity_level, case))
    activity_by_age = summarize_by(activity_rows, "age")
    activity_by_level = summarize_by(activity_rows, "cell_level")

    delta_rows: List[dict] = []
    for case in args.cases:
        delta_rows.extend(last_growth_delta_rows(args.activity_level, case))
    delta_by_node_level = summarize_metrics_by(delta_rows, "node_level", ["g_delta_abs", "relative_delta"])
    delta_by_distance = summarize_metrics_by(delta_rows, "distance_to_active_parent_frontier", ["g_delta_abs", "relative_delta"])

    patch_rows, patch_summaries, primary_patch = patch_scaling_rows(args.levels, args.patch_sizes, args.cases, args)

    write_csv(args.out / "cell_activity_rows.csv", activity_rows)
    write_csv(args.out / "activity_by_age.csv", activity_by_age)
    write_csv(args.out / "activity_by_cell_level.csv", activity_by_level)
    write_csv(args.out / "last_growth_delta_rows.csv", delta_rows)
    write_csv(args.out / "last_growth_delta_by_node_level.csv", delta_by_node_level)
    write_csv(args.out / "last_growth_delta_by_distance.csv", delta_by_distance)
    write_csv(args.out / "patch_rows_all.csv", patch_rows)
    write_csv(args.out / "patch_summary_table_all.csv", patch_summaries)
    write_csv(args.out / "primary_same_suffix_patch_summary.csv", primary_patch)
    write_results(args.out / "RESULTS_activity_gradient_patch_scaling.md", args, activity_by_age, activity_by_level, delta_by_node_level, delta_by_distance, primary_patch)
    summary = {
        "elapsed_seconds": time.time() - t0,
        "activity_rows": len(activity_rows),
        "last_growth_delta_rows": len(delta_rows),
        "patch_rows": len(patch_rows),
        "patch_summaries": len(patch_summaries),
        "primary_patch_rows": len(primary_patch),
        "activity_level": args.activity_level,
        "levels": args.levels,
        "patch_sizes": args.patch_sizes,
        "degree": args.degree,
        "cases": args.cases,
    }
    (args.out / "SUMMARY.txt").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    zip_path = args.out.with_suffix(".zip")
    with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as z:
        for p in args.out.rglob("*"):
            z.write(p, p.relative_to(args.out.parent))
    print(json.dumps({"out": str(args.out), "zip": str(zip_path), **summary}, indent=2))


if __name__ == "__main__":
    main()
