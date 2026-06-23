"""
CNNA / growing real complement network
Patch-size / level scaling diagnostic for local real operator-system closure.

Purpose
-------
Previous diagnostics showed:
- single double-history pairs are too small for stable DtN-adjoint closure;
- same-suffix local patches improve over random controls but still do not close;
- adding a few obvious history intertwiners does not solve the problem.

This script tests whether the missing *-closure behaves like a finite-patch
artifact or a local-net/limit phenomenon. It runs the same patch operator-system
closure gate across:
- levels L,
- same-suffix patch sizes k,
- word degrees d,
- growth controls.

No J, no i, no complex phase, and no C*-norm are inserted. The test remains a
pre-* / pre-J diagnostic.
"""
from __future__ import annotations

import argparse
import csv
import json
import time
import zipfile
from pathlib import Path
from typing import Dict, List, Any

import numpy as np

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


def finite_mean(xs):
    vals = [float(x) for x in xs if x is not None and np.isfinite(float(x))]
    return float(np.mean(vals)) if vals else float("nan")


def finite_min(xs):
    vals = [float(x) for x in xs if x is not None and np.isfinite(float(x))]
    return float(np.min(vals)) if vals else float("nan")


def finite_max(xs):
    vals = [float(x) for x in xs if x is not None and np.isfinite(float(x))]
    return float(np.max(vals)) if vals else float("nan")


def trend_slope(points: List[tuple[float, float]]) -> float:
    pts = [(float(x), float(y)) for x, y in points if np.isfinite(float(y))]
    if len(pts) < 2:
        return float("nan")
    x = np.array([p[0] for p in pts], dtype=float)
    y = np.array([p[1] for p in pts], dtype=float)
    if np.std(x) < 1e-12:
        return float("nan")
    return float(np.polyfit(x, y, 1)[0])


def run_one(args: argparse.Namespace, level: int, patch_size: int, degree: int, tier: str, cases: List[str]) -> tuple[List[dict], List[dict]]:
    ns = argparse.Namespace(
        levels=[level],
        cases=cases,
        tiers=[tier],
        patch_size=patch_size,
        max_patches_per_control=args.max_patches_per_control,
        degree=degree,
        word_cap=args.word_cap,
        mult_sample_cap=args.mult_sample_cap,
        operator_mode=args.operator_mode,
        longitudinal_mode=args.longitudinal_mode,
        metric_source=args.metric_source,
        ridge=args.ridge,
        tol=args.tol,
        seed=args.seed + 1000 * level + 100 * patch_size + 10 * degree,
        out=args.out / f"tmp_L{level}_p{patch_size}_d{degree}_{tier}",
        include_projectors=args.include_projectors,
        include_diag_powers=args.include_diag_powers,
    )
    rows, summaries = core.run(ns)
    for r in rows:
        r["scan_level"] = level
        r["scan_patch_size"] = patch_size
        r["scan_degree"] = degree
        r["scan_tier"] = tier
    for s in summaries:
        s["scan_level"] = level
        s["scan_patch_size"] = patch_size
        s["scan_degree"] = degree
        s["scan_tier"] = tier
    return rows, summaries


def extract_primary_summaries(summaries: List[dict]) -> List[dict]:
    out = []
    for s in summaries:
        label = str(s.get("label", ""))
        if ":same_suffix_multi_history_patch:" in label:
            out.append(s)
    return out


def compact_primary_rows(primary: List[dict]) -> List[dict]:
    rows = []
    for s in primary:
        rows.append({
            "label": s.get("label"),
            "level": s.get("scan_level"),
            "patch_size": s.get("scan_patch_size"),
            "degree": s.get("scan_degree"),
            "tier": s.get("scan_tier"),
            "count": s.get("count"),
            "word_dim_mean": s.get("word_space_dim_mean"),
            "word_frac_mean": s.get("word_space_fraction_full_mean"),
            "star_basis_rel_mean": s.get("star_basis_rel_mean_mean"),
            "star_basis_rel_p95": s.get("star_basis_rel_p95_mean"),
            "star_basis_rel_lt_0p25": s.get("star_basis_rel_lt_0p25_fraction_mean"),
            "star_seed_rel_mean": s.get("star_seed_rel_mean_mean"),
            "star_seed_rel_lt_0p25": s.get("star_seed_rel_lt_0p25_fraction_mean"),
            "mult_left_rel_mean": s.get("mult_left_rel_mean_mean"),
            "mult_left_rel_p95": s.get("mult_left_rel_p95_mean"),
            "mult_left_rel_lt_0p25": s.get("mult_left_rel_lt_0p25_fraction_mean"),
            "comm_rel_mean": s.get("patch_comm_C_rel_mean_mean"),
            "metric_cond_mean": s.get("patch_metric_global_cond_mean"),
        })
    return rows


def build_trend_rows(primary_rows: List[dict]) -> List[dict]:
    trend_rows: List[dict] = []
    tiers = sorted({r["tier"] for r in primary_rows})
    degrees = sorted({int(r["degree"]) for r in primary_rows})
    patch_sizes = sorted({int(r["patch_size"]) for r in primary_rows})
    levels = sorted({int(r["level"]) for r in primary_rows})
    metrics = ["star_basis_rel_mean", "star_seed_rel_mean", "mult_left_rel_mean"]
    for tier in tiers:
        for degree in degrees:
            for metric in metrics:
                # level trend per patch size
                for ps in patch_sizes:
                    pts = [(r["level"], r[metric]) for r in primary_rows if r["tier"] == tier and int(r["degree"]) == degree and int(r["patch_size"]) == ps]
                    trend_rows.append({
                        "tier": tier,
                        "degree": degree,
                        "patch_size": ps,
                        "trend_axis": "level",
                        "metric": metric,
                        "slope": trend_slope(pts),
                        "values": json.dumps(pts),
                    })
                # patch-size trend per level
                for L in levels:
                    pts = [(r["patch_size"], r[metric]) for r in primary_rows if r["tier"] == tier and int(r["degree"]) == degree and int(r["level"]) == L]
                    trend_rows.append({
                        "tier": tier,
                        "degree": degree,
                        "level": L,
                        "trend_axis": "patch_size",
                        "metric": metric,
                        "slope": trend_slope(pts),
                        "values": json.dumps(pts),
                    })
    return trend_rows


def write_results(path: Path, args: argparse.Namespace, primary_rows: List[dict], trend_rows: List[dict]) -> None:
    def find(level, patch_size, degree, tier):
        for r in primary_rows:
            if int(r["level"]) == level and int(r["patch_size"]) == patch_size and int(r["degree"]) == degree and r["tier"] == tier:
                return r
        return None

    lines: List[str] = []
    lines.append("# RESULTS: patch-size / level scaling star-closure diagnostic")
    lines.append("")
    lines.append("## Gate")
    lines.append("")
    lines.append("This test asks whether the local real operator-system closure improves systematically with patch size, growth level, or word degree. It remains pre-J and pre-C*. No complex structure is inserted.")
    lines.append("")
    lines.append("## Parameters")
    lines.append("")
    lines.append("```json")
    lines.append(json.dumps({
        "levels": args.levels,
        "patch_sizes": args.patch_sizes,
        "degrees": args.degrees,
        "tiers": args.tiers,
        "cases": args.cases,
        "max_patches_per_control": args.max_patches_per_control,
        "word_cap": args.word_cap,
        "operator_mode": args.operator_mode,
        "longitudinal_mode": args.longitudinal_mode,
        "metric_source": args.metric_source,
    }, indent=2))
    lines.append("```")
    lines.append("")
    lines.append("## Primary same-suffix growth summaries")
    lines.append("")
    for tier in args.tiers:
        lines.append(f"### tier = `{tier}`")
        lines.append("")
        for degree in args.degrees:
            lines.append(f"#### degree = {degree}")
            lines.append("")
            lines.append("```text")
            header = "L  k  count  word_dim  star_basis  star_seed  mult  comm_rel"
            lines.append(header)
            for L in args.levels:
                for k in args.patch_sizes:
                    r = find(L, k, degree, tier)
                    if not r:
                        lines.append(f"{L:<2} {k:<2} -- no same-suffix patch")
                    else:
                        lines.append(
                            f"{L:<2} {k:<2} {int(r['count']):<6} "
                            f"{float(r['word_dim_mean']):<8.1f} "
                            f"{float(r['star_basis_rel_mean']):<10.4f} "
                            f"{float(r['star_seed_rel_mean']):<9.4f} "
                            f"{float(r['mult_left_rel_mean']):<7.4f} "
                            f"{float(r['comm_rel_mean']):<8.4f}"
                        )
            lines.append("```")
            lines.append("")
    lines.append("## Trend diagnosis")
    lines.append("")
    def slope_for(axis, metric, tier, degree, **kwargs):
        for t in trend_rows:
            if t.get("trend_axis") == axis and t.get("metric") == metric and t.get("tier") == tier and int(t.get("degree")) == degree:
                ok = True
                for k, v in kwargs.items():
                    if k not in t or int(t[k]) != int(v):
                        ok = False
                        break
                if ok:
                    return t.get("slope")
        return float("nan")
    lines.append("Negative slopes would indicate improving closure as the axis grows. Mixed or positive slopes mean no clear limiting trend at this scale.")
    lines.append("")
    for tier in args.tiers:
        for degree in args.degrees:
            lines.append(f"### trend tier={tier}, degree={degree}")
            lines.append("")
            lines.append("```text")
            for k in args.patch_sizes:
                sb = slope_for("level", "star_basis_rel_mean", tier, degree, patch_size=k)
                ss = slope_for("level", "star_seed_rel_mean", tier, degree, patch_size=k)
                ml = slope_for("level", "mult_left_rel_mean", tier, degree, patch_size=k)
                lines.append(f"level-slope k={k}: star_basis={sb:.5g}, star_seed={ss:.5g}, mult={ml:.5g}")
            for L in args.levels:
                sb = slope_for("patch_size", "star_basis_rel_mean", tier, degree, level=L)
                ss = slope_for("patch_size", "star_seed_rel_mean", tier, degree, level=L)
                ml = slope_for("patch_size", "mult_left_rel_mean", tier, degree, level=L)
                lines.append(f"patch-size-slope L={L}: star_basis={sb:.5g}, star_seed={ss:.5g}, mult={ml:.5g}")
            lines.append("```")
            lines.append("")
    lines.append("## Interpretation")
    lines.append("")
    lines.append("At this sampled scale, a genuine local-net/*-limit signal would require star and multiplication residuals to decrease consistently with patch size and/or level. If only selected entries improve while multiplication closure stays high, the current finite generated operator family is still too small or not the right local algebra carrier.")
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--levels", nargs="+", type=int, default=[4, 5, 6])
    ap.add_argument("--patch-sizes", nargs="+", type=int, default=[2, 3])
    ap.add_argument("--degrees", nargs="+", type=int, default=[2, 3])
    ap.add_argument("--tiers", nargs="+", default=["base_patch"])
    ap.add_argument("--cases", nargs="+", default=["real_growth"])
    ap.add_argument("--max-patches-per-control", type=int, default=2)
    ap.add_argument("--word-cap", type=int, default=220)
    ap.add_argument("--mult-sample-cap", type=int, default=120)
    ap.add_argument("--operator-mode", default="triangular_handoff")
    ap.add_argument("--longitudinal-mode", default="triangular_record_live")
    ap.add_argument("--metric-source", default="record_live_block")
    ap.add_argument("--ridge", type=float, default=1e-9)
    ap.add_argument("--tol", type=float, default=1e-9)
    ap.add_argument("--seed", type=int, default=20260621)
    ap.add_argument("--out", type=Path, default=Path("patch_size_level_scaling_star_closure_out"))
    ap.add_argument("--include-projectors", action="store_true", default=True)
    ap.add_argument("--include-diag-powers", action="store_true", default=True)
    args = ap.parse_args()

    t0 = time.time()
    args.out.mkdir(parents=True, exist_ok=True)
    all_rows: List[dict] = []
    all_summaries: List[dict] = []
    for level in args.levels:
        for patch_size in args.patch_sizes:
            for degree in args.degrees:
                for tier in args.tiers:
                    rows, summaries = run_one(args, level, patch_size, degree, tier, args.cases)
                    all_rows.extend(rows)
                    all_summaries.extend(summaries)
    primary = extract_primary_summaries(all_summaries)
    primary_rows = compact_primary_rows(primary)
    trend_rows = build_trend_rows(primary_rows)

    write_csv(args.out / "patch_rows_all.csv", all_rows)
    write_csv(args.out / "summary_table_all.csv", all_summaries)
    write_csv(args.out / "primary_same_suffix_summary.csv", primary_rows)
    write_csv(args.out / "trend_table.csv", trend_rows)
    write_results(args.out / "RESULTS_patch_size_level_scaling_star_closure.md", args, primary_rows, trend_rows)
    summary = {
        "elapsed_seconds": time.time() - t0,
        "patch_rows": len(all_rows),
        "summaries": len(all_summaries),
        "primary_rows": len(primary_rows),
        "trend_rows": len(trend_rows),
        "levels": args.levels,
        "patch_sizes": args.patch_sizes,
        "degrees": args.degrees,
        "tiers": args.tiers,
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
