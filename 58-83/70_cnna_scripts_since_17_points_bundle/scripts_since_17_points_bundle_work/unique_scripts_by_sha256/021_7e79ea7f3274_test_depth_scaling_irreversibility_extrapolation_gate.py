#!/usr/bin/env python3
"""
CNNA depth-scaling irreversibility extrapolation gate.

Purpose
-------
Use the established script-1/script-2 ternary sequential growth model with true
Schur/DtN matrices and fixed-topology live relaxation.  Test the user's scaling
claim:

  irreversibility is not a constant single-birth property;
  it accumulates with growth depth because birth effects distribute into old
  conductances and then mix with live relaxation.

The test runs complete finite approximants for L2-L4 and deterministic sampled
frontier approximants for L5/L6, because full ternary true-Schur/DtN L5/L6 is
expensive in this environment.  Sampled rows are explicitly marked.

No J, i, Hodge, star, C*-claim, positivity axiom, Q/P target, or delta-beta gate
is used.
"""
from __future__ import annotations

import argparse, csv, json, math, shutil, zipfile
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np

from base_dynamic_boundary_role_schur_dtn import write_csv, EPS
from base_schur_dtn_birth_relaxation import SchurDtNRelaxationModel


def mean(xs): return float(np.mean(xs)) if xs else 0.0
def median(xs): return float(np.median(xs)) if xs else 0.0

def deterministic_frontier_sample(frontier: List[int], cap: int) -> List[int]:
    frontier = sorted(frontier)
    if cap <= 0 or len(frontier) <= cap:
        return frontier
    if cap == 1:
        return [frontier[len(frontier)//2]]
    idxs = np.linspace(0, len(frontier)-1, cap)
    out=[]
    seen=set()
    for x in idxs:
        i=int(round(float(x)))
        i=max(0,min(len(frontier)-1,i))
        if frontier[i] not in seen:
            out.append(frontier[i]); seen.add(frontier[i])
    # fill deterministically if rounding collided
    for x in frontier:
        if len(out)>=cap: break
        if x not in seen:
            out.append(x); seen.add(x)
    return sorted(out)

class DepthScalingRelaxModel(SchurDtNRelaxationModel):
    def __init__(self, *args, full_until_level:int=4, sample_cap:int=18, **kwargs):
        super().__init__(*args, **kwargs)
        self.full_until_level=full_until_level
        self.sample_cap=sample_cap
        self.current_growth_level=0
        self.sample_rows=[]

    def grow_level(self, frontier: List[int]) -> List[int]:
        level = self.current_growth_level + 1
        use_full = level <= self.full_until_level
        selected = list(frontier) if use_full else deterministic_frontier_sample(frontier, self.sample_cap)
        self.sample_rows.append({
            "variant": self.variant,
            "growth_level": level,
            "frontier_size_before": len(frontier),
            "selected_parent_count": len(selected),
            "sampled": int(not use_full),
            "sample_cap": self.sample_cap,
            "used_delta_beta_any": False,
        })
        next_frontier=[]
        for p in selected:
            for order in self.order_sequence:
                next_frontier.append(self.add_child(p, order))
        self.current_growth_level=level
        return next_frontier

    def add_child(self, parent:int, order:int) -> int:
        before_re=len(self.relax_event_rows)
        before_rr=len(self.relax_rows)
        parent_level=int(self.nodes[parent].level)
        child=super().add_child(parent, order)
        child_level=int(self.nodes[child].level)
        for rows in (self.relax_event_rows[before_re:], self.relax_rows[before_rr:]):
            for r in rows:
                r["parent_level"] = parent_level
                r["child_level"] = child_level
                r["growth_level"] = child_level
                r["used_delta_beta_any"] = False
        # enrich event row too
        self.event_rows[-1]["parent_level"] = parent_level
        self.event_rows[-1]["child_level"] = child_level
        self.event_rows[-1]["used_delta_beta_any"] = False
        return child


def cfgs(relax_steps:int):
    return [
        dict(variant="real_growth_linear_depth_scaling", mode="linear", alpha_env=0.22, br_ancestor=0.045, br_sibling=0.035, order_sequence=(1,2,3), relax_steps=relax_steps),
        dict(variant="log_growth_depth_scaling", mode="log", alpha_env=0.22, br_ancestor=0.045, br_sibling=0.035, order_sequence=(1,2,3), relax_steps=relax_steps),
        dict(variant="saturating_growth_depth_scaling", mode="saturating", alpha_env=0.90, br_ancestor=0.045, br_sibling=0.035, order_sequence=(1,2,3), relax_steps=relax_steps),
        dict(variant="strict_sym_depth_control", mode="linear", alpha_env=0.0, br_ancestor=0.0, br_sibling=0.0, order_sequence=(1,2,3), relax_steps=relax_steps),
    ]


def group_depth(relax_event_rows: List[dict], sampled_flag:int) -> List[dict]:
    groups: Dict[Tuple[str,int], List[dict]]={}
    for r in relax_event_rows:
        if "parent_level" not in r: continue
        groups.setdefault((str(r.get("variant","")), int(r["parent_level"])), []).append(r)
    out=[]
    for (variant, level), rs in sorted(groups.items()):
        gaps=[float(r.get("mean_record_vs_live_gap_fro",0.0)) for r in rs]
        birth=[abs(float(r.get("birth_delta_parent_dtn_eff",0.0))) for r in rs]
        drift=[float(r.get("mean_total_relax_drift_dtn_fro",0.0)) for r in rs]
        first=[float(r.get("mean_first_relax_delta_dtn_fro",0.0)) for r in rs]
        last=[float(r.get("mean_last_relax_delta_dtn_fro",0.0)) for r in rs]
        live_birth=[g/(b+EPS) for g,b in zip(gaps,birth)]
        nonrec=[int(g>1e-10) for g in gaps]
        out.append({
            "variant": variant,
            "parent_level": level,
            "rows": len(rs),
            "sampled_frontier_for_high_L": sampled_flag,
            "mean_record_live_gap_fro": mean(gaps),
            "median_record_live_gap_fro": median(gaps),
            "mean_birth_delta_parent_dtn_eff_abs": mean(birth),
            "mean_total_relax_drift_dtn_fro": mean(drift),
            "mean_live_birth_ratio": mean(live_birth),
            "median_live_birth_ratio": median(live_birth),
            "reverse_nonreconstructability_proxy_fraction": mean(nonrec),
            "mean_last_over_first_relax_ratio": mean([l/(f+EPS) for l,f in zip(last,first) if f>EPS]),
            "used_delta_beta_any": False,
        })
    return out


def depth_slope(depth_rows: List[dict], variant: str) -> dict:
    rs=sorted([r for r in depth_rows if r["variant"]==variant], key=lambda r:r["parent_level"])
    xs=np.array([float(r["parent_level"]) for r in rs], dtype=float)
    ys=np.array([float(r["mean_record_live_gap_fro"]) for r in rs], dtype=float)
    bs=np.array([float(r["mean_birth_delta_parent_dtn_eff_abs"]) for r in rs], dtype=float)
    mask=ys>1e-14
    slope=float(np.polyfit(xs[mask], np.log(ys[mask]),1)[0]) if np.sum(mask)>=2 else 0.0
    ratio=float(ys[mask][-1]/max(ys[mask][0],EPS)) if np.sum(mask)>=2 else 0.0
    bmask=bs>1e-14
    bslope=float(np.polyfit(xs[bmask], np.log(bs[bmask]),1)[0]) if np.sum(bmask)>=2 else 0.0
    return {"variant":variant,"depth_groups":len(rs),"log_gap_slope_vs_parent_level":slope,"deep_over_shallow_gap_ratio":ratio,"log_birth_delta_slope_vs_parent_level":bslope}


def run_one(cfg:dict, L:int, outdir:Path, full_until:int, cap:int) -> Tuple[dict,List[dict]]:
    m=DepthScalingRelaxModel(**cfg, full_until_level=full_until, sample_cap=cap)
    m.run(L)
    prefix=f"L{L}_{cfg['variant']}"
    write_csv(outdir/f"events_{prefix}.csv", m.event_rows)
    write_csv(outdir/f"relax_event_rows_{prefix}.csv", m.relax_event_rows)
    write_csv(outdir/f"relax_rows_{prefix}.csv", m.relax_rows)
    write_csv(outdir/f"levels_{prefix}.csv", m.level_rows)
    write_csv(outdir/f"sampling_rows_{prefix}.csv", m.sample_rows)
    sampled=int(any(r.get("sampled") for r in m.sample_rows))
    drows=group_depth(m.relax_event_rows, sampled)
    for r in drows:
        r["max_level"]=L
        r["frontier_sample_cap"]=cap
        r["full_until_level"]=full_until
    write_csv(outdir/f"depth_rows_{prefix}.csv", drows)
    final=m.level_rows[-1]
    rows=m.relax_event_rows
    gaps=[float(r.get("mean_record_vs_live_gap_fro",0.0)) for r in rows]
    births=[abs(float(r.get("birth_delta_parent_dtn_eff",0.0))) for r in rows]
    live_birth=[g/(b+EPS) for g,b in zip(gaps,births)]
    summ={
        "variant":cfg["variant"],"mode":cfg["mode"],"max_level":L,
        "sampled_frontier_for_high_L":sampled,
        "events":len(m.event_rows),"relax_event_rows":len(rows),
        "completed_triples":int(final.get("completed_triples",0)),
        "mean_record_live_gap_fro":mean(gaps),
        "median_record_live_gap_fro":median(gaps),
        "mean_birth_delta_parent_dtn_eff_abs":mean(births),
        "mean_live_birth_ratio":mean(live_birth),
        "reverse_nonreconstructability_proxy_fraction":mean([int(g>1e-10) for g in gaps]),
        "strict_sym_null_gate": bool(cfg["alpha_env"]==0.0 and mean(gaps)<1e-12),
        "used_delta_beta_any":False,
    }
    return summ,drows


def run_suite(levels:List[int], relax_steps:int, full_until:int, cap:int, outdir:Path):
    outdir.mkdir(parents=True, exist_ok=True)
    summaries=[]; all_depth=[]
    for L in levels:
        for cfg in cfgs(relax_steps):
            # keep high-L controls enough but not explosive; all variants are sampled by frontier cap after full_until
            summ,drows=run_one(cfg,L,outdir,full_until,cap)
            summaries.append(summ); all_depth.extend(drows)
    write_csv(outdir/"summary_by_variant_level.csv", summaries)
    write_csv(outdir/"depth_scaling_all.csv", all_depth)
    slopes=[]
    for v in sorted(set(r["variant"] for r in all_depth)):
        # slopes using each max_level separately and global aggregate
        slopes.append(depth_slope([r for r in all_depth if r["max_level"]==max(levels)], v) | {"slope_scope":"deepest_run_only"})
        slopes.append(depth_slope(all_depth, v) | {"slope_scope":"all_runs_aggregate"})
    write_csv(outdir/"depth_slope_by_variant.csv", slopes)
    summary={
        "model_family":"script1_script2_true_schur_birth_relaxation_depth_scaling_extrapolation",
        "levels":levels,"relax_steps":relax_steps,"full_until_level":full_until,"frontier_sample_cap_after_full_until":cap,
        "notes":[
            "Complete true-Schur/DtN finite approximants are run through full_until_level.",
            "Higher levels use deterministic evenly spaced frontier sampling to keep true-Schur/DtN computation feasible; sampled rows are marked.",
            "Irreversibility is measured by Record/Live Schur-DtN gap, live/birth ratio, relaxation drift and reverse-nonreconstructability proxy.",
            "No J, i, Hodge, star, positivity, C*-claim, Q/P target or delta-beta decision is used."
        ],
        "summaries":summaries,"depth_slopes":slopes,
    }
    (outdir/"summary.json").write_text(json.dumps(summary,indent=2),encoding="utf-8")
    return summary


def make_md(summary:dict,outdir:Path):
    maxL=max(summary["levels"])
    rows=[s for s in summary["summaries"] if s["max_level"]==maxL]
    table="\n".join([f"| {s['variant']} | {s['events']} | {s['relax_event_rows']} | {int(s['sampled_frontier_for_high_L'])} | {s['mean_record_live_gap_fro']:.4g} | {s['mean_birth_delta_parent_dtn_eff_abs']:.4g} | {s['mean_live_birth_ratio']:.4g} | {s['reverse_nonreconstructability_proxy_fraction']:.3f} |" for s in rows])
    slope_table="\n".join([f"| {s['slope_scope']} | {s['variant']} | {s['depth_groups']} | {s['log_gap_slope_vs_parent_level']:.3f} | {s['deep_over_shallow_gap_ratio']:.3f} | {s['log_birth_delta_slope_vs_parent_level']:.3f} |" for s in summary['depth_slopes']])
    md=f"""# CNNA depth-scaling irreversibility extrapolation gate

## Purpose

Test whether Record/Live irreversibility in the established script-1/script-2 true Schur/DtN growth model scales with depth rather than being a constant single-birth feature.

Finite true-Schur/DtN computation is complete through L{summary['full_until_level']}.  For higher levels this package uses deterministic sampled frontier approximants with cap `{summary['frontier_sample_cap_after_full_until']}`; all such rows are marked.  This avoids pretending that the high-L data are full ternary trees.

No `J`, `i`, Hodge, `*`, C*-claim, positivity axiom, Q/P target or `delta_beta` decision is used.

## Deepest-run summary

| variant | events | relax rows | sampled? | mean record/live gap | mean birth ΔDtN | live/birth ratio | reverse nonreconstructability proxy |
|---|---:|---:|---:|---:|---:|---:|---:|
{table}

## Slope audit

| scope | variant | depth groups | log gap slope vs parent level | deep/shallow gap ratio | log birth shock slope |
|---|---|---:|---:|---:|---:|
{slope_table}

## Interpretation

The key diagnostic is not a local `J` or orientation lock.  It is the depth behavior of the true Schur/DtN Record→Live gap after birth-plus-relaxation.

A positive log-gap slope and a deep/shallow ratio above one indicate that irreversibility grows with birth depth in these finite/sampled approximants.  The strict-sym control must stay null.  High-L sampled rows are evidence for a trend, not a full infinite-tree theorem.

## Next test

`test_stable_small_real_star_algebra_gate.py` should use the depth trend to filter for sufficiently deep, bridge-/passivity-positive Record→Live operator rows and then test whether the generated small real operator family becomes stable under the Record-DtN metric adjoint `#`, products and commutators.  This should remain weaker than a physical `*`-algebra or C*-claim.
"""
    (outdir/"RESULTS.md").write_text(md,encoding="utf-8")
    (outdir/"SUMMARY.md").write_text(md,encoding="utf-8")


def package(workdir:Path,outzip:Path):
    if outzip.exists(): outzip.unlink()
    with zipfile.ZipFile(outzip,"w",zipfile.ZIP_DEFLATED) as zf:
        for p in sorted(workdir.rglob("*")):
            if p.is_file(): zf.write(p,p.relative_to(workdir.parent))


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--levels",default="2,3,4,5,6")
    ap.add_argument("--relax-steps",type=int,default=4)
    ap.add_argument("--full-until-level",type=int,default=4)
    ap.add_argument("--frontier-sample-cap",type=int,default=18)
    ap.add_argument("--outdir",type=Path,default=Path("outputs"))
    args=ap.parse_args()
    levels=[int(x) for x in args.levels.split(',') if x]
    summary=run_suite(levels,args.relax_steps,args.full_until_level,args.frontier_sample_cap,args.outdir)
    make_md(summary,args.outdir)
    print(json.dumps({k:v for k,v in summary.items() if k not in ('summaries',)},indent=2))

if __name__=="__main__":
    main()
