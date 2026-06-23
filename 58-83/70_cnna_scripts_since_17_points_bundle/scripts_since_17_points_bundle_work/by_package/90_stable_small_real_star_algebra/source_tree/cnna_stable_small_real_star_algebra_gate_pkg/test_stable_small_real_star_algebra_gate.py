#!/usr/bin/env python3
"""
CNNA stable small real #/*-algebra-like operator-family gate.

Purpose
-------
After the depth-scaling irreversibility audit, test the next requested step:
not merely whether a metric adjoint # exists, but whether bridge-/passivity-positive
Record->Live relaxation operators form a small stable real operator family under:

  - record-DtN metric adjoint #_G;
  - products;
  - commutators;
  - depth filtering / stability across levels.

This is not a physical *-algebra or C*-algebra claim.  It is a derived-only,
finite, real, DtN-record/live operator-family diagnostic.

No J, i, Hodge, physical Hilbert adjoint, C*-norm, positivity axiom, Q/P target,
or delta-beta is used.
"""
from __future__ import annotations

import argparse, csv, json, zipfile
from pathlib import Path
from typing import Dict, List, Tuple
import numpy as np

from base_dynamic_boundary_role_schur_dtn import write_csv, EPS
from base_star_candidate import StarFamilyCandidateModel, configs, mean, median


def deterministic_frontier_sample(frontier: List[int], cap: int) -> List[int]:
    frontier=sorted(frontier)
    if cap<=0 or len(frontier)<=cap: return frontier
    idxs=np.linspace(0,len(frontier)-1,cap)
    out=[]; seen=set()
    for x in idxs:
        i=int(round(float(x))); i=max(0,min(len(frontier)-1,i))
        if frontier[i] not in seen:
            out.append(frontier[i]); seen.add(frontier[i])
    for x in frontier:
        if len(out)>=cap: break
        if x not in seen: out.append(x); seen.add(x)
    return sorted(out)

class SampledStarFamilyModel(StarFamilyCandidateModel):
    def __init__(self,*args,full_until_level:int=3,sample_cap:int=8,**kwargs):
        super().__init__(*args,**kwargs)
        self.full_until_level=full_until_level
        self.sample_cap=sample_cap
        self.current_growth_level=0
        self.sample_rows=[]

    def grow_level(self, frontier: List[int]) -> List[int]:
        level=self.current_growth_level+1
        use_full=level<=self.full_until_level
        selected=list(frontier) if use_full else deterministic_frontier_sample(frontier,self.sample_cap)
        self.sample_rows.append({"variant":self.variant,"growth_level":level,"frontier_size_before":len(frontier),"selected_parent_count":len(selected),"sampled":int(not use_full),"sample_cap":self.sample_cap,"used_delta_beta_any":False})
        nxt=[]
        for p in selected:
            for order in self.order_sequence:
                nxt.append(self.add_child(p,order))
        self.current_growth_level=level
        return nxt

    def add_child(self,parent:int,order:int)->int:
        before=len(self.star_rows)
        parent_level=int(self.nodes[parent].level)
        child=super().add_child(parent,order)
        for r in self.star_rows[before:]:
            r["parent_level"]=parent_level
            r["child_level"]=int(self.nodes[child].level)
            r["used_delta_beta_any"]=False
            r["stable_strong_star_algebra_gate"]=stable_row_gate(r)
            r["stable_weak_star_algebra_gate"]=weak_row_gate(r)
        return child

def weak_row_gate(r:dict)->bool:
    return bool(
        r.get("valid") and
        r.get("candidate_positive_gate") and
        float(r.get("star_span_residual_mean",1.0)) < 1e-8 and
        float(r.get("star_involutive_residual_mean",1.0)) < 1e-6 and
        float(r.get("anti_multiplicative_residual_mean",1.0)) < 1e-6 and
        float(r.get("product_closure_residual_mean",1.0)) < 0.35 and
        float(r.get("commutator_closure_residual_mean",1.0)) < 0.35
    )

def stable_row_gate(r:dict)->bool:
    return bool(
        r.get("valid") and
        r.get("candidate_positive_gate") and
        r.get("passivity_surrogate_gate") and
        float(r.get("star_span_residual_mean",1.0)) < 1e-8 and
        float(r.get("star_involutive_residual_mean",1.0)) < 1e-6 and
        float(r.get("anti_multiplicative_residual_mean",1.0)) < 1e-6 and
        float(r.get("product_closure_residual_mean",1.0)) < 0.12 and
        float(r.get("product_closure_residual_max",1.0)) < 0.25 and
        float(r.get("commutator_closure_residual_mean",1.0)) < 0.18 and
        float(r.get("commutator_closure_residual_max",1.0)) < 0.40 and
        int(r.get("operator_family_rank",0)) >= 2
    )

def summarize_rows(rows:List[dict], variant:str, L:int, sampled:int)->dict:
    valid=[r for r in rows if r.get("valid")]
    cand=[r for r in valid if r.get("candidate_positive_gate")]
    weak=[r for r in valid if r.get("stable_weak_star_algebra_gate")]
    strong=[r for r in valid if r.get("stable_strong_star_algebra_gate")]
    return {
        "variant":variant,"max_level":L,"sampled_frontier_for_high_L":sampled,
        "valid_rows":len(valid),"candidate_positive_rows":len(cand),
        "candidate_positive_fraction":len(cand)/len(valid) if valid else 0.0,
        "weak_star_family_pass_fraction":len(weak)/len(valid) if valid else 0.0,
        "strong_stable_star_algebra_pass_fraction":len(strong)/len(valid) if valid else 0.0,
        "mean_product_closure_residual_candidate_rows":mean([float(r.get("product_closure_residual_mean",0.0)) for r in cand]),
        "mean_commutator_closure_residual_candidate_rows":mean([float(r.get("commutator_closure_residual_mean",0.0)) for r in cand]),
        "mean_commutator_relative_norm_candidate_rows":mean([float(r.get("commutator_relative_norm_mean",0.0)) for r in cand]),
        "mean_record_live_gap_candidate_rows":mean([float(r.get("delta_fro",0.0)) for r in cand]),
        "mean_operator_family_rank_candidate_rows":mean([float(r.get("operator_family_rank",0.0)) for r in cand]),
        "strict_sym_null_gate": variant.startswith("strict_sym") and len(cand)==0 and len(strong)==0,
        "used_delta_beta_any":False,
    }

def depth_rows(rows:List[dict], sampled:int)->List[dict]:
    groups:Dict[Tuple[str,int],List[dict]]={}
    for r in rows:
        if not r.get("valid") or "parent_level" not in r: continue
        groups.setdefault((str(r.get("variant","")),int(r["parent_level"])),[]).append(r)
    out=[]
    for (v,pl),rs in sorted(groups.items()):
        valid=rs; cand=[r for r in rs if r.get("candidate_positive_gate")]; strong=[r for r in rs if r.get("stable_strong_star_algebra_gate")]
        out.append({"variant":v,"parent_level":pl,"rows":len(rs),"sampled_frontier_for_high_L":sampled,
                    "candidate_positive_fraction":len(cand)/len(valid) if valid else 0.0,
                    "strong_stable_star_algebra_pass_fraction":len(strong)/len(valid) if valid else 0.0,
                    "mean_gap_candidate_rows":mean([float(r.get("delta_fro",0.0)) for r in cand]),
                    "mean_product_resid_candidate_rows":mean([float(r.get("product_closure_residual_mean",0.0)) for r in cand]),
                    "mean_comm_resid_candidate_rows":mean([float(r.get("commutator_closure_residual_mean",0.0)) for r in cand]),
                    "used_delta_beta_any":False})
    return out

def filtered_configs(relax_steps:int):
    # Keep runtime bounded.  Real+saturating+strict are the decisive lines here.
    out=[]
    for c in configs(relax_steps):
        if c["variant"] in {"real_growth_linear_star_candidate","saturating_growth_star_candidate","strict_symmetrized_response_star_control"}:
            out.append(c)
    return out

def run_one(cfg:dict,L:int,outdir:Path,full_until:int,cap:int):
    m=SampledStarFamilyModel(**cfg,full_until_level=full_until,sample_cap=cap)
    m.run(L)
    sampled=int(any(r.get("sampled") for r in m.sample_rows))
    prefix=f"L{L}_{cfg['variant']}_stable_star"
    write_csv(outdir/f"events_{prefix}.csv",m.event_rows)
    write_csv(outdir/f"star_operator_rows_{prefix}.csv",m.star_rows)
    write_csv(outdir/f"depth_star_rows_{prefix}.csv",depth_rows(m.star_rows,sampled))
    write_csv(outdir/f"sampling_rows_{prefix}.csv",m.sample_rows)
    return summarize_rows(m.star_rows,cfg["variant"],L,sampled), depth_rows(m.star_rows,sampled)

def run_suite(levels:List[int],relax_steps:int,full_until:int,cap:int,outdir:Path):
    outdir.mkdir(parents=True,exist_ok=True)
    summaries=[]; all_depth=[]
    for L in levels:
        for cfg in filtered_configs(relax_steps):
            s,d=run_one(cfg,L,outdir,full_until,cap)
            summaries.append(s); 
            for x in d:
                x["max_level"]=L; all_depth.append(x)
    write_csv(outdir/"summary_by_variant_level.csv",summaries)
    write_csv(outdir/"depth_star_rows_all.csv",all_depth)
    summary={"model_family":"script1_script2_true_schur_live_semigroup_stable_small_real_star_algebra_gate", "levels":levels,"relax_steps":relax_steps,"full_until_level":full_until,"frontier_sample_cap_after_full_until":cap,
             "notes":["Uses established true-Schur/DtN Record->Live relaxation model.","Higher levels use deterministic sampled frontier approximants.","A row passes only if bridge/passivity positive and stable under #, product and commutator thresholds.","No J, i, Hodge, physical Hilbert adjoint, C*-norm, Q/P target or delta-beta gate is used."],"summaries":summaries}
    (outdir/"summary.json").write_text(json.dumps(summary,indent=2),encoding="utf-8")
    return summary

def make_md(summary:dict,outdir:Path):
    maxL=max(summary["levels"])
    rows=[s for s in summary["summaries"] if s["max_level"]==maxL]
    table="\n".join([f"| {s['variant']} | {s['valid_rows']} | {s['candidate_positive_fraction']:.3f} | {s['weak_star_family_pass_fraction']:.3f} | {s['strong_stable_star_algebra_pass_fraction']:.3f} | {s['mean_product_closure_residual_candidate_rows']:.3f} | {s['mean_commutator_closure_residual_candidate_rows']:.3f} | {s['mean_commutator_relative_norm_candidate_rows']:.3f} |" for s in rows])
    md=f"""# CNNA stable small real #/*-algebra-like operator-family gate

## Purpose

Use bridge-/passivity-positive Record→Live Schur/DtN relaxation rows to ask whether a **small real operator family** stabilizes under the record-DtN metric adjoint `#`, products and commutators.

This is deliberately weaker than a physical `*`-algebra and much weaker than a C*-algebra.  The metric adjoint `#_G` exists by linear algebra; it is not counted as a result unless the generated finite operator family is also approximately closed under products and commutators.

High-L rows are deterministic sampled frontier approximants after L{summary['full_until_level']}.

## Deepest-run summary

| variant | valid rows | candidate+ | weak # family pass | strong stable pass | product resid | commutator resid | commutator norm |
|---|---:|---:|---:|---:|---:|---:|---:|
{table}

## Interpretation

- `candidate+` means the Record→Live row is bridge-/passivity-positive and bounded.
- `weak # family pass` allows broad closure residuals.
- `strong stable pass` requires stricter simultaneous `#`, product and commutator closure.
- A positive `#` operation alone is not success; closure of the generated family is the actual gate.

## Main conclusion

The live semigroup supplies a real adjunction/passivity precursor, especially in the saturating response variant, but a robust small real `*`-algebra-like family is not yet uniformly established.  If strong-pass rows concentrate only in saturating/deep samples, the next step is not to claim a C*-structure but to identify the missing closure term or the correct coarse-grained family.
"""
    (outdir/"RESULTS.md").write_text(md,encoding="utf-8")
    (outdir/"SUMMARY.md").write_text(md,encoding="utf-8")

def main():
    ap=argparse.ArgumentParser(); ap.add_argument("--levels",default="3,4,5,6"); ap.add_argument("--relax-steps",type=int,default=4); ap.add_argument("--full-until-level",type=int,default=3); ap.add_argument("--frontier-sample-cap",type=int,default=8); ap.add_argument("--outdir",type=Path,default=Path("outputs"))
    args=ap.parse_args(); levels=[int(x) for x in args.levels.split(',') if x]
    summary=run_suite(levels,args.relax_steps,args.full_until_level,args.frontier_sample_cap,args.outdir); make_md(summary,args.outdir)
    print(json.dumps({k:v for k,v in summary.items() if k!='summaries'},indent=2))
if __name__=="__main__": main()
