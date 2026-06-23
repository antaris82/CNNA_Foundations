#!/usr/bin/env python3
"""
CNNA commutator-generator completion gate.

Purpose
-------
Continue the established script-1/script-2 true-Schur/DtN Record->Live
operator-family line after the coarse-grained live #/* closure test.

The preceding package showed that coarse-graining, especially mean_plus_principal,
improves product and commutator closure but does not yet give a stable small real
#/*-algebra-like family.  This test asks a narrower question:

  If we take the best coarse-grained real operator family and add ONLY derived
  commutator modes [A_i,A_j], does closure stabilize?

No arbitrary generators are added.  No J, i, Hodge, physical Hilbert adjoint,
C*-norm, Q/P target, positivity axiom, or delta-beta decision is used.

The metric adjoint #_G is still only diagnostic.  A pass requires simultaneous
# closure, product closure, commutator closure, and strict_sym null.
"""
from __future__ import annotations

import argparse, json, zipfile, shutil
from pathlib import Path
from typing import Dict, List, Tuple
import numpy as np

from base_dynamic_boundary_role_schur_dtn import DynamicSchurDtNRoleModel, write_csv, EPS
from base_schur_dtn_birth_relaxation import SchurDtNRelaxationModel
from base_live_semigroup_polarity_to_operator_bridge import zero_sum_basis, project_zero_sum, matrix_metrics, mean, median, sym, rel_fro
from base_star_candidate import configs, metric_adjoint, span_project_residual, numerical_rank_of_family, operator_sequence_family_metrics


def _norm(A: np.ndarray) -> float:
    return float(np.linalg.norm(A, ord='fro')) if A.size else 0.0


def deterministic_frontier_sample(frontier: List[int], cap: int) -> List[int]:
    frontier = sorted(frontier)
    if cap <= 0 or len(frontier) <= cap:
        return frontier
    idxs = np.linspace(0, len(frontier)-1, cap)
    out, seen = [], set()
    for x in idxs:
        i = int(round(float(x)))
        i = max(0, min(len(frontier)-1, i))
        if frontier[i] not in seen:
            out.append(frontier[i]); seen.add(frontier[i])
    for x in frontier:
        if len(out) >= cap: break
        if x not in seen:
            out.append(x); seen.add(x)
    return sorted(out)


class MatrixCollectingModel(SchurDtNRelaxationModel):
    def __init__(self, *args, axis_sign:int=1, relax_steps:int=4, full_until_level:int=3, sample_cap:int=8, **kwargs):
        super().__init__(*args, **kwargs)
        self.axis_sign = int(axis_sign)
        self.relax_steps = int(relax_steps)
        self.full_until_level = int(full_until_level)
        self.sample_cap = int(sample_cap)
        self.current_growth_level = 0
        self.matrix_rows: List[dict] = []
        self.metric_rows: List[dict] = []
        self.sample_rows: List[dict] = []

    def grow_level(self, frontier: List[int]) -> List[int]:
        level = self.current_growth_level + 1
        selected = list(frontier) if level <= self.full_until_level else deterministic_frontier_sample(frontier, self.sample_cap)
        self.sample_rows.append({"growth_level": level, "frontier_size_before": len(frontier), "selected_parent_count": len(selected), "sampled": int(level > self.full_until_level), "used_delta_beta_any": False})
        nxt=[]
        for p in selected:
            for order in self.order_sequence:
                nxt.append(self.add_child(p, order))
        self.current_growth_level = level
        return nxt

    def add_child(self, parent:int, order:int)->int:
        parent_line_before = self.parent_line(parent)
        older_before = list(self.nodes[parent].children)
        parent_level = int(self.nodes[parent].level)
        child = DynamicSchurDtNRoleModel.add_child(self, parent, order)
        self._register_new_base_edges()
        birth_event = self.event_rows[-1]
        birth_id = int(birth_event["birth_id"])
        after_nodes, after_edges = self.snapshot()
        fixed=[]
        for cut in parent_line_before:
            rec = self.dtn_for_boundary(after_nodes, after_edges, cut)
            fixed.append({"cut":cut, "role_kind":"ancestor_parent_line_cut", "record_dtn":rec, "seq":[]})
        for cut in older_before:
            rec = self.dtn_for_boundary(after_nodes, after_edges, cut)
            fixed.append({"cut":cut, "role_kind":"older_sibling_backreaction_cut", "record_dtn":rec, "seq":[]})
        for step in range(1, self.relax_steps+1):
            self.live_relax_step()
            live_nodes, live_edges = self.snapshot()
            for fc in fixed:
                live = self.fixed_boundary_dtn(live_nodes, live_edges, fc["cut"], list(fc["record_dtn"].boundary))
                fc["seq"].append(live)
        cand_count=0; valid_count=0; gap_vals=[]
        for fc in fixed:
            record = fc["record_dtn"]; seq=fc["seq"]
            base = {
                "variant": self.variant, "mode": self.mode, "birth_id": birth_id, "birth_order": order,
                "parent": parent, "child": child, "parent_level": parent_level, "child_level": int(self.nodes[child].level),
                "cut_node": fc["cut"], "cut_depth_from_parent": parent_line_before.index(fc["cut"]) if fc["cut"] in parent_line_before else -1,
                "role_kind": fc["role_kind"], "boundary": " ".join(map(str, record.boundary)), "used_delta_beta_any": False,
            }
            met = operator_sequence_family_metrics(record, seq, fc["cut"], self.axis_sign)
            mrow = dict(base); mrow.update({k:v for k,v in met.items() if not isinstance(v, np.ndarray)})
            mrow["coarse_level"] = "raw_per_cut"
            self.metric_rows.append(mrow)
            if met.get("valid"):
                valid_count += 1
                if met.get("candidate_positive_gate"): cand_count += 1
                gap_vals.append(float(met.get("delta_fro",0.0)))
            bundle = build_operator_bundle(record, seq, fc["cut"], self.axis_sign)
            if bundle.get("valid"):
                brow = dict(base); brow.update(bundle)
                self.matrix_rows.append(brow)
        birth_event.update({"matrix_valid_rows": valid_count, "matrix_candidate_rows": cand_count, "mean_gap": mean(gap_vals), "used_delta_beta_any": False})
        return child


def build_operator_bundle(record, seq, cut:int, axis_sign:int)->dict:
    mm = matrix_metrics(record, seq, cut, axis_sign)
    if not mm.get("valid"):
        return {"valid": False, "reason": mm.get("reason", "invalid")}
    G = project_zero_sum(sym(np.asarray(record.lam, dtype=float)))
    d = G.shape[0]
    if d <= 0:
        return {"valid": False, "reason": "zero_projected_dim"}
    evals = np.linalg.eigvalsh(sym(G))
    if evals.size == 0 or float(np.min(evals)) < -1e-8:
        return {"valid": False, "reason": "record_metric_not_psd"}
    reg = max(1e-10, 1e-8 * max(float(np.max(evals)), 1.0))
    Gp = np.linalg.pinv(G + reg*np.eye(d), rcond=1e-10)
    As=[]; prevP=G.copy(); inc=[]
    for item in seq:
        P = project_zero_sum(sym(np.asarray(item.lam, dtype=float)))
        if P.shape != G.shape: continue
        Delta = sym(P - G)
        A = Gp @ Delta
        As.append(A)
        inc.append(_norm(P - prevP)); prevP=P
    if not As:
        return {"valid": False, "reason": "no_sequence_operators"}
    candidate_positive = bool(
        mm.get("valid")
        and mm.get("bounded_record_metric_gate", False)
        and mm.get("passivity_surrogate_gate", False)
        and mm.get("bounded_decay_gate", False)
        and float(mm.get("delta_fro", 0.0)) > 1e-10
    )
    out = dict(mm)
    out.update({"valid": True, "G": G, "As": As, "projected_dim": d, "operator_count": len(As), "candidate_positive_gate": candidate_positive})
    return out


def normalize_family(mats: List[np.ndarray]) -> List[np.ndarray]:
    out=[]
    for A in mats:
        A=np.asarray(A,dtype=float)
        n=_norm(A)
        if n > 1e-12:
            out.append(A/n)
    return out


def independent_extend(basis: List[np.ndarray], candidates: List[np.ndarray], max_new:int=8, tol:float=1e-8) -> Tuple[List[np.ndarray], int, int]:
    """Add independent normalized candidates to basis by span residual."""
    out=list(basis)
    added=0; seen=0
    for C in candidates:
        n=_norm(C)
        if n <= 1e-12: continue
        seen += 1
        Cn=C/n
        if not out or span_project_residual(Cn,out) > tol:
            out.append(Cn); added += 1
            if added >= max_new: break
    return out, added, seen


def derived_commutators(gens: List[np.ndarray]) -> List[np.ndarray]:
    cs=[]
    for i,A in enumerate(gens):
        for j,B in enumerate(gens):
            if j <= i: continue
            C=A@B-B@A
            if _norm(C) > 1e-12:
                cs.append(C)
    return cs


def completion_generators(initial: List[np.ndarray], mode:str, max_new:int=8) -> Tuple[List[np.ndarray], dict]:
    base=normalize_family(initial)
    meta={"initial_generator_count": len(base), "commutator_candidates_seen":0, "commutator_modes_added":0, "second_order_candidates_seen":0, "second_order_modes_added":0}
    if mode == "base":
        return base, meta
    C1=derived_commutators(base)
    ext, added, seen=independent_extend(base,C1,max_new=max_new)
    meta["commutator_candidates_seen"]=seen; meta["commutator_modes_added"]=added
    if mode == "commutator_completed_1":
        return ext, meta
    if mode == "commutator_completed_2":
        # second-order commutators with the just-extended family, still derived from previous modes.
        C2=derived_commutators(ext)
        ext2, added2, seen2=independent_extend(ext,C2,max_new=max_new)
        meta["second_order_candidates_seen"]=seen2; meta["second_order_modes_added"]=added2
        return ext2, meta
    return ext, meta


def closure_metrics(G: np.ndarray, generators: List[np.ndarray], label:dict)->dict:
    if G.size == 0 or not generators:
        return {**label, "valid": False, "reason":"empty_generators"}
    d = G.shape[0]
    I = np.eye(d)
    raw=[I]+[np.asarray(A,dtype=float) for A in generators if np.asarray(A).shape == G.shape]
    if len(raw) <= 1:
        return {**label, "valid": False, "reason":"no_matching_generators"}
    basis=normalize_family(raw)
    fam_rank = numerical_rank_of_family(basis)
    sample = raw[1:min(len(raw),6)]
    star_res=[]; invol_res=[]; anti_res=[]; prod_res=[]; comm_res=[]; comm_norm=[]
    for A in sample:
        Ah = metric_adjoint(A,G)
        star_res.append(span_project_residual(Ah,basis))
        invol_res.append(rel_fro(metric_adjoint(Ah,G),A))
    for A in sample:
        for B in sample:
            AB=A@B; BA=B@A; C=AB-BA
            prod_res.append(span_project_residual(AB,basis))
            comm_res.append(span_project_residual(C,basis))
            comm_norm.append(_norm(C)/(_norm(AB)+_norm(BA)+EPS))
            anti_res.append(rel_fro(metric_adjoint(AB,G), metric_adjoint(B,G)@metric_adjoint(A,G)))
    out = dict(label)
    out.update({
        "valid": True,
        "operator_family_dim": d,
        "operator_family_rank": fam_rank,
        "generator_count": len(sample),
        "span_basis_count": len(basis),
        "star_span_residual_mean": float(np.mean(star_res)) if star_res else 0.0,
        "star_involutive_residual_mean": float(np.mean(invol_res)) if invol_res else 0.0,
        "anti_multiplicative_residual_mean": float(np.mean(anti_res)) if anti_res else 0.0,
        "product_closure_residual_mean": float(np.mean(prod_res)) if prod_res else 0.0,
        "product_closure_residual_max": float(np.max(prod_res)) if prod_res else 0.0,
        "commutator_closure_residual_mean": float(np.mean(comm_res)) if comm_res else 0.0,
        "commutator_closure_residual_max": float(np.max(comm_res)) if comm_res else 0.0,
        "commutator_relative_norm_mean": float(np.mean(comm_norm)) if comm_norm else 0.0,
        "used_delta_beta_any": False,
    })
    weak = (
        out["star_span_residual_mean"] < 1e-8 and out["star_involutive_residual_mean"] < 1e-6 and
        out["anti_multiplicative_residual_mean"] < 1e-6 and out["product_closure_residual_mean"] < 0.35 and
        out["commutator_closure_residual_mean"] < 0.35 and fam_rank >= 2
    )
    strong = (
        out["star_span_residual_mean"] < 1e-8 and out["star_involutive_residual_mean"] < 1e-6 and
        out["anti_multiplicative_residual_mean"] < 1e-6 and out["product_closure_residual_mean"] < 0.12 and
        out["product_closure_residual_max"] < 0.30 and out["commutator_closure_residual_mean"] < 0.18 and
        out["commutator_closure_residual_max"] < 0.45 and fam_rank >= 2
    )
    out["weak_star_family_pass"] = bool(weak)
    out["strong_stable_star_family_pass"] = bool(strong)
    return out


def average_matrices(mats: List[np.ndarray]) -> np.ndarray:
    return np.mean(np.stack(mats, axis=0), axis=0)


def group_key(row:dict, level_mode:str)->Tuple:
    if level_mode == "depth":
        return (row["variant"], row["parent_level"], row["role_kind"], int(row["projected_dim"]))
    if level_mode == "depth_boundarydim":
        return (row["variant"], row["parent_level"], int(row["projected_dim"]))
    return (row["variant"], int(row["projected_dim"]))


def initial_coarse_gens(rs:List[dict], mode_name:str, svd_k:int=3)->Tuple[np.ndarray,List[np.ndarray],dict]:
    d=int(rs[0]["projected_dim"])
    G=average_matrices([r["G"] for r in rs])
    minlen=min(len(r["As"]) for r in rs)
    if minlen <= 0:
        return G, [], {"min_sequence_len":0, "svd_singular_values":""}
    meta={"min_sequence_len":minlen, "svd_singular_values":""}
    gens=[]
    if mode_name in ("depth_group_average", "mean_plus_principal"):
        for k in range(minlen):
            gens.append(average_matrices([r["As"][k] for r in rs]))
    if mode_name in ("principal_stable_modes", "mean_plus_principal"):
        vecs=[]
        for r in rs:
            for A in r["As"][:minlen]:
                vecs.append(A.reshape(-1))
        if vecs:
            M=np.stack(vecs,axis=1)
            U,S,Vt=np.linalg.svd(M,full_matrices=False)
            meta["svd_singular_values"]=" ".join(f"{float(x):.6g}" for x in S[:min(6,len(S))])
            for j in range(min(svd_k, U.shape[1])):
                if S[j] > 1e-12:
                    gens.append(U[:,j].reshape(d,d))
    return G, gens, meta


def coarse_group_completion_metrics(rows: List[dict], level_mode:str, initial_mode:str, completion_mode:str, svd_k:int=3, max_comm_new:int=8)->List[dict]:
    groups:Dict[Tuple,List[dict]]={}
    for r in rows:
        if not r.get("valid") or not r.get("candidate_positive_gate"):
            continue
        groups.setdefault(group_key(r, level_mode), []).append(r)
    out=[]
    for key, rs in sorted(groups.items(), key=lambda kv: str(kv[0])):
        if len(rs) < 3: continue
        d=int(rs[0]["projected_dim"])
        G, gens0, imeta = initial_coarse_gens(rs, initial_mode, svd_k=svd_k)
        if not gens0: continue
        gens, cmeta = completion_generators(gens0, completion_mode, max_new=max_comm_new)
        label={
            "coarse_level": initial_mode,
            "completion_level": completion_mode,
            "group_mode": level_mode,
            "group_key": str(key),
            "variant": rs[0]["variant"],
            "parent_level": key[1] if len(key)>2 and isinstance(key[1], int) else -1,
            "rows_in_group": len(rs),
            "projected_dim": d,
            "mean_delta_fro": mean([float(r.get("delta_fro",0.0)) for r in rs]),
            "mean_parent_level": mean([float(r.get("parent_level",0.0)) for r in rs]),
            **imeta,
            **cmeta,
        }
        out.append(closure_metrics(G, gens, label))
    return out


def summarize(rows:List[dict], variant:str, coarse_level:str, completion_level:str)->dict:
    rs=[r for r in rows if r.get("valid") and r.get("variant")==variant and r.get("coarse_level")==coarse_level and r.get("completion_level")==completion_level]
    return {
        "variant": variant, "coarse_level": coarse_level, "completion_level": completion_level, "rows": len(rs),
        "weak_pass_fraction": mean([int(r.get("weak_star_family_pass",False)) for r in rs]),
        "strong_pass_fraction": mean([int(r.get("strong_stable_star_family_pass",False)) for r in rs]),
        "mean_product_resid": mean([float(r.get("product_closure_residual_mean",0.0)) for r in rs]),
        "mean_commutator_resid": mean([float(r.get("commutator_closure_residual_mean",0.0)) for r in rs]),
        "mean_commutator_max": mean([float(r.get("commutator_closure_residual_max",0.0)) for r in rs]),
        "mean_commutator_norm": mean([float(r.get("commutator_relative_norm_mean",0.0)) for r in rs]),
        "mean_rank": mean([float(r.get("operator_family_rank",0.0)) for r in rs]),
        "mean_comm_modes_added": mean([float(r.get("commutator_modes_added",0.0)) for r in rs]),
        "mean_second_modes_added": mean([float(r.get("second_order_modes_added",0.0)) for r in rs]),
        "used_delta_beta_any": False,
    }


def selected_configs(relax_steps:int)->List[dict]:
    keep={"real_growth_linear_star_candidate","saturating_growth_star_candidate","strict_symmetrized_response_star_control"}
    return [c for c in configs(relax_steps) if c["variant"] in keep]


def run_variant(cfg:dict,L:int,relax_steps:int,full_until:int,cap:int,outdir:Path)->Tuple[List[dict],List[dict]]:
    c=dict(cfg); c["relax_steps"]=relax_steps
    m=MatrixCollectingModel(**c, full_until_level=full_until, sample_cap=cap)
    m.run(L)
    rows=[]
    for gm in ["depth", "depth_boundarydim", "global_dim"]:
        for initial in ["depth_group_average", "principal_stable_modes", "mean_plus_principal"]:
            for comp in ["base", "commutator_completed_1", "commutator_completed_2"]:
                rows.extend(coarse_group_completion_metrics(m.matrix_rows, gm, initial, comp))
    prefix=f"L{L}_{cfg['variant']}_comm_completion"
    serial=[]
    for r in rows:
        serial.append({k:v for k,v in r.items() if not isinstance(v,(np.ndarray,list,dict))})
    write_csv(outdir/f"completion_rows_{prefix}.csv", serial)
    write_csv(outdir/f"events_{prefix}.csv", m.event_rows)
    write_csv(outdir/f"sampling_{prefix}.csv", m.sample_rows)
    return serial, m.sample_rows


def run_suite(levels:List[int], relax_steps:int, full_until:int, cap:int, outdir:Path)->dict:
    outdir.mkdir(parents=True, exist_ok=True)
    all_rows=[]; all_sampling=[]
    for L in levels:
        for cfg in selected_configs(relax_steps):
            rows,samp=run_variant(cfg,L,relax_steps,full_until,cap,outdir)
            for r in rows:
                r["max_level"]=L
                all_rows.append(r)
            for s in samp:
                s["variant"]=cfg["variant"]; s["max_level"]=L; all_sampling.append(s)
    write_csv(outdir/"all_completion_rows.csv", all_rows)
    write_csv(outdir/"sampling_rows_all.csv", all_sampling)
    summaries=[]
    for L in levels:
        variants=sorted({r.get("variant","") for r in all_rows if r.get("max_level")==L})
        coarse_levels=sorted({r.get("coarse_level","") for r in all_rows if r.get("max_level")==L})
        completion_levels=sorted({r.get("completion_level","") for r in all_rows if r.get("max_level")==L})
        for v in variants:
            for cl in coarse_levels:
                for cm in completion_levels:
                    rs=[r for r in all_rows if r.get("max_level")==L and r.get("variant")==v and r.get("coarse_level")==cl and r.get("completion_level")==cm]
                    if not rs: continue
                    summaries.append({"max_level":L, **summarize(rs,v,cl,cm)})
    write_csv(outdir/"summary_by_level_variant_completion.csv", summaries)
    summary={
        "model_family":"script1_script2_true_schur_live_commutator_generator_completion_gate",
        "levels":levels,
        "relax_steps":relax_steps,
        "full_until_level":full_until,
        "frontier_sample_cap_after_full_until":cap,
        "notes":[
            "Uses actual Record->Live true Schur/DtN operators, not scalar residual CSVs.",
            "Starts from coarse-grained families and adds only derived commutator modes [A_i,A_j].",
            "Tests base, one-round, and two-round commutator completion.",
            "Metric adjoint #_G is not counted alone; simultaneous product and commutator closure is required.",
            "No J, i, Hodge, physical Hilbert adjoint, C*-norm, Q/P target, or delta-beta gate is used."
        ],
        "summaries":summaries,
    }
    (outdir/"summary.json").write_text(json.dumps(summary,indent=2),encoding="utf-8")
    return summary


def make_md(summary:dict,outdir:Path):
    maxL=max(summary["levels"])
    rows=[s for s in summary["summaries"] if s["max_level"]==maxL]
    table="\n".join([
        f"| {s['variant']} | {s['coarse_level']} | {s['completion_level']} | {s['rows']} | {s['weak_pass_fraction']:.3f} | {s['strong_pass_fraction']:.3f} | {s['mean_product_resid']:.3f} | {s['mean_commutator_resid']:.3f} | {s['mean_commutator_max']:.3f} | {s['mean_comm_modes_added']:.2f} | {s['mean_rank']:.2f} |"
        for s in rows
    ])
    md=f"""# CNNA commutator-generator completion gate

## Purpose

This package tests whether the missing commutator closure of the coarse-grained live Record→Live Schur/DtN operator family can be repaired by adding **only derived commutator modes**.

It starts from coarse-grained families:

```text
raw not used as target; tested coarse families are:
depth_group_average
principal_stable_modes
mean_plus_principal
```

and compares:

```text
base
commutator_completed_1
commutator_completed_2
```

No arbitrary generators are added.  No `J`, `i`, Hodge, physical Hilbert adjoint, C*-norm, Q/P target, positivity axiom, or delta-beta decision is used.

## Deepest-level summary

| variant | coarse family | completion | rows | weak pass | strong pass | product resid | commutator resid | commutator max | comm modes added | rank |
|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
{table}

## Interpretation rule

A real success would require that commutator completion lowers the commutator residual **without** destroying product closure and # closure, and that strict-sym remains null.

If completion improves product closure but leaves commutator residual high, then the family is not yet a stable small real *-algebra-like operator family.  It may still be a real Lie-/operator precursor, but not a closed algebraic structure.

## Next methodological step

If the best row is still only weak/partial, the next test should not claim C*-structure.  It should inspect whether the residual commutator modes are depth-coherent, central/noise-like, or require a larger coarse-graining window.
"""
    (outdir/"RESULTS.md").write_text(md,encoding="utf-8")
    (outdir/"SUMMARY.md").write_text(md,encoding="utf-8")


def package_dir(pkgdir:Path,outzip:Path):
    if outzip.exists(): outzip.unlink()
    with zipfile.ZipFile(outzip,"w",zipfile.ZIP_DEFLATED) as zf:
        for p in sorted(pkgdir.rglob("*")):
            if p.is_dir(): continue
            zf.write(p,p.relative_to(pkgdir.parent))


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--levels",default="3,4")
    ap.add_argument("--relax-steps",type=int,default=4)
    ap.add_argument("--full-until-level",type=int,default=3)
    ap.add_argument("--frontier-sample-cap",type=int,default=8)
    ap.add_argument("--outdir",type=Path,default=Path("outputs"))
    args=ap.parse_args()
    levels=[int(x) for x in args.levels.split(',') if x]
    summary=run_suite(levels,args.relax_steps,args.full_until_level,args.frontier_sample_cap,args.outdir)
    make_md(summary,args.outdir)
    print(json.dumps({k:v for k,v in summary.items() if k!='summaries'},indent=2))

if __name__=="__main__":
    main()
