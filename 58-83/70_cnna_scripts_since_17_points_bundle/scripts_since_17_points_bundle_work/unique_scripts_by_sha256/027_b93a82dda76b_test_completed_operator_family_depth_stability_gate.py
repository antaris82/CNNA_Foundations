#!/usr/bin/env python3
"""
CNNA completed operator-family depth stability gate.

Purpose
-------
Follow-up to the commutator-generator completion gate.  The previous test showed
that a coarse-grained Record->Live true-Schur/DtN operator family becomes much
more stable after adding only derived commutator modes.  This test asks whether
that completed family is stable across depth/generation and finite approximant
level.

No new algebraic structure is introduced.  It reuses the established script-1/2
true-Schur/DtN Record->Live growth model and the previous completion machinery.
No J, i, Hodge, physical Hilbert adjoint, C*-norm, Q/P target, positivity axiom,
or delta-beta gate is used.

Gate dimensions
---------------
1. level stability: L2,L3,L4 sampled/full-until-3 comparison
2. generation stability: parent_level grouped residuals/pass fractions
3. generator role stability: family rank and commutator modes added
4. strict_sym null

Important limitation
--------------------
Full true-Schur/DtN L4+ complete ternary growth is expensive in this environment.
The default uses full_until_level=3 and deterministic frontier sampling after
that, exactly like the preceding completion package.  The script marks sampling
rows explicitly and does not claim a full infinite-tree theorem.
"""
from __future__ import annotations

import argparse, json, zipfile, shutil
from pathlib import Path
from typing import Dict, List, Tuple
import numpy as np

from base_commutator_completion import run_suite, mean, median, write_csv


def slope(xs, ys):
    xs = np.asarray(xs, dtype=float); ys = np.asarray(ys, dtype=float)
    mask = np.isfinite(xs) & np.isfinite(ys)
    xs = xs[mask]; ys = ys[mask]
    if xs.size < 2 or np.allclose(xs, xs[0]):
        return 0.0
    return float(np.polyfit(xs, ys, 1)[0])


def load_rows(outdir: Path):
    import csv
    p = outdir / "all_completion_rows.csv"
    rows=[]
    if not p.exists(): return rows
    with p.open(newline='', encoding='utf-8') as f:
        for r in csv.DictReader(f):
            out={}
            for k,v in r.items():
                if v in ("", None): out[k]=v; continue
                if v in ("True","False"):
                    out[k]=(v=="True"); continue
                try:
                    if any(c in v for c in ['.','e','E']): out[k]=float(v)
                    else: out[k]=int(v)
                except Exception:
                    out[k]=v
            rows.append(out)
    return rows


def load_sampling(outdir: Path):
    import csv
    p = outdir / "sampling_rows_all.csv"
    rows=[]
    if not p.exists(): return rows
    with p.open(newline='', encoding='utf-8') as f:
        for r in csv.DictReader(f):
            out={}
            for k,v in r.items():
                if v in ("", None): out[k]=v; continue
                if v in ("True","False"):
                    out[k]=(v=="True"); continue
                try:
                    if any(c in v for c in ['.','e','E']): out[k]=float(v)
                    else: out[k]=int(v)
                except Exception:
                    out[k]=v
            rows.append(out)
    return rows


def grouped(rows: List[dict], keys: List[str]):
    d={}
    for r in rows:
        k=tuple(r.get(x) for x in keys)
        d.setdefault(k,[]).append(r)
    return d


def summarize_level(rows: List[dict], key: Tuple):
    # key: variant, coarse, completion, max_level
    return {
        "variant": key[0],
        "coarse_level": key[1],
        "completion_level": key[2],
        "max_level": key[3],
        "rows": len(rows),
        "weak_pass_fraction": mean([int(bool(r.get('weak_star_family_pass', False))) for r in rows]),
        "strong_pass_fraction": mean([int(bool(r.get('strong_stable_star_family_pass', False))) for r in rows]),
        "mean_product_resid": mean([float(r.get('product_closure_residual_mean', 0.0)) for r in rows]),
        "mean_product_max": mean([float(r.get('product_closure_residual_max', 0.0)) for r in rows]),
        "mean_commutator_resid": mean([float(r.get('commutator_closure_residual_mean', 0.0)) for r in rows]),
        "mean_commutator_max": mean([float(r.get('commutator_closure_residual_max', 0.0)) for r in rows]),
        "mean_rank": mean([float(r.get('operator_family_rank', 0.0)) for r in rows]),
        "mean_span_basis_count": mean([float(r.get('span_basis_count', 0.0)) for r in rows]),
        "mean_comm_modes_added": mean([float(r.get('commutator_modes_added', 0.0)) for r in rows]),
        "mean_second_modes_added": mean([float(r.get('second_order_modes_added', 0.0)) for r in rows]),
        "used_delta_beta_any": any(bool(r.get('used_delta_beta_any', False)) for r in rows),
    }


def level_stability(level_summaries: List[dict]):
    groups = grouped(level_summaries, ["variant", "coarse_level", "completion_level"])
    out=[]
    for (variant, coarse, completion), rs in sorted(groups.items(), key=lambda kv: str(kv[0])):
        rs=sorted(rs, key=lambda r: r["max_level"])
        if len(rs)<2: continue
        levels=[r["max_level"] for r in rs]
        strong=[float(r["strong_pass_fraction"]) for r in rs]
        weak=[float(r["weak_pass_fraction"]) for r in rs]
        prod=[float(r["mean_product_resid"]) for r in rs]
        comm=[float(r["mean_commutator_resid"]) for r in rs]
        rank=[float(r["mean_rank"]) for r in rs]
        last=rs[-1]
        prev=rs[-2]
        stable = (
            last["rows"] > 0 and
            last["strong_pass_fraction"] >= 0.45 and
            last["mean_product_resid"] < 0.02 and
            last["mean_commutator_resid"] < 0.02 and
            (last["strong_pass_fraction"] + 0.20 >= prev["strong_pass_fraction"]) and
            abs(slope(levels, comm)) < 0.03 and
            not last["used_delta_beta_any"]
        )
        out.append({
            "variant": variant,
            "coarse_level": coarse,
            "completion_level": completion,
            "levels": " ".join(map(str, levels)),
            "strong_pass_by_level": " ".join(f"{x:.4f}" for x in strong),
            "weak_pass_by_level": " ".join(f"{x:.4f}" for x in weak),
            "product_resid_by_level": " ".join(f"{x:.6g}" for x in prod),
            "comm_resid_by_level": " ".join(f"{x:.6g}" for x in comm),
            "rank_by_level": " ".join(f"{x:.4f}" for x in rank),
            "strong_pass_slope": slope(levels, strong),
            "comm_resid_slope": slope(levels, comm),
            "rank_slope": slope(levels, rank),
            "last_level": last["max_level"],
            "last_rows": last["rows"],
            "last_strong_pass_fraction": last["strong_pass_fraction"],
            "last_product_resid": last["mean_product_resid"],
            "last_commutator_resid": last["mean_commutator_resid"],
            "last_rank": last["mean_rank"],
            "level_stability_gate": bool(stable),
            "used_delta_beta_any": False,
        })
    return out


def generation_stability(rows: List[dict]):
    # evaluate parent_level spread at deepest level for each variant/coarse/completion
    maxL = max([int(r.get('max_level',0)) for r in rows], default=0)
    deep=[r for r in rows if int(r.get('max_level',0))==maxL]
    groups = grouped(deep, ["variant", "coarse_level", "completion_level", "parent_level"])
    parent_summ=[]
    for key, rs in sorted(groups.items(), key=lambda kv: str(kv[0])):
        if key[3] in (-1, "-1", ""): continue
        parent_summ.append({
            "variant": key[0], "coarse_level": key[1], "completion_level": key[2], "parent_level": key[3],
            "rows": len(rs),
            "strong_pass_fraction": mean([int(bool(r.get('strong_stable_star_family_pass', False))) for r in rs]),
            "weak_pass_fraction": mean([int(bool(r.get('weak_star_family_pass', False))) for r in rs]),
            "mean_product_resid": mean([float(r.get('product_closure_residual_mean', 0.0)) for r in rs]),
            "mean_commutator_resid": mean([float(r.get('commutator_closure_residual_mean', 0.0)) for r in rs]),
            "mean_rank": mean([float(r.get('operator_family_rank', 0.0)) for r in rs]),
            "used_delta_beta_any": False,
        })
    fam_groups=grouped(parent_summ,["variant","coarse_level","completion_level"])
    fam_out=[]
    for key, rs in sorted(fam_groups.items(), key=lambda kv: str(kv[0])):
        if len(rs)<2: continue
        pls=[float(r["parent_level"]) for r in rs]
        strong=[float(r["strong_pass_fraction"]) for r in rs]
        comm=[float(r["mean_commutator_resid"]) for r in rs]
        prod=[float(r["mean_product_resid"]) for r in rs]
        gate=(
            median(strong) >= 0.45 and
            max(comm) < 0.05 and
            max(prod) < 0.05 and
            not any(bool(r.get('used_delta_beta_any',False)) for r in rs)
        )
        fam_out.append({
            "variant": key[0], "coarse_level": key[1], "completion_level": key[2],
            "max_level": maxL,
            "parent_levels": " ".join(str(int(x)) for x in sorted(pls)),
            "parent_strong_pass": " ".join(f"{x:.4f}" for _,x in sorted(zip(pls,strong))),
            "parent_comm_resid": " ".join(f"{x:.6g}" for _,x in sorted(zip(pls,comm))),
            "parent_product_resid": " ".join(f"{x:.6g}" for _,x in sorted(zip(pls,prod))),
            "median_parent_strong_pass": median(strong),
            "max_parent_comm_resid": max(comm) if comm else 0.0,
            "max_parent_product_resid": max(prod) if prod else 0.0,
            "parent_strong_slope": slope(pls,strong),
            "parent_comm_slope": slope(pls,comm),
            "generation_stability_gate": bool(gate),
            "used_delta_beta_any": False,
        })
    return parent_summ, fam_out


def strict_sym_audit(level_summaries: List[dict]):
    strict=[s for s in level_summaries if 'strict_sym' in str(s.get('variant',''))]
    return {
        "strict_sym_rows": sum(int(s.get('rows',0)) for s in strict),
        "strict_sym_any_strong_pass": any(float(s.get('strong_pass_fraction',0.0))>0 for s in strict),
        "strict_sym_any_weak_pass": any(float(s.get('weak_pass_fraction',0.0))>0 for s in strict),
        "strict_sym_max_product_resid": max([float(s.get('mean_product_resid',0.0)) for s in strict], default=0.0),
        "strict_sym_max_comm_resid": max([float(s.get('mean_commutator_resid',0.0)) for s in strict], default=0.0),
        "strict_sym_null_gate": not any(float(s.get('strong_pass_fraction',0.0))>0 or float(s.get('weak_pass_fraction',0.0))>0 for s in strict),
        "used_delta_beta_any": False,
    }


def make_markdown(summary: dict, outdir: Path):
    stab=[r for r in summary["level_stability_rows"] if r.get("level_stability_gate")]
    gen=[r for r in summary["generation_stability_rows"] if r.get("generation_stability_gate")]
    # best rows by last comm residual among real growth and saturating
    ls=summary["level_stability_rows"]
    best=sorted([r for r in ls if 'strict_sym' not in r['variant']], key=lambda r: (r['last_commutator_resid'], -r['last_strong_pass_fraction']))[:10]
    best_table="\n".join(
        f"| {r['variant']} | {r['coarse_level']} | {r['completion_level']} | {r['last_level']} | {r['last_strong_pass_fraction']:.3f} | {r['last_product_resid']:.4g} | {r['last_commutator_resid']:.4g} | {r['rank_by_level']} | {r['level_stability_gate']} |"
        for r in best
    )
    gate_table="\n".join(
        f"| {r['variant']} | {r['coarse_level']} | {r['completion_level']} | {r['last_strong_pass_fraction']:.3f} | {r['last_product_resid']:.4g} | {r['last_commutator_resid']:.4g} |"
        for r in stab[:20]
    ) or "| none | | | | | |"
    gen_table="\n".join(
        f"| {r['variant']} | {r['coarse_level']} | {r['completion_level']} | {r['median_parent_strong_pass']:.3f} | {r['max_parent_product_resid']:.4g} | {r['max_parent_comm_resid']:.4g} | {r['parent_levels']} |"
        for r in gen[:20]
    ) or "| none | | | | | | |"
    md=f"""# CNNA completed operator-family depth stability gate

## Purpose

This package tests whether the commutator-completed, coarse-grained Record→Live true-Schur/DtN operator family remains stable across finite approximant depth and parent-generation.

It does **not** introduce `J`, `i`, Hodge, physical Hilbert adjoint, C*-norm, Q/P targets, positivity axioms, or delta-beta decisions.

## Model and limitation

The test reuses the established script-1/script-2 true-Schur/DtN growth line.  Full true-Schur/DtN growth is run up to `full_until_level = {summary['full_until_level']}`.  Above that level, frontier parents are selected by deterministic sampling with cap `{summary['frontier_sample_cap_after_full_until']}`.  Therefore L4 here is a sampled/frontier approximant if `full_until_level < 4`; this package does **not** claim a full infinite-tree theorem.

## Best deepest-level rows

| variant | coarse family | completion | level | strong pass | product resid | commutator resid | rank by level | stable? |
|---|---|---|---:|---:|---:|---:|---|---|
{best_table}

## Level-stability passes

| variant | coarse family | completion | last strong | last product | last commutator |
|---|---|---|---:|---:|---:|
{gate_table}

## Generation-stability passes at deepest level

| variant | coarse family | completion | median parent strong | max parent product | max parent comm | parent levels |
|---|---|---|---:|---:|---:|---|
{gen_table}

## Strict-sym audit

```json
{json.dumps(summary['strict_sym_audit'], indent=2)}
```

## Interpretation

A positive row means only this:

```text
The completed coarse-grained real operator family is stable under #_G, products,
and commutators across the finite/sampled depth audit.
```

It does **not** mean:

```text
C*-algebra, Hilbert positivity, J, i, or a physical quantum algebra.
```

If the best pass rows are present but rely on sampled L4, the correct next step is to harden the depth evidence with either a more efficient Schur/DtN implementation or a narrower full-L4 subset chosen by a predeclared rule.
"""
    (outdir/"RESULTS.md").write_text(md, encoding="utf-8")
    (outdir/"SUMMARY.md").write_text(md, encoding="utf-8")


def package(pkgdir: Path, outzip: Path):
    if outzip.exists(): outzip.unlink()
    with zipfile.ZipFile(outzip, 'w', zipfile.ZIP_DEFLATED) as zf:
        for p in sorted(pkgdir.rglob('*')):
            if p.is_dir(): continue
            zf.write(p, p.relative_to(pkgdir.parent))


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--levels', default='2,3,4')
    ap.add_argument('--relax-steps', type=int, default=4)
    ap.add_argument('--full-until-level', type=int, default=3)
    ap.add_argument('--frontier-sample-cap', type=int, default=8)
    ap.add_argument('--outdir', type=Path, default=Path('outputs'))
    ap.add_argument('--zip', type=Path, default=None)
    args=ap.parse_args()
    levels=[int(x) for x in args.levels.split(',') if x]
    args.outdir.mkdir(parents=True, exist_ok=True)
    # Run underlying suite into a subdir to avoid overwriting analysis products.
    rawdir=args.outdir/'raw_completion_suite'
    raw_summary=run_suite(levels, args.relax_steps, args.full_until_level, args.frontier_sample_cap, rawdir)
    raw_rows=load_rows(rawdir)
    sampling=load_sampling(rawdir)
    level_summ=[]
    for key, rs in grouped(raw_rows, ['variant','coarse_level','completion_level','max_level']).items():
        level_summ.append(summarize_level(rs,key))
    level_stab=level_stability(level_summ)
    parent_summ, gen_stab=generation_stability(raw_rows)
    strict=strict_sym_audit(level_summ)
    write_csv(args.outdir/'level_summary_rows.csv', level_summ)
    write_csv(args.outdir/'level_stability_rows.csv', level_stab)
    write_csv(args.outdir/'parent_level_summary_rows.csv', parent_summ)
    write_csv(args.outdir/'generation_stability_rows.csv', gen_stab)
    write_csv(args.outdir/'sampling_rows.csv', sampling)
    summary={
        'model_family':'script1_script2_true_schur_completed_operator_family_depth_stability_gate',
        'levels':levels,
        'relax_steps':args.relax_steps,
        'full_until_level':args.full_until_level,
        'frontier_sample_cap_after_full_until':args.frontier_sample_cap,
        'level_stability_gate_count':sum(int(r.get('level_stability_gate',False)) for r in level_stab),
        'generation_stability_gate_count':sum(int(r.get('generation_stability_gate',False)) for r in gen_stab),
        'strict_sym_audit':strict,
        'level_stability_rows':level_stab,
        'generation_stability_rows':gen_stab,
        'notes':[
            'Uses actual Record->Live true Schur/DtN operators through the previous completion machinery.',
            'Tests stability across L2,L3,L4 with full_until_level sampling explicitly logged.',
            'No J, i, Hodge, physical Hilbert adjoint, C*-norm, Q/P target, positivity axiom, or delta-beta gate is used.',
            'A pass is finite/sampled evidence for a small real #/*-like operator-family precursor only, not a C*-algebra claim.'
        ]
    }
    (args.outdir/'summary.json').write_text(json.dumps(summary, indent=2), encoding='utf-8')
    make_markdown(summary, args.outdir)
    print(json.dumps({k:v for k,v in summary.items() if k not in ('level_stability_rows','generation_stability_rows')}, indent=2))
    if args.zip:
        package(args.outdir.parent, args.zip)

if __name__=='__main__':
    main()
