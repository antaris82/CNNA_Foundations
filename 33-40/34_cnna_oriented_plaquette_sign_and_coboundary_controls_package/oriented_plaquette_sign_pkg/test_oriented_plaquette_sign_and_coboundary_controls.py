#!/usr/bin/env python3
"""
CNNA Growth Test: oriented plaquette sign and coboundary controls.

Purpose
-------
This is the immediate audit test after the first simplicial DtN plaquette
frustration diagnostic.  It does not extract J.  It checks whether the
face-level signal behaves like a genuine oriented plaquette effect rather than
like a first-order height/coboundary artefact.

Threads joined
--------------
- Tests 1/2: event-resolved sequential growth and backreaction.
- Test 3: provenance-generated 2-simplices.
- Test 9: geometry/minimal-axis holonomy alone is flat.
- Tests 22/23: dynamic Record/Live DtN operators.
- Tests 24-32: node/edge provenance operator closure alone remains too
  gradient-like.
- Previous bridge test: full DtN matrices on faces have nonzero second-order
  plaquette commutators while diagonal/trace controls collapse.

Gate questions
--------------
1. Orientation: Does reversing an oriented triangle send the commutator/skew
   signal to its negative, not merely preserve a norm?
2. Coboundary: Since A_ab = S_b-S_a, the first-order circulation is forced to
   vanish.  Is the remaining second-order signal explained away by scalar,
   diagonal, or common-commuting-potential ablations?
3. Face separation: Do sibling triangles and parent-fan triangles carry
   different signal profiles?

No i, no J, no Hilbert space, no C*-norm, and no *-algebra are introduced.
"""

from __future__ import annotations

import argparse
import csv
import math
import random
import shutil
import zipfile
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

import numpy as np

from record_live_block_base import RealGrowth, LocalCell, build_cells, write_csv

EPS = 1e-12


def fro(M: np.ndarray) -> float:
    return float(np.linalg.norm(M, ord="fro"))


def mean(xs: Iterable[float]) -> float:
    vals = [float(x) for x in xs if np.isfinite(float(x))]
    return float(np.mean(vals)) if vals else float("nan")


def perc(xs: Iterable[float], q: float) -> float:
    vals = [float(x) for x in xs if np.isfinite(float(x))]
    return float(np.percentile(vals, q)) if vals else float("nan")


def skew(M: np.ndarray) -> np.ndarray:
    return 0.5 * (M - M.T)


def sym(M: np.ndarray) -> np.ndarray:
    return 0.5 * (M + M.T)


def sym_expm(M: np.ndarray, step: float) -> np.ndarray:
    A = sym(np.asarray(M, dtype=float))
    vals, vecs = np.linalg.eigh(A)
    vals = np.clip(step * vals, -30.0, 30.0)
    return (vecs * np.exp(vals)) @ vecs.T


def max_imag_eig(M: np.ndarray) -> float:
    vals = np.linalg.eigvals(M)
    return float(np.max(np.abs(np.imag(vals))))


def inv_safe(M: np.ndarray) -> np.ndarray:
    return np.linalg.pinv(M, rcond=1e-10)


def address_tuple(model: RealGrowth, node_id: int) -> Tuple[int, ...]:
    out: List[int] = []
    cur = int(node_id)
    while model.nodes[cur].parent is not None:
        out.append(int(model.nodes[cur].birth_order))
        cur = int(model.nodes[cur].parent)
    return tuple(reversed(out))


@dataclass(frozen=True)
class OrientedFacePair:
    mode: str
    level: int
    base_parent: int
    forward: Tuple[int, int, int]
    reverse: Tuple[int, int, int]


def build_face_pairs(model: RealGrowth, cells: Dict[int, LocalCell], *, include_random: int, seed: int) -> List[OrientedFacePair]:
    pairs: List[OrientedFacePair] = []
    cell_ids = set(cells.keys())
    for p in sorted(model.completed_parent_ids()):
        children = model.child_ids_ordered(p)
        if len(children) != 3:
            continue
        c1, c2, c3 = children
        if c1 in cell_ids and c2 in cell_ids and c3 in cell_ids:
            level = int(model.nodes[c1].level)
            pairs.append(OrientedFacePair("sibling_triangle", level, p, (c1, c2, c3), (c1, c3, c2)))
        if p in cell_ids:
            for a, b in [(c1, c2), (c2, c3), (c3, c1)]:
                if a in cell_ids and b in cell_ids:
                    level = int(model.nodes[p].level)
                    pairs.append(OrientedFacePair("parent_fan_triangle", level, p, (p, a, b), (p, b, a)))

    if include_random > 0:
        rng = random.Random(seed)
        by_level: Dict[int, List[int]] = {}
        for cid, c in cells.items():
            by_level.setdefault(c.level, []).append(cid)
        levels = [lv for lv, xs in by_level.items() if len(xs) >= 3]
        for _ in range(include_random):
            lv = rng.choice(levels)
            tri = tuple(rng.sample(by_level[lv], 3))
            pairs.append(OrientedFacePair("random_same_level_triangle", lv, -1, tri, (tri[0], tri[2], tri[1])))
    return pairs


def raw_vertex_operator(cell: LocalCell, source_mode: str) -> np.ndarray:
    return sym(cell.matrix_for_source(source_mode))


def reduce_ops(Sa: np.ndarray, Sb: np.ndarray, Sc: np.ndarray, reduction: str) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    if reduction == "full":
        return Sa, Sb, Sc
    if reduction == "diagonal":
        return np.diag(np.diag(Sa)), np.diag(np.diag(Sb)), np.diag(np.diag(Sc))
    if reduction == "trace_scalar":
        return (
            np.eye(3) * (float(np.trace(Sa)) / 3.0),
            np.eye(3) * (float(np.trace(Sb)) / 3.0),
            np.eye(3) * (float(np.trace(Sc)) / 3.0),
        )
    if reduction == "common_mean_diagonal":
        M = sym((Sa + Sb + Sc) / 3.0)
        vals, Q = np.linalg.eigh(M)
        def project(S: np.ndarray) -> np.ndarray:
            D = np.diag(np.diag(Q.T @ S @ Q))
            return Q @ D @ Q.T
        return project(Sa), project(Sb), project(Sc)
    raise ValueError(reduction)


def face_core(vertices: Tuple[int, int, int], cells: Dict[int, LocalCell], source_mode: str, reduction: str, step: float) -> dict:
    a, b, c = vertices
    Sa0 = raw_vertex_operator(cells[a], source_mode)
    Sb0 = raw_vertex_operator(cells[b], source_mode)
    Sc0 = raw_vertex_operator(cells[c], source_mode)
    Sa, Sb, Sc = reduce_ops(Sa0, Sb0, Sc0, reduction)
    Aab = Sb - Sa
    Abc = Sc - Sb
    Aca = Sa - Sc
    first = Aab + Abc + Aca
    comm = Aab @ Abc - Abc @ Aab
    Uab = sym_expm(Aab, step)
    Ubc = sym_expm(Abc, step)
    Uca = sym_expm(Aca, step)
    P = Uca @ Ubc @ Uab
    I = np.eye(3)
    PmI = P - I
    pair_comm = mean([fro(Sa @ Sb - Sb @ Sa), fro(Sb @ Sc - Sc @ Sb), fro(Sc @ Sa - Sa @ Sc)])
    return {
        "S": (Sa, Sb, Sc),
        "A": (Aab, Abc, Aca),
        "P": P,
        "comm": comm,
        "skewP": skew(P),
        "first_norm": fro(first),
        "edge_mean_norm": mean([fro(Aab), fro(Abc), fro(Aca)]),
        "comm_norm": fro(comm),
        "plaquette_norm": fro(PmI),
        "skew_norm": fro(skew(P)),
        "skew_fraction": fro(skew(P)) / (fro(PmI) + EPS),
        "max_imag_eig": max_imag_eig(P),
        "has_complex_pair": int(max_imag_eig(P) > 1e-8),
        "pairwise_vertex_comm_norm": pair_comm,
    }


def pair_metrics(pair: OrientedFacePair, cells: Dict[int, LocalCell], source_mode: str, reduction: str, step: float) -> dict:
    f = face_core(pair.forward, cells, source_mode, reduction, step)
    r = face_core(pair.reverse, cells, source_mode, reduction, step)
    Csign = fro(f["comm"] + r["comm"]) / (fro(f["comm"]) + fro(r["comm"]) + EPS)
    Ssign = fro(f["skewP"] + r["skewP"]) / (fro(f["skewP"]) + fro(r["skewP"]) + EPS)
    try:
        invP = inv_safe(f["P"])
        inv_res = fro(r["P"] - invP) / (fro(invP - np.eye(3)) + fro(r["P"] - np.eye(3)) + EPS)
    except Exception:
        inv_res = float("nan")
    # sign alignment of axial-vector representatives of skew matrices
    # K -> (K32, K13, K21).  Reversal should produce negative vector.
    def axial(K: np.ndarray) -> np.ndarray:
        return np.array([K[2, 1], K[0, 2], K[1, 0]], dtype=float)
    af = axial(f["comm"])
    ar = axial(r["comm"])
    axial_cos_reversed = float(np.dot(af, -ar) / ((np.linalg.norm(af) * np.linalg.norm(ar)) + EPS))
    sf = axial(f["skewP"])
    sr = axial(r["skewP"])
    skew_axial_cos_reversed = float(np.dot(sf, -sr) / ((np.linalg.norm(sf) * np.linalg.norm(sr)) + EPS))
    return {
        "mode": pair.mode,
        "face_level": pair.level,
        "base_parent": pair.base_parent,
        "forward_vertices": " ".join(map(str, pair.forward)),
        "reverse_vertices": " ".join(map(str, pair.reverse)),
        "forward_addr": " | ".join(".".join(map(str, address_tuple(cells[v].model, v))) for v in pair.forward),
        "source_mode": source_mode,
        "reduction": reduction,
        "first_order_closure_norm_f": f["first_norm"],
        "first_order_closure_norm_r": r["first_norm"],
        "edge_generator_mean_norm": f["edge_mean_norm"],
        "commutator_norm_f": f["comm_norm"],
        "commutator_norm_r": r["comm_norm"],
        "commutator_sign_residual": Csign,
        "commutator_axial_cos_reversal": axial_cos_reversed,
        "plaquette_norm_f": f["plaquette_norm"],
        "plaquette_norm_r": r["plaquette_norm"],
        "plaquette_skew_norm_f": f["skew_norm"],
        "plaquette_skew_norm_r": r["skew_norm"],
        "plaquette_skew_sign_residual": Ssign,
        "plaquette_skew_axial_cos_reversal": skew_axial_cos_reversed,
        "plaquette_inverse_residual": inv_res,
        "skew_fraction_f": f["skew_fraction"],
        "max_imag_eig_f": f["max_imag_eig"],
        "has_complex_pair_f": f["has_complex_pair"],
        "pairwise_vertex_comm_norm": f["pairwise_vertex_comm_norm"],
    }


def summarize(rows: List[dict], keys: List[str]) -> List[dict]:
    groups: Dict[Tuple, List[dict]] = {}
    for r in rows:
        groups.setdefault(tuple(r[k] for k in keys), []).append(r)
    out: List[dict] = []
    for k, rs in sorted(groups.items(), key=lambda kv: kv[0]):
        d = {keys[i]: k[i] for i in range(len(keys))}
        d.update(
            count=len(rs),
            mean_first_order_closure=mean(r["first_order_closure_norm_f"] for r in rs),
            mean_edge_generator_norm=mean(r["edge_generator_mean_norm"] for r in rs),
            mean_pairwise_vertex_comm_norm=mean(r["pairwise_vertex_comm_norm"] for r in rs),
            mean_commutator_norm=mean(r["commutator_norm_f"] for r in rs),
            p95_commutator_norm=perc((r["commutator_norm_f"] for r in rs), 95),
            mean_commutator_sign_residual=mean(r["commutator_sign_residual"] for r in rs),
            mean_commutator_axial_cos_reversal=mean(r["commutator_axial_cos_reversal"] for r in rs),
            mean_plaquette_norm=mean(r["plaquette_norm_f"] for r in rs),
            p95_plaquette_norm=perc((r["plaquette_norm_f"] for r in rs), 95),
            mean_skew_norm=mean(r["plaquette_skew_norm_f"] for r in rs),
            mean_skew_sign_residual=mean(r["plaquette_skew_sign_residual"] for r in rs),
            mean_skew_axial_cos_reversal=mean(r["plaquette_skew_axial_cos_reversal"] for r in rs),
            mean_inverse_residual=mean(r["plaquette_inverse_residual"] for r in rs),
            mean_max_imag_eig=mean(r["max_imag_eig_f"] for r in rs),
            frac_complex_pair=mean(r["has_complex_pair_f"] for r in rs),
        )
        out.append(d)
    return out


def build_model(control: str, max_level: int, mode: str) -> RealGrowth:
    if control == "real_growth":
        model = RealGrowth(mode=mode)
    elif control == "symmetrized_birth":
        model = RealGrowth(mode=mode, growth_rule="symmetrized_birth")
    elif control == "no_backreaction":
        model = RealGrowth(mode=mode, br_ancestor=0.0, br_sibling=0.0)
    else:
        raise ValueError(control)
    model.grow(max_level)
    return model


def control_drop_rows(rows: List[dict]) -> List[dict]:
    # Compare full vs ablated reductions within same control/mode/source.
    groups: Dict[Tuple, Dict[str, List[dict]]] = {}
    for r in rows:
        key = (r["control"], r["mode"], r["source_mode"])
        groups.setdefault(key, {}).setdefault(r["reduction"], []).append(r)
    out: List[dict] = []
    for (control, mode, source), byred in sorted(groups.items()):
        full = byred.get("full", [])
        if not full:
            continue
        full_comm = mean(r["commutator_norm_f"] for r in full)
        full_skew = mean(r["plaquette_skew_norm_f"] for r in full)
        for red in ["diagonal", "trace_scalar", "common_mean_diagonal"]:
            rr = byred.get(red, [])
            if not rr:
                continue
            comm = mean(r["commutator_norm_f"] for r in rr)
            skewn = mean(r["plaquette_skew_norm_f"] for r in rr)
            out.append({
                "control": control,
                "mode": mode,
                "source_mode": source,
                "ablation": red,
                "full_commutator_norm": full_comm,
                "ablation_commutator_norm": comm,
                "commutator_remaining_fraction": comm / (full_comm + EPS),
                "full_skew_norm": full_skew,
                "ablation_skew_norm": skewn,
                "skew_remaining_fraction": skewn / (full_skew + EPS),
            })
    return out


def run(args: argparse.Namespace) -> str:
    outdir: Path = args.outdir
    if outdir.exists():
        shutil.rmtree(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    controls = ["real_growth", "symmetrized_birth", "no_backreaction"]
    source_modes = ["record", "live", "handoff", "aging"]
    reductions = ["full", "diagonal", "trace_scalar", "common_mean_diagonal"]
    rows: List[dict] = []
    model_rows: List[dict] = []

    for control in controls:
        model = build_model(control, args.max_level, args.mode)
        cells = build_cells(model, args.max_level)
        pairs = build_face_pairs(model, cells, include_random=args.random_faces, seed=args.seed)
        model_rows.append({
            "control": control,
            "nodes": len(model.nodes),
            "completed_cells": len(cells),
            "face_pairs_total": len(pairs),
            "sibling_pairs": sum(1 for p in pairs if p.mode == "sibling_triangle"),
            "parent_fan_pairs": sum(1 for p in pairs if p.mode == "parent_fan_triangle"),
            "random_pairs": sum(1 for p in pairs if p.mode == "random_same_level_triangle"),
        })
        for source_mode in source_modes:
            for reduction in reductions:
                for pair in pairs:
                    r = pair_metrics(pair, cells, source_mode, reduction, args.step)
                    r["control"] = control
                    r["max_level"] = args.max_level
                    r["step"] = args.step
                    rows.append(r)

    by_main = summarize(rows, ["control", "mode", "source_mode", "reduction"])
    by_face_level = summarize(rows, ["control", "mode", "face_level", "source_mode", "reduction"])
    by_source = summarize(rows, ["control", "source_mode", "reduction"])
    ablation_rows = control_drop_rows(rows)

    write_csv(outdir / "oriented_pair_rows.csv", rows)
    write_csv(outdir / "orientation_summary_main.csv", by_main)
    write_csv(outdir / "orientation_summary_by_face_level.csv", by_face_level)
    write_csv(outdir / "orientation_summary_by_source.csv", by_source)
    write_csv(outdir / "coboundary_ablation_summary.csv", ablation_rows)
    write_csv(outdir / "model_summaries.csv", model_rows)

    def fmt(d: dict) -> str:
        return (
            f"{d['control']} | {d.get('mode','*')} | {d.get('source_mode','*')} | {d.get('reduction','*')}: "
            f"count={d['count']}, comm={d['mean_commutator_norm']:.6g}, "
            f"comm_sign={d['mean_commutator_sign_residual']:.3g}, "
            f"comm_cos={d['mean_commutator_axial_cos_reversal']:.3f}, "
            f"plaq={d['mean_plaquette_norm']:.6g}, skew={d['mean_skew_norm']:.6g}, "
            f"skew_sign={d['mean_skew_sign_residual']:.3g}, "
            f"inv={d['mean_inverse_residual']:.3g}, complex_frac={d['frac_complex_pair']:.3f}, "
            f"first={d['mean_first_order_closure']:.2e}"
        )

    lines = [
        "CNNA ORIENTED PLAQUETTE SIGN AND COBOUNDARY CONTROLS",
        f"max_level={args.max_level}, mode={args.mode}, step={args.step}",
        "",
        "MODEL SUMMARIES",
    ]
    for m in model_rows:
        lines.append(f"  {m['control']}: nodes={m['nodes']}, completed={m['completed_cells']}, sibling={m['sibling_pairs']}, parent_fan={m['parent_fan_pairs']}, random={m['random_pairs']}")
    lines.append("")
    lines.append("SELECTED FULL-DTN FACE SUMMARIES")
    for d in by_main:
        if d["reduction"] == "full" and d["source_mode"] in ("live", "handoff", "aging"):
            lines.append("  " + fmt(d))
    lines.append("")
    lines.append("REAL-GROWTH LIVE COBOUNDARY/COMMUTING ABLATIONS")
    for d in ablation_rows:
        if d["control"] == "real_growth" and d["source_mode"] == "live":
            lines.append(
                f"  {d['mode']} | {d['ablation']}: comm_remain={d['commutator_remaining_fraction']:.3g}, "
                f"skew_remain={d['skew_remaining_fraction']:.3g}, "
                f"full_comm={d['full_commutator_norm']:.6g}, abl_comm={d['ablation_commutator_norm']:.6g}"
            )
    lines.append("")
    lines.append("INTERPRETATION")
    lines.extend([
        "  Good orientation behavior means commutator_sign_residual and skew_sign_residual are near 0,",
        "  and axial_cos_reversal is near +1 when comparing forward with negative reversed orientation.",
        "  Diagonal/trace/common-basis ablations test whether the signal is merely scalar/commuting coboundary data.",
        "  first_order_closure is forced to telescope and should remain numerical zero.",
        "  This is still not a J-test; it audits the face-level noncommutative carrier before any J extraction.",
    ])
    summary = "\n".join(lines)
    (outdir / "SUMMARY.txt").write_text(summary, encoding="utf-8")

    results = f"""# Results: oriented plaquette sign and coboundary controls

This test audits the first simplicial DtN plaquette-frustration signal.

It still does **not** set `i`, `J`, a Hilbert space, a C*-norm, or a *-algebra.

## Gate questions

1. Orientation reversal: `(a,b,c)` versus `(a,c,b)` should flip the commutator/skew sign.
2. Coboundary controls: first-order edge circulation is exact and telescopes; diagonal, trace-scalar, and common-mean-diagonal ablations should remove any scalar/commuting explanation.
3. Face separation: sibling triangles and parent-fan triangles are kept as distinct geometric carriers.

## Current run

```text
{summary}
```

## Output files

- `oriented_pair_rows.csv`
- `orientation_summary_main.csv`
- `orientation_summary_by_face_level.csv`
- `orientation_summary_by_source.csv`
- `coboundary_ablation_summary.csv`
- `model_summaries.csv`
- `SUMMARY.txt`

## Reading rule

A useful positive sign audit is:

```text
commutator_sign_residual ~ 0
plaquette_skew_sign_residual ~ 0
commutator_axial_cos_reversal ~ 1
first_order_closure ~ 0
full DtN signal >> diagonal/trace/common-commuting ablations
```

If these fail, the prior plaquette signal is not yet reliable enough for any J-sector projection.
"""
    (outdir / "RESULTS_oriented_plaquette_sign_and_coboundary_controls.md").write_text(results, encoding="utf-8")
    return summary


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--max-level", type=int, default=6)
    ap.add_argument("--mode", choices=["linear", "log", "saturating"], default="linear")
    ap.add_argument("--step", type=float, default=0.04)
    ap.add_argument("--random-faces", type=int, default=80)
    ap.add_argument("--seed", type=int, default=1731)
    ap.add_argument("--outdir", type=Path, default=Path("oriented_plaquette_sign_out_L6"))
    args = ap.parse_args()
    print(run(args))


if __name__ == "__main__":
    main()
