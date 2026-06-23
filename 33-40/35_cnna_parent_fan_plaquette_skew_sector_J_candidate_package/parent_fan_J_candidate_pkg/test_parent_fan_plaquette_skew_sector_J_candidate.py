#!/usr/bin/env python3
"""
CNNA Growth Test: parent-fan plaquette skew-sector J candidate.

Purpose
-------
This is the first deliberately J-near test after the simplicial DtN
plaquette-frustration and oriented-sign audits.

It does not introduce a complex Hilbert space, a C*-norm, a GNS
representation, or an AQFT net.  It asks a narrower local question:

    Does the robust parent-fan plaquette commutator
        K_abc = [A_ab, A_bc]
    carry a stable rank-2 skew sector that behaves like a real local
    rotation plane?

For symmetric full DtN vertex operators S_a,S_b,S_c, the edge differences
A_ab=S_b-S_a, A_bc=S_c-S_b are symmetric; hence K=[A_ab,A_bc] is skew.
A nonzero 3x3 skew operator has a 1D kernel and a 2D image plane.  On that
plane, after normalization by its nonzero singular value omega, it should
satisfy J^2=-I_plane.

This is still a diagnostic.  A nonzero skew matrix automatically determines a
formal plane rotation.  The nontrivial gates are therefore:

1. The skew signal must come from the robust real-growth parent-fan full-DtN
   carrier, not from diagonal/scalar/commuting controls.
2. Orientation reversal must send J -> -J on the same plane.
3. Symmetrized-birth should suppress the stable J candidate.
4. Sibling triangles are tracked separately.
5. No-backreaction is separated into birth-environment/live contribution and
   handoff/aging contribution.
6. The rank-2 plane/axis should be nondegenerate and measurable across faces.

No free J is set; J is only normalized from the measured commutator K.
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


def max_imag_eig(M: np.ndarray) -> float:
    vals = np.linalg.eigvals(M)
    return float(np.max(np.abs(np.imag(vals))))


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


def commutator_for_vertices(vertices: Tuple[int, int, int], cells: Dict[int, LocalCell], source_mode: str, reduction: str) -> dict:
    a, b, c = vertices
    Sa0 = raw_vertex_operator(cells[a], source_mode)
    Sb0 = raw_vertex_operator(cells[b], source_mode)
    Sc0 = raw_vertex_operator(cells[c], source_mode)
    Sa, Sb, Sc = reduce_ops(Sa0, Sb0, Sc0, reduction)
    Aab = Sb - Sa
    Abc = Sc - Sb
    Aca = Sa - Sc
    K = skew(Aab @ Abc - Abc @ Aab)
    return {
        "S": (Sa, Sb, Sc),
        "A": (Aab, Abc, Aca),
        "K": K,
        "first_norm": fro(Aab + Abc + Aca),
        "edge_mean_norm": mean([fro(Aab), fro(Abc), fro(Aca)]),
        "vertex_pair_comm_norm": mean([fro(Sa @ Sb - Sb @ Sa), fro(Sb @ Sc - Sc @ Sb), fro(Sc @ Sa - Sa @ Sc)]),
    }


def axial_from_skew(K: np.ndarray) -> np.ndarray:
    # For skew matrix [[0,-z,y],[z,0,-x],[-y,x,0]], axial=(x,y,z).
    return np.array([K[2, 1], K[0, 2], K[1, 0]], dtype=float)


def j_sector_metrics(K: np.ndarray, *, min_omega: float) -> dict:
    K = skew(np.asarray(K, dtype=float))
    normK = fro(K)
    vals, vecs = np.linalg.eigh(-K @ K)
    order = np.argsort(vals)[::-1]
    vals = np.maximum(vals[order], 0.0)
    vecs = vecs[:, order]
    omega1 = math.sqrt(float(vals[0])) if vals.size else 0.0
    omega2 = math.sqrt(float(vals[1])) if vals.size > 1 else 0.0
    omega3 = math.sqrt(float(vals[2])) if vals.size > 2 else 0.0
    omega = 0.5 * (omega1 + omega2)
    rank2_balance = omega2 / (omega1 + EPS)
    kernel_gap = omega3 / (omega2 + EPS)
    if omega <= min_omega or omega2 <= min_omega:
        return {
            "K_norm": normK,
            "omega": omega,
            "omega1": omega1,
            "omega2": omega2,
            "omega3": omega3,
            "rank2_balance": rank2_balance,
            "kernel_gap": kernel_gap,
            "J_valid": 0,
            "J2_plane_residual": float("nan"),
            "J_leakage": float("nan"),
            "plane_projector_trace": 0.0,
            "axis_norm": float(np.linalg.norm(axial_from_skew(K))),
            "max_imag_eig_J": float("nan"),
        }
    U2 = vecs[:, :2]
    P = U2 @ U2.T
    I = np.eye(K.shape[0])
    J = K / (omega + EPS)
    # On the extracted plane, J^2 should be -P and J should preserve the plane.
    J2_res = fro(P @ (J @ J + P) @ P) / (fro(P) + EPS)
    leakage = fro((I - P) @ J @ P) / (fro(J @ P) + EPS)
    max_imag = max_imag_eig(P @ J @ P)
    return {
        "K_norm": normK,
        "omega": omega,
        "omega1": omega1,
        "omega2": omega2,
        "omega3": omega3,
        "rank2_balance": rank2_balance,
        "kernel_gap": kernel_gap,
        "J_valid": 1,
        "J2_plane_residual": J2_res,
        "J_leakage": leakage,
        "plane_projector_trace": float(np.trace(P)),
        "axis_norm": float(np.linalg.norm(axial_from_skew(K))),
        "max_imag_eig_J": max_imag,
    }


def plane_projector_from_K(K: np.ndarray, min_omega: float) -> np.ndarray | None:
    vals, vecs = np.linalg.eigh(-skew(K) @ skew(K))
    order = np.argsort(vals)[::-1]
    vals = np.maximum(vals[order], 0.0)
    vecs = vecs[:, order]
    if len(vals) < 2 or math.sqrt(float(vals[1])) <= min_omega:
        return None
    U2 = vecs[:, :2]
    return U2 @ U2.T


def oriented_pair_row(pair: OrientedFacePair, cells: Dict[int, LocalCell], source_mode: str, reduction: str, min_omega: float) -> dict:
    f = commutator_for_vertices(pair.forward, cells, source_mode, reduction)
    r = commutator_for_vertices(pair.reverse, cells, source_mode, reduction)
    Kf = f["K"]
    Kr = r["K"]
    mf = j_sector_metrics(Kf, min_omega=min_omega)
    mr = j_sector_metrics(Kr, min_omega=min_omega)
    af = axial_from_skew(Kf)
    ar = axial_from_skew(Kr)
    axial_cos_reversal = float(np.dot(af, -ar) / ((np.linalg.norm(af) * np.linalg.norm(ar)) + EPS))
    sign_res = fro(Kf + Kr) / (fro(Kf) + fro(Kr) + EPS)
    # Compare forward plane with reversed plane.  They should be same plane, opposite orientation.
    Pf = plane_projector_from_K(Kf, min_omega)
    Pr = plane_projector_from_K(Kr, min_omega)
    if Pf is not None and Pr is not None:
        plane_projector_res = fro(Pf - Pr) / (fro(Pf) + fro(Pr) + EPS)
    else:
        plane_projector_res = float("nan")
    # J sign residual after normalizing by each omega.
    if mf["J_valid"] and mr["J_valid"]:
        Jf = Kf / (mf["omega"] + EPS)
        Jr = Kr / (mr["omega"] + EPS)
        J_sign_res = fro(Jf + Jr) / (fro(Jf) + fro(Jr) + EPS)
    else:
        J_sign_res = float("nan")
    return {
        "mode": pair.mode,
        "face_level": pair.level,
        "base_parent": pair.base_parent,
        "forward_vertices": " ".join(map(str, pair.forward)),
        "reverse_vertices": " ".join(map(str, pair.reverse)),
        "forward_addr": " | ".join(".".join(map(str, address_tuple(cells[v].model, v))) for v in pair.forward),
        "source_mode": source_mode,
        "reduction": reduction,
        "first_order_closure_norm": f["first_norm"],
        "edge_generator_mean_norm": f["edge_mean_norm"],
        "vertex_pair_comm_norm": f["vertex_pair_comm_norm"],
        "K_norm": mf["K_norm"],
        "omega": mf["omega"],
        "omega1": mf["omega1"],
        "omega2": mf["omega2"],
        "omega3": mf["omega3"],
        "rank2_balance": mf["rank2_balance"],
        "kernel_gap": mf["kernel_gap"],
        "J_valid": mf["J_valid"],
        "J2_plane_residual": mf["J2_plane_residual"],
        "J_leakage": mf["J_leakage"],
        "plane_projector_trace": mf["plane_projector_trace"],
        "max_imag_eig_J": mf["max_imag_eig_J"],
        "orientation_K_sign_residual": sign_res,
        "orientation_J_sign_residual": J_sign_res,
        "orientation_axial_cos_reversal": axial_cos_reversal,
        "orientation_plane_projector_residual": plane_projector_res,
        "reverse_K_norm": mr["K_norm"],
        "reverse_J2_plane_residual": mr["J2_plane_residual"],
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
            mean_first_order_closure=mean(r["first_order_closure_norm"] for r in rs),
            mean_edge_generator_norm=mean(r["edge_generator_mean_norm"] for r in rs),
            mean_vertex_pair_comm_norm=mean(r["vertex_pair_comm_norm"] for r in rs),
            mean_K_norm=mean(r["K_norm"] for r in rs),
            p95_K_norm=perc((r["K_norm"] for r in rs), 95),
            mean_omega=mean(r["omega"] for r in rs),
            mean_rank2_balance=mean(r["rank2_balance"] for r in rs),
            mean_kernel_gap=mean(r["kernel_gap"] for r in rs),
            frac_J_valid=mean(r["J_valid"] for r in rs),
            mean_J2_plane_residual=mean(r["J2_plane_residual"] for r in rs),
            p95_J2_plane_residual=perc((r["J2_plane_residual"] for r in rs), 95),
            mean_J_leakage=mean(r["J_leakage"] for r in rs),
            mean_max_imag_eig_J=mean(r["max_imag_eig_J"] for r in rs),
            mean_orientation_K_sign_residual=mean(r["orientation_K_sign_residual"] for r in rs),
            mean_orientation_J_sign_residual=mean(r["orientation_J_sign_residual"] for r in rs),
            mean_orientation_axial_cos_reversal=mean(r["orientation_axial_cos_reversal"] for r in rs),
            mean_orientation_plane_projector_residual=mean(r["orientation_plane_projector_residual"] for r in rs),
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


def ablation_rows(rows: List[dict]) -> List[dict]:
    groups: Dict[Tuple, Dict[str, List[dict]]] = {}
    for r in rows:
        key = (r["control"], r["mode"], r["source_mode"])
        groups.setdefault(key, {}).setdefault(r["reduction"], []).append(r)
    out: List[dict] = []
    for (control, mode, source), byred in sorted(groups.items()):
        full = byred.get("full", [])
        if not full:
            continue
        full_K = mean(r["K_norm"] for r in full)
        full_valid = mean(r["J_valid"] for r in full)
        for red in ["diagonal", "trace_scalar", "common_mean_diagonal"]:
            rr = byred.get(red, [])
            if not rr:
                continue
            k = mean(r["K_norm"] for r in rr)
            out.append({
                "control": control,
                "mode": mode,
                "source_mode": source,
                "ablation": red,
                "full_K_norm": full_K,
                "ablation_K_norm": k,
                "K_remaining_fraction": k / (full_K + EPS),
                "full_frac_J_valid": full_valid,
                "ablation_frac_J_valid": mean(r["J_valid"] for r in rr),
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
                    r = oriented_pair_row(pair, cells, source_mode, reduction, args.min_omega)
                    r["control"] = control
                    r["max_level"] = args.max_level
                    r["min_omega"] = args.min_omega
                    rows.append(r)

    by_main = summarize(rows, ["control", "mode", "source_mode", "reduction"])
    by_face_level = summarize(rows, ["control", "mode", "face_level", "source_mode", "reduction"])
    by_source = summarize(rows, ["control", "source_mode", "reduction"])
    abls = ablation_rows(rows)

    write_csv(outdir / "J_candidate_face_rows.csv", rows)
    write_csv(outdir / "J_candidate_summary_main.csv", by_main)
    write_csv(outdir / "J_candidate_summary_by_face_level.csv", by_face_level)
    write_csv(outdir / "J_candidate_summary_by_source.csv", by_source)
    write_csv(outdir / "J_candidate_ablation_summary.csv", abls)
    write_csv(outdir / "model_summaries.csv", model_rows)

    def fmt(d: dict) -> str:
        return (
            f"{d['control']} | {d.get('mode','*')} | {d.get('source_mode','*')} | {d.get('reduction','*')}: "
            f"count={d['count']}, K={d['mean_K_norm']:.6g}, omega={d['mean_omega']:.6g}, "
            f"J_valid={d['frac_J_valid']:.3f}, J2={d['mean_J2_plane_residual']:.3g}, "
            f"leak={d['mean_J_leakage']:.3g}, rank2={d['mean_rank2_balance']:.6g}, "
            f"kgap={d['mean_kernel_gap']:.3g}, Jsign={d['mean_orientation_J_sign_residual']:.3g}, "
            f"Ksign={d['mean_orientation_K_sign_residual']:.3g}, axis_cos={d['mean_orientation_axial_cos_reversal']:.3f}, "
            f"plane_rev={d['mean_orientation_plane_projector_residual']:.3g}, first={d['mean_first_order_closure']:.2e}"
        )

    lines = [
        "CNNA PARENT-FAN PLAQUETTE SKEW-SECTOR J-CANDIDATE DIAGNOSTIC",
        f"max_level={args.max_level}, mode={args.mode}, min_omega={args.min_omega}",
        "",
        "MODEL SUMMARIES",
    ]
    for m in model_rows:
        lines.append(f"  {m['control']}: nodes={m['nodes']}, completed={m['completed_cells']}, sibling={m['sibling_pairs']}, parent_fan={m['parent_fan_pairs']}, random={m['random_pairs']}")
    lines.append("")
    lines.append("SELECTED FULL-DTN J-SECTOR SUMMARIES")
    for d in by_main:
        if d["reduction"] == "full" and d["source_mode"] in ("live", "handoff", "aging"):
            lines.append("  " + fmt(d))
    lines.append("")
    lines.append("REAL-GROWTH LIVE FULL VS COMMUTING ABLATIONS")
    for d in abls:
        if d["control"] == "real_growth" and d["source_mode"] == "live":
            lines.append(
                f"  {d['mode']} | {d['ablation']}: K_remain={d['K_remaining_fraction']:.3g}, "
                f"Jvalid_full={d['full_frac_J_valid']:.3f}, Jvalid_abl={d['ablation_frac_J_valid']:.3f}, "
                f"full_K={d['full_K_norm']:.6g}, abl_K={d['ablation_K_norm']:.6g}"
            )
    lines.append("")
    lines.append("INTERPRETATION RULE")
    lines.extend([
        "  A useful local J-candidate carrier requires K_norm and omega above threshold,",
        "  rank2_balance near 1, kernel_gap near 0, J2_plane_residual near 0,",
        "  orientation_J_sign_residual near 0 under reversed orientation, and collapse under",
        "  diagonal/trace/common-commuting reductions.  These are necessary diagnostic gates,",
        "  not yet a derived complex algebra or global J field.",
    ])
    summary = "\n".join(lines)
    (outdir / "SUMMARY.txt").write_text(summary, encoding="utf-8")

    results = f"""# Results: parent-fan plaquette skew-sector J-candidate diagnostic

This is the first J-near test after the face-level plaquette-frustration and
orientation/coboundary audits.

It still does **not** introduce a complex Hilbert space, a C*-norm, GNS, AQFT,
or a global complex algebra.  The only candidate is extracted from the measured
real skew commutator

```text
K_abc = [A_ab, A_bc]
```

on oriented provenance-generated 2-simplices.

For each nonzero 3x3 skew `K`, the test extracts its rank-2 image plane and
normalizes `J = K / omega` on that plane.  The nontrivial issue is not the
linear algebra fact itself, but whether the robust parent-fan/live/full-DtN
carrier survives controls and has the right orientation behavior.

## Current run

```text
{summary}
```

## Output files

- `J_candidate_face_rows.csv`
- `J_candidate_summary_main.csv`
- `J_candidate_summary_by_face_level.csv`
- `J_candidate_summary_by_source.csv`
- `J_candidate_ablation_summary.csv`
- `model_summaries.csv`
- `SUMMARY.txt`

## Reading rule

Positive diagnostic signs are:

```text
real_growth / parent_fan_triangle / live / full:
  K_norm, omega nonzero
  frac_J_valid near 1
  rank2_balance near 1
  kernel_gap near 0
  J2_plane_residual near 0
  orientation_J_sign_residual near 0
  diagonal/trace/common-commuting ablations collapse the signal
```

Even if these hold, the output is still a local skew-sector J-candidate, not a
reelle *-Algebra, not a C*-completion, and not a GNS/AQFT net.
"""
    (outdir / "RESULTS_parent_fan_plaquette_skew_sector_J_candidate.md").write_text(results, encoding="utf-8")
    return summary


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--max-level", type=int, default=6)
    ap.add_argument("--mode", choices=["linear", "log", "saturating"], default="linear")
    ap.add_argument("--random-faces", type=int, default=80)
    ap.add_argument("--seed", type=int, default=1731)
    ap.add_argument("--min-omega", type=float, default=1e-10)
    ap.add_argument("--outdir", type=Path, default=Path("parent_fan_plaquette_J_candidate_out_L6"))
    args = ap.parse_args()
    print(run(args))


if __name__ == "__main__":
    main()
