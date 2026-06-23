#!/usr/bin/env python3
"""
CNNA Growth Test: double-history parent-fan plaquette J-holonomy.

Purpose
-------
This is the first test that combines the two previously separated carriers:

1. parent-fan 2-simplices, where full DtN plaquette commutators produce a
   local rank-2 skew sector and hence a local diagnostic J_F;
2. double-history suffix identification, where different first-root histories
   reach the same suffix and expose non-telescoping DtN handoff defects.

Previous tests showed that minimal SO(3) transport between extracted J axes is
flat/signless.  This test does not use that synthetic transport.  Instead it
uses the provenance quotient itself:

    same parent-level + same suffix after forgetting the first root sector,
    with directed root-sector order 1 -> 2 -> 3 -> 1.

For each quotient class containing all three root sectors, it compares the
three parent-fan J-nets under several label-gluing laws.  The primary law is
"directed_cyclic_root_shift": an edge from root sector r to s shifts the
face-label index by (s-r mod 3).  This is directed by birth/root order, not by a
best-fit SO(3) rotation.

The gate is whether the local Z3 parent-fan sectors remain nontrivial under
this provenance gluing, while controls vanish or collapse:

- identical-history control must vanish;
- diagonal/trace/common-commuting reductions must kill J;
- symmetrized birth should suppress K/J amplitudes;
- random quotient cycles should not mimic the same structured result.

No physical i, no global J-field, no C*-norm, no GNS representation and no AQFT
net are introduced here.  This is a finite diagnostic for whether the local
parent-fan J sectors have a non-flat provenance-identification holonomy.
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
from typing import Dict, Iterable, List, Optional, Tuple

import numpy as np

from record_live_block_base import RealGrowth, LocalCell, build_cells, write_csv

EPS = 1e-12
LABELS = ("F12", "F23", "F31")


def fro(M: np.ndarray) -> float:
    return float(np.linalg.norm(M, ord="fro"))


def nrm(v: np.ndarray) -> float:
    return float(np.linalg.norm(v))


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


def address_tuple(model: RealGrowth, node_id: int) -> Tuple[int, ...]:
    out: List[int] = []
    cur = int(node_id)
    while model.nodes[cur].parent is not None:
        out.append(int(model.nodes[cur].birth_order))
        cur = int(model.nodes[cur].parent)
    return tuple(reversed(out))


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


def axial_from_skew(K: np.ndarray) -> np.ndarray:
    return np.array([K[2, 1], K[0, 2], K[1, 0]], dtype=float)


@dataclass(frozen=True)
class Face:
    label: str
    parent: int
    level: int
    vertices: Tuple[int, int, int]


@dataclass
class JFace:
    face: Face
    K: np.ndarray
    J: Optional[np.ndarray]
    P: Optional[np.ndarray]
    axis: np.ndarray
    omega: float
    K_norm: float
    valid: bool
    J2_res: float
    leakage: float
    first_norm: float
    edge_norm: float


@dataclass
class JNet:
    parent: int
    root_sector: int
    suffix: Tuple[int, ...]
    level: int
    address: Tuple[int, ...]
    faces: Dict[str, JFace]


def local_commutator(vertices: Tuple[int, int, int], cells: Dict[int, LocalCell], source_mode: str, reduction: str) -> Tuple[np.ndarray, float, float]:
    a, b, c = vertices
    Sa0 = raw_vertex_operator(cells[a], source_mode)
    Sb0 = raw_vertex_operator(cells[b], source_mode)
    Sc0 = raw_vertex_operator(cells[c], source_mode)
    Sa, Sb, Sc = reduce_ops(Sa0, Sb0, Sc0, reduction)
    Aab = Sb - Sa
    Abc = Sc - Sb
    Aca = Sa - Sc
    K = skew(Aab @ Abc - Abc @ Aab)
    first = fro(Aab + Abc + Aca)
    edge = mean([fro(Aab), fro(Abc), fro(Aca)])
    return K, first, edge


def j_from_face(face: Face, cells: Dict[int, LocalCell], source_mode: str, reduction: str, min_omega: float) -> JFace:
    K, first, edge = local_commutator(face.vertices, cells, source_mode, reduction)
    K = skew(K)
    K_norm = fro(K)
    vals, vecs = np.linalg.eigh(-K @ K)
    order = np.argsort(vals)[::-1]
    vals = np.maximum(vals[order], 0.0)
    vecs = vecs[:, order]
    omega1 = math.sqrt(float(vals[0])) if vals.size > 0 else 0.0
    omega2 = math.sqrt(float(vals[1])) if vals.size > 1 else 0.0
    omega = 0.5 * (omega1 + omega2)
    axis = axial_from_skew(K)
    if omega <= min_omega or omega2 <= min_omega:
        return JFace(face, K, None, None, axis, omega, K_norm, False, float("nan"), float("nan"), first, edge)
    U2 = vecs[:, :2]
    P = U2 @ U2.T
    J = K / (omega + EPS)
    I = np.eye(3)
    J2_res = fro(P @ (J @ J + P) @ P) / (fro(P) + EPS)
    leakage = fro((I - P) @ J @ P) / (fro(J @ P) + EPS)
    return JFace(face, K, J, P, axis, omega, K_norm, True, J2_res, leakage, first, edge)


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


def build_parent_faces(model: RealGrowth, cells: Dict[int, LocalCell], parent: int) -> Optional[List[Face]]:
    cell_ids = set(cells.keys())
    if parent not in cell_ids:
        return None
    ch = model.child_ids_ordered(parent)
    if len(ch) != 3:
        return None
    c1, c2, c3 = ch
    if not (c1 in cell_ids and c2 in cell_ids and c3 in cell_ids):
        return None
    lvl = int(model.nodes[parent].level)
    return [
        Face("F12", parent, lvl, (parent, c1, c2)),
        Face("F23", parent, lvl, (parent, c2, c3)),
        Face("F31", parent, lvl, (parent, c3, c1)),
    ]


def build_jnets(model: RealGrowth, cells: Dict[int, LocalCell], source_mode: str, reduction: str, min_omega: float) -> Dict[int, JNet]:
    out: Dict[int, JNet] = {}
    for p in sorted(model.completed_parent_ids()):
        addr = address_tuple(model, p)
        if len(addr) < 1:
            continue
        faces = build_parent_faces(model, cells, p)
        if faces is None:
            continue
        jfaces = {f.label: j_from_face(f, cells, source_mode, reduction, min_omega) for f in faces}
        root_sector = int(addr[0])
        suffix = tuple(addr[1:])
        out[p] = JNet(p, root_sector, suffix, int(model.nodes[p].level), addr, jfaces)
    return out


def group_double_history(nets: Dict[int, JNet]) -> Dict[Tuple[int, Tuple[int, ...]], List[JNet]]:
    groups: Dict[Tuple[int, Tuple[int, ...]], List[JNet]] = {}
    for n in nets.values():
        groups.setdefault((n.level, n.suffix), []).append(n)
    return {k: sorted(v, key=lambda x: x.root_sector) for k, v in groups.items() if len({x.root_sector for x in v}) >= 2}


def root_index(root_sector: int) -> int:
    return (int(root_sector) - 1) % 3


def shift_for_law(src_root: int, dst_root: int, law: str, rng: Optional[random.Random] = None) -> int:
    delta = (root_index(dst_root) - root_index(src_root)) % 3
    if law == "directed_cyclic_root_shift":
        return delta
    if law == "same_label":
        return 0
    if law == "reverse_cyclic_root_shift":
        return (-delta) % 3
    if law == "random_label_shift":
        if rng is None:
            return 0
        return rng.randrange(3)
    raise ValueError(law)


def J_res(A: Optional[np.ndarray], B: Optional[np.ndarray], signed: bool) -> float:
    if A is None or B is None:
        return float("nan")
    if signed:
        return fro(A - B) / (fro(A) + fro(B) + EPS)
    return min(fro(A - B), fro(A + B)) / (fro(A) + fro(B) + EPS)


def best_sign(A: Optional[np.ndarray], B: Optional[np.ndarray]) -> int:
    if A is None or B is None:
        return 0
    return 1 if fro(A - B) <= fro(A + B) else -1


def axis_cos(a: np.ndarray, b: np.ndarray, signless: bool = False) -> float:
    na = nrm(a)
    nb = nrm(b)
    if na < EPS or nb < EPS:
        return float("nan")
    c = float(np.dot(a, b) / (na * nb))
    return abs(c) if signless else c


def label_index(label: str) -> int:
    return LABELS.index(label)


def label_shift(label: str, shift: int) -> str:
    return LABELS[(label_index(label) + shift) % 3]


def compare_edge(src: JNet, dst: JNet, law: str, rng: Optional[random.Random]) -> List[dict]:
    sh = shift_for_law(src.root_sector, dst.root_sector, law, rng)
    rows: List[dict] = []
    for lab in LABELS:
        dlab = label_shift(lab, sh)
        a = src.faces[lab]
        b = dst.faces[dlab]
        rows.append({
            "src_parent": src.parent,
            "dst_parent": dst.parent,
            "src_root": src.root_sector,
            "dst_root": dst.root_sector,
            "src_label": lab,
            "dst_label": dlab,
            "shift": sh,
            "valid_src": int(a.valid),
            "valid_dst": int(b.valid),
            "both_valid": int(a.valid and b.valid),
            "K_src": a.K_norm,
            "K_dst": b.K_norm,
            "J_fixed_res": J_res(a.J, b.J, signed=True),
            "J_signless_res": J_res(a.J, b.J, signed=False),
            "best_sign": best_sign(a.J, b.J),
            "axis_cos": axis_cos(a.axis, b.axis, False),
            "axis_abs_cos": axis_cos(a.axis, b.axis, True),
            "J2_src": a.J2_res,
            "J2_dst": b.J2_res,
        })
    return rows


def three_sector_cycle(nets: List[JNet]) -> Optional[List[JNet]]:
    by_root: Dict[int, JNet] = {n.root_sector: n for n in nets}
    if not all(r in by_root for r in (1, 2, 3)):
        return None
    return [by_root[1], by_root[2], by_root[3]]


def summarize_cycle(cycle: List[JNet], law: str, control: str, source_mode: str, reduction: str, class_key: Tuple[int, Tuple[int, ...]], rng: Optional[random.Random], mode: str) -> Tuple[dict, List[dict]]:
    edges = [(cycle[0], cycle[1]), (cycle[1], cycle[2]), (cycle[2], cycle[0])]
    edge_rows: List[dict] = []
    for ei, (a, b) in enumerate(edges):
        rs = compare_edge(a, b, law, rng)
        for r in rs:
            r.update(edge_index=ei, mode=mode, control=control, source_mode=source_mode, reduction=reduction, gluing_law=law, level=class_key[0], suffix=".".join(map(str, class_key[1])))
        edge_rows.extend(rs)

    fixed = [r["J_fixed_res"] for r in edge_rows]
    signless = [r["J_signless_res"] for r in edge_rows]
    both_valid = [r["both_valid"] for r in edge_rows]
    # Track Z2 sign products around the 1->2->3->1 cycle for each starting face label.
    z2_products = []
    z2_obstruct = []
    track_fixed = []
    track_signless = []
    for start in LABELS:
        lab = start
        prod = 1
        ok = True
        tf = []
        ts = []
        for (src, dst) in edges:
            sh = shift_for_law(src.root_sector, dst.root_sector, law, rng if law == "random_label_shift" else None)
            dlab = label_shift(lab, sh)
            a = src.faces[lab]
            b = dst.faces[dlab]
            s = best_sign(a.J, b.J)
            if s == 0 or not (a.valid and b.valid):
                ok = False
            else:
                prod *= s
                tf.append(J_res(a.J, b.J, signed=True))
                ts.append(J_res(a.J, b.J, signed=False))
            lab = dlab
        if ok:
            z2_products.append(prod)
            z2_obstruct.append(1 if prod < 0 else 0)
            track_fixed.append(mean(tf))
            track_signless.append(mean(ts))
    row = {
        "mode": mode,
        "control": control,
        "source_mode": source_mode,
        "reduction": reduction,
        "gluing_law": law,
        "level": class_key[0],
        "suffix": ".".join(map(str, class_key[1])),
        "parents": " ".join(str(n.parent) for n in cycle),
        "addresses": " | ".join(".".join(map(str, n.address)) for n in cycle),
        "valid_edge_frac": mean(both_valid),
        "mean_K": mean([r["K_src"] for r in edge_rows] + [r["K_dst"] for r in edge_rows]),
        "mean_J_fixed_res": mean(fixed),
        "mean_J_signless_res": mean(signless),
        "p95_J_signless_res": perc(signless, 95),
        "mean_axis_abs_cos": mean(r["axis_abs_cos"] for r in edge_rows),
        "mean_axis_cos": mean(r["axis_cos"] for r in edge_rows),
        "z2_tracks": len(z2_products),
        "z2_obstruction_frac": mean(z2_obstruct),
        "z2_product_mean": mean(z2_products),
        "track_fixed_res": mean(track_fixed),
        "track_signless_res": mean(track_signless),
    }
    return row, edge_rows


def summarize(rows: List[dict], keys: List[str]) -> List[dict]:
    groups: Dict[Tuple, List[dict]] = {}
    for r in rows:
        groups.setdefault(tuple(r[k] for k in keys), []).append(r)
    out: List[dict] = []
    for key, rs in sorted(groups.items(), key=lambda kv: kv[0]):
        d = {keys[i]: key[i] for i in range(len(keys))}
        d.update(
            count=len(rs),
            mean_valid_edge_frac=mean(r["valid_edge_frac"] for r in rs),
            mean_K=mean(r["mean_K"] for r in rs),
            mean_J_fixed_res=mean(r["mean_J_fixed_res"] for r in rs),
            mean_J_signless_res=mean(r["mean_J_signless_res"] for r in rs),
            p95_J_signless_res=perc((r["mean_J_signless_res"] for r in rs), 95),
            mean_axis_abs_cos=mean(r["mean_axis_abs_cos"] for r in rs),
            mean_axis_cos=mean(r["mean_axis_cos"] for r in rs),
            mean_z2_obstruction_frac=mean(r["z2_obstruction_frac"] for r in rs),
            mean_z2_product=mean(r["z2_product_mean"] for r in rs),
            mean_track_fixed_res=mean(r["track_fixed_res"] for r in rs),
            mean_track_signless_res=mean(r["track_signless_res"] for r in rs),
        )
        out.append(d)
    return out


def random_cycles(nets: Dict[int, JNet], count: int, seed: int) -> List[List[JNet]]:
    rng = random.Random(seed)
    buckets: Dict[int, List[JNet]] = {}
    for n in nets.values():
        buckets.setdefault(n.level, []).append(n)
    levels = [lv for lv, xs in buckets.items() if len(xs) >= 3]
    out: List[List[JNet]] = []
    attempts = 0
    while len(out) < count and attempts < max(2000, 40 * count):
        attempts += 1
        lv = rng.choice(levels)
        tri = rng.sample(buckets[lv], 3)
        # assign artificial sorted root sectors for shift law by sorting existing roots/address
        tri = sorted(tri, key=lambda n: (n.root_sector, n.parent))
        if len({n.suffix for n in tri}) == 1:
            continue
        out.append(tri)
    return out


def identical_cycles(nets: Dict[int, JNet], count: int) -> List[List[JNet]]:
    out = []
    for n in sorted(nets.values(), key=lambda x: (x.level, x.parent)):
        out.append([n, n, n])
        if len(out) >= count:
            break
    return out


def run(args: argparse.Namespace) -> str:
    outdir: Path = args.outdir
    if outdir.exists():
        shutil.rmtree(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    rng = random.Random(args.seed)

    controls = ["real_growth", "symmetrized_birth", "no_backreaction"]
    source_modes = ["live", "handoff", "aging"]
    reductions = ["full", "diagonal", "trace_scalar", "common_mean_diagonal"]
    laws = ["directed_cyclic_root_shift", "same_label", "reverse_cyclic_root_shift", "random_label_shift"]

    cycle_rows: List[dict] = []
    edge_rows_all: List[dict] = []
    model_rows: List[dict] = []

    for control in controls:
        model = build_model(control, args.max_level, args.mode)
        cells = build_cells(model, args.max_level)
        for source_mode in source_modes:
            for reduction in reductions:
                nets = build_jnets(model, cells, source_mode, reduction, args.min_omega)
                groups = group_double_history(nets)
                dh_cycles: List[Tuple[Tuple[int, Tuple[int, ...]], List[JNet]]] = []
                for key, gs in sorted(groups.items(), key=lambda kv: (kv[0][0], kv[0][1])):
                    cyc = three_sector_cycle(gs)
                    if cyc is not None:
                        dh_cycles.append((key, cyc))
                if source_mode == "live" and reduction == "full":
                    model_rows.append({
                        "control": control,
                        "nodes": len(model.nodes),
                        "cells": len(cells),
                        "jnets": len(nets),
                        "double_history_classes": len(groups),
                        "three_sector_cycles": len(dh_cycles),
                    })
                rand = random_cycles(nets, args.random_cycles, args.seed + 17) if args.random_cycles else []
                ident = identical_cycles(nets, min(args.identical_cycles, len(nets)))
                for law in laws:
                    for key, cyc in dh_cycles:
                        row, edges = summarize_cycle(cyc, law, control, source_mode, reduction, key, rng, "double_history_parent_fan_cycle")
                        cycle_rows.append(row)
                        edge_rows_all.extend(edges)
                    # Controls only for primary law and same_label to avoid huge CSVs.
                    if law in ("directed_cyclic_root_shift", "same_label"):
                        for i, cyc in enumerate(ident):
                            key = (cyc[0].level, cyc[0].suffix)
                            row, edges = summarize_cycle(cyc, law, control, source_mode, reduction, key, rng, "identical_history_control")
                            row["suffix"] = f"identity:{cyc[0].parent}"
                            cycle_rows.append(row)
                            edge_rows_all.extend(edges)
                        for i, cyc in enumerate(rand):
                            key = (cyc[0].level, tuple([-999, i]))
                            row, edges = summarize_cycle(cyc, law, control, source_mode, reduction, key, rng, "random_same_level_cycle_baseline")
                            cycle_rows.append(row)
                            edge_rows_all.extend(edges)

    main = summarize(cycle_rows, ["mode", "control", "source_mode", "reduction", "gluing_law"])
    by_level = summarize(cycle_rows, ["mode", "control", "level", "source_mode", "reduction", "gluing_law"])

    write_csv(outdir / "double_history_J_cycle_rows.csv", cycle_rows)
    write_csv(outdir / "double_history_J_edge_rows.csv", edge_rows_all)
    write_csv(outdir / "double_history_J_summary_main.csv", main)
    write_csv(outdir / "double_history_J_summary_by_level.csv", by_level)
    write_csv(outdir / "model_summaries.csv", model_rows)

    def findrow(mode: str, control: str, source: str, red: str, law: str) -> Optional[dict]:
        for d in main:
            if d["mode"] == mode and d["control"] == control and d["source_mode"] == source and d["reduction"] == red and d["gluing_law"] == law:
                return d
        return None

    selected_keys = [
        ("double_history_parent_fan_cycle", "real_growth", "live", "full", "directed_cyclic_root_shift"),
        ("double_history_parent_fan_cycle", "real_growth", "live", "full", "same_label"),
        ("double_history_parent_fan_cycle", "real_growth", "live", "full", "reverse_cyclic_root_shift"),
        ("double_history_parent_fan_cycle", "symmetrized_birth", "live", "full", "directed_cyclic_root_shift"),
        ("double_history_parent_fan_cycle", "real_growth", "handoff", "full", "directed_cyclic_root_shift"),
        ("identical_history_control", "real_growth", "live", "full", "directed_cyclic_root_shift"),
        ("random_same_level_cycle_baseline", "real_growth", "live", "full", "directed_cyclic_root_shift"),
        ("double_history_parent_fan_cycle", "real_growth", "live", "diagonal", "directed_cyclic_root_shift"),
        ("double_history_parent_fan_cycle", "real_growth", "live", "trace_scalar", "directed_cyclic_root_shift"),
    ]

    def fmt(d: Optional[dict]) -> str:
        if d is None:
            return "missing"
        return (
            f"{d['mode']} | {d['control']} | {d['source_mode']} | {d['reduction']} | {d['gluing_law']}: "
            f"count={d['count']}, valid={d['mean_valid_edge_frac']:.3f}, K={d['mean_K']:.6g}, "
            f"Jfixed={d['mean_J_fixed_res']:.4g}, J±={d['mean_J_signless_res']:.4g}, "
            f"axis_abs={d['mean_axis_abs_cos']:.4f}, z2obs={d['mean_z2_obstruction_frac']:.4g}, "
            f"z2prod={d['mean_z2_product']:.4g}, track±={d['mean_track_signless_res']:.4g}"
        )

    lines = [
        "CNNA DOUBLE-HISTORY PARENT-FAN PLAQUETTE J-HOLONOMY DIAGNOSTIC",
        f"max_level={args.max_level}, mode={args.mode}, min_omega={args.min_omega}",
        "",
        "MODEL SUMMARIES",
    ]
    for m in model_rows:
        lines.append(f"  {m['control']}: nodes={m['nodes']}, cells={m['cells']}, jnets={m['jnets']}, double_classes={m['double_history_classes']}, three_sector_cycles={m['three_sector_cycles']}")
    lines.append("")
    lines.append("SELECTED SUMMARIES")
    for key in selected_keys:
        lines.append("  " + fmt(findrow(*key)))
    lines.append("")
    lines.append("READING RULE")
    lines.extend([
        "  directed_cyclic_root_shift uses root-sector birth order 1->2->3->1 and shifts",
        "  parent-fan face labels by the directed sector difference. This is the primary",
        "  provenance-identification test; it is not a minimal SO(3) best-fit transport.",
        "  identical_history_control should have zero/near-zero residuals and no obstruction.",
        "  diagonal/trace/common-commuting reductions should invalidate/kill K.",
        "  symmetrized_birth should suppress the K amplitude and any robust holonomy signal.",
        "  A nonzero z2obs would indicate a sign-lift obstruction around the root-sector",
        "  double-history cycle; a small J± with high axis_abs but large Jfixed indicates",
        "  signless plane gluing without fixed-sign locking.",
    ])
    summary = "\n".join(lines)
    (outdir / "SUMMARY.txt").write_text(summary, encoding="utf-8")

    results = f"""# Results: double-history parent-fan plaquette J-holonomy

This test combines two previously separate structures:

1. parent-fan plaquette J-candidates from full DtN face commutators;
2. double-history suffix identification across root sectors.

It deliberately avoids minimal SO(3) axis transport.  The primary gluing law is
provenance-directed: root sector 1 -> 2 -> 3 -> 1, with face labels shifted by
the directed root-sector difference.

## Current run

```text
{summary}
```

## Output files

- `double_history_J_cycle_rows.csv`
- `double_history_J_edge_rows.csv`
- `double_history_J_summary_main.csv`
- `double_history_J_summary_by_level.csv`
- `model_summaries.csv`
- `SUMMARY.txt`

## Interpretation status

This is still not a derivation of a global complex structure.  It tests whether
the local parent-fan Z3/J sectors remain nontrivial when glued by the same
suffix-forgetting double-history quotient that previously exposed DtN handoff
defects.

A positive Stufe-4 candidate would require more than local J-validity: it would
need a double-history-specific sign/holonomy obstruction that vanishes for
identical history, collapses under symmetrized birth and commuting reductions,
and is not reproduced by random same-level cycles.
"""
    (outdir / "RESULTS_double_history_parent_fan_plaquette_J_holonomy.md").write_text(results, encoding="utf-8")

    # Package full directory.
    zip_path = outdir.with_suffix(".zip")
    if zip_path.exists():
        zip_path.unlink()
    with zipfile.ZipFile(zip_path, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        for p in outdir.rglob("*"):
            if p.is_file():
                zf.write(p, p.relative_to(outdir.parent))
    return summary


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--max-level", type=int, default=6)
    ap.add_argument("--mode", default="linear")
    ap.add_argument("--min-omega", type=float, default=1e-10)
    ap.add_argument("--random-cycles", type=int, default=121)
    ap.add_argument("--identical-cycles", type=int, default=121)
    ap.add_argument("--seed", type=int, default=20260621)
    ap.add_argument("--outdir", type=Path, default=Path("double_history_parent_fan_J_holonomy_out_L6"))
    args = ap.parse_args()
    print(run(args))


if __name__ == "__main__":
    main()
