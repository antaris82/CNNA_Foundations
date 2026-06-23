#!/usr/bin/env python3
"""
CNNA Growth Test: parent-fan J-field gluing and holonomy.

Purpose
-------
This test follows the parent-fan plaquette skew-sector J-candidate diagnostic.
The previous test showed that the robust carrier

    K_abc = [A_ab, A_bc]

on provenance-generated parent-fan 2-simplices has a stable rank-2 skew sector
and yields a local diagnostic J = K/omega with J^2 ~= -I on the extracted
plane.

This test asks the next question:

    Do the local parent-fan J-candidates glue coherently across the face net?

For each completed parent p with ordered children c1,c2,c3, the oriented
parent-fan faces are

    F12=(p,c1,c2), F23=(p,c2,c3), F31=(p,c3,c1).

Each face yields a local skew K_F and normalized J_F on a 2D plane in the
common 3-port DtN coordinate space.  The test measures:

1. Per-face J validity and J^2 residual.
2. Pairwise plane/projector agreement between F12,F23,F31.
3. Pairwise J agreement up to sign and with fixed cyclic sign.
4. Pairwise J-commutators [J_F,J_G].
5. Minimal-axis connection holonomy around F12 -> F23 -> F31 -> F12.
6. Controls: diagonal/trace/common-commuting reductions, symmetrized birth,
   no-backreaction, sibling faces, and random same-level faces.

No global J, no complex Hilbert space, no C*-norm, no GNS representation and no
AQFT net are introduced here.  This is a finite diagnostic for whether the
local face-level J candidates behave like a coherent growth-defined J field.
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


def fro(M: np.ndarray) -> float:
    return float(np.linalg.norm(M, ord="fro"))


def nrm(v: np.ndarray) -> float:
    return float(np.linalg.norm(v))


def unit(v: np.ndarray) -> np.ndarray:
    n = nrm(v)
    if n < EPS:
        return np.zeros_like(v, dtype=float)
    return np.asarray(v, dtype=float) / n


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
    # K = [[0,-z,y],[z,0,-x],[-y,x,0]] corresponds to cross-product by (x,y,z).
    return np.array([K[2, 1], K[0, 2], K[1, 0]], dtype=float)


def skew_from_axis(n: np.ndarray) -> np.ndarray:
    x, y, z = np.asarray(n, dtype=float)
    return np.array([[0.0, -z, y], [z, 0.0, -x], [-y, x, 0.0]], dtype=float)


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


@dataclass(frozen=True)
class Face:
    kind: str
    parent: int
    level: int
    label: str
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
    rank2_balance: float
    kernel_gap: float
    first_norm: float
    edge_norm: float


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
    omega3 = math.sqrt(float(vals[2])) if vals.size > 2 else 0.0
    omega = 0.5 * (omega1 + omega2)
    rank2_balance = omega2 / (omega1 + EPS)
    kernel_gap = omega3 / (omega2 + EPS)
    axis = axial_from_skew(K)
    if omega <= min_omega or omega2 <= min_omega:
        return JFace(face, K, None, None, axis, omega, K_norm, False, float("nan"), float("nan"), rank2_balance, kernel_gap, first, edge)
    U2 = vecs[:, :2]
    P = U2 @ U2.T
    J = K / (omega + EPS)
    I = np.eye(3)
    J2_res = fro(P @ (J @ J + P) @ P) / (fro(P) + EPS)
    leakage = fro((I - P) @ J @ P) / (fro(J @ P) + EPS)
    return JFace(face, K, J, P, axis, omega, K_norm, True, J2_res, leakage, rank2_balance, kernel_gap, first, edge)


def plane_res(P: Optional[np.ndarray], Q: Optional[np.ndarray]) -> float:
    if P is None or Q is None:
        return float("nan")
    return fro(P - Q) / (fro(P) + fro(Q) + EPS)


def J_res(A: Optional[np.ndarray], B: Optional[np.ndarray], signed: bool) -> float:
    if A is None or B is None:
        return float("nan")
    if signed:
        return fro(A - B) / (fro(A) + fro(B) + EPS)
    return min(fro(A - B), fro(A + B)) / (fro(A) + fro(B) + EPS)


def J_comm(A: Optional[np.ndarray], B: Optional[np.ndarray]) -> float:
    if A is None or B is None:
        return float("nan")
    return fro(A @ B - B @ A)


def axis_cos(a: np.ndarray, b: np.ndarray, signless: bool = False) -> float:
    na = nrm(a)
    nb = nrm(b)
    if na < EPS or nb < EPS:
        return float("nan")
    c = float(np.dot(a, b) / (na * nb))
    return abs(c) if signless else c


def minimal_rotation(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    a = unit(a)
    b = unit(b)
    if nrm(a) < EPS or nrm(b) < EPS:
        return np.eye(3)
    v = np.cross(a, b)
    s = nrm(v)
    c = float(np.dot(a, b))
    if s < 1e-10:
        if c > 0:
            return np.eye(3)
        # 180-degree rotation around an arbitrary axis perpendicular to a.
        tmp = np.array([1.0, 0.0, 0.0])
        if abs(float(np.dot(tmp, a))) > 0.9:
            tmp = np.array([0.0, 1.0, 0.0])
        u = unit(np.cross(a, tmp))
        U = skew_from_axis(u)
        return np.eye(3) + 2.0 * (U @ U)
    vx = skew_from_axis(v)
    return np.eye(3) + vx + vx @ vx * ((1.0 - c) / (s * s + EPS))


def rotation_angle(R: np.ndarray) -> float:
    x = (float(np.trace(R)) - 1.0) / 2.0
    x = max(-1.0, min(1.0, x))
    return float(math.acos(x))


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


def build_parent_face_groups(model: RealGrowth, cells: Dict[int, LocalCell]) -> Dict[int, List[Face]]:
    out: Dict[int, List[Face]] = {}
    cell_ids = set(cells.keys())
    for p in sorted(model.completed_parent_ids()):
        children = model.child_ids_ordered(p)
        if len(children) != 3 or p not in cell_ids:
            continue
        c1, c2, c3 = children
        if not (c1 in cell_ids and c2 in cell_ids and c3 in cell_ids):
            continue
        lvl = int(model.nodes[p].level)
        out[p] = [
            Face("parent_fan_triangle", p, lvl, "F12", (p, c1, c2)),
            Face("parent_fan_triangle", p, lvl, "F23", (p, c2, c3)),
            Face("parent_fan_triangle", p, lvl, "F31", (p, c3, c1)),
        ]
    return out


def build_sibling_faces(model: RealGrowth, cells: Dict[int, LocalCell]) -> List[Face]:
    out: List[Face] = []
    cell_ids = set(cells.keys())
    for p in sorted(model.completed_parent_ids()):
        ch = model.child_ids_ordered(p)
        if len(ch) != 3:
            continue
        c1, c2, c3 = ch
        if c1 in cell_ids and c2 in cell_ids and c3 in cell_ids:
            out.append(Face("sibling_triangle", p, int(model.nodes[c1].level), "S123", (c1, c2, c3)))
    return out


def random_face_groups(model: RealGrowth, cells: Dict[int, LocalCell], count: int, seed: int) -> Dict[int, List[Face]]:
    rng = random.Random(seed)
    by_level: Dict[int, List[int]] = {}
    for cid, c in cells.items():
        by_level.setdefault(c.level, []).append(cid)
    levels = [lv for lv, xs in by_level.items() if len(xs) >= 9]
    out: Dict[int, List[Face]] = {}
    for i in range(count):
        lv = rng.choice(levels)
        verts = rng.sample(by_level[lv], 9)
        out[-100000 - i] = [
            Face("random_same_level_face_net", -100000 - i, lv, "R1", tuple(verts[0:3])),
            Face("random_same_level_face_net", -100000 - i, lv, "R2", tuple(verts[3:6])),
            Face("random_same_level_face_net", -100000 - i, lv, "R3", tuple(verts[6:9])),
        ]
    return out


def summarize_face_net(parent: int, faces: List[JFace], model: RealGrowth, source_mode: str, reduction: str, control: str) -> dict:
    labels = [f.face.label for f in faces]
    valid = [f for f in faces if f.valid]
    pairs = [(faces[0], faces[1]), (faces[1], faces[2]), (faces[2], faces[0])] if len(faces) == 3 else []
    pair_plane = [plane_res(a.P, b.P) for a, b in pairs]
    pair_J_fixed = [J_res(a.J, b.J, signed=True) for a, b in pairs]
    pair_J_signless = [J_res(a.J, b.J, signed=False) for a, b in pairs]
    pair_axis_cos = [axis_cos(a.axis, b.axis, False) for a, b in pairs]
    pair_axis_abs = [axis_cos(a.axis, b.axis, True) for a, b in pairs]
    pair_comm = [J_comm(a.J, b.J) for a, b in pairs]

    hol = float("nan")
    hol_abs = float("nan")
    if len(faces) == 3 and all(nrm(f.axis) > EPS for f in faces):
        n1, n2, n3 = (unit(f.axis) for f in faces)
        R12 = minimal_rotation(n1, n2)
        R23 = minimal_rotation(n2, n3)
        R31 = minimal_rotation(n3, n1)
        H = R31 @ R23 @ R12
        hol = rotation_angle(H)
        # Signless axes identify J and -J planes; use absolute alignment by flipping to previous axis.
        m1 = n1
        m2 = n2 if np.dot(m1, n2) >= 0 else -n2
        m3 = n3 if np.dot(m2, n3) >= 0 else -n3
        Q12 = minimal_rotation(m1, m2)
        Q23 = minimal_rotation(m2, m3)
        Q31 = minimal_rotation(m3, m1)
        hol_abs = rotation_angle(Q31 @ Q23 @ Q12)

    addr = ""
    if parent >= 0:
        addr = ".".join(map(str, address_tuple(model, parent)))
    return {
        "control": control,
        "source_mode": source_mode,
        "reduction": reduction,
        "face_net_kind": faces[0].face.kind if faces else "none",
        "parent": parent,
        "parent_addr": addr,
        "parent_level": int(model.nodes[parent].level) if parent >= 0 else int(faces[0].face.level),
        "face_labels": " ".join(labels),
        "faces": len(faces),
        "valid_faces": len(valid),
        "frac_valid_faces": len(valid) / max(len(faces), 1),
        "mean_K_norm": mean(f.K_norm for f in faces),
        "min_K_norm": min((f.K_norm for f in faces), default=float("nan")),
        "mean_omega": mean(f.omega for f in faces),
        "mean_J2_res": mean(f.J2_res for f in faces),
        "mean_leakage": mean(f.leakage for f in faces),
        "mean_rank2_balance": mean(f.rank2_balance for f in faces),
        "mean_kernel_gap": mean(f.kernel_gap for f in faces),
        "mean_first_norm": mean(f.first_norm for f in faces),
        "mean_edge_norm": mean(f.edge_norm for f in faces),
        "mean_plane_res": mean(pair_plane),
        "p95_plane_res": perc(pair_plane, 95),
        "mean_J_fixed_res": mean(pair_J_fixed),
        "mean_J_signless_res": mean(pair_J_signless),
        "mean_axis_cos": mean(pair_axis_cos),
        "mean_axis_abs_cos": mean(pair_axis_abs),
        "mean_J_pair_comm": mean(pair_comm),
        "holonomy_angle": hol,
        "holonomy_angle_signless": hol_abs,
    }


def summarize(rows: List[dict], keys: List[str]) -> List[dict]:
    groups: Dict[Tuple, List[dict]] = {}
    for r in rows:
        groups.setdefault(tuple(r[k] for k in keys), []).append(r)
    out: List[dict] = []
    for key, rs in sorted(groups.items(), key=lambda kv: kv[0]):
        d = {keys[i]: key[i] for i in range(len(keys))}
        d.update(
            count=len(rs),
            mean_valid=mean(r["frac_valid_faces"] for r in rs),
            mean_K_norm=mean(r["mean_K_norm"] for r in rs),
            mean_omega=mean(r["mean_omega"] for r in rs),
            mean_J2_res=mean(r["mean_J2_res"] for r in rs),
            mean_rank2_balance=mean(r["mean_rank2_balance"] for r in rs),
            mean_kernel_gap=mean(r["mean_kernel_gap"] for r in rs),
            mean_plane_res=mean(r["mean_plane_res"] for r in rs),
            p95_plane_res=perc((r["mean_plane_res"] for r in rs), 95),
            mean_J_fixed_res=mean(r["mean_J_fixed_res"] for r in rs),
            mean_J_signless_res=mean(r["mean_J_signless_res"] for r in rs),
            mean_axis_cos=mean(r["mean_axis_cos"] for r in rs),
            mean_axis_abs_cos=mean(r["mean_axis_abs_cos"] for r in rs),
            mean_J_pair_comm=mean(r["mean_J_pair_comm"] for r in rs),
            mean_holonomy_angle=mean(r["holonomy_angle"] for r in rs),
            mean_holonomy_angle_signless=mean(r["holonomy_angle_signless"] for r in rs),
            mean_first_norm=mean(r["mean_first_norm"] for r in rs),
        )
        out.append(d)
    return out


def run(args: argparse.Namespace) -> str:
    outdir: Path = args.outdir
    if outdir.exists():
        shutil.rmtree(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    controls = ["real_growth", "symmetrized_birth", "no_backreaction"]
    source_modes = ["live", "handoff", "aging"]
    reductions = ["full", "diagonal", "trace_scalar", "common_mean_diagonal"]
    face_rows: List[dict] = []
    net_rows: List[dict] = []
    model_rows: List[dict] = []

    for control in controls:
        model = build_model(control, args.max_level, args.mode)
        cells = build_cells(model, args.max_level)
        parent_groups = build_parent_face_groups(model, cells)
        rand_groups = random_face_groups(model, cells, args.random_nets, args.seed) if args.random_nets else {}
        sibling_faces = build_sibling_faces(model, cells)
        # Convert sibling faces into pseudo three-face nets by grouping consecutive sibling faces at same level.
        sibling_groups: Dict[int, List[Face]] = {}
        for i in range(0, min(len(sibling_faces) - 2, args.random_nets * 3), 3):
            sibling_groups[-200000 - i] = sibling_faces[i:i+3]

        model_rows.append({
            "control": control,
            "nodes": len(model.nodes),
            "completed_cells": len(cells),
            "parent_fan_nets": len(parent_groups),
            "sibling_pseudo_nets": len(sibling_groups),
            "random_nets": len(rand_groups),
        })
        all_groups = []
        for p, fs in parent_groups.items():
            all_groups.append(("parent_fan_net", p, fs))
        for p, fs in sibling_groups.items():
            all_groups.append(("sibling_pseudo_net", p, fs))
        for p, fs in rand_groups.items():
            all_groups.append(("random_same_level_net", p, fs))

        for source_mode in source_modes:
            for reduction in reductions:
                for net_kind, p, faces in all_groups:
                    jfaces = [j_from_face(f, cells, source_mode, reduction, args.min_omega) for f in faces]
                    for jf in jfaces:
                        face_rows.append({
                            "control": control,
                            "source_mode": source_mode,
                            "reduction": reduction,
                            "net_kind": net_kind,
                            "face_kind": jf.face.kind,
                            "parent": jf.face.parent,
                            "face_label": jf.face.label,
                            "face_level": jf.face.level,
                            "vertices": " ".join(map(str, jf.face.vertices)),
                            "K_norm": jf.K_norm,
                            "omega": jf.omega,
                            "valid": int(jf.valid),
                            "J2_res": jf.J2_res,
                            "leakage": jf.leakage,
                            "rank2_balance": jf.rank2_balance,
                            "kernel_gap": jf.kernel_gap,
                            "first_norm": jf.first_norm,
                            "edge_norm": jf.edge_norm,
                            "axis_x": float(jf.axis[0]),
                            "axis_y": float(jf.axis[1]),
                            "axis_z": float(jf.axis[2]),
                        })
                    r = summarize_face_net(p, jfaces, model, source_mode, reduction, control)
                    r["net_kind"] = net_kind
                    net_rows.append(r)

    main = summarize(net_rows, ["control", "net_kind", "source_mode", "reduction"])
    by_level = summarize(net_rows, ["control", "net_kind", "parent_level", "source_mode", "reduction"])
    write_csv(outdir / "J_field_face_rows.csv", face_rows)
    write_csv(outdir / "J_field_net_rows.csv", net_rows)
    write_csv(outdir / "J_field_summary_main.csv", main)
    write_csv(outdir / "J_field_summary_by_level.csv", by_level)
    write_csv(outdir / "model_summaries.csv", model_rows)

    def fmt(d: dict) -> str:
        return (
            f"{d['control']} | {d['net_kind']} | {d['source_mode']} | {d['reduction']}: "
            f"count={d['count']}, valid={d['mean_valid']:.3f}, K={d['mean_K_norm']:.6g}, "
            f"J2={d['mean_J2_res']:.3g}, plane={d['mean_plane_res']:.3g}, "
            f"J±={d['mean_J_signless_res']:.3g}, Jfixed={d['mean_J_fixed_res']:.3g}, "
            f"axis_abs={d['mean_axis_abs_cos']:.3f}, Jcomm={d['mean_J_pair_comm']:.3g}, "
            f"hol={d['mean_holonomy_angle']:.3g}, hol±={d['mean_holonomy_angle_signless']:.3g}"
        )

    lines = [
        "CNNA PARENT-FAN J-FIELD GLUING AND HOLONOMY DIAGNOSTIC",
        f"max_level={args.max_level}, mode={args.mode}, min_omega={args.min_omega}",
        "",
        "MODEL SUMMARIES",
    ]
    for m in model_rows:
        lines.append(f"  {m['control']}: nodes={m['nodes']}, completed={m['completed_cells']}, parent_nets={m['parent_fan_nets']}, sibling_pseudo={m['sibling_pseudo_nets']}, random={m['random_nets']}")
    lines.append("")
    lines.append("SELECTED FULL-DTN GLUING SUMMARIES")
    for d in main:
        if d["reduction"] == "full" and d["source_mode"] in ("live", "handoff", "aging"):
            lines.append("  " + fmt(d))
    lines.append("")
    lines.append("REAL-GROWTH LIVE ABLATIONS")
    for d in main:
        if d["control"] == "real_growth" and d["source_mode"] == "live" and d["net_kind"] == "parent_fan_net":
            lines.append("  " + fmt(d))
    lines.append("")
    lines.append("READING RULE")
    lines.extend([
        "  Local J validity was already checked per face in the previous test.",
        "  This test asks if same-parent face triples glue: small plane residuals,",
        "  small signless J residuals, high absolute axis cosine, and small J-commutators.",
        "  Holonomy here uses the minimal SO(3) connection between extracted axes;",
        "  small holonomy means coherent/flat local J-field gluing, not absence of face-level K.",
    ])
    summary = "\n".join(lines)
    (outdir / "SUMMARY.txt").write_text(summary, encoding="utf-8")

    results = f"""# Results: parent-fan J-field gluing and holonomy

This diagnostic follows the local parent-fan plaquette skew-sector J-candidate.
It does **not** derive a global complex structure.  It tests whether the local
face-level candidates glue coherently across the parent-fan face net.

## Current run

```text
{summary}
```

## Output files

- `J_field_face_rows.csv`
- `J_field_net_rows.csv`
- `J_field_summary_main.csv`
- `J_field_summary_by_level.csv`
- `model_summaries.csv`
- `SUMMARY.txt`

## Interpretation guide

A coherent local J-field carrier should show, for
`real_growth / parent_fan_net / live / full`:

```text
valid faces near 1
per-face J2 residual near 0
small plane residual across F12,F23,F31
small signless J residual across F12,F23,F31
high absolute axis cosine
small J-pair commutator
collapse or invalidation under diagonal/trace/common-commuting reductions
```

A small minimal-axis holonomy means the local J-field glues flatly across the
three parent-fan faces.  It does not erase the nonzero face-level commutator K;
it says the extracted J-planes themselves are mutually coherent.
"""
    (outdir / "RESULTS_parent_fan_J_field_gluing_and_holonomy.md").write_text(results, encoding="utf-8")

    with zipfile.ZipFile(outdir.with_suffix(".zip"), "w", compression=zipfile.ZIP_DEFLATED) as zf:
        for p in outdir.rglob("*"):
            if p.is_file():
                zf.write(p, p.relative_to(outdir.parent))

    return summary


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--max-level", type=int, default=6)
    ap.add_argument("--mode", default="linear")
    ap.add_argument("--min-omega", type=float, default=1e-10)
    ap.add_argument("--random-nets", type=int, default=60)
    ap.add_argument("--seed", type=int, default=91031)
    ap.add_argument("--outdir", type=Path, default=Path("parent_fan_J_field_gluing_out_L6"))
    args = ap.parse_args()
    print(run(args))


if __name__ == "__main__":
    main()
