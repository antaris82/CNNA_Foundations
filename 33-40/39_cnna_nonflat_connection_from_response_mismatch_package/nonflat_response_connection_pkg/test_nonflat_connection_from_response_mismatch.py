#!/usr/bin/env python3
"""
CNNA Growth Test: nonflat connection from same-label response mismatch.

Purpose
-------
The previous boundary-resolved correspondence test inferred identity/same-label
face gluing between double-history parent-fan nets. This test accepts that
finding and asks the next question:

    If same-label gluing is the natural correspondence, does the remaining
    DtN/response mismatch define a nonflat connection around double-history
    cycles?

The test does not impose a cyclic face-label shift and does not use the minimal
SO(3) axis transport from earlier J-field tests. It works in the canonical
three-port DtN boundary basis and builds a response-metric transport from the
symmetric DtN data on the same parent-fan face label.

For each double-history suffix class with root sectors 1,2,3 and each face
label F12,F23,F31, the loop is

    root 1 --identity face label--> root 2 --identity--> root 3 --identity--> root 1.

For each face label, it measures:

- raw same-label mismatch: amplitude, plane tilt, signed/signless J mismatch,
  K mismatch and metric mismatch;
- response-metric transport T_ij built from face DtN metrics G_i,G_j;
- noncommuting metric content [G_i,G_j];
- loop holonomy H = T_31 T_23 T_12;
- whether H transports the local J sector nontrivially around the loop.

Controls include identical-history clones, symmetrized birth, no-backreaction,
random same-level cycles, and diagonal/trace/common-mean reductions. Diagonal
or trace reductions should eliminate the noncommutative K sector; identical
history should have zero mismatch/holonomy.

No physical i, no global J, no *-algebra, no C*-norm, no GNS representation
and no AQFT net are introduced here.
"""

from __future__ import annotations

import argparse
import math
import random
import shutil
import zipfile
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

import numpy as np
from scipy.linalg import expm

from record_live_block_base import RealGrowth, LocalCell, build_cells, write_csv, offdiag_norm

EPS = 1e-12
LABELS = ("F12", "F23", "F31")
LOOP_STEP = 0.04


def fro(M: np.ndarray) -> float:
    return float(np.linalg.norm(M, ord="fro"))


def nrm(v: np.ndarray) -> float:
    return float(np.linalg.norm(np.asarray(v, dtype=float).reshape(-1)))


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


def upper(M: np.ndarray) -> np.ndarray:
    return np.asarray(M, dtype=float)[np.triu_indices(M.shape[0])]


def rel_dist(A: np.ndarray, B: np.ndarray, signless: bool = False) -> float:
    A = np.asarray(A, dtype=float)
    B = np.asarray(B, dtype=float)
    den = fro(A) + fro(B) + EPS
    if signless:
        return min(fro(A - B), fro(A + B)) / den
    return fro(A - B) / den


def axial_from_skew(K: np.ndarray) -> np.ndarray:
    return np.array([K[2, 1], K[0, 2], K[1, 0]], dtype=float)


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
        _, Q = np.linalg.eigh(M)
        def project(S: np.ndarray) -> np.ndarray:
            return Q @ np.diag(np.diag(Q.T @ S @ Q)) @ Q.T
        return project(Sa), project(Sb), project(Sc)
    raise ValueError(reduction)



def expm_step(A: np.ndarray, step: float) -> np.ndarray:
    A = np.asarray(A, dtype=float)
    scale = max(fro(A), EPS)
    # Keep the exponential diagnostic in the small-loop regime while preserving
    # relative comparisons across controls.
    M = (step / scale) * A if scale > 1.0 else step * A
    return np.asarray(expm(M), dtype=float)


def wilson_loop(edge_generators: List[np.ndarray], step: float) -> dict:
    # edge_generators are ordered 1->2, 2->3, 3->1.
    U12 = expm_step(edge_generators[0], step)
    U23 = expm_step(edge_generators[1], step)
    U31 = expm_step(edge_generators[2], step)
    W = U31 @ U23 @ U12
    I = np.eye(W.shape[0])
    vals = np.linalg.eigvals(W)
    return {
        "first": fro(edge_generators[0] + edge_generators[1] + edge_generators[2]),
        "edge": mean(fro(A) for A in edge_generators),
        "norm": fro(W - I) / (fro(I) + EPS),
        "skew_frac": fro(skew(W)) / (fro(W) + EPS),
        "sym_res": fro(sym(W) - I) / (fro(I) + EPS),
        "eig_imag_max": float(np.max(np.abs(np.imag(vals)))),
        "complex": int(float(np.max(np.abs(np.imag(vals)))) > 1e-8),
        "det": float(np.linalg.det(W)),
    }

def spd_project(G: np.ndarray, ridge: float) -> Tuple[np.ndarray, float, float]:
    G = sym(G)
    vals, Q = np.linalg.eigh(G)
    lam_max = max(float(np.max(np.abs(vals))), 1.0)
    floor = max(ridge, ridge * lam_max)
    vals_clip = np.maximum(vals, floor)
    out = Q @ np.diag(vals_clip) @ Q.T
    cond = float(np.max(vals_clip) / max(np.min(vals_clip), EPS))
    clipped = float(np.mean(vals < floor))
    return sym(out), cond, clipped


def invsqrt_and_sqrt(G: np.ndarray, ridge: float) -> Tuple[np.ndarray, np.ndarray, float, float]:
    Gp, cond, clipped = spd_project(G, ridge)
    vals, Q = np.linalg.eigh(Gp)
    vals = np.maximum(vals, ridge)
    S = Q @ np.diag(np.sqrt(vals)) @ Q.T
    IS = Q @ np.diag(1.0 / np.sqrt(vals)) @ Q.T
    return sym(S), sym(IS), cond, clipped


def metric_transport(G_src: np.ndarray, G_dst: np.ndarray, ridge: float) -> Tuple[np.ndarray, float, float]:
    S_src, _, cond_src, clip_src = invsqrt_and_sqrt(G_src, ridge)
    _, IS_dst, cond_dst, clip_dst = invsqrt_and_sqrt(G_dst, ridge)
    T = IS_dst @ S_src
    return T, max(cond_src, cond_dst), max(clip_src, clip_dst)


def matrix_safe_inv(T: np.ndarray) -> np.ndarray:
    return np.linalg.pinv(T, rcond=1e-10)


@dataclass(frozen=True)
class Face:
    label: str
    label_index: int
    ports: Tuple[int, int]
    parent: int
    level: int
    vertices: Tuple[int, int, int]


@dataclass
class JFace:
    face: Face
    Sa: np.ndarray
    Sb: np.ndarray
    Sc: np.ndarray
    G: np.ndarray
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
    metric_cond: float
    metric_clipped: float


@dataclass
class JNet:
    parent: int
    root_sector: int
    suffix: Tuple[int, ...]
    level: int
    address: Tuple[int, ...]
    faces: Dict[str, JFace]


def face_metric(Sa: np.ndarray, Sb: np.ndarray, Sc: np.ndarray, metric_mode: str, ridge: float) -> Tuple[np.ndarray, float, float]:
    if metric_mode == "mean_face":
        G0 = (Sa + Sb + Sc) / 3.0
    elif metric_mode == "parent_only":
        G0 = Sa
    elif metric_mode == "children_mean":
        G0 = 0.5 * (Sb + Sc)
    elif metric_mode == "edge_energy":
        Aab = Sb - Sa
        Abc = Sc - Sb
        Aca = Sa - Sc
        G0 = (Aab.T @ Aab + Abc.T @ Abc + Aca.T @ Aca) / 3.0
    else:
        raise ValueError(metric_mode)
    return spd_project(G0, ridge)


def local_face(vertices: Tuple[int, int, int], cells: Dict[int, LocalCell], source_mode: str, reduction: str, metric_mode: str, ridge: float, min_omega: float, face: Face) -> JFace:
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
    G, cond, clipped = face_metric(Sa, Sb, Sc, metric_mode, ridge)

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
        return JFace(face, Sa, Sb, Sc, G, K, None, None, axis, omega, K_norm, False, float("nan"), float("nan"), first, edge, cond, clipped)
    U2 = vecs[:, :2]
    P = U2 @ U2.T
    J = K / (omega + EPS)
    I = np.eye(3)
    J2_res = fro(P @ (J @ J + P) @ P) / (fro(P) + EPS)
    leakage = fro((I - P) @ J @ P) / (fro(J @ P) + EPS)
    return JFace(face, Sa, Sb, Sc, G, K, J, P, axis, omega, K_norm, True, J2_res, leakage, first, edge, cond, clipped)


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
        Face("F12", 0, (1, 2), parent, lvl, (parent, c1, c2)),
        Face("F23", 1, (2, 3), parent, lvl, (parent, c2, c3)),
        Face("F31", 2, (3, 1), parent, lvl, (parent, c3, c1)),
    ]


def build_jnets(model: RealGrowth, cells: Dict[int, LocalCell], source_mode: str, reduction: str, metric_mode: str, ridge: float, min_omega: float) -> Dict[int, JNet]:
    out: Dict[int, JNet] = {}
    for p in sorted(model.completed_parent_ids()):
        addr = address_tuple(model, p)
        if len(addr) < 1:
            continue
        faces = build_parent_faces(model, cells, p)
        if faces is None:
            continue
        jfaces = {f.label: local_face(f.vertices, cells, source_mode, reduction, metric_mode, ridge, min_omega, f) for f in faces}
        out[p] = JNet(p, int(addr[0]), tuple(addr[1:]), int(model.nodes[p].level), addr, jfaces)
    return out


def group_double_history(nets: Dict[int, JNet]) -> Dict[Tuple[int, Tuple[int, ...]], List[JNet]]:
    groups: Dict[Tuple[int, Tuple[int, ...]], List[JNet]] = {}
    for n in nets.values():
        groups.setdefault((n.level, n.suffix), []).append(n)
    return {k: sorted(v, key=lambda x: x.root_sector) for k, v in groups.items() if len({x.root_sector for x in v}) >= 2}


def three_sector_cycle(nets: List[JNet]) -> Optional[List[JNet]]:
    by_root = {n.root_sector: n for n in nets}
    if not all(r in by_root for r in (1, 2, 3)):
        return None
    return [by_root[1], by_root[2], by_root[3]]


def clone_with_root(net: JNet, root_sector: int) -> JNet:
    return JNet(net.parent, int(root_sector), net.suffix, net.level, (int(root_sector),) + tuple(net.suffix), net.faces)


def identical_cycles(nets: Dict[int, JNet], count: int) -> List[List[JNet]]:
    out: List[List[JNet]] = []
    for n in sorted(nets.values(), key=lambda x: (x.level, x.parent)):
        out.append([clone_with_root(n, 1), clone_with_root(n, 2), clone_with_root(n, 3)])
        if len(out) >= count:
            break
    return out


def random_cycles(nets: Dict[int, JNet], count: int, seed: int) -> List[List[JNet]]:
    rng = random.Random(seed)
    buckets: Dict[int, List[JNet]] = {}
    for n in nets.values():
        buckets.setdefault(n.level, []).append(n)
    levels = [lv for lv, xs in buckets.items() if len(xs) >= 3]
    out: List[List[JNet]] = []
    attempts = 0
    while levels and len(out) < count and attempts < max(3000, 50 * count):
        attempts += 1
        lv = rng.choice(levels)
        tri0 = rng.sample(buckets[lv], 3)
        if len({n.suffix for n in tri0}) == 1:
            continue
        tri0 = sorted(tri0, key=lambda n: (n.root_sector, n.parent))
        out.append([clone_with_root(tri0[i], i + 1) for i in range(3)])
    return out


def edge_stats(src: JFace, dst: JFace, ridge: float) -> dict:
    T, cond, clipped = metric_transport(src.G, dst.G, ridge)
    Ti = matrix_safe_inv(T)
    if src.J is not None and dst.J is not None:
        Jt = T @ src.J @ Ti
        J_signed = rel_dist(Jt, dst.J, signless=False)
        J_pm = rel_dist(Jt, dst.J, signless=True)
    else:
        J_signed = float("nan")
        J_pm = float("nan")
    Kt = T @ src.K @ Ti
    K_signed = rel_dist(Kt, dst.K, signless=False)
    K_pm = rel_dist(Kt, dst.K, signless=True)
    raw_J_signed = rel_dist(src.J, dst.J, False) if src.J is not None and dst.J is not None else float("nan")
    raw_J_pm = rel_dist(src.J, dst.J, True) if src.J is not None and dst.J is not None else float("nan")
    raw_K_signed = rel_dist(src.K, dst.K, False)
    raw_K_pm = rel_dist(src.K, dst.K, True)
    plane = rel_dist(src.P, dst.P, False) if src.P is not None and dst.P is not None else float("nan")
    amp = abs(src.omega - dst.omega) / (abs(src.omega) + abs(dst.omega) + EPS)
    metric = rel_dist(src.G, dst.G, False)
    metric_comm = fro(src.G @ dst.G - dst.G @ src.G) / ((fro(src.G) * fro(dst.G)) + EPS)
    return {
        "T": T,
        "T_norm": fro(T - np.eye(3)) / (fro(np.eye(3)) + EPS),
        "transport_cond": cond,
        "transport_clipped": clipped,
        "raw_J_signed": raw_J_signed,
        "raw_J_pm": raw_J_pm,
        "raw_K_signed": raw_K_signed,
        "raw_K_pm": raw_K_pm,
        "transport_J_signed": J_signed,
        "transport_J_pm": J_pm,
        "transport_K_signed": K_signed,
        "transport_K_pm": K_pm,
        "plane_tilt": plane,
        "amplitude_mismatch": amp,
        "metric_mismatch": metric,
        "metric_comm": metric_comm,
    }


def loop_face_summary(cycle: List[JNet], label: str, mode: str, control: str, source_mode: str, reduction: str, metric_mode: str, ridge: float, class_key: Tuple[int, Tuple[int, ...]]) -> Tuple[dict, List[dict]]:
    n1, n2, n3 = cycle
    fs = [n.faces[label] for n in cycle]
    edges = [(0, 1), (1, 2), (2, 0)]
    edge_rows: List[dict] = []
    edge_dicts: List[dict] = []
    for ei, (a, b) in enumerate(edges):
        d = edge_stats(fs[a], fs[b], ridge)
        edge_dicts.append(d)
        row = {k: v for k, v in d.items() if k != "T"}
        row.update(
            mode=mode,
            control=control,
            source_mode=source_mode,
            reduction=reduction,
            metric_mode=metric_mode,
            level=class_key[0],
            suffix=".".join(map(str, class_key[1])),
            face_label=label,
            edge_index=ei,
            src_root=cycle[a].root_sector,
            dst_root=cycle[b].root_sector,
            src_parent=cycle[a].parent,
            dst_parent=cycle[b].parent,
            src_addr=".".join(map(str, cycle[a].address)),
            dst_addr=".".join(map(str, cycle[b].address)),
        )
        edge_rows.append(row)

    H = edge_dicts[2]["T"] @ edge_dicts[1]["T"] @ edge_dicts[0]["T"]
    I = np.eye(3)
    hol = fro(H - I) / (fro(I) + EPS)
    hol_skew = fro(skew(H)) / (fro(H) + EPS)
    hol_sym = fro(sym(H) - I) / (fro(I) + EPS)
    detH = float(np.linalg.det(H))
    eigH = np.linalg.eigvals(H)
    eig_imag_max = float(np.max(np.abs(np.imag(eigH))))
    eig_complex = int(eig_imag_max > 1e-8)

    f0 = fs[0]
    Hi = matrix_safe_inv(H)
    if f0.J is not None:
        J_loop = H @ f0.J @ Hi
        loop_J_signed = rel_dist(J_loop, f0.J, False)
        loop_J_pm = rel_dist(J_loop, f0.J, True)
    else:
        loop_J_signed = float("nan")
        loop_J_pm = float("nan")
    K_loop = H @ f0.K @ Hi
    loop_K_signed = rel_dist(K_loop, f0.K, False)
    loop_K_pm = rel_dist(K_loop, f0.K, True)

    Gcomm_edges = [
        fs[0].G @ fs[1].G - fs[1].G @ fs[0].G,
        fs[1].G @ fs[2].G - fs[2].G @ fs[1].G,
        fs[2].G @ fs[0].G - fs[0].G @ fs[2].G,
    ]
    Gcomm_cycle = mean([
        fro(Gcomm_edges[0]) / ((fro(fs[0].G) * fro(fs[1].G)) + EPS),
        fro(Gcomm_edges[1]) / ((fro(fs[1].G) * fro(fs[2].G)) + EPS),
        fro(Gcomm_edges[2]) / ((fro(fs[2].G) * fro(fs[0].G)) + EPS),
    ])

    # Wilson-loop diagnostics from response mismatch itself.  These are not exact
    # metric-identifying maps, so they can detect second-order noncommutativity
    # that the pure metric-matching transport necessarily gauges away.
    metric_delta_edges = [fs[1].G - fs[0].G, fs[2].G - fs[1].G, fs[0].G - fs[2].G]
    K_delta_edges = [fs[1].K - fs[0].K, fs[2].K - fs[1].K, fs[0].K - fs[2].K]
    metric_comm_edges = [skew(Gcomm_edges[0]), skew(Gcomm_edges[1]), skew(Gcomm_edges[2])]
    mixed_edges = [skew(fs[0].G @ fs[1].K - fs[1].K @ fs[0].G),
                   skew(fs[1].G @ fs[2].K - fs[2].K @ fs[1].G),
                   skew(fs[2].G @ fs[0].K - fs[0].K @ fs[2].G)]
    W_metric = wilson_loop(metric_delta_edges, LOOP_STEP)
    W_K = wilson_loop(K_delta_edges, LOOP_STEP)
    W_Gcomm = wilson_loop(metric_comm_edges, LOOP_STEP)
    W_mixed = wilson_loop(mixed_edges, LOOP_STEP)

    row = {
        "mode": mode,
        "control": control,
        "source_mode": source_mode,
        "reduction": reduction,
        "metric_mode": metric_mode,
        "level": class_key[0],
        "suffix": ".".join(map(str, class_key[1])),
        "parents": " ".join(str(n.parent) for n in cycle),
        "addresses": " | ".join(".".join(map(str, n.address)) for n in cycle),
        "face_label": label,
        "valid_frac": mean(int(f.valid) for f in fs),
        "mean_K": mean(f.K_norm for f in fs),
        "mean_omega": mean(f.omega for f in fs),
        "mean_J2": mean(f.J2_res for f in fs),
        "mean_first": mean(f.first_norm for f in fs),
        "mean_edge": mean(f.edge_norm for f in fs),
        "mean_metric_cond": mean(f.metric_cond for f in fs),
        "mean_metric_clipped": mean(f.metric_clipped for f in fs),
        "raw_J_signed": mean(d["raw_J_signed"] for d in edge_dicts),
        "raw_J_pm": mean(d["raw_J_pm"] for d in edge_dicts),
        "raw_K_signed": mean(d["raw_K_signed"] for d in edge_dicts),
        "raw_K_pm": mean(d["raw_K_pm"] for d in edge_dicts),
        "transport_J_signed": mean(d["transport_J_signed"] for d in edge_dicts),
        "transport_J_pm": mean(d["transport_J_pm"] for d in edge_dicts),
        "transport_K_signed": mean(d["transport_K_signed"] for d in edge_dicts),
        "transport_K_pm": mean(d["transport_K_pm"] for d in edge_dicts),
        "plane_tilt": mean(d["plane_tilt"] for d in edge_dicts),
        "amplitude_mismatch": mean(d["amplitude_mismatch"] for d in edge_dicts),
        "metric_mismatch": mean(d["metric_mismatch"] for d in edge_dicts),
        "metric_comm": mean(d["metric_comm"] for d in edge_dicts),
        "metric_comm_cycle": Gcomm_cycle,
        "T_norm": mean(d["T_norm"] for d in edge_dicts),
        "transport_cond": mean(d["transport_cond"] for d in edge_dicts),
        "transport_clipped": mean(d["transport_clipped"] for d in edge_dicts),
        "wil_metric_first": W_metric["first"],
        "wil_metric_edge": W_metric["edge"],
        "wil_metric_norm": W_metric["norm"],
        "wil_metric_skew": W_metric["skew_frac"],
        "wil_metric_complex": W_metric["complex"],
        "wil_K_first": W_K["first"],
        "wil_K_edge": W_K["edge"],
        "wil_K_norm": W_K["norm"],
        "wil_K_skew": W_K["skew_frac"],
        "wil_K_complex": W_K["complex"],
        "wil_Gcomm_edge": W_Gcomm["edge"],
        "wil_Gcomm_norm": W_Gcomm["norm"],
        "wil_Gcomm_skew": W_Gcomm["skew_frac"],
        "wil_Gcomm_complex": W_Gcomm["complex"],
        "wil_mixed_edge": W_mixed["edge"],
        "wil_mixed_norm": W_mixed["norm"],
        "wil_mixed_skew": W_mixed["skew_frac"],
        "wil_mixed_complex": W_mixed["complex"],
        "holonomy_norm": hol,
        "holonomy_skew_frac": hol_skew,
        "holonomy_sym_res": hol_sym,
        "holonomy_det": detH,
        "holonomy_eig_imag_max": eig_imag_max,
        "holonomy_complex": eig_complex,
        "loop_J_signed": loop_J_signed,
        "loop_J_pm": loop_J_pm,
        "loop_K_signed": loop_K_signed,
        "loop_K_pm": loop_K_pm,
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
            valid=mean(r["valid_frac"] for r in rs),
            K=mean(r["mean_K"] for r in rs),
            omega=mean(r["mean_omega"] for r in rs),
            J2=mean(r["mean_J2"] for r in rs),
            raw_J_pm=mean(r["raw_J_pm"] for r in rs),
            raw_J_signed=mean(r["raw_J_signed"] for r in rs),
            plane_tilt=mean(r["plane_tilt"] for r in rs),
            amp=mean(r["amplitude_mismatch"] for r in rs),
            metric_mismatch=mean(r["metric_mismatch"] for r in rs),
            metric_comm=mean(r["metric_comm"] for r in rs),
            metric_comm_cycle=mean(r["metric_comm_cycle"] for r in rs),
            T_norm=mean(r["T_norm"] for r in rs),
            wil_metric=mean(r["wil_metric_norm"] for r in rs),
            wil_metric_skew=mean(r["wil_metric_skew"] for r in rs),
            wil_metric_complex=mean(r["wil_metric_complex"] for r in rs),
            wil_K=mean(r["wil_K_norm"] for r in rs),
            wil_K_skew=mean(r["wil_K_skew"] for r in rs),
            wil_K_complex=mean(r["wil_K_complex"] for r in rs),
            wil_Gcomm=mean(r["wil_Gcomm_norm"] for r in rs),
            wil_Gcomm_skew=mean(r["wil_Gcomm_skew"] for r in rs),
            wil_Gcomm_complex=mean(r["wil_Gcomm_complex"] for r in rs),
            wil_mixed=mean(r["wil_mixed_norm"] for r in rs),
            wil_mixed_skew=mean(r["wil_mixed_skew"] for r in rs),
            wil_mixed_complex=mean(r["wil_mixed_complex"] for r in rs),
            hol=mean(r["holonomy_norm"] for r in rs),
            hol_p95=perc([r["holonomy_norm"] for r in rs], 95),
            hol_skew=mean(r["holonomy_skew_frac"] for r in rs),
            hol_complex=mean(r["holonomy_complex"] for r in rs),
            loop_J_pm=mean(r["loop_J_pm"] for r in rs),
            loop_J_signed=mean(r["loop_J_signed"] for r in rs),
            loop_K_pm=mean(r["loop_K_pm"] for r in rs),
            transport_J_pm=mean(r["transport_J_pm"] for r in rs),
            transport_J_signed=mean(r["transport_J_signed"] for r in rs),
            transport_K_pm=mean(r["transport_K_pm"] for r in rs),
            cond=mean(r["mean_metric_cond"] for r in rs),
            clipped=mean(r["mean_metric_clipped"] for r in rs),
        )
        out.append(d)
    return out


def run(args: argparse.Namespace) -> str:
    outdir: Path = args.outdir
    if outdir.exists():
        shutil.rmtree(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    controls = ["real_growth", "symmetrized_birth", "no_backreaction"]
    source_modes = ["live"] if args.quick_suite else ["live", "handoff", "aging"]
    reductions = ["full", "diagonal", "trace_scalar"] if args.quick_suite else ["full", "diagonal", "trace_scalar", "common_mean_diagonal"]
    metric_modes = args.metric_modes.split(",") if args.metric_modes else (["mean_face", "parent_only", "children_mean", "edge_energy"] if not args.quick_suite else ["mean_face", "edge_energy"])

    loop_rows: List[dict] = []
    edge_rows: List[dict] = []
    model_rows: List[dict] = []

    for control in controls:
        model = build_model(control, args.max_level, args.mode)
        cells = build_cells(model, args.max_level)
        for source_mode in source_modes:
            for reduction in reductions:
                for metric_mode in metric_modes:
                    nets = build_jnets(model, cells, source_mode, reduction, metric_mode, args.ridge, args.min_omega)
                    groups = group_double_history(nets)
                    cycles: List[Tuple[Tuple[int, Tuple[int, ...]], List[JNet]]] = []
                    for key, gs in sorted(groups.items(), key=lambda kv: (kv[0][0], kv[0][1])):
                        cyc = three_sector_cycle(gs)
                        if cyc is not None:
                            cycles.append((key, cyc))
                    if source_mode == "live" and reduction == "full" and metric_mode == metric_modes[0]:
                        model_rows.append({
                            "control": control,
                            "nodes": len(model.nodes),
                            "cells": len(cells),
                            "jnets": len(nets),
                            "double_history_classes": len(groups),
                            "three_sector_cycles": len(cycles),
                        })
                    ident = identical_cycles(nets, min(args.identical_cycles, len(nets)))
                    rand = random_cycles(nets, args.random_cycles, args.seed + 131) if args.random_cycles else []
                    for key, cyc in cycles:
                        for lab in LABELS:
                            r, ers = loop_face_summary(cyc, lab, "double_history_suffix_cycle", control, source_mode, reduction, metric_mode, args.ridge, key)
                            loop_rows.append(r)
                            edge_rows.extend(ers)
                    for cyc in ident:
                        key = (cyc[0].level, tuple([-111, cyc[0].parent]))
                        for lab in LABELS:
                            r, ers = loop_face_summary(cyc, lab, "identical_history_clone_control", control, source_mode, reduction, metric_mode, args.ridge, key)
                            loop_rows.append(r)
                            edge_rows.extend(ers)
                    for i, cyc in enumerate(rand):
                        key = (cyc[0].level, tuple([-999, i]))
                        for lab in LABELS:
                            r, ers = loop_face_summary(cyc, lab, "random_same_level_cycle_baseline", control, source_mode, reduction, metric_mode, args.ridge, key)
                            loop_rows.append(r)
                            edge_rows.extend(ers)

    main = summarize(loop_rows, ["mode", "control", "source_mode", "reduction", "metric_mode"])
    by_label = summarize(loop_rows, ["mode", "control", "source_mode", "reduction", "metric_mode", "face_label"])
    by_level = summarize(loop_rows, ["mode", "control", "level", "source_mode", "reduction", "metric_mode"])

    write_csv(outdir / "response_connection_loop_rows.csv", loop_rows)
    write_csv(outdir / "response_connection_edge_rows.csv", edge_rows)
    write_csv(outdir / "response_connection_summary_main.csv", main)
    write_csv(outdir / "response_connection_summary_by_label.csv", by_label)
    write_csv(outdir / "response_connection_summary_by_level.csv", by_level)
    write_csv(outdir / "model_summaries.csv", model_rows)

    def findrow(mode: str, control: str, source: str, red: str, metric: str) -> Optional[dict]:
        for d in main:
            if d["mode"] == mode and d["control"] == control and d["source_mode"] == source and d["reduction"] == red and d["metric_mode"] == metric:
                return d
        return None

    selected = [
        ("double_history_suffix_cycle", "real_growth", "live", "full", "mean_face"),
        ("double_history_suffix_cycle", "real_growth", "live", "full", "edge_energy"),
        ("double_history_suffix_cycle", "symmetrized_birth", "live", "full", "mean_face"),
        ("double_history_suffix_cycle", "no_backreaction", "live", "full", "mean_face"),
        ("identical_history_clone_control", "real_growth", "live", "full", "mean_face"),
        ("random_same_level_cycle_baseline", "real_growth", "live", "full", "mean_face"),
        ("double_history_suffix_cycle", "real_growth", "live", "diagonal", "mean_face"),
        ("double_history_suffix_cycle", "real_growth", "live", "trace_scalar", "mean_face"),
        ("double_history_suffix_cycle", "real_growth", "handoff", "full", "mean_face"),
        ("double_history_suffix_cycle", "real_growth", "aging", "full", "mean_face"),
    ]

    def fmt(d: Optional[dict]) -> str:
        if d is None:
            return "missing"
        return (
            f"{d['mode']} | {d['control']} | {d['source_mode']} | {d['reduction']} | {d['metric_mode']}: "
            f"count={d['count']}, valid={d['valid']:.3f}, K={d['K']:.6g}, "
            f"rawJ±={d['raw_J_pm']:.4g}, plane={d['plane_tilt']:.4g}, amp={d['amp']:.4g}, "
            f"metric={d['metric_mismatch']:.4g}, Gcomm={d['metric_comm']:.4g}, "
            f"T={d['T_norm']:.4g}, metricWilson={d['wil_metric']:.4g}, KWilson={d['wil_K']:.4g}, "
            f"GcommWilson={d['wil_Gcomm']:.4g}, mixedWilson={d['wil_mixed']:.4g}, "
            f"hol={d['hol']:.4g}, hol_skew={d['hol_skew']:.4g}, hol_complex={d['hol_complex']:.3f}, "
            f"loopJ±={d['loop_J_pm']:.4g}, transJ±={d['transport_J_pm']:.4g}, cond={d['cond']:.4g}"
        )

    lines = [
        "CNNA NONFLAT CONNECTION FROM SAME-LABEL RESPONSE MISMATCH",
        f"max_level={args.max_level}, mode={args.mode}, min_omega={args.min_omega}, ridge={args.ridge}",
        "",
        "MODEL SUMMARIES",
    ]
    for m in model_rows:
        lines.append(
            f"  {m['control']}: nodes={m['nodes']}, cells={m['cells']}, jnets={m['jnets']}, "
            f"double_classes={m['double_history_classes']}, three_sector_cycles={m['three_sector_cycles']}"
        )
    lines.append("")
    lines.append("SELECTED SUMMARIES")
    for key in selected:
        lines.append("  " + fmt(findrow(*key)))
    lines.append("")
    lines.append("READING RULE")
    lines.extend([
        "  Identity/same-label face gluing is taken as the boundary-resolved result",
        "  from the previous test.  This script does not search over face labels and",
        "  does not impose a cyclic Z3 shift.",
        "  The response connection uses only symmetric DtN face metrics G and response",
        "  mismatch generators on same-label faces.  The exact metric transport",
        "  T_ij = G_j^{-1/2} G_i^{1/2} is included only as a gauge audit, because it",
        "  can telescope by construction.  The Wilson columns metricWilson, KWilson,",
        "  GcommWilson and mixedWilson test the response-mismatch generators themselves.",
        "  A positive nonflat-connection gate would need real_growth/full/live to have",
        "  Wilson curvature, transported-J loop residual or holonomy clearly above identical-history,",
        "  diagonal/trace and random controls, and preferably stronger than",
        "  symmetrized_birth.  Otherwise the local parent-fan J sector remains a local",
        "  Stufe-2/3 phenomenon without Stufe-4 locking.",
    ])
    summary = "\n".join(lines)
    (outdir / "SUMMARY.txt").write_text(summary, encoding="utf-8")

    results = f"""# Results: nonflat connection from same-label response mismatch

This test accepts identity/same-label face gluing and asks whether the remaining
DtN/response mismatch defines a nonflat connection around double-history cycles.

```text
{summary}
```

## Output files

- `response_connection_loop_rows.csv`
- `response_connection_edge_rows.csv`
- `response_connection_summary_main.csv`
- `response_connection_summary_by_label.csv`
- `response_connection_summary_by_level.csv`
- `model_summaries.csv`
- `SUMMARY.txt`

## Status

The primary row is `double_history_suffix_cycle | real_growth | live | full |
mean_face`.  The strongest anti-smuggling controls are identical-history,
diagonal/trace reduction and random same-level cycles.
"""
    (outdir / "RESULTS_nonflat_connection_from_response_mismatch.md").write_text(results, encoding="utf-8")

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
    ap.add_argument("--ridge", type=float, default=1e-9)
    ap.add_argument("--random-cycles", type=int, default=121)
    ap.add_argument("--identical-cycles", type=int, default=121)
    ap.add_argument("--seed", type=int, default=20260621)
    ap.add_argument("--metric-modes", default="")
    ap.add_argument("--quick-suite", action="store_true")
    ap.add_argument("--outdir", type=Path, default=Path("nonflat_response_connection_out_L6"))
    args = ap.parse_args()
    print(run(args))


if __name__ == "__main__":
    main()
