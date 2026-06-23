#!/usr/bin/env python3
"""
CNNA Growth Test: boundary-resolved face-correspondence transport.

Purpose
-------
The preceding Z3-label test showed that a prescribed directed cyclic label
shift is not inferred by the parent-fan J/Plaquette data.  This diagnostic now
removes face labels from the identification step.

For each double-history suffix class with root sectors 1,2,3, every parent-fan
net has three faces.  Instead of comparing F12 to F12 or imposing a cyclic
shift, the script builds boundary-resolved fingerprints for each face and
infers the best face correspondence between nets by minimizing the total cost
over all six S3 permutations.

The anti-smuggling rule is central: the primary fingerprint `boundary_ports`
uses DtN/Record-Live boundary and port data, but not J, K or complex eigenvalue
information.  K/J-based fingerprints are included only as audits.

The gate is whether the boundary data themselves force a nontrivial face
correspondence across double-history gluing.  If the best permutation remains
identity/same-label, then the current Double-History quotient carries the local
parent-fan J sector but does not transport it nontrivially.

No physical i, no global J-field, no C*-norm, no GNS representation and no AQFT
net are introduced here.
"""

from __future__ import annotations

import argparse
import itertools
import math
import random
import shutil
import zipfile
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

import numpy as np

from record_live_block_base import RealGrowth, LocalCell, build_cells, write_csv, fro as base_fro, offdiag_norm

EPS = 1e-12
LABELS = ("F12", "F23", "F31")
PERMS = list(itertools.permutations(range(3)))
PERM_NAMES = {p: "".join(str(i + 1) for i in p) for p in PERMS}
IDENTITY = (0, 1, 2)
CYCLIC_PLUS = (1, 2, 0)
CYCLIC_MINUS = (2, 0, 1)
REVERSALS = {p for p in PERMS if p not in {IDENTITY, CYCLIC_PLUS, CYCLIC_MINUS}}


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


def upper(M: np.ndarray) -> np.ndarray:
    idx = np.triu_indices(M.shape[0])
    return np.asarray(M, dtype=float)[idx]


def row_unit(v: np.ndarray) -> np.ndarray:
    v = np.asarray(v, dtype=float)
    s = nrm(v)
    if s < EPS:
        return np.zeros_like(v)
    return v / s


def standardize_vector(v: np.ndarray, scale_mode: str) -> np.ndarray:
    v = np.asarray(v, dtype=float).reshape(-1)
    if scale_mode == "raw":
        return v
    if scale_mode == "unit":
        s = nrm(v)
        return v / (s + EPS)
    if scale_mode == "zscore":
        mu = float(np.mean(v)) if v.size else 0.0
        sd = float(np.std(v)) if v.size else 0.0
        return (v - mu) / (sd + EPS)
    raise ValueError(scale_mode)


def vector_distance(a: np.ndarray, b: np.ndarray, signless: bool = False) -> float:
    a = np.asarray(a, dtype=float).reshape(-1)
    b = np.asarray(b, dtype=float).reshape(-1)
    if a.size != b.size:
        raise ValueError("fingerprints have different sizes")
    if nrm(a) < EPS and nrm(b) < EPS:
        return 0.0
    den = nrm(a) + nrm(b) + EPS
    if signless:
        return min(nrm(a - b), nrm(a + b)) / den
    return nrm(a - b) / den


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
    label_index: int
    ports: Tuple[int, int]
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
        Face("F12", 0, (1, 2), parent, lvl, (parent, c1, c2)),
        Face("F23", 1, (2, 3), parent, lvl, (parent, c2, c3)),
        Face("F31", 2, (3, 1), parent, lvl, (parent, c3, c1)),
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


def face_matrices(jf: JFace, cells: Dict[int, LocalCell], source_mode: str, reduction: str) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    a, b, c = jf.face.vertices
    Sa0 = raw_vertex_operator(cells[a], source_mode)
    Sb0 = raw_vertex_operator(cells[b], source_mode)
    Sc0 = raw_vertex_operator(cells[c], source_mode)
    return reduce_ops(Sa0, Sb0, Sc0, reduction)


def cell_profile(cell: LocalCell, source_mode: str) -> np.ndarray:
    S = sym(cell.matrix_for_source(source_mode))
    vals = np.linalg.eigvalsh(S)
    parts = [
        upper(S),
        np.diag(S),
        vals,
        np.array([float(np.trace(S)), fro(S), offdiag_norm(S), float(cell.age), float(cell.level)], dtype=float),
    ]
    # Add DtN decomposed boundary rows when available. These are scalar real
    # boundary profiles, not J/K data.
    if source_mode == "live":
        d = cell.live1
    elif source_mode == "record":
        d = cell.record
    else:
        d = cell.live1
    for key in ("env_vec", "uv_vec", "port_shunt"):
        if key in d:
            parts.append(np.asarray(d[key], dtype=float).reshape(-1))
    return np.concatenate(parts)


def port_boundary_profile(cell: LocalCell, ports: Tuple[int, int], source_mode: str) -> np.ndarray:
    S = sym(cell.matrix_for_source(source_mode))
    i, j = ports
    ri = S[i - 1, :]
    rj = S[j - 1, :]
    return np.concatenate([
        ri,
        rj,
        ri - rj,
        np.array([S[i - 1, i - 1], S[j - 1, j - 1], S[i - 1, j - 1], fro(S), float(np.trace(S))], dtype=float),
    ])


def face_fingerprint(jf: JFace, cells: Dict[int, LocalCell], source_mode: str, reduction: str, fingerprint: str, scale_mode: str) -> np.ndarray:
    face = jf.face
    parent, child_i, child_j = face.vertices
    ci = cells[parent]
    cb = cells[child_i]
    cc = cells[child_j]
    Sa, Sb, Sc = face_matrices(jf, cells, source_mode, reduction)
    Aab = Sb - Sa
    Abc = Sc - Sb
    Aca = Sa - Sc
    i, j = face.ports

    if fingerprint == "boundary_ports":
        # Primary anti-smuggling fingerprint: uses boundary/port DtN profiles and
        # vertex response summaries, but not K/J/axis data.
        v = np.concatenate([
            port_boundary_profile(ci, (i, j), source_mode),
            cell_profile(cb, source_mode),
            cell_profile(cc, source_mode),
            upper(Sa), upper(Sb), upper(Sc),
            upper(Aab), upper(Abc), upper(Aca),
            np.array([fro(Aab), fro(Abc), fro(Aca), jf.edge_norm, jf.first_norm], dtype=float),
        ])
        return standardize_vector(v, scale_mode)

    if fingerprint == "boundary_ports_no_edges":
        v = np.concatenate([
            port_boundary_profile(ci, (i, j), source_mode),
            cell_profile(cb, source_mode),
            cell_profile(cc, source_mode),
            upper(Sa), upper(Sb), upper(Sc),
        ])
        return standardize_vector(v, scale_mode)

    if fingerprint == "record_live_boundary":
        # Uses record/live/handoff profiles simultaneously; still no K/J.
        parts: List[np.ndarray] = []
        for sm in ("record", "live", "handoff", "aging"):
            parts.append(port_boundary_profile(ci, (i, j), sm))
            parts.append(cell_profile(cb, sm if sm in {"record", "live"} else source_mode))
            parts.append(cell_profile(cc, sm if sm in {"record", "live"} else source_mode))
        return standardize_vector(np.concatenate(parts), scale_mode)

    if fingerprint == "K_signed":
        return standardize_vector(jf.K.reshape(-1), scale_mode)

    if fingerprint == "K_signless":
        return standardize_vector(jf.K.reshape(-1), scale_mode)

    if fingerprint == "axis_signed":
        return standardize_vector(jf.axis, scale_mode)

    if fingerprint == "axis_abs":
        return standardize_vector(np.abs(jf.axis), scale_mode)

    if fingerprint == "amplitude":
        v = np.array([jf.K_norm, jf.omega, jf.edge_norm, jf.first_norm, jf.J2_res if np.isfinite(jf.J2_res) else 0.0], dtype=float)
        return standardize_vector(v, scale_mode)

    if fingerprint == "face_dtn_full":
        v = np.concatenate([upper(Sa), upper(Sb), upper(Sc), upper(Aab), upper(Abc), upper(Aca), cell_profile(ci, source_mode)])
        return standardize_vector(v, scale_mode)

    raise ValueError(fingerprint)


def fingerprint_signless(fingerprint: str) -> bool:
    return fingerprint in {"K_signless", "axis_abs"}


def perm_cost(src: JNet, dst: JNet, cells: Dict[int, LocalCell], source_mode: str, reduction: str, fingerprint: str, scale_mode: str, perm: Tuple[int, int, int]) -> Tuple[float, float]:
    ds: List[float] = []
    valid = 0
    src_faces = [src.faces[lab] for lab in LABELS]
    dst_faces = [dst.faces[lab] for lab in LABELS]
    for si, di in enumerate(perm):
        a = src_faces[si]
        b = dst_faces[di]
        try:
            va = face_fingerprint(a, cells, source_mode, reduction, fingerprint, scale_mode)
            vb = face_fingerprint(b, cells, source_mode, reduction, fingerprint, scale_mode)
            d = vector_distance(va, vb, signless=fingerprint_signless(fingerprint))
        except Exception:
            d = float("nan")
        if np.isfinite(d):
            ds.append(d)
            valid += 1
    return mean(ds), valid / 3.0


def directed_perm(src_root: int, dst_root: int) -> Tuple[int, int, int]:
    # root sectors are 1,2,3. The directed label shift q=(dst-src) mod 3 maps
    # source face index i -> i+q. As a permutation tuple, perm[si]=di.
    q = ((int(dst_root) - 1) - (int(src_root) - 1)) % 3
    return tuple((i + q) % 3 for i in range(3))


def inverse_perm(p: Tuple[int, int, int]) -> Tuple[int, int, int]:
    inv = [0, 0, 0]
    for i, j in enumerate(p):
        inv[j] = i
    return tuple(inv)


def compose_perm(p: Tuple[int, int, int], q: Tuple[int, int, int]) -> Tuple[int, int, int]:
    # Apply p, then q.
    return tuple(q[p[i]] for i in range(3))


def perm_category(p: Tuple[int, int, int]) -> str:
    if p == IDENTITY:
        return "identity"
    if p == CYCLIC_PLUS:
        return "cyclic_plus"
    if p == CYCLIC_MINUS:
        return "cyclic_minus"
    return "reflection"


def infer_perm(src: JNet, dst: JNet, cells: Dict[int, LocalCell], source_mode: str, reduction: str, fingerprint: str, scale_mode: str) -> dict:
    costs: Dict[Tuple[int, int, int], float] = {}
    valids: Dict[Tuple[int, int, int], float] = {}
    for p in PERMS:
        c, vf = perm_cost(src, dst, cells, source_mode, reduction, fingerprint, scale_mode, p)
        costs[p] = c
        valids[p] = vf
    finite = [(p, c) for p, c in costs.items() if np.isfinite(c)]
    if not finite:
        best = IDENTITY
        best_cost = float("nan")
        second_cost = float("nan")
    else:
        finite.sort(key=lambda x: x[1])
        best = finite[0][0]
        best_cost = finite[0][1]
        second_cost = finite[1][1] if len(finite) > 1 else float("nan")
    dperm = directed_perm(src.root_sector, dst.root_sector)
    return {
        "src_parent": src.parent,
        "dst_parent": dst.parent,
        "src_root": src.root_sector,
        "dst_root": dst.root_sector,
        "best_perm": PERM_NAMES[best],
        "best_category": perm_category(best),
        "best_cost": best_cost,
        "second_cost": second_cost,
        "gap": second_cost - best_cost if np.isfinite(second_cost) and np.isfinite(best_cost) else float("nan"),
        "identity_cost": costs[IDENTITY],
        "cyclic_plus_cost": costs[CYCLIC_PLUS],
        "cyclic_minus_cost": costs[CYCLIC_MINUS],
        "directed_perm": PERM_NAMES[dperm],
        "directed_cost": costs[dperm],
        "valid_frac": valids[best],
        "identity_match": int(best == IDENTITY),
        "directed_match": int(best == dperm),
        "cyclic_match": int(best in {CYCLIC_PLUS, CYCLIC_MINUS}),
        "reflection_match": int(best in REVERSALS),
        "adv_directed_vs_identity": costs[IDENTITY] - costs[dperm] if np.isfinite(costs[IDENTITY]) and np.isfinite(costs[dperm]) else float("nan"),
    }


def cycle_summary(cycle: List[JNet], cells: Dict[int, LocalCell], mode: str, control: str, source_mode: str, reduction: str, fingerprint: str, scale_mode: str, class_key: Tuple[int, Tuple[int, ...]]) -> Tuple[dict, List[dict]]:
    edges = [(cycle[0], cycle[1]), (cycle[1], cycle[2]), (cycle[2], cycle[0])]
    edge_rows: List[dict] = []
    inferred: List[Tuple[int, int, int]] = []
    directed: List[Tuple[int, int, int]] = []
    for ei, (src, dst) in enumerate(edges):
        er = infer_perm(src, dst, cells, source_mode, reduction, fingerprint, scale_mode)
        p = tuple(int(ch) - 1 for ch in er["best_perm"])
        d = tuple(int(ch) - 1 for ch in er["directed_perm"])
        inferred.append(p)
        directed.append(d)
        er.update(
            edge_index=ei,
            mode=mode,
            control=control,
            source_mode=source_mode,
            reduction=reduction,
            fingerprint=fingerprint,
            scale_mode=scale_mode,
            level=class_key[0],
            suffix=".".join(map(str, class_key[1])),
            src_addr=".".join(map(str, src.address)),
            dst_addr=".".join(map(str, dst.address)),
        )
        edge_rows.append(er)
    product = compose_perm(compose_perm(inferred[0], inferred[1]), inferred[2]) if len(inferred) == 3 else IDENTITY
    directed_product = compose_perm(compose_perm(directed[0], directed[1]), directed[2]) if len(directed) == 3 else IDENTITY
    faces = [n.faces[lab] for n in cycle for lab in LABELS]
    row = {
        "mode": mode,
        "control": control,
        "source_mode": source_mode,
        "reduction": reduction,
        "fingerprint": fingerprint,
        "scale_mode": scale_mode,
        "level": class_key[0],
        "suffix": ".".join(map(str, class_key[1])),
        "parents": " ".join(str(n.parent) for n in cycle),
        "addresses": " | ".join(".".join(map(str, n.address)) for n in cycle),
        "cycle_valid_face_frac": mean(int(f.valid) for f in faces),
        "mean_K": mean(f.K_norm for f in faces),
        "mean_J2": mean(f.J2_res for f in faces),
        "mean_edge_valid_frac": mean(er["valid_frac"] for er in edge_rows),
        "mean_best_cost": mean(er["best_cost"] for er in edge_rows),
        "mean_identity_cost": mean(er["identity_cost"] for er in edge_rows),
        "mean_directed_cost": mean(er["directed_cost"] for er in edge_rows),
        "mean_gap": mean(er["gap"] for er in edge_rows),
        "identity_match_frac": mean(er["identity_match"] for er in edge_rows),
        "directed_match_frac": mean(er["directed_match"] for er in edge_rows),
        "cyclic_match_frac": mean(er["cyclic_match"] for er in edge_rows),
        "reflection_match_frac": mean(er["reflection_match"] for er in edge_rows),
        "adv_directed_vs_identity": mean(er["adv_directed_vs_identity"] for er in edge_rows),
        "inferred_cycle_product": PERM_NAMES[product],
        "inferred_cycle_category": perm_category(product),
        "directed_cycle_product": PERM_NAMES[directed_product],
        "cycle_identity_product": int(product == IDENTITY),
        "cycle_nontrivial_product": int(product != IDENTITY),
        "cycle_matches_directed_product": int(product == directed_product),
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
            valid=mean(r["cycle_valid_face_frac"] for r in rs),
            K=mean(r["mean_K"] for r in rs),
            J2=mean(r["mean_J2"] for r in rs),
            edge_valid=mean(r["mean_edge_valid_frac"] for r in rs),
            best_cost=mean(r["mean_best_cost"] for r in rs),
            identity_cost=mean(r["mean_identity_cost"] for r in rs),
            directed_cost=mean(r["mean_directed_cost"] for r in rs),
            gap=mean(r["mean_gap"] for r in rs),
            identity_match=mean(r["identity_match_frac"] for r in rs),
            directed_match=mean(r["directed_match_frac"] for r in rs),
            cyclic_match=mean(r["cyclic_match_frac"] for r in rs),
            reflection_match=mean(r["reflection_match_frac"] for r in rs),
            adv_directed_vs_identity=mean(r["adv_directed_vs_identity"] for r in rs),
            cycle_identity=mean(r["cycle_identity_product"] for r in rs),
            cycle_nontrivial=mean(r["cycle_nontrivial_product"] for r in rs),
            cycle_matches_directed=mean(r["cycle_matches_directed_product"] for r in rs),
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
    fingerprints = (["boundary_ports", "boundary_ports_no_edges", "record_live_boundary", "face_dtn_full", "K_signed", "K_signless", "amplitude"] if args.quick_suite else [
        "boundary_ports",
        "boundary_ports_no_edges",
        "record_live_boundary",
        "face_dtn_full",
        "K_signed",
        "K_signless",
        "axis_abs",
        "amplitude",
    ])
    if args.fingerprints:
        fingerprints = args.fingerprints.split(",")

    cycle_rows: List[dict] = []
    edge_rows: List[dict] = []
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
                ident = identical_cycles(nets, min(args.identical_cycles, len(nets)))
                rand = random_cycles(nets, args.random_cycles, args.seed + 77) if args.random_cycles else []
                for fingerprint in fingerprints:
                    for key, cyc in dh_cycles:
                        row, ers = cycle_summary(cyc, cells, "double_history_suffix_cycle", control, source_mode, reduction, fingerprint, args.scale_mode, key)
                        cycle_rows.append(row)
                        edge_rows.extend(ers)
                    for cyc in ident:
                        key = (cyc[0].level, tuple([-111, cyc[0].parent]))
                        row, ers = cycle_summary(cyc, cells, "identical_history_clone_control", control, source_mode, reduction, fingerprint, args.scale_mode, key)
                        cycle_rows.append(row)
                        edge_rows.extend(ers)
                    for i, cyc in enumerate(rand):
                        key = (cyc[0].level, tuple([-999, i]))
                        row, ers = cycle_summary(cyc, cells, "random_same_level_cycle_baseline", control, source_mode, reduction, fingerprint, args.scale_mode, key)
                        cycle_rows.append(row)
                        edge_rows.extend(ers)

    main = summarize(cycle_rows, ["mode", "control", "source_mode", "reduction", "fingerprint", "scale_mode"])
    by_level = summarize(cycle_rows, ["mode", "control", "level", "source_mode", "reduction", "fingerprint", "scale_mode"])

    write_csv(outdir / "boundary_face_cycle_rows.csv", cycle_rows)
    write_csv(outdir / "boundary_face_edge_rows.csv", edge_rows)
    write_csv(outdir / "boundary_face_summary_main.csv", main)
    write_csv(outdir / "boundary_face_summary_by_level.csv", by_level)
    write_csv(outdir / "model_summaries.csv", model_rows)

    def findrow(mode: str, control: str, source: str, red: str, fp: str) -> Optional[dict]:
        for d in main:
            if d["mode"] == mode and d["control"] == control and d["source_mode"] == source and d["reduction"] == red and d["fingerprint"] == fp:
                return d
        return None

    selected = [
        ("double_history_suffix_cycle", "real_growth", "live", "full", "boundary_ports"),
        ("double_history_suffix_cycle", "real_growth", "live", "full", "boundary_ports_no_edges"),
        ("double_history_suffix_cycle", "real_growth", "live", "full", "record_live_boundary"),
        ("double_history_suffix_cycle", "real_growth", "live", "full", "face_dtn_full"),
        ("double_history_suffix_cycle", "real_growth", "live", "full", "K_signed"),
        ("double_history_suffix_cycle", "real_growth", "live", "full", "K_signless"),
        ("double_history_suffix_cycle", "symmetrized_birth", "live", "full", "boundary_ports"),
        ("double_history_suffix_cycle", "no_backreaction", "live", "full", "boundary_ports"),
        ("identical_history_clone_control", "real_growth", "live", "full", "boundary_ports"),
        ("random_same_level_cycle_baseline", "real_growth", "live", "full", "boundary_ports"),
        ("double_history_suffix_cycle", "real_growth", "live", "diagonal", "boundary_ports"),
        ("double_history_suffix_cycle", "real_growth", "live", "trace_scalar", "boundary_ports"),
        ("double_history_suffix_cycle", "real_growth", "handoff", "full", "boundary_ports"),
        ("double_history_suffix_cycle", "real_growth", "aging", "full", "boundary_ports"),
    ]

    def fmt(d: Optional[dict]) -> str:
        if d is None:
            return "missing"
        return (
            f"{d['mode']} | {d['control']} | {d['source_mode']} | {d['reduction']} | {d['fingerprint']}: "
            f"count={d['count']}, valid={d['valid']:.3f}, K={d['K']:.6g}, "
            f"id_match={d['identity_match']:.3f}, dir_match={d['directed_match']:.3f}, "
            f"cyclic={d['cyclic_match']:.3f}, refl={d['reflection_match']:.3f}, "
            f"best={d['best_cost']:.4g}, id_cost={d['identity_cost']:.4g}, dir_cost={d['directed_cost']:.4g}, "
            f"gap={d['gap']:.4g}, adv_dir_vs_id={d['adv_directed_vs_identity']:.4g}, "
            f"cycle_id={d['cycle_identity']:.3f}, cycle_nontriv={d['cycle_nontrivial']:.3f}"
        )

    lines = [
        "CNNA BOUNDARY-RESOLVED FACE-CORRESPONDENCE TRANSPORT",
        f"max_level={args.max_level}, mode={args.mode}, min_omega={args.min_omega}, scale={args.scale_mode}",
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
        "  The primary anti-smuggling fingerprint is boundary_ports: it uses DtN",
        "  port/vertex boundary profiles and edge-response magnitudes, but not J/K",
        "  as the matching feature. K/J fingerprints are audits only.",
        "  The script does not restrict correspondences to Z3 label shifts. It tests",
        "  all six S3 permutations and classifies the best one as identity, cyclic",
        "  plus/minus or reflection.",
        "  A positive boundary-resolved nontrivial transport would need real_growth",
        "  double-history cycles to prefer cyclic/reflection permutations over identity,",
        "  with identical-history and random controls not doing the same. Diagonal/trace",
        "  reductions should kill the noncommutative K-sector even if boundary labels",
        "  still have trivial identity structure.",
    ])
    summary = "\n".join(lines)
    (outdir / "SUMMARY.txt").write_text(summary, encoding="utf-8")

    results = f"""# Results: boundary-resolved face-correspondence transport

This test infers the face correspondence between double-history parent-fan nets
from boundary/port fingerprints rather than from face labels or a prescribed Z3
shift.

```text
{summary}
```

## Output files

- `boundary_face_cycle_rows.csv`
- `boundary_face_edge_rows.csv`
- `boundary_face_summary_main.csv`
- `boundary_face_summary_by_level.csv`
- `model_summaries.csv`
- `SUMMARY.txt`

## Status

The primary row is `double_history_suffix_cycle | real_growth | live | full |
boundary_ports`.  It is the cleanest anti-smuggling gate because it excludes
J/K from the correspondence feature.
"""
    (outdir / "RESULTS_boundary_resolved_face_correspondence_transport.md").write_text(results, encoding="utf-8")

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
    ap.add_argument("--scale-mode", default="unit", choices=["unit", "raw", "zscore"])
    ap.add_argument("--fingerprints", default="")
    ap.add_argument("--quick-suite", action="store_true")
    ap.add_argument("--outdir", type=Path, default=Path("boundary_resolved_face_correspondence_out_L6"))
    args = ap.parse_args()
    print(run(args))


if __name__ == "__main__":
    main()
