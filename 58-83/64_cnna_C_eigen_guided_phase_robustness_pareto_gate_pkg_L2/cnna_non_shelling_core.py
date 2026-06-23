#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import math
import shutil
import zipfile
from collections import defaultdict
from dataclasses import dataclass, field
from itertools import combinations, permutations
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

import numpy as np

EPS = 1e-12
Face = Tuple[int, int, int]
Tet = Tuple[int, int, int, int]
Simplex = Tuple[int, ...]


def unit(v: np.ndarray) -> np.ndarray:
    v = np.asarray(v, dtype=float).reshape(3)
    n = float(np.linalg.norm(v))
    if n < EPS:
        return np.array([0.0, 0.0, 1.0])
    return v / n


def frame_from_radial(radial: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    r = unit(radial)
    seed = np.array([0.0, 0.0, 1.0])
    if abs(float(np.dot(seed, r))) > 0.92:
        seed = np.array([1.0, 0.0, 0.0])
    e1 = unit(np.cross(r, seed))
    e2 = unit(np.cross(r, e1))
    return r, e1, e2


def fro(M: np.ndarray) -> float:
    return float(np.linalg.norm(np.asarray(M, dtype=float), ord="fro"))


def skew(M: np.ndarray) -> np.ndarray:
    M = np.asarray(M, dtype=float)
    return 0.5 * (M - M.T)


def sym(M: np.ndarray) -> np.ndarray:
    M = np.asarray(M, dtype=float)
    return 0.5 * (M + M.T)


def mean(xs: Iterable[float]) -> float:
    vals = [float(x) for x in xs if np.isfinite(float(x))]
    return float(np.mean(vals)) if vals else 0.0


def perc(xs: Iterable[float], q: float) -> float:
    vals = [float(x) for x in xs if np.isfinite(float(x))]
    return float(np.percentile(vals, q)) if vals else 0.0


@dataclass
class Node:
    id: int
    parent: Optional[int]
    level: int
    birth_order: int
    birth_time: int
    birth_g: float
    g: float
    pos: np.ndarray
    radial: np.ndarray
    e1: np.ndarray
    e2: np.ndarray
    children: List[int] = field(default_factory=list)


class DynamicProvenanceGrowth:
    def __init__(
        self,
        mode: str = "linear",
        growth_rule: str = "real_growth",
        branching: int = 3,
        base: float = 1.0,
        alpha_env: float = 0.22,
        ancestor_decay: float = 0.55,
        br_ancestor: float = 0.045,
        br_sibling: float = 0.035,
        transverse_amp: float = 0.42,
        radial_step: float = 1.0,
    ):
        if branching != 3:
            raise ValueError("This diagnostic uses ternary sibling fans.")
        if growth_rule not in {"real_growth", "symmetrized_birth", "no_backreaction"}:
            raise ValueError(growth_rule)
        self.mode = mode
        self.growth_rule = growth_rule
        self.branching = branching
        self.base = base
        self.alpha_env = alpha_env
        self.ancestor_decay = ancestor_decay
        self.br_ancestor = 0.0 if growth_rule == "no_backreaction" else br_ancestor
        self.br_sibling = 0.0 if growth_rule == "no_backreaction" else br_sibling
        self.transverse_amp = transverse_amp
        self.radial_step = radial_step
        self.nodes: Dict[int, Node] = {}
        self.t = 0
        self.next_id = 0
        self.directed_edges: Dict[Tuple[int, int], float] = defaultdict(float)
        self.birth_events: List[dict] = []
        r, e1, e2 = frame_from_radial(np.array([0.0, 0.0, 1.0]))
        root = self._new_node(None, 0, 0, 1.0, np.zeros(3), r, e1, e2)
        self.root = root.id

    def _new_node(
        self,
        parent: Optional[int],
        level: int,
        birth_order: int,
        birth_g: float,
        pos: np.ndarray,
        radial: np.ndarray,
        e1: np.ndarray,
        e2: np.ndarray,
    ) -> Node:
        n = Node(self.next_id, parent, level, birth_order, self.t, birth_g, birth_g, pos, radial, e1, e2)
        self.nodes[n.id] = n
        self.next_id += 1
        if parent is not None:
            self.nodes[parent].children.append(n.id)
        return n

    def parent_line(self, parent: int) -> List[int]:
        line: List[int] = []
        cur: Optional[int] = parent
        while cur is not None:
            line.append(cur)
            cur = self.nodes[cur].parent
        return line

    def address_tuple(self, node: int) -> Tuple[int, ...]:
        out: List[int] = []
        cur: Optional[int] = node
        while cur is not None and self.nodes[cur].parent is not None:
            out.append(self.nodes[cur].birth_order)
            cur = self.nodes[cur].parent
        return tuple(reversed(out))

    def birth_environment_load(self, parent: int, older_siblings: List[int]) -> float:
        env = 0.0
        for d, a in enumerate(self.parent_line(parent), start=1):
            env += self.nodes[a].g * (self.ancestor_decay ** (d - 1))
        if self.growth_rule != "symmetrized_birth":
            for s in older_siblings:
                env += self.nodes[s].g
        return env

    def child_conductance_from_env(self, env_load: float) -> float:
        if self.mode == "linear":
            return self.base + self.alpha_env * env_load
        if self.mode == "log":
            return self.base + self.alpha_env * math.log1p(env_load)
        if self.mode == "saturating":
            return self.base + self.alpha_env * (env_load / (1.0 + env_load))
        raise ValueError(self.mode)

    def child_position(self, parent: int, order: int, older_siblings: List[int]) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        p = self.nodes[parent]
        theta = 2.0 * math.pi * (order - 1) / 3.0
        twist = 0.37 * sum((i + 1) * x for i, x in enumerate(self.address_tuple(parent)))
        transverse = math.cos(theta + twist) * p.e1 + math.sin(theta + twist) * p.e2
        older_push = np.zeros(3)
        for s in older_siblings:
            older_push += unit(p.pos - self.nodes[s].pos)
        direction = unit(p.radial + self.transverse_amp * transverse + 0.08 * older_push)
        step = self.radial_step * (1.0 + 0.06 * (order - 2))
        pos = p.pos + step * direction
        r, e1, e2 = frame_from_radial(pos if np.linalg.norm(pos) > EPS else direction)
        return pos, r, e1, e2

    def add_child(self, parent: int, order: int) -> int:
        self.t += 1
        older = list(self.nodes[parent].children)
        env_load = self.birth_environment_load(parent, older)
        birth_g = self.child_conductance_from_env(env_load)
        pos, r, e1, e2 = self.child_position(parent, order, older)
        child = self._new_node(parent, self.nodes[parent].level + 1, order, birth_g, pos, r, e1, e2)
        c = child.id
        total_env = env_load + EPS
        for d, a in enumerate(self.parent_line(parent), start=1):
            contrib = self.nodes[a].g * (self.ancestor_decay ** (d - 1))
            self.directed_edges[(a, c)] += self.alpha_env * contrib / total_env * birth_g
        if self.growth_rule != "symmetrized_birth":
            for s in older:
                contrib = self.nodes[s].g
                self.directed_edges[(s, c)] += self.alpha_env * contrib / total_env * birth_g
        for d, a in enumerate(self.parent_line(parent), start=1):
            delta = self.br_ancestor * birth_g / (d * d)
            self.nodes[a].g += delta
            self.directed_edges[(c, a)] += delta
        for s in older:
            delta = self.br_sibling * birth_g
            self.nodes[s].g += delta
            self.directed_edges[(c, s)] += delta
        self.birth_events.append({
            "t": self.t,
            "parent": parent,
            "child": c,
            "order": order,
            "older_siblings": list(older),
            "env_load": env_load,
            "birth_g": birth_g,
            "level": child.level,
        })
        return c

    def grow_level(self, frontier: List[int]) -> List[int]:
        nxt: List[int] = []
        for p in frontier:
            for k in range(1, 4):
                nxt.append(self.add_child(p, k))
        return nxt

    def grow(self, max_level: int) -> None:
        frontier = [self.root]
        for _ in range(max_level):
            frontier = self.grow_level(frontier)

    def completed_parent_ids(self) -> List[int]:
        return [i for i, n in self.nodes.items() if len(n.children) == 3]

    def child_ids_ordered(self, parent: int) -> List[int]:
        return sorted(self.nodes[parent].children, key=lambda c: self.nodes[c].birth_order)


@dataclass
class SimplicialComplex:
    name: str
    vertices: set[int] = field(default_factory=set)
    tets: List[Tet] = field(default_factory=list)
    face_birth: Dict[Face, int] = field(default_factory=dict)

    def clone(self, name: Optional[str] = None) -> "SimplicialComplex":
        return SimplicialComplex(name or self.name, set(self.vertices), list(self.tets), dict(self.face_birth))

    def add_tet(self, tet: Iterable[int], birth_time: int = 0) -> bool:
        tt = tuple(sorted(set(int(x) for x in tet)))
        if len(tt) != 4:
            return False
        if tt in self.tets:
            return False
        self.tets.append(tt)
        self.vertices.update(tt)
        for f in faces_of_tet(tt):
            self.face_birth.setdefault(f, birth_time)
        return True

    def faces(self) -> List[Face]:
        fs = set()
        for t in self.tets:
            fs.update(faces_of_tet(t))
        return sorted(fs)

    def edges(self) -> List[Tuple[int, int]]:
        es = set()
        for f in self.faces():
            a, b, c = f
            es.add(tuple(sorted((a, b))))
            es.add(tuple(sorted((a, c))))
            es.add(tuple(sorted((b, c))))
        return sorted(es)

    def face_occupancy(self) -> Dict[Face, int]:
        occ: Dict[Face, int] = defaultdict(int)
        for t in self.tets:
            for f in faces_of_tet(t):
                occ[f] += 1
        return dict(occ)


def faces_of_tet(t: Tet) -> List[Face]:
    a, b, c, d = tuple(sorted(t))
    return [tuple(sorted(x)) for x in [(b, c, d), (a, c, d), (a, b, d), (a, b, c)]]


def update_face_maps(K: SimplicialComplex, faces_by_vertex: Dict[int, set[Face]], occ: Dict[Face, int], tet: Tet) -> None:
    for f in faces_of_tet(tet):
        occ[f] = occ.get(f, 0) + 1
        for v in f:
            faces_by_vertex.setdefault(v, set()).add(f)


def choose_boundary_face_for_parent(model: DynamicProvenanceGrowth, faces_by_vertex: Dict[int, set[Face]], occ: Dict[Face, int], parent: int) -> Optional[Face]:
    candidates = [f for f in faces_by_vertex.get(parent, set()) if occ.get(f, 0) == 1]
    if not candidates:
        return None
    pr = model.nodes[parent].radial
    def score(f: Face) -> Tuple[float, int]:
        centroid = sum((model.nodes[v].pos for v in f), np.zeros(3)) / 3.0
        outward = unit(centroid - model.nodes[parent].pos)
        return (float(np.dot(outward, pr)), -sum(model.nodes[v].birth_time for v in f))
    return max(candidates, key=score)


def build_dynamic_outward_ngf_complex(model: DynamicProvenanceGrowth) -> SimplicialComplex:
    K = SimplicialComplex("dynamic_outward_ngf_base")
    occ: Dict[Face, int] = {}
    faces_by_vertex: Dict[int, set[Face]] = defaultdict(set)
    root_seeded = False
    for ev in model.birth_events:
        parent = int(ev["parent"])
        child = int(ev["child"])
        if not root_seeded and len(model.nodes[model.root].children) == 3:
            ch = model.child_ids_ordered(model.root)
            tet = tuple(sorted((model.root, ch[0], ch[1], ch[2])))
            if K.add_tet(tet, birth_time=max(model.nodes[c].birth_time for c in ch)):
                update_face_maps(K, faces_by_vertex, occ, tet)
            root_seeded = True
        if child in K.vertices:
            continue
        face = choose_boundary_face_for_parent(model, faces_by_vertex, occ, parent)
        if face is None:
            continue
        tet = tuple(sorted((*face, child)))
        if any(occ.get(ff, 0) >= 2 for ff in faces_of_tet(tet)):
            continue
        if K.add_tet(tet, birth_time=int(ev["t"])):
            update_face_maps(K, faces_by_vertex, occ, tet)
    return K


def gf2_rank(A: np.ndarray) -> int:
    A = (np.asarray(A, dtype=np.uint8) & 1).copy()
    m, n = A.shape
    rank = 0
    row = 0
    for col in range(n):
        piv = None
        for r in range(row, m):
            if A[r, col]:
                piv = r
                break
        if piv is None:
            continue
        if piv != row:
            A[[row, piv]] = A[[piv, row]]
        for r in range(m):
            if r != row and A[r, col]:
                A[r, :] ^= A[row, :]
        rank += 1
        row += 1
        if row == m:
            break
    return rank


def boundary_matrix_mod2(high: List[Simplex], low: List[Simplex]) -> np.ndarray:
    idx = {s: i for i, s in enumerate(low)}
    B = np.zeros((len(low), len(high)), dtype=np.uint8)
    for j, s in enumerate(high):
        for k in range(len(s)):
            face = tuple(sorted(s[:k] + s[k + 1:]))
            if face in idx:
                B[idx[face], j] ^= 1
    return B


def boundary_matrix_real(high: List[Simplex], low: List[Simplex]) -> np.ndarray:
    idx = {s: i for i, s in enumerate(low)}
    B = np.zeros((len(low), len(high)), dtype=float)
    for j, s in enumerate(high):
        for k in range(len(s)):
            face = tuple(sorted(s[:k] + s[k + 1:]))
            if face in idx:
                B[idx[face], j] += -1.0 if k % 2 else 1.0
    return B


def topology(K: SimplicialComplex) -> dict:
    verts = sorted(K.vertices)
    edges = K.edges()
    faces = K.faces()
    tets = sorted(K.tets)
    C = [len(verts), len(edges), len(faces), len(tets)]
    B1 = boundary_matrix_mod2([tuple(e) for e in edges], [(v,) for v in verts]) if edges else np.zeros((len(verts), 0), dtype=np.uint8)
    B2 = boundary_matrix_mod2([tuple(f) for f in faces], [tuple(e) for e in edges]) if faces else np.zeros((len(edges), 0), dtype=np.uint8)
    B3 = boundary_matrix_mod2([tuple(t) for t in tets], [tuple(f) for f in faces]) if tets else np.zeros((len(faces), 0), dtype=np.uint8)
    r1, r2, r3 = gf2_rank(B1), gf2_rank(B2), gf2_rank(B3)
    occ = K.face_occupancy()
    boundary_faces = sum(1 for n in occ.values() if n == 1)
    saturated_faces = sum(1 for n in occ.values() if n == 2)
    overfull_faces = sum(1 for n in occ.values() if n > 2)
    return {
        "vertices": C[0],
        "edges": C[1],
        "faces": C[2],
        "tets": C[3],
        "euler": C[0] - C[1] + C[2] - C[3],
        "beta0": C[0] - r1,
        "beta1": C[1] - r1 - r2,
        "beta2": C[2] - r2 - r3,
        "beta3": C[3] - r3,
        "boundary_faces": boundary_faces,
        "saturated_faces": saturated_faces,
        "overfull_faces": overfull_faces,
        "boundary_fraction": boundary_faces / (len(occ) + EPS),
        "saturated_fraction": saturated_faces / (len(occ) + EPS),
        "manifold_face_fraction": (boundary_faces + saturated_faces) / (len(occ) + EPS),
    }


def vertex_operator(model: DynamicProvenanceGrowth, node: int, source: str = "live") -> np.ndarray:
    n = model.nodes[node]
    r, e1, e2 = n.radial, n.e1, n.e2
    order_phase = 2.0 * math.pi * (n.birth_order - 1) / 3.0 if n.parent is not None else 0.0
    q = math.cos(order_phase) * e1 + math.sin(order_phase) * e2
    h = unit(0.7 * r + 0.3 * q)
    birth = n.birth_g
    live = n.g
    aging = max(0.0, live - birth)
    if source == "record":
        a, b, c = birth, 0.22 * birth, 0.08 * birth
    elif source == "live":
        a, b, c = live, 0.25 * birth + 0.55 * aging, 0.12 * live
    elif source == "handoff":
        inc = 0.0
        if n.parent is not None:
            inc = model.directed_edges.get((n.parent, node), 0.0) + model.directed_edges.get((node, n.parent), 0.0)
        a, b, c = birth + inc, 0.18 * live + inc, 0.15 * inc + 0.05 * birth
    elif source == "aging":
        a, b, c = aging + 0.1 * birth, 0.6 * aging + 0.03 * birth, 0.3 * aging
    else:
        raise ValueError(source)
    M = a * np.outer(r, r) + b * np.outer(q, q) + c * np.outer(h, h) + 0.04 * birth * np.eye(3)
    return sym(M)


def face_K(model: DynamicProvenanceGrowth, face: Face, source: str = "live") -> np.ndarray:
    a, b, c = face
    Sa = vertex_operator(model, a, source)
    Sb = vertex_operator(model, b, source)
    Sc = vertex_operator(model, c, source)
    Aab = Sb - Sa
    Abc = Sc - Sb
    return skew(Aab @ Abc - Abc @ Aab)


def exactness_metrics(model: DynamicProvenanceGrowth, K: SimplicialComplex, source: str = "live") -> dict:
    faces = K.faces()
    edges = K.edges()
    tets = sorted(K.tets)
    if not faces:
        return {"K_mean": 0.0, "closed_ratio": 0.0, "exact_ratio": 0.0, "harmonic_ratio": 0.0}
    weights = np.array([fro(face_K(model, f, source)) for f in faces], dtype=float)
    norm = float(np.linalg.norm(weights)) + EPS
    B2 = boundary_matrix_real([tuple(f) for f in faces], [tuple(e) for e in edges]) if edges else np.zeros((0, len(faces)))
    B3 = boundary_matrix_real([tuple(t) for t in tets], [tuple(f) for f in faces]) if tets else np.zeros((len(faces), 0))
    d2 = B3.T
    closed_res = float(np.linalg.norm(d2 @ weights)) if d2.size else 0.0
    A_exact = B2.T if B2.size else np.zeros((len(faces), 0))
    if A_exact.size and A_exact.shape[1] > 0:
        sol, *_ = np.linalg.lstsq(A_exact, weights, rcond=None)
        exact = A_exact @ sol
    else:
        exact = np.zeros_like(weights)
    exact_res = float(np.linalg.norm(weights - exact))
    harmonic_ratio = 0.0
    if B2.size and B3.size:
        L2 = B2.T @ B2 + B3 @ B3.T
        vals, vecs = np.linalg.eigh(L2)
        mask = vals < 1e-9
        if np.any(mask):
            H = vecs[:, mask]
            hproj = H @ (H.T @ weights)
            harmonic_ratio = float(np.linalg.norm(hproj) / norm)
    return {
        "K_mean": mean(weights),
        "K_p95": perc(weights, 95),
        "closed_ratio": closed_res / norm,
        "exact_residual_ratio": exact_res / norm,
        "harmonic_ratio": harmonic_ratio,
    }


def attachment_face_count(K: SimplicialComplex, tet: Tet) -> Tuple[int, int, List[Face]]:
    occ = K.face_occupancy()
    boundary = []
    existing = 0
    for f in faces_of_tet(tet):
        if occ.get(f, 0) > 0:
            existing += 1
        if occ.get(f, 0) == 1:
            boundary.append(f)
    return len(boundary), existing, boundary


def add_tets_result(K: SimplicialComplex, tets: List[Tet], name: str = "candidate") -> Tuple[Optional[SimplicialComplex], str]:
    L = K.clone(name)
    occ = L.face_occupancy()
    for tet in tets:
        tt = tuple(sorted(set(tet)))
        if len(tt) != 4:
            return None, "degenerate_tet"
        if tt in L.tets:
            return None, "duplicate_tet"
        for f in faces_of_tet(tt):
            if occ.get(f, 0) >= 2:
                return None, "face_occupancy_over_2"
        if not L.add_tet(tt, birth_time=0):
            return None, "add_failed"
        for f in faces_of_tet(tt):
            occ[f] = occ.get(f, 0) + 1
    return L, "ok"


def face_centroid(model: DynamicProvenanceGrowth, f: Face) -> np.ndarray:
    return sum((model.nodes[v].pos for v in f), np.zeros(3)) / 3.0


def face_normal_proxy(model: DynamicProvenanceGrowth, f: Face) -> np.ndarray:
    a, b, c = [model.nodes[v].pos for v in f]
    return unit(np.cross(b - a, c - a))


def address_similarity(model: DynamicProvenanceGrowth, f: Face, g: Face) -> float:
    def parent_multiset(face: Face) -> List[int]:
        return [model.nodes[v].parent if model.nodes[v].parent is not None else -1 for v in face]
    pf, pg = parent_multiset(f), parent_multiset(g)
    same_parent = len(set(pf).intersection(pg))
    def suffix(face: Face) -> set[Tuple[int, ...]]:
        return {model.address_tuple(v)[-2:] for v in face}
    overlap = len(suffix(f).intersection(suffix(g)))
    return 0.4 * same_parent + 0.25 * overlap


def transverse_complementarity(model: DynamicProvenanceGrowth, f: Face, g: Face) -> float:
    nf = face_normal_proxy(model, f)
    ng = face_normal_proxy(model, g)
    cf = unit(face_centroid(model, f))
    cg = unit(face_centroid(model, g))
    normal_opposition = max(0.0, -float(np.dot(nf, ng)))
    radial_alignment = max(0.0, float(np.dot(cf, cg)))
    return normal_opposition + 0.35 * radial_alignment


def directed_face_coupling(model: DynamicProvenanceGrowth, f: Face, g: Face) -> Tuple[float, float]:
    fg = 0.0
    gf = 0.0
    for a in f:
        for b in g:
            fg += model.directed_edges.get((a, b), 0.0)
            gf += model.directed_edges.get((b, a), 0.0)
    imbalance = abs(fg - gf)
    total = fg + gf
    return total, imbalance


def face_pair_response_score(model: DynamicProvenanceGrowth, f: Face, g: Face, source: str = "live") -> dict:
    Kf = fro(face_K(model, f, source))
    Kg = fro(face_K(model, g, source))
    total_coupling, imbalance = directed_face_coupling(model, f, g)
    cf = face_centroid(model, f)
    cg = face_centroid(model, g)
    distance = float(np.linalg.norm(cf - cg))
    trans_comp = transverse_complementarity(model, f, g)
    addr = address_similarity(model, f, g)
    radial_shell_match = 1.0 / (1.0 + abs(float(np.linalg.norm(cf) - np.linalg.norm(cg))))
    score = (
        1.00 * (Kf + Kg)
        + 2.00 * total_coupling
        + 1.50 * imbalance
        + 0.75 * trans_comp
        + 0.50 * addr
        + 0.40 * radial_shell_match
        - 0.03 * distance
    )
    return {
        "response_score": float(score),
        "K_pair_norm": float(Kf + Kg),
        "directed_coupling": float(total_coupling),
        "directed_imbalance": float(imbalance),
        "transverse_complementarity": float(trans_comp),
        "address_similarity": float(addr),
        "radial_shell_match": float(radial_shell_match),
        "centroid_distance": float(distance),
    }


def single_face_vertex_score(model: DynamicProvenanceGrowth, f: Face, x: int, source: str = "live") -> dict:
    fx = tuple(sorted((*f, x)))
    face_score = fro(face_K(model, f, source))
    gvals = [model.nodes[v].g for v in f] + [model.nodes[x].g]
    spread = max(gvals) - min(gvals)
    parent_coupling = 0.0
    for v in f:
        parent_coupling += model.directed_edges.get((v, x), 0.0) + model.directed_edges.get((x, v), 0.0)
    dist = float(np.linalg.norm(model.nodes[x].pos - face_centroid(model, f)))
    score = face_score + 2.0 * parent_coupling + 0.45 * spread - 0.04 * dist
    return {
        "response_score": float(score),
        "K_pair_norm": float(face_score),
        "directed_coupling": float(parent_coupling),
        "directed_imbalance": 0.0,
        "transverse_complementarity": 0.0,
        "address_similarity": 0.0,
        "radial_shell_match": 0.0,
        "centroid_distance": dist,
    }


def quotient_complex(K: SimplicialComplex, mapping: Dict[int, int]) -> Tuple[Optional[SimplicialComplex], str]:
    parent: Dict[int, int] = {v: v for v in K.vertices}
    def find(x: int) -> int:
        parent.setdefault(x, x)
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    def union(a: int, b: int) -> None:
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[max(ra, rb)] = min(ra, rb)
    for a, b in mapping.items():
        union(a, b)
    L = SimplicialComplex("quotient_candidate")
    collapsed = 0
    for t in K.tets:
        qt = tuple(sorted({find(v) for v in t}))
        if len(qt) < 4:
            collapsed += 1
            continue
        L.add_tet(qt, birth_time=0)
    if not L.tets:
        return None, "all_tets_collapsed"
    occ = L.face_occupancy()
    if any(n > 2 for n in occ.values()):
        return L, "overfull_after_quotient"
    if collapsed > max(2, 0.15 * len(K.tets)):
        return L, "too_many_collapsed_tets"
    return L, "ok"


def parity_of_permutation(p: Tuple[int, int, int]) -> int:
    inv = 0
    for i in range(len(p)):
        for j in range(i + 1, len(p)):
            if p[i] > p[j]:
                inv += 1
    return -1 if inv % 2 else 1


def summarize_delta(topo0: dict, topo1: dict) -> dict:
    keys = ["beta0", "beta1", "beta2", "beta3", "boundary_faces", "saturated_faces", "overfull_faces", "euler"]
    return {"delta_" + k: topo1.get(k, 0) - topo0.get(k, 0) for k in keys}


def enumerate_moves(
    model: DynamicProvenanceGrowth,
    K: SimplicialComplex,
    *,
    source: str = "live",
    max_boundary_faces: int = 90,
    max_single_vertices: int = 12,
    max_pair_candidates: int = 2500,
    max_rows: int = 5000,
) -> Tuple[List[dict], dict]:
    topo0 = topology(K)
    occ = K.face_occupancy()
    boundary = [f for f, n in occ.items() if n == 1]
    def boundary_base_score(f: Face) -> float:
        c = face_centroid(model, f)
        radial = float(np.linalg.norm(c))
        return fro(face_K(model, f, source)) + 0.03 * radial + 0.01 * sum(model.nodes[v].birth_time for v in f)
    boundary = sorted(boundary, key=boundary_base_score, reverse=True)[:max_boundary_faces]
    rows: List[dict] = []
    cid = 0

    for f in boundary:
        fset = set(f)
        c = face_centroid(model, f)
        vertices = sorted(K.vertices, key=lambda x: np.linalg.norm(model.nodes[x].pos - c))[:max_single_vertices]
        for x in vertices:
            if x in fset:
                continue
            tet = tuple(sorted((*f, x)))
            if tet in K.tets:
                continue
            newfaces = faces_of_tet(tet)
            if any(occ.get(ff, 0) >= 2 for ff in newfaces):
                cls = "illegal"
                reason = "face_occupancy_over_2"
                L = None
            else:
                attach_count, existing_count, _ = attachment_face_count(K, tet)
                if attach_count == 0:
                    cls = "illegal"
                    reason = "detached_tet"
                    L = None
                elif attach_count == 1:
                    cls = "shelling_disk_move"
                    reason = "single_boundary_face_attachment"
                    L, reason = add_tets_result(K, [tet], "shelling_candidate")
                elif attach_count in {2, 3, 4}:
                    cls = "cap_move"
                    reason = f"connected_{attach_count}_face_boundary_cap"
                    L, reason = add_tets_result(K, [tet], "cap_candidate")
                else:
                    cls = "illegal"
                    reason = "unclassified_single_tet"
                    L = None
            score = single_face_vertex_score(model, f, x, source)
            topo1 = topology(L) if L is not None else topo0
            delta = summarize_delta(topo0, topo1)
            rows.append({
                "candidate_id": cid,
                "move_class": cls,
                "status": reason,
                "face_a": " ".join(map(str, f)),
                "face_b": "",
                "new_tets": "|".join(" ".join(map(str, t)) for t in ([tet] if cls != "illegal" else [])),
                "attach_components": 1 if cls in {"shelling_disk_move", "cap_move"} else 0,
                "attach_boundary_faces": attachment_face_count(K, tet)[0] if len(set(tet)) == 4 else 0,
                **score,
                **delta,
                "new_beta1": topo1["beta1"],
                "new_beta2": topo1["beta2"],
                "new_boundary_fraction": topo1["boundary_fraction"],
            })
            cid += 1
            if len(rows) >= max_rows:
                break
        if len(rows) >= max_rows:
            break

    pair_pre: List[Tuple[float, Face, Face, dict]] = []
    for f, g in combinations(boundary, 2):
        if set(f).intersection(g):
            # shared-vertex pairs are mostly local caps, not global complement closures.
            continue
        score = face_pair_response_score(model, f, g, source)
        pair_pre.append((score["response_score"], f, g, score))
    pair_pre.sort(reverse=True, key=lambda z: z[0])
    pair_pre = pair_pre[:max_pair_candidates]

    for _, f, g, score in pair_pre:
        # Handle bridge: attach a triangular-prism-like 3-tet bridge to two disconnected boundary faces.
        # This is not a CNNA claim; it audits whether the current boundary even contains non-shelling bridge candidates.
        best_bridge: Optional[Tuple[List[Tet], Optional[SimplicialComplex], str, dict]] = None
        for perm in permutations(g):
            d, e, h = perm
            a, b, c = f
            tets = [tuple(sorted((a, b, c, d))), tuple(sorted((b, c, d, e))), tuple(sorted((c, d, e, h)))]
            L, reason = add_tets_result(K, tets, "handle_bridge_candidate")
            if L is None:
                topo1 = topo0
            else:
                topo1 = topology(L)
            delta = summarize_delta(topo0, topo1)
            topo_gain = 2.0 * max(0, delta["delta_beta1"]) + 2.5 * max(0, delta["delta_beta2"]) - 0.02 * max(0, delta["delta_boundary_faces"])
            cand = (tets, L, reason, {**delta, "topology_gain": topo_gain})
            if best_bridge is None or cand[3]["topology_gain"] > best_bridge[3]["topology_gain"]:
                best_bridge = cand
        if best_bridge is not None:
            tets, L, reason, delta = best_bridge
            topo1 = topology(L) if L is not None else topo0
            cls = "handle_candidate" if L is not None and reason == "ok" else "illegal"
            rows.append({
                "candidate_id": cid,
                "move_class": cls,
                "status": reason,
                "face_a": " ".join(map(str, f)),
                "face_b": " ".join(map(str, g)),
                "new_tets": "|".join(" ".join(map(str, t)) for t in tets),
                "attach_components": 2,
                "attach_boundary_faces": 2,
                **score,
                **{k: v for k, v in delta.items() if k.startswith("delta_")},
                "new_beta1": topo1["beta1"],
                "new_beta2": topo1["beta2"],
                "new_boundary_fraction": topo1["boundary_fraction"],
            })
            cid += 1

        # Quotient: identify two boundary faces with opposite orientation. Test all odd permutations.
        best_quot: Optional[Tuple[Tuple[int, int, int], Optional[SimplicialComplex], str, dict]] = None
        for p in permutations((0, 1, 2)):
            if parity_of_permutation(p) != -1:
                continue
            mapping = {g[p[i]]: f[i] for i in range(3)}
            L, reason = quotient_complex(K, mapping)
            topo1 = topology(L) if L is not None else topo0
            delta = summarize_delta(topo0, topo1)
            topo_gain = 2.0 * max(0, delta["delta_beta1"]) + 2.5 * max(0, delta["delta_beta2"]) + 0.04 * max(0, -delta["delta_boundary_faces"]) - 4.0 * max(0, topo1.get("overfull_faces", 0))
            cand = (p, L, reason, {**delta, "topology_gain": topo_gain})
            if best_quot is None or cand[3]["topology_gain"] > best_quot[3]["topology_gain"]:
                best_quot = cand
        if best_quot is not None:
            p, L, reason, delta = best_quot
            topo1 = topology(L) if L is not None else topo0
            legal = L is not None and reason == "ok"
            rows.append({
                "candidate_id": cid,
                "move_class": "quotient_candidate" if legal else "illegal",
                "status": reason,
                "face_a": " ".join(map(str, f)),
                "face_b": " ".join(map(str, g)) + f" perm={p}",
                "new_tets": "quotient_face_pairing",
                "attach_components": 2,
                "attach_boundary_faces": 2,
                **score,
                **{k: v for k, v in delta.items() if k.startswith("delta_")},
                "new_beta1": topo1["beta1"],
                "new_beta2": topo1["beta2"],
                "new_boundary_fraction": topo1["boundary_fraction"],
            })
            cid += 1
        if len(rows) >= max_rows:
            break

    legal_rows = [r for r in rows if r["move_class"] != "illegal"]
    ranked = sorted(legal_rows, key=lambda r: r["response_score"], reverse=True)
    for rank, row in enumerate(ranked, start=1):
        row["response_rank_legal"] = rank
    rank_map = {row["candidate_id"]: row.get("response_rank_legal", "") for row in ranked}
    for row in rows:
        row["response_rank_legal"] = rank_map.get(row["candidate_id"], "")
    counts: Dict[str, int] = defaultdict(int)
    for row in rows:
        counts[row["move_class"]] += 1
    best_by_class = {}
    for cls in sorted(counts):
        sub = [r for r in rows if r["move_class"] == cls]
        sub.sort(key=lambda r: r["response_score"], reverse=True)
        if sub:
            best_by_class[cls] = sub[0]
    top_handle_rank = min([int(r["response_rank_legal"]) for r in rows if r["move_class"] == "handle_candidate" and r["response_rank_legal"] != ""], default=None)
    top_quot_rank = min([int(r["response_rank_legal"]) for r in rows if r["move_class"] == "quotient_candidate" and r["response_rank_legal"] != ""], default=None)
    topological_effective = [r for r in rows if r["move_class"] in {"handle_candidate", "quotient_candidate"} and (r["delta_beta1"] > 0 or r["delta_beta2"] > 0 or r["delta_boundary_faces"] < 0)]
    audit = {
        "base_topology": topo0,
        "base_exactness": exactness_metrics(model, K, source),
        "candidate_counts": dict(counts),
        "legal_candidate_count": len(legal_rows),
        "top_handle_rank": top_handle_rank,
        "top_quotient_rank": top_quot_rank,
        "non_shelling_candidates_exist": counts.get("handle_candidate", 0) + counts.get("quotient_candidate", 0) > 0,
        "topologically_effective_non_shelling_count": len(topological_effective),
        "best_by_class": {k: compact_candidate(v) for k, v in best_by_class.items()},
        "interpretation_flags": {
            "current_growth_is_shelling_like": topo0["beta1"] == 0 and topo0["beta2"] == 0,
            "handle_or_quotient_exists": counts.get("handle_candidate", 0) + counts.get("quotient_candidate", 0) > 0,
            "response_ranks_non_shelling_top10": any((r is not None and r <= 10) for r in [top_handle_rank, top_quot_rank]),
            "missing_operation_if_not_applied": counts.get("handle_candidate", 0) + counts.get("quotient_candidate", 0) > 0,
        },
    }
    return rows, audit


def compact_candidate(row: dict) -> dict:
    keep = [
        "candidate_id", "move_class", "status", "response_score", "response_rank_legal",
        "face_a", "face_b", "attach_components", "delta_beta1", "delta_beta2",
        "delta_boundary_faces", "new_beta1", "new_beta2", "new_boundary_fraction",
        "K_pair_norm", "directed_coupling", "directed_imbalance", "transverse_complementarity",
        "address_similarity", "radial_shell_match", "centroid_distance",
    ]
    return {k: row.get(k, "") for k in keep}


def write_csv(path: Path, rows: List[dict]) -> None:
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    fields = list(rows[0].keys())
    for r in rows[1:]:
        for k in r.keys():
            if k not in fields:
                fields.append(k)
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for r in rows:
            w.writerow(r)



import re


def parse_face_string(s: str) -> Face:
    nums = [int(x) for x in re.findall(r"-?\d+", s.split("perm=")[0])[:3]]
    if len(nums) != 3:
        raise ValueError(f"cannot parse face: {s}")
    return tuple(sorted(nums))


def parse_tets_string(s: str) -> List[Tet]:
    out: List[Tet] = []
    if not s or s == "quotient_face_pairing":
        return out
    for part in s.split("|"):
        nums = [int(x) for x in part.split()]
        if len(nums) == 4:
            out.append(tuple(sorted(nums)))
    return out


def parse_perm_string(s: str) -> Tuple[int, int, int]:
    m = re.search(r"perm=\(([^)]*)\)", s)
    if not m:
        raise ValueError(f"cannot parse perm: {s}")
    nums = [int(x.strip()) for x in m.group(1).split(',')]
    if len(nums) != 3:
        raise ValueError(f"bad perm: {s}")
    return tuple(nums)  # type: ignore[return-value]


def edge_link_cycle_metrics(K: SimplicialComplex) -> dict:
    tets = sorted(K.tets)
    edges = K.edges()
    cycle_count = 0
    candidate_count = 0
    bad_link_count = 0
    for e in edges:
        a, b = e
        link_vertices = set()
        link_edges = set()
        for t in tets:
            if a in t and b in t:
                others = [x for x in t if x not in {a, b}]
                if len(others) == 2:
                    u, v = sorted(others)
                    link_vertices.add(u)
                    link_vertices.add(v)
                    link_edges.add((u, v))
        if len(link_vertices) >= 3:
            candidate_count += 1
            adj: Dict[int, set[int]] = {v: set() for v in link_vertices}
            for u, v in link_edges:
                adj.setdefault(u, set()).add(v)
                adj.setdefault(v, set()).add(u)
            seen = set()
            stack = [next(iter(link_vertices))]
            while stack:
                x = stack.pop()
                if x in seen:
                    continue
                seen.add(x)
                stack.extend(sorted(adj.get(x, set()) - seen))
            connected = len(seen) == len(link_vertices)
            all_degree_two = all(len(adj.get(v, set())) == 2 for v in link_vertices)
            has_cycle = connected and len(link_edges) >= len(link_vertices)
            if has_cycle:
                cycle_count += 1
            if not all_degree_two and has_cycle:
                bad_link_count += 1
    return {
        "edge_link_cycle_count": cycle_count,
        "edge_link_candidate_count": candidate_count,
        "edge_link_cycle_fraction": cycle_count / (candidate_count + EPS),
        "edge_link_noncircle_cycle_count": bad_link_count,
    }


def full_metrics(model: DynamicProvenanceGrowth, K: SimplicialComplex, source: str) -> dict:
    topo = topology(K)
    exact = exactness_metrics(model, K, source)
    links = edge_link_cycle_metrics(K)
    return {**topo, **links, **exact}


def apply_candidate_row(K: SimplicialComplex, row: dict) -> Tuple[Optional[SimplicialComplex], str, str]:
    cls = row.get("move_class", "")
    if cls in {"shelling_disk_move", "cap_move", "handle_candidate"}:
        tets = parse_tets_string(str(row.get("new_tets", "")))
        if not tets:
            return None, "no_tets", ""
        L, reason = add_tets_result(K, tets, f"applied_{cls}")
        return L, reason, "|".join(" ".join(map(str, t)) for t in tets)
    if cls == "quotient_candidate":
        f = parse_face_string(str(row.get("face_a", "")))
        g = parse_face_string(str(row.get("face_b", "")))
        p = parse_perm_string(str(row.get("face_b", "")))
        mapping = {g[p[i]]: f[i] for i in range(3)}
        L, reason = quotient_complex(K, mapping)
        return L, reason, json.dumps({str(k): v for k, v in mapping.items()}, sort_keys=True)
    return None, "unsupported_move_class", ""


def move_topology_gain(row: dict) -> float:
    return (
        4.0 * max(0.0, float(row.get("delta_beta2", 0)))
        + 3.0 * max(0.0, float(row.get("delta_beta1", 0)))
        + 0.03 * max(0.0, -float(row.get("delta_boundary_faces", 0)))
        + 0.002 * float(row.get("response_score", 0.0))
    )


def pick_rows(rows: List[dict]) -> List[Tuple[str, str, dict]]:
    picked: List[Tuple[str, str, dict]] = []
    for cls in ["shelling_disk_move", "cap_move", "handle_candidate", "quotient_candidate"]:
        legal = [r for r in rows if r.get("move_class") == cls and r.get("status") == "ok"]
        if not legal:
            continue
        by_response = sorted(legal, key=lambda r: int(r.get("response_rank_legal") or 10**9))[0]
        picked.append((cls, "top_response", by_response))
        by_topology = sorted(legal, key=move_topology_gain, reverse=True)[0]
        if by_topology.get("candidate_id") != by_response.get("candidate_id"):
            picked.append((cls, "top_topology_gain", by_topology))
    # also apply the best legal non-shelling candidate by response and by topology, even if already selected
    nonshell = [r for r in rows if r.get("move_class") in {"handle_candidate", "quotient_candidate"} and r.get("status") == "ok"]
    if nonshell:
        picked.append(("non_shelling", "top_response", sorted(nonshell, key=lambda r: int(r.get("response_rank_legal") or 10**9))[0]))
        picked.append(("non_shelling", "top_topology_gain", sorted(nonshell, key=move_topology_gain, reverse=True)[0]))
    seen = set()
    unique: List[Tuple[str, str, dict]] = []
    for cls, kind, row in picked:
        key = (cls, kind, row.get("candidate_id"))
        if key in seen:
            continue
        seen.add(key)
        unique.append((cls, kind, row))
    return unique


def audit_and_apply(growth_rule: str, args: argparse.Namespace, out: Path) -> dict:
    model = DynamicProvenanceGrowth(mode=args.mode, growth_rule=growth_rule, transverse_amp=args.transverse_amp)
    model.grow(args.max_level)
    K = build_dynamic_outward_ngf_complex(model)
    rows, audit = enumerate_moves(
        model,
        K,
        source=args.source,
        max_boundary_faces=args.max_boundary_faces,
        max_single_vertices=args.max_single_vertices,
        max_pair_candidates=args.max_pair_candidates,
        max_rows=args.max_rows,
    )
    sub = out / growth_rule
    sub.mkdir(parents=True, exist_ok=True)
    write_csv(sub / "move_candidates.csv", sorted(rows, key=lambda r: (r.get("response_rank_legal") == "", r.get("response_rank_legal") or 10**9)))
    base_metrics = full_metrics(model, K, args.source)
    applied_rows: List[dict] = []
    for selection_class, selection_kind, row in pick_rows(rows):
        L, reason, payload = apply_candidate_row(K, row)
        if L is None:
            new_metrics = base_metrics
            applied_ok = False
        else:
            new_metrics = full_metrics(model, L, args.source)
            applied_ok = reason == "ok"
        applied_rows.append({
            "growth_rule": growth_rule,
            "selection_class": selection_class,
            "selection_kind": selection_kind,
            "candidate_id": row.get("candidate_id"),
            "move_class": row.get("move_class"),
            "apply_status": reason,
            "applied_ok": applied_ok,
            "payload": payload,
            "response_rank_legal": row.get("response_rank_legal"),
            "response_score": row.get("response_score"),
            "candidate_delta_beta1": row.get("delta_beta1"),
            "candidate_delta_beta2": row.get("delta_beta2"),
            "candidate_delta_boundary_faces": row.get("delta_boundary_faces"),
            "base_beta0": base_metrics["beta0"],
            "base_beta1": base_metrics["beta1"],
            "base_beta2": base_metrics["beta2"],
            "base_beta3": base_metrics["beta3"],
            "new_beta0": new_metrics["beta0"],
            "new_beta1": new_metrics["beta1"],
            "new_beta2": new_metrics["beta2"],
            "new_beta3": new_metrics["beta3"],
            "delta_beta1_actual": new_metrics["beta1"] - base_metrics["beta1"],
            "delta_beta2_actual": new_metrics["beta2"] - base_metrics["beta2"],
            "base_boundary_fraction": base_metrics["boundary_fraction"],
            "new_boundary_fraction": new_metrics["boundary_fraction"],
            "delta_boundary_fraction": new_metrics["boundary_fraction"] - base_metrics["boundary_fraction"],
            "base_edge_link_cycle_fraction": base_metrics["edge_link_cycle_fraction"],
            "new_edge_link_cycle_fraction": new_metrics["edge_link_cycle_fraction"],
            "base_harmonic_ratio": base_metrics["harmonic_ratio"],
            "new_harmonic_ratio": new_metrics["harmonic_ratio"],
            "base_exact_residual_ratio": base_metrics["exact_residual_ratio"],
            "new_exact_residual_ratio": new_metrics["exact_residual_ratio"],
            "base_K_mean": base_metrics["K_mean"],
            "new_K_mean": new_metrics["K_mean"],
            "face_a": row.get("face_a"),
            "face_b": row.get("face_b"),
        })
    write_csv(sub / "applied_move_results.csv", applied_rows)
    payload = {
        "growth_rule": growth_rule,
        "base_metrics": base_metrics,
        "candidate_counts": audit["candidate_counts"],
        "best_by_class": audit["best_by_class"],
        "applied_moves": applied_rows,
    }
    (sub / "summary.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
    return payload


def run_apply(args: argparse.Namespace) -> dict:
    out = Path(args.out)
    if out.exists():
        shutil.rmtree(out)
    out.mkdir(parents=True, exist_ok=True)
    growth_rules = args.growth_rules.split(',')
    payloads = [audit_and_apply(gr.strip(), args, out) for gr in growth_rules if gr.strip()]
    comparative_rows: List[dict] = []
    for p in payloads:
        comparative_rows.extend(p["applied_moves"])
    write_csv(out / "comparative_applied_move_results.csv", comparative_rows)
    summary = {
        "script": "test_apply_top_ranked_non_shelling_move_and_reaudit.py",
        "max_level": args.max_level,
        "mode": args.mode,
        "source": args.source,
        "growth_rules": growth_rules,
        "payloads": payloads,
    }
    (out / "comparative_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    return summary


def make_markdown(summary: dict) -> Tuple[str, str, str]:
    real = next((p for p in summary["payloads"] if p["growth_rule"] == "real_growth"), summary["payloads"][0])
    base = real["base_metrics"]
    applied = real["applied_moves"]
    def fmt_rows(rows: List[dict]) -> str:
        lines = ["| selection | move | rank | Δβ1 | Δβ2 | boundary | harmonic | exact_res | status |", "|---|---:|---:|---:|---:|---:|---:|---:|---|"]
        for r in rows:
            lines.append(
                f"| {r['selection_class']}:{r['selection_kind']} | {r['move_class']} | {r['response_rank_legal']} | "
                f"{r['delta_beta1_actual']} | {r['delta_beta2_actual']} | {r['new_boundary_fraction']:.4g} | "
                f"{r['new_harmonic_ratio']:.4g} | {r['new_exact_residual_ratio']:.4g} | {r['apply_status']} |"
            )
        return "\n".join(lines)
    table = fmt_rows(applied)
    results = f"""# RESULTS

## What was tested

This test applies the previously audited top-ranked moves instead of only listing them.
It compares shelling/cap moves with non-shelling handle and quotient candidates, then recomputes Betti numbers and K exactness/harmonic projection.

## Real-growth base complex

```json
{json.dumps({k: base[k] for k in ['vertices','edges','faces','tets','euler','beta0','beta1','beta2','beta3','boundary_fraction','saturated_fraction','edge_link_cycle_fraction','K_mean','exact_residual_ratio','harmonic_ratio']}, indent=2)}
```

## Applied real-growth moves

{table}

## Main interpretation

A shelling/cap move may locally reduce boundary or alter local links, but it is not expected to create a stable global cohomology carrier. A handle/quotient move is the first tested operation that can change the topology class of the carrier.

If a non-shelling move gives delta beta1 or beta2 greater than zero and a nonzero harmonic projection, then the missing ingredient is not the local K-sector but the growth operation class. If it changes Betti numbers but harmonic remains near zero, the topology carrier is opened but the actual K-field still does not land in the harmonic sector. If it does not change Betti numbers, the candidate was only superficially non-shelling.
"""
    summary_md = f"""# SUMMARY

Package: `apply_top_ranked_non_shelling_move_and_reaudit`

Parameters:

```json
{json.dumps({k: summary[k] for k in ['max_level','mode','source','growth_rules']}, indent=2)}
```

Core result for real growth:

- base beta: ({base['beta0']}, {base['beta1']}, {base['beta2']}, {base['beta3']})
- base boundary fraction: {base['boundary_fraction']:.6g}
- base harmonic ratio: {base['harmonic_ratio']:.6g}

The package contains per-growth-rule folders with candidate lists and applied-move reaudit tables, plus a comparative CSV/JSON.
"""
    readme = """# Apply top-ranked non-shelling move and reaudit

Run:

```bash
python3 test_apply_top_ranked_non_shelling_move_and_reaudit.py --max-level 3 --out apply_out
```

Outputs:

- comparative_summary.json
- comparative_applied_move_results.csv
- per-growth-rule move_candidates.csv
- per-growth-rule applied_move_results.csv
- RESULTS.md
- SUMMARY.md

This is a numerical topology/operator diagnostic, not a proof and not yet a CNNA-derived closure law.
"""
    return summary_md, results, readme


def package(out: Path, script_path: Path, zip_path: Path) -> None:
    with zipfile.ZipFile(zip_path, "w", compression=zipfile.ZIP_DEFLATED) as z:
        z.write(script_path, arcname=script_path.name)
        for p in sorted(out.rglob("*")):
            if p.is_file():
                z.write(p, arcname=p.relative_to(out.parent))


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--max-level", type=int, default=3)
    p.add_argument("--mode", choices=["linear", "log", "saturating"], default="linear")
    p.add_argument("--growth-rules", default="real_growth,symmetrized_birth,no_backreaction")
    p.add_argument("--source", choices=["record", "live", "handoff", "aging"], default="live")
    p.add_argument("--transverse-amp", type=float, default=0.42)
    p.add_argument("--max-boundary-faces", type=int, default=90)
    p.add_argument("--max-single-vertices", type=int, default=12)
    p.add_argument("--max-pair-candidates", type=int, default=2500)
    p.add_argument("--max-rows", type=int, default=5000)
    p.add_argument("--out", type=str, default="apply_top_ranked_move_out")
    p.add_argument("--make-zip", action="store_true")
    args = p.parse_args()
    summary = run_apply(args)
    out = Path(args.out)
    summary_md, results_md, readme = make_markdown(summary)
    (out / "SUMMARY.md").write_text(summary_md, encoding="utf-8")
    (out / "RESULTS.md").write_text(results_md, encoding="utf-8")
    (out / "README.md").write_text(readme, encoding="utf-8")
    source_audit = """# SOURCE_AUDIT_1_40

This package continues the non-shelling pairing audit.
It carries forward:

- Script 1/2 dynamic birth/backreaction and transverse sibling offset.
- Script 12 shell-normalized inverse-square kernel as a response-kernel/locality result, not as a shelling/topology result.
- Script 35 K_abc=[A_ab,A_bc] local noncommutative sector.
- Script 40 obstruction of immediate local tetrahedral closure.

It tests whether applying handle/quotient candidate moves changes Betti numbers and K harmonic projection.
"""
    (out / "SOURCE_AUDIT_1_40.md").write_text(source_audit, encoding="utf-8")
    print(json.dumps({"out": str(out), "real_base": summary["payloads"][0]["base_metrics"], "real_applied": summary["payloads"][0]["applied_moves"][:4]}, indent=2))

if __name__ == "__main__":
    main()
