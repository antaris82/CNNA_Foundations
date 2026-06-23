#!/usr/bin/env python3
"""
CNNA / growing real complement network
Real local operator calculus from longitudinal Env/UV planes and transverse DtN handoff defects.

Purpose
-------
This diagnostic deliberately stops before a *-algebra and before J.
It tests whether real deterministic growth generates a composable local operator
calculus:

    local longitudinal role-plane E_p = span_R{Env-response, UV-tail-response}
    L_p : E_p -> E_p              (longitudinal DtN/response operator)
    T_pq : E_p -> E_q             (transverse double-history handoff operator)

The tested obstruction/candidate is the mixed composition defect

    C_pq = L_q T_pq - T_pq L_p.

Nonzero C_pq means the longitudinal and transverse updates are not captured by a
single scalar/radial response.  It is only a real operator-calculus gate, not a
claim that J^2 = -I has been derived.

Growth rule
-----------
The growth is event-resolved and follows the earlier dynamic birth tests:
- ternary children are born sequentially;
- child k senses parent-line + older siblings;
- the newborn immediately backreacts on parent-line and older siblings;
- descendant shell-loads keep aging older local cells.

The DtN layer follows the later dynamic-DtN tests:
- every completed local parent stores immutable record/completion data;
- live DtN data are replayed later from current conductances, ancestor load,
  descendant shell-loads, and directed local response.

Limitations
-----------
Numerical diagnostic only.  No physical i, no J, no *-algebra, no C*-norm, no
AQFT additivity, and no Type-III claim are proved here.
"""

from __future__ import annotations

import argparse
import csv
import itertools
import math
import random
import time
import zipfile
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

import numpy as np

EPS = 1e-12


def mean(xs: Iterable[float]) -> float:
    vals = [float(x) for x in xs if np.isfinite(float(x))]
    return float(np.mean(vals)) if vals else float("nan")


def std(xs: Iterable[float]) -> float:
    vals = [float(x) for x in xs if np.isfinite(float(x))]
    return float(np.std(vals)) if vals else float("nan")


def perc(xs: Iterable[float], q: float) -> float:
    vals = [float(x) for x in xs if np.isfinite(float(x))]
    return float(np.percentile(vals, q)) if vals else float("nan")


def fro(M: np.ndarray) -> float:
    return float(np.linalg.norm(M, ord="fro"))


def norm(v: np.ndarray) -> float:
    return float(np.linalg.norm(v))


def unit(v: np.ndarray) -> np.ndarray:
    n = norm(v)
    if n < EPS:
        return np.zeros_like(v, dtype=float)
    return np.asarray(v, dtype=float) / n


def offdiag_norm(M: np.ndarray) -> float:
    A = np.asarray(M, dtype=float).copy()
    np.fill_diagonal(A, 0.0)
    return fro(A)


def write_csv(path: Path, rows: List[dict]) -> None:
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    keys = sorted({k for r in rows for k in r.keys()})
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)


@dataclass
class Node:
    id: int
    parent: Optional[int]
    level: int
    birth_order: int
    birth_time: int
    birth_g: float
    g: float
    children: List[int] = field(default_factory=list)


class RealGrowth:
    def __init__(
        self,
        *,
        branching: int = 3,
        mode: str = "linear",
        growth_rule: str = "sequential",
        base: float = 1.0,
        alpha_env: float = 0.22,
        ancestor_env_decay: float = 0.55,
        br_ancestor: float = 0.045,
        br_sibling: float = 0.035,
        shell_normalized: bool = True,
    ):
        if branching != 3:
            raise ValueError("this diagnostic currently assumes ternary local cells")
        if mode not in {"linear", "log", "saturating"}:
            raise ValueError("mode must be linear, log, or saturating")
        if growth_rule not in {"sequential", "symmetrized_birth"}:
            raise ValueError("growth_rule must be sequential or symmetrized_birth")
        self.b = branching
        self.mode = mode
        self.growth_rule = growth_rule
        self.base = base
        self.alpha_env = alpha_env
        self.ancestor_env_decay = ancestor_env_decay
        self.br_ancestor = br_ancestor
        self.br_sibling = br_sibling
        self.shell_normalized = shell_normalized

        self.nodes: Dict[int, Node] = {}
        self.local_w: Dict[int, Dict[Tuple[int, int], float]] = defaultdict(lambda: defaultdict(float))
        self.record_completion_matrix: Dict[int, np.ndarray] = {}
        self.live_completion_matrix: Dict[int, np.ndarray] = {}
        self.completion_level: Dict[int, int] = {}
        self.completion_child_g: Dict[int, np.ndarray] = {}
        self.completion_ancestor_env: Dict[int, float] = {}
        self.completion_desc_loads: Dict[int, np.ndarray] = {}
        self.desc_shell_load_by_node: Dict[int, float] = defaultdict(float)

        self.next_id = 0
        self.t = 0
        root = self._new_node(None, 0, 0, 1.0)
        self.root = root.id

    def _new_node(self, parent: Optional[int], level: int, birth_order: int, birth_g: float) -> Node:
        n = Node(self.next_id, parent, level, birth_order, self.t, birth_g, birth_g)
        self.nodes[n.id] = n
        self.next_id += 1
        if parent is not None:
            self.nodes[parent].children.append(n.id)
        return n

    def kernel_value(self, d: int) -> float:
        if not self.shell_normalized:
            return 1.0 / (d * d)
        return 1.0 / ((self.b ** (d - 1)) * d * d)

    def parent_line(self, parent: int) -> List[int]:
        line: List[int] = []
        cur: Optional[int] = parent
        while cur is not None:
            line.append(cur)
            cur = self.nodes[cur].parent
        return line

    def ancestor_env_load(self, parent: int) -> float:
        env = 0.0
        for d, a in enumerate(self.parent_line(parent), start=1):
            env += self.nodes[a].g * (self.ancestor_env_decay ** (d - 1))
        return env

    def birth_env_load(self, parent: int, older: List[int]) -> float:
        if self.growth_rule == "symmetrized_birth":
            return self.ancestor_env_load(parent)
        env = self.ancestor_env_load(parent)
        for s in older:
            env += self.nodes[s].g
        return env

    def child_g_from_env(self, env: float) -> float:
        if self.mode == "linear":
            return self.base + self.alpha_env * env
        if self.mode == "log":
            return self.base + self.alpha_env * math.log1p(env)
        return self.base + self.alpha_env * (env / (1.0 + env))

    def add_child(self, parent: int, order: int) -> int:
        self.t += 1
        older = list(self.nodes[parent].children)
        env = self.birth_env_load(parent, older)
        bg = self.child_g_from_env(env)
        child = self._new_node(parent, self.nodes[parent].level + 1, order, bg)

        if self.growth_rule == "sequential":
            for s in older:
                i = self.nodes[s].birth_order
                j = order
                self.local_w[parent][(i, j)] += self.alpha_env * self.nodes[s].g / (env + EPS) * bg
            for s in older:
                i = order
                j = self.nodes[s].birth_order
                delta = self.br_sibling * bg
                self.nodes[s].g += delta
                self.local_w[parent][(i, j)] += delta
        else:
            if len(self.nodes[parent].children) == self.b:
                ch = self.child_ids_ordered(parent)
                gmean = float(np.mean([self.nodes[c].g for c in ch]))
                w = self.alpha_env * gmean / (self.ancestor_env_load(parent) + EPS) * gmean
                for i in range(1, self.b + 1):
                    for j in range(1, self.b + 1):
                        if i != j:
                            self.local_w[parent][(i, j)] = w

        for d, a in enumerate(self.parent_line(parent), start=1):
            kd = self.kernel_value(d)
            self.nodes[a].g += self.br_ancestor * bg * kd
            self.desc_shell_load_by_node[a] += bg * kd

        if len(self.nodes[parent].children) == self.b:
            self._snapshot_completion(parent)
        return child.id

    def grow_one_level(self, frontier: List[int]) -> List[int]:
        new_frontier: List[int] = []
        for p in frontier:
            for k in range(1, self.b + 1):
                new_frontier.append(self.add_child(p, k))
        return new_frontier

    def grow(self, max_level: int) -> None:
        frontier = [self.root]
        for _ in range(1, max_level + 1):
            frontier = self.grow_one_level(frontier)

    def completed_parent_ids(self) -> List[int]:
        return [n.id for n in self.nodes.values() if len(n.children) == self.b]

    def child_ids_ordered(self, parent: int) -> List[int]:
        ch = list(self.nodes[parent].children)
        ch.sort(key=lambda c: self.nodes[c].birth_order)
        return ch

    def child_g_vector(self, parent: int) -> np.ndarray:
        return np.array([self.nodes[c].g for c in self.child_ids_ordered(parent)], dtype=float)

    def child_descendant_loads(self, parent: int) -> np.ndarray:
        return np.array([self.desc_shell_load_by_node.get(c, 0.0) for c in self.child_ids_ordered(parent)], dtype=float)

    def local_record_matrix(self, parent: int) -> Optional[np.ndarray]:
        if len(self.nodes[parent].children) != self.b:
            return None
        M = np.zeros((3, 3), dtype=float)
        w = self.local_w[parent]
        for i in range(1, 4):
            for j in range(1, 4):
                if i != j:
                    M[j - 1, i - 1] = w.get((i, j), 0.0)
        return M

    def live_replay_matrix(self, parent: int) -> Optional[np.ndarray]:
        if len(self.nodes[parent].children) != self.b:
            return None
        g = self.child_g_vector(parent)
        ancestor_env = self.ancestor_env_load(parent)
        M = np.zeros((3, 3), dtype=float)
        if self.growth_rule == "symmetrized_birth":
            w = self.alpha_env * float(np.mean(g)) / (ancestor_env + EPS) * float(np.mean(g))
            for i in range(3):
                for j in range(3):
                    if i != j:
                        M[j, i] = w
            return M
        for i in range(1, 4):
            for j in range(1, 4):
                if i == j:
                    continue
                source = g[i - 1]
                target = g[j - 1]
                if i < j:
                    prefix = ancestor_env + float(np.sum(g[: j - 1]))
                    M[j - 1, i - 1] = self.alpha_env * source / (prefix + EPS) * target
                else:
                    M[j - 1, i - 1] = self.br_sibling * source
        return M

    def _snapshot_completion(self, parent: int) -> None:
        if parent in self.completion_level:
            return
        R = self.local_record_matrix(parent)
        L = self.live_replay_matrix(parent)
        if R is None or L is None:
            return
        self.record_completion_matrix[parent] = R.copy()
        self.live_completion_matrix[parent] = L.copy()
        self.completion_level[parent] = self.nodes[parent].level + 1
        self.completion_child_g[parent] = self.child_g_vector(parent).copy()
        self.completion_ancestor_env[parent] = self.ancestor_env_load(parent)
        self.completion_desc_loads[parent] = self.child_descendant_loads(parent).copy()


def sym_response(M: np.ndarray) -> np.ndarray:
    C = 0.5 * (M + M.T)
    C = np.maximum(C, 0.0)
    np.fill_diagonal(C, 0.0)
    return C


def dtn_decomposition(
    child_g: np.ndarray,
    ancestor_env: float,
    directed_response: np.ndarray,
    descendant_loads: np.ndarray,
    *,
    eta_sibling: float = 0.35,
    eta_descendant_shunt: float = 0.08,
    eta_env_shunt: float = 0.18,
) -> dict:
    g = np.maximum(np.asarray(child_g, dtype=float), EPS)
    desc = np.maximum(np.asarray(descendant_loads, dtype=float), 0.0)
    C = eta_sibling * sym_response(directed_response)
    port_shunt = eta_descendant_shunt * desc
    core_shunt = max(eta_env_shunt * float(ancestor_env), EPS)

    D_env = np.diag(g) - np.outer(g, g) / (float(np.sum(g)) + core_shunt)
    D_sib = np.zeros((3, 3), dtype=float)
    for i in range(3):
        for j in range(i + 1, 3):
            c = max(float(C[i, j]), 0.0)
            D_sib[i, i] += c
            D_sib[j, j] += c
            D_sib[i, j] -= c
            D_sib[j, i] -= c
    D_uv = np.diag(port_shunt)
    D_total = D_env + D_sib + D_uv
    uv_port_mass = g + port_shunt
    env_row = np.sum(D_env, axis=1)
    return {
        "D_total": 0.5 * (D_total + D_total.T),
        "D_env": 0.5 * (D_env + D_env.T),
        "D_sib": 0.5 * (D_sib + D_sib.T),
        "D_uv": D_uv,
        "env_vec": env_row,
        "uv_vec": uv_port_mass,
        "core_shunt": core_shunt,
        "port_shunt": port_shunt,
    }


def address(model: RealGrowth, node_id: int) -> Tuple[int, ...]:
    out: List[int] = []
    cur = node_id
    while model.nodes[cur].parent is not None:
        out.append(int(model.nodes[cur].birth_order))
        cur = model.nodes[cur].parent
    return tuple(reversed(out))


class LocalCell:
    def __init__(self, model: RealGrowth, parent: int, max_level: int):
        self.model = model
        self.parent = int(parent)
        self.level = int(model.nodes[parent].level)
        self.completion_level = int(model.completion_level[parent])
        self.age = int(max_level - self.completion_level)
        self.address = address(model, parent)
        self.root_sector = int(self.address[0]) if self.address else 0
        self.suffix = tuple(self.address[1:]) if self.address else tuple()

        R = model.record_completion_matrix[parent]
        L0 = model.live_completion_matrix[parent]
        L1 = model.live_replay_matrix(parent)
        if L1 is None:
            raise ValueError("missing live replay")
        self.R_dir = R.copy()
        self.L0_dir = L0.copy()
        self.L1_dir = L1.copy()

        self.record = dtn_decomposition(
            model.completion_child_g[parent],
            model.completion_ancestor_env[parent],
            self.R_dir,
            model.completion_desc_loads[parent],
        )
        self.live0 = dtn_decomposition(
            model.completion_child_g[parent],
            model.completion_ancestor_env[parent],
            self.L0_dir,
            model.completion_desc_loads[parent],
        )
        self.live1 = dtn_decomposition(
            model.child_g_vector(parent),
            model.ancestor_env_load(parent),
            self.L1_dir,
            model.child_descendant_loads(parent),
        )
        self.aging = self.live1["D_total"] - self.live0["D_total"]
        self.handoff = self.live1["D_total"] - self.record["D_total"]
        self.L_record = role_longitudinal_operator(self.record)
        self.L_live = role_longitudinal_operator(self.live1)
        self.L_aging = role_aging_operator(self.live0, self.live1)

    def port_vector(self, port: int, source_mode: str) -> np.ndarray:
        k = port - 1
        if source_mode == "record":
            return self.record["D_total"][k, :].copy()
        if source_mode == "live":
            return self.live1["D_total"][k, :].copy()
        if source_mode == "aging":
            return self.aging[k, :].copy()
        if source_mode == "handoff":
            return self.handoff[k, :].copy()
        raise ValueError(f"unknown source_mode {source_mode}")

    def longitudinal_operator(self, source_mode: str) -> np.ndarray:
        if source_mode == "record":
            return self.L_record.copy()
        if source_mode == "aging":
            return self.L_aging.copy()
        return self.L_live.copy()


def role_longitudinal_operator(dec: dict) -> np.ndarray:
    env_vec = np.asarray(dec["env_vec"], dtype=float)
    uv_vec = np.asarray(dec["uv_vec"], dtype=float)
    D_env = np.asarray(dec["D_env"], dtype=float)
    D_uv_role = np.diag(np.maximum(uv_vec, 0.0))
    D_cross = np.asarray(dec["D_total"], dtype=float) - D_env - np.asarray(dec["D_uv"], dtype=float)
    a = norm(env_vec)
    d = norm(uv_vec)
    b = offdiag_norm(D_cross) + abs(float(np.dot(unit(env_vec), unit(uv_vec)))) * math.sqrt(a * d + EPS)
    return np.array([[a, b], [b, d]], dtype=float)


def role_aging_operator(dec0: dict, dec1: dict) -> np.ndarray:
    return role_longitudinal_operator(dec1) - role_longitudinal_operator(dec0)


def transverse_operator(a: LocalCell, b: LocalCell, port: int, source_mode: str) -> Tuple[np.ndarray, dict]:
    delta = b.port_vector(port, source_mode) - a.port_vector(port, source_mode)
    env_ref = unit(a.live1["env_vec"] + b.live1["env_vec"])
    uv_ref = unit(a.live1["uv_vec"] + b.live1["uv_vec"])
    e = float(np.dot(delta, env_ref))
    u = float(np.dot(delta, uv_ref))
    proj = e * env_ref + u * uv_ref
    residual = delta - proj
    T = np.array([[0.0, u], [e, 0.0]], dtype=float)
    return T, {
        "delta_norm": norm(delta),
        "delta_env_component": e,
        "delta_uv_component": u,
        "delta_residual_norm": norm(residual),
        "T_norm": fro(T),
        "T_skew_norm": fro(0.5 * (T - T.T)),
        "T_sym_norm": fro(0.5 * (T + T.T)),
    }


def pair_row(a: LocalCell, b: LocalCell, port: int, source_mode: str, control: str, pair_index: int) -> dict:
    T, tstats = transverse_operator(a, b, port, source_mode)
    La = a.longitudinal_operator(source_mode)
    Lb = b.longitudinal_operator(source_mode)
    left = Lb @ T
    right = T @ La
    C = left - right
    denom = fro(left) + fro(right) + EPS
    row = {
        "control": control,
        "source_mode": source_mode,
        "pair_index": pair_index,
        "parent_a": a.parent,
        "parent_b": b.parent,
        "port": int(port),
        "level": a.level,
        "age_a": a.age,
        "age_b": b.age,
        "root_sector_a": a.root_sector,
        "root_sector_b": b.root_sector,
        "address_a": ".".join(map(str, a.address)),
        "address_b": ".".join(map(str, b.address)),
        "suffix": ".".join(map(str, a.suffix)),
        "L_a_norm": fro(La),
        "L_b_norm": fro(Lb),
        "compose_left_norm": fro(left),
        "compose_right_norm": fro(right),
        "commutator_norm": fro(C),
        "commutator_rel": fro(C) / denom,
        "commutator_skew_norm": fro(0.5 * (C - C.T)),
        "commutator_sym_norm": fro(0.5 * (C + C.T)),
        "plane_drift_L_norm": fro(Lb - La),
    }
    row.update(tstats)
    return row


def build_cells(model: RealGrowth, max_level: int) -> Dict[int, LocalCell]:
    cells: Dict[int, LocalCell] = {}
    for p in model.completed_parent_ids():
        if p in model.completion_level:
            cells[p] = LocalCell(model, p, max_level)
    return cells


def double_history_pairs(cells: Dict[int, LocalCell]) -> List[Tuple[LocalCell, LocalCell, int]]:
    groups: Dict[Tuple[int, Tuple[int, ...], int], List[LocalCell]] = defaultdict(list)
    for c in cells.values():
        if c.level < 1:
            continue
        for port in (1, 2, 3):
            groups[(c.level, c.suffix, port)].append(c)
    pairs = []
    for _key, gs in sorted(groups.items(), key=lambda kv: (kv[0][0], kv[0][2], kv[0][1])):
        roots = {g.root_sector for g in gs}
        if len(gs) < 2 or len(roots) < 2:
            continue
        port = int(_key[2])
        for a, b in itertools.combinations(gs, 2):
            if a.root_sector != b.root_sector:
                pairs.append((a, b, port))
    return pairs


def identical_pairs(cells: Dict[int, LocalCell]) -> List[Tuple[LocalCell, LocalCell, int]]:
    out = []
    for c in cells.values():
        if c.level < 1:
            continue
        for port in (1, 2, 3):
            out.append((c, c, port))
    return out


def random_pairs(cells: Dict[int, LocalCell], n: int, seed: int) -> List[Tuple[LocalCell, LocalCell, int]]:
    rng = random.Random(seed)
    buckets: Dict[Tuple[int, int], List[LocalCell]] = defaultdict(list)
    for c in cells.values():
        if c.level >= 1:
            for port in (1, 2, 3):
                buckets[(c.level, port)].append(c)
    keys = [k for k, v in buckets.items() if len(v) >= 2]
    out = []
    attempts = 0
    while len(out) < n and attempts < max(2000, 40 * n):
        attempts += 1
        level, port = rng.choice(keys)
        a, b = rng.sample(buckets[(level, port)], 2)
        if a.suffix == b.suffix:
            continue
        out.append((a, b, port))
    return out


def summarize(rows: List[dict], label: str) -> dict:
    out = {"label": label, "count": len(rows)}
    keys = [
        "T_norm", "delta_norm", "delta_residual_norm", "commutator_norm",
        "commutator_rel", "commutator_skew_norm", "commutator_sym_norm",
        "plane_drift_L_norm",
    ]
    for k in keys:
        vals = [r[k] for r in rows if k in r]
        out[f"{k}_mean"] = mean(vals)
        out[f"{k}_std"] = std(vals)
        out[f"{k}_p50"] = perc(vals, 50)
        out[f"{k}_p95"] = perc(vals, 95)
        out[f"{k}_max"] = max([float(v) for v in vals], default=float("nan"))
    out["fraction_nonzero_T"] = mean([float(r["T_norm"] > 1e-10) for r in rows]) if rows else 0.0
    out["fraction_nonzero_commutator"] = mean([float(r["commutator_norm"] > 1e-10) for r in rows]) if rows else 0.0
    return out


def by_level(rows: List[dict], label: str) -> List[dict]:
    buckets: Dict[Tuple[int, str], List[dict]] = defaultdict(list)
    for r in rows:
        buckets[(int(r["level"]), str(r["source_mode"]))].append(r)
    out = []
    for (level, mode), rs in sorted(buckets.items()):
        s = summarize(rs, f"{label}_L{level}_{mode}")
        s["level"] = level
        s["source_mode"] = mode
        out.append(s)
    return out


def run_case(case_name: str, model: RealGrowth, max_level: int, source_modes: List[str], random_n: int, seed: int) -> Tuple[List[dict], List[dict]]:
    cells = build_cells(model, max_level)
    dh = double_history_pairs(cells)
    ident = identical_pairs(cells)
    rand = random_pairs(cells, min(random_n, max(1, len(dh))), seed)
    all_rows: List[dict] = []
    summaries: List[dict] = []
    pair_sets = [
        ("double_history_suffix_quotient", dh),
        ("identical_history_control", ident),
        ("random_same_level_port_baseline", rand),
    ]
    for source_mode in source_modes:
        for control, pairs in pair_sets:
            rows = [pair_row(a, b, port, source_mode, f"{case_name}:{control}", i) for i, (a, b, port) in enumerate(pairs)]
            all_rows.extend(rows)
            summaries.append(summarize(rows, f"{case_name}:{control}:{source_mode}"))
            summaries.extend(by_level(rows, f"{case_name}:{control}"))
    meta = {
        "label": f"{case_name}:meta",
        "count": len(cells),
        "completed_cells": len(cells),
        "double_history_pairs": len(dh),
        "identical_pairs": len(ident),
        "random_pairs": len(rand),
        "nodes": len(model.nodes),
    }
    summaries.insert(0, meta)
    return all_rows, summaries


def write_results_md(path: Path, args: argparse.Namespace, summaries: List[dict]) -> None:
    def find(label: str) -> Optional[dict]:
        for s in summaries:
            if s.get("label") == label:
                return s
        return None

    primary = find("real_growth:double_history_suffix_quotient:handoff")
    ident = find("real_growth:identical_history_control:handoff")
    nobr = find("no_backreaction:double_history_suffix_quotient:handoff")
    sym = find("symmetrized_birth:double_history_suffix_quotient:handoff")
    rand = find("real_growth:random_same_level_port_baseline:handoff")

    lines = []
    lines.append("# RESULTS: real local operator calculus, longitudinal/transverse")
    lines.append("")
    lines.append("## Status")
    lines.append("")
    lines.append("This is a real operator-calculus diagnostic. It does not derive J and does not test J²=-I.")
    lines.append("It checks whether true event-resolved growth produces composable longitudinal and transverse operators.")
    lines.append("")
    lines.append("## Parameters")
    lines.append("")
    lines.append(f"- max_level: `{args.max_level}`")
    lines.append(f"- conductance mode: `{args.mode}`")
    lines.append(f"- source modes: `{', '.join(args.source_modes)}`")
    lines.append(f"- random seed: `{args.seed}`")
    lines.append("")
    lines.append("## Primary handoff gate")
    lines.append("")
    if primary:
        lines.append(f"- double-history handoff pairs: `{int(primary['count'])}`")
        lines.append(f"- mean transverse T norm: `{primary['T_norm_mean']:.12g}`")
        lines.append(f"- mean composition defect ||L_q T - T L_p||: `{primary['commutator_norm_mean']:.12g}`")
        lines.append(f"- mean relative composition defect: `{primary['commutator_rel_mean']:.12g}`")
        lines.append(f"- nonzero T fraction: `{primary['fraction_nonzero_T']:.6g}`")
        lines.append(f"- nonzero composition-defect fraction: `{primary['fraction_nonzero_commutator']:.6g}`")
    lines.append("")
    lines.append("## Controls")
    lines.append("")
    for name, s in [
        ("identical-history", ident),
        ("real-growth random same-level/port baseline", rand),
        ("no-backreaction double-history", nobr),
        ("symmetrized-birth double-history", sym),
    ]:
        if not s:
            continue
        lines.append(f"### {name}")
        lines.append(f"- rows: `{int(s['count'])}`")
        lines.append(f"- mean T norm: `{s['T_norm_mean']:.12g}`")
        lines.append(f"- mean composition defect: `{s['commutator_norm_mean']:.12g}`")
        lines.append(f"- mean relative composition defect: `{s['commutator_rel_mean']:.12g}`")
        lines.append("")
    lines.append("## Interpretation")
    lines.append("")
    lines.append("Positive result criterion for this gate: double-history handoff rows have nonzero transverse T and nonzero mixed composition defect, while identical-history controls vanish. This would mean that Script 23's DtN defect can be promoted to a real longitudinal/transverse operator-calculus candidate.")
    lines.append("")
    lines.append("Negative result criterion: if T or the composition defect vanishes under real growth, then Script 23 remains a scalar/vector handoff obstruction, not yet an operator-calculus object.")
    lines.append("")
    lines.append("Next gate after a positive result: define a canonical pairing/energy form and test whether the generated real operator system admits a stable adjoint operation A -> A*.")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--max-level", type=int, default=6)
    ap.add_argument("--mode", choices=["linear", "log", "saturating"], default="linear")
    ap.add_argument("--out", type=Path, default=Path("real_local_operator_calculus_out_L6"))
    ap.add_argument("--seed", type=int, default=20260621)
    ap.add_argument("--random-n", type=int, default=3000)
    ap.add_argument("--source-modes", nargs="+", default=["handoff", "aging", "live", "record"], choices=["handoff", "aging", "live", "record"])
    args = ap.parse_args()

    t0 = time.time()
    args.out.mkdir(parents=True, exist_ok=True)

    cases = []
    real = RealGrowth(mode=args.mode, growth_rule="sequential")
    real.grow(args.max_level)
    cases.append(("real_growth", real))

    nobr = RealGrowth(mode=args.mode, growth_rule="sequential", br_ancestor=0.0, br_sibling=0.0)
    nobr.grow(args.max_level)
    cases.append(("no_backreaction", nobr))

    sym = RealGrowth(mode=args.mode, growth_rule="symmetrized_birth")
    sym.grow(args.max_level)
    cases.append(("symmetrized_birth", sym))

    all_rows: List[dict] = []
    all_summaries: List[dict] = []
    for idx, (name, model) in enumerate(cases):
        rows, summaries = run_case(name, model, args.max_level, args.source_modes, args.random_n, args.seed + idx)
        all_rows.extend(rows)
        all_summaries.extend(summaries)
        write_csv(args.out / f"pair_table_{name}.csv", rows)
        write_csv(args.out / f"summary_table_{name}.csv", summaries)

    write_csv(args.out / "pair_table_all.csv", all_rows)
    write_csv(args.out / "summary_table_all.csv", all_summaries)
    write_results_md(args.out / "RESULTS_real_local_operator_calculus.md", args, all_summaries)

    summary = [
        "CNNA real local operator calculus diagnostic",
        f"max_level={args.max_level}",
        f"mode={args.mode}",
        f"rows={len(all_rows)} summaries={len(all_summaries)}",
        f"elapsed_sec={time.time() - t0:.3f}",
        "Primary file: RESULTS_real_local_operator_calculus.md",
    ]
    (args.out / "SUMMARY.txt").write_text("\n".join(summary) + "\n", encoding="utf-8")

    zip_path = args.out.with_suffix(".zip")
    with zipfile.ZipFile(zip_path, "w", compression=zipfile.ZIP_DEFLATED) as z:
        z.write(Path(__file__), Path(__file__).name)
        for p in args.out.rglob("*"):
            if p.is_file():
                z.write(p, p.relative_to(args.out.parent))
    print("wrote", args.out)
    print("wrote", zip_path)


if __name__ == "__main__":
    main()
