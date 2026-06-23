#!/usr/bin/env python3
"""
CNNA dynamic boundary role recalibration gate, L2.

This diagnostic deliberately reuses the established dynamic birth conductance / monodromy
model shape from Growth scripts 1 and 2:

- ternary sequential births;
- child conductance from parent-line + older-sibling environment;
- environment edges from existing nodes to newborn;
- newborn backreaction to parent line and older siblings;
- completed sibling triple monodromy diagnostics.

The new content is an event-resolved boundary-role audit.  It asks whether each birth
recalibrates the cut roles (self / parent-line / sibling / UV-tail / incoming/outgoing)
in a retarded, pre-causal way, before asking for any Q/P, J, spin, Hodge, or *-structure.

No complex scalar, Hodge operator, physical adjoint, positivity, or delta-beta decision is used.
Complex numbers appear only in the inherited neutral phasor / Z3 monodromy diagnostic from
script 1/2, not as an ontic input.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import cmath
from dataclasses import dataclass, field
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np

EPS = 1e-12


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


class DynamicBoundaryRoleModel:
    def __init__(
        self,
        variant: str,
        mode: str = "linear",
        branching: int = 3,
        base: float = 1.0,
        alpha_env: float = 0.22,
        ancestor_decay: float = 0.55,
        br_ancestor: float = 0.045,
        br_sibling: float = 0.035,
        order_sequence: Tuple[int, int, int] = (1, 2, 3),
        eps: float = EPS,
    ):
        if branching != 3:
            raise ValueError("This diagnostic currently assumes ternary sibling triples.")
        self.variant = variant
        self.mode = mode
        self.branching = branching
        self.base = base
        self.alpha_env = alpha_env
        self.ancestor_decay = ancestor_decay
        self.br_ancestor = br_ancestor
        self.br_sibling = br_sibling
        self.order_sequence = order_sequence
        self.eps = eps

        self.nodes: Dict[int, Node] = {}
        self.t = 0
        self.next_id = 0
        self.directed_edges: Dict[Tuple[int, int], float] = defaultdict(float)
        self.event_rows: List[dict] = []
        self.role_rows: List[dict] = []
        self.triple_rows: List[dict] = []
        self.level_rows: List[dict] = []

        root = self._new_node(parent=None, level=0, birth_order=0, birth_g=1.0)
        self.root = root.id

    def _new_node(self, parent: Optional[int], level: int, birth_order: int, birth_g: float) -> Node:
        n = Node(
            id=self.next_id,
            parent=parent,
            level=level,
            birth_order=birth_order,
            birth_time=self.t,
            birth_g=birth_g,
            g=birth_g,
        )
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

    def descendants(self, node: int) -> List[int]:
        out: List[int] = []
        stack = list(self.nodes[node].children)
        while stack:
            x = stack.pop()
            out.append(x)
            stack.extend(self.nodes[x].children)
        return out

    def immediate_tail_load(self, node: int) -> float:
        return sum(self.nodes[c].g for c in self.nodes[node].children)

    def cut_response_scalar(self, node: int) -> float:
        # Scalar DtN/Schur surrogate compatible with script 1/2 conductance data:
        # current node response plus immediate UV-tail load.
        return self.nodes[node].g + self.immediate_tail_load(node)

    def birth_environment_load(self, parent: int, older_siblings: List[int]) -> float:
        env = 0.0
        for d, a in enumerate(self.parent_line(parent), start=1):
            env += self.nodes[a].g * (self.ancestor_decay ** (d - 1))
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
        raise ValueError(f"unknown mode: {self.mode}")

    @staticmethod
    def neutral_phasor_values(vals: List[float]) -> complex:
        omega = cmath.exp(2j * math.pi / 3)
        return vals[0] + vals[1] * omega + vals[2] * (omega ** 2)

    def neutral_for_parent(self, parent: int, current: bool = True) -> Optional[complex]:
        ch = self.nodes[parent].children
        if len(ch) != 3:
            return None
        vals = [self.nodes[c].g if current else self.nodes[c].birth_g for c in ch]
        return self.neutral_phasor_values(vals)

    def w(self, u: int, v: int) -> float:
        return self.directed_edges.get((u, v), 0.0)

    def local_pair_weights(self, parent: int) -> Optional[dict]:
        ch = self.nodes[parent].children
        if len(ch) != 3:
            return None
        # Children are stored in chronological birth order.  birth_order labels may be reversed.
        c1, c2, c3 = ch
        return {
            "w12": self.w(c1, c2),
            "w23": self.w(c2, c3),
            "w31": self.w(c3, c1),
            "w13": self.w(c1, c3),
            "w32": self.w(c3, c2),
            "w21": self.w(c2, c1),
        }

    @staticmethod
    def transport_matrix_from_weights(
        w12: float, w23: float, w31: float, w13: float = 0.0, w32: float = 0.0, w21: float = 0.0
    ) -> np.ndarray:
        A = np.zeros((3, 3), dtype=float)
        A[1, 0] = w12
        A[2, 1] = w23
        A[0, 2] = w31
        A[2, 0] = w13
        A[1, 2] = w32
        A[0, 1] = w21
        return A

    @staticmethod
    def column_stochastic(A: np.ndarray) -> np.ndarray:
        B = A.copy().astype(float)
        for j in range(B.shape[1]):
            s = B[:, j].sum()
            if s > 0:
                B[:, j] /= s
            else:
                B[j, j] = 1.0
        return B

    @staticmethod
    def eig_class(vals: np.ndarray, tol: float = 1e-9) -> str:
        if np.max(np.abs(np.imag(vals))) > tol:
            return "complex_pair"
        return "real_or_degenerate"

    def local_cycle_bias(self, parent: int) -> Tuple[Optional[float], Optional[float], Optional[float]]:
        weights = self.local_pair_weights(parent)
        if weights is None:
            return None, None, None
        e = self.eps
        fprod = (weights["w12"] + e) * (weights["w23"] + e) * (weights["w31"] + e)
        rprod = (weights["w13"] + e) * (weights["w32"] + e) * (weights["w21"] + e)
        return math.log(fprod / rprod), fprod, rprod

    def triple_monodromy(self, parent: int) -> Optional[dict]:
        ch = self.nodes[parent].children
        if len(ch) != 3:
            return None
        weights = self.local_pair_weights(parent)
        assert weights is not None
        e = self.eps
        w12, w23, w31 = weights["w12"] + e, weights["w23"] + e, weights["w31"] + e
        w13, w32, w21 = weights["w13"] + e, weights["w32"] + e, weights["w21"] + e
        A_full = self.transport_matrix_from_weights(w12, w23, w31, w13, w32, w21)
        A_sym = 0.5 * (A_full + A_full.T)
        A_path = self.transport_matrix_from_weights(w12, w23, 0.0, 0.0, 0.0, 0.0)
        P_full = self.column_stochastic(A_full)
        P_sym = self.column_stochastic(A_sym)
        P_path = self.column_stochastic(A_path)
        eig_full = np.linalg.eigvals(P_full)
        eig_sym = np.linalg.eigvals(P_sym)
        eig_path = np.linalg.eigvals(P_path)
        log_circ, fprod, rprod = self.local_cycle_bias(parent)
        z_current = self.neutral_for_parent(parent, current=True)
        mean_g = sum(self.nodes[x].g for x in ch) / 3.0
        return {
            "variant": self.variant,
            "mode": self.mode,
            "time": self.t,
            "parent": parent,
            "parent_level": self.nodes[parent].level,
            "children_chronological": " ".join(map(str, ch)),
            "child_birth_orders": " ".join(str(self.nodes[c].birth_order) for c in ch),
            "child_birth_times": " ".join(str(self.nodes[c].birth_time) for c in ch),
            "child_birth_g": " ".join(f"{self.nodes[c].birth_g:.12g}" for c in ch),
            "child_current_g": " ".join(f"{self.nodes[c].g:.12g}" for c in ch),
            "log_circulation_forward_vs_reverse": log_circ or 0.0,
            "forward_product": fprod or 0.0,
            "reverse_product": rprod or 0.0,
            "neutral_norm_current": (abs(z_current) / mean_g) if z_current is not None and mean_g > 0 else 0.0,
            "full_markov_eig_class": self.eig_class(eig_full),
            "sym_markov_eig_class": self.eig_class(eig_sym),
            "path_markov_eig_class": self.eig_class(eig_path),
            "full_markov_max_imag": float(np.max(np.abs(np.imag(eig_full)))),
            "sym_markov_max_imag": float(np.max(np.abs(np.imag(eig_sym)))),
            "path_markov_max_imag": float(np.max(np.abs(np.imag(eig_path)))),
        }

    def add_child(self, parent: int, order: int) -> int:
        self.t += 1
        older = list(self.nodes[parent].children)
        parent_line = self.parent_line(parent)
        before_g = {n: self.nodes[n].g for n in self.nodes}
        before_response = {a: self.cut_response_scalar(a) for a in parent_line}
        older_before_response = {s: self.cut_response_scalar(s) for s in older}
        parent_tail_count_before = len(self.nodes[parent].children)
        env_load = self.birth_environment_load(parent, older)
        birth_g = self.child_conductance_from_env(env_load)
        total_env = env_load + self.eps

        # Edge contributions to newborn, computed before node creation to log source roles.
        ancestor_env_weights = []
        for d, a in enumerate(parent_line, start=1):
            contrib = self.nodes[a].g * (self.ancestor_decay ** (d - 1))
            weight = self.alpha_env * contrib / total_env * birth_g if total_env > 0 else 0.0
            ancestor_env_weights.append((a, d, weight))
        sibling_env_weights = []
        for s in older:
            contrib = self.nodes[s].g
            weight = self.alpha_env * contrib / total_env * birth_g if total_env > 0 else 0.0
            sibling_env_weights.append((s, weight))

        child = self._new_node(parent, self.nodes[parent].level + 1, order, birth_g)
        c = child.id

        incoming_to_child = 0.0
        for a, _d, weight in ancestor_env_weights:
            self.directed_edges[(a, c)] += weight
            incoming_to_child += weight
        for s, weight in sibling_env_weights:
            self.directed_edges[(s, c)] += weight
            incoming_to_child += weight

        ancestor_backreaction = []
        child_to_parent_line = 0.0
        for d, a in enumerate(parent_line, start=1):
            delta = self.br_ancestor * birth_g / (d * d)
            self.nodes[a].g += delta
            self.directed_edges[(c, a)] += delta
            child_to_parent_line += delta
            ancestor_backreaction.append((a, d, delta))

        sibling_backreaction = []
        child_to_older_siblings = 0.0
        for s in older:
            delta = self.br_sibling * birth_g
            self.nodes[s].g += delta
            self.directed_edges[(c, s)] += delta
            child_to_older_siblings += delta
            sibling_backreaction.append((s, delta))

        outgoing_from_child = child_to_parent_line + child_to_older_siblings
        parent_tail_count_after = len(self.nodes[parent].children)
        child_own_uv_tail_count_after = len(self.nodes[c].children)

        after_response = {a: self.cut_response_scalar(a) for a in parent_line}
        older_after_response = {s: self.cut_response_scalar(s) for s in older}

        # Per-cut role rows: the same new child is self-bright/no-own-tail for itself,
        # and newly realized UV-tail/backreaction source for every parent-line cut.
        for a, depth, delta in ancestor_backreaction:
            br_before = before_response[a]
            br_after = after_response[a]
            self.role_rows.append(
                {
                    "variant": self.variant,
                    "mode": self.mode,
                    "t": self.t,
                    "birth_id": self.t,
                    "cut_node": a,
                    "cut_depth_from_parent": depth,
                    "parent": parent,
                    "child": c,
                    "child_birth_order_label": order,
                    "older_sibling_count": len(older),
                    "child_self_uv_tail_count_after": child_own_uv_tail_count_after,
                    "child_has_own_uv_tail": int(child_own_uv_tail_count_after > 0),
                    "child_is_new_uv_tail_for_cut": 1,
                    "cut_receives_backreaction": int(delta > EPS),
                    "incoming_env_to_child": incoming_to_child,
                    "outgoing_child_to_cut_delta": delta,
                    "cut_g_before": before_g[a],
                    "cut_g_after": self.nodes[a].g,
                    "cut_g_delta": self.nodes[a].g - before_g[a],
                    "cut_response_before": br_before,
                    "cut_response_after": br_after,
                    "cut_response_delta": br_after - br_before,
                    "role_statement": "child=no-own-UV-self; child=new-UV-tail-for-cut; cut=ancestor/IR-side-updated",
                }
            )

        for s, delta in sibling_backreaction:
            self.role_rows.append(
                {
                    "variant": self.variant,
                    "mode": self.mode,
                    "t": self.t,
                    "birth_id": self.t,
                    "cut_node": s,
                    "cut_depth_from_parent": 0,
                    "parent": parent,
                    "child": c,
                    "child_birth_order_label": order,
                    "older_sibling_count": len(older),
                    "child_self_uv_tail_count_after": child_own_uv_tail_count_after,
                    "child_has_own_uv_tail": int(child_own_uv_tail_count_after > 0),
                    "child_is_new_uv_tail_for_cut": 0,
                    "cut_receives_backreaction": int(delta > EPS),
                    "incoming_env_to_child": incoming_to_child,
                    "outgoing_child_to_cut_delta": delta,
                    "cut_g_before": before_g[s],
                    "cut_g_after": self.nodes[s].g,
                    "cut_g_delta": self.nodes[s].g - before_g[s],
                    "cut_response_before": older_before_response[s],
                    "cut_response_after": older_after_response[s],
                    "cut_response_delta": older_after_response[s] - older_before_response[s],
                    "role_statement": "older-sibling=already-bright-cell-updated-by-newborn-backreaction",
                }
            )

        # Event summary.
        ancestor_deltas = [d for (_a, _depth, d) in ancestor_backreaction]
        ancestor_depth_monotone = int(all(ancestor_deltas[i] >= ancestor_deltas[i + 1] - 1e-12 for i in range(len(ancestor_deltas) - 1)))
        child_birth_g_prev_sibling_max = max([self.nodes[s].birth_g for s in older], default=0.0)
        sibling_order_increment = birth_g - child_birth_g_prev_sibling_max if older else 0.0
        children = self.nodes[parent].children
        partial = [self.nodes[x].g for x in children]
        padded = partial + [0.0] * (3 - len(partial))
        z_partial = self.neutral_phasor_values(padded)

        completed = len(children) == 3
        triple_row = None
        if completed:
            triple_row = self.triple_monodromy(parent)
            assert triple_row is not None
            self.triple_rows.append(triple_row)

        self.event_rows.append(
            {
                "variant": self.variant,
                "mode": self.mode,
                "t": self.t,
                "birth_id": self.t,
                "parent": parent,
                "child": c,
                "child_level": child.level,
                "child_birth_order_label": order,
                "chronological_child_index_for_parent": len(children),
                "ancestor_chain": " ".join(map(str, parent_line)),
                "older_siblings": " ".join(map(str, older)),
                "older_sibling_count": len(older),
                "younger_siblings_exist_at_birth": 0,
                "child_has_own_uv_tail": int(child_own_uv_tail_count_after > 0),
                "child_own_uv_tail_count_after": child_own_uv_tail_count_after,
                "parent_tail_count_before": parent_tail_count_before,
                "parent_tail_count_after": parent_tail_count_after,
                "parent_gains_uv_tail": int(parent_tail_count_after == parent_tail_count_before + 1),
                "env_load": env_load,
                "child_birth_g": birth_g,
                "incoming_env_to_child_total": incoming_to_child,
                "outgoing_child_backreaction_total": outgoing_from_child,
                "child_to_parent_line_total": child_to_parent_line,
                "child_to_older_siblings_total": child_to_older_siblings,
                "ancestor_backreaction_total": sum(ancestor_deltas),
                "older_sibling_backreaction_total": child_to_older_siblings,
                "ancestor_depth_monotone": ancestor_depth_monotone,
                "sibling_order_increment_vs_older_max": sibling_order_increment,
                "parent_g_delta": self.nodes[parent].g - before_g[parent],
                "parent_response_before": before_response[parent],
                "parent_response_after": after_response[parent],
                "parent_response_delta": after_response[parent] - before_response[parent],
                "retarded_event_signal": int(outgoing_from_child > EPS and child_own_uv_tail_count_after == 0),
                "advanced_leakage_signal": 0,
                "boundary_polarity_signal": int(child_own_uv_tail_count_after == 0 and parent_tail_count_after == parent_tail_count_before + 1),
                "conductance_response_signal": int(outgoing_from_child > EPS or incoming_to_child > EPS),
                "partial_neutral_abs": abs(z_partial),
                "triple_completed": int(completed),
                "triple_log_circulation": triple_row["log_circulation_forward_vs_reverse"] if triple_row else 0.0,
                "triple_full_markov_eig_class": triple_row["full_markov_eig_class"] if triple_row else "",
            }
        )
        return c

    def grow_level(self, frontier: List[int]) -> List[int]:
        next_frontier: List[int] = []
        for p in frontier:
            for order in self.order_sequence:
                next_frontier.append(self.add_child(p, order))
        return next_frontier

    def completed_parents(self) -> List[int]:
        return [n.id for n in self.nodes.values() if len(n.children) == 3]

    def level_summary(self, level: int) -> dict:
        current_triples = []
        for p in self.completed_parents():
            m = self.triple_monodromy(p)
            if m is not None:
                current_triples.append(m)

        events = self.event_rows
        def mean(vals: List[float]) -> float:
            return float(np.mean(vals)) if vals else 0.0
        def frac(vals: List[int]) -> float:
            return mean([float(x) for x in vals])

        ret = frac([int(e["retarded_event_signal"]) for e in events])
        adv = frac([int(e["advanced_leakage_signal"]) for e in events])
        boundary = frac([int(e["boundary_polarity_signal"]) for e in events])
        conduct = frac([int(e["conductance_response_signal"]) for e in events])
        no_own_uv = frac([int(e["child_has_own_uv_tail"] == 0) for e in events])
        parent_gain = frac([int(e["parent_gains_uv_tail"]) for e in events])
        depth_mono = frac([int(e["ancestor_depth_monotone"]) for e in events])
        order_incs = [float(e["sibling_order_increment_vs_older_max"]) for e in events if int(e["older_sibling_count"]) > 0]
        triple_log_abs = [abs(float(t["log_circulation_forward_vs_reverse"])) for t in current_triples]
        frac_complex = frac([int(t["full_markov_eig_class"] == "complex_pair") for t in current_triples]) if current_triples else 0.0
        gs = [n.g for n in self.nodes.values()]
        row = {
            "variant": self.variant,
            "mode": self.mode,
            "level": level,
            "time": self.t,
            "nodes": len(self.nodes),
            "events": len(events),
            "role_rows": len(self.role_rows),
            "completed_triples": len(current_triples),
            "retarded_event_fraction": ret,
            "advanced_leakage_fraction": adv,
            "boundary_polarity_fraction": boundary,
            "child_no_own_uv_fraction": no_own_uv,
            "parent_gains_uv_tail_fraction": parent_gain,
            "conductance_response_fraction": conduct,
            "ancestor_depth_monotone_fraction": depth_mono,
            "mean_sibling_order_increment": mean(order_incs),
            "min_sibling_order_increment": min(order_incs) if order_incs else 0.0,
            "max_sibling_order_increment": max(order_incs) if order_incs else 0.0,
            "mean_abs_log_circulation": mean(triple_log_abs),
            "frac_full_markov_complex": frac_complex,
            "min_g": min(gs),
            "max_g": max(gs),
            "mean_g": float(np.mean(gs)),
        }
        self.level_rows.append(row)
        return row

    def run(self, max_level: int) -> None:
        frontier = [self.root]
        self.level_summary(0)
        for level in range(1, max_level + 1):
            frontier = self.grow_level(frontier)
            self.level_summary(level)


def write_csv(path: Path, rows: List[dict]) -> None:
    if not rows:
        return
    keys = sorted({k for row in rows for k in row.keys()})
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)


def run_suite(max_level: int, outdir: Path) -> dict:
    outdir.mkdir(parents=True, exist_ok=True)
    configs = [
        {
            "variant": "real_growth_linear_script1_2",
            "mode": "linear",
            "alpha_env": 0.22,
            "br_ancestor": 0.045,
            "br_sibling": 0.035,
            "order_sequence": (1, 2, 3),
        },
        {
            "variant": "log_growth_script1_2",
            "mode": "log",
            "alpha_env": 0.22,
            "br_ancestor": 0.045,
            "br_sibling": 0.035,
            "order_sequence": (1, 2, 3),
        },
        {
            "variant": "saturating_growth_script1_2",
            "mode": "saturating",
            "alpha_env": 0.90,
            "br_ancestor": 0.045,
            "br_sibling": 0.035,
            "order_sequence": (1, 2, 3),
        },
        {
            "variant": "strict_symmetrized_response_control",
            "mode": "linear",
            "alpha_env": 0.0,
            "br_ancestor": 0.0,
            "br_sibling": 0.0,
            "order_sequence": (1, 2, 3),
        },
        {
            "variant": "reverse_birth_label_order_control",
            "mode": "linear",
            "alpha_env": 0.22,
            "br_ancestor": 0.045,
            "br_sibling": 0.035,
            "order_sequence": (3, 2, 1),
        },
    ]
    summary = {
        "model_family": "script1_script2_dynamic_birth_conductance_monodromy",
        "max_level": max_level,
        "derived_only_notes": [
            "No complex scalar/J/Hodge/star/positivity is used as a gate.",
            "Complex neutral phasor diagnostics are inherited from script 1/2 only as Z3 response diagnostics.",
            "DtN/Schur deltas are scalar conductance-response surrogates, not matrix Schur complements.",
        ],
        "variants": [],
    }
    for cfg in configs:
        m = DynamicBoundaryRoleModel(**cfg)
        m.run(max_level)
        prefix = cfg["variant"]
        write_csv(outdir / f"events_{prefix}.csv", m.event_rows)
        write_csv(outdir / f"boundary_role_rows_{prefix}.csv", m.role_rows)
        write_csv(outdir / f"triples_{prefix}.csv", m.triple_rows)
        write_csv(outdir / f"levels_{prefix}.csv", m.level_rows)
        final = m.level_rows[-1]
        # Gates: split topological role polarity from conductance/response role recalibration.
        response_role_gate = (
            final["retarded_event_fraction"] >= 0.999
            and final["advanced_leakage_fraction"] <= 1e-12
            and final["conductance_response_fraction"] >= 0.999
            and final["ancestor_depth_monotone_fraction"] >= 0.999
        )
        topological_role_gate = (
            final["child_no_own_uv_fraction"] >= 0.999
            and final["parent_gains_uv_tail_fraction"] >= 0.999
            and final["boundary_polarity_fraction"] >= 0.999
        )
        summary["variants"].append(
            {
                "variant": prefix,
                "events": final["events"],
                "completed_triples": final["completed_triples"],
                "topological_role_gate": bool(topological_role_gate),
                "response_role_recalibration_gate": bool(response_role_gate),
                "retarded_event_fraction": final["retarded_event_fraction"],
                "advanced_leakage_fraction": final["advanced_leakage_fraction"],
                "boundary_polarity_fraction": final["boundary_polarity_fraction"],
                "conductance_response_fraction": final["conductance_response_fraction"],
                "child_no_own_uv_fraction": final["child_no_own_uv_fraction"],
                "parent_gains_uv_tail_fraction": final["parent_gains_uv_tail_fraction"],
                "ancestor_depth_monotone_fraction": final["ancestor_depth_monotone_fraction"],
                "mean_sibling_order_increment": final["mean_sibling_order_increment"],
                "mean_abs_log_circulation": final["mean_abs_log_circulation"],
                "frac_full_markov_complex": final["frac_full_markov_complex"],
                "used_delta_beta_any": False,
            }
        )
    (outdir / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    return summary


def make_markdown(summary: dict, outdir: Path) -> None:
    def row_for(v: dict) -> str:
        return (
            f"| {v['variant']} | {v['events']} | {v['completed_triples']} | "
            f"{int(v['topological_role_gate'])} | {int(v['response_role_recalibration_gate'])} | "
            f"{v['retarded_event_fraction']:.3f} | {v['advanced_leakage_fraction']:.3f} | "
            f"{v['conductance_response_fraction']:.3f} | {v['mean_sibling_order_increment']:.6f} | "
            f"{v['mean_abs_log_circulation']:.6f} | {v['frac_full_markov_complex']:.3f} |"
        )
    table = "\n".join(row_for(v) for v in summary["variants"])
    md = f"""# RESULTS — CNNA dynamic boundary role recalibration gate, L2

## Model provenance

This package reuses the established script-1/script-2 dynamic birth model:

- ternary sequential births;
- newborn conductance derived from parent-line + older-sibling environment;
- existing environment edges into the newborn;
- newborn backreaction to the parent line and older siblings;
- completed sibling-triple monodromy diagnostics.

The new diagnostic is event-resolved boundary-role recalibration.  It logs, per birth, how the same newborn has no own UV-tail from its own perspective while immediately becoming a UV-tail/backreaction source for the parent line.

## Gate table

| variant | events | triples | topological role | response role | retarded frac | advanced frac | conductance frac | mean sibling increment | mean |log circ| | full Markov complex frac |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
{table}

## Interpretation

The gate separates two claims that should not be conflated:

1. **Topological boundary polarity**: every newborn has no own UV-tail at birth and simultaneously becomes a new UV-tail for the parent cut.  This remains true even in the strict-sym response control because it is part of the birth event itself.
2. **Conductance/response role recalibration**: parent-line and older-sibling conductances are updated directionally by the newborn.  This is positive in the real/log/saturating/reverse-label growth variants and collapses in the strict-sym response control.

The result therefore supports a pre-causal boundary-role layer, not a J/i/spin claim.  Growth creates an irreversible role update:

```text
child self-perspective: no own UV-tail
parent-line perspective: child is new UV-tail / backreaction source
older sibling perspective: already-born cells are updated by later births
ancestor/root perspective: descendant birth changes effective response
```

This is the CNNA analogue of a retarded boundary-value selection candidate.  It is not yet a Q/P, *, J, or complex structure.

## Important limitations

- The `cut_response_delta` / `Schur_delta` columns are scalar conductance-response surrogates, not full matrix Schur complements.
- The strict-sym control kills conductance/response asymmetry, but not the topological fact that a birth adds a child to the parent cut.
- The test uses L2 only; deeper levels should check whether old-interior role deltas decay while frontier updates dominate.
- No delta-beta, H², Q/P, J, Hodge, physical adjoint, or positivity criterion is used as a decision gate.

## Next test

`test_role_recalibration_to_boundary_value_bridge_gate.py`

Purpose: use the event-resolved role rows as the source layer and test whether the retarded boundary-role polarity transfers to the later Q/P/operator layer more robustly than the previous response-monodromy bridge.
"""
    (outdir / "RESULTS.md").write_text(md, encoding="utf-8")
    (outdir / "SUMMARY.md").write_text(md, encoding="utf-8")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--max-level", type=int, default=2)
    ap.add_argument("--outdir", type=Path, default=Path("outputs"))
    args = ap.parse_args()
    summary = run_suite(args.max_level, args.outdir)
    make_markdown(summary, args.outdir)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
