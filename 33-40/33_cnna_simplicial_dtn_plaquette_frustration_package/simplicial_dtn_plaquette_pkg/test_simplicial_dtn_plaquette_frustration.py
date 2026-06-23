#!/usr/bin/env python3
"""
CNNA Growth Test: simplicial DtN plaquette frustration.

Purpose
-------
This is the first bridge test after the pure provenance operator-closure
obstruction.  It couples all earlier strands:

1. real event-resolved growth from tests 1/2,
2. provenance-generated 2-simplices from test 3,
3. the obstruction from test 9 that minimal axis geometry alone is flat,
4. dynamic record/live DtN response germs from tests 22/23,
5. the post-23 lesson that node/edge provenance alone remains gradient-like.

The test does NOT set i, J, a complex phase, a Hilbert space, a C*-norm, or a
*-algebra.  It asks only whether the geometry glued from provenance supplies a
2-dimensional plaquette carrier on which DtN/Record-Live response operators have
non-telescoping, noncommuting curvature.

Core construction
-----------------
For every completed local cell p, build its children c1,c2,c3.  Once those
children themselves are completed local cells, the provenance generates oriented
2-simplices:

  sibling face:       (c1,c2,c3)
  parent-fan faces:   (p,c1,c2), (p,c2,c3), (p,c3,c1)

Each vertex v of a face carries a 3-port DtN operator S_v, selected from one of

  record, live, handoff = live-record, aging = live-live_at_completion.

For an oriented edge a->b define the response-derived edge generator

  A_ab = S_b - S_a.

This is intentionally a coboundary at first order, hence

  A_ab + A_bc + A_ca = 0.

Therefore any nonzero plaquette signal is not allowed to come from first-order
height/gradient circulation.  The measured second-order core is the commutator

  [A_ab, A_bc] = A_ab A_bc - A_bc A_ab.

The finite plaquette holonomy surrogate is

  U_ab = exp(eps A_ab),
  P_abc = U_ca U_bc U_ab.

Nonzero skew(P_abc), complex eigenpairs, and nonzero commutator norm diagnose a
real operatorial plaquette-frustration candidate.  A positive result here is
still not a J proof; it only restores the missing 2-simplex/face carrier.

Controls
--------
- sequential real growth
- symmetrized birth
- no backreaction
- diagonal-only DtN operators
- trace-scalar DtN operators
- random same-level triples
- reversed face orientation sign check
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


def address_tuple(model: RealGrowth, node_id: int) -> Tuple[int, ...]:
    out: List[int] = []
    cur = int(node_id)
    while model.nodes[cur].parent is not None:
        out.append(int(model.nodes[cur].birth_order))
        cur = int(model.nodes[cur].parent)
    return tuple(reversed(out))


@dataclass(frozen=True)
class Face:
    mode: str
    orientation: str
    level: int
    base_parent: int
    vertices: Tuple[int, int, int]


def build_faces(model: RealGrowth, cells: Dict[int, LocalCell], *, include_random: int, seed: int) -> List[Face]:
    faces: List[Face] = []
    cell_ids = set(cells.keys())
    for p in sorted(model.completed_parent_ids()):
        children = model.child_ids_ordered(p)
        if len(children) != 3:
            continue
        c1, c2, c3 = children
        if c1 in cell_ids and c2 in cell_ids and c3 in cell_ids:
            level = int(model.nodes[c1].level)
            faces.append(Face("sibling_triangle", "forward", level, p, (c1, c2, c3)))
            faces.append(Face("sibling_triangle", "reverse", level, p, (c1, c3, c2)))
        if p in cell_ids:
            fan = [(c1, c2), (c2, c3), (c3, c1)]
            for a, b in fan:
                if a in cell_ids and b in cell_ids:
                    level = int(model.nodes[p].level)
                    faces.append(Face("parent_fan_triangle", "forward", level, p, (p, a, b)))
                    faces.append(Face("parent_fan_triangle", "reverse", level, p, (p, b, a)))

    if include_random > 0:
        rng = random.Random(seed)
        by_level: Dict[int, List[int]] = {}
        for cid, c in cells.items():
            by_level.setdefault(c.level, []).append(cid)
        levels = [lv for lv, xs in by_level.items() if len(xs) >= 3]
        for k in range(include_random):
            lv = rng.choice(levels)
            tri = tuple(rng.sample(by_level[lv], 3))
            faces.append(Face("random_same_level_triangle", "forward", lv, -1, tri))
    return faces


def vertex_operator(cell: LocalCell, source_mode: str, reduction: str) -> np.ndarray:
    S = cell.matrix_for_source(source_mode)
    S = sym(S)
    if reduction == "full":
        return S
    if reduction == "diagonal":
        return np.diag(np.diag(S))
    if reduction == "trace_scalar":
        return np.eye(3) * (float(np.trace(S)) / 3.0)
    raise ValueError(f"unknown reduction {reduction}")


def face_metrics(face: Face, cells: Dict[int, LocalCell], source_mode: str, reduction: str, step: float) -> dict:
    a, b, c = face.vertices
    Sa = vertex_operator(cells[a], source_mode, reduction)
    Sb = vertex_operator(cells[b], source_mode, reduction)
    Sc = vertex_operator(cells[c], source_mode, reduction)
    Aab = Sb - Sa
    Abc = Sc - Sb
    Aca = Sa - Sc
    first = Aab + Abc + Aca
    comm = Aab @ Abc - Abc @ Aab
    comm_alt = Abc @ Aca - Aca @ Abc
    comm_cyclic = comm + comm_alt + (Aca @ Aab - Aab @ Aca)
    Uab = sym_expm(Aab, step)
    Ubc = sym_expm(Abc, step)
    Uca = sym_expm(Aca, step)
    P = Uca @ Ubc @ Uab
    Pinv_order = Uab @ Ubc @ Uca
    I = np.eye(3)
    PmI = P - I
    sk = skew(P)
    sy = sym(PmI)
    vals = np.linalg.eigvals(P)
    detP = float(np.linalg.det(P))
    return {
        "mode": face.mode,
        "orientation": face.orientation,
        "face_level": face.level,
        "base_parent": face.base_parent,
        "vertices": " ".join(map(str, face.vertices)),
        "addr0": ".".join(map(str, address_tuple(cells[a].model, a))),
        "addr1": ".".join(map(str, address_tuple(cells[b].model, b))),
        "addr2": ".".join(map(str, address_tuple(cells[c].model, c))),
        "source_mode": source_mode,
        "reduction": reduction,
        "first_order_closure_norm": fro(first),
        "edge_generator_mean_norm": mean([fro(Aab), fro(Abc), fro(Aca)]),
        "commutator_norm": fro(comm),
        "cyclic_commutator_norm": fro(comm_cyclic),
        "plaquette_norm": fro(PmI),
        "plaquette_skew_norm": fro(sk),
        "plaquette_sym_norm": fro(sy),
        "skew_fraction": fro(sk) / (fro(PmI) + EPS),
        "max_imag_eig": max_imag_eig(P),
        "has_complex_pair": int(max_imag_eig(P) > 1e-8),
        "det_minus_one_abs": abs(detP - 1.0),
        "reverse_order_mismatch": fro(P - Pinv_order),
        "eigvals_real": " ".join(f"{float(np.real(x)):.9g}" for x in vals),
        "eigvals_imag": " ".join(f"{float(np.imag(x)):.9g}" for x in vals),
    }


def summarize(rows: List[dict], keys: List[str]) -> List[dict]:
    groups: Dict[Tuple, List[dict]] = {}
    for r in rows:
        k = tuple(r[x] for x in keys)
        groups.setdefault(k, []).append(r)
    out: List[dict] = []
    for k, rs in sorted(groups.items(), key=lambda kv: kv[0]):
        d = {keys[i]: k[i] for i in range(len(keys))}
        d.update(
            count=len(rs),
            mean_first_order_closure_norm=mean(r["first_order_closure_norm"] for r in rs),
            mean_edge_generator_norm=mean(r["edge_generator_mean_norm"] for r in rs),
            mean_commutator_norm=mean(r["commutator_norm"] for r in rs),
            p95_commutator_norm=perc((r["commutator_norm"] for r in rs), 95),
            mean_cyclic_commutator_norm=mean(r["cyclic_commutator_norm"] for r in rs),
            mean_plaquette_norm=mean(r["plaquette_norm"] for r in rs),
            p95_plaquette_norm=perc((r["plaquette_norm"] for r in rs), 95),
            mean_skew_norm=mean(r["plaquette_skew_norm"] for r in rs),
            mean_skew_fraction=mean(r["skew_fraction"] for r in rs),
            mean_max_imag_eig=mean(r["max_imag_eig"] for r in rs),
            frac_complex_pair=mean(r["has_complex_pair"] for r in rs),
            mean_reverse_order_mismatch=mean(r["reverse_order_mismatch"] for r in rs),
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


def run(args: argparse.Namespace) -> str:
    outdir: Path = args.outdir
    if outdir.exists():
        shutil.rmtree(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    rows: List[dict] = []
    model_summaries: List[dict] = []
    controls = ["real_growth", "symmetrized_birth", "no_backreaction"]
    source_modes = ["record", "live", "handoff", "aging"]
    reductions = ["full", "diagonal", "trace_scalar"]

    for control in controls:
        model = build_model(control, args.max_level, args.mode)
        cells = build_cells(model, args.max_level)
        faces = build_faces(model, cells, include_random=args.random_faces, seed=args.seed)
        model_summaries.append({
            "control": control,
            "nodes": len(model.nodes),
            "completed_cells": len(cells),
            "faces_total": len(faces),
            "sibling_faces": sum(1 for f in faces if f.mode == "sibling_triangle" and f.orientation == "forward"),
            "parent_fan_faces": sum(1 for f in faces if f.mode == "parent_fan_triangle" and f.orientation == "forward"),
            "random_faces": sum(1 for f in faces if f.mode == "random_same_level_triangle"),
        })
        for source_mode in source_modes:
            for reduction in reductions:
                for face in faces:
                    r = face_metrics(face, cells, source_mode, reduction, args.step)
                    r["control"] = control
                    r["max_level"] = args.max_level
                    r["step"] = args.step
                    rows.append(r)

    by_main = summarize(rows, ["control", "mode", "orientation", "source_mode", "reduction"])
    by_source = summarize(rows, ["control", "source_mode", "reduction"])
    by_face = summarize(rows, ["control", "mode", "source_mode", "reduction"])
    by_level = summarize(rows, ["control", "mode", "face_level", "source_mode", "reduction"])

    write_csv(outdir / "plaquette_face_rows.csv", rows)
    write_csv(outdir / "plaquette_summary_main.csv", by_main)
    write_csv(outdir / "plaquette_summary_by_source.csv", by_source)
    write_csv(outdir / "plaquette_summary_by_face.csv", by_face)
    write_csv(outdir / "plaquette_summary_by_level.csv", by_level)
    write_csv(outdir / "model_summaries.csv", model_summaries)

    def line_for(d: dict) -> str:
        return (
            f"{d['control']} | {d.get('mode','*')} | {d.get('source_mode','*')} | {d.get('reduction','*')}: "
            f"count={d['count']}, comm={d['mean_commutator_norm']:.6g}, "
            f"plaq={d['mean_plaquette_norm']:.6g}, skew={d['mean_skew_norm']:.6g}, "
            f"skew_frac={d['mean_skew_fraction']:.4f}, complex_frac={d['frac_complex_pair']:.3f}, "
            f"first={d['mean_first_order_closure_norm']:.2e}"
        )

    lines = [
        "CNNA SIMPLICIAL DTN PLAQUETTE FRUSTRATION TEST",
        f"max_level={args.max_level}, mode={args.mode}, step={args.step}",
        "",
        "MODEL SUMMARIES",
    ]
    for m in model_summaries:
        lines.append(
            f"  {m['control']}: nodes={m['nodes']}, completed={m['completed_cells']}, "
            f"sibling_faces={m['sibling_faces']}, parent_fan_faces={m['parent_fan_faces']}, random={m['random_faces']}"
        )
    lines += ["", "SELECTED FULL-DTN FORWARD FACE SUMMARIES"]
    for r in by_main:
        if r["orientation"] == "forward" and r["reduction"] == "full" and r["source_mode"] in ("live", "handoff", "aging") and r["mode"] in ("sibling_triangle", "parent_fan_triangle", "random_same_level_triangle"):
            lines.append("  " + line_for(r))
    lines += ["", "DIAGONAL AND TRACE-SCALAR CONTROLS, REAL GROWTH / LIVE"]
    for r in by_main:
        if r["control"] == "real_growth" and r["orientation"] == "forward" and r["source_mode"] == "live" and r["reduction"] in ("diagonal", "trace_scalar"):
            lines.append("  " + line_for(r))
    lines += ["", "INTERPRETATION GUARDRAILS"]
    lines += [
        "  first_order_closure_norm should be numerical zero by construction: A_ab+A_bc+A_ca telescopes.",
        "  nonzero commutator/skew/complex-pair diagnostics are therefore second-order operatorial plaquette effects, not height-gradient circulation.",
        "  positive plaquette frustration is not a J proof; it is the missing face-level carrier needed before any J extraction.",
    ]
    summary = "\n".join(lines)
    (outdir / "SUMMARY.txt").write_text(summary, encoding="utf-8")

    results = f"""# Results: simplicial DtN plaquette frustration

This diagnostic couples real event-resolved growth, provenance-generated 2-simplices, and dynamic DtN/Record-Live response operators.

It deliberately does **not** set `i`, `J`, a complex Hilbert space, a `C*` norm, or a `*`-algebra.

## Core gate

For each oriented face `(a,b,c)` and selected DtN source `S_v`, define

```text
A_ab = S_b - S_a
P_abc = exp(eps A_ca) exp(eps A_bc) exp(eps A_ab)
```

Since `A_ab + A_bc + A_ca = 0`, first-order circulation is forced to telescope.  The tested signal is therefore the second-order/noncommuting plaquette part:

```text
[A_ab,A_bc] and skew(P_abc)
```

## Output files

- `plaquette_face_rows.csv`
- `plaquette_summary_main.csv`
- `plaquette_summary_by_source.csv`
- `plaquette_summary_by_face.csv`
- `plaquette_summary_by_level.csv`
- `model_summaries.csv`
- `SUMMARY.txt`

## Current run

```text
{summary}
```

## Reading rule

A useful positive result is:

```text
full DtN live/handoff plaquette commutator > diagonal/trace controls,
real growth > symmetrized/no-backreaction where appropriate,
first-order closure ~ 0,
orientation reversal changes the skew-sign-related rows consistently.
```

A useful negative result is:

```text
full DtN plaquette signal collapses to diagonal/trace controls,
or all signals are explained by random same-level faces.
```

Either way, this is the correct first test after the pure-provenance operator-closure obstruction: it asks whether the noncommutative core lives on provenance-glued 2-simplices.
"""
    (outdir / "RESULTS_simplicial_dtn_plaquette_frustration.md").write_text(results, encoding="utf-8")
    return summary


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--max-level", type=int, default=6)
    ap.add_argument("--mode", choices=["linear", "log", "saturating"], default="linear")
    ap.add_argument("--step", type=float, default=0.04)
    ap.add_argument("--random-faces", type=int, default=80)
    ap.add_argument("--seed", type=int, default=1729)
    ap.add_argument("--outdir", type=Path, default=Path("simplicial_dtn_plaquette_out_L6"))
    args = ap.parse_args()
    print(run(args))


if __name__ == "__main__":
    main()
