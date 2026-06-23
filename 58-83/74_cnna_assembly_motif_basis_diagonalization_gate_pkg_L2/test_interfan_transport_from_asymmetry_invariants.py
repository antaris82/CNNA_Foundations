#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import math
import shutil
import zipfile
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

import numpy as np

import cnna_non_shelling_core as core

EPS = 1e-12
Face = Tuple[int, int, int]


def write_csv(path: Path, rows: List[dict]) -> None:
    if not rows:
        path.write_text('', encoding='utf-8')
        return
    fields: List[str] = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    with path.open('w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for row in rows:
            w.writerow(row)


class StrictSymmetrizedGrowth(core.DynamicProvenanceGrowth):
    """A stricter symmetrized control than the historical symmetrized_birth switch.

    The historical generated core still allows sequential ancestor updates during a
    sibling fan and can therefore retain residual order traces.  This control keeps
    the same outward geometry labels, but suppresses sibling-directed and
    backreaction-produced asymmetry as far as this lightweight diagnostic can do
    without replacing the whole growth engine.
    """

    def __init__(self, *args, **kwargs):
        kwargs['growth_rule'] = 'symmetrized_birth'
        kwargs['br_sibling'] = 0.0
        kwargs['br_ancestor'] = 0.0
        super().__init__(*args, **kwargs)

    def add_child(self, parent: int, order: int) -> int:
        # same geometric birth labels as before, but no older-sibling environment,
        # no sibling backreaction, no ancestor backreaction update.
        self.t += 1
        older: List[int] = []
        env_load = self.birth_environment_load(parent, older)
        birth_g = self.child_conductance_from_env(env_load)
        pos, r, e1, e2 = self.child_position(parent, order, older)
        child = self._new_node(parent, self.nodes[parent].level + 1, order, birth_g, pos, r, e1, e2)
        c = child.id
        total_env = env_load + EPS
        for d, a in enumerate(self.parent_line(parent), start=1):
            contrib = self.nodes[a].g * (self.ancestor_decay ** (d - 1))
            # symmetric feed-forward environment only; no reciprocal birth-order residue.
            self.directed_edges[(a, c)] += self.alpha_env * contrib / total_env * birth_g
        self.birth_events.append({
            't': self.t,
            'parent': parent,
            'child': c,
            'order': order,
            'older_siblings': [],
            'env_load': env_load,
            'birth_g': birth_g,
            'level': child.level,
            'strict_sym_control': True,
        })
        return c


def ancestors(model: core.DynamicProvenanceGrowth, node: int) -> List[int]:
    out: List[int] = []
    cur = model.nodes[node].parent
    while cur is not None:
        out.append(cur)
        cur = model.nodes[cur].parent
    return out


def is_descendant(model: core.DynamicProvenanceGrowth, child: int, parent: int) -> bool:
    cur = model.nodes[child].parent
    while cur is not None:
        if cur == parent:
            return True
        cur = model.nodes[cur].parent
    return False


def descendant_shell_loads(model: core.DynamicProvenanceGrowth) -> Dict[int, float]:
    loads = {i: 0.0 for i in model.nodes}
    for v, nv in model.nodes.items():
        cur = nv.parent
        dist = 1
        while cur is not None:
            # same semantic family as shell-normalized inverse-square: near UV matters more;
            # distant descendants are present but damped.
            k = 1.0 / ((model.branching ** max(0, dist - 1)) * (dist * dist))
            loads[cur] += nv.birth_g * k
            cur = model.nodes[cur].parent
            dist += 1
    return loads


def normalized_spread(xs: Iterable[float]) -> float:
    vals = [float(x) for x in xs]
    if not vals:
        return 0.0
    m = max(abs(x) for x in vals) + EPS
    return (max(vals) - min(vals)) / m


def fan_invariants(model: core.DynamicProvenanceGrowth, shell_load: Dict[int, float]) -> Dict[int, dict]:
    out: Dict[int, dict] = {}
    for parent in model.completed_parent_ids():
        ch = model.child_ids_ordered(parent)
        raw_imb = 0.0
        total = 0.0
        for i, a in enumerate(ch):
            for b in ch[i + 1:]:
                wab = model.directed_edges.get((a, b), 0.0)
                wba = model.directed_edges.get((b, a), 0.0)
                raw_imb += abs(wab - wba)
                total += wab + wba
        g_updates = [model.nodes[c].g - model.nodes[c].birth_g for c in ch]
        g_values = [model.nodes[c].g for c in ch]
        b_values = [model.nodes[c].birth_g for c in ch]
        shell_values = [shell_load.get(c, 0.0) for c in ch]
        parent_line_updates = [model.nodes[a].g - model.nodes[a].birth_g for a in model.parent_line(parent)]
        age_grad = sum(abs(x) / (i + 1) for i, x in enumerate(parent_line_updates))
        inv = {
            'parent': parent,
            'fan_directed_sibling_imbalance_raw': raw_imb,
            'fan_directed_sibling_imbalance_norm': raw_imb / (total + EPS),
            'fan_live_record_spread': normalized_spread(g_updates),
            'fan_conductance_spread': normalized_spread(g_values),
            'fan_birth_conductance_spread': normalized_spread(b_values),
            'fan_desc_shell_spread': normalized_spread(shell_values),
            'fan_parent_line_aging_gradient': age_grad,
        }
        inv['fan_asymmetry_invariant'] = (
            1.40 * inv['fan_directed_sibling_imbalance_norm']
            + 0.90 * inv['fan_live_record_spread']
            + 0.70 * inv['fan_desc_shell_spread']
            + 0.45 * inv['fan_conductance_spread']
            + 0.35 * inv['fan_parent_line_aging_gradient']
        )
        out[parent] = inv
    return out


def face_parent_ids(model: core.DynamicProvenanceGrowth, f: Face) -> List[int]:
    ps: List[int] = []
    for v in f:
        p = model.nodes[v].parent
        if p is not None:
            ps.append(p)
    return ps


def mean_value(xs: Iterable[float]) -> float:
    vals = [float(x) for x in xs]
    return sum(vals) / len(vals) if vals else 0.0


def face_mean(model: core.DynamicProvenanceGrowth, f: Face, values: Dict[int, float]) -> float:
    return mean_value(values.get(v, 0.0) for v in f)


def nonreciprocal_residue(model: core.DynamicProvenanceGrowth, f: Face, g: Face) -> Tuple[float, float]:
    raw = 0.0
    total = 0.0
    for a in f:
        for b in g:
            wab = model.directed_edges.get((a, b), 0.0)
            wba = model.directed_edges.get((b, a), 0.0)
            raw += abs(wab - wba)
            total += wab + wba
    return raw, raw / (total + EPS)


def candidate_asymmetry_features(
    model: core.DynamicProvenanceGrowth,
    row: dict,
    fan_inv: Dict[int, dict],
    shell_load: Dict[int, float],
) -> dict:
    if not row.get('face_b') or 'perm=' in str(row.get('face_b', '')):
        fb_txt = str(row.get('face_b', '')).split('perm=')[0].strip()
    else:
        fb_txt = str(row.get('face_b', '')).strip()
    try:
        f = core.parse_face_string(str(row.get('face_a', '')))
        g = core.parse_face_string(fb_txt)
    except Exception:
        return {
            'asymmetry_candidate': False,
            'A_invariant': 0.0,
            'A_gate': False,
            'A_rank_score': float(row.get('response_score', 0.0)),
        }

    parents = set(face_parent_ids(model, f) + face_parent_ids(model, g))
    fan_vals = [fan_inv[p]['fan_asymmetry_invariant'] for p in parents if p in fan_inv]
    fan_directed = [fan_inv[p]['fan_directed_sibling_imbalance_norm'] for p in parents if p in fan_inv]
    fan_shell = [fan_inv[p]['fan_desc_shell_spread'] for p in parents if p in fan_inv]
    fan_aging = [fan_inv[p]['fan_live_record_spread'] for p in parents if p in fan_inv]
    fan_parent_grad = [fan_inv[p]['fan_parent_line_aging_gradient'] for p in parents if p in fan_inv]
    nr_raw, nr_norm = nonreciprocal_residue(model, f, g)
    aging_by_node = {i: n.g - n.birth_g for i, n in model.nodes.items()}
    live_record_gap = abs(face_mean(model, f, aging_by_node) - face_mean(model, g, aging_by_node))
    shell_gap = abs(face_mean(model, f, shell_load) - face_mean(model, g, shell_load))
    conductance_gap = abs(mean_value(model.nodes[v].g for v in f) - mean_value(model.nodes[v].g for v in g))
    directed_imbalance = float(row.get('directed_imbalance', 0.0) or 0.0)
    trans = float(row.get('transverse_complementarity', 0.0) or 0.0)
    addr = float(row.get('address_similarity', 0.0) or 0.0)
    centroid_distance = float(row.get('centroid_distance', 0.0) or 0.0)
    A = (
        1.30 * directed_imbalance
        + 1.10 * mean_value(fan_directed)
        + 0.95 * nr_norm
        + 0.70 * mean_value(fan_vals)
        + 0.45 * live_record_gap
        + 0.38 * shell_gap
        + 0.30 * mean_value(fan_shell)
        + 0.24 * mean_value(fan_aging)
        + 0.18 * conductance_gap
        + 0.12 * mean_value(fan_parent_grad)
    )
    complement_gate = trans >= 0.18 and centroid_distance >= 0.55 and directed_imbalance >= 0.035
    invariant_gate = A >= 0.18 and mean_value(fan_directed) >= 0.001
    A_gate = bool(complement_gate and invariant_gate)
    # Ranking score deliberately does NOT include delta_beta; this tests whether the
    # real response/provenance asymmetry selects topology-effective moves without
    # looking at topology after the fact.
    A_rank_score = (
        1.00 * float(row.get('response_score', 0.0) or 0.0)
        + 3.50 * A
        + 1.20 * trans
        + 0.35 * addr
        - 0.03 * centroid_distance
    ) if A_gate else (
        0.25 * float(row.get('response_score', 0.0) or 0.0)
        + 0.50 * A
        + 0.15 * trans
    )
    return {
        'asymmetry_candidate': True,
        'A_invariant': A,
        'A_gate': A_gate,
        'A_rank_score': A_rank_score,
        'A_fan_mean': mean_value(fan_vals),
        'A_fan_directed_mean': mean_value(fan_directed),
        'A_nonreciprocal_raw': nr_raw,
        'A_nonreciprocal_norm': nr_norm,
        'A_live_record_gap': live_record_gap,
        'A_shell_gap': shell_gap,
        'A_conductance_gap': conductance_gap,
        'A_fan_shell_mean': mean_value(fan_shell),
        'A_fan_aging_mean': mean_value(fan_aging),
        'A_parent_grad_mean': mean_value(fan_parent_grad),
        'A_complement_gate': complement_gate,
        'A_invariant_gate': invariant_gate,
    }


def rerank_by_asymmetry(rows: List[dict]) -> None:
    legal = [r for r in rows if r.get('move_class') != 'illegal' and r.get('status') == 'ok']
    ranked_all = sorted(legal, key=lambda r: float(r.get('A_rank_score', 0.0)), reverse=True)
    for i, r in enumerate(ranked_all, start=1):
        r['A_rank_legal'] = i
    gated = [r for r in legal if r.get('A_gate')]
    ranked_gated = sorted(gated, key=lambda r: float(r.get('A_rank_score', 0.0)), reverse=True)
    for i, r in enumerate(ranked_gated, start=1):
        r['A_rank_gated'] = i
    for r in rows:
        r.setdefault('A_rank_legal', '')
        r.setdefault('A_rank_gated', '')


def top_rows(rows: List[dict], predicate, n: int = 5, key: str = 'A_rank_score') -> List[dict]:
    sub = [r for r in rows if predicate(r)]
    return sorted(sub, key=lambda r: float(r.get(key, 0.0)), reverse=True)[:n]


def compact_row(r: Optional[dict]) -> dict:
    if not r:
        return {}
    keys = [
        'candidate_id', 'move_class', 'response_rank_legal', 'A_rank_legal', 'A_rank_gated',
        'response_score', 'A_rank_score', 'A_invariant', 'A_gate',
        'delta_beta1', 'delta_beta2', 'new_beta1', 'new_beta2', 'delta_boundary_faces',
        'directed_imbalance', 'transverse_complementarity', 'address_similarity',
        'A_fan_directed_mean', 'A_nonreciprocal_norm', 'A_live_record_gap', 'A_shell_gap',
        'face_a', 'face_b',
    ]
    return {k: r.get(k, '') for k in keys}


def apply_and_measure(model: core.DynamicProvenanceGrowth, K: core.SimplicialComplex, row: Optional[dict], source: str) -> dict:
    if row is None:
        return {'status': 'missing'}
    L, reason, encoded = core.apply_candidate_row(K, row)
    base = core.full_metrics(model, K, source)
    if L is None:
        out = compact_row(row)
        out.update({'status': reason, 'encoded_move': encoded})
        return out
    after = core.full_metrics(model, L, source)
    out = compact_row(row)
    out.update({
        'status': reason,
        'encoded_move': encoded,
        'after_beta0': after['beta0'],
        'after_beta1': after['beta1'],
        'after_beta2': after['beta2'],
        'after_beta3': after['beta3'],
        'after_boundary_fraction': after['boundary_fraction'],
        'after_edge_link_cycle_fraction': after['edge_link_cycle_fraction'],
        'after_K_mean': after['K_mean'],
        'after_harmonic_ratio': after['harmonic_ratio'],
        'after_exact_residual_ratio': after['exact_residual_ratio'],
        'measured_delta_beta1': after['beta1'] - base['beta1'],
        'measured_delta_beta2': after['beta2'] - base['beta2'],
        'measured_delta_beta3': after['beta3'] - base['beta3'],
    })
    return out


def build_model(variant: str, args: argparse.Namespace) -> core.DynamicProvenanceGrowth:
    if variant == 'strict_symmetrized_control':
        return StrictSymmetrizedGrowth(mode=args.mode, transverse_amp=args.transverse_amp)
    growth_rule = {
        'real_growth': 'real_growth',
        'historical_symmetrized_birth': 'symmetrized_birth',
        'no_backreaction': 'no_backreaction',
    }[variant]
    return core.DynamicProvenanceGrowth(mode=args.mode, growth_rule=growth_rule, transverse_amp=args.transverse_amp)


def run_variant(variant: str, args: argparse.Namespace, out: Path) -> dict:
    model = build_model(variant, args)
    model.grow(args.max_level)
    K = core.build_dynamic_outward_ngf_complex(model)
    rows, audit = core.enumerate_moves(
        model,
        K,
        source=args.source,
        max_boundary_faces=args.max_boundary_faces,
        max_single_vertices=args.max_single_vertices,
        max_pair_candidates=args.max_pair_candidates,
        max_rows=args.max_rows,
    )
    shell_load = descendant_shell_loads(model)
    fan_inv = fan_invariants(model, shell_load)
    for r in rows:
        if r.get('move_class') in {'handle_candidate', 'quotient_candidate'}:
            r.update(candidate_asymmetry_features(model, r, fan_inv, shell_load))
        else:
            r.update({
                'asymmetry_candidate': False,
                'A_invariant': 0.0,
                'A_gate': False,
                'A_rank_score': 0.25 * float(r.get('response_score', 0.0) or 0.0),
                'A_fan_mean': 0.0,
                'A_fan_directed_mean': 0.0,
                'A_nonreciprocal_raw': 0.0,
                'A_nonreciprocal_norm': 0.0,
                'A_live_record_gap': 0.0,
                'A_shell_gap': 0.0,
                'A_conductance_gap': 0.0,
                'A_fan_shell_mean': 0.0,
                'A_fan_aging_mean': 0.0,
                'A_parent_grad_mean': 0.0,
                'A_complement_gate': False,
                'A_invariant_gate': False,
            })
    rerank_by_asymmetry(rows)

    legal = [r for r in rows if r.get('move_class') != 'illegal' and r.get('status') == 'ok']
    gated = [r for r in legal if r.get('A_gate')]
    topology_effective = [r for r in legal if r.get('move_class') in {'handle_candidate', 'quotient_candidate'} and (int(r.get('delta_beta1', 0)) > 0 or int(r.get('delta_beta2', 0)) > 0)]
    gated_topology_effective = [r for r in topology_effective if r.get('A_gate')]
    best_A = top_rows(rows, lambda r: r.get('status') == 'ok' and r.get('move_class') != 'illegal', n=1)[0] if legal else None
    best_A_gated = top_rows(rows, lambda r: r.get('status') == 'ok' and r.get('A_gate'), n=1)[0] if gated else None
    best_topology_effective_A = top_rows(rows, lambda r: r in topology_effective, n=1)[0] if topology_effective else None
    best_delta_beta2 = sorted(topology_effective, key=lambda r: (int(r.get('delta_beta2', 0)), float(r.get('A_rank_score', 0.0))), reverse=True)[0] if topology_effective else None

    applied = [
        apply_and_measure(model, K, best_A, args.source),
        apply_and_measure(model, K, best_A_gated, args.source),
        apply_and_measure(model, K, best_topology_effective_A, args.source),
        apply_and_measure(model, K, best_delta_beta2, args.source),
    ]
    for label, row in zip(['top_A_any', 'top_A_gated', 'top_A_topology_effective', 'max_delta_beta2_then_A'], applied):
        row['selection'] = label

    vout = out / variant
    vout.mkdir(parents=True, exist_ok=True)
    write_csv(vout / 'move_candidates_with_asymmetry_invariants.csv', rows)
    write_csv(vout / 'fan_asymmetry_invariants.csv', list(fan_inv.values()))
    write_csv(vout / 'applied_asymmetry_selected_moves.csv', applied)

    metrics = core.full_metrics(model, K, args.source)
    fan_values = [x['fan_asymmetry_invariant'] for x in fan_inv.values()]
    variant_summary = {
        'variant': variant,
        'max_level': args.max_level,
        'base_metrics': metrics,
        'move_counts': audit.get('move_counts', audit.get('candidate_counts', {})),
        'legal_count': len(legal),
        'A_gated_count': len(gated),
        'topology_effective_count': len(topology_effective),
        'A_gated_topology_effective_count': len(gated_topology_effective),
        'fan_asymmetry_mean': mean_value(fan_values),
        'fan_asymmetry_max': max(fan_values) if fan_values else 0.0,
        'best_A_any': compact_row(best_A),
        'best_A_gated': compact_row(best_A_gated),
        'best_A_topology_effective': compact_row(best_topology_effective_A),
        'best_delta_beta2': compact_row(best_delta_beta2),
        'applied': applied,
    }
    (vout / 'variant_summary.json').write_text(json.dumps(variant_summary, indent=2), encoding='utf-8')
    return variant_summary


def make_docs(summary: dict) -> Tuple[str, str, str]:
    rows = summary['variant_rows']
    def table() -> str:
        lines = ['| variant | A_gated | topo_eff | gated_topo_eff | fan_A_mean | best_A_class | best_A_delta_b2 | best_A_rank | after_b2 | harmonic |',
                 '|---|---:|---:|---:|---:|---|---:|---:|---:|---:|']
        for r in rows:
            lines.append(
                f"| {r['variant']} | {r['A_gated_count']} | {r['topology_effective_count']} | {r['A_gated_topology_effective_count']} | "
                f"{r['fan_asymmetry_mean']:.6g} | {r['best_A_any'].get('move_class','')} | {r['best_A_any'].get('delta_beta2','')} | "
                f"{r['best_A_any'].get('A_rank_legal','')} | {r['applied'][0].get('after_beta2','')} | {r['applied'][0].get('after_harmonic_ratio','')} |"
            )
        return '\n'.join(lines)
    summary_md = f"""# SUMMARY — Inter-fan transport from asymmetry invariants

This package tests the sharper Claude/CNNA hypothesis:

```text
not every face/child needs a complement;
only distant, transversely complementary face pairs with measurable directed birth/backreaction asymmetry should become complement-pairing candidates.
```

The ranking score is deliberately based on local/provenance data only: directed imbalance, fan sibling nonreciprocity, live-record aging, descendant-shell UV-tail difference, conductance asymmetry, and transverse complementarity. It does not use `delta_beta2` for ranking.

{table()}

Main interpretation: if `real_growth` selects a topology-effective handle/quotient at high asymmetry rank while strict symmetrized/no-backreaction controls lose that selection, then selective complement pairing is a plausible derived-growth gate. If controls still select it, the operation is still too externally offered.
"""
    results_md = f"""# RESULTS — Asymmetry-invariant complement pairing gate

## What was tested

The previous inter-fan test transported Z3/birth-order labels between neighboring parent-fans. That was too generic because symmetrized controls still carry labels. This test instead gates complement pairing by invariants that should collapse when sequential birth/backreaction is removed:

```text
A_fan =
  directed sibling imbalance
+ parent-line live-record aging gradient
+ descendant-shell-load difference
+ nonreciprocal backreaction residue
+ conductance update asymmetry
```

For a face pair to pass the asymmetry gate, it must also be a nonlocal transverse-complement pair. Shelling/cap moves are retained as controls but do not pass the pair gate.

## Comparative result

{table()}

## Critical reading

The strongest possible result would be:

```text
real_growth: high-ranked gated handle/quotient with delta_beta2 > 0
strict_symmetrized_control: A_gate collapses or the beta2-positive move falls in rank
no_backreaction: A_gate strongly reduced
```

If `historical_symmetrized_birth` still shows a strong signal, that indicates the older symmetrized control is not strict enough; it still preserves sequential residues in this lightweight growth core. For that reason the package includes both `historical_symmetrized_birth` and `strict_symmetrized_control`.

## Model limitation

A handle/quotient candidate is still an offered move. The test does not yet prove that the growth dynamics must execute it. It tests whether the real sequential provenance data select exactly the non-shelling move that opens beta2, without looking at beta2 during ranking.
"""
    audit_md = """# SOURCE_AUDIT_1_40

Integrated non-obstructed threads:

- Script 1/2: sequential birth, parent-line plus older-sibling environment, newborn backreaction, local transverse sibling offset.
- Script 12: shell-normalized inverse-square should be read as a response/backreaction kernel, not as a claim that shelling-ball geometry is sufficient.
- Scripts 33-35: the relevant local operator is the real plaquette commutator `K_abc=[A_ab,A_bc]`, not a synthetic face field.
- Script 40: immediate parent+three-children tetrahedral closure is an obstruction if it only creates local boundary filling.
- New correction: the missing move class is selective non-shelling complement-pairing of boundary faces, gated by invariant asymmetry rather than by raw birth labels.
"""
    return summary_md, results_md, audit_md


def package(out: Path, script_path: Path, zip_path: Path) -> None:
    if zip_path.exists():
        zip_path.unlink()
    with zipfile.ZipFile(zip_path, 'w', compression=zipfile.ZIP_DEFLATED) as z:
        z.write(script_path, script_path.name)
        z.write(Path(__file__).with_name('cnna_non_shelling_core.py'), 'cnna_non_shelling_core.py')
        for p in sorted(out.rglob('*')):
            if p.is_file():
                z.write(p, str(p.relative_to(out.parent)))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument('--max-level', type=int, default=2)
    ap.add_argument('--mode', default='linear', choices=['linear', 'log', 'saturating'])
    ap.add_argument('--source', default='live', choices=['record', 'live', 'full'])
    ap.add_argument('--transverse-amp', type=float, default=0.42)
    ap.add_argument('--max-boundary-faces', type=int, default=90)
    ap.add_argument('--max-single-vertices', type=int, default=12)
    ap.add_argument('--max-pair-candidates', type=int, default=2500)
    ap.add_argument('--max-rows', type=int, default=5000)
    ap.add_argument('--out', default='asymmetry_invariant_pairing_out_L2')
    ap.add_argument('--zip', default='cnna_interfan_transport_from_asymmetry_invariants_pkg_L2.zip')
    args = ap.parse_args()

    out = Path(args.out)
    if out.exists():
        shutil.rmtree(out)
    out.mkdir(parents=True)
    variants = ['real_growth', 'historical_symmetrized_birth', 'strict_symmetrized_control', 'no_backreaction']
    summaries = [run_variant(v, args, out) for v in variants]
    variant_rows = []
    for s in summaries:
        variant_rows.append({
            'variant': s['variant'],
            'A_gated_count': s['A_gated_count'],
            'topology_effective_count': s['topology_effective_count'],
            'A_gated_topology_effective_count': s['A_gated_topology_effective_count'],
            'fan_asymmetry_mean': s['fan_asymmetry_mean'],
            'fan_asymmetry_max': s['fan_asymmetry_max'],
            'best_A_any': s['best_A_any'],
            'best_A_gated': s['best_A_gated'],
            'best_A_topology_effective': s['best_A_topology_effective'],
            'best_delta_beta2': s['best_delta_beta2'],
            'applied': s['applied'],
        })
    summary = {'args': vars(args), 'variant_rows': variant_rows, 'variants': summaries}
    (out / 'comparative_summary.json').write_text(json.dumps(summary, indent=2), encoding='utf-8')
    flat_rows = []
    for s in summaries:
        for a in s['applied']:
            row = {'variant': s['variant']}
            row.update(a)
            flat_rows.append(row)
    write_csv(out / 'comparative_applied_asymmetry_selected_moves.csv', flat_rows)
    comp_rows = []
    for s in summaries:
        comp_rows.append({
            'variant': s['variant'],
            'A_gated_count': s['A_gated_count'],
            'topology_effective_count': s['topology_effective_count'],
            'A_gated_topology_effective_count': s['A_gated_topology_effective_count'],
            'fan_asymmetry_mean': s['fan_asymmetry_mean'],
            'fan_asymmetry_max': s['fan_asymmetry_max'],
            'best_A_move_class': s['best_A_any'].get('move_class', ''),
            'best_A_delta_beta2': s['best_A_any'].get('delta_beta2', ''),
            'best_A_rank': s['best_A_any'].get('A_rank_legal', ''),
            'best_A_response_rank': s['best_A_any'].get('response_rank_legal', ''),
            'best_A_A_invariant': s['best_A_any'].get('A_invariant', ''),
            'top_A_after_beta2': s['applied'][0].get('after_beta2', ''),
            'top_A_after_harmonic_ratio': s['applied'][0].get('after_harmonic_ratio', ''),
        })
    write_csv(out / 'comparative_asymmetry_summary.csv', comp_rows)
    summary_md, results_md, audit_md = make_docs(summary)
    (out / 'SUMMARY.md').write_text(summary_md, encoding='utf-8')
    (out / 'RESULTS.md').write_text(results_md, encoding='utf-8')
    (out / 'SOURCE_AUDIT_1_40.md').write_text(audit_md, encoding='utf-8')
    (out / 'README.md').write_text('Run: python test_interfan_transport_from_asymmetry_invariants.py\n', encoding='utf-8')
    package(out, Path(__file__), Path(args.zip))


if __name__ == '__main__':
    main()
