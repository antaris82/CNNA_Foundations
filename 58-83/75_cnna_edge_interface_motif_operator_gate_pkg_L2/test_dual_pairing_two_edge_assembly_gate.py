#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import itertools
import json
import math
import shutil
import zipfile
from pathlib import Path
from typing import Dict, Iterable, List, Tuple, Set

EPS = 1e-12
Face = Tuple[int, int, int]
Edge = Tuple[int, int]


def read_csv(path: Path) -> List[dict]:
    if not path.exists() or path.stat().st_size == 0:
        return []
    with path.open(newline='', encoding='utf-8') as f:
        return list(csv.DictReader(f))


def write_csv(path: Path, rows: List[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text('', encoding='utf-8')
        return
    keys = sorted({k for r in rows for k in r.keys()})
    with path.open('w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        for r in rows:
            w.writerow(r)


def ffloat(x, default=0.0) -> float:
    try:
        if x is None or x == '':
            return default
        if isinstance(x, str) and x.lower() in {'true', 'false'}:
            return 1.0 if x.lower() == 'true' else 0.0
        return float(x)
    except Exception:
        return default


def fint(x, default=0) -> int:
    try:
        if x is None or x == '':
            return default
        return int(float(x))
    except Exception:
        return default


def fbool(x) -> bool:
    if isinstance(x, bool):
        return x
    if x is None:
        return False
    s = str(x).strip().lower()
    if s in {'true','1','yes','y'}:
        return True
    if s in {'false','0','no','n',''}:
        return False
    try:
        return float(s) != 0.0
    except Exception:
        return False


def parse_pair_key(key: str) -> Tuple[str, Face, Face]:
    parts = str(key).split('|')
    if len(parts) < 3:
        return str(key), (), ()  # type: ignore
    move = parts[0]
    def parse_face(txt: str) -> Face:
        vals = [int(v) for v in txt.replace('[',' ').replace(']',' ').replace(',',' ').split()]
        if len(vals) != 3:
            raise ValueError(f'bad face: {txt}')
        return tuple(sorted(vals))  # type: ignore
    return move, parse_face(parts[1]), parse_face(parts[2])


def face_edges(f: Face) -> Set[Edge]:
    a,b,c = f
    return {tuple(sorted((a,b))), tuple(sorted((a,c))), tuple(sorted((b,c)))}  # type: ignore


def candidate_geometry(row: dict) -> dict:
    move, f1, f2 = parse_pair_key(row.get('candidate_pair_key',''))
    faces = {f1, f2}
    verts = set(f1) | set(f2)
    edges = face_edges(f1) | face_edges(f2)
    return {'move': move, 'faces': faces, 'verts': verts, 'edges': edges, 'f1': f1, 'f2': f2}


def context_between(ga: dict, gb: dict) -> dict:
    face_overlap = len(ga['faces'] & gb['faces'])
    edge_overlap = len(ga['edges'] & gb['edges'])
    vertex_overlap = len(ga['verts'] & gb['verts'])
    disjoint = vertex_overlap == 0
    # Connected if candidate-pair supports share at least one vertex. Strong if share edge or face.
    context = 'disjoint'
    if face_overlap:
        context = 'shared_face'
    elif edge_overlap:
        context = 'shared_edge'
    elif vertex_overlap:
        context = 'shared_vertex'
    return {
        'face_overlap': face_overlap,
        'edge_overlap': edge_overlap,
        'vertex_overlap': vertex_overlap,
        'context': context,
        'is_connected_context': not disjoint,
        'is_strong_context': face_overlap > 0 or edge_overlap > 0,
    }


def add_role_flags(row: dict) -> dict:
    r = dict(row)
    r['_A_gate'] = fbool(r.get('passes_A_gate_both'))
    r['_beta2'] = fbool(r.get('passes_beta2_audit_both'))
    r['_QP'] = fbool(r.get('passes_QP_proxy'))
    r['_C'] = fbool(r.get('passes_C_lock_worst'))
    r['_K'] = fbool(r.get('passes_kappa_signed_flip'))
    r['_candidate_id'] = str(r.get('candidate_pair_key',''))
    r['_scan_id'] = str(r.get('scan_id',''))
    r['_event_t'] = str(r.get('event_t',''))
    r['_cascade_index'] = str(r.get('cascade_index',''))
    r['_C_worst'] = ffloat(r.get('C_lock_worst'), 999.0)
    r['_flip_abs'] = ffloat(r.get('signed_flip_abs'), 999.0)
    r['_QP_balance_min'] = ffloat(r.get('QP_balance_min'), 0.0)
    r['_comm_area'] = ffloat(r.get('comm_abs_area_avg'), 0.0)
    r['_signed_amp_min'] = ffloat(r.get('signed_amplitude_min'), 0.0)
    r['_beta_delta_min'] = fint(r.get('beta_delta_min'), 0)
    r['_directed'] = ffloat(r.get('directed_imbalance_avg'), 0.0)
    r['_transverse'] = ffloat(r.get('transverse_complementarity_avg'), 0.0)
    r['_transport_cosine'] = ffloat(r.get('transport_cosine_avg'), 0.0)
    r['_selected'] = fbool(r.get('selected_A_identity')) or fbool(r.get('selected_A_kappa'))
    # Roles for two-edge split. beta2 is audit-only; this script does not choose moves.
    r['_role_A_beta_QP'] = r['_A_gate'] and r['_beta2'] and r['_QP']
    r['_role_B_C_kappa'] = r['_A_gate'] and r['_C'] and r['_K']
    r['_role_B_C_only'] = r['_A_gate'] and r['_C']
    r['_role_B_kappa_only'] = r['_A_gate'] and r['_K']
    return r


def assembly_score(a: dict, b: dict, ctx: dict) -> float:
    # Smaller is better. Audit score only; not a growth decision.
    c = a['_C_worst'] if a['_C'] else b['_C_worst']
    flip = a['_flip_abs'] if a['_K'] else b['_flip_abs']
    qp = min(a['_QP_balance_min'], b['_QP_balance_min'])
    amp = max(a['_signed_amp_min'], b['_signed_amp_min'])
    context_bonus = {'shared_face': 0.0, 'shared_edge': 0.1, 'shared_vertex': 0.25, 'disjoint': 0.6}[ctx['context']]
    return float(c + flip + 0.5 * (1.0 - min(1.0, qp)) + 0.2 * (1.0 - min(1.0, amp)) + context_bonus)


def analyze_variant(variant: str, rows: List[dict], out: Path, args: argparse.Namespace) -> dict:
    vout = out / variant
    vout.mkdir(parents=True, exist_ok=True)
    rows = [add_role_flags(r) for r in rows]
    geoms = {i: candidate_geometry(r) for i, r in enumerate(rows)}

    # In a simultaneous assembly, require same scan_id by default.
    by_scan: Dict[str, List[int]] = {}
    for i, r in enumerate(rows):
        by_scan.setdefault(r['_scan_id'], []).append(i)

    assemblies: List[dict] = []
    strict_assemblies: List[dict] = []
    relaxed_assemblies: List[dict] = []
    roleA_count = sum(1 for r in rows if r['_role_A_beta_QP'])
    roleB_count = sum(1 for r in rows if r['_role_B_C_kappa'])

    for scan, idxs in by_scan.items():
        # Keep bounded but exhaustive for this L2 candidate set.
        for ia, ib in itertools.permutations(idxs, 2):
            if ia == ib:
                continue
            a, b = rows[ia], rows[ib]
            if not (a['_role_A_beta_QP'] and b['_role_B_C_kappa']):
                continue
            ctx = context_between(geoms[ia], geoms[ib])
            rec = {
                'variant': variant,
                'scan_id': scan,
                'event_t_A': a['_event_t'],
                'event_t_B': b['_event_t'],
                'candidate_A_key': a.get('candidate_pair_key',''),
                'candidate_B_key': b.get('candidate_pair_key',''),
                'A_beta_delta_min': a['_beta_delta_min'],
                'A_QP_balance_min': a['_QP_balance_min'],
                'A_directed_imbalance': a['_directed'],
                'A_transverse_complementarity': a['_transverse'],
                'A_transport_cosine': a['_transport_cosine'],
                'A_selected_by_old_rule': a['_selected'],
                'B_C_lock_worst': b['_C_worst'],
                'B_signed_flip_abs': b['_flip_abs'],
                'B_signed_amp_min': b['_signed_amp_min'],
                'B_QP_balance_min': b['_QP_balance_min'],
                'B_directed_imbalance': b['_directed'],
                'B_transverse_complementarity': b['_transverse'],
                'B_transport_cosine': b['_transport_cosine'],
                'B_selected_by_old_rule': b['_selected'],
                'context': ctx['context'],
                'face_overlap': ctx['face_overlap'],
                'edge_overlap': ctx['edge_overlap'],
                'vertex_overlap': ctx['vertex_overlap'],
                'is_connected_context': ctx['is_connected_context'],
                'is_strong_context': ctx['is_strong_context'],
            }
            rec['assembly_score_low_is_good'] = assembly_score(a, b, ctx)
            assemblies.append(rec)
            if ctx['is_connected_context']:
                relaxed_assemblies.append(rec)
            if ctx['is_strong_context']:
                strict_assemblies.append(rec)

    # Also audit partial assemblies: beta/QP role with C-only or kappa-only to see missing ingredient.
    partials: List[dict] = []
    for scan, idxs in by_scan.items():
        for ia, ib in itertools.permutations(idxs, 2):
            if ia == ib:
                continue
            a, b = rows[ia], rows[ib]
            if not a['_role_A_beta_QP']:
                continue
            if not (b['_role_B_C_only'] or b['_role_B_kappa_only']):
                continue
            ctx = context_between(geoms[ia], geoms[ib])
            if not ctx['is_connected_context']:
                continue
            partials.append({
                'variant': variant,
                'scan_id': scan,
                'candidate_A_key': a.get('candidate_pair_key',''),
                'candidate_B_key': b.get('candidate_pair_key',''),
                'B_C_only': b['_role_B_C_only'],
                'B_kappa_only': b['_role_B_kappa_only'],
                'B_C_and_kappa': b['_role_B_C_kappa'],
                'A_beta_delta_min': a['_beta_delta_min'],
                'A_QP_balance_min': a['_QP_balance_min'],
                'B_C_lock_worst': b['_C_worst'],
                'B_signed_flip_abs': b['_flip_abs'],
                'context': ctx['context'],
                'face_overlap': ctx['face_overlap'],
                'edge_overlap': ctx['edge_overlap'],
                'vertex_overlap': ctx['vertex_overlap'],
            })

    assemblies.sort(key=lambda r: r['assembly_score_low_is_good'])
    relaxed_assemblies.sort(key=lambda r: r['assembly_score_low_is_good'])
    strict_assemblies.sort(key=lambda r: r['assembly_score_low_is_good'])
    partials.sort(key=lambda r: (not r['B_C_and_kappa'], r['B_C_lock_worst'] + r['B_signed_flip_abs']))

    write_csv(vout / 'dual_pair_assemblies_all.csv', assemblies)
    write_csv(vout / 'dual_pair_assemblies_connected.csv', relaxed_assemblies)
    write_csv(vout / 'dual_pair_assemblies_strong_context.csv', strict_assemblies)
    write_csv(vout / 'partial_role_split_assemblies_connected.csv', partials[:args.keep_top])
    write_csv(vout / 'top_dual_pair_assemblies.csv', assemblies[:args.keep_top])

    context_counts = {}
    for r in assemblies:
        context_counts[r['context']] = context_counts.get(r['context'], 0) + 1

    # stronger gate: same scan, connected context, role A, role B.
    summary = {
        'variant': variant,
        'candidate_count': len(rows),
        'scan_count': len(by_scan),
        'role_A_beta_QP_count': roleA_count,
        'role_B_C_kappa_count': roleB_count,
        'role_A_and_role_B_same_single_pair_count': sum(1 for r in rows if r['_role_A_beta_QP'] and r['_role_B_C_kappa']),
        'dual_pair_assemblies_total_same_scan': len(assemblies),
        'dual_pair_assemblies_connected_context': len(relaxed_assemblies),
        'dual_pair_assemblies_strong_context': len(strict_assemblies),
        'assembly_context_counts': context_counts,
        'partial_connected_assemblies_count': len(partials),
        'best_assembly': assemblies[0] if assemblies else None,
        'best_connected_assembly': relaxed_assemblies[0] if relaxed_assemblies else None,
        'best_strong_context_assembly': strict_assemblies[0] if strict_assemblies else None,
        'decision_used_delta_beta_any': any(fbool(r.get('decision_used_delta_beta')) for r in rows),
    }
    (vout / 'variant_dual_pair_assembly_summary.json').write_text(json.dumps(summary, indent=2), encoding='utf-8')
    return summary


def discover_input_rows(input_root: Path, variants: List[str]) -> Dict[str, Path]:
    paths = {}
    # Prefer pair_property output.
    for v in variants:
        p = input_root / 'pair_property_tradeoff_obstruction_out_L2' / v / 'prepared_candidate_rows.csv'
        if p.exists():
            paths[v] = p
    # fallback if user points directly to output folder
    for v in variants:
        if v not in paths:
            p = input_root / v / 'prepared_candidate_rows.csv'
            if p.exists():
                paths[v] = p
    return paths


def make_docs(out: Path, summaries: List[dict]) -> None:
    table_lines = [
        '| variant | candidates | role A beta/QP | role B C/kappa | single-pair all | two-edge same-scan | connected | strong | used Δβ? |',
        '|---|---:|---:|---:|---:|---:|---:|---:|---:|'
    ]
    for s in summaries:
        table_lines.append(f"| {s['variant']} | {s['candidate_count']} | {s['role_A_beta_QP_count']} | {s['role_B_C_kappa_count']} | {s['role_A_and_role_B_same_single_pair_count']} | {s['dual_pair_assemblies_total_same_scan']} | {s['dual_pair_assemblies_connected_context']} | {s['dual_pair_assemblies_strong_context']} | {s['decision_used_delta_beta_any']} |")
    table = '\n'.join(table_lines)
    results = f"""# RESULTS — dual pairing two-edge assembly gate

## Comparative table

{table}

## Interpretation

This audit tests whether the roles that did not coincide on a single pair can be distributed across two same-scan pair candidates:

```text
Pair A: beta2/QP carrier
Pair B: C-lock/kappa-flip carrier
Assembly: same scan plus geometric context
```

The script does not introduce a new growth rule and does not use beta2 as a move decision.  `delta_beta2` is an audit label inherited from the candidate evaluation.

A positive two-edge signal requires nonzero same-scan assemblies; a stronger signal requires connected or strong face/edge context.  If connected assemblies exist, the single-pair obstruction is not the end of the line; the relevant object is at least a two-edge assembly.
"""
    summary = f"""# SUMMARY — dual pairing two-edge assembly gate

This package follows the tradeoff obstruction gate.  The previous result showed that beta2-opening, Q/P support, C-eigen lock and kappa-signed flip exist in the candidate space but do not coincide on one candidate.

This package checks whether they can be split across two coupled pairings in the same scan.

{table}

Conservative reading: a nonzero connected two-edge assembly count is not a derivation of `i`, `J`, `*`, positivity or a C*-structure.  It only means the single-pair-local obstruction was too strict and the next object is a multi-edge assembly.
"""
    audit = """# SOURCE AUDIT

Input source: `pair_property_tradeoff_obstruction_out_L2/*/prepared_candidate_rows.csv` from the previous package.

Methodological constraints:

- derived-only audit;
- no imported `i`, global `J`, Hodge star, positivity, norm axiom or physical adjoint;
- no new growth rule;
- no beta/H2/kappa target used as a move decision;
- beta2 is used only as an audit role label for Pair A.

This is not a complete new dynamic run.  It is a compatibility audit over an already generated L2 candidate space.
"""
    readme = """# Dual pairing two-edge assembly gate

Run:

```bash
python3 test_dual_pairing_two_edge_assembly_gate.py \
  --input-root . \
  --out dual_pairing_two_edge_assembly_out_L2 \
  --zip cnna_dual_pairing_two_edge_assembly_gate_pkg_L2.zip
```

The expected input folder is `pair_property_tradeoff_obstruction_out_L2` from the previous tradeoff package.
"""
    (out / 'RESULTS.md').write_text(results, encoding='utf-8')
    (out / 'SUMMARY.md').write_text(summary, encoding='utf-8')
    (out / 'SOURCE_AUDIT.md').write_text(audit, encoding='utf-8')
    (out / 'README.md').write_text(readme, encoding='utf-8')


def package(out: Path, zip_path: Path, extra_files: List[Path]) -> None:
    if zip_path.exists():
        zip_path.unlink()
    with zipfile.ZipFile(zip_path, 'w', zipfile.ZIP_DEFLATED) as z:
        z.write(Path(__file__), Path(__file__).name)
        for f in extra_files:
            if f.exists():
                z.write(f, f.name)
        for p in sorted(out.rglob('*')):
            if p.is_file():
                z.write(p, p.relative_to(out.parent))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument('--input-root', default='.')
    ap.add_argument('--out', default='dual_pairing_two_edge_assembly_out_L2')
    ap.add_argument('--zip', default='cnna_dual_pairing_two_edge_assembly_gate_pkg_L2.zip')
    ap.add_argument('--variants', nargs='*', default=['real_growth','strict_symmetrized_control','no_backreaction'])
    ap.add_argument('--keep-top', type=int, default=80)
    args = ap.parse_args()

    out = Path(args.out)
    if out.exists():
        shutil.rmtree(out)
    out.mkdir(parents=True)

    input_root = Path(args.input_root)
    paths = discover_input_rows(input_root, args.variants)
    summaries = []
    copied_inputs = out / 'input_prepared_candidate_rows'
    copied_inputs.mkdir(parents=True, exist_ok=True)
    for v in args.variants:
        p = paths.get(v)
        rows = read_csv(p) if p else []
        if p and p.exists():
            shutil.copy2(p, copied_inputs / f'{v}_prepared_candidate_rows.csv')
        summaries.append(analyze_variant(v, rows, out, args))

    (out / 'comparative_summary.json').write_text(json.dumps({'variant_rows': summaries, 'args': vars(args)}, indent=2), encoding='utf-8')
    flat = []
    for s in summaries:
        flat.append({k:v for k,v in s.items() if not isinstance(v, (dict, list)) and k not in {'best_assembly','best_connected_assembly','best_strong_context_assembly'}})
    write_csv(out / 'comparative_dual_pair_assembly_summary.csv', flat)
    # best assembly excerpt table
    excerpts = []
    for s in summaries:
        for key in ['best_assembly','best_connected_assembly','best_strong_context_assembly']:
            rec = s.get(key)
            if rec:
                q = dict(rec); q['assembly_kind'] = key; excerpts.append(q)
    write_csv(out / 'comparative_best_assembly_excerpts.csv', excerpts)
    make_docs(out, summaries)
    package(out, Path(args.zip), [])
    print(json.dumps({'zip': args.zip, 'out': args.out, 'summaries': summaries}, indent=2))


if __name__ == '__main__':
    main()
