#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import shutil
import zipfile
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import cnna_non_shelling_core as core
from test_interfan_transport_from_asymmetry_invariants import (
    StrictSymmetrizedGrowth,
    build_model,
    candidate_asymmetry_features,
    descendant_shell_loads,
    fan_invariants,
    rerank_by_asymmetry,
    compact_row,
    write_csv,
)

Face = Tuple[int, int, int]
Tet = Tuple[int, int, int, int]


def rebuild_face_maps(K: core.SimplicialComplex):
    occ: Dict[Face, int] = {}
    faces_by_vertex: Dict[int, set[Face]] = {}
    for tet in K.tets:
        for f in core.faces_of_tet(tet):
            occ[f] = occ.get(f, 0) + 1
            for v in f:
                faces_by_vertex.setdefault(v, set()).add(f)
    return faces_by_vertex, occ


def ordinary_outward_step(model: core.DynamicProvenanceGrowth, K: core.SimplicialComplex, ev: dict, root_seeded: bool):
    faces_by_vertex, occ = rebuild_face_maps(K)
    parent = int(ev['parent'])
    child = int(ev['child'])
    added = False
    encoded = ''
    # Seed the inner root fan as the first primal tetrahedral cell when the root fan is complete.
    if not root_seeded and len(model.nodes[model.root].children) == 3:
        ch = model.child_ids_ordered(model.root)
        tet = tuple(sorted((model.root, ch[0], ch[1], ch[2])))
        if K.add_tet(tet, birth_time=max(model.nodes[c].birth_time for c in ch)):
            added = True
            encoded = 'root_seed:' + ' '.join(map(str, tet))
        root_seeded = True
        faces_by_vertex, occ = rebuild_face_maps(K)
    if child in K.vertices:
        return root_seeded, added, encoded
    face = core.choose_boundary_face_for_parent(model, faces_by_vertex, occ, parent)
    if face is None:
        return root_seeded, added, encoded
    tet = tuple(sorted((*face, child)))
    if any(occ.get(ff, 0) >= 2 for ff in core.faces_of_tet(tet)):
        return root_seeded, added, encoded
    if K.add_tet(tet, birth_time=int(ev['t'])):
        added = True
        encoded = (encoded + ';' if encoded else '') + 'ordinary:' + ' '.join(map(str, tet))
    return root_seeded, added, encoded


def enrich_and_rerank(model: core.DynamicProvenanceGrowth, rows: List[dict], source: str) -> List[dict]:
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
    return rows


def pick_automatic_pairing(rows: List[dict], prefer: str, allow_quotient: bool) -> Optional[dict]:
    classes = []
    if prefer == 'handle_first':
        classes = ['handle_candidate'] + (['quotient_candidate'] if allow_quotient else [])
    elif prefer == 'quotient_first':
        classes = (['quotient_candidate'] if allow_quotient else []) + ['handle_candidate']
    else:
        classes = ['handle_candidate'] + (['quotient_candidate'] if allow_quotient else [])
    for cls in classes:
        sub = [r for r in rows if r.get('status') == 'ok' and r.get('A_gate') and r.get('move_class') == cls]
        sub.sort(key=lambda r: float(r.get('A_rank_score', 0.0)), reverse=True)
        if sub:
            return sub[0]
    return None


def scan_and_apply_pairing(
    model: core.DynamicProvenanceGrowth,
    K: core.SimplicialComplex,
    args: argparse.Namespace,
    event_t: int,
    pair_index: int,
) -> Tuple[core.SimplicialComplex, Optional[dict], List[dict]]:
    rows, audit = core.enumerate_moves(
        model,
        K,
        source=args.source,
        max_boundary_faces=args.max_boundary_faces,
        max_single_vertices=args.max_single_vertices,
        max_pair_candidates=args.max_pair_candidates,
        max_rows=args.max_rows,
    )
    rows = enrich_and_rerank(model, rows, args.source)
    chosen = pick_automatic_pairing(rows, args.prefer_pairing, args.allow_quotient)
    if chosen is None:
        return K, None, rows
    # Decision gate deliberately does not inspect delta_beta.  It uses only A_gate, move class,
    # manifold legality/status, and A_rank_score.  delta_* remains in the log for audit.
    L, reason, encoded = core.apply_candidate_row(K, chosen)
    if L is None or reason != 'ok':
        log = compact_row(chosen)
        log.update({'event_t': event_t, 'pair_index': pair_index, 'applied': False, 'apply_reason': reason, 'encoded_move': encoded})
        return K, log, rows
    before = core.full_metrics(model, K, args.source)
    after = core.full_metrics(model, L, args.source)
    log = compact_row(chosen)
    log.update({
        'event_t': event_t,
        'pair_index': pair_index,
        'applied': True,
        'apply_reason': reason,
        'encoded_move': encoded,
        'decision_used_delta_beta': False,
        'before_beta0': before['beta0'],
        'before_beta1': before['beta1'],
        'before_beta2': before['beta2'],
        'before_beta3': before['beta3'],
        'after_beta0': after['beta0'],
        'after_beta1': after['beta1'],
        'after_beta2': after['beta2'],
        'after_beta3': after['beta3'],
        'after_boundary_fraction': after['boundary_fraction'],
        'after_edge_link_cycle_fraction': after['edge_link_cycle_fraction'],
        'after_K_mean': after['K_mean'],
        'after_harmonic_ratio': after['harmonic_ratio'],
        'after_exact_residual_ratio': after['exact_residual_ratio'],
        'measured_delta_beta1': after['beta1'] - before['beta1'],
        'measured_delta_beta2': after['beta2'] - before['beta2'],
        'measured_delta_beta3': after['beta3'] - before['beta3'],
    })
    return L, log, rows


def build_auto_pairing_complex(model: core.DynamicProvenanceGrowth, args: argparse.Namespace, out: Path, variant: str):
    K = core.SimplicialComplex(f'{variant}_auto_asymmetry_pairing')
    root_seeded = False
    birth_log: List[dict] = []
    pairing_log: List[dict] = []
    scans_with_gate = 0
    scans = 0
    pair_count = 0
    candidate_sample: List[dict] = []
    for ev in sorted(model.birth_events, key=lambda x: int(x['t'])):
        root_seeded, added, encoded = ordinary_outward_step(model, K, ev, root_seeded)
        birth_log.append({
            't': int(ev['t']),
            'parent': int(ev['parent']),
            'child': int(ev['child']),
            'level': int(ev['level']),
            'ordinary_added': added,
            'ordinary_encoded': encoded,
            'pair_count_before': pair_count,
            'tet_count_after_ordinary': len(K.tets),
        })
        if not added:
            continue
        if pair_count >= args.max_auto_pairings:
            continue
        if len(K.tets) < args.min_tets_before_pairing:
            continue
        if int(ev['t']) < args.min_birth_time_before_pairing:
            continue
        K2, log, rows = scan_and_apply_pairing(model, K, args, int(ev['t']), pair_count + 1)
        scans += 1
        if any(r.get('A_gate') and r.get('status') == 'ok' for r in rows):
            scans_with_gate += 1
        if not candidate_sample:
            # keep first scan rows, compact enough for audit CSV
            for r in sorted(rows, key=lambda rr: float(rr.get('A_rank_score', 0.0)), reverse=True)[:min(80, len(rows))]:
                candidate_sample.append(compact_row(r))
        if log is not None:
            pairing_log.append(log)
            if log.get('applied'):
                K = K2
                pair_count += 1
    return K, birth_log, pairing_log, {'scans': scans, 'scans_with_gate': scans_with_gate, 'candidate_sample': candidate_sample}


def run_variant(variant: str, args: argparse.Namespace, out: Path) -> dict:
    model = build_model(variant, args)
    model.grow(args.max_level)
    baseline_K = core.build_dynamic_outward_ngf_complex(model)
    baseline_metrics = core.full_metrics(model, baseline_K, args.source)
    vout = out / variant
    vout.mkdir(parents=True, exist_ok=True)
    auto_K, birth_log, pairing_log, scan_info = build_auto_pairing_complex(model, args, vout, variant)
    auto_metrics = core.full_metrics(model, auto_K, args.source)
    write_csv(vout / 'birth_geometry_log.csv', birth_log)
    write_csv(vout / 'automatic_pairing_log.csv', pairing_log)
    write_csv(vout / 'first_scan_top_candidates.csv', scan_info['candidate_sample'])
    summary = {
        'variant': variant,
        'max_level': args.max_level,
        'baseline_metrics': baseline_metrics,
        'auto_metrics': auto_metrics,
        'ordinary_tets_baseline': len(baseline_K.tets),
        'auto_tets_or_quotient_tets': len(auto_K.tets),
        'automatic_pairings_applied': sum(1 for x in pairing_log if x.get('applied')),
        'automatic_pairing_attempts_logged': len(pairing_log),
        'scans': scan_info['scans'],
        'scans_with_A_gate': scan_info['scans_with_gate'],
        'pairing_log_compact': pairing_log,
        'interpretation_flags': {
            'auto_growth_opened_beta2': auto_metrics['beta2'] > baseline_metrics['beta2'],
            'auto_growth_nonexact_K_positive': auto_metrics['harmonic_ratio'] > args.harmonic_positive_threshold,
            'baseline_ball_like': baseline_metrics['beta1'] == 0 and baseline_metrics['beta2'] == 0,
            'pairing_was_automatic_not_posthoc': True,
            'decision_used_delta_beta': False,
        },
    }
    (vout / 'variant_summary.json').write_text(json.dumps(summary, indent=2), encoding='utf-8')
    return summary


def make_docs(summary: dict):
    rows = summary['variant_rows']
    def md_table():
        lines = ['| variant | baseline beta | auto beta | pairings | auto harmonic | opened beta2 | nonexact K |', '|---|---:|---:|---:|---:|---:|---:|']
        for r in rows:
            b = r['baseline_metrics']; a = r['auto_metrics']
            lines.append(f"| {r['variant']} | ({b['beta0']},{b['beta1']},{b['beta2']},{b['beta3']}) | ({a['beta0']},{a['beta1']},{a['beta2']},{a['beta3']}) | {r['automatic_pairings_applied']} | {a['harmonic_ratio']:.6g} | {r['interpretation_flags']['auto_growth_opened_beta2']} | {r['interpretation_flags']['auto_growth_nonexact_K_positive']} |")
        return '\n'.join(lines)
    result = md_table()
    summary_md = f"""# SUMMARY — growth with asymmetry-gated complement pairing

This package tests the previously missing step directly inside primal geometry growth:

```text
ordinary outward birth
-> boundary scan
-> if a nonlocal boundary-face pair passes the asymmetry-invariant complement gate
   and the resulting move is manifold-legal:
      apply the non-shelling complement pairing immediately
-> otherwise continue ordinary outward NGF/CQNM-like attachment
```

The decision gate uses A-invariant provenance/response quantities and transverse complementarity. It deliberately does **not** use delta beta in the decision. Delta beta is logged only after application for audit.

{result}

Interpretation: a positive result requires beta2 to open and the K-field to acquire a nonzero harmonic component during the automatic growth run itself. If beta2 opens only in real growth but not in strict symmetry control, the selective complement-pairing mechanism is no longer merely post-hoc candidate application. If strict controls also open beta2, the move class is still too generic.
"""
    results_md = f"""# RESULTS — automatic asymmetry-gated pairing

{result}

## Key decision rule

The automatic rule prefers handle candidates by default. Quotient candidates are available behind `--allow-quotient` because quotienting during later births requires a persistent vertex-map; without that map later birth events can reintroduce already identified vertices. The present run therefore tests the safer handle realization of complement pairing.

## Critical audit point

The automatic decision used:

```text
A_gate
move_class in handle_candidate (default)
status == ok
A_rank_score
```

It did not use `delta_beta1`, `delta_beta2`, or any after-the-fact topology gain to choose the move. Those fields are written only for post-run audit.

## Variant summaries

```json
{json.dumps(rows, indent=2)[:12000]}
```
"""
    audit = """# SOURCE_AUDIT_1_40

Non-obstructed threads carried forward:

- Script 1/2: sequential birth, parent-line plus older-sibling environment, directed backreaction, newborn as UV-tail/backreaction for ancestors, transverse sibling offset.
- Script 12: shell-normalized inverse-square response/backreaction weighting. This is not the same as ball-like shelling topology.
- Script 35: real operator sector uses K_abc=[A_ab,A_bc], not a synthetic K-field.
- Script 40 and later: local parent-fan tetrahedral closure is a Korand trap if treated as the full geometry.
- Recent move audits: non-shelling handle/quotient moves can open beta2, but prior tests applied them post-hoc. This package moves the pairing into the growth loop.
"""
    readme = """# Automatic asymmetry-gated complement pairing package

Run:

```bash
python3 test_growth_with_asymmetry_gated_complement_pairing.py --max-level 2 --max-auto-pairings 1 --out auto_pairing_out_L2
```

The default run is intentionally small and auditable. Increase `--max-level` and `--max-auto-pairings` only after inspecting the logs.
"""
    return summary_md, results_md, audit, readme


def write_comparative(out: Path, summaries: List[dict]):
    rows = []
    for s in summaries:
        b = s['baseline_metrics']; a = s['auto_metrics']
        rows.append({
            'variant': s['variant'],
            'baseline_beta0': b['beta0'], 'baseline_beta1': b['beta1'], 'baseline_beta2': b['beta2'], 'baseline_beta3': b['beta3'],
            'auto_beta0': a['beta0'], 'auto_beta1': a['beta1'], 'auto_beta2': a['beta2'], 'auto_beta3': a['beta3'],
            'auto_boundary_fraction': a['boundary_fraction'],
            'auto_edge_link_cycle_fraction': a['edge_link_cycle_fraction'],
            'auto_K_mean': a['K_mean'],
            'auto_harmonic_ratio': a['harmonic_ratio'],
            'auto_exact_residual_ratio': a['exact_residual_ratio'],
            'automatic_pairings_applied': s['automatic_pairings_applied'],
            'scans': s['scans'],
            'scans_with_A_gate': s['scans_with_A_gate'],
            'auto_growth_opened_beta2': s['interpretation_flags']['auto_growth_opened_beta2'],
            'auto_growth_nonexact_K_positive': s['interpretation_flags']['auto_growth_nonexact_K_positive'],
        })
    write_csv(out / 'comparative_auto_pairing_summary.csv', rows)


def package(out: Path, script: Path, zip_path: Path):
    with zipfile.ZipFile(zip_path, 'w', zipfile.ZIP_DEFLATED) as z:
        z.write(script, arcname=script.name)
        z.write(Path(__file__).with_name('cnna_non_shelling_core.py'), arcname='cnna_non_shelling_core.py')
        z.write(Path(__file__).with_name('test_interfan_transport_from_asymmetry_invariants.py'), arcname='test_interfan_transport_from_asymmetry_invariants.py')
        for p in out.rglob('*'):
            if p.is_file():
                z.write(p, arcname=str(p.relative_to(out.parent)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--max-level', type=int, default=2)
    ap.add_argument('--mode', default='linear')
    ap.add_argument('--source', default='live')
    ap.add_argument('--transverse-amp', type=float, default=0.42)
    ap.add_argument('--max-boundary-faces', type=int, default=90)
    ap.add_argument('--max-single-vertices', type=int, default=12)
    ap.add_argument('--max-pair-candidates', type=int, default=1800)
    ap.add_argument('--max-rows', type=int, default=3600)
    ap.add_argument('--max-auto-pairings', type=int, default=1)
    ap.add_argument('--min-tets-before-pairing', type=int, default=4)
    ap.add_argument('--min-birth-time-before-pairing', type=int, default=4)
    ap.add_argument('--prefer-pairing', choices=['handle_first','quotient_first'], default='handle_first')
    ap.add_argument('--allow-quotient', action='store_true')
    ap.add_argument('--harmonic-positive-threshold', type=float, default=1e-4)
    ap.add_argument('--variants', nargs='*', default=['real_growth','historical_symmetrized_birth','strict_symmetrized_control','no_backreaction'])
    ap.add_argument('--out', default='auto_asymmetry_pairing_out_L2')
    ap.add_argument('--zip', default='cnna_growth_with_asymmetry_gated_complement_pairing_pkg_L2.zip')
    args = ap.parse_args()
    out = Path(args.out)
    if out.exists():
        shutil.rmtree(out)
    out.mkdir(parents=True)
    summaries = [run_variant(v, args, out) for v in args.variants]
    summary = {'args': vars(args), 'variant_rows': summaries}
    (out / 'comparative_summary.json').write_text(json.dumps(summary, indent=2), encoding='utf-8')
    write_comparative(out, summaries)
    smd, rmd, amd, readme = make_docs(summary)
    (out / 'SUMMARY.md').write_text(smd, encoding='utf-8')
    (out / 'RESULTS.md').write_text(rmd, encoding='utf-8')
    (out / 'SOURCE_AUDIT_1_40.md').write_text(amd, encoding='utf-8')
    (out / 'README.md').write_text(readme, encoding='utf-8')
    package(out, Path(__file__), Path(args.zip))
    print(json.dumps({
        'zip': args.zip,
        'out': args.out,
        'summary': [
            {
                'variant': s['variant'],
                'baseline_beta': [s['baseline_metrics'][f'beta{i}'] for i in range(4)],
                'auto_beta': [s['auto_metrics'][f'beta{i}'] for i in range(4)],
                'pairings': s['automatic_pairings_applied'],
                'harmonic_ratio': s['auto_metrics']['harmonic_ratio'],
            } for s in summaries
        ]
    }, indent=2))


if __name__ == '__main__':
    main()
