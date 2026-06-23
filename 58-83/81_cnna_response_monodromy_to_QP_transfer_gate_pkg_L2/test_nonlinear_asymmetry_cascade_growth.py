#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import math
import shutil
import zipfile
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np

import cnna_non_shelling_core as core
from test_interfan_transport_from_asymmetry_invariants import (
    StrictSymmetrizedGrowth,
    candidate_asymmetry_features,
    compact_row,
    descendant_shell_loads,
    fan_invariants,
    rerank_by_asymmetry,
    write_csv,
)
from test_growth_with_asymmetry_gated_complement_pairing import ordinary_outward_step

Face = Tuple[int, int, int]
EPS = 1e-12


class NonlinearConductanceGrowth(core.DynamicProvenanceGrowth):
    def __init__(self, *args, response_mode: str = 'linear', nonlinear_gamma: float = 1.65, nonlinear_threshold: float = 1.8, **kwargs):
        self.response_mode = response_mode
        self.nonlinear_gamma = nonlinear_gamma
        self.nonlinear_threshold = nonlinear_threshold
        mode = response_mode if response_mode in {'linear', 'log', 'saturating'} else 'linear'
        super().__init__(*args, mode=mode, **kwargs)

    def child_conductance_from_env(self, env_load: float) -> float:
        if self.response_mode in {'linear', 'log', 'saturating'}:
            return super().child_conductance_from_env(env_load)
        if self.response_mode == 'power_saturating':
            x = max(0.0, float(env_load))
            y = (x ** self.nonlinear_gamma) / (1.0 + x ** (self.nonlinear_gamma - 1.0))
            return self.base + self.alpha_env * y
        if self.response_mode == 'threshold_power':
            x = max(0.0, float(env_load) - self.nonlinear_threshold)
            return self.base + self.alpha_env * (x ** self.nonlinear_gamma)
        raise ValueError(self.response_mode)


class StrictSymmetrizedNonlinear(StrictSymmetrizedGrowth):
    def __init__(self, *args, response_mode: str = 'linear', nonlinear_gamma: float = 1.65, nonlinear_threshold: float = 1.8, **kwargs):
        self.response_mode = response_mode
        self.nonlinear_gamma = nonlinear_gamma
        self.nonlinear_threshold = nonlinear_threshold
        mode = response_mode if response_mode in {'linear', 'log', 'saturating'} else 'linear'
        super().__init__(*args, mode=mode, **kwargs)

    def child_conductance_from_env(self, env_load: float) -> float:
        if self.response_mode in {'linear', 'log', 'saturating'}:
            return core.DynamicProvenanceGrowth.child_conductance_from_env(self, env_load)
        if self.response_mode == 'power_saturating':
            x = max(0.0, float(env_load))
            y = (x ** self.nonlinear_gamma) / (1.0 + x ** (self.nonlinear_gamma - 1.0))
            return self.base + self.alpha_env * y
        if self.response_mode == 'threshold_power':
            x = max(0.0, float(env_load) - self.nonlinear_threshold)
            return self.base + self.alpha_env * (x ** self.nonlinear_gamma)
        raise ValueError(self.response_mode)


def build_model(variant: str, args: argparse.Namespace) -> core.DynamicProvenanceGrowth:
    kwargs = {
        'response_mode': args.response_mode,
        'nonlinear_gamma': args.nonlinear_gamma,
        'nonlinear_threshold': args.nonlinear_threshold,
        'transverse_amp': args.transverse_amp,
    }
    if variant == 'strict_symmetrized_control':
        return StrictSymmetrizedNonlinear(**kwargs)
    rule = {
        'real_growth': 'real_growth',
        'historical_symmetrized_birth': 'symmetrized_birth',
        'no_backreaction': 'no_backreaction',
    }[variant]
    return NonlinearConductanceGrowth(growth_rule=rule, **kwargs)


def enrich_and_rerank(model: core.DynamicProvenanceGrowth, rows: List[dict]) -> List[dict]:
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


def nonlinear_score(row: dict, args: argparse.Namespace, cascade_index: int) -> float:
    A = max(0.0, float(row.get('A_invariant', 0.0) or 0.0))
    trans = max(0.0, float(row.get('transverse_complementarity', 0.0) or 0.0))
    directed = max(0.0, float(row.get('directed_imbalance', 0.0) or 0.0))
    base = max(0.0, float(row.get('A_rank_score', 0.0) or 0.0))
    # Deliberately no delta_beta: nonlinear reinforcement is based on provenance asymmetry only.
    threshold = max(0.0, A - args.cascade_A_threshold)
    reinforcement = (1.0 + threshold) ** args.cascade_gamma
    transverse_boost = (1.0 + args.transverse_nonlinear_weight * trans) ** args.cascade_transverse_gamma
    directed_boost = 1.0 + args.directed_nonlinear_weight * directed
    fatigue = 1.0 / (1.0 + args.cascade_fatigue * max(0, cascade_index - 1))
    return base * reinforcement * transverse_boost * directed_boost * fatigue


def pick_nonlinear_pair(rows: List[dict], args: argparse.Namespace, cascade_index: int, used_faces: set[Face]) -> Optional[dict]:
    allowed_classes = {'handle_candidate'}
    if args.allow_quotient:
        allowed_classes.add('quotient_candidate')
    sub: List[dict] = []
    for r in rows:
        if r.get('status') != 'ok' or not r.get('A_gate'):
            continue
        if r.get('move_class') not in allowed_classes:
            continue
        try:
            fa = core.parse_face_string(str(r.get('face_a', '')))
            fb_txt = str(r.get('face_b', '')).split('perm=')[0].strip()
            fb = core.parse_face_string(fb_txt)
        except Exception:
            fa = fb = tuple()
        if not args.allow_reuse_faces and (fa in used_faces or fb in used_faces):
            continue
        rr = dict(r)
        rr['nonlinear_cascade_score'] = nonlinear_score(rr, args, cascade_index)
        sub.append(rr)
    sub.sort(key=lambda r: float(r.get('nonlinear_cascade_score', 0.0)), reverse=True)
    return sub[0] if sub else None


def scan_rows(model: core.DynamicProvenanceGrowth, K: core.SimplicialComplex, args: argparse.Namespace) -> List[dict]:
    rows, _audit = core.enumerate_moves(
        model,
        K,
        source=args.source,
        max_boundary_faces=args.max_boundary_faces,
        max_single_vertices=args.max_single_vertices,
        max_pair_candidates=args.max_pair_candidates,
        max_rows=args.max_rows,
    )
    return enrich_and_rerank(model, rows)


def apply_pair(model: core.DynamicProvenanceGrowth, K: core.SimplicialComplex, chosen: dict, args: argparse.Namespace, event_t: int, cascade_index: int) -> Tuple[core.SimplicialComplex, dict, bool]:
    before = core.full_metrics(model, K, args.source)
    L, reason, encoded = core.apply_candidate_row(K, chosen)
    log = compact_row(chosen)
    log.update({
        'event_t': event_t,
        'cascade_index': cascade_index,
        'nonlinear_cascade_score': chosen.get('nonlinear_cascade_score', ''),
        'apply_reason': reason,
        'encoded_move': encoded,
        'decision_used_delta_beta': False,
    })
    if L is None or reason != 'ok':
        log.update({'applied': False})
        return K, log, False
    after = core.full_metrics(model, L, args.source)
    log.update({
        'applied': True,
        'before_beta0': before['beta0'], 'before_beta1': before['beta1'], 'before_beta2': before['beta2'], 'before_beta3': before['beta3'],
        'after_beta0': after['beta0'], 'after_beta1': after['beta1'], 'after_beta2': after['beta2'], 'after_beta3': after['beta3'],
        'after_boundary_fraction': after['boundary_fraction'],
        'after_edge_link_cycle_fraction': after['edge_link_cycle_fraction'],
        'after_K_mean': after['K_mean'],
        'after_harmonic_ratio': after['harmonic_ratio'],
        'after_exact_residual_ratio': after['exact_residual_ratio'],
        'measured_delta_beta1': after['beta1'] - before['beta1'],
        'measured_delta_beta2': after['beta2'] - before['beta2'],
        'measured_delta_beta3': after['beta3'] - before['beta3'],
    })
    return L, log, True


def nonlinear_cascade_after_birth(model: core.DynamicProvenanceGrowth, K: core.SimplicialComplex, args: argparse.Namespace, event_t: int, used_faces: set[Face], global_pair_count: int) -> Tuple[core.SimplicialComplex, List[dict], List[dict], int]:
    logs: List[dict] = []
    first_scan_sample: List[dict] = []
    cascade_index = 1
    while cascade_index <= args.max_cascade_per_birth and global_pair_count < args.max_auto_pairings:
        if len(K.tets) < args.min_tets_before_pairing:
            break
        rows = scan_rows(model, K, args)
        if not first_scan_sample:
            top = sorted(rows, key=lambda r: float(r.get('A_rank_score', 0.0) or 0.0), reverse=True)[:args.keep_top_candidates]
            first_scan_sample.extend(compact_row(r) | {'nonlinear_cascade_score': nonlinear_score(r, args, cascade_index)} for r in top)
        chosen = pick_nonlinear_pair(rows, args, cascade_index, used_faces)
        if chosen is None:
            break
        # optional nonlinear trigger: require score large enough after first pair.
        if float(chosen.get('nonlinear_cascade_score', 0.0) or 0.0) < args.min_nonlinear_score:
            break
        K2, log, applied = apply_pair(model, K, chosen, args, event_t, cascade_index)
        logs.append(log)
        if not applied:
            break
        try:
            fa = core.parse_face_string(str(chosen.get('face_a', '')))
            fb_txt = str(chosen.get('face_b', '')).split('perm=')[0].strip()
            fb = core.parse_face_string(fb_txt)
            used_faces.add(fa); used_faces.add(fb)
        except Exception:
            pass
        K = K2
        global_pair_count += 1
        cascade_index += 1
        if not args.cascade_rescan:
            break
    return K, logs, first_scan_sample, global_pair_count


def build_nonlinear_auto_complex(model: core.DynamicProvenanceGrowth, args: argparse.Namespace, out: Path, variant: str):
    K = core.SimplicialComplex(f'{variant}_nonlinear_cascade')
    root_seeded = False
    birth_log: List[dict] = []
    pairing_log: List[dict] = []
    candidate_sample: List[dict] = []
    used_faces: set[Face] = set()
    global_pair_count = 0
    scans_triggered = 0
    for ev in sorted(model.birth_events, key=lambda x: int(x['t'])):
        root_seeded, added, encoded = ordinary_outward_step(model, K, ev, root_seeded)
        birth_log.append({
            't': int(ev['t']),
            'parent': int(ev['parent']),
            'child': int(ev['child']),
            'level': int(ev['level']),
            'ordinary_added': added,
            'ordinary_encoded': encoded,
            'pair_count_before': global_pair_count,
            'tet_count_after_ordinary': len(K.tets),
        })
        if not added:
            continue
        if int(ev['t']) < args.min_birth_time_before_pairing:
            continue
        if global_pair_count >= args.max_auto_pairings:
            continue
        K, logs, sample, global_pair_count = nonlinear_cascade_after_birth(model, K, args, int(ev['t']), used_faces, global_pair_count)
        if sample and not candidate_sample:
            candidate_sample.extend(sample)
        if logs:
            scans_triggered += 1
            pairing_log.extend(logs)
    return K, birth_log, pairing_log, candidate_sample, scans_triggered


def run_variant(variant: str, args: argparse.Namespace, out: Path) -> dict:
    model = build_model(variant, args)
    model.grow(args.max_level)
    baseline_K = core.build_dynamic_outward_ngf_complex(model)
    baseline_metrics = core.full_metrics(model, baseline_K, args.source)
    vout = out / variant
    vout.mkdir(parents=True, exist_ok=True)
    auto_K, birth_log, pairing_log, candidate_sample, scans_triggered = build_nonlinear_auto_complex(model, args, vout, variant)
    auto_metrics = core.full_metrics(model, auto_K, args.source)
    write_csv(vout / 'birth_geometry_log.csv', birth_log)
    write_csv(vout / 'nonlinear_pairing_cascade_log.csv', pairing_log)
    write_csv(vout / 'first_scan_top_candidates.csv', candidate_sample)
    summary = {
        'variant': variant,
        'max_level': args.max_level,
        'response_mode': args.response_mode,
        'nonlinear_gamma': args.nonlinear_gamma,
        'baseline_metrics': baseline_metrics,
        'auto_metrics': auto_metrics,
        'automatic_pairings_applied': sum(1 for x in pairing_log if x.get('applied')),
        'automatic_pairing_attempts_logged': len(pairing_log),
        'births_with_cascade_logs': scans_triggered,
        'cascade_log_compact': pairing_log,
        'ordinary_tets_baseline': len(baseline_K.tets),
        'auto_tets': len(auto_K.tets),
        'interpretation_flags': {
            'nonlinear_growth_opened_beta2': auto_metrics['beta2'] > baseline_metrics['beta2'],
            'nonlinear_growth_nonexact_K_positive': auto_metrics['harmonic_ratio'] > args.harmonic_positive_threshold,
            'strict_decision_used_delta_beta': False,
            'cascade_not_linear_level_scaling': True,
            'baseline_ball_like': baseline_metrics['beta1'] == 0 and baseline_metrics['beta2'] == 0,
        },
    }
    (vout / 'variant_summary.json').write_text(json.dumps(summary, indent=2), encoding='utf-8')
    return summary


def make_docs(summary: dict) -> Tuple[str, str, str, str]:
    rows = summary['variant_rows']
    def table() -> str:
        lines = ['| variant | response_mode | baseline beta | nonlinear beta | pairings | harmonic | opened beta2 | nonexact K |', '|---|---|---:|---:|---:|---:|---:|---:|']
        for r in rows:
            b = r['baseline_metrics']; a = r['auto_metrics']
            lines.append(f"| {r['variant']} | {r['response_mode']} | ({b['beta0']},{b['beta1']},{b['beta2']},{b['beta3']}) | ({a['beta0']},{a['beta1']},{a['beta2']},{a['beta3']}) | {r['automatic_pairings_applied']} | {a['harmonic_ratio']:.6g} | {r['interpretation_flags']['nonlinear_growth_opened_beta2']} | {r['interpretation_flags']['nonlinear_growth_nonexact_K_positive']} |")
        return '\n'.join(lines)
    summary_md = f"""# SUMMARY — nonlinear asymmetry-cascade growth

This package implements the correction that growth need not proceed linearly as "one birth, one local closure, next level".

The tested loop is event-driven and nonlinear:

```text
ordinary outward birth
-> local/nonlocal boundary scan
-> asymmetry-gated complement pairing
-> if a pairing fires, rescan the updated complex immediately
-> allow a bounded cascade before the next ordinary birth
```

The decision rule uses provenance/response invariants and a nonlinear reinforcement of the A-gate. It does not inspect delta beta when selecting a move.

{table()}

Interpretation: a positive result means the complex opens beta2 and the K-sector acquires a nonzero harmonic component during nonlinear growth itself, not after a post-hoc candidate application. A strict symmetry control must remain negative for the mechanism to be selective rather than merely an offered topological move.
"""
    results_md = f"""# RESULTS — nonlinear asymmetry-cascade growth

{table()}

## Nonlinear components

The package differs from the previous automatic L2 test in three ways:

1. A birth can trigger a cascade of complement-pairing moves before the next ordinary birth.
2. Candidate ranking uses a nonlinear score based on A-invariant excess, transverse complementarity, and directed imbalance.
3. The conductance response can be nonlinear via `power_saturating` or `threshold_power` modes.

No decision uses `delta_beta`. Topological changes are measured only after applying the selected move.

## Critical reading

If real growth opens beta2 but `strict_symmetrized_control` remains beta2=0, the mechanism is selective for sequential provenance asymmetry. If `no_backreaction` also opens beta2, the topological opening is driven by older-sibling/environment asymmetry rather than by backreaction alone. That is not a failure; it localizes the generative source.

## Raw compact summary

```json
{json.dumps(rows, indent=2)[:18000]}
```
"""
    audit_md = """# SOURCE_AUDIT_1_40

Threads carried forward:

- Script 1/2: sequential birth, parent-line + older siblings, directed sibling asymmetry, newborn as later UV-tail/backreaction for ancestors.
- Script 12: shell-normalized inverse-square is a response/backreaction weighting; it is not an argument for shelling-ball topology.
- Script 35: use real `K_abc=[A_ab,A_bc]`-style operator content inherited through the core metrics, not a synthetic face field.
- Script 40: immediate parent-fan tetrahedral closure remains a local Korand trap if no non-shelling complement pairing is allowed.
- New correction: growth can be event-driven/nonlinear; a birth may trigger a complement-pairing cascade rather than one linear local attachment.
"""
    readme = """# Nonlinear asymmetry-cascade growth package

Run default:

```bash
python3 test_nonlinear_asymmetry_cascade_growth.py --max-level 2 --response-mode power_saturating --max-auto-pairings 2 --max-cascade-per-birth 2
```

For a faster audit, reduce `--max-pair-candidates` and `--max-rows`. For scaling, use local-indexing before increasing level too aggressively.
"""
    return summary_md, results_md, audit_md, readme


def write_comparative(out: Path, summaries: List[dict]) -> None:
    rows = []
    for s in summaries:
        b = s['baseline_metrics']; a = s['auto_metrics']
        rows.append({
            'variant': s['variant'],
            'response_mode': s['response_mode'],
            'baseline_beta0': b['beta0'], 'baseline_beta1': b['beta1'], 'baseline_beta2': b['beta2'], 'baseline_beta3': b['beta3'],
            'nonlinear_beta0': a['beta0'], 'nonlinear_beta1': a['beta1'], 'nonlinear_beta2': a['beta2'], 'nonlinear_beta3': a['beta3'],
            'nonlinear_boundary_fraction': a['boundary_fraction'],
            'nonlinear_edge_link_cycle_fraction': a['edge_link_cycle_fraction'],
            'nonlinear_K_mean': a['K_mean'],
            'nonlinear_harmonic_ratio': a['harmonic_ratio'],
            'nonlinear_exact_residual_ratio': a['exact_residual_ratio'],
            'pairings_applied': s['automatic_pairings_applied'],
            'births_with_cascade_logs': s['births_with_cascade_logs'],
            'opened_beta2': s['interpretation_flags']['nonlinear_growth_opened_beta2'],
            'nonexact_K': s['interpretation_flags']['nonlinear_growth_nonexact_K_positive'],
        })
    write_csv(out / 'comparative_nonlinear_cascade_summary.csv', rows)


def package(out: Path, script_path: Path, zip_path: Path) -> None:
    if zip_path.exists():
        zip_path.unlink()
    with zipfile.ZipFile(zip_path, 'w', zipfile.ZIP_DEFLATED) as z:
        z.write(script_path, script_path.name)
        z.write(Path(__file__).with_name('cnna_non_shelling_core.py'), 'cnna_non_shelling_core.py')
        z.write(Path(__file__).with_name('test_interfan_transport_from_asymmetry_invariants.py'), 'test_interfan_transport_from_asymmetry_invariants.py')
        z.write(Path(__file__).with_name('test_growth_with_asymmetry_gated_complement_pairing.py'), 'test_growth_with_asymmetry_gated_complement_pairing.py')
        for p in out.rglob('*'):
            if p.is_file():
                z.write(p, str(p.relative_to(out.parent)))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument('--max-level', type=int, default=2)
    ap.add_argument('--response-mode', choices=['linear', 'log', 'saturating', 'power_saturating', 'threshold_power'], default='power_saturating')
    ap.add_argument('--source', default='live', choices=['record', 'live', 'full'])
    ap.add_argument('--transverse-amp', type=float, default=0.42)
    ap.add_argument('--nonlinear-gamma', type=float, default=1.65)
    ap.add_argument('--nonlinear-threshold', type=float, default=1.8)
    ap.add_argument('--cascade-A-threshold', type=float, default=0.18)
    ap.add_argument('--cascade-gamma', type=float, default=1.75)
    ap.add_argument('--cascade-transverse-gamma', type=float, default=1.25)
    ap.add_argument('--transverse-nonlinear-weight', type=float, default=1.4)
    ap.add_argument('--directed-nonlinear-weight', type=float, default=1.1)
    ap.add_argument('--cascade-fatigue', type=float, default=0.25)
    ap.add_argument('--cascade-rescan', action='store_true', default=True)
    ap.add_argument('--allow-reuse-faces', action='store_true')
    ap.add_argument('--allow-quotient', action='store_true')
    ap.add_argument('--max-boundary-faces', type=int, default=90)
    ap.add_argument('--max-single-vertices', type=int, default=12)
    ap.add_argument('--max-pair-candidates', type=int, default=2200)
    ap.add_argument('--max-rows', type=int, default=4400)
    ap.add_argument('--max-auto-pairings', type=int, default=2)
    ap.add_argument('--max-cascade-per-birth', type=int, default=2)
    ap.add_argument('--min-tets-before-pairing', type=int, default=4)
    ap.add_argument('--min-birth-time-before-pairing', type=int, default=4)
    ap.add_argument('--min-nonlinear-score', type=float, default=0.0)
    ap.add_argument('--keep-top-candidates', type=int, default=80)
    ap.add_argument('--harmonic-positive-threshold', type=float, default=1e-4)
    ap.add_argument('--variants', nargs='*', default=['real_growth','strict_symmetrized_control','no_backreaction'])
    ap.add_argument('--out', default='nonlinear_asymmetry_cascade_out_L2')
    ap.add_argument('--zip', default='cnna_nonlinear_asymmetry_cascade_growth_pkg_L2.zip')
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
                'nonlinear_beta': [s['auto_metrics'][f'beta{i}'] for i in range(4)],
                'pairings': s['automatic_pairings_applied'],
                'harmonic_ratio': s['auto_metrics']['harmonic_ratio'],
            } for s in summaries
        ]
    }, indent=2))


if __name__ == '__main__':
    main()
