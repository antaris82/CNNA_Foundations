#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import shutil
import zipfile
from pathlib import Path
from typing import List, Tuple, Optional

import numpy as np

import cnna_non_shelling_core as core
import test_nonlinear_asymmetry_cascade_growth as nl
from test_growth_with_asymmetry_gated_complement_pairing import ordinary_outward_step
import test_pairing_transport_antisym_birth_coherence_gate as p56
import test_signed_quadrature_area_kappa_gate as p59
import test_pair_J_alignment_search_gate as p61
import test_C_eigen_quadrature_refinement_gate as p62

EPS = 1e-12
Face = Tuple[int, int, int]


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


def fval(row: dict, key: str, default: float = 0.0) -> float:
    try:
        return float(row.get(key, default) or default)
    except Exception:
        return default


def parse_faces(row: dict) -> Tuple[Optional[Face], Optional[Face]]:
    return p62.parse_candidate_faces(row)


def faces_used(row: dict) -> Tuple[Optional[Face], Optional[Face]]:
    fa, fb = parse_faces(row)
    return fa, fb


def allowed_row(row: dict, args: argparse.Namespace, used_faces: set[Face]) -> bool:
    if row.get('status') != 'ok' or not bool(row.get('A_gate')):
        return False
    allowed_classes = {'handle_candidate'}
    if args.allow_quotient:
        allowed_classes.add('quotient_candidate')
    if row.get('move_class') not in allowed_classes:
        return False
    fa, fb = faces_used(row)
    if fa is None or fb is None:
        return False
    if not args.allow_reuse_faces and (fa in used_faces or fb in used_faces):
        return False
    return True


def eval_candidate(model, K, row: dict, args: argparse.Namespace, cascade_index: int, rank_A: int) -> dict:
    rr = dict(row)
    rr['nonlinear_cascade_score'] = nl.nonlinear_score(rr, args, cascade_index)
    out = {
        'candidate_id': rr.get('candidate_id', ''),
        'move_class': rr.get('move_class', ''),
        'face_a': rr.get('face_a', ''),
        'face_b': rr.get('face_b', ''),
        'A_gate': rr.get('A_gate', ''),
        'A_rank_score': rr.get('A_rank_score', ''),
        'A_invariant': rr.get('A_invariant', ''),
        'directed_imbalance': rr.get('directed_imbalance', ''),
        'transverse_complementarity': rr.get('transverse_complementarity', ''),
        'nonlinear_cascade_score': rr.get('nonlinear_cascade_score', ''),
        'rank_A_score': rank_A,
        'delta_beta1_audit_only': rr.get('delta_beta1', ''),
        'delta_beta2_audit_only': rr.get('delta_beta2', ''),
        'delta_beta3_audit_only': rr.get('delta_beta3', ''),
    }
    out.update(p62.candidate_alignment(model, K, rr, args))
    d2 = int(fval(rr, 'delta_beta2', 0.0))
    out['beta2_opening_audit_only'] = d2 > 0
    out['C_lock_lt_threshold'] = fval(out, 'best_C_eigen_J_lock_mean_resid', 1.0) < args.lock_residual_threshold
    out['C_lock_max_lt_threshold'] = fval(out, 'best_C_eigen_J_lock_max_resid', 1.0) < args.lock_max_threshold
    out['decision_used_delta_beta'] = False
    return out


def choose_current_rule(rows: List[dict], args: argparse.Namespace, cascade_index: int, used_faces: set[Face]) -> Optional[dict]:
    return nl.pick_nonlinear_pair(rows, args, cascade_index, used_faces)


def choose_C_guided(rows: List[dict], eval_rows: List[dict], args: argparse.Namespace) -> Optional[dict]:
    # Ranking is only among provenance A-gated legal candidates and uses C-eigen J-lock residual.
    # delta_beta* columns are audit-only and are deliberately not referenced here.
    ok = [r for r in eval_rows if r.get('candidate_eval_status') == 'ok']
    if not ok:
        return None
    ok.sort(key=lambda r: (
        fval(r, 'best_C_eigen_J_lock_mean_resid', 1.0),
        fval(r, 'best_C_eigen_J_lock_max_resid', 1.0),
        -fval(r, 'A_rank_score', 0.0),
        -fval(r, 'nonlinear_cascade_score', 0.0),
        str(r.get('face_a', '')),
        str(r.get('face_b', '')),
    ))
    best = ok[0]
    best_id = best.get('candidate_id', '')
    best_key = (best.get('face_a', ''), best.get('face_b', ''), best.get('move_class', ''))
    for row in rows:
        rr = dict(row)
        rr['nonlinear_cascade_score'] = nl.nonlinear_score(rr, args, int(best.get('cascade_index', 1) or 1))
        key = (rr.get('face_a', ''), rr.get('face_b', ''), rr.get('move_class', ''))
        if rr.get('candidate_id', '') == best_id or key == best_key:
            return rr
    return None


def scan_apply_rule(model, args: argparse.Namespace, variant: str, rule: str):
    K = core.SimplicialComplex(f'{variant}_{rule}_C_eigen_guided')
    root_seeded = False
    birth_log: List[dict] = []
    pairing_log: List[dict] = []
    candidate_rows: List[dict] = []
    scan_rows_log: List[dict] = []
    used_faces: set[Face] = set()
    global_pair_count = 0
    scan_id = 0

    for ev in sorted(model.birth_events, key=lambda x: int(x['t'])):
        root_seeded, added, encoded = ordinary_outward_step(model, K, ev, root_seeded)
        event_t = int(ev['t'])
        birth_log.append({
            't': event_t, 'parent': int(ev['parent']), 'child': int(ev['child']), 'level': int(ev['level']),
            'ordinary_added': added, 'ordinary_encoded': encoded, 'pair_count_before': global_pair_count,
            'tet_count_after_ordinary': len(K.tets),
        })
        if not added or event_t < args.min_birth_time_before_pairing or global_pair_count >= args.max_auto_pairings:
            continue
        cascade_index = 1
        while cascade_index <= args.max_cascade_per_birth and global_pair_count < args.max_auto_pairings:
            if len(K.tets) < args.min_tets_before_pairing:
                break
            rows = nl.scan_rows(model, K, args)
            legal = [r for r in rows if allowed_row(r, args, used_faces)]
            legal_by_A = sorted(legal, key=lambda r: fval(r, 'A_rank_score', 0.0), reverse=True)
            if args.max_eval_candidates > 0:
                legal_by_A = legal_by_A[:args.max_eval_candidates]
            eval_rows: List[dict] = []
            for rank, r in enumerate(legal_by_A, start=1):
                er = eval_candidate(model, K, r, args, cascade_index, rank)
                er.update({'variant': variant, 'rule': rule, 'phase_sign': args.phase_sign, 'event_t': event_t, 'scan_id': scan_id, 'cascade_index': cascade_index})
                eval_rows.append(er)
            if rule == 'A_rank_rule':
                chosen = choose_current_rule(rows, args, cascade_index, used_faces)
            elif rule == 'C_eigen_guided_rule':
                # Attach cascade_index to eval rows for deterministic score reconstruction in choose_C_guided.
                for er in eval_rows:
                    er['cascade_index'] = cascade_index
                chosen = choose_C_guided(rows, eval_rows, args)
            else:
                raise ValueError(f'unknown rule: {rule}')

            chosen_eval = None
            if chosen is not None:
                cfa, cfb = parse_faces(chosen)
                ckey = (str(chosen.get('candidate_id','')), str(chosen.get('face_a','')), str(chosen.get('face_b','')), str(chosen.get('move_class','')))
                for er in eval_rows:
                    ekey = (str(er.get('candidate_id','')), str(er.get('face_a','')), str(er.get('face_b','')), str(er.get('move_class','')))
                    if ekey == ckey:
                        er[f'selected_by_{rule}'] = True
                        chosen_eval = er
                    else:
                        er[f'selected_by_{rule}'] = False
            for er in eval_rows:
                candidate_rows.append(er)
            beta2_candidates = [r for r in eval_rows if r.get('beta2_opening_audit_only')]
            best_C = min([r for r in eval_rows if r.get('candidate_eval_status') == 'ok'], key=lambda r: (fval(r,'best_C_eigen_J_lock_mean_resid',1.0), fval(r,'best_C_eigen_J_lock_max_resid',1.0)), default=None)
            scan_rows_log.append({
                'variant': variant, 'rule': rule, 'phase_sign': args.phase_sign,
                'event_t': event_t, 'scan_id': scan_id, 'cascade_index': cascade_index,
                'legal_candidate_count': len(eval_rows),
                'beta2_opening_candidate_count_audit_only': len(beta2_candidates),
                'chosen_candidate_id': chosen.get('candidate_id','') if chosen else '',
                'chosen_face_a': chosen.get('face_a','') if chosen else '',
                'chosen_face_b': chosen.get('face_b','') if chosen else '',
                'chosen_C_lock_mean': chosen_eval.get('best_C_eigen_J_lock_mean_resid','') if chosen_eval else '',
                'chosen_C_lock_max': chosen_eval.get('best_C_eigen_J_lock_max_resid','') if chosen_eval else '',
                'chosen_delta_beta2_audit_only': chosen.get('delta_beta2','') if chosen else '',
                'best_C_candidate_id': best_C.get('candidate_id','') if best_C else '',
                'best_C_lock_mean': best_C.get('best_C_eigen_J_lock_mean_resid','') if best_C else '',
                'best_C_delta_beta2_audit_only': best_C.get('delta_beta2_audit_only','') if best_C else '',
            })
            scan_id += 1
            if chosen is None:
                break
            if fval(chosen, 'nonlinear_cascade_score', 0.0) < args.min_nonlinear_score:
                break
            K2, log, applied = nl.apply_pair(model, K, chosen, args, event_t, cascade_index)
            log['selection_rule'] = rule
            if chosen_eval:
                log['selected_C_eigen_J_lock_mean_resid'] = chosen_eval.get('best_C_eigen_J_lock_mean_resid','')
                log['selected_C_eigen_J_lock_max_resid'] = chosen_eval.get('best_C_eigen_J_lock_max_resid','')
                log['selected_best_C_eigen_candidate'] = chosen_eval.get('best_C_eigen_candidate','')
                log['selected_transport_cosine_ka_kb_reversed'] = chosen_eval.get('transport_cosine_ka_kb_reversed','')
                log['selected_comm_signed_birth_over_abs'] = chosen_eval.get('comm_signed_birth_over_abs','')
            pairing_log.append(log)
            if not applied:
                break
            fa, fb = parse_faces(chosen)
            if fa is not None and fb is not None:
                used_faces.add(fa); used_faces.add(fb)
            K = K2
            global_pair_count += 1
            cascade_index += 1
            if not args.cascade_rescan:
                break
    return K, birth_log, pairing_log, candidate_rows, scan_rows_log


def mean_vals(rows: List[dict], key: str) -> float:
    vals = []
    for r in rows:
        try:
            vals.append(float(r.get(key, 0.0) or 0.0))
        except Exception:
            pass
    return float(np.mean(vals)) if vals else 0.0


def selected_summary(pairing_log: List[dict]) -> dict:
    applied = [r for r in pairing_log if r.get('applied')]
    return {
        'selected_pair_count': len(applied),
        'selected_C_eigen_J_lock_mean_avg': mean_vals(applied, 'selected_C_eigen_J_lock_mean_resid'),
        'selected_C_eigen_J_lock_max_avg': mean_vals(applied, 'selected_C_eigen_J_lock_max_resid'),
        'selected_transport_cosine_avg': mean_vals(applied, 'selected_transport_cosine_ka_kb_reversed'),
        'selected_comm_signed_birth_over_abs_avg': mean_vals(applied, 'selected_comm_signed_birth_over_abs'),
        'decision_used_delta_beta_any': any(str(r.get('decision_used_delta_beta','')).lower() == 'true' for r in applied),
        'measured_delta_beta2_sum': sum(int(fval(r, 'measured_delta_beta2', 0.0)) for r in applied),
    }


def run_variant_rule(variant: str, phase_sign: int, rule: str, args: argparse.Namespace, out: Path) -> dict:
    local_args = argparse.Namespace(**vars(args))
    local_args.phase_sign = phase_sign
    vname = f'{variant}_{rule}_phase{phase_sign:+d}'.replace('+','plus').replace('-','minus')
    vout = out / vname
    vout.mkdir(parents=True, exist_ok=True)
    model = nl.build_model(variant, local_args)
    model.grow(local_args.max_level)
    baseline_K = core.build_dynamic_outward_ngf_complex(model)
    baseline_metrics = core.full_metrics(model, baseline_K, local_args.source)
    K, birth_log, pairing_log, candidate_rows, scan_rows_log = scan_apply_rule(model, local_args, variant, rule)
    auto_metrics = core.full_metrics(model, K, local_args.source)
    directed_metrics, directed_pair_rows, directed_face_rows, three_rows = p56.directed_metrics(model, K, pairing_log, local_args)
    signed_metrics, signed_pair_rows, signed_face_rows = p59.signed_quadrature_rows(model, K, pairing_log, local_args)
    alignment_metrics, pair_rows, candidate_summary, alignment_candidate_rows = p61.alignment_search_metrics(model, K, pairing_log, local_args)
    sel = selected_summary(pairing_log)
    write_csv(vout / 'birth_geometry_log.csv', birth_log)
    write_csv(vout / 'pairing_cascade_log.csv', pairing_log)
    write_csv(vout / 'candidate_eval_rows.csv', candidate_rows)
    write_csv(vout / 'scan_summaries.csv', scan_rows_log)
    write_csv(vout / 'directed_pair_rows.csv', directed_pair_rows)
    write_csv(vout / 'signed_quadrature_pair_rows.csv', signed_pair_rows)
    write_csv(vout / 'alignment_pair_rows.csv', pair_rows)
    write_csv(vout / 'alignment_candidate_summary.csv', candidate_summary)
    summary = {
        'variant': variant,
        'selection_rule': rule,
        'phase_sign': phase_sign,
        'variant_rule_phase': vname,
        'baseline_metrics': baseline_metrics,
        'auto_metrics': auto_metrics,
        'automatic_pairings_applied': sum(1 for r in pairing_log if r.get('applied')),
        'selection_metrics': sel,
        'directed_antisym_metrics': directed_metrics,
        'signed_quadrature_metrics': signed_metrics,
        'alignment_metrics': alignment_metrics,
        'interpretation_flags': {
            'beta2_opened': auto_metrics['beta2'] > baseline_metrics['beta2'],
            'pair_transport_harmonic_positive': directed_metrics['pair_transport_harmonic_ratio'] > args.harmonic_positive_threshold,
            'Q_harmonic_positive': alignment_metrics['Q_even_harmonic_ratio'] > args.harmonic_positive_threshold,
            'P_harmonic_positive': alignment_metrics['P_odd_harmonic_ratio'] > args.harmonic_positive_threshold,
            'C_eigen_lock_pass_mean': alignment_metrics['best_per_pair_mean_J_lock_resid'] < args.lock_residual_threshold,
            'C_eigen_lock_pass_max': alignment_metrics['best_per_pair_max_J_lock_resid'] < args.lock_max_threshold,
            'decision_used_delta_beta_any': sel['decision_used_delta_beta_any'],
        }
    }
    (vout / 'variant_rule_summary.json').write_text(json.dumps(summary, indent=2), encoding='utf-8')
    return summary


def compare_phase(rows: List[dict]) -> List[dict]:
    out = []
    keys = {}
    for r in rows:
        keys.setdefault((r['variant'], r['selection_rule']), {})[r['phase_sign']] = r
    for (variant, rule), d in sorted(keys.items()):
        if 1 not in d or -1 not in d:
            continue
        p = d[1]; m = d[-1]
        sp = p['signed_quadrature_metrics']['signed_birth_over_abs_sum_ratio']
        sm = m['signed_quadrature_metrics']['signed_birth_over_abs_sum_ratio']
        out.append({
            'variant': variant,
            'selection_rule': rule,
            'signed_birth_plus': sp,
            'signed_birth_minus': sm,
            'flip_score_zero_if_perfect': (sp + sm) / (abs(sp) + abs(sm) + EPS),
            'pair_transport_harm_plus': p['directed_antisym_metrics']['pair_transport_harmonic_ratio'],
            'pair_transport_harm_minus': m['directed_antisym_metrics']['pair_transport_harmonic_ratio'],
            'best_J_lock_plus': p['alignment_metrics']['best_per_pair_mean_J_lock_resid'],
            'best_J_lock_minus': m['alignment_metrics']['best_per_pair_mean_J_lock_resid'],
        })
    return out


def write_comparative(out: Path, rows: List[dict], phase_rows: List[dict]) -> None:
    flat = []
    for r in rows:
        a = r['auto_metrics']; dm = r['directed_antisym_metrics']; sm = r['signed_quadrature_metrics']; am = r['alignment_metrics']; sel = r['selection_metrics']
        flat.append({
            'variant_rule_phase': r['variant_rule_phase'],
            'variant': r['variant'],
            'selection_rule': r['selection_rule'],
            'phase_sign': r['phase_sign'],
            'beta0': a['beta0'], 'beta1': a['beta1'], 'beta2': a['beta2'], 'beta3': a['beta3'],
            'pairings': r['automatic_pairings_applied'],
            'selected_C_eigen_lock_mean_avg': sel['selected_C_eigen_J_lock_mean_avg'],
            'selected_C_eigen_lock_max_avg': sel['selected_C_eigen_J_lock_max_avg'],
            'selected_transport_cosine_avg': sel['selected_transport_cosine_avg'],
            'selected_comm_signed_birth_over_abs_avg': sel['selected_comm_signed_birth_over_abs_avg'],
            'pair_transport_harmonic_ratio': dm['pair_transport_harmonic_ratio'],
            'pair_kappa_orientation_ratio': dm['pair_kappa_orientation_ratio'],
            'pair_orientation_coherence': dm['pair_orientation_coherence'],
            'Q_harmonic_ratio': am['Q_even_harmonic_ratio'],
            'P_harmonic_ratio': am['P_odd_harmonic_ratio'],
            'best_global_candidate': am['best_global_candidate'],
            'best_per_pair_mean_J_lock_resid': am['best_per_pair_mean_J_lock_resid'],
            'best_per_pair_max_J_lock_resid': am['best_per_pair_max_J_lock_resid'],
            'signed_birth_over_abs_sum_ratio': sm['signed_birth_over_abs_sum_ratio'],
            'abs_area_mean_ratio': sm['abs_area_mean_ratio'],
            'decision_used_delta_beta_any': sel['decision_used_delta_beta_any'],
            'measured_delta_beta2_sum': sel['measured_delta_beta2_sum'],
        })
    write_csv(out / 'comparative_C_eigen_guided_pairing_summary.csv', flat)
    write_csv(out / 'phase_flip_summary.csv', phase_rows)


def make_docs(summary: dict, phase_rows: List[dict]) -> Tuple[str,str,str,str]:
    rows = summary['variant_rows']
    lines = [
        '| variant/rule/phase | beta | pairs | selected C-lock | pair harm | Q harm | P harm | best J-lock | signed birth | used dBeta? |',
        '|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|',
    ]
    for r in rows:
        a = r['auto_metrics']; dm = r['directed_antisym_metrics']; am = r['alignment_metrics']; sm = r['signed_quadrature_metrics']; sel = r['selection_metrics']
        lines.append(f"| {r['variant_rule_phase']} | ({a['beta0']},{a['beta1']},{a['beta2']},{a['beta3']}) | {r['automatic_pairings_applied']} | {sel['selected_C_eigen_J_lock_mean_avg']:.6g} | {dm['pair_transport_harmonic_ratio']:.6g} | {am['Q_even_harmonic_ratio']:.6g} | {am['P_odd_harmonic_ratio']:.6g} | {am['best_per_pair_mean_J_lock_resid']:.6g} | {sm['signed_birth_over_abs_sum_ratio']:.6g} | {sel['decision_used_delta_beta_any']} |")
    table = '\n'.join(lines)
    plines = ['| variant/rule | signed + | signed - | flip score | pair harm + | pair harm - | J-lock + | J-lock - |','|---|---:|---:|---:|---:|---:|---:|---:|']
    for p in phase_rows:
        plines.append(f"| {p['variant']} / {p['selection_rule']} | {p['signed_birth_plus']:.6g} | {p['signed_birth_minus']:.6g} | {p['flip_score_zero_if_perfect']:.6g} | {p['pair_transport_harm_plus']:.6g} | {p['pair_transport_harm_minus']:.6g} | {p['best_J_lock_plus']:.6g} | {p['best_J_lock_minus']:.6g} |")
    ptable = '\n'.join(plines)
    smd = f"""# SUMMARY — C-eigen guided pairing rule gate

Model label:
CNNA growing primal simplicial complex with deterministic sequential provenance growth, nonlinear asymmetry-gated complement-pair candidate space, directed antisymmetric birth-transport operators, and local C/J pair algebra.

This package compares two selection rules:

```text
A_rank_rule:
  inherited nonlinear cascade ranking.

C_eigen_guided_rule:
  select only from legal A-gated candidates,
  but rank by native C-eigen J-lock residual.
```

The guided rule does not use delta_beta/H2 as a decision input.  Delta-beta columns are audit-only after scan enumeration.

{table}

## Phase flip comparison

{ptable}
"""
    rmd = f"""# RESULTS — C-eigen guided pairing rule gate

## Comparative table

{table}

## Phase flip comparison

{ptable}

## Interpretation protocol

A constructive success would require:

```text
1. beta2 still opens,
2. pair_transport_harmonic_ratio remains positive,
3. Q and P harmonic channels remain positive,
4. C-eigen J-lock residual improves relative to A_rank_rule,
5. strict_sym stays killed,
6. decision_used_delta_beta_any remains false.
```

A strong J/i claim is still forbidden.  This test only checks whether the previous obstruction came from selecting the wrong pairings.
"""
    audit = """# SOURCE AUDIT

This package follows the refinement gate that found beta2-opening A-gated candidates with better C-eigen J-lock residuals than the currently selected pairs.

Anti-smuggling conditions:

- no complex scalars, no i;
- no Hodge star, no positivity, no physical adjoint;
- no final sym(M) in the directed birth-transport vertex operator;
- no arbitrary fitted rotation;
- no delta_beta/H2/kappa used in the C-eigen-guided selection rule.
"""
    readme = """# C-eigen guided pairing rule gate

Run:

```bash
python3 test_C_eigen_guided_pairing_rule_gate.py
```

Outputs include comparative CSV/JSON, per-rule candidate scans, RESULTS.md, SUMMARY.md, SOURCE_AUDIT.md, and the ZIP package.
"""
    return smd, rmd, audit, readme


def package(out: Path, zip_path: Path) -> None:
    files = [
        Path(__file__).name,
        'test_C_eigen_quadrature_refinement_gate.py',
        'test_pair_J_alignment_search_gate.py',
        'test_pairing_quadrature_adjoint_pairing_gate.py',
        'test_signed_quadrature_area_kappa_gate.py',
        'test_pairing_quadrature_split_symplectic_defect_gate.py',
        'test_pairing_transport_antisym_birth_coherence_gate.py',
        'test_pairing_transport_harmonic_kappa_gate.py',
        'test_nonlinear_asymmetry_cascade_growth.py',
        'test_harmonic_k_orientation_kappa_gate.py',
        'cnna_non_shelling_core.py',
        'test_interfan_transport_from_asymmetry_invariants.py',
        'test_growth_with_asymmetry_gated_complement_pairing.py',
    ]
    if zip_path.exists():
        zip_path.unlink()
    with zipfile.ZipFile(zip_path, 'w', zipfile.ZIP_DEFLATED) as z:
        for f in files:
            if Path(f).exists():
                z.write(f, f)
        for p in sorted(out.rglob('*')):
            if p.is_file():
                z.write(p, p.resolve().relative_to(Path.cwd()))


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
    ap.add_argument('--keep-top-candidates', type=int, default=120)
    ap.add_argument('--keep-top-faces', type=int, default=80)
    ap.add_argument('--max-eval-candidates', type=int, default=0, help='0 means evaluate all legal A-gated pair candidates')
    ap.add_argument('--harmonic-positive-threshold', type=float, default=1e-4)
    ap.add_argument('--antisym-eta', type=float, default=1.0)
    ap.add_argument('--erase-phase-for-strict-sym', action='store_true', default=True)
    ap.add_argument('--lock-residual-threshold', type=float, default=0.20)
    ap.add_argument('--lock-max-threshold', type=float, default=0.30)
    ap.add_argument('--variants', nargs='*', default=['real_growth','strict_symmetrized_control','no_backreaction'])
    ap.add_argument('--phase-signs', nargs='*', type=int, default=[1, -1])
    ap.add_argument('--rules', nargs='*', default=['A_rank_rule','C_eigen_guided_rule'])
    ap.add_argument('--out', default='C_eigen_guided_pairing_rule_out_L2')
    ap.add_argument('--zip', default='cnna_C_eigen_guided_pairing_rule_gate_pkg_L2.zip')
    args = ap.parse_args()

    out = Path(args.out)
    if out.exists():
        shutil.rmtree(out)
    out.mkdir(parents=True)
    rows = []
    for variant in args.variants:
        for phase in args.phase_signs:
            for rule in args.rules:
                rows.append(run_variant_rule(variant, int(phase), str(rule), args, out))
    phase_rows = compare_phase(rows)
    summary = {'args': vars(args), 'variant_rows': rows, 'phase_flip_rows': phase_rows}
    (out / 'comparative_summary.json').write_text(json.dumps(summary, indent=2), encoding='utf-8')
    write_comparative(out, rows, phase_rows)
    smd, rmd, audit, readme = make_docs(summary, phase_rows)
    (out / 'SUMMARY.md').write_text(smd, encoding='utf-8')
    (out / 'RESULTS.md').write_text(rmd, encoding='utf-8')
    (out / 'SOURCE_AUDIT.md').write_text(audit, encoding='utf-8')
    (out / 'README.md').write_text(readme, encoding='utf-8')
    package(out, Path(args.zip))
    print(json.dumps({
        'zip': args.zip,
        'out': args.out,
        'summary': [
            {
                'variant_rule_phase': r['variant_rule_phase'],
                'beta': [r['auto_metrics'][f'beta{i}'] for i in range(4)],
                'pairings': r['automatic_pairings_applied'],
                'selected_C_lock': r['selection_metrics']['selected_C_eigen_J_lock_mean_avg'],
                'pair_harm': r['directed_antisym_metrics']['pair_transport_harmonic_ratio'],
                'Q_harm': r['alignment_metrics']['Q_even_harmonic_ratio'],
                'P_harm': r['alignment_metrics']['P_odd_harmonic_ratio'],
                'best_J_lock': r['alignment_metrics']['best_per_pair_mean_J_lock_resid'],
                'signed_birth': r['signed_quadrature_metrics']['signed_birth_over_abs_sum_ratio'],
                'used_delta_beta': r['selection_metrics']['decision_used_delta_beta_any'],
            } for r in rows
        ]
    }, indent=2))


if __name__ == '__main__':
    main()
