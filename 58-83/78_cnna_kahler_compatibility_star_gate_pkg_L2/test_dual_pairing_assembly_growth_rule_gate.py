#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import math
import shutil
import zipfile
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Set

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
Edge = Tuple[int, int]


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
        v = row.get(key, default)
        if v is None or v == '':
            return default
        if isinstance(v, str) and v.lower() in {'true','false'}:
            return 1.0 if v.lower() == 'true' else 0.0
        return float(v)
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


def parse_faces(row: dict) -> Tuple[Optional[Face], Optional[Face]]:
    return p62.parse_candidate_faces(row)


def face_edges(f: Face) -> Set[Edge]:
    a,b,c = f
    return {tuple(sorted((a,b))), tuple(sorted((a,c))), tuple(sorted((b,c)))}  # type: ignore


def support_context(row_a: dict, row_b: dict) -> dict:
    fa1, fa2 = parse_faces(row_a)
    fb1, fb2 = parse_faces(row_b)
    if fa1 is None or fa2 is None or fb1 is None or fb2 is None:
        return {'context': 'bad_parse', 'face_overlap':0, 'edge_overlap':0, 'vertex_overlap':0, 'is_connected_context':False, 'is_strong_context':False}
    faces_a = {fa1, fa2}; faces_b = {fb1, fb2}
    verts_a = set(fa1) | set(fa2); verts_b = set(fb1) | set(fb2)
    edges_a = face_edges(fa1) | face_edges(fa2); edges_b = face_edges(fb1) | face_edges(fb2)
    fo = len(faces_a & faces_b)
    eo = len(edges_a & edges_b)
    vo = len(verts_a & verts_b)
    ctx = 'disjoint'
    if fo:
        ctx = 'shared_face'
    elif eo:
        ctx = 'shared_edge'
    elif vo:
        ctx = 'shared_vertex'
    return {'context':ctx, 'face_overlap':fo, 'edge_overlap':eo, 'vertex_overlap':vo, 'is_connected_context':vo>0, 'is_strong_context':fo>0 or eo>0}


def allowed_row(row: dict, args: argparse.Namespace, used_faces: set[Face]) -> bool:
    if row.get('status') != 'ok' or not fbool(row.get('A_gate')):
        return False
    allowed = {'handle_candidate'}
    if args.allow_quotient:
        allowed.add('quotient_candidate')
    if row.get('move_class') not in allowed:
        return False
    fa, fb = parse_faces(row)
    if fa is None or fb is None:
        return False
    if not args.allow_reuse_faces and (fa in used_faces or fb in used_faces):
        return False
    return True


def eval_candidate(model, K, row: dict, args: argparse.Namespace, cascade_index: int, rank_A: int, kappa: bool = True) -> dict:
    rr = dict(row)
    rr['nonlinear_cascade_score'] = nl.nonlinear_score(rr, args, cascade_index)
    out = {
        'candidate_id': rr.get('candidate_id',''),
        'move_class': rr.get('move_class',''),
        'status': rr.get('status',''),
        'face_a': rr.get('face_a',''),
        'face_b': rr.get('face_b',''),
        'new_tets': rr.get('new_tets',''),
        'A_gate': rr.get('A_gate',''),
        'A_rank_score': rr.get('A_rank_score',''),
        'A_invariant': rr.get('A_invariant',''),
        'directed_imbalance': rr.get('directed_imbalance',''),
        'transverse_complementarity': rr.get('transverse_complementarity',''),
        'nonlinear_cascade_score': rr.get('nonlinear_cascade_score',''),
        'rank_A_score': rank_A,
        'delta_beta1_audit_only': rr.get('delta_beta1',''),
        'delta_beta2_audit_only': rr.get('delta_beta2',''),
        'delta_beta3_audit_only': rr.get('delta_beta3',''),
        'decision_used_delta_beta': False,
    }
    out.update(p62.candidate_alignment(model, K, rr, args))
    if kappa and out.get('candidate_eval_status') == 'ok':
        signed_id = fval(out, 'comm_signed_birth_over_abs', 0.0)
        flip = candidate_alignment_under_kappa(model, K, rr, args)
        signed_k = fval(flip, 'comm_signed_birth_over_abs', 0.0)
        out['kappa_comm_signed_birth_over_abs'] = signed_k
        out['kappa_signed_flip_abs'] = abs(signed_id + signed_k) / (abs(signed_id) + abs(signed_k) + EPS)
        out['kappa_signed_amp_min'] = min(abs(signed_id), abs(signed_k))
    else:
        out['kappa_comm_signed_birth_over_abs'] = ''
        out['kappa_signed_flip_abs'] = 999.0
        out['kappa_signed_amp_min'] = 0.0
    return out


def mirror_birth_orders(model) -> Dict[int, int]:
    saved: Dict[int, int] = {}
    for vid, node in model.nodes.items():
        saved[vid] = int(node.birth_order)
        if int(node.birth_order) == 1:
            node.birth_order = 3
        elif int(node.birth_order) == 3:
            node.birth_order = 1
        else:
            node.birth_order = int(node.birth_order)
    return saved


def restore_birth_orders(model, saved: Dict[int, int]) -> None:
    for vid, bo in saved.items():
        model.nodes[vid].birth_order = bo


def candidate_alignment_under_kappa(model, K, row: dict, args: argparse.Namespace) -> dict:
    saved = mirror_birth_orders(model)
    try:
        return p62.candidate_alignment(model, K, row, args)
    finally:
        restore_birth_orders(model, saved)


def qp_balance(row: dict) -> float:
    q = fval(row, 'Q_norm', 0.0)
    p = fval(row, 'P_norm', 0.0)
    return min(q, p) / (max(q, p) + EPS) if max(q,p) > EPS else 0.0


def row_key(row: dict) -> Tuple[str,str,str,str]:
    return (str(row.get('candidate_id','')), str(row.get('face_a','')), str(row.get('face_b','')), str(row.get('move_class','')))


def find_original_row(rows: List[dict], erow: dict) -> Optional[dict]:
    key = row_key(erow)
    for r in rows:
        if row_key(r) == key:
            rr = dict(r)
            rr['nonlinear_cascade_score'] = erow.get('nonlinear_cascade_score','')
            return rr
    # fallback by face/move only
    fkey = (str(erow.get('face_a','')), str(erow.get('face_b','')), str(erow.get('move_class','')))
    for r in rows:
        if (str(r.get('face_a','')), str(r.get('face_b','')), str(r.get('move_class',''))) == fkey:
            rr = dict(r); rr['nonlinear_cascade_score'] = erow.get('nonlinear_cascade_score','')
            return rr
    return None


def proxy_A_score(row: dict) -> float:
    # Pair A should provide a provenance/QP carrier, without using delta_beta/H2.
    return (
        1.00 * qp_balance(row)
        + 0.35 * min(1.0, fval(row, 'comm_abs_area', 0.0))
        + 0.25 * min(1.0, fval(row, 'directed_imbalance', 0.0))
        + 0.25 * min(1.0, fval(row, 'transverse_complementarity', 0.0))
        + 0.15 * min(1.0, fval(row, 'nonlinear_cascade_score', 0.0) / 10.0)
        + 0.10 * min(1.0, fval(row, 'A_rank_score', 0.0))
    )


def proxy_B_score(row: dict) -> float:
    # Pair B should supply C/J-local lock plus genuine kappa flip.
    c = fval(row, 'best_C_eigen_J_lock_mean_resid', 9.0)
    cm = fval(row, 'best_C_eigen_J_lock_max_resid', 9.0)
    flip = fval(row, 'kappa_signed_flip_abs', 9.0)
    amp = fval(row, 'kappa_signed_amp_min', 0.0)
    return (
        -1.50 * min(2.0, c)
        -0.75 * min(2.0, cm)
        -1.25 * min(2.0, flip)
        +0.45 * min(1.0, amp)
        +0.25 * qp_balance(row)
        +0.15 * min(1.0, fval(row, 'directed_imbalance', 0.0))
    )


def context_bonus(ctx: dict) -> float:
    return {'shared_face':0.45, 'shared_edge':0.30, 'shared_vertex':0.12, 'disjoint':-0.75, 'bad_parse':-2.0}.get(ctx['context'], -1.0)


def choose_pair_A(eval_rows: List[dict], args: argparse.Namespace) -> Optional[dict]:
    ok = [r for r in eval_rows if r.get('candidate_eval_status') == 'ok']
    if not ok:
        return None
    ok.sort(key=lambda r: (proxy_A_score(r), fval(r,'A_rank_score',0.0), str(r.get('face_a','')), str(r.get('face_b',''))), reverse=True)
    return ok[0]


def choose_pair_B(eval_rows: List[dict], pair_A: dict, args: argparse.Namespace) -> Optional[dict]:
    candidates = []
    for r in eval_rows:
        if r is pair_A or row_key(r) == row_key(pair_A):
            continue
        if r.get('candidate_eval_status') != 'ok':
            continue
        ctx = support_context(pair_A, r)
        if args.require_connected_assembly and not ctx['is_connected_context']:
            continue
        if args.require_strong_assembly_context and not ctx['is_strong_context']:
            continue
        score = proxy_B_score(r) + context_bonus(ctx)
        rr = dict(r)
        rr['_assembly_context'] = ctx
        rr['_assembly_score'] = score
        candidates.append(rr)
    if not candidates:
        return None
    candidates.sort(key=lambda r: (r['_assembly_score'], -fval(r,'best_C_eigen_J_lock_mean_resid',9.0), -fval(r,'kappa_signed_flip_abs',9.0)), reverse=True)
    return candidates[0]


def eval_legal_rows(model, K, rows: List[dict], args: argparse.Namespace, used_faces: set[Face], cascade_index: int) -> List[dict]:
    legal = [r for r in rows if allowed_row(r, args, used_faces)]
    legal_by_A = sorted(legal, key=lambda r: fval(r, 'A_rank_score', 0.0), reverse=True)
    if args.max_eval_candidates > 0:
        legal_by_A = legal_by_A[:args.max_eval_candidates]
    return [eval_candidate(model, K, r, args, cascade_index, rank, kappa=args.eval_kappa) for rank, r in enumerate(legal_by_A, start=1)]


def apply_logged_pair(model, K, original_row: dict, args, event_t: int, cascade_index: int, role: str) -> Tuple[object, dict, bool]:
    K2, log, applied = nl.apply_pair(model, K, original_row, args, event_t, cascade_index)
    log['assembly_role'] = role
    log['selection_rule'] = 'dual_pairing_assembly_growth_rule'
    return K2, log, applied


def add_used_faces_from(row: dict, used_faces: set[Face]) -> None:
    fa, fb = parse_faces(row)
    if fa is not None:
        used_faces.add(fa)
    if fb is not None:
        used_faces.add(fb)


def build_model(variant: str, args: argparse.Namespace):
    return nl.build_model(variant, args)


def build_dual_assembly_complex(model, args: argparse.Namespace, variant: str, out: Path):
    K = core.SimplicialComplex(f'{variant}_dual_pairing_assembly_growth')
    root_seeded = False
    birth_log: List[dict] = []
    pairing_log: List[dict] = []
    assembly_log: List[dict] = []
    candidate_rows: List[dict] = []
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
            # Step A scan on current K.
            rows_A = nl.scan_rows(model, K, args)
            eval_A = eval_legal_rows(model, K, rows_A, args, used_faces, cascade_index)
            for er in eval_A:
                er.update({'variant':variant, 'event_t':event_t, 'scan_id':scan_id, 'cascade_index':cascade_index, 'pre_or_post':'pre_A', 'role_selected':''})
                candidate_rows.append(er)
            pair_A_eval = choose_pair_A(eval_A, args)
            if pair_A_eval is None:
                break
            pair_A_orig = find_original_row(rows_A, pair_A_eval)
            if pair_A_orig is None:
                break
            if fval(pair_A_orig, 'nonlinear_cascade_score', 0.0) < args.min_nonlinear_score:
                break
            before_A = core.full_metrics(model, K, args.source)
            K1, logA, appliedA = apply_logged_pair(model, K, pair_A_orig, args, event_t, cascade_index, 'A_beta_QP_proxy')
            logA['proxy_A_score'] = proxy_A_score(pair_A_eval)
            logA['C_eigen_J_lock_mean_resid_eval'] = pair_A_eval.get('best_C_eigen_J_lock_mean_resid','')
            logA['kappa_signed_flip_abs_eval'] = pair_A_eval.get('kappa_signed_flip_abs','')
            pairing_log.append(logA)
            if not appliedA:
                assembly_log.append({'variant':variant, 'event_t':event_t, 'scan_id':scan_id, 'cascade_index':cascade_index, 'assembly_applied':False, 'failure':'A_apply_failed', 'A_reason':logA.get('apply_reason','')})
                break
            add_used_faces_from(pair_A_orig, used_faces)
            global_pair_count += 1
            K = K1

            # Step B rescan on K after A.  This is the dynamic part; B is not stale from the old complex.
            if global_pair_count >= args.max_auto_pairings:
                assembly_log.append({'variant':variant, 'event_t':event_t, 'scan_id':scan_id, 'cascade_index':cascade_index, 'assembly_applied':False, 'B_attempted':False, 'reason':'max_pair_count_after_A'})
                break
            rows_B = nl.scan_rows(model, K, args)
            eval_B = eval_legal_rows(model, K, rows_B, args, used_faces if not args.allow_B_reuse_A_faces else set(), cascade_index+1)
            # Use original A-eval as context, but B rows live on post-A K.
            for er in eval_B:
                er.update({'variant':variant, 'event_t':event_t, 'scan_id':scan_id, 'cascade_index':cascade_index, 'pre_or_post':'post_A_pre_B', 'role_selected':''})
                candidate_rows.append(er)
            pair_B_eval = choose_pair_B(eval_B, pair_A_eval, args)
            ctx = pair_B_eval.get('_assembly_context', {}) if pair_B_eval else {}
            pair_B_orig = find_original_row(rows_B, pair_B_eval) if pair_B_eval else None
            if pair_B_eval is None or pair_B_orig is None:
                after_A = core.full_metrics(model, K, args.source)
                assembly_log.append({
                    'variant':variant, 'event_t':event_t, 'scan_id':scan_id, 'cascade_index':cascade_index,
                    'assembly_applied':False, 'A_applied':True, 'B_attempted':False, 'failure':'no_B_candidate',
                    'A_face_a': pair_A_eval.get('face_a',''), 'A_face_b': pair_A_eval.get('face_b',''),
                    'after_A_beta2': after_A['beta2'], 'A_measured_delta_beta2': logA.get('measured_delta_beta2',''),
                })
                cascade_index += 1
                scan_id += 1
                if not args.cascade_rescan:
                    break
                continue
            before_B = core.full_metrics(model, K, args.source)
            K2, logB, appliedB = apply_logged_pair(model, K, pair_B_orig, args, event_t, cascade_index+1, 'B_C_kappa_context')
            logB['proxy_B_score'] = proxy_B_score(pair_B_eval)
            logB['C_eigen_J_lock_mean_resid_eval'] = pair_B_eval.get('best_C_eigen_J_lock_mean_resid','')
            logB['kappa_signed_flip_abs_eval'] = pair_B_eval.get('kappa_signed_flip_abs','')
            logB['assembly_context'] = ctx.get('context','')
            pairing_log.append(logB)
            if appliedB:
                add_used_faces_from(pair_B_orig, used_faces)
                global_pair_count += 1
                K = K2
            after = core.full_metrics(model, K, args.source)
            assembly_log.append({
                'variant':variant, 'event_t':event_t, 'scan_id':scan_id, 'cascade_index':cascade_index,
                'assembly_applied': bool(appliedB), 'A_applied': True, 'B_attempted': True, 'B_applied': bool(appliedB),
                'A_face_a': pair_A_eval.get('face_a',''), 'A_face_b': pair_A_eval.get('face_b',''),
                'B_face_a': pair_B_eval.get('face_a',''), 'B_face_b': pair_B_eval.get('face_b',''),
                'context': ctx.get('context',''), 'face_overlap': ctx.get('face_overlap',''), 'edge_overlap': ctx.get('edge_overlap',''), 'vertex_overlap': ctx.get('vertex_overlap',''),
                'A_proxy_score': proxy_A_score(pair_A_eval), 'A_QP_balance': qp_balance(pair_A_eval),
                'A_delta_beta2_audit_only': pair_A_eval.get('delta_beta2_audit_only',''), 'A_measured_delta_beta2': logA.get('measured_delta_beta2',''),
                'B_proxy_score': proxy_B_score(pair_B_eval), 'B_C_lock': pair_B_eval.get('best_C_eigen_J_lock_mean_resid',''), 'B_kappa_flip_abs': pair_B_eval.get('kappa_signed_flip_abs',''),
                'B_delta_beta2_audit_only': pair_B_eval.get('delta_beta2_audit_only',''), 'B_measured_delta_beta2': logB.get('measured_delta_beta2',''),
                'after_beta0': after['beta0'], 'after_beta1': after['beta1'], 'after_beta2': after['beta2'], 'after_beta3': after['beta3'],
                'decision_used_delta_beta': False,
            })
            scan_id += 1
            cascade_index += 2
            if not args.cascade_rescan:
                break
    return K, birth_log, pairing_log, assembly_log, candidate_rows


def summarize_selection(pairing_log: List[dict]) -> dict:
    vals_c = [fval(r, 'C_eigen_J_lock_mean_resid_eval', math.nan) for r in pairing_log if str(r.get('C_eigen_J_lock_mean_resid_eval','')).strip() != '']
    vals_k = [fval(r, 'kappa_signed_flip_abs_eval', math.nan) for r in pairing_log if str(r.get('kappa_signed_flip_abs_eval','')).strip() != '']
    return {
        'selected_C_eigen_J_lock_mean_avg': float(np.mean(vals_c)) if vals_c else 0.0,
        'selected_kappa_flip_abs_avg': float(np.mean(vals_k)) if vals_k else 0.0,
        'decision_used_delta_beta_any': any(fbool(r.get('decision_used_delta_beta')) for r in pairing_log),
        'measured_delta_beta2_sum': sum(int(fval(r, 'measured_delta_beta2', 0.0)) for r in pairing_log),
        'pair_A_count': sum(1 for r in pairing_log if r.get('assembly_role') == 'A_beta_QP_proxy' and fbool(r.get('applied'))),
        'pair_B_count': sum(1 for r in pairing_log if r.get('assembly_role') == 'B_C_kappa_context' and fbool(r.get('applied'))),
    }


def run_variant(variant: str, args: argparse.Namespace, out: Path) -> dict:
    model = build_model(variant, args)
    model.grow(args.max_level)
    baseline_K = core.build_dynamic_outward_ngf_complex(model)
    baseline_metrics = core.full_metrics(model, baseline_K, args.source)
    vout = out / variant
    vout.mkdir(parents=True, exist_ok=True)
    K, birth_log, pairing_log, assembly_log, candidate_rows = build_dual_assembly_complex(model, args, variant, vout)
    auto_metrics = core.full_metrics(model, K, args.source)
    dm, pair_rows, top_rows, three_rows = p56.directed_metrics(model, K, pairing_log, args)
    sm, signed_rows, signed_face_rows = p59.signed_quadrature_rows(model, K, pairing_log, args)
    am, align_pair_rows, align_candidate_rows, align_candidate_summary = p61.alignment_search_metrics(model, K, pairing_log, args)
    sel = summarize_selection(pairing_log)
    write_csv(vout / 'birth_geometry_log.csv', birth_log)
    write_csv(vout / 'dual_assembly_pairing_log.csv', pairing_log)
    write_csv(vout / 'dual_assembly_log.csv', assembly_log)
    write_csv(vout / 'candidate_eval_rows.csv', candidate_rows)
    write_csv(vout / 'directed_pair_rows.csv', pair_rows)
    write_csv(vout / 'signed_quadrature_rows.csv', signed_rows)
    write_csv(vout / 'alignment_pair_rows.csv', align_pair_rows)
    write_csv(vout / 'alignment_candidate_rows.csv', align_candidate_rows)
    write_csv(vout / 'alignment_candidate_summary.csv', align_candidate_summary)
    summary = {
        'variant': variant,
        'max_level': args.max_level,
        'phase_sign': args.phase_sign,
        'baseline_metrics': baseline_metrics,
        'auto_metrics': auto_metrics,
        'directed_metrics': dm,
        'signed_quadrature_metrics': sm,
        'alignment_metrics': am,
        'selection_metrics': sel,
        'automatic_pairings_applied': sum(1 for x in pairing_log if fbool(x.get('applied'))),
        'assemblies_applied': sum(1 for x in assembly_log if fbool(x.get('B_applied'))),
        'assemblies_attempted': len(assembly_log),
        'interpretation_flags': {
            'beta2_opened_vs_baseline': auto_metrics['beta2'] > baseline_metrics['beta2'],
            'Q_and_P_harmonic_positive': am['Q_even_harmonic_ratio'] > args.harmonic_positive_threshold and am['P_odd_harmonic_ratio'] > args.harmonic_positive_threshold,
            'pair_transport_harmonic_positive': dm['pair_transport_harmonic_ratio'] > args.harmonic_positive_threshold,
            'J_lock_below_threshold': am['best_per_pair_mean_J_lock_resid'] < args.lock_residual_threshold,
            'strict_decision_used_delta_beta': sel['decision_used_delta_beta_any'],
        },
    }
    (vout / 'variant_dual_assembly_growth_summary.json').write_text(json.dumps(summary, indent=2), encoding='utf-8')
    return summary


def write_comparative(out: Path, rows: List[dict]) -> None:
    flat = []
    for r in rows:
        a = r['auto_metrics']; dm = r['directed_metrics']; sm = r['signed_quadrature_metrics']; am = r['alignment_metrics']; sel = r['selection_metrics']
        flat.append({
            'variant': r['variant'],
            'beta0': a['beta0'], 'beta1': a['beta1'], 'beta2': a['beta2'], 'beta3': a['beta3'],
            'pairings': r['automatic_pairings_applied'],
            'assemblies_applied': r['assemblies_applied'],
            'pair_A_count': sel['pair_A_count'],
            'pair_B_count': sel['pair_B_count'],
            'pair_transport_harmonic_ratio': dm['pair_transport_harmonic_ratio'],
            'pair_kappa_orientation_ratio': dm['pair_kappa_orientation_ratio'],
            'pair_H_coherence': dm['pair_orientation_coherence'],
            'Q_harmonic_ratio': am['Q_even_harmonic_ratio'],
            'P_harmonic_ratio': am['P_odd_harmonic_ratio'],
            'best_global_candidate': am['best_global_candidate'],
            'best_per_pair_mean_J_lock_resid': am['best_per_pair_mean_J_lock_resid'],
            'best_per_pair_max_J_lock_resid': am['best_per_pair_max_J_lock_resid'],
            'signed_birth_over_abs_sum_ratio': sm['signed_birth_over_abs_sum_ratio'],
            'selected_C_eigen_J_lock_mean_avg': sel['selected_C_eigen_J_lock_mean_avg'],
            'selected_kappa_flip_abs_avg': sel['selected_kappa_flip_abs_avg'],
            'decision_used_delta_beta_any': sel['decision_used_delta_beta_any'],
            'measured_delta_beta2_sum': sel['measured_delta_beta2_sum'],
        })
    write_csv(out / 'comparative_dual_assembly_growth_summary.csv', flat)


def make_docs(summary: dict) -> Tuple[str,str,str,str]:
    rows = summary['variant_rows']
    lines = [
        '| variant | beta | pairs | assemblies | pair harm | Q harm | P harm | best J-lock | signed birth | selected C-lock | selected kappa-flip | used dBeta? |',
        '|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|',
    ]
    for r in rows:
        a = r['auto_metrics']; dm = r['directed_metrics']; am = r['alignment_metrics']; sm = r['signed_quadrature_metrics']; sel = r['selection_metrics']
        lines.append(f"| {r['variant']} | ({a['beta0']},{a['beta1']},{a['beta2']},{a['beta3']}) | {r['automatic_pairings_applied']} | {r['assemblies_applied']} | {dm['pair_transport_harmonic_ratio']:.6g} | {am['Q_even_harmonic_ratio']:.6g} | {am['P_odd_harmonic_ratio']:.6g} | {am['best_per_pair_mean_J_lock_resid']:.6g} | {sm['signed_birth_over_abs_sum_ratio']:.6g} | {sel['selected_C_eigen_J_lock_mean_avg']:.6g} | {sel['selected_kappa_flip_abs_avg']:.6g} | {sel['decision_used_delta_beta_any']} |")
    table = '\n'.join(lines)
    smd = f"""# SUMMARY — dual-pairing assembly growth rule gate

Model label:
CNNA growing primal simplicial complex with deterministic sequential provenance growth, directed antisymmetric birth-transport operators, and a dynamic two-step pair assembly rule.

This package tests the motif found by the previous audit:

```text
Pair A: provenance/QP carrier proxy, selected without delta-beta.
Pair B: C-lock/kappa context proxy, selected after Pair A by a rescan.
```

No decision uses delta_beta, H2, complex scalars, Hodge, positivity, or a physical adjoint.  The vertex operator uses the antisymmetric birth-transport term and no final sym(M).

{table}
"""
    rmd = f"""# RESULTS — dual-pairing assembly growth rule gate

## Comparative table

{table}

## Interpretation protocol

A constructive success would require:

```text
1. strict_sym remains zero/killed;
2. real_growth opens beta2;
3. Q/P and pair-transport harmonic channels remain positive;
4. J-lock improves or remains competitive;
5. signed orientation improves without phase/i/Hodge/* input;
6. decision_used_delta_beta_any remains false.
```

The rule is deliberately two-step: Pair B is chosen after applying Pair A and rescanning the updated complex.  This avoids pretending that a stale same-scan shared-face candidate is automatically still legal after Pair A.
"""
    audit = """# SOURCE AUDIT

Carried forward:

- Pair-property tradeoff gate showed single-pair locality is obstructed: beta2/QP, C-lock, and kappa-flip split across different candidates.
- Dual-pair assembly audit showed same-scan connected two-pair motifs exist.

This package turns that audit into a dynamic growth-rule test.

Anti-smuggling constraints:

- Pair A is selected using provenance/QP proxies, not delta_beta.
- Pair B is selected using C-lock/kappa/context proxies after a rescan, not delta_beta.
- delta_beta, H2 and harmonic diagnostics are audit-only after application.
- no i, no global J, no Hodge star, no positivity, no C*-norm, no final sym(M).

Limitations:

- This is L2 only.
- The kappa diagnostic is label-kappa on a fixed grown model, not a fully regrown reverse-sibling universe.
- If Pair B cannot be found after Pair A, that is an actual dynamic-legality obstruction, not a code failure.
"""
    readme = """# Dual-pairing assembly growth rule gate

Run:

```bash
python3 test_dual_pairing_assembly_growth_rule_gate.py
```

Outputs include per-variant candidate scans, assembly logs, pairing logs, comparative CSV/JSON, RESULTS.md, SUMMARY.md, SOURCE_AUDIT.md, and the ZIP package.
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
    ap.add_argument('--allow-B-reuse-A-faces', action='store_true', default=True)
    ap.add_argument('--allow-quotient', action='store_true')
    ap.add_argument('--max-boundary-faces', type=int, default=90)
    ap.add_argument('--max-single-vertices', type=int, default=12)
    ap.add_argument('--max-pair-candidates', type=int, default=2200)
    ap.add_argument('--max-rows', type=int, default=4400)
    ap.add_argument('--max-auto-pairings', type=int, default=4)
    ap.add_argument('--max-cascade-per-birth', type=int, default=4)
    ap.add_argument('--min-tets-before-pairing', type=int, default=4)
    ap.add_argument('--min-birth-time-before-pairing', type=int, default=4)
    ap.add_argument('--min-nonlinear-score', type=float, default=0.0)
    ap.add_argument('--keep-top-candidates', type=int, default=120)
    ap.add_argument('--keep-top-faces', type=int, default=80)
    ap.add_argument('--max-eval-candidates', type=int, default=0)
    ap.add_argument('--harmonic-positive-threshold', type=float, default=1e-4)
    ap.add_argument('--antisym-eta', type=float, default=1.0)
    ap.add_argument('--erase-phase-for-strict-sym', action='store_true', default=True)
    ap.add_argument('--eval-kappa', action='store_true', default=True)
    ap.add_argument('--lock-residual-threshold', type=float, default=0.20)
    ap.add_argument('--lock-max-threshold', type=float, default=0.30)
    ap.add_argument('--require-connected-assembly', action='store_true', default=True)
    ap.add_argument('--require-strong-assembly-context', action='store_true', default=False)
    ap.add_argument('--variants', nargs='*', default=['real_growth','strict_symmetrized_control','no_backreaction'])
    ap.add_argument('--phase-sign', type=int, default=1)
    ap.add_argument('--out', default='dual_pairing_assembly_growth_rule_out_L2')
    ap.add_argument('--zip', default='cnna_dual_pairing_assembly_growth_rule_gate_pkg_L2.zip')
    args = ap.parse_args()

    out = Path(args.out)
    if out.exists():
        shutil.rmtree(out)
    out.mkdir(parents=True)
    rows = [run_variant(v, args, out) for v in args.variants]
    summary = {'args': vars(args), 'variant_rows': rows}
    (out / 'comparative_summary.json').write_text(json.dumps(summary, indent=2), encoding='utf-8')
    write_comparative(out, rows)
    smd, rmd, audit, readme = make_docs(summary)
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
                'variant': r['variant'],
                'beta': [r['auto_metrics'][f'beta{i}'] for i in range(4)],
                'pairings': r['automatic_pairings_applied'],
                'assemblies': r['assemblies_applied'],
                'pair_harm': r['directed_metrics']['pair_transport_harmonic_ratio'],
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
