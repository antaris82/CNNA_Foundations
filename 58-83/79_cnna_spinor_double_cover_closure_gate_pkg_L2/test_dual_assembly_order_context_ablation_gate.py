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
import test_dual_pairing_assembly_growth_rule_gate as p68

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
    return p68.fval(row, key, default)


def fbool(x) -> bool:
    return p68.fbool(x)


def parse_faces(row: dict):
    return p68.parse_faces(row)


def qp_balance(row: dict) -> float:
    return p68.qp_balance(row)


def row_key(row: dict):
    return p68.row_key(row)


def add_used_faces_from(row: dict, used_faces: set[Face]) -> None:
    return p68.add_used_faces_from(row, used_faces)


def context_bonus(ctx: dict) -> float:
    return p68.context_bonus(ctx)


def choose_pair_B_alone(eval_rows: List[dict], args: argparse.Namespace) -> Optional[dict]:
    ok = [r for r in eval_rows if r.get('candidate_eval_status') == 'ok']
    if not ok:
        return None
    ok.sort(key=lambda r: (
        p68.proxy_B_score(r),
        -fval(r, 'best_C_eigen_J_lock_mean_resid', 9.0),
        -fval(r, 'kappa_signed_flip_abs', 9.0),
        qp_balance(r),
        str(r.get('face_a','')),
        str(r.get('face_b','')),
    ), reverse=True)
    return ok[0]


def choose_pair_A_with_context(eval_rows: List[dict], pair_B: dict, args: argparse.Namespace) -> Optional[dict]:
    candidates = []
    for r in eval_rows:
        if row_key(r) == row_key(pair_B):
            continue
        if r.get('candidate_eval_status') != 'ok':
            continue
        ctx = p68.support_context(pair_B, r)
        if args.require_connected_assembly and not ctx['is_connected_context']:
            continue
        if args.require_strong_assembly_context and not ctx['is_strong_context']:
            continue
        rr = dict(r)
        rr['_assembly_context'] = ctx
        rr['_assembly_score'] = p68.proxy_A_score(r) + context_bonus(ctx)
        candidates.append(rr)
    if not candidates:
        return None
    candidates.sort(key=lambda r: (
        r['_assembly_score'],
        p68.proxy_A_score(r),
        fval(r,'A_rank_score',0.0),
        str(r.get('face_a','')),
        str(r.get('face_b','')),
    ), reverse=True)
    return candidates[0]


def choose_stale_AB(eval_rows: List[dict], args: argparse.Namespace) -> Tuple[Optional[dict], Optional[dict]]:
    A = p68.choose_pair_A(eval_rows, args)
    if A is None:
        return None, None
    B = p68.choose_pair_B(eval_rows, A, args)
    return A, B


def original_or_none(rows: List[dict], erow: Optional[dict]) -> Optional[dict]:
    if erow is None:
        return None
    return p68.find_original_row(rows, erow)


def apply_pair(model, K, orig: dict, args, event_t: int, cascade_index: int, role: str, rule: str):
    K2, log, applied = nl.apply_pair(model, K, orig, args, event_t, cascade_index)
    log['assembly_role'] = role
    log['selection_rule'] = rule
    return K2, log, applied


def eval_rows_for_scan(model, K, rows: List[dict], args, used_faces: set[Face], cascade_index: int) -> List[dict]:
    return p68.eval_legal_rows(model, K, rows, args, used_faces, cascade_index)


def option_tag(args: argparse.Namespace) -> str:
    ctx = 'strong' if args.require_strong_assembly_context else ('connected' if args.require_connected_assembly else 'anyctx')
    reuse = 'reuseB' if args.allow_B_reuse_A_faces else 'noReuseB'
    return f"{args.assembly_order}_{ctx}_{reuse}"


def build_ablation_complex(model, args: argparse.Namespace, variant: str, out: Path):
    K = core.SimplicialComplex(f'{variant}_{option_tag(args)}')
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

            rows0 = nl.scan_rows(model, K, args)
            eval0 = eval_rows_for_scan(model, K, rows0, args, used_faces, cascade_index)
            for er in eval0:
                er.update({'variant':variant, 'option':option_tag(args), 'event_t':event_t, 'scan_id':scan_id, 'cascade_index':cascade_index, 'stage':'initial'})
                candidate_rows.append(er)

            if args.assembly_order == 'A_to_B_rescan':
                A_eval = p68.choose_pair_A(eval0, args)
                A_orig = original_or_none(rows0, A_eval)
                if A_eval is None or A_orig is None:
                    break
                K1, logA, appliedA = apply_pair(model, K, A_orig, args, event_t, cascade_index, 'A_beta_QP_proxy', option_tag(args))
                logA['proxy_A_score'] = p68.proxy_A_score(A_eval)
                logA['C_eigen_J_lock_mean_resid_eval'] = A_eval.get('best_C_eigen_J_lock_mean_resid','')
                logA['kappa_signed_flip_abs_eval'] = A_eval.get('kappa_signed_flip_abs','')
                pairing_log.append(logA)
                if not appliedA:
                    assembly_log.append({'variant':variant, 'option':option_tag(args), 'event_t':event_t, 'scan_id':scan_id, 'assembly_applied':False, 'failure':'A_apply_failed'})
                    break
                add_used_faces_from(A_orig, used_faces)
                global_pair_count += 1
                K = K1
                if global_pair_count >= args.max_auto_pairings:
                    break

                rows1 = nl.scan_rows(model, K, args)
                used_for_B = set() if args.allow_B_reuse_A_faces else used_faces
                eval1 = eval_rows_for_scan(model, K, rows1, args, used_for_B, cascade_index+1)
                for er in eval1:
                    er.update({'variant':variant, 'option':option_tag(args), 'event_t':event_t, 'scan_id':scan_id, 'cascade_index':cascade_index+1, 'stage':'post_A'})
                    candidate_rows.append(er)
                B_eval = p68.choose_pair_B(eval1, A_eval, args)
                B_orig = original_or_none(rows1, B_eval)
                ctx = B_eval.get('_assembly_context',{}) if B_eval else {}
                if B_eval is None or B_orig is None:
                    assembly_log.append({'variant':variant, 'option':option_tag(args), 'event_t':event_t, 'scan_id':scan_id, 'assembly_applied':False, 'A_applied':True, 'B_attempted':False, 'failure':'no_B_candidate'})
                    cascade_index += 1; scan_id += 1
                    if not args.cascade_rescan: break
                    continue
                K2, logB, appliedB = apply_pair(model, K, B_orig, args, event_t, cascade_index+1, 'B_C_kappa_context', option_tag(args))
                logB['proxy_B_score'] = p68.proxy_B_score(B_eval)
                logB['C_eigen_J_lock_mean_resid_eval'] = B_eval.get('best_C_eigen_J_lock_mean_resid','')
                logB['kappa_signed_flip_abs_eval'] = B_eval.get('kappa_signed_flip_abs','')
                logB['assembly_context'] = ctx.get('context','')
                pairing_log.append(logB)
                if appliedB:
                    add_used_faces_from(B_orig, used_faces)
                    global_pair_count += 1
                    K = K2
                after = core.full_metrics(model, K, args.source)
                assembly_log.append(make_assembly_record(variant,args,event_t,scan_id,cascade_index,A_eval,B_eval,ctx,appliedB,after))
                cascade_index += 2; scan_id += 1
                if not args.cascade_rescan: break

            elif args.assembly_order == 'B_to_A_rescan':
                B_eval = choose_pair_B_alone(eval0, args)
                B_orig = original_or_none(rows0, B_eval)
                if B_eval is None or B_orig is None:
                    break
                K1, logB, appliedB = apply_pair(model, K, B_orig, args, event_t, cascade_index, 'B_C_kappa_context_first', option_tag(args))
                logB['proxy_B_score'] = p68.proxy_B_score(B_eval)
                logB['C_eigen_J_lock_mean_resid_eval'] = B_eval.get('best_C_eigen_J_lock_mean_resid','')
                logB['kappa_signed_flip_abs_eval'] = B_eval.get('kappa_signed_flip_abs','')
                pairing_log.append(logB)
                if not appliedB:
                    assembly_log.append({'variant':variant, 'option':option_tag(args), 'event_t':event_t, 'scan_id':scan_id, 'assembly_applied':False, 'failure':'B_apply_failed'})
                    break
                add_used_faces_from(B_orig, used_faces)
                global_pair_count += 1
                K = K1
                if global_pair_count >= args.max_auto_pairings:
                    break
                rows1 = nl.scan_rows(model, K, args)
                used_for_A = set() if args.allow_B_reuse_A_faces else used_faces
                eval1 = eval_rows_for_scan(model, K, rows1, args, used_for_A, cascade_index+1)
                for er in eval1:
                    er.update({'variant':variant, 'option':option_tag(args), 'event_t':event_t, 'scan_id':scan_id, 'cascade_index':cascade_index+1, 'stage':'post_B'})
                    candidate_rows.append(er)
                A_eval = choose_pair_A_with_context(eval1, B_eval, args)
                A_orig = original_or_none(rows1, A_eval)
                ctx = A_eval.get('_assembly_context',{}) if A_eval else {}
                if A_eval is None or A_orig is None:
                    assembly_log.append({'variant':variant, 'option':option_tag(args), 'event_t':event_t, 'scan_id':scan_id, 'assembly_applied':False, 'B_applied':True, 'A_attempted':False, 'failure':'no_A_candidate'})
                    cascade_index += 1; scan_id += 1
                    if not args.cascade_rescan: break
                    continue
                K2, logA, appliedA = apply_pair(model, K, A_orig, args, event_t, cascade_index+1, 'A_beta_QP_proxy_second', option_tag(args))
                logA['proxy_A_score'] = p68.proxy_A_score(A_eval)
                logA['C_eigen_J_lock_mean_resid_eval'] = A_eval.get('best_C_eigen_J_lock_mean_resid','')
                logA['kappa_signed_flip_abs_eval'] = A_eval.get('kappa_signed_flip_abs','')
                logA['assembly_context'] = ctx.get('context','')
                pairing_log.append(logA)
                if appliedA:
                    add_used_faces_from(A_orig, used_faces)
                    global_pair_count += 1
                    K = K2
                after = core.full_metrics(model, K, args.source)
                assembly_log.append(make_assembly_record(variant,args,event_t,scan_id,cascade_index,A_eval,B_eval,ctx,appliedA,after, order='B_to_A_rescan'))
                cascade_index += 2; scan_id += 1
                if not args.cascade_rescan: break

            elif args.assembly_order == 'stale_same_scan':
                A_eval, B_eval = choose_stale_AB(eval0, args)
                A_orig = original_or_none(rows0, A_eval)
                B_orig = original_or_none(rows0, B_eval)
                ctx = p68.support_context(A_eval, B_eval) if A_eval and B_eval else {}
                if A_eval is None or B_eval is None or A_orig is None or B_orig is None:
                    break
                K1, logA, appliedA = apply_pair(model, K, A_orig, args, event_t, cascade_index, 'A_beta_QP_proxy_stale', option_tag(args))
                logA['proxy_A_score'] = p68.proxy_A_score(A_eval)
                logA['C_eigen_J_lock_mean_resid_eval'] = A_eval.get('best_C_eigen_J_lock_mean_resid','')
                logA['kappa_signed_flip_abs_eval'] = A_eval.get('kappa_signed_flip_abs','')
                pairing_log.append(logA)
                if not appliedA:
                    assembly_log.append({'variant':variant, 'option':option_tag(args), 'event_t':event_t, 'scan_id':scan_id, 'assembly_applied':False, 'failure':'stale_A_apply_failed'})
                    break
                add_used_faces_from(A_orig, used_faces)
                global_pair_count += 1
                K = K1
                if global_pair_count >= args.max_auto_pairings:
                    break
                # Apply B from stale pre-A scan.  It may fail; that is the stale-legality diagnostic.
                K2, logB, appliedB = apply_pair(model, K, B_orig, args, event_t, cascade_index+1, 'B_C_kappa_context_stale', option_tag(args))
                logB['proxy_B_score'] = p68.proxy_B_score(B_eval)
                logB['C_eigen_J_lock_mean_resid_eval'] = B_eval.get('best_C_eigen_J_lock_mean_resid','')
                logB['kappa_signed_flip_abs_eval'] = B_eval.get('kappa_signed_flip_abs','')
                logB['assembly_context'] = ctx.get('context','')
                logB['stale_same_scan_B'] = True
                pairing_log.append(logB)
                if appliedB:
                    add_used_faces_from(B_orig, used_faces)
                    global_pair_count += 1
                    K = K2
                after = core.full_metrics(model, K, args.source)
                rec = make_assembly_record(variant,args,event_t,scan_id,cascade_index,A_eval,B_eval,ctx,appliedB,after, order='stale_same_scan')
                rec['stale_B_apply_reason'] = logB.get('apply_reason','')
                assembly_log.append(rec)
                cascade_index += 2; scan_id += 1
                if not args.cascade_rescan: break

            else:
                raise ValueError(f'Unknown assembly_order {args.assembly_order}')

    return K, birth_log, pairing_log, assembly_log, candidate_rows


def make_assembly_record(variant, args, event_t, scan_id, cascade_index, A_eval, B_eval, ctx, applied, after, order=None):
    order = order or args.assembly_order
    return {
        'variant': variant,
        'option': option_tag(args),
        'assembly_order': order,
        'event_t': event_t,
        'scan_id': scan_id,
        'cascade_index': cascade_index,
        'assembly_applied': bool(applied),
        'A_face_a': A_eval.get('face_a','') if A_eval else '',
        'A_face_b': A_eval.get('face_b','') if A_eval else '',
        'B_face_a': B_eval.get('face_a','') if B_eval else '',
        'B_face_b': B_eval.get('face_b','') if B_eval else '',
        'context': ctx.get('context','') if ctx else '',
        'face_overlap': ctx.get('face_overlap','') if ctx else '',
        'edge_overlap': ctx.get('edge_overlap','') if ctx else '',
        'vertex_overlap': ctx.get('vertex_overlap','') if ctx else '',
        'A_proxy_score': p68.proxy_A_score(A_eval) if A_eval else '',
        'A_QP_balance': qp_balance(A_eval) if A_eval else '',
        'A_C_lock': A_eval.get('best_C_eigen_J_lock_mean_resid','') if A_eval else '',
        'A_kappa_flip_abs': A_eval.get('kappa_signed_flip_abs','') if A_eval else '',
        'A_delta_beta2_audit_only': A_eval.get('delta_beta2_audit_only','') if A_eval else '',
        'B_proxy_score': p68.proxy_B_score(B_eval) if B_eval else '',
        'B_QP_balance': qp_balance(B_eval) if B_eval else '',
        'B_C_lock': B_eval.get('best_C_eigen_J_lock_mean_resid','') if B_eval else '',
        'B_kappa_flip_abs': B_eval.get('kappa_signed_flip_abs','') if B_eval else '',
        'B_delta_beta2_audit_only': B_eval.get('delta_beta2_audit_only','') if B_eval else '',
        'after_beta0': after['beta0'], 'after_beta1': after['beta1'], 'after_beta2': after['beta2'], 'after_beta3': after['beta3'],
        'decision_used_delta_beta': False,
    }


def summarize_selection(pairing_log: List[dict]) -> dict:
    return p68.summarize_selection(pairing_log)


def run_variant_option(variant: str, args: argparse.Namespace, out: Path) -> dict:
    model = nl.build_model(variant, args)
    model.grow(args.max_level)
    baseline_K = core.build_dynamic_outward_ngf_complex(model)
    baseline_metrics = core.full_metrics(model, baseline_K, args.source)
    tag = option_tag(args)
    vout = out / tag / variant
    vout.mkdir(parents=True, exist_ok=True)
    K, birth_log, pairing_log, assembly_log, candidate_rows = build_ablation_complex(model, args, variant, vout)
    auto_metrics = core.full_metrics(model, K, args.source)
    dm, pair_rows, top_rows, three_rows = p56.directed_metrics(model, K, pairing_log, args)
    sm, signed_rows, signed_face_rows = p59.signed_quadrature_rows(model, K, pairing_log, args)
    am, align_pair_rows, align_candidate_rows, align_candidate_summary = p61.alignment_search_metrics(model, K, pairing_log, args)
    sel = summarize_selection(pairing_log)
    write_csv(vout / 'birth_geometry_log.csv', birth_log)
    write_csv(vout / 'assembly_pairing_log.csv', pairing_log)
    write_csv(vout / 'assembly_ablation_log.csv', assembly_log)
    write_csv(vout / 'candidate_eval_rows.csv', candidate_rows)
    write_csv(vout / 'directed_pair_rows.csv', pair_rows)
    write_csv(vout / 'signed_quadrature_rows.csv', signed_rows)
    write_csv(vout / 'alignment_pair_rows.csv', align_pair_rows)
    write_csv(vout / 'alignment_candidate_rows.csv', align_candidate_rows)
    write_csv(vout / 'alignment_candidate_summary.csv', align_candidate_summary)
    summary = {
        'variant': variant,
        'option': tag,
        'assembly_order': args.assembly_order,
        'require_connected_assembly': args.require_connected_assembly,
        'require_strong_assembly_context': args.require_strong_assembly_context,
        'allow_B_reuse_A_faces': args.allow_B_reuse_A_faces,
        'baseline_metrics': baseline_metrics,
        'auto_metrics': auto_metrics,
        'directed_metrics': dm,
        'signed_quadrature_metrics': sm,
        'alignment_metrics': am,
        'selection_metrics': sel,
        'automatic_pairings_applied': sum(1 for x in pairing_log if fbool(x.get('applied'))),
        'assemblies_applied': sum(1 for x in assembly_log if fbool(x.get('assembly_applied'))),
        'assemblies_attempted': len(assembly_log),
        'decision_used_delta_beta_any': sel['decision_used_delta_beta_any'],
    }
    (vout / 'variant_option_summary.json').write_text(json.dumps(summary, indent=2), encoding='utf-8')
    return summary


def clone_args(args: argparse.Namespace, **updates) -> argparse.Namespace:
    d = vars(args).copy()
    d.update(updates)
    return argparse.Namespace(**d)


def run_experiment(args: argparse.Namespace, out: Path) -> List[dict]:
    options = []
    for order in args.orders:
        for ctx in args.context_modes:
            for reuse in args.reuse_modes:
                updates = {'assembly_order': order, 'allow_B_reuse_A_faces': (reuse == 'reuseB')}
                if ctx == 'anyctx':
                    updates.update({'require_connected_assembly': False, 'require_strong_assembly_context': False})
                elif ctx == 'connected':
                    updates.update({'require_connected_assembly': True, 'require_strong_assembly_context': False})
                elif ctx == 'strong':
                    updates.update({'require_connected_assembly': True, 'require_strong_assembly_context': True})
                else:
                    raise ValueError(ctx)
                options.append(clone_args(args, **updates))
    summaries: List[dict] = []
    for opt in options:
        for variant in args.variants:
            summaries.append(run_variant_option(variant, opt, out))
    return summaries


def write_comparative(out: Path, rows: List[dict]) -> None:
    flat = []
    for r in rows:
        a = r['auto_metrics']; dm = r['directed_metrics']; sm = r['signed_quadrature_metrics']; am = r['alignment_metrics']; sel = r['selection_metrics']
        flat.append({
            'option': r['option'],
            'variant': r['variant'],
            'order': r['assembly_order'],
            'context_connected': r['require_connected_assembly'],
            'context_strong': r['require_strong_assembly_context'],
            'allow_B_reuse_A_faces': r['allow_B_reuse_A_faces'],
            'beta0': a['beta0'], 'beta1': a['beta1'], 'beta2': a['beta2'], 'beta3': a['beta3'],
            'pairings': r['automatic_pairings_applied'],
            'assemblies_applied': r['assemblies_applied'],
            'assemblies_attempted': r['assemblies_attempted'],
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
        })
    write_csv(out / 'comparative_dual_assembly_order_context_ablation_summary.csv', flat)


def summarize_best(rows: List[dict]) -> dict:
    nonstrict = [r for r in rows if r['variant'] != 'strict_symmetrized_control']
    def key_signed(r): return abs(r['signed_quadrature_metrics']['signed_birth_over_abs_sum_ratio'])
    def key_j(r): return -r['alignment_metrics']['best_per_pair_mean_J_lock_resid']
    def key_beta(r): return r['auto_metrics']['beta2']
    def key_harm(r): return r['directed_metrics']['pair_transport_harmonic_ratio']
    best = {}
    if nonstrict:
        best['best_signed'] = slim(max(nonstrict, key=key_signed))
        best['best_J_lock'] = slim(max(nonstrict, key=key_j))
        best['best_beta2'] = slim(max(nonstrict, key=key_beta))
        best['best_pair_harm'] = slim(max(nonstrict, key=key_harm))
    return best


def slim(r: dict) -> dict:
    return {
        'option': r['option'],
        'variant': r['variant'],
        'beta2': r['auto_metrics']['beta2'],
        'pairings': r['automatic_pairings_applied'],
        'assemblies': r['assemblies_applied'],
        'pair_harm': r['directed_metrics']['pair_transport_harmonic_ratio'],
        'Q_harm': r['alignment_metrics']['Q_even_harmonic_ratio'],
        'P_harm': r['alignment_metrics']['P_odd_harmonic_ratio'],
        'J_lock': r['alignment_metrics']['best_per_pair_mean_J_lock_resid'],
        'signed_birth': r['signed_quadrature_metrics']['signed_birth_over_abs_sum_ratio'],
        'used_delta_beta': r['decision_used_delta_beta_any'],
    }


def make_docs(summary: dict) -> Tuple[str,str,str,str]:
    rows = summary['variant_rows']
    lines = [
        '| option | variant | beta | pairs | asm | pair harm | Q harm | P harm | J-lock | signed birth | used dBeta? |',
        '|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|',
    ]
    for r in rows:
        a = r['auto_metrics']; dm = r['directed_metrics']; am = r['alignment_metrics']; sm = r['signed_quadrature_metrics']
        lines.append(f"| {r['option']} | {r['variant']} | ({a['beta0']},{a['beta1']},{a['beta2']},{a['beta3']}) | {r['automatic_pairings_applied']} | {r['assemblies_applied']} | {dm['pair_transport_harmonic_ratio']:.6g} | {am['Q_even_harmonic_ratio']:.6g} | {am['P_odd_harmonic_ratio']:.6g} | {am['best_per_pair_mean_J_lock_resid']:.6g} | {sm['signed_birth_over_abs_sum_ratio']:.6g} | {r['decision_used_delta_beta_any']} |")
    table = '\n'.join(lines)
    best = json.dumps(summary.get('best_summary',{}), indent=2)
    smd = f"""# SUMMARY — dual assembly order/context ablation gate

Model label:
CNNA growing primal simplicial complex with deterministic sequential provenance growth, directed antisymmetric birth-transport operators, and dynamic/stale two-pair assembly ablations.

This package tests whether the two-pair assembly result depends on:

```text
A->B after rescan
B->A after rescan
stale same-scan A/B
connected versus strong context
allowing B to reuse A faces
```

No decision uses delta_beta, H2, complex scalars, Hodge, positivity, a physical adjoint, or final sym(M).

{table}

## Best-option summary

```json
{best}
```
"""
    rmd = f"""# RESULTS — dual assembly order/context ablation gate

## Comparative table

{table}

## Interpretation protocol

A robust constructive result would show one option with:

```text
strict_sym killed,
beta2 open,
Q/P harmonic positive,
signed_birth strong,
J-lock low,
decision_used_delta_beta_any false.
```

The ablation asks whether the signed/J-lock tension seen in the dynamic two-pair assembly is caused by order, stale legality, context restriction, or face reuse.
"""
    audit = """# SOURCE AUDIT

Carried forward:

- Single-pair property tradeoff gate: beta2/QP, C-lock, and kappa-flip split across candidates.
- Dual-pair assembly audit: same-scan two-pair motifs exist.
- Dual-pair dynamic growth gate: dynamic A->B assemblies open larger beta2 and Q/P channels but signed orientation and J-lock do not stabilize together.

This package does not introduce a new final rule. It ablates the dynamic assembly mechanics:

- A_to_B_rescan
- B_to_A_rescan
- stale_same_scan
- connected/strong/any context
- allow_B_reuse_A_faces on/off

Anti-smuggling constraints:

- delta_beta/H2/harmonic data remain audit-only.
- no i, no global J, no Hodge star, no positivity, no C*-norm, no physical adjoint, no final sym(M).
- stale same-scan is explicitly marked as a legality diagnostic, not as a claimed growth law.
"""
    readme = """# Dual assembly order/context ablation gate

Run:

```bash
python3 test_dual_assembly_order_context_ablation_gate.py
```

Default run uses a reduced ablation set:
A_to_B_rescan, B_to_A_rescan, stale_same_scan
with connected/strong contexts and B face reuse on/off.

Outputs include per-option/per-variant logs, comparative CSV/JSON, RESULTS.md, SUMMARY.md, SOURCE_AUDIT.md, and a ZIP package.
"""
    return smd, rmd, audit, readme


def package(out: Path, zip_path: Path) -> None:
    files = [
        Path(__file__).name,
        'test_dual_pairing_assembly_growth_rule_gate.py',
        'test_dual_pairing_two_edge_assembly_gate.py',
        'test_pair_property_tradeoff_obstruction_gate.py',
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
    ap.add_argument('--assembly-order', choices=['A_to_B_rescan','B_to_A_rescan','stale_same_scan'], default='A_to_B_rescan')
    ap.add_argument('--orders', nargs='*', default=['A_to_B_rescan','B_to_A_rescan','stale_same_scan'])
    ap.add_argument('--context-modes', nargs='*', default=['connected','strong'])
    ap.add_argument('--reuse-modes', nargs='*', default=['reuseB','noReuseB'])
    ap.add_argument('--variants', nargs='*', default=['real_growth','strict_symmetrized_control','no_backreaction'])
    ap.add_argument('--phase-sign', type=int, default=1)
    ap.add_argument('--out', default='dual_assembly_order_context_ablation_out_L2')
    ap.add_argument('--zip', default='cnna_dual_assembly_order_context_ablation_gate_pkg_L2.zip')
    args = ap.parse_args()

    out = Path(args.out)
    if out.exists():
        shutil.rmtree(out)
    out.mkdir(parents=True)
    rows = run_experiment(args, out)
    summary = {'args': vars(args), 'variant_rows': rows}
    summary['best_summary'] = summarize_best(rows)
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
        'best_summary': summary['best_summary'],
        'rows': [
            {
                'option': r['option'],
                'variant': r['variant'],
                'beta': [r['auto_metrics'][f'beta{i}'] for i in range(4)],
                'pairings': r['automatic_pairings_applied'],
                'assemblies': r['assemblies_applied'],
                'pair_harm': r['directed_metrics']['pair_transport_harmonic_ratio'],
                'Q_harm': r['alignment_metrics']['Q_even_harmonic_ratio'],
                'P_harm': r['alignment_metrics']['P_odd_harmonic_ratio'],
                'J_lock': r['alignment_metrics']['best_per_pair_mean_J_lock_resid'],
                'signed_birth': r['signed_quadrature_metrics']['signed_birth_over_abs_sum_ratio'],
                'used_delta_beta': r['decision_used_delta_beta_any'],
            } for r in rows
        ]
    }, indent=2))


if __name__ == '__main__':
    main()
