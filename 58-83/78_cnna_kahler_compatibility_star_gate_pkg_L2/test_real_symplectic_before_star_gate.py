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
import test_nonlinear_asymmetry_cascade_growth as nl
import test_pairing_transport_antisym_birth_coherence_gate as p56
import test_pair_J_alignment_search_gate as p61
import test_dual_assembly_order_context_ablation_gate as p69
import test_signed_Jlock_role_coupling_gate as p70
import test_assembly_motif_basis_diagonalization_gate as p71
import test_edge_interface_motif_operator_gate as p72
import test_shared_edge_link_cycle_operator_gate as p73

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


def fbool(x) -> bool:
    return p69.fbool(x)


def norm(x: np.ndarray) -> float:
    return float(np.linalg.norm(np.asarray(x, dtype=float)))


def orth(cols: List[np.ndarray], tol: float = 1e-10) -> np.ndarray:
    return p71.orth_basis(cols, tol=tol)


def mat_rank_svals(M: np.ndarray, tol: float = 1e-9) -> Tuple[int, List[float]]:
    if M.size == 0:
        return 0, []
    s = np.linalg.svd(M, compute_uv=False)
    mx = float(s[0]) if len(s) else 0.0
    r = int(np.sum(s > tol * max(mx, 1.0)))
    return r, [float(x) for x in s]


def skew_resid_matrix(M: np.ndarray) -> float:
    return norm(M + M.T) / (norm(M) + EPS)


def projected_form(S: np.ndarray, cols: List[np.ndarray]) -> Tuple[np.ndarray, np.ndarray]:
    U = orth(cols)
    if U.shape[1] == 0:
        return U, np.zeros((0, 0), dtype=float)
    return U, U.T @ S @ U


def symplectic_metrics(prefix: str, S: np.ndarray, qcols: List[np.ndarray], pcols: List[np.ndarray], args: argparse.Namespace) -> dict:
    allcols = qcols + pcols
    U, M = projected_form(S, allcols)
    dim = int(U.shape[1])
    rank, svals = mat_rank_svals(M, tol=args.singular_tol)
    max_sv = max(svals) if svals else 0.0
    min_nonzero = min([x for x in svals if x > args.singular_tol * max(max_sv, 1.0)], default=0.0)
    min_all = min(svals) if svals else 0.0
    nondeg_ratio = float(min_all / (max_sv + EPS)) if dim > 0 else 0.0
    effective_nonzero_ratio = float(min_nonzero / (max_sv + EPS)) if max_sv > 0 else 0.0

    Q = orth(qcols)
    P = orth(pcols)
    q_rank, p_rank = int(Q.shape[1]), int(P.shape[1])
    QQ = Q.T @ S @ Q if q_rank else np.zeros((0,0))
    PP = P.T @ S @ P if p_rank else np.zeros((0,0))
    QP = Q.T @ S @ P if q_rank and p_rank else np.zeros((0,0))
    qp_rank, qp_svals = mat_rank_svals(QP, tol=args.singular_tol)
    qp_max = max(qp_svals) if qp_svals else 0.0
    qp_min = min(qp_svals) if qp_svals else 0.0
    qp_ratio = float(qp_min / (qp_max + EPS)) if qp_svals else 0.0
    total_form_norm = norm(M)
    q_iso = norm(QQ) / (total_form_norm + EPS)
    p_iso = norm(PP) / (total_form_norm + EPS)
    qp_norm = norm(QP)
    qp_dominance = qp_norm / (total_form_norm + EPS)
    full_even = bool(dim > 0 and dim % 2 == 0)
    nondegenerate = bool(full_even and rank == dim and nondeg_ratio >= args.symplectic_ratio_threshold and skew_resid_matrix(M) <= args.skew_threshold)
    lagrangian_pairing = bool(
        q_rank == p_rank and q_rank > 0 and qp_rank == q_rank
        and qp_ratio >= args.qp_ratio_threshold
        and q_iso <= args.isotropic_threshold
        and p_iso <= args.isotropic_threshold
        and skew_resid_matrix(M) <= args.skew_threshold
    )
    # A form can be symplectic without Q/P being a clean Lagrangian split; both are logged separately.
    gate_pass = bool(nondegenerate and lagrangian_pairing)
    return {
        f'{prefix}_basis_dim': dim,
        f'{prefix}_form_rank': rank,
        f'{prefix}_skew_resid': skew_resid_matrix(M),
        f'{prefix}_max_sv': float(max_sv),
        f'{prefix}_min_sv': float(min_all),
        f'{prefix}_nondeg_ratio': float(nondeg_ratio),
        f'{prefix}_effective_nonzero_ratio': float(effective_nonzero_ratio),
        f'{prefix}_Q_rank': q_rank,
        f'{prefix}_P_rank': p_rank,
        f'{prefix}_QP_rank': qp_rank,
        f'{prefix}_QP_ratio': float(qp_ratio),
        f'{prefix}_Q_isotropic_resid': float(q_iso),
        f'{prefix}_P_isotropic_resid': float(p_iso),
        f'{prefix}_QP_dominance': float(qp_dominance),
        f'{prefix}_nondegenerate_gate': nondegenerate,
        f'{prefix}_lagrangian_QP_gate': lagrangian_pairing,
        f'{prefix}_symplectic_gate_pass': gate_pass,
        f'{prefix}_svals_json': json.dumps(svals),
        f'{prefix}_QP_svals_json': json.dumps(qp_svals),
    }


def data_wedge_form(qcols: List[np.ndarray], pcols: List[np.ndarray]) -> np.ndarray:
    # Diagnostic upper bound: uses Q/P data themselves, so it is NOT counted as a derived primary candidate.
    n = len(qcols[0]) if qcols else (len(pcols[0]) if pcols else 0)
    S = np.zeros((n, n), dtype=float)
    for q, p in zip(qcols, pcols):
        S += np.outer(q, p) - np.outer(p, q)
    return S


def pair_level_rows(p: dict, label: str, args: argparse.Namespace) -> List[dict]:
    q = np.concatenate([p['Q_a'], p['Q_b']])
    pp = np.concatenate([p['P_a'], p['P_b']])
    S_pair = p['J']
    S_data = data_wedge_form([q], [pp])
    rows = []
    for name, S, primary in [('pair_exchange', S_pair, True), ('data_wedge_control', S_data, False)]:
        m = symplectic_metrics(name, S, [q], [pp], args)
        rows.append({
            'level': 'pair', 'pair_label': label, 'omega_candidate': name, 'primary_candidate': primary,
            'face_a': str(list(p['fa'])), 'face_b': str(list(p['fb'])), **m,
        })
    return rows


def union_vecs(faces: List[Face], pA: dict, pB: dict) -> Dict[str, np.ndarray]:
    return {
        'A_Q': p71.put_pair_vec_union(faces, pA, 'Q'),
        'A_P': p71.put_pair_vec_union(faces, pA, 'P'),
        'B_Q': p71.put_pair_vec_union(faces, pB, 'Q'),
        'B_P': p71.put_pair_vec_union(faces, pB, 'P'),
    }


def motif_symplectic_rows(model, K, pA: dict, pB: dict, args: argparse.Namespace) -> List[dict]:
    rows: List[dict] = []
    faces = p72.union_faces_from_pairs(pA, pB)
    vecs = union_vecs(faces, pA, pB)
    qcols = [vecs['A_Q'], vecs['B_Q']]
    pcols = [vecs['A_P'], vecs['B_P']]
    Jbase, Cbase = p71.union_JC(faces, [pA, pB])
    candidates: List[Tuple[str, np.ndarray, bool, dict]] = [('union_pair_exchange', Jbase, True, {})]
    Sdata = data_wedge_form(qcols, pcols)
    candidates.append(('data_wedge_control', Sdata, False, {'control_reason':'uses_QP_data'}))
    links = p72.interface_links_for_assembly(model, pA, pB)
    for mode in args.interface_modes:
        Jint, Cint, istat = p72.build_interface_operator(model, faces, links, mode)
        candidates.append((f'edge_interface_only_{mode}', Jint, True, istat))
        for sc in args.interface_scales:
            Jeff, Ceff, lam = p72.scaled_add_base_interface(Jbase, Cbase, Jint, Cint, sc)
            candidates.append((f'union_plus_edge_interface_{mode}_{sc}', Jeff, True, {**istat, 'lambda_used': lam}))
    # Shared-edge link-cycle candidates from package 73, if shared edges exist.
    for edge in p73.shared_edges_between_pairs(pA, pB):
        for carrier_mode in ['motif_only']:
            lfaces = p73.carrier_faces_for_mode(K, faces, edge, carrier_mode)
            # lift motif vectors to link carrier when extra link faces are requested; motif_only means same faces here.
            lvecs = union_vecs(lfaces, pA, pB)
            lqcols = [lvecs['A_Q'], lvecs['B_Q']]
            lpcols = [lvecs['A_P'], lvecs['B_P']]
            Jb, Cb = p71.union_JC(lfaces, [pA, pB])
            for order_mode in args.link_order_modes:
                flip = p73.kappa_flip_stats(model, K, edge, order_mode, args)
                for block in args.link_block_modes:
                    Jcyc, Ccyc, cstat = p73.build_link_cycle_operator(model, K, lfaces, edge, order_mode, block)
                    candidates2 = [(f'link_cycle_only_{order_mode}_{block}_{edge}', Jcyc, True, {**cstat, **flip, 'shared_edge': str(list(edge))})]
                    for scale in args.link_scales:
                        Jeff, Ceff, lam = p73.scaled_add(Jb, Cb, Jcyc, Ccyc, scale)
                        candidates2.append((f'union_plus_link_cycle_{order_mode}_{block}_{scale}_{edge}', Jeff, True, {**cstat, **flip, 'shared_edge': str(list(edge)), 'lambda_used': lam}))
                    for cname, S, prim, extra in candidates2:
                        m = symplectic_metrics('omega', S, lqcols, lpcols, args)
                        rows.append({'level':'motif', 'carrier_mode': carrier_mode, 'omega_candidate': cname, 'primary_candidate': prim, **extra, **m})
    for cname, S, prim, extra in candidates:
        m = symplectic_metrics('omega', S, qcols, pcols, args)
        rows.append({'level':'motif', 'carrier_mode':'motif_only', 'omega_candidate': cname, 'primary_candidate': prim, **extra, **m})
    return rows


def parse_assembly_pairs(model, assembly_row: dict, args: argparse.Namespace) -> Optional[Tuple[dict, dict]]:
    A1, A2 = p71.parse_faces_from_assembly(assembly_row, 'A')
    B1, B2 = p71.parse_faces_from_assembly(assembly_row, 'B')
    if A1 is None or A2 is None or B1 is None or B2 is None:
        return None
    pA = p71.pair_fields(model, A1, A2, args)
    pB = p71.pair_fields(model, B1, B2, args)
    if pA is None or pB is None:
        return None
    return pA, pB


def symplectic_rows(model, K, assembly_log: List[dict], args: argparse.Namespace) -> Tuple[List[dict], dict]:
    rows: List[dict] = []
    for i, a in enumerate(assembly_log):
        if not fbool(a.get('assembly_applied')):
            continue
        parsed = parse_assembly_pairs(model, a, args)
        if parsed is None:
            continue
        pA, pB = parsed
        base_context = {'assembly_index': i, 'context': a.get('context',''), 'A_face_a': str(list(pA['fa'])), 'A_face_b': str(list(pA['fb'])), 'B_face_a': str(list(pB['fa'])), 'B_face_b': str(list(pB['fb']))}
        for r in pair_level_rows(pA, 'A', args) + pair_level_rows(pB, 'B', args):
            rows.append({**base_context, **r})
        for r in motif_symplectic_rows(model, K, pA, pB, args):
            rows.append({**base_context, **r})
    if not rows:
        return rows, {'row_count':0, 'primary_symplectic_pass_count':0, 'primary_motif_symplectic_pass_count':0, 'best_primary_nondeg_ratio':0.0, 'local_negative_marker': True}
    primary = [r for r in rows if fbool(r.get('primary_candidate'))]
    primary_motif = [r for r in primary if r.get('level') == 'motif']
    control = [r for r in rows if not fbool(r.get('primary_candidate'))]
    def best(rows, key, reverse=False, default=0.0):
        xs = [float(r.get(key,0.0)) for r in rows]
        if not xs: return default
        return (max(xs) if reverse else min(xs))
    def count_gate(rows, gate):
        return sum(1 for r in rows if fbool(r.get(gate)))
    best_primary = max(primary, key=lambda r: (fbool(r.get('omega_symplectic_gate_pass')), float(r.get('omega_nondeg_ratio',0.0)), float(r.get('omega_QP_ratio',0.0))), default=None)
    best_primary_motif = max(primary_motif, key=lambda r: (fbool(r.get('omega_symplectic_gate_pass')), float(r.get('omega_nondeg_ratio',0.0)), float(r.get('omega_QP_ratio',0.0))), default=None)
    summary = {
        'row_count': len(rows),
        'primary_row_count': len(primary),
        'control_row_count': len(control),
        'primary_symplectic_pass_count': count_gate(primary, 'omega_symplectic_gate_pass'),
        'primary_nondegenerate_count': count_gate(primary, 'omega_nondegenerate_gate'),
        'primary_lagrangian_QP_count': count_gate(primary, 'omega_lagrangian_QP_gate'),
        'primary_motif_symplectic_pass_count': count_gate(primary_motif, 'omega_symplectic_gate_pass'),
        'primary_motif_nondegenerate_count': count_gate(primary_motif, 'omega_nondegenerate_gate'),
        'primary_motif_lagrangian_QP_count': count_gate(primary_motif, 'omega_lagrangian_QP_gate'),
        'control_symplectic_pass_count': count_gate(control, 'omega_symplectic_gate_pass'),
        'best_primary_nondeg_ratio': float(best_primary.get('omega_nondeg_ratio',0.0)) if best_primary else 0.0,
        'best_primary_QP_ratio': float(best_primary.get('omega_QP_ratio',0.0)) if best_primary else 0.0,
        'best_primary_Q_isotropic': float(best_primary.get('omega_Q_isotropic_resid',0.0)) if best_primary else 0.0,
        'best_primary_P_isotropic': float(best_primary.get('omega_P_isotropic_resid',0.0)) if best_primary else 0.0,
        'best_primary_motif_nondeg_ratio': float(best_primary_motif.get('omega_nondeg_ratio',0.0)) if best_primary_motif else 0.0,
        'best_primary_motif_QP_ratio': float(best_primary_motif.get('omega_QP_ratio',0.0)) if best_primary_motif else 0.0,
        'best_primary_candidate': best_primary.get('omega_candidate','') if best_primary else '',
        'best_primary_level': best_primary.get('level','') if best_primary else '',
        'best_primary_dim': int(best_primary.get('omega_basis_dim',0)) if best_primary else 0,
        'best_primary_rank': int(best_primary.get('omega_form_rank',0)) if best_primary else 0,
        'best_primary_gate_pass': fbool(best_primary.get('omega_symplectic_gate_pass')) if best_primary else False,
        'best_primary_motif_candidate': best_primary_motif.get('omega_candidate','') if best_primary_motif else '',
        'best_primary_motif_dim': int(best_primary_motif.get('omega_basis_dim',0)) if best_primary_motif else 0,
        'best_primary_motif_rank': int(best_primary_motif.get('omega_form_rank',0)) if best_primary_motif else 0,
        'best_primary_motif_gate_pass': fbool(best_primary_motif.get('omega_symplectic_gate_pass')) if best_primary_motif else False,
        'union_pair_exchange_motif_pass_count': sum(1 for r in primary_motif if r.get('omega_candidate') == 'union_pair_exchange' and fbool(r.get('omega_symplectic_gate_pass'))),
        'local_symplectic_positive_marker': bool(count_gate(primary_motif, 'omega_symplectic_gate_pass') > 0 or count_gate(primary, 'omega_symplectic_gate_pass') > 0),
        'local_negative_marker': bool(count_gate(primary, 'omega_symplectic_gate_pass') == 0),
    }
    return rows, summary


def option_tag(args: argparse.Namespace) -> str:
    return p69.option_tag(args)


def clone_args(args: argparse.Namespace, **updates) -> argparse.Namespace:
    d = vars(args).copy(); d.update(updates); return argparse.Namespace(**d)


def run_variant(variant: str, args: argparse.Namespace, out: Path) -> dict:
    model = nl.build_model(variant, args)
    model.grow(args.max_level)
    baseline_K = core.build_dynamic_outward_ngf_complex(model)
    baseline_metrics = core.full_metrics(model, baseline_K, args.source)
    tag = option_tag(args)
    vout = out / tag / variant
    vout.mkdir(parents=True, exist_ok=True)
    K, birth_log, pairing_log, assembly_log, candidate_rows = p69.build_ablation_complex(model, args, variant, vout)
    auto_metrics = core.full_metrics(model, K, args.source)
    dm, pair_rows, top_rows, three_rows = p56.directed_metrics(model, K, pairing_log, args)
    sm, signed_rows, signed_face_rows = p69.p59.signed_quadrature_rows(model, K, pairing_log, args)
    am, align_pair_rows, align_candidate_rows, align_candidate_summary = p61.alignment_search_metrics(model, K, pairing_log, args)
    sel = p69.summarize_selection(pairing_log)
    motif_rows, motif_summary = p71.assembly_motif_rows(model, K, pairing_log, assembly_log, args)
    symp_rows, symp_summary = symplectic_rows(model, K, assembly_log, args)
    write_csv(vout / 'birth_geometry_log.csv', birth_log)
    write_csv(vout / 'assembly_pairing_log.csv', pairing_log)
    write_csv(vout / 'assembly_ablation_log.csv', assembly_log)
    write_csv(vout / 'candidate_eval_rows.csv', candidate_rows)
    write_csv(vout / 'directed_pair_rows.csv', pair_rows)
    write_csv(vout / 'signed_quadrature_rows.csv', signed_rows)
    write_csv(vout / 'alignment_pair_rows.csv', align_pair_rows)
    write_csv(vout / 'alignment_candidate_rows.csv', align_candidate_rows)
    write_csv(vout / 'alignment_candidate_summary.csv', align_candidate_summary)
    write_csv(vout / 'assembly_motif_basis_rows.csv', motif_rows)
    write_csv(vout / 'real_symplectic_before_star_rows.csv', symp_rows)
    summary = {
        'variant': variant,
        'option': tag,
        'baseline_metrics': baseline_metrics,
        'auto_metrics': auto_metrics,
        'directed_metrics': dm,
        'signed_quadrature_metrics': sm,
        'alignment_metrics': am,
        'selection_metrics': sel,
        'motif_basis_metrics': motif_summary,
        'real_symplectic_metrics': symp_summary,
        'automatic_pairings_applied': sum(1 for x in pairing_log if fbool(x.get('applied'))),
        'assemblies_applied': sum(1 for x in assembly_log if fbool(x.get('assembly_applied'))),
        'assemblies_attempted': len(assembly_log),
        'decision_used_delta_beta_any': sel['decision_used_delta_beta_any'],
    }
    (vout / 'variant_real_symplectic_before_star_summary.json').write_text(json.dumps(summary, indent=2), encoding='utf-8')
    return summary


def run_experiment(args: argparse.Namespace, out: Path) -> List[dict]:
    p70.patch_modules()
    rows = []
    options = []
    for order in args.orders:
        for ctx in args.context_modes:
            for reuse in args.reuse_modes:
                updates = {'assembly_order': order, 'allow_B_reuse_A_faces': (reuse == 'reuseB')}
                if ctx == 'connected':
                    updates.update({'require_connected_assembly': True, 'require_strong_assembly_context': False})
                elif ctx == 'strong':
                    updates.update({'require_connected_assembly': True, 'require_strong_assembly_context': True})
                elif ctx == 'anyctx':
                    updates.update({'require_connected_assembly': False, 'require_strong_assembly_context': False})
                else:
                    raise ValueError(ctx)
                options.append(clone_args(args, **updates))
    for opt in options:
        for variant in opt.variants:
            rows.append(run_variant(variant, opt, out))
    return rows


def slim(r: dict) -> dict:
    a = r['auto_metrics']; dm = r['directed_metrics']; sm = r['signed_quadrature_metrics']; am = r['alignment_metrics']; mm = r['motif_basis_metrics']; sy = r['real_symplectic_metrics']
    return {
        'option': r['option'], 'variant': r['variant'],
        'beta': [a['beta0'], a['beta1'], a['beta2'], a['beta3']],
        'pairings': r['automatic_pairings_applied'], 'assemblies': r['assemblies_applied'],
        'pair_harm': dm['pair_transport_harmonic_ratio'],
        'Q_harm': am['Q_even_harmonic_ratio'], 'P_harm': am['P_odd_harmonic_ratio'],
        'pair_local_J_lock': am['best_per_pair_mean_J_lock_resid'],
        'motif_union_J_lock': mm.get('union_sum_best_mean_resid',0.0),
        'signed_birth': sm['signed_birth_over_abs_sum_ratio'],
        'symp_primary_pass': sy.get('primary_symplectic_pass_count',0),
        'symp_primary_motif_pass': sy.get('primary_motif_symplectic_pass_count',0),
        'symp_primary_nondeg_count': sy.get('primary_nondegenerate_count',0),
        'symp_primary_lagrangian_count': sy.get('primary_lagrangian_QP_count',0),
        'best_symp_candidate': sy.get('best_primary_candidate',''),
        'best_symp_level': sy.get('best_primary_level',''),
        'best_symp_dim': sy.get('best_primary_dim',0),
        'best_symp_rank': sy.get('best_primary_rank',0),
        'best_symp_nondeg_ratio': sy.get('best_primary_nondeg_ratio',0.0),
        'best_symp_QP_ratio': sy.get('best_primary_QP_ratio',0.0),
        'best_motif_symp_candidate': sy.get('best_primary_motif_candidate',''),
        'best_motif_symp_dim': sy.get('best_primary_motif_dim',0),
        'best_motif_symp_rank': sy.get('best_primary_motif_rank',0),
        'best_motif_symp_nondeg_ratio': sy.get('best_primary_motif_nondeg_ratio',0.0),
        'best_motif_symp_QP_ratio': sy.get('best_primary_motif_QP_ratio',0.0),
        'used_delta_beta': r['decision_used_delta_beta_any'],
    }


def write_comparative(out: Path, rows: List[dict]) -> None:
    flat = []
    for r in rows:
        s = slim(r); sy = r['real_symplectic_metrics']
        flat.append({**s,
            'beta0': s['beta'][0], 'beta1': s['beta'][1], 'beta2': s['beta'][2], 'beta3': s['beta'][3],
            **{k: sy.get(k, '') for k in sorted(sy.keys()) if not isinstance(sy.get(k), (dict, list))}
        })
    write_csv(out / 'comparative_real_symplectic_before_star_summary.csv', flat)


def make_docs(summary: dict) -> Tuple[str, str, str, str]:
    rows = summary['variant_rows']
    lines = ['| option | variant | beta | pairs | asm | pair harm | Q harm | P harm | pair J-lock | union J-lock | signed | primary symp pass | motif symp pass | best primary Ω | level | dim/rank | Ω nondeg ratio | Ω QP ratio | used dβ? |',
             '|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|---|---:|---:|---:|---:|']
    for r in rows:
        s = slim(r)
        lines.append(f"| {s['option']} | {s['variant']} | ({s['beta'][0]},{s['beta'][1]},{s['beta'][2]},{s['beta'][3]}) | {s['pairings']} | {s['assemblies']} | {s['pair_harm']:.6g} | {s['Q_harm']:.6g} | {s['P_harm']:.6g} | {s['pair_local_J_lock']:.6g} | {s['motif_union_J_lock']:.6g} | {s['signed_birth']:.6g} | {s['symp_primary_pass']} | {s['symp_primary_motif_pass']} | {s['best_symp_candidate']} | {s['best_symp_level']} | {s['best_symp_dim']}/{s['best_symp_rank']} | {s['best_symp_nondeg_ratio']:.6g} | {s['best_symp_QP_ratio']:.6g} | {s['used_delta_beta']} |")
    table = '\n'.join(lines)
    nonstrict = [r for r in rows if r['variant'] != 'strict_symmetrized_control']
    any_primary_pass = any(r['real_symplectic_metrics'].get('primary_symplectic_pass_count',0) > 0 for r in nonstrict)
    any_motif_pass = any(r['real_symplectic_metrics'].get('primary_motif_symplectic_pass_count',0) > 0 for r in nonstrict)
    smd = f"""# SUMMARY — real symplectic before star gate

Model label:
CNNA growing primal simplicial complex with deterministic sequential provenance growth, directed antisymmetric birth-transport operators, signed-Jlock two-pair assembly motifs, and a real symplectic-first diagnostic.

Corrected ladder order:

```text
real Q/P carrier
→ real symplectic form Ω
→ possible derived #/*-structure
→ only then J / complex orientation
```

This package therefore does **not** use `J²≈-I` as the primary gate.  It asks first whether any already-derived skew operator candidate restricts to a nondegenerate real two-form on the actual Q/P or A/B motif space.

Primary Ω candidates:

```text
pair_exchange / union_pair_exchange
edge_interface_only
union_plus_edge_interface
link_cycle_only
union_plus_link_cycle
```

The `data_wedge_control` is logged only as a tautological upper-bound control because it uses Q/P data directly.

Positive primary symplectic gate requires:

```text
projected Ω skew residual <= {summary['args']['skew_threshold']}
projected Ω full rank on an even-dimensional Q/P span
nondegeneracy ratio >= {summary['args']['symplectic_ratio_threshold']}
Q and P have equal positive rank
Ω(Q,P) is full rank with ratio >= {summary['args']['qp_ratio_threshold']}
Q and P are approximately Ω-isotropic <= {summary['args']['isotropic_threshold']}
strict_sym remains null
used_delta_beta remains false
```

{table}

## High-level decision

```json
{{
  "any_non_strict_primary_symplectic_pass": {str(any_primary_pass).lower()},
  "any_non_strict_primary_motif_symplectic_pass": {str(any_motif_pass).lower()}
}}
```
"""
    rmd = f"""# RESULTS — real symplectic before star gate

## Comparative table

{table}

## Interpretation rule

This package separates three things that earlier tests tended to entangle:

```text
1. Q/P carrier existence,
2. real symplectic nondegeneracy of Ω on Q/P,
3. J-lock / complex orientation.
```

A positive Ω result does not mean `i` or `J` has been derived.  It means the ladder can advance to a derived real #/* search.  Important audit detail: `union_pair_exchange` already passes on the motif in the non-strict variants; the link-cycle variants are not needed to make Ω nondegenerate, even when they tie or duplicate the best score.  A negative Ω result means the obstruction is even earlier than # or J: the Q/P carrier exists, but no tested derived real symplectic form has been found on it.

## Stop/continue rule

- If a primary Ω candidate passes in non-strict variants and strict_sym is null, the next test should be `test_real_star_from_symplectic_gate.py`.
- If only `data_wedge_control` passes, the Q/P plane is independent but the symplectic form is not derived; do not count it as success.
- If no primary Ω passes, document a symplectic-level obstruction before returning to J-search.
"""
    audit = """# SOURCE AUDIT

No i, global J, Hodge star, physical adjoint, positivity, C*-norm, final sym(M), or delta-beta/H2 decision is introduced.

The tested primary two-forms are the skew parts already present in earlier derived operators: pair exchange, edge-interface, and shared-edge link-cycle variants.  The data-wedge control is included explicitly as a non-primary tautological diagnostic only.

This is a Python diagnostic on the L2 model, not a Lean theorem.
"""
    readme = """# Real symplectic before star gate

Run:

```bash
python3 test_real_symplectic_before_star_gate.py
```

The test checks whether the already-derived Q/P carrier admits a primary real symplectic form before any #/* or J claim is made.
"""
    return smd, rmd, audit, readme


def package(out: Path, zip_path: Path) -> None:
    files = [
        Path(__file__).name,
        'test_shared_edge_link_cycle_operator_gate.py',
        'test_edge_interface_motif_operator_gate.py',
        'test_assembly_motif_basis_diagonalization_gate.py',
        'test_signed_Jlock_role_coupling_gate.py',
        'test_dual_assembly_order_context_ablation_gate.py',
        'test_dual_pairing_assembly_growth_rule_gate.py',
        'test_dual_pairing_two_edge_assembly_gate.py',
        'test_pair_property_tradeoff_obstruction_gate.py',
        'test_C_eigen_quadrature_refinement_gate.py',
        'test_C_eigen_guided_pairing_rule_gate.py',
        'test_pair_J_alignment_search_gate.py',
        'test_pairing_quadrature_adjoint_pairing_gate.py',
        'test_signed_quadrature_area_kappa_gate.py',
        'test_pairing_quadrature_split_symplectic_defect_gate.py',
        'test_pairing_transport_antisym_birth_coherence_gate.py',
        'test_pairing_transport_harmonic_kappa_gate.py',
        'test_nonlinear_asymmetry_cascade_growth.py',
        'test_harmonic_k_orientation_kappa_gate.py',
        'test_interfan_transport_from_asymmetry_invariants.py',
        'test_growth_with_asymmetry_gated_complement_pairing.py',
        'cnna_non_shelling_core.py',
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
    ap.add_argument('--require-strong-assembly-context', action='store_true', default=True)
    ap.add_argument('--assembly-order', choices=['A_to_B_rescan','B_to_A_rescan','stale_same_scan'], default='A_to_B_rescan')
    ap.add_argument('--orders', nargs='*', default=['A_to_B_rescan'])
    ap.add_argument('--context-modes', nargs='*', default=['strong'])
    ap.add_argument('--reuse-modes', nargs='*', default=['reuseB'])
    ap.add_argument('--variants', nargs='*', default=['real_growth','strict_symmetrized_control','no_backreaction'])
    ap.add_argument('--phase-sign', type=int, default=1)
    ap.add_argument('--signed-comm-threshold', type=float, default=0.10)
    ap.add_argument('--interface-modes', nargs='*', default=['incidence_identity','edge_projector','edge_complement_projector'])
    ap.add_argument('--interface-scales', nargs='*', default=['unit','half_base_norm','base_norm'])
    ap.add_argument('--link-order-modes', nargs='*', default=['birth_order','address_order','geometric_angle'])
    ap.add_argument('--link-block-modes', nargs='*', default=['identity','edge_projector','edge_complement'])
    ap.add_argument('--link-scales', nargs='*', default=['unit','quarter_base_norm','half_base_norm','base_norm'])
    ap.add_argument('--link-circulation-threshold', type=float, default=1e-6)
    ap.add_argument('--kappa-flip-threshold', type=float, default=0.20)
    ap.add_argument('--singular-tol', type=float, default=1e-9)
    ap.add_argument('--skew-threshold', type=float, default=1e-8)
    ap.add_argument('--symplectic-ratio-threshold', type=float, default=1e-3)
    ap.add_argument('--qp-ratio-threshold', type=float, default=1e-3)
    ap.add_argument('--isotropic-threshold', type=float, default=0.35)
    ap.add_argument('--out', default='real_symplectic_before_star_out_L2')
    ap.add_argument('--zip', default='cnna_real_symplectic_before_star_gate_pkg_L2.zip')
    args = ap.parse_args()
    out = Path(args.out)
    if out.exists():
        shutil.rmtree(out)
    out.mkdir(parents=True)
    rows = run_experiment(args, out)
    write_comparative(out, rows)
    summary = {'args': vars(args), 'variant_rows': rows}
    (out / 'comparative_summary.json').write_text(json.dumps(summary, indent=2), encoding='utf-8')
    smd, rmd, audit, readme = make_docs(summary)
    (out / 'SUMMARY.md').write_text(smd, encoding='utf-8')
    (out / 'RESULTS.md').write_text(rmd, encoding='utf-8')
    (out / 'SOURCE_AUDIT.md').write_text(audit, encoding='utf-8')
    (out / 'README.md').write_text(readme, encoding='utf-8')
    package(out, Path(args.zip))
    print(json.dumps({'zip': args.zip, 'out': args.out, 'rows': [slim(r) for r in rows]}, indent=2))


if __name__ == '__main__':
    main()
