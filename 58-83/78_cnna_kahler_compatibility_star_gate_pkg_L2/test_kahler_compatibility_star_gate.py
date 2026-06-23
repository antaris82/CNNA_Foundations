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
import test_real_symplectic_before_star_gate as p74

EPS = 1e-12
Face = Tuple[int, int, int]


def fbool(x) -> bool:
    if isinstance(x, str):
        return x.lower() in {'true', '1', 'yes'}
    return bool(x)


def norm(x: np.ndarray) -> float:
    return float(np.linalg.norm(np.asarray(x, dtype=float)))


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


def orth(cols: List[np.ndarray], tol: float = 1e-10) -> np.ndarray:
    return p71.orth_basis(cols, tol=tol)


def projector(U: np.ndarray) -> np.ndarray:
    if U.size == 0 or U.shape[1] == 0:
        return np.zeros((U.shape[0], U.shape[0]), dtype=float)
    return U @ U.T


def mat_rank_svals(M: np.ndarray, tol: float = 1e-9) -> Tuple[int, List[float]]:
    if M.size == 0:
        return 0, []
    s = np.linalg.svd(M, compute_uv=False)
    mx = float(s[0]) if len(s) else 0.0
    r = int(np.sum(s > tol * max(mx, 1.0)))
    return r, [float(x) for x in s]


def nondeg_ratio(M: np.ndarray, tol: float = 1e-9) -> float:
    if M.size == 0:
        return 0.0
    s = np.linalg.svd(M, compute_uv=False)
    if len(s) == 0:
        return 0.0
    return float(s[-1] / (s[0] + EPS))


def skew_resid(M: np.ndarray) -> float:
    return norm(M + M.T) / (norm(M) + EPS)


def sym_resid(M: np.ndarray) -> float:
    return norm(M - M.T) / (norm(M) + EPS)


def subspace_image_residual(J: np.ndarray, source_cols: List[np.ndarray], target_cols: List[np.ndarray]) -> float:
    S = orth(source_cols)
    T = orth(target_cols)
    if S.shape[1] == 0:
        return 0.0
    JS = J @ S
    den = norm(JS)
    if den < EPS:
        return 0.0
    PT = projector(T)
    return norm((np.eye(J.shape[0]) - PT) @ JS) / (den + EPS)


def span_residual_to_ops(A: np.ndarray, ops: List[np.ndarray]) -> float:
    if not ops:
        return 1.0
    target = A.reshape(-1)
    B = np.column_stack([O.reshape(-1) for O in ops])
    if norm(target) < EPS:
        return 0.0
    try:
        coeff, *_ = np.linalg.lstsq(B, target, rcond=None)
        pred = B @ coeff
        return norm(target - pred) / (norm(target) + EPS)
    except Exception:
        return 1.0


def metric_adjoint(A: np.ndarray, G: np.ndarray) -> np.ndarray:
    return np.linalg.pinv(G, rcond=1e-10) @ A.T @ G


def projected_data(faces: List[Face], pA: dict, pB: dict) -> Tuple[Dict[str, np.ndarray], List[np.ndarray], List[np.ndarray], List[np.ndarray], np.ndarray]:
    vecs = p74.union_vecs(faces, pA, pB)
    qcols = [vecs['A_Q'], vecs['B_Q']]
    pcols = [vecs['A_P'], vecs['B_P']]
    allcols = qcols + pcols
    U = orth(allcols)
    return vecs, qcols, pcols, allcols, U


def metric_candidates(n: int, U: np.ndarray, faces: List[Face], pA: dict, pB: dict, Jbase: np.ndarray, Cbase: np.ndarray, qcols: List[np.ndarray], pcols: List[np.ndarray], args: argparse.Namespace) -> List[Tuple[str, np.ndarray, bool, str]]:
    # Primary candidates are structural metrics already present in the motif/cochain representation.
    # Controls deliberately use Q/P data or omega-derived polar algebra and are not counted as primary.
    I = np.eye(n)
    Csym = 0.5 * (Cbase + Cbase.T)
    Csq = Cbase.T @ Cbase
    # Data Gram control on the Q/P carrier, lifted with a tiny identity regularizer only for numerical inversion.
    D = np.zeros((n, n), dtype=float)
    for v in qcols + pcols:
        D += np.outer(v, v)
    if norm(D) > EPS:
        D = D / (norm(D) + EPS)
    return [
        ('cochain_identity_metric', I, True, 'ambient coefficient metric already used by residuals'),
        ('union_C_symmetric_metric', Csym, True, 'symmetric pair-conjugation operator C on motif union'),
        ('union_C_square_metric', Csq, True, 'C^T C structural metric diagnostic'),
        ('QP_data_gram_control', D + args.metric_regularizer * I, False, 'uses Q/P data directly'),
    ]


def omega_candidates(model, K, faces: List[Face], pA: dict, pB: dict, args: argparse.Namespace) -> List[Tuple[str, np.ndarray, bool, str, dict]]:
    Jbase, Cbase = p71.union_JC(faces, [pA, pB])
    out: List[Tuple[str, np.ndarray, bool, str, dict]] = []
    # These are the already-tested structural skew operators.  They are NOT independent of earlier J-like pair exchange;
    # the audit reports that explicitly.
    out.append(('union_pair_exchange_skew', 0.5 * (Jbase - Jbase.T), True, 'skew part of union pair-exchange J-like operator; not independent of prior J tests', {}))
    links = p72.interface_links_for_assembly(model, pA, pB)
    for mode in args.interface_modes:
        Jint, Cint, istat = p72.build_interface_operator(model, faces, links, mode)
        out.append((f'edge_interface_skew_{mode}', 0.5 * (Jint - Jint.T), True, 'skew part of Face->shared-edge->Face interface operator', istat))
        for sc in args.interface_scales:
            Jeff, Ceff, lam = p72.scaled_add_base_interface(Jbase, Cbase, Jint, Cint, sc)
            out.append((f'union_plus_edge_interface_skew_{mode}_{sc}', 0.5 * (Jeff - Jeff.T), True, 'skew part of union plus edge interface', {**istat, 'lambda_used': lam}))
    for edge in p73.shared_edges_between_pairs(pA, pB):
        for order_mode in args.link_order_modes:
            flip = p73.kappa_flip_stats(model, K, edge, order_mode, args)
            for block in args.link_block_modes:
                Jcyc, Ccyc, cstat = p73.build_link_cycle_operator(model, K, faces, edge, order_mode, block)
                out.append((f'link_cycle_skew_{order_mode}_{block}_{edge}', 0.5 * (Jcyc - Jcyc.T), True, 'skew part of shared-edge link-cycle operator', {**cstat, **flip, 'shared_edge': str(list(edge))}))
                for scale in args.link_scales:
                    Jeff, Ceff, lam = p73.scaled_add(Jbase, Cbase, Jcyc, Ccyc, scale)
                    out.append((f'union_plus_link_cycle_skew_{order_mode}_{block}_{scale}_{edge}', 0.5 * (Jeff - Jeff.T), True, 'skew part of union plus shared-edge link-cycle', {**cstat, **flip, 'shared_edge': str(list(edge)), 'lambda_used': lam}))
    # Explicit tautological data wedge control: not primary.
    vecs, qcols, pcols, allcols, U = projected_data(faces, pA, pB)
    Sdata = p74.data_wedge_form(qcols, pcols)
    out.append(('QP_data_wedge_control', Sdata, False, 'uses Q/P vectors directly; tautological upper bound', {}))
    return out


def compatibility_metrics(prefix: str, Omega: np.ndarray, G: np.ndarray, Cbase: np.ndarray, qcols: List[np.ndarray], pcols: List[np.ndarray], allcols: List[np.ndarray], args: argparse.Namespace) -> dict:
    U = orth(allcols)
    dim = int(U.shape[1])
    if dim == 0:
        return {f'{prefix}_basis_dim': 0, f'{prefix}_gate_pass': False}
    Op = U.T @ Omega @ U
    Gp = U.T @ G @ U
    Cp = U.T @ Cbase @ U
    qproj = [U.T @ q for q in qcols]
    pproj = [U.T @ p for p in pcols]
    Ir = np.eye(dim)
    omega_rank, omega_s = mat_rank_svals(Op, args.singular_tol)
    g_rank, g_s = mat_rank_svals(Gp, args.singular_tol)
    omega_nd = nondeg_ratio(Op, args.singular_tol)
    g_nd = nondeg_ratio(Gp, args.singular_tol)
    valid = bool(dim > 0 and dim % 2 == 0 and omega_rank == dim and g_rank == dim and omega_nd >= args.nondeg_threshold and g_nd >= args.nondeg_threshold and skew_resid(Op) <= args.skew_threshold and sym_resid(Gp) <= args.sym_threshold)
    if not valid:
        J = np.zeros_like(Op)
    else:
        J = np.linalg.pinv(Gp, rcond=1e-10) @ Op
    J2 = norm(J @ J + Ir) / (norm(Ir) + EPS) if dim else 0.0
    J_metric_orth = norm(J.T @ Gp @ J - Gp) / (norm(Gp) + EPS) if dim else 0.0
    GJ_minus_O = norm(Gp @ J - Op) / (norm(Op) + EPS) if norm(Op) > EPS else 0.0
    J_g_skew = norm(J.T @ Gp + Gp @ J) / (norm(Gp @ J) + EPS) if dim and norm(Gp @ J) > EPS else 0.0
    q_to_p = subspace_image_residual(J, qproj, pproj)
    p_to_q = subspace_image_residual(J, pproj, qproj)
    mean_lock = 0.5 * (q_to_p + p_to_q)
    max_lock = max(q_to_p, p_to_q)
    J_hash = metric_adjoint(J, Gp) if valid else np.zeros_like(J)
    C_hash = metric_adjoint(Cp, Gp) if valid else np.zeros_like(Cp)
    J_anti_hash = norm(J_hash + J) / (norm(J) + EPS) if norm(J) > EPS else 0.0
    C_self_hash = norm(C_hash - Cp) / (norm(Cp) + EPS) if norm(Cp) > EPS else 0.0
    ops = [Ir, Cp, J, Op]
    star_span = max(span_residual_to_ops(metric_adjoint(A, Gp), ops) for A in ops) if valid else 1.0
    # Metric adjoint identities are tautological for invertible G; logged but not a success condition.
    star_involution = 0.0
    star_antimult = 0.0
    if valid:
        A, B = Cp, J
        A2 = metric_adjoint(metric_adjoint(A, Gp), Gp)
        star_involution = norm(A2 - A) / (norm(A) + EPS) if norm(A) > EPS else 0.0
        lhs = metric_adjoint(A @ B, Gp)
        rhs = metric_adjoint(B, Gp) @ metric_adjoint(A, Gp)
        star_antimult = norm(lhs - rhs) / (norm(lhs) + EPS) if norm(lhs) > EPS else 0.0
    gate = bool(
        valid
        and J2 <= args.compat_J2_threshold
        and mean_lock <= args.compat_lock_mean_threshold
        and max_lock <= args.compat_lock_max_threshold
        and J_metric_orth <= args.metric_orth_threshold
        and J_anti_hash <= args.hash_anti_threshold
        and star_span <= args.star_span_threshold
    )
    return {
        f'{prefix}_basis_dim': dim,
        f'{prefix}_omega_rank': omega_rank,
        f'{prefix}_g_rank': g_rank,
        f'{prefix}_omega_nondeg_ratio': float(omega_nd),
        f'{prefix}_g_nondeg_ratio': float(g_nd),
        f'{prefix}_omega_skew_resid': skew_resid(Op),
        f'{prefix}_g_sym_resid': sym_resid(Gp),
        f'{prefix}_valid_omega_g': valid,
        f'{prefix}_J2_plus_I_resid': float(J2),
        f'{prefix}_J_metric_orth_resid': float(J_metric_orth),
        f'{prefix}_GJ_minus_Omega_resid': float(GJ_minus_O),
        f'{prefix}_J_g_skew_resid': float(J_g_skew),
        f'{prefix}_J_Q_to_P_resid': float(q_to_p),
        f'{prefix}_J_P_to_Q_resid': float(p_to_q),
        f'{prefix}_J_QP_mean_lock_resid': float(mean_lock),
        f'{prefix}_J_QP_max_lock_resid': float(max_lock),
        f'{prefix}_J_hash_anti_self_resid': float(J_anti_hash),
        f'{prefix}_C_hash_self_resid': float(C_self_hash),
        f'{prefix}_star_family_span_resid': float(star_span),
        f'{prefix}_star_involution_resid_diagnostic': float(star_involution),
        f'{prefix}_star_antimult_resid_diagnostic': float(star_antimult),
        f'{prefix}_compatibility_gate_pass': gate,
        f'{prefix}_omega_svals_json': json.dumps(omega_s),
        f'{prefix}_g_svals_json': json.dumps(g_s),
    }


def compatibility_rows_for_assembly(model, K, assembly_row: dict, args: argparse.Namespace) -> List[dict]:
    parsed = p74.parse_assembly_pairs(model, assembly_row, args)
    if parsed is None:
        return []
    pA, pB = parsed
    faces = p72.union_faces_from_pairs(pA, pB)
    Jbase, Cbase = p71.union_JC(faces, [pA, pB])
    vecs, qcols, pcols, allcols, U = projected_data(faces, pA, pB)
    n = len(faces) * 3
    rows: List[dict] = []
    omegas = omega_candidates(model, K, faces, pA, pB, args)
    metrics = metric_candidates(n, U, faces, pA, pB, Jbase, Cbase, qcols, pcols, args)
    base = {
        'context': assembly_row.get('context',''),
        'A_face_a': str(list(pA['fa'])), 'A_face_b': str(list(pA['fb'])),
        'B_face_a': str(list(pB['fa'])), 'B_face_b': str(list(pB['fb'])),
        'union_faces': str([list(f) for f in faces]),
    }
    for oname, O, oprim, oreason, oextra in omegas:
        for gname, G, gprim, greason in metrics:
            m = compatibility_metrics('compat', O, G, Cbase, qcols, pcols, allcols, args)
            rows.append({
                **base,
                'omega_candidate': oname,
                'omega_primary': oprim,
                'omega_source_note': oreason,
                'metric_candidate': gname,
                'metric_primary': gprim,
                'metric_source_note': greason,
                'primary_pair': bool(oprim and gprim),
                **oextra,
                **m,
            })
    return rows


def summarize_compat(rows: List[dict]) -> dict:
    if not rows:
        return {
            'compat_row_count': 0,
            'primary_compat_row_count': 0,
            'valid_omega_g_count': 0,
            'primary_valid_omega_g_count': 0,
            'compatibility_gate_pass_count': 0,
            'primary_compatibility_gate_pass_count': 0,
            'best_primary_J2_resid': 0.0,
            'best_primary_lock_mean_resid': 0.0,
            'best_primary_metric_orth_resid': 0.0,
            'best_primary_candidate': '',
            'best_primary_metric': '',
        }
    primary = [r for r in rows if fbool(r.get('primary_pair'))]
    def count(rs, key): return sum(1 for r in rs if fbool(r.get(key)))
    def score(r):
        valid_bonus = 0 if fbool(r.get('compat_valid_omega_g')) else 100
        return (
            valid_bonus,
            float(r.get('compat_J2_plus_I_resid', 99.0)),
            float(r.get('compat_J_QP_mean_lock_resid', 99.0)),
            float(r.get('compat_J_metric_orth_resid', 99.0)),
        )
    best_primary = min(primary, key=score, default=None)
    best_any = min(rows, key=score, default=None)
    def pick(r, k, default=0.0):
        return float(r.get(k, default)) if r else default
    return {
        'compat_row_count': len(rows),
        'primary_compat_row_count': len(primary),
        'valid_omega_g_count': count(rows, 'compat_valid_omega_g'),
        'primary_valid_omega_g_count': count(primary, 'compat_valid_omega_g'),
        'compatibility_gate_pass_count': count(rows, 'compat_compatibility_gate_pass'),
        'primary_compatibility_gate_pass_count': count(primary, 'compat_compatibility_gate_pass'),
        'best_primary_J2_resid': pick(best_primary, 'compat_J2_plus_I_resid'),
        'best_primary_lock_mean_resid': pick(best_primary, 'compat_J_QP_mean_lock_resid'),
        'best_primary_lock_max_resid': pick(best_primary, 'compat_J_QP_max_lock_resid'),
        'best_primary_metric_orth_resid': pick(best_primary, 'compat_J_metric_orth_resid'),
        'best_primary_hash_anti_resid': pick(best_primary, 'compat_J_hash_anti_self_resid'),
        'best_primary_star_span_resid': pick(best_primary, 'compat_star_family_span_resid'),
        'best_primary_omega_nondeg_ratio': pick(best_primary, 'compat_omega_nondeg_ratio'),
        'best_primary_g_nondeg_ratio': pick(best_primary, 'compat_g_nondeg_ratio'),
        'best_primary_candidate': best_primary.get('omega_candidate','') if best_primary else '',
        'best_primary_metric': best_primary.get('metric_candidate','') if best_primary else '',
        'best_any_J2_resid': pick(best_any, 'compat_J2_plus_I_resid'),
        'best_any_lock_mean_resid': pick(best_any, 'compat_J_QP_mean_lock_resid'),
        'best_any_candidate': best_any.get('omega_candidate','') if best_any else '',
        'best_any_metric': best_any.get('metric_candidate','') if best_any else '',
    }


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
    symp_rows, symp_summary = p74.symplectic_rows(model, K, assembly_log, args)
    compat_rows: List[dict] = []
    for i, a in enumerate(assembly_log):
        if not fbool(a.get('assembly_applied')):
            continue
        for r in compatibility_rows_for_assembly(model, K, a, args):
            r['assembly_index'] = i
            compat_rows.append(r)
    compat_summary = summarize_compat(compat_rows)
    write_csv(vout / 'birth_geometry_log.csv', birth_log)
    write_csv(vout / 'assembly_pairing_log.csv', pairing_log)
    write_csv(vout / 'assembly_ablation_log.csv', assembly_log)
    write_csv(vout / 'candidate_eval_rows.csv', candidate_rows)
    write_csv(vout / 'directed_pair_rows.csv', pair_rows)
    write_csv(vout / 'signed_quadrature_rows.csv', signed_rows)
    write_csv(vout / 'alignment_pair_rows.csv', align_pair_rows)
    write_csv(vout / 'assembly_motif_basis_rows.csv', motif_rows)
    write_csv(vout / 'real_symplectic_before_star_rows.csv', symp_rows)
    write_csv(vout / 'kahler_compatibility_star_rows.csv', compat_rows)
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
        'kahler_compatibility_metrics': compat_summary,
        'automatic_pairings_applied': sum(1 for x in pairing_log if fbool(x.get('applied'))),
        'assemblies_applied': sum(1 for x in assembly_log if fbool(x.get('assembly_applied'))),
        'assemblies_attempted': len(assembly_log),
        'decision_used_delta_beta_any': sel['decision_used_delta_beta_any'],
    }
    (vout / 'variant_kahler_compatibility_star_summary.json').write_text(json.dumps(summary, indent=2), encoding='utf-8')
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
    a = r['auto_metrics']; dm = r['directed_metrics']; sm = r['signed_quadrature_metrics']; am = r['alignment_metrics']; mm = r['motif_basis_metrics']; sy = r['real_symplectic_metrics']; kc = r['kahler_compatibility_metrics']
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
        'symp_best_candidate': sy.get('best_primary_candidate',''),
        'symp_best_nondeg_ratio': sy.get('best_primary_nondeg_ratio',0.0),
        'compat_primary_pass': kc.get('primary_compatibility_gate_pass_count',0),
        'compat_valid_omega_g': kc.get('primary_valid_omega_g_count',0),
        'compat_best_J2': kc.get('best_primary_J2_resid',0.0),
        'compat_best_lock_mean': kc.get('best_primary_lock_mean_resid',0.0),
        'compat_best_lock_max': kc.get('best_primary_lock_max_resid',0.0),
        'compat_best_metric_orth': kc.get('best_primary_metric_orth_resid',0.0),
        'compat_best_hash_anti': kc.get('best_primary_hash_anti_resid',0.0),
        'compat_best_star_span': kc.get('best_primary_star_span_resid',0.0),
        'compat_best_omega': kc.get('best_primary_candidate',''),
        'compat_best_metric': kc.get('best_primary_metric',''),
        'used_delta_beta': r['decision_used_delta_beta_any'],
    }


def write_comparative(out: Path, rows: List[dict]) -> None:
    flat = []
    for r in rows:
        s = slim(r)
        kc = r['kahler_compatibility_metrics']
        flat.append({**s,
            'beta0': s['beta'][0], 'beta1': s['beta'][1], 'beta2': s['beta'][2], 'beta3': s['beta'][3],
            **{k: kc.get(k, '') for k in sorted(kc.keys()) if not isinstance(kc.get(k), (dict, list))}
        })
    write_csv(out / 'comparative_kahler_compatibility_star_summary.csv', flat)


def make_docs(summary: dict) -> Tuple[str, str, str, str]:
    rows = summary['variant_rows']
    lines = ['| option | variant | beta | pairs | asm | pair harm | Q harm | P harm | symp pass | best Ω | Ω ratio | compat pass | valid Ω-g | best Ω/g | J2 resid | QP lock | metric orth | #J anti | star span | used dβ? |',
             '|---|---|---:|---:|---:|---:|---:|---:|---:|---|---:|---:|---:|---|---:|---:|---:|---:|---:|---:|']
    for r in rows:
        s = slim(r)
        lines.append(f"| {s['option']} | {s['variant']} | ({s['beta'][0]},{s['beta'][1]},{s['beta'][2]},{s['beta'][3]}) | {s['pairings']} | {s['assemblies']} | {s['pair_harm']:.6g} | {s['Q_harm']:.6g} | {s['P_harm']:.6g} | {s['symp_primary_pass']} | {s['symp_best_candidate']} | {s['symp_best_nondeg_ratio']:.6g} | {s['compat_primary_pass']} | {s['compat_valid_omega_g']} | {s['compat_best_omega']} / {s['compat_best_metric']} | {s['compat_best_J2']:.6g} | {s['compat_best_lock_mean']:.6g} | {s['compat_best_metric_orth']:.6g} | {s['compat_best_hash_anti']:.6g} | {s['compat_best_star_span']:.6g} | {s['used_delta_beta']} |")
    table = '\n'.join(lines)
    nonstrict = [r for r in rows if r['variant'] != 'strict_symmetrized_control']
    any_pass = any(r['kahler_compatibility_metrics'].get('primary_compatibility_gate_pass_count',0) > 0 for r in nonstrict)
    smd = f"""# SUMMARY — Kähler compatibility / star-from-symplectic gate

Model tag: `CQNM/s=-1 saturated geometry reference, provenance-growth L2 diagnostic`.

This package corrects the previous symplectic test: a nondegenerate skew form alone is not counted as a breakthrough.  The gate now asks whether a primary Ω candidate and a primary g candidate are simultaneously compatible enough that

```text
J = g^-1 Ω
```

has small `J²+I`, maps the derived Q-subspace to the derived P-subspace, is metric-compatible, and yields a stable metric-adjoint `#` on the tested operator family.

{table}

Decision:

```json
{{
  "any_non_strict_primary_compatibility_pass": {str(any_pass).lower()}
}}
```
"""
    rmd = f"""# RESULTS — Kähler compatibility / star-from-symplectic gate

## Comparative table

{table}

## Gate definition

A row only passes if all of the following hold:

```text
Ω is skew and nondegenerate on the actual Q/P motif span,
g is symmetric and nondegenerate on the same span,
J = g^-1 Ω satisfies J² ≈ -I,
J maps span(Q_motif) to span(P_motif),
J is g-orthogonal enough,
J# ≈ -J for the metric-adjoint #,
the tested operator family is # stable,
strict_sym remains null,
used_delta_beta remains false.
```

The test deliberately logs, but does not reward, the tautological facts that a metric adjoint is involutive and anti-multiplicative for invertible g.  The success condition is simultaneous compatibility of Ω, g, #, and the Q/P split.

## Interpretation rule

- If `compat pass > 0`, the ladder may proceed from real Ω to a candidate real #/*-algebra and then to J.
- If Ω passes but compatibility fails, the result is a Kähler-like compatibility obstruction: the ingredients exist separately but do not align.
- If strict_sym is nonzero, the result is invalid as a provenance asymmetry diagnostic.
"""
    audit = """# SOURCE AUDIT

No i, global J, Hodge star, physical adjoint, positivity, C*-norm, final sym(M), or delta-beta/H2 decision is introduced.

Important limitation: the primary Ω candidates include skew parts of earlier J-like pair/union/interface/link operators.  They are therefore not independent derivations of Ω.  This package treats them as structural candidates and tests the stronger simultaneous compatibility with g and #, precisely to avoid overclaiming from isolated skew nondegeneracy.

Metric candidates:
- `cochain_identity_metric`: ambient coefficient metric already used by residual diagnostics.
- `union_C_symmetric_metric`: symmetric part of the derived pair-conjugation C on the motif union.
- `union_C_square_metric`: C^T C structural diagnostic.
- `QP_data_gram_control`: non-primary control using Q/P data directly.

This is a Python L2 diagnostic, not a Lean theorem.
"""
    readme = """# Kähler compatibility / star-from-symplectic gate

Run:

```bash
python3 test_kahler_compatibility_star_gate.py
```

The test checks simultaneous compatibility of Ω, g, and a metric-adjoint # before any J/i claim is made.
"""
    return smd, rmd, audit, readme


def package(out: Path, zip_path: Path) -> None:
    files = [
        Path(__file__).name,
        'test_real_symplectic_before_star_gate.py',
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
    ap.add_argument('--sym-threshold', type=float, default=1e-8)
    ap.add_argument('--nondeg-threshold', type=float, default=1e-6)
    ap.add_argument('--metric-regularizer', type=float, default=1e-9)
    ap.add_argument('--compat-J2-threshold', type=float, default=0.20)
    ap.add_argument('--compat-lock-mean-threshold', type=float, default=0.20)
    ap.add_argument('--compat-lock-max-threshold', type=float, default=0.30)
    ap.add_argument('--metric-orth-threshold', type=float, default=0.25)
    ap.add_argument('--hash-anti-threshold', type=float, default=0.25)
    ap.add_argument('--star-span-threshold', type=float, default=0.25)
    # Args required by p74.symplectic_rows.
    ap.add_argument('--symplectic-ratio-threshold', type=float, default=1e-3)
    ap.add_argument('--qp-ratio-threshold', type=float, default=1e-3)
    ap.add_argument('--isotropic-threshold', type=float, default=0.35)
    ap.add_argument('--out', default='kahler_compatibility_star_out_L2')
    ap.add_argument('--zip', default='cnna_kahler_compatibility_star_gate_pkg_L2.zip')
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
