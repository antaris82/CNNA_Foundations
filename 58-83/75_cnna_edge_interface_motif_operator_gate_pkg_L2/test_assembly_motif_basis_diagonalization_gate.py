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
import test_harmonic_k_orientation_kappa_gate as hk
import test_pairing_transport_antisym_birth_coherence_gate as p56
import test_pairing_quadrature_split_symplectic_defect_gate as p58
import test_pairing_quadrature_adjoint_pairing_gate as p60
import test_pair_J_alignment_search_gate as p61
import test_dual_pairing_assembly_growth_rule_gate as p68
import test_dual_assembly_order_context_ablation_gate as p69
import test_signed_Jlock_role_coupling_gate as p70

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
    return p68.fbool(x)


def norm(x: np.ndarray) -> float:
    return float(np.linalg.norm(x))


def orth_basis(cols: List[np.ndarray], tol: float = 1e-10) -> np.ndarray:
    if not cols:
        return np.zeros((0, 0), dtype=float)
    M = np.column_stack([np.asarray(c, dtype=float).reshape(-1) for c in cols])
    if M.size == 0:
        return np.zeros((M.shape[0], 0), dtype=float)
    U, S, _ = np.linalg.svd(M, full_matrices=False)
    keep = S > tol * max(1.0, float(S[0]) if S.size else 1.0)
    return U[:, keep]


def projector(U: np.ndarray) -> np.ndarray:
    if U.size == 0 or U.shape[1] == 0:
        return np.zeros((U.shape[0], U.shape[0]), dtype=float)
    return U @ U.T


def subspace_image_residual(J: np.ndarray, source_cols: List[np.ndarray], target_cols: List[np.ndarray]) -> float:
    S = orth_basis(source_cols)
    T = orth_basis(target_cols)
    if S.shape[1] == 0:
        return 0.0
    JS = J @ S
    den = norm(JS)
    if den < EPS:
        return 0.0
    PT = projector(T)
    return norm((np.eye(J.shape[0]) - PT) @ JS) / (den + EPS)


def span_leakage(J: np.ndarray, cols: List[np.ndarray]) -> float:
    U = orth_basis(cols)
    if U.shape[1] == 0:
        return 0.0
    JU = J @ U
    den = norm(JU)
    if den < EPS:
        return 0.0
    PU = projector(U)
    return norm((np.eye(J.shape[0]) - PU) @ JU) / (den + EPS)


def projected_J2_residual(J: np.ndarray, cols: List[np.ndarray]) -> float:
    U = orth_basis(cols)
    r = U.shape[1]
    if r == 0:
        return 0.0
    Jp = U.T @ J @ U
    return norm(Jp @ Jp + np.eye(r)) / (norm(np.eye(r)) + EPS)


def projected_anticomm_residual(C: np.ndarray, J: np.ndarray, cols: List[np.ndarray]) -> float:
    U = orth_basis(cols)
    r = U.shape[1]
    if r == 0:
        return 0.0
    Cp = U.T @ C @ U
    Jp = U.T @ J @ U
    return norm(Cp @ Jp @ Cp + Jp) / (norm(Jp) + EPS)


def skew_residual(J: np.ndarray, cols: List[np.ndarray]) -> float:
    U = orth_basis(cols)
    if U.shape[1] == 0:
        return 0.0
    Jp = U.T @ J @ U
    return norm(Jp + Jp.T) / (norm(Jp) + EPS)


def parse_face_from_field(x) -> Optional[Face]:
    try:
        f = core.parse_face_string(str(x))
        if len(f) == 3:
            return tuple(sorted(f))
    except Exception:
        return None
    return None


def parse_faces_from_assembly(row: dict, prefix: str) -> Tuple[Optional[Face], Optional[Face]]:
    fa = parse_face_from_field(row.get(f'{prefix}_face_a', ''))
    fb = parse_face_from_field(row.get(f'{prefix}_face_b', ''))
    return fa, fb


def pair_fields(model, fa: Face, fb: Face, args: argparse.Namespace) -> Optional[dict]:
    try:
        ka = p56.axial(p56.face_K_directed(model, fa, args.source, args.phase_sign, args.antisym_eta, args.erase_phase_for_strict_sym))
        kb = p56.axial(p56.face_K_directed(model, fb, args.source, args.phase_sign, args.antisym_eta, args.erase_phase_for_strict_sym))
        na_out = hk.face_normal(model, fa, 'outward')
        nb_out = hk.face_normal(model, fb, 'outward')
        na_birth = hk.face_normal(model, fa, 'birth_order')
        R = p56.rotation_from_to(nb_out, -na_out)
        kb_r = R @ kb
        Q_a = ka + kb_r
        P_a = ka - kb_r
        Q_b = -(R.T @ Q_a)
        P_b = +(R.T @ P_a)
        C_a = np.cross(Q_a, P_a)
        C_b = -(R.T @ C_a)
        return {
            'fa': fa, 'fb': fb, 'ka': ka, 'kb': kb, 'kb_r': kb_r, 'R': R,
            'Q_a': Q_a, 'P_a': P_a, 'Q_b': Q_b, 'P_b': P_b,
            'C_a': C_a, 'C_b': C_b,
            'J': p60.block_J(R), 'C': p60.block_C(R),
            'na_birth': na_birth, 'na_out': na_out,
        }
    except Exception:
        return None


def put_pair_vec_union(faces: List[Face], p: dict, kind: str) -> np.ndarray:
    idx = {f: i for i, f in enumerate(faces)}
    v = np.zeros((len(faces), 3), dtype=float)
    if kind == 'Q':
        va, vb = p['Q_a'], p['Q_b']
    elif kind == 'P':
        va, vb = p['P_a'], p['P_b']
    elif kind == 'C':
        va, vb = p['C_a'], p['C_b']
    elif kind == 'raw':
        va, vb = p['ka'], p['kb']
    else:
        raise ValueError(kind)
    v[idx[p['fa']]] += va
    v[idx[p['fb']]] += vb
    return v.reshape(-1)


def union_JC(faces: List[Face], pairs: List[dict]) -> Tuple[np.ndarray, np.ndarray]:
    n = 3 * len(faces)
    J = np.zeros((n, n), dtype=float)
    C = np.zeros((n, n), dtype=float)
    idx = {f: i for i, f in enumerate(faces)}
    def sl(f: Face):
        i = idx[f]
        return slice(3*i, 3*i+3)
    for p in pairs:
        fa, fb, R = p['fa'], p['fb'], p['R']
        sa, sb = sl(fa), sl(fb)
        J[sa, sb] += R
        J[sb, sa] += -R.T
        C[sa, sb] += R
        C[sb, sa] += R.T
    return J, C


def direct_sum_motif(pA: dict, pB: dict) -> Tuple[np.ndarray, np.ndarray, Dict[str, np.ndarray]]:
    J = np.zeros((12, 12), dtype=float)
    C = np.zeros((12, 12), dtype=float)
    J[:6,:6] = pA['J']; J[6:,6:] = pB['J']
    C[:6,:6] = pA['C']; C[6:,6:] = pB['C']
    z = np.zeros(6, dtype=float)
    v = {
        'A_Q': np.concatenate([np.concatenate([pA['Q_a'], pA['Q_b']]), z]),
        'A_P': np.concatenate([np.concatenate([pA['P_a'], pA['P_b']]), z]),
        'B_Q': np.concatenate([z, np.concatenate([pB['Q_a'], pB['Q_b']])]),
        'B_P': np.concatenate([z, np.concatenate([pB['P_a'], pB['P_b']])]),
    }
    return J, C, v


def motif_metrics(prefix: str, J: np.ndarray, C: np.ndarray, vecs: Dict[str, np.ndarray]) -> dict:
    Qcols = [vecs['A_Q'], vecs['B_Q']]
    Pcols = [vecs['A_P'], vecs['B_P']]
    allcols = Qcols + Pcols
    q_to_p = subspace_image_residual(J, Qcols, Pcols)
    p_to_q = subspace_image_residual(J, Pcols, Qcols)
    leak = span_leakage(J, allcols)
    j2 = projected_J2_residual(J, allcols)
    cjc = projected_anticomm_residual(C, J, allcols)
    skew = skew_residual(J, allcols)
    U = orth_basis(allcols)
    Q = orth_basis(Qcols); P = orth_basis(Pcols)
    mean_lock = 0.5 * (q_to_p + p_to_q)
    max_lock = max(q_to_p, p_to_q)
    return {
        f'{prefix}_basis_rank': int(U.shape[1]),
        f'{prefix}_Q_rank': int(Q.shape[1]),
        f'{prefix}_P_rank': int(P.shape[1]),
        f'{prefix}_J_Q_to_P_subspace_resid': q_to_p,
        f'{prefix}_J_P_to_Q_subspace_resid': p_to_q,
        f'{prefix}_J_QP_subspace_mean_resid': mean_lock,
        f'{prefix}_J_QP_subspace_max_resid': max_lock,
        f'{prefix}_J_span_leakage': leak,
        f'{prefix}_projected_J2_plus_I_resid': j2,
        f'{prefix}_projected_CJC_plus_J_resid': cjc,
        f'{prefix}_projected_J_skew_resid': skew,
        f'{prefix}_Q_norm': norm(np.column_stack(Qcols)) if Qcols else 0.0,
        f'{prefix}_P_norm': norm(np.column_stack(Pcols)) if Pcols else 0.0,
        f'{prefix}_gate_pass': bool(mean_lock < 0.20 and max_lock < 0.30 and leak < 0.30 and j2 < 0.40),
    }


def signed_motif_stats(pA: dict, pB: dict) -> dict:
    rows = {}
    for label, p in [('A', pA), ('B', pB)]:
        cross = np.cross(p['Q_a'], p['P_a'])
        a = norm(cross)
        s = float(np.dot(cross, p['na_birth']))
        rows[f'{label}_signed_birth_over_abs'] = s / (a + EPS)
        rows[f'{label}_area_abs'] = a
    total_abs = rows['A_area_abs'] + rows['B_area_abs'] + EPS
    rows['motif_signed_birth_over_abs_weighted'] = (
        rows['A_signed_birth_over_abs'] * rows['A_area_abs'] + rows['B_signed_birth_over_abs'] * rows['B_area_abs']
    ) / total_abs
    return rows


def assembly_motif_rows(model, K, pairing_log: List[dict], assembly_log: List[dict], args: argparse.Namespace) -> Tuple[List[dict], dict]:
    rows: List[dict] = []
    for i, a in enumerate(assembly_log):
        if not fbool(a.get('assembly_applied')):
            continue
        A1, A2 = parse_faces_from_assembly(a, 'A')
        B1, B2 = parse_faces_from_assembly(a, 'B')
        if A1 is None or A2 is None or B1 is None or B2 is None:
            continue
        pA = pair_fields(model, A1, A2, args)
        pB = pair_fields(model, B1, B2, args)
        if pA is None or pB is None:
            continue
        Jds, Cds, Vds = direct_sum_motif(pA, pB)
        ds = motif_metrics('direct_sum', Jds, Cds, Vds)
        uf = sorted({A1, A2, B1, B2})
        Ju, Cu = union_JC(uf, [pA, pB])
        Vu = {
            'A_Q': put_pair_vec_union(uf, pA, 'Q'),
            'A_P': put_pair_vec_union(uf, pA, 'P'),
            'B_Q': put_pair_vec_union(uf, pB, 'Q'),
            'B_P': put_pair_vec_union(uf, pB, 'P'),
        }
        um = motif_metrics('union_sum', Ju, Cu, Vu)
        singA = p61.J_lock_residual(pA['J'], np.concatenate([pA['Q_a'], pA['Q_b']]), np.concatenate([pA['P_a'], pA['P_b']]))
        singB = p61.J_lock_residual(pB['J'], np.concatenate([pB['Q_a'], pB['Q_b']]), np.concatenate([pB['P_a'], pB['P_b']]))
        shared_faces = len({A1, A2} & {B1, B2})
        shared_vertices = len((set(A1)|set(A2)) & (set(B1)|set(B2)))
        row = {
            'assembly_index': i,
            'context': a.get('context',''),
            'shared_faces': shared_faces,
            'shared_vertices': shared_vertices,
            'A_face_a': str(list(A1)), 'A_face_b': str(list(A2)),
            'B_face_a': str(list(B1)), 'B_face_b': str(list(B2)),
            'A_pair_J_lock_mean_raw_QP': singA['J_lock_mean_resid'],
            'B_pair_J_lock_mean_raw_QP': singB['J_lock_mean_resid'],
            'pair_local_mean_J_lock_raw_QP': 0.5*(singA['J_lock_mean_resid'] + singB['J_lock_mean_resid']),
            'pair_local_max_J_lock_raw_QP': max(singA['J_lock_max_resid'], singB['J_lock_max_resid']),
            **signed_motif_stats(pA, pB),
            **ds,
            **um,
        }
        row['motif_improves_over_pairlocal_direct_sum'] = bool(row['direct_sum_J_QP_subspace_mean_resid'] < row['pair_local_mean_J_lock_raw_QP'])
        row['motif_improves_over_pairlocal_union_sum'] = bool(row['union_sum_J_QP_subspace_mean_resid'] < row['pair_local_mean_J_lock_raw_QP'])
        rows.append(row)
    if not rows:
        return rows, {
            'assembly_count': 0,
            'direct_sum_best_mean_resid': 0.0,
            'union_sum_best_mean_resid': 0.0,
            'direct_sum_gate_pass_count': 0,
            'union_sum_gate_pass_count': 0,
        }
    def avg(k): return float(np.mean([float(r[k]) for r in rows]))
    def mn(k): return float(np.min([float(r[k]) for r in rows]))
    summary = {
        'assembly_count': len(rows),
        'direct_sum_best_mean_resid': mn('direct_sum_J_QP_subspace_mean_resid'),
        'direct_sum_avg_mean_resid': avg('direct_sum_J_QP_subspace_mean_resid'),
        'direct_sum_best_max_resid': mn('direct_sum_J_QP_subspace_max_resid'),
        'direct_sum_avg_span_leakage': avg('direct_sum_J_span_leakage'),
        'direct_sum_avg_J2_resid': avg('direct_sum_projected_J2_plus_I_resid'),
        'direct_sum_gate_pass_count': sum(1 for r in rows if fbool(r['direct_sum_gate_pass'])),
        'union_sum_best_mean_resid': mn('union_sum_J_QP_subspace_mean_resid'),
        'union_sum_avg_mean_resid': avg('union_sum_J_QP_subspace_mean_resid'),
        'union_sum_best_max_resid': mn('union_sum_J_QP_subspace_max_resid'),
        'union_sum_avg_span_leakage': avg('union_sum_J_span_leakage'),
        'union_sum_avg_J2_resid': avg('union_sum_projected_J2_plus_I_resid'),
        'union_sum_gate_pass_count': sum(1 for r in rows if fbool(r['union_sum_gate_pass'])),
        'pair_local_avg_mean_J_lock_raw_QP': avg('pair_local_mean_J_lock_raw_QP'),
        'motif_improves_direct_sum_count': sum(1 for r in rows if fbool(r['motif_improves_over_pairlocal_direct_sum'])),
        'motif_improves_union_sum_count': sum(1 for r in rows if fbool(r['motif_improves_over_pairlocal_union_sum'])),
        'signed_birth_weighted_avg': avg('motif_signed_birth_over_abs_weighted'),
    }
    return rows, summary


def run_variant_option_with_motif(variant: str, args: argparse.Namespace, out: Path) -> dict:
    model = nl.build_model(variant, args)
    model.grow(args.max_level)
    baseline_K = core.build_dynamic_outward_ngf_complex(model)
    baseline_metrics = core.full_metrics(model, baseline_K, args.source)
    tag = p69.option_tag(args)
    vout = out / tag / variant
    vout.mkdir(parents=True, exist_ok=True)
    K, birth_log, pairing_log, assembly_log, candidate_rows = p69.build_ablation_complex(model, args, variant, vout)
    auto_metrics = core.full_metrics(model, K, args.source)
    dm, pair_rows, top_rows, three_rows = p56.directed_metrics(model, K, pairing_log, args)
    sm, signed_rows, signed_face_rows = p69.p59.signed_quadrature_rows(model, K, pairing_log, args)
    am, align_pair_rows, align_candidate_rows, align_candidate_summary = p61.alignment_search_metrics(model, K, pairing_log, args)
    sel = p69.summarize_selection(pairing_log)
    motif_rows, motif_summary = assembly_motif_rows(model, K, pairing_log, assembly_log, args)
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
        'motif_basis_metrics': motif_summary,
        'automatic_pairings_applied': sum(1 for x in pairing_log if fbool(x.get('applied'))),
        'assemblies_applied': sum(1 for x in assembly_log if fbool(x.get('assembly_applied'))),
        'assemblies_attempted': len(assembly_log),
        'decision_used_delta_beta_any': sel['decision_used_delta_beta_any'],
    }
    (vout / 'variant_option_motif_summary.json').write_text(json.dumps(summary, indent=2), encoding='utf-8')
    return summary


def clone_args(args: argparse.Namespace, **updates) -> argparse.Namespace:
    d = vars(args).copy(); d.update(updates); return argparse.Namespace(**d)


def run_experiment(args: argparse.Namespace, out: Path) -> List[dict]:
    p70.patch_modules()
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
    rows = []
    for opt in options:
        for variant in args.variants:
            rows.append(run_variant_option_with_motif(variant, opt, out))
    return rows


def slim(r: dict) -> dict:
    a = r['auto_metrics']; dm = r['directed_metrics']; sm = r['signed_quadrature_metrics']; am = r['alignment_metrics']; mm = r['motif_basis_metrics']
    return {
        'option': r['option'], 'variant': r['variant'],
        'beta': [a['beta0'], a['beta1'], a['beta2'], a['beta3']],
        'pairings': r['automatic_pairings_applied'],
        'assemblies': r['assemblies_applied'],
        'pair_harm': dm['pair_transport_harmonic_ratio'],
        'Q_harm': am['Q_even_harmonic_ratio'],
        'P_harm': am['P_odd_harmonic_ratio'],
        'pair_local_J_lock': am['best_per_pair_mean_J_lock_resid'],
        'signed_birth': sm['signed_birth_over_abs_sum_ratio'],
        'motif_count': mm.get('assembly_count',0),
        'direct_sum_motif_lock': mm.get('direct_sum_best_mean_resid',0.0),
        'union_sum_motif_lock': mm.get('union_sum_best_mean_resid',0.0),
        'direct_sum_gate_pass': mm.get('direct_sum_gate_pass_count',0),
        'union_sum_gate_pass': mm.get('union_sum_gate_pass_count',0),
        'used_delta_beta': r['decision_used_delta_beta_any'],
    }


def write_comparative(out: Path, rows: List[dict]) -> None:
    flat = []
    for r in rows:
        s = slim(r); mm = r['motif_basis_metrics']
        flat.append({
            **s,
            'beta0': s['beta'][0], 'beta1': s['beta'][1], 'beta2': s['beta'][2], 'beta3': s['beta'][3],
            'direct_sum_avg_mean_resid': mm.get('direct_sum_avg_mean_resid',0.0),
            'direct_sum_avg_span_leakage': mm.get('direct_sum_avg_span_leakage',0.0),
            'direct_sum_avg_J2_resid': mm.get('direct_sum_avg_J2_resid',0.0),
            'union_sum_avg_mean_resid': mm.get('union_sum_avg_mean_resid',0.0),
            'union_sum_avg_span_leakage': mm.get('union_sum_avg_span_leakage',0.0),
            'union_sum_avg_J2_resid': mm.get('union_sum_avg_J2_resid',0.0),
            'pair_local_avg_mean_J_lock_raw_QP': mm.get('pair_local_avg_mean_J_lock_raw_QP',0.0),
            'motif_improves_direct_sum_count': mm.get('motif_improves_direct_sum_count',0),
            'motif_improves_union_sum_count': mm.get('motif_improves_union_sum_count',0),
        })
    write_csv(out / 'comparative_assembly_motif_basis_diagonalization_summary.csv', flat)


def make_docs(summary: dict) -> Tuple[str,str,str,str]:
    rows = summary['variant_rows']
    lines = ['| option | variant | beta | pairs | asm | pair harm | Q harm | P harm | pair J-lock | signed | motif n | direct motif lock | union motif lock | direct pass | union pass | used dBeta? |',
             '|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|']
    for r in rows:
        s = slim(r)
        lines.append(f"| {s['option']} | {s['variant']} | ({s['beta'][0]},{s['beta'][1]},{s['beta'][2]},{s['beta'][3]}) | {s['pairings']} | {s['assemblies']} | {s['pair_harm']:.6g} | {s['Q_harm']:.6g} | {s['P_harm']:.6g} | {s['pair_local_J_lock']:.6g} | {s['signed_birth']:.6g} | {s['motif_count']} | {s['direct_sum_motif_lock']:.6g} | {s['union_sum_motif_lock']:.6g} | {s['direct_sum_gate_pass']} | {s['union_sum_gate_pass']} | {s['used_delta_beta']} |")
    table = '\n'.join(lines)
    nonstrict = [r for r in rows if r['variant'] != 'strict_symmetrized_control']
    best_direct = min(nonstrict, key=lambda r: r['motif_basis_metrics'].get('direct_sum_best_mean_resid',9.0), default=None)
    best_union = min(nonstrict, key=lambda r: r['motif_basis_metrics'].get('union_sum_best_mean_resid',9.0), default=None)
    best = {'best_direct_sum': slim(best_direct) if best_direct else {}, 'best_union_sum': slim(best_union) if best_union else {}}
    smd = f"""# SUMMARY — assembly motif basis diagonalization gate

Model label:
CNNA growing primal simplicial complex with deterministic sequential provenance growth, directed antisymmetric birth-transport operators, and signed-Jlock two-pair assembly motifs.

This test asks whether the previous J-lock obstruction is a measurement-basis artifact.  Instead of testing each pair separately, each complete A/B assembly is represented by the native motif basis

```text
A_Q, A_P, B_Q, B_P
```

Two derived motif operators are tested:

```text
direct_sum: pair-local block-diagonal J/C on A and B separately
union_sum: shared-face motif operator obtained by summing the two pair-exchange maps on the actual union of faces
```

A positive result would require J to map the motif Q-plane span(A_Q,B_Q) into the motif P-plane span(A_P,B_P), with low span leakage and reasonable projected J^2 = -I behavior.  No i, global J, Hodge, *, positivity, C*-norm, final sym(M), or delta-beta decision is introduced.

{table}

## Best rows

```json
{json.dumps(best, indent=2)}
```
"""
    rmd = f"""# RESULTS — assembly motif basis diagonalization gate

## Comparative table

{table}

## Gate criterion

The motif basis gate passes only if the combined A/B motif, not merely an individual pair, shows:

```text
J(span(A_Q,B_Q)) approximately subset span(A_P,B_P),
J(span(A_P,B_P)) approximately subset span(A_Q,B_Q),
low J-span leakage,
projected J^2 approximately -I,
strict_sym killed,
used_delta_beta = false.
```

This is a subspace-lock diagnostic.  It permits derived A/B mixing inside the motif planes, but it does not fit an arbitrary rotation and does not define P as JQ.
"""
    audit = """# SOURCE AUDIT

Carried forward:

- Single-pair gates found Q/P channels and local C/J pair algebra, but no robust dynamic J_pair(Q)=P lock.
- Tradeoff and kappa-permutation gates showed beta2, C-lock, Q/P, and signed kappa flip split across single-pair candidates.
- Dual-pair assembly gates showed beta2/QP/signed features can be carried by a two-pair motif.
- signed-Jlock role-coupling strengthened beta2/QP/signed magnitudes but still failed pair-local J-lock.

This package changes only the measurement basis: it tests complete A/B motifs.  It does not set i, J, Hodge, *, positivity, a C*-norm, or final sym(M).  Delta-beta/H2/harmonic data are measured after the fact only.

Caveat: the union_sum operator is a derived diagnostic obtained by summing local pair-exchange maps on the shared-face motif.  It is not yet a proven CNNA dynamic law; if it succeeds, it would still need formal derived-only justification.
"""
    readme = """# Assembly motif basis diagonalization gate

Run:

```bash
python3 test_assembly_motif_basis_diagonalization_gate.py
```

The script evaluates the complete two-pair motif basis A_Q, A_P, B_Q, B_P for dynamic A/B assemblies and writes CSV/JSON summaries plus this ZIP package.
"""
    return smd, rmd, audit, readme


def package(out: Path, zip_path: Path) -> None:
    files = [
        Path(__file__).name,
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
    ap.add_argument('--out', default='assembly_motif_basis_diagonalization_out_L2')
    ap.add_argument('--zip', default='cnna_assembly_motif_basis_diagonalization_gate_pkg_L2.zip')
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
