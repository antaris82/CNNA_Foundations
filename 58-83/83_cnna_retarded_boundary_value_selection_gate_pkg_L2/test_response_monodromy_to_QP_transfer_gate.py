#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import math
import shutil
import zipfile
from pathlib import Path
from typing import Dict, List, Tuple, Optional

import numpy as np

import cnna_non_shelling_core as core
import test_nonlinear_asymmetry_cascade_growth as nl
import test_dual_assembly_order_context_ablation_gate as p69
import test_signed_Jlock_role_coupling_gate as p70
import test_assembly_motif_basis_diagonalization_gate as p71
import test_edge_interface_motif_operator_gate as p72
import test_real_symplectic_before_star_gate as p74
import test_kahler_compatibility_star_gate as p75
import test_pairing_transport_antisym_birth_coherence_gate as p56
import test_signed_quadrature_area_kappa_gate as p59
import test_pair_J_alignment_search_gate as p61

EPS = 1e-12
Face = Tuple[int, int, int]


def fbool(x) -> bool:
    if isinstance(x, str):
        return x.lower() in {'true','1','yes'}
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


def unit(v: np.ndarray) -> np.ndarray:
    v = np.asarray(v, dtype=float).reshape(3)
    n = norm(v)
    return v / (n + EPS) if n > EPS else np.zeros(3)


def birth_q_vector(model: core.DynamicProvenanceGrowth, node_id: int) -> np.ndarray:
    n = model.nodes[int(node_id)]
    phase = 2.0 * math.pi * (int(n.birth_order) - 1) / 3.0 if n.parent is not None else 0.0
    return unit(math.cos(phase) * n.e1 + math.sin(phase) * n.e2)


def birth_h_vector(model: core.DynamicProvenanceGrowth, node_id: int) -> np.ndarray:
    n = model.nodes[int(node_id)]
    return unit(0.7 * n.radial + 0.3 * birth_q_vector(model, node_id))


def orth(cols: List[np.ndarray], tol: float = 1e-10) -> np.ndarray:
    return p71.orth_basis(cols, tol=tol)


def projector(U: np.ndarray) -> np.ndarray:
    if U.size == 0 or U.shape[1] == 0:
        return np.zeros((U.shape[0], U.shape[0]), dtype=float)
    return U @ U.T


def subspace_image_residual(T: np.ndarray, source_cols: List[np.ndarray], target_cols: List[np.ndarray]) -> float:
    S = orth(source_cols)
    R = orth(target_cols)
    if S.shape[1] == 0:
        return 0.0
    TS = T @ S
    den = norm(TS)
    if den < EPS:
        return 0.0
    PR = projector(R)
    return norm((np.eye(T.shape[0]) - PR) @ TS) / (den + EPS)


def projected_J2_residual(T: np.ndarray, cols: List[np.ndarray]) -> float:
    U = orth(cols)
    if U.shape[1] == 0:
        return 0.0
    A = U.T @ T @ U
    I = np.eye(A.shape[0])
    return norm(A @ A + I) / (norm(I) + EPS)


def residual_to_sign(A: np.ndarray, sign: int) -> float:
    if A.size == 0:
        return 1.0
    I = np.eye(A.shape[0])
    return norm(A - sign * I) / (norm(A) + norm(I) + EPS)


def matrix_logsafe(x: float) -> float:
    return math.log(max(float(x), EPS))


def fan_weights(model: core.DynamicProvenanceGrowth, parent: int) -> Optional[dict]:
    children = model.child_ids_ordered(parent)
    if len(children) != 3:
        return None
    W = np.zeros((3,3), dtype=float)
    for i, a in enumerate(children):
        for j, b in enumerate(children):
            if i == j:
                continue
            W[i,j] = float(model.directed_edges.get((a,b), 0.0))
    # Add tiny EPS only for logs/normalization; raw norm remains reported separately.
    Wp = W + EPS * (np.ones_like(W) - np.eye(3))
    circ = (matrix_logsafe(Wp[0,1]) + matrix_logsafe(Wp[1,2]) + matrix_logsafe(Wp[2,0])
            - matrix_logsafe(Wp[1,0]) - matrix_logsafe(Wp[2,1]) - matrix_logsafe(Wp[0,2]))
    row = Wp.copy()
    rs = row.sum(axis=1, keepdims=True) + EPS
    P = row / rs
    A = 0.5 * (P - P.T)
    # weighted forward/reverse cycle permutation; columns are source coordinates.
    F = np.zeros((3,3), dtype=float)
    R = np.zeros((3,3), dtype=float)
    # Normalize edge weights by geometric mean to keep scale diagnostic but not explosive.
    fw = [Wp[0,1], Wp[1,2], Wp[2,0]]
    rv = [Wp[0,2], Wp[2,1], Wp[1,0]]
    gf = math.exp(sum(matrix_logsafe(x) for x in fw) / 3.0)
    gr = math.exp(sum(matrix_logsafe(x) for x in rv) / 3.0)
    F[1,0] = Wp[0,1] / (gf + EPS); F[2,1] = Wp[1,2] / (gf + EPS); F[0,2] = Wp[2,0] / (gf + EPS)
    R[2,0] = Wp[0,2] / (gr + EPS); R[1,2] = Wp[2,1] / (gr + EPS); R[0,1] = Wp[1,0] / (gr + EPS)
    # orientation-only C3 skew generator; magnitude from tanh(circ), no imported J.
    C = np.array([[0.0, -1.0, 1.0], [1.0, 0.0, -1.0], [-1.0, 1.0, 0.0]], dtype=float)
    C = math.tanh(0.25 * circ) * C / math.sqrt(3.0)
    eigP = np.linalg.eigvals(P)
    return {
        'parent': int(parent), 'children': children, 'W': W, 'P': P, 'A': A, 'F': F, 'R': R,
        'cycle_difference': 0.5 * (F - R), 'cycle_sum': 0.5 * (F + R), 'signed_C3_skew': C,
        'circulation_log_ratio': float(circ),
        'raw_weight_norm': norm(W),
        'complex_eig_imag_max': float(max(abs(np.imag(z)) for z in eigP)) if eigP.size else 0.0,
        'complex_eig_real_spread': float(np.ptp(np.real(eigP))) if eigP.size else 0.0,
        'markov_eigvals_json': json.dumps([[float(np.real(z)), float(np.imag(z))] for z in eigP]),
    }


def response_op_from_fan(fan: dict, mode: str) -> np.ndarray:
    return {
        'markov': fan['P'].T,  # acts on column child-coordinate states
        'skew_markov': fan['A'],
        'forward_cycle': fan['F'],
        'reverse_cycle': fan['R'],
        'cycle_difference': fan['cycle_difference'],
        'cycle_sum': fan['cycle_sum'],
        'signed_C3_skew': fan['signed_C3_skew'],
    }[mode]


def motif_parent_fans(model: core.DynamicProvenanceGrowth, faces: List[Face], args: argparse.Namespace) -> List[int]:
    verts = sorted({v for f in faces for v in f})
    by_parent: Dict[int, set[int]] = {}
    for v in verts:
        p = model.nodes[int(v)].parent
        if p is None:
            continue
        if len(model.nodes[p].children) != 3:
            continue
        by_parent.setdefault(int(p), set()).add(int(v))
    out = []
    for p, seen in sorted(by_parent.items()):
        min_seen = args.min_children_in_motif_for_fan
        if len(seen) >= min_seen:
            out.append(p)
    return out


def build_bridge_matrix(model: core.DynamicProvenanceGrowth, faces: List[Face], parents: List[int], mode: str) -> Tuple[np.ndarray, List[Tuple[int,int,int]]]:
    # B maps motif face-vector coefficients (3 per face) -> sibling child coordinates.
    rows: List[Tuple[int,int,int]] = []
    for p in parents:
        for c in model.child_ids_ordered(p):
            rows.append((p, c, int(model.nodes[c].birth_order)))
    row_index = {(p,c): i for i,(p,c,_) in enumerate(rows)}
    B = np.zeros((len(rows), 3 * len(faces)), dtype=float)
    for fi, f in enumerate(faces):
        for v in f:
            p = model.nodes[int(v)].parent
            if p is None or (int(p), int(v)) not in row_index:
                continue
            n = model.nodes[int(v)]
            if mode == 'radial':
                d = unit(n.radial)
                weight = 1.0
            elif mode == 'birth_q':
                d = birth_q_vector(model, int(v))
                weight = 1.0
            elif mode == 'birth_h':
                d = birth_h_vector(model, int(v))
                weight = 1.0
            elif mode == 'conductance_birth_q':
                d = birth_q_vector(model, int(v))
                weight = float(n.g)
            elif mode == 'record_birth_q':
                d = birth_q_vector(model, int(v))
                weight = float(n.birth_g)
            else:
                raise ValueError(mode)
            B[row_index[(int(p), int(v))], 3*fi:3*fi+3] += weight * d / 3.0
    return B, rows


def block_response_operator(model: core.DynamicProvenanceGrowth, parents: List[int], op_mode: str) -> Tuple[np.ndarray, dict]:
    n = 3 * len(parents)
    R = np.zeros((n,n), dtype=float)
    circs = []
    imags = []
    norms = []
    fan_json = []
    for bi, p in enumerate(parents):
        fan = fan_weights(model, p)
        if fan is None:
            continue
        Op = response_op_from_fan(fan, op_mode)
        sl = slice(3*bi, 3*bi+3)
        R[sl, sl] = Op
        circs.append(float(fan['circulation_log_ratio']))
        imags.append(float(fan['complex_eig_imag_max']))
        norms.append(float(fan['raw_weight_norm']))
        fan_json.append({k: fan[k] for k in ('parent','children','circulation_log_ratio','raw_weight_norm','complex_eig_imag_max')})
    stats = {
        'fan_count': len(parents),
        'mean_abs_response_circulation': float(np.mean(np.abs(circs))) if circs else 0.0,
        'max_abs_response_circulation': float(np.max(np.abs(circs))) if circs else 0.0,
        'mean_response_complex_eig_imag': float(np.mean(imags)) if imags else 0.0,
        'max_response_complex_eig_imag': float(np.max(imags)) if imags else 0.0,
        'mean_raw_response_weight_norm': float(np.mean(norms)) if norms else 0.0,
        'fans_json': json.dumps(fan_json),
    }
    return R, stats


def lift_response_to_motif(B: np.ndarray, R: np.ndarray, args: argparse.Namespace) -> Tuple[np.ndarray, dict]:
    if B.size == 0 or R.size == 0 or norm(B) < EPS:
        return np.zeros((B.shape[1] if B.ndim == 2 else 0, B.shape[1] if B.ndim == 2 else 0)), {
            'bridge_rank': 0, 'bridge_min_sv': 0.0, 'bridge_cond': 0.0, 'bridge_norm': 0.0, 'response_transfer_norm': 0.0,
        }
    s = np.linalg.svd(B, compute_uv=False)
    rank = int(np.sum(s > args.bridge_rank_tol * max(float(s[0]) if len(s) else 0.0, 1.0)))
    cond = float(s[0] / (s[-1] + EPS)) if len(s) else 0.0
    Binv = np.linalg.pinv(B, rcond=args.bridge_pinv_rcond)
    T = Binv @ R @ B
    return T, {
        'bridge_rank': rank,
        'bridge_min_sv': float(s[-1]) if len(s) else 0.0,
        'bridge_max_sv': float(s[0]) if len(s) else 0.0,
        'bridge_cond': cond,
        'bridge_norm': norm(B),
        'response_transfer_norm': norm(T),
    }


def projected_spectrum_metrics(T: np.ndarray, cols: List[np.ndarray], prefix: str) -> dict:
    U = orth(cols)
    if U.shape[1] == 0 or T.size == 0:
        return {
            f'{prefix}_rank': 0, f'{prefix}_norm': 0.0, f'{prefix}_eig_imag_max': 0.0,
            f'{prefix}_eig_real_max_abs': 0.0, f'{prefix}_skew_ratio': 0.0, f'{prefix}_sym_ratio': 0.0,
            f'{prefix}_identity_resid': 1.0, f'{prefix}_minus_identity_resid': 1.0,
            f'{prefix}_square_minus_identity_resid': 1.0, f'{prefix}_cube_identity_resid': 1.0,
            f'{prefix}_eigvals_json': '[]', f'{prefix}_singular_values_json': '[]',
        }
    A = U.T @ T @ U
    eig = np.linalg.eigvals(A)
    s = np.linalg.svd(A, compute_uv=False)
    return {
        f'{prefix}_rank': int(U.shape[1]),
        f'{prefix}_norm': norm(A),
        f'{prefix}_eig_imag_max': float(max(abs(np.imag(z)) for z in eig)) if eig.size else 0.0,
        f'{prefix}_eig_real_max_abs': float(max(abs(np.real(z)) for z in eig)) if eig.size else 0.0,
        f'{prefix}_skew_ratio': norm(0.5*(A-A.T)) / (norm(A)+EPS),
        f'{prefix}_sym_ratio': norm(0.5*(A+A.T)) / (norm(A)+EPS),
        f'{prefix}_identity_resid': residual_to_sign(A, +1),
        f'{prefix}_minus_identity_resid': residual_to_sign(A, -1),
        f'{prefix}_square_minus_identity_resid': residual_to_sign(A @ A, -1),
        f'{prefix}_cube_identity_resid': residual_to_sign(A @ A @ A, +1),
        f'{prefix}_eigvals_json': json.dumps([[float(np.real(z)), float(np.imag(z))] for z in eig]),
        f'{prefix}_singular_values_json': json.dumps([float(x) for x in s]),
    }


def transfer_metrics_for_motif(model, K, assembly_row: dict, args: argparse.Namespace) -> List[dict]:
    parsed = p74.parse_assembly_pairs(model, assembly_row, args)
    if parsed is None:
        return []
    pA, pB = parsed
    faces = p72.union_faces_from_pairs(pA, pB)
    vecs, qcols, pcols, allcols, U = p75.projected_data(faces, pA, pB)
    parents = motif_parent_fans(model, faces, args)
    base = {
        'event_t': assembly_row.get('event_t',''), 'scan_id': assembly_row.get('scan_id',''),
        'context': assembly_row.get('context',''), 'face_overlap': assembly_row.get('face_overlap',''),
        'edge_overlap': assembly_row.get('edge_overlap',''), 'vertex_overlap': assembly_row.get('vertex_overlap',''),
        'A_face_a': str(list(pA['fa'])), 'A_face_b': str(list(pA['fb'])),
        'B_face_a': str(list(pB['fa'])), 'B_face_b': str(list(pB['fb'])),
        'union_faces': str([list(f) for f in faces]),
        'motif_basis_rank': int(U.shape[1]), 'Q_rank': int(orth(qcols).shape[1]), 'P_rank': int(orth(pcols).shape[1]),
        'candidate_parent_fans': str(parents), 'parent_fan_count': len(parents),
    }
    if not parents or U.shape[1] == 0:
        return [{**base, 'transfer_status': 'no_parent_fan_or_basis'}]
    rows = []
    for bridge_mode in args.bridge_modes:
        B, brow = build_bridge_matrix(model, faces, parents, bridge_mode)
        for op_mode in args.response_ops:
            R, rstats = block_response_operator(model, parents, op_mode)
            T, bstats = lift_response_to_motif(B, R, args)
            q_to_p = subspace_image_residual(T, qcols, pcols)
            p_to_q = subspace_image_residual(T, pcols, qcols)
            mean_lock = 0.5 * (q_to_p + p_to_q)
            leak = p71.span_leakage(T, allcols) if T.size else 1.0
            j2 = projected_J2_residual(T, allcols)
            qp_cols = qcols + pcols
            spec = projected_spectrum_metrics(T, qp_cols, 'projected_QP_transfer')
            transfer_active = bool(bstats['bridge_rank'] > 0 and bstats['response_transfer_norm'] > args.min_transfer_norm and rstats['max_abs_response_circulation'] > args.min_response_circulation)
            qp_lock_pass = bool(transfer_active and mean_lock <= args.qp_lock_threshold and max(q_to_p,p_to_q) <= args.qp_lock_max_threshold)
            complex_like = bool(transfer_active and spec['projected_QP_transfer_eig_imag_max'] >= args.complex_imag_threshold)
            c3_like = bool(transfer_active and spec['projected_QP_transfer_cube_identity_resid'] <= args.c3_identity_threshold and spec['projected_QP_transfer_identity_resid'] >= args.not_identity_threshold)
            double_cover_like = bool(transfer_active and spec['projected_QP_transfer_minus_identity_resid'] <= args.minus_identity_threshold)
            rows.append({
                **base, 'transfer_status': 'ok', 'bridge_mode': bridge_mode, 'response_op': op_mode,
                'bridge_rows': str(brow),
                **rstats, **bstats,
                'Q_to_P_transfer_resid': q_to_p,
                'P_to_Q_transfer_resid': p_to_q,
                'QP_transfer_mean_resid': mean_lock,
                'QP_transfer_max_resid': max(q_to_p, p_to_q),
                'motif_span_leakage': leak,
                'projected_J2_plus_I_resid': j2,
                'transfer_active': transfer_active,
                'QP_lock_transfer_pass': qp_lock_pass,
                'complex_like_transfer_pass': complex_like,
                'C3_like_transfer_pass': c3_like,
                'double_cover_like_transfer_pass': double_cover_like,
                'result_open_transfer_gate_pass': bool(qp_lock_pass or complex_like or c3_like or double_cover_like),
                **spec,
            })
    return rows


def summarize_transfer(rows: List[dict]) -> dict:
    ok = [r for r in rows if r.get('transfer_status') == 'ok']
    def count(k): return sum(1 for r in ok if fbool(r.get(k)))
    def mn(k, default=0.0):
        vals = [float(r.get(k, math.nan)) for r in ok if str(r.get(k,'')) not in ('','nan') and np.isfinite(float(r.get(k, math.nan)))]
        return float(min(vals)) if vals else default
    def mx(k, default=0.0):
        vals = [float(r.get(k, math.nan)) for r in ok if str(r.get(k,'')) not in ('','nan') and np.isfinite(float(r.get(k, math.nan)))]
        return float(max(vals)) if vals else default
    best_lock = min(ok, key=lambda r: float(r.get('QP_transfer_mean_resid', 99.0)), default=None)
    best_imag = max(ok, key=lambda r: float(r.get('projected_QP_transfer_eig_imag_max', 0.0)), default=None)
    best_c3 = min(ok, key=lambda r: float(r.get('projected_QP_transfer_cube_identity_resid', 99.0)), default=None)
    return {
        'transfer_row_count': len(rows),
        'ok_transfer_row_count': len(ok),
        'transfer_active_count': count('transfer_active'),
        'QP_lock_transfer_pass_count': count('QP_lock_transfer_pass'),
        'complex_like_transfer_pass_count': count('complex_like_transfer_pass'),
        'C3_like_transfer_pass_count': count('C3_like_transfer_pass'),
        'double_cover_like_transfer_pass_count': count('double_cover_like_transfer_pass'),
        'result_open_transfer_gate_pass_count': count('result_open_transfer_gate_pass'),
        'best_QP_transfer_mean_resid': mn('QP_transfer_mean_resid'),
        'best_QP_transfer_max_resid': mn('QP_transfer_max_resid'),
        'best_projected_J2_plus_I_resid': mn('projected_J2_plus_I_resid'),
        'max_projected_QP_transfer_eig_imag': mx('projected_QP_transfer_eig_imag_max'),
        'best_C3_cube_identity_resid': mn('projected_QP_transfer_cube_identity_resid'),
        'best_minus_identity_resid': mn('projected_QP_transfer_minus_identity_resid'),
        'max_response_circulation': mx('max_abs_response_circulation'),
        'max_response_complex_eig_imag': mx('max_response_complex_eig_imag'),
        'best_lock_row': slim_row(best_lock),
        'best_complex_row': slim_row(best_imag),
        'best_C3_row': slim_row(best_c3),
    }


def slim_row(r: Optional[dict]) -> dict:
    if not r:
        return {}
    keys = ['bridge_mode','response_op','context','QP_transfer_mean_resid','QP_transfer_max_resid','projected_QP_transfer_eig_imag_max','projected_QP_transfer_cube_identity_resid','projected_QP_transfer_minus_identity_resid','projected_J2_plus_I_resid','bridge_rank','parent_fan_count','max_abs_response_circulation','max_response_complex_eig_imag','transfer_active','QP_lock_transfer_pass','complex_like_transfer_pass','C3_like_transfer_pass','double_cover_like_transfer_pass']
    return {k: r.get(k) for k in keys}


def run_variant(variant: str, args: argparse.Namespace, out: Path) -> dict:
    model = nl.build_model(variant, args)
    model.grow(args.max_level)
    baseline_K = core.build_dynamic_outward_ngf_complex(model)
    baseline_metrics = core.full_metrics(model, baseline_K, args.source)
    vout = out / variant
    vout.mkdir(parents=True, exist_ok=True)
    K, birth_log, pairing_log, assembly_log, candidate_rows = p69.build_ablation_complex(model, args, variant, vout)
    auto_metrics = core.full_metrics(model, K, args.source)
    dm, pair_rows, top_rows, three_rows = p56.directed_metrics(model, K, pairing_log, args)
    sm, signed_rows, signed_face_rows = p59.signed_quadrature_rows(model, K, pairing_log, args)
    am, align_pair_rows, align_candidate_rows, align_candidate_summary = p61.alignment_search_metrics(model, K, pairing_log, args)
    transfer_rows: List[dict] = []
    for a in assembly_log:
        if fbool(a.get('assembly_applied')):
            transfer_rows.extend(transfer_metrics_for_motif(model, K, a, args))
    transfer_summary = summarize_transfer(transfer_rows)
    write_csv(vout / 'birth_geometry_log.csv', birth_log)
    write_csv(vout / 'assembly_pairing_log.csv', pairing_log)
    write_csv(vout / 'assembly_log.csv', assembly_log)
    write_csv(vout / 'candidate_eval_rows.csv', candidate_rows)
    write_csv(vout / 'response_to_QP_transfer_rows.csv', transfer_rows)
    write_csv(vout / 'directed_pair_rows.csv', pair_rows)
    write_csv(vout / 'signed_quadrature_rows.csv', signed_rows)
    write_csv(vout / 'alignment_pair_rows.csv', align_pair_rows)
    summary = {
        'variant': variant,
        'baseline_metrics': baseline_metrics,
        'auto_metrics': auto_metrics,
        'directed_metrics': dm,
        'signed_quadrature_metrics': sm,
        'alignment_metrics': am,
        'automatic_pairings_applied': sum(1 for x in pairing_log if fbool(x.get('applied'))),
        'assemblies_applied': sum(1 for x in assembly_log if fbool(x.get('assembly_applied'))),
        'decision_used_delta_beta_any': any(fbool(x.get('decision_used_delta_beta')) for x in pairing_log),
        'transfer_summary': transfer_summary,
    }
    (vout / 'variant_transfer_summary.json').write_text(json.dumps(summary, indent=2), encoding='utf-8')
    return summary


def make_docs(summaries: List[dict]) -> Tuple[str,str,str,str]:
    lines = ['| variant | beta | assemblies | pair harm | Q harm | P harm | base J-lock | response transfer active | QP-lock pass | complex-like pass | C3-like pass | double-cover pass | best QP transfer | max imag | best C3 resid | used Δβ? |', '|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|']
    for r in summaries:
        a = r['auto_metrics']; dm = r['directed_metrics']; am = r['alignment_metrics']; ts = r['transfer_summary']
        lines.append(f"| {r['variant']} | ({a['beta0']},{a['beta1']},{a['beta2']},{a['beta3']}) | {r['assemblies_applied']} | {dm['pair_transport_harmonic_ratio']:.6g} | {am['Q_even_harmonic_ratio']:.6g} | {am['P_odd_harmonic_ratio']:.6g} | {am['best_per_pair_mean_J_lock_resid']:.6g} | {ts['transfer_active_count']} | {ts['QP_lock_transfer_pass_count']} | {ts['complex_like_transfer_pass_count']} | {ts['C3_like_transfer_pass_count']} | {ts['double_cover_like_transfer_pass_count']} | {ts['best_QP_transfer_mean_resid']:.6g} | {ts['max_projected_QP_transfer_eig_imag']:.6g} | {ts['best_C3_cube_identity_resid']:.6g} | {r['decision_used_delta_beta_any']} |")
    table = '\n'.join(lines)
    smd = f"""# SUMMARY — response monodromy to Q/P transfer gate

This package is deliberately result-open.  It does not ask only for a Spin-1/2/double-cover signature.  It asks whether the lower sibling-response monodromy/circulation layer can be transferred to the later Q/P assembly motif layer in any clear form:

- Q/P-lock transfer,
- complex-eigen/rotational transfer,
- C3-like monodromy transfer,
- double-cover-like sign transfer,
- identity/contraction/leakage diagnostics.

{table}

The test uses actual response data from the directed sibling fans (`directed_edges`) and bridges it to the Q/P motif carrier through face/node incidence and birth-frame projections.  It does not use matrix powers of a prebuilt Q/P operator as a substitute for monodromy.
"""
    rmd = f"""# RESULTS — response monodromy to Q/P transfer gate

{table}

## Method

For each complete A/B assembly motif, the script:

1. extracts completed sibling fans touched by the motif faces;
2. builds response-layer operators from directed sibling conductance/backreaction data:
   `markov`, `skew_markov`, `forward_cycle`, `reverse_cycle`, `cycle_difference`, `cycle_sum`, and `signed_C3_skew`;
3. builds several derived bridge maps from motif face-vector data to sibling fan coordinates:
   `radial`, `birth_q`, `birth_h`, `conductance_birth_q`, `record_birth_q`;
4. lifts the response operator back to the motif carrier by least-squares reconstruction;
5. audits the resulting operator on the Q/P motif basis without imposing a target interpretation.

## Interpretation rule

This is not a Spin-only test.  A positive output could be:

- a good Q/P transfer lock,
- a robust complex-eigen rotation,
- a C3-like cyclic closure,
- a double-cover-like sign closure,
- or another nontrivial transfer mode visible in the CSV rows.

A negative output means: the lower response monodromy exists, but the tested natural incidence/birth-frame bridges do not transfer it cleanly to the later Q/P motif carrier.

## Compact JSON

```json
{json.dumps(summaries, indent=2)[:26000]}
```
"""
    audit = """# SOURCE AUDIT

This package corrects the previous alpha/double-cover shortcut.

- It does not infer spinor structure from `T^2 ≈ -alpha I` on one operator.
- It uses response-layer data from directed sibling fans, following the old dynamic birth monodromy idea: forward/reverse cyclic response asymmetry is lower-layer provenance information.
- It keeps the later Q/P assembly layer separate and asks whether a bridge transfers that response information upward.
- It does not set `i`, global `J`, Hodge, positivity, C*-norm, or a physical adjoint.
- Delta-beta/H2 are measured after growth and are not used as selection inputs.

Limitations:

- The bridge maps are candidate derived bridges from incidence plus birth-frame geometry; they are not yet theorems.
- A failure here does not prove no transfer exists; it says these local bridge candidates do not produce a clear Q/P transfer mode.
"""
    readme = """# Response monodromy to Q/P transfer gate

Run:

```bash
python3 test_response_monodromy_to_QP_transfer_gate.py
```

The output package contains per-variant summaries and full `response_to_QP_transfer_rows.csv` files.
"""
    return smd, rmd, audit, readme


def package(out: Path, zip_path: Path) -> None:
    files = [
        Path(__file__).name,
        'test_signed_Jlock_role_coupling_gate.py',
        'test_dual_assembly_order_context_ablation_gate.py',
        'test_dual_pairing_assembly_growth_rule_gate.py',
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
    ap.add_argument('--response-mode', choices=['linear','log','saturating','power_saturating','threshold_power'], default='power_saturating')
    ap.add_argument('--source', default='live', choices=['record','live','full'])
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
    ap.add_argument('--phase-sign', type=int, default=1)
    ap.add_argument('--variants', nargs='*', default=['real_growth','strict_symmetrized_control','no_backreaction'])
    ap.add_argument('--bridge-modes', nargs='*', default=['radial','birth_q','birth_h','conductance_birth_q','record_birth_q'])
    ap.add_argument('--response-ops', nargs='*', default=['markov','skew_markov','forward_cycle','reverse_cycle','cycle_difference','cycle_sum','signed_C3_skew'])
    ap.add_argument('--min-children-in-motif-for-fan', type=int, default=2)
    ap.add_argument('--bridge-pinv-rcond', type=float, default=1e-10)
    ap.add_argument('--bridge-rank-tol', type=float, default=1e-9)
    ap.add_argument('--min-transfer-norm', type=float, default=1e-8)
    ap.add_argument('--min-response-circulation', type=float, default=1e-6)
    ap.add_argument('--qp-lock-threshold', type=float, default=0.25)
    ap.add_argument('--qp-lock-max-threshold', type=float, default=0.35)
    ap.add_argument('--complex-imag-threshold', type=float, default=0.05)
    ap.add_argument('--c3-identity-threshold', type=float, default=0.35)
    ap.add_argument('--not-identity-threshold', type=float, default=0.35)
    ap.add_argument('--minus-identity-threshold', type=float, default=0.35)
    ap.add_argument('--out', default='response_monodromy_to_QP_transfer_out_L2')
    ap.add_argument('--zip', default='cnna_response_monodromy_to_QP_transfer_gate_pkg_L2.zip')
    args = ap.parse_args()

    p70.patch_modules()
    out = Path(args.out)
    if out.exists(): shutil.rmtree(out)
    out.mkdir(parents=True)
    summaries = []
    for v in args.variants:
        summaries.append(run_variant(v, args, out))
    write_csv(out / 'comparative_response_monodromy_to_QP_transfer_summary.csv', [
        {
            'variant': s['variant'],
            'beta0': s['auto_metrics']['beta0'], 'beta1': s['auto_metrics']['beta1'], 'beta2': s['auto_metrics']['beta2'], 'beta3': s['auto_metrics']['beta3'],
            'assemblies_applied': s['assemblies_applied'],
            'pair_harmonic': s['directed_metrics']['pair_transport_harmonic_ratio'],
            'Q_harmonic': s['alignment_metrics']['Q_even_harmonic_ratio'],
            'P_harmonic': s['alignment_metrics']['P_odd_harmonic_ratio'],
            **{k: v for k,v in s['transfer_summary'].items() if not isinstance(v, dict)},
            'decision_used_delta_beta_any': s['decision_used_delta_beta_any'],
        } for s in summaries
    ])
    smd, rmd, audit, readme = make_docs(summaries)
    (out / 'SUMMARY.md').write_text(smd, encoding='utf-8')
    (out / 'RESULTS.md').write_text(rmd, encoding='utf-8')
    (out / 'SOURCE_AUDIT.md').write_text(audit, encoding='utf-8')
    (out / 'README.md').write_text(readme, encoding='utf-8')
    (out / 'summary.json').write_text(json.dumps({'summaries': summaries}, indent=2), encoding='utf-8')
    package(out, Path(args.zip))
    print(json.dumps({'out': str(out), 'zip': str(Path(args.zip).resolve()), 'summaries': summaries}, indent=2)[:6000])


if __name__ == '__main__':
    main()
