#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import math
import shutil
import zipfile
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np

import cnna_non_shelling_core as core
import test_nonlinear_asymmetry_cascade_growth as nl
import test_harmonic_k_orientation_kappa_gate as hk

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


def unit(v: np.ndarray) -> np.ndarray:
    n = float(np.linalg.norm(v))
    if n < EPS:
        return np.zeros_like(v, dtype=float)
    return np.array(v, dtype=float) / n


def skew(M: np.ndarray) -> np.ndarray:
    return 0.5 * (M - M.T)


def axial(A: np.ndarray) -> np.ndarray:
    # A x = w × x convention for skew matrix [[0,-wz,wy],[wz,0,-wx],[-wy,wx,0]]
    return np.array([A[2, 1], A[0, 2], A[1, 0]], dtype=float)


def vertex_operator_directed(
    model: core.DynamicProvenanceGrowth,
    node: int,
    source: str,
    phase_sign: int,
    antisym_eta: float,
    erase_phase_for_strict_sym: bool,
) -> np.ndarray:
    """Directed real birth-transport operator.

    Main anti-smuggling rule: no final sym(M).  The antisymmetric part is not free;
    it is derived from the already present ternary birth order via q and h.
    """
    n = model.nodes[node]
    r, e1, e2 = n.radial, n.e1, n.e2
    if n.parent is None or n.birth_order == 0:
        order_phase = 0.0
    else:
        if erase_phase_for_strict_sym and getattr(model, 'growth_rule', '') == 'symmetrized_birth':
            effective_order = 2
        else:
            effective_order = n.birth_order
        order_phase = float(phase_sign) * 2.0 * math.pi * (effective_order - 1) / 3.0
    q = math.cos(order_phase) * e1 + math.sin(order_phase) * e2
    h = core.unit(0.7 * r + 0.3 * q)
    birth = n.birth_g
    live = n.g
    aging = max(0.0, live - birth)
    if source == 'record':
        a, b, c = birth, 0.22 * birth, 0.08 * birth
    elif source == 'live':
        a, b, c = live, 0.25 * birth + 0.55 * aging, 0.12 * live
    elif source == 'full':
        a, b, c = 0.5 * (birth + live), 0.235 * birth + 0.275 * aging, 0.1 * live
    elif source == 'handoff':
        inc = 0.0
        if n.parent is not None:
            inc = model.directed_edges.get((n.parent, node), 0.0) + model.directed_edges.get((node, n.parent), 0.0)
        a, b, c = birth + inc, 0.18 * live + inc, 0.15 * inc + 0.05 * birth
    elif source == 'aging':
        a, b, c = aging + 0.1 * birth, 0.6 * aging + 0.03 * birth, 0.3 * aging
    else:
        raise ValueError(source)
    metric_part = a * np.outer(r, r) + b * np.outer(q, q) + c * np.outer(h, h) + 0.04 * birth * np.eye(3)
    transport_part = antisym_eta * (0.5 * b + c + 0.07 * birth) * (np.outer(q, h) - np.outer(h, q))
    return metric_part + transport_part


def face_K_directed(
    model: core.DynamicProvenanceGrowth,
    face: Face,
    source: str,
    phase_sign: int,
    antisym_eta: float,
    erase_phase_for_strict_sym: bool,
) -> np.ndarray:
    a, b, c = face
    Sa = vertex_operator_directed(model, a, source, phase_sign, antisym_eta, erase_phase_for_strict_sym)
    Sb = vertex_operator_directed(model, b, source, phase_sign, antisym_eta, erase_phase_for_strict_sym)
    Sc = vertex_operator_directed(model, c, source, phase_sign, antisym_eta, erase_phase_for_strict_sym)
    Aab = Sb - Sa
    Abc = Sc - Sb
    return skew(Aab @ Abc - Abc @ Aab)


def rotation_from_to(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    a = np.array(a, dtype=float)
    b = np.array(b, dtype=float)
    na = float(np.linalg.norm(a))
    nb = float(np.linalg.norm(b))
    if na < EPS or nb < EPS:
        return np.eye(3)
    a = a / na
    b = b / nb
    c = float(np.clip(np.dot(a, b), -1.0, 1.0))
    if c > 1.0 - 1e-10:
        return np.eye(3)
    if c < -1.0 + 1e-10:
        axis = np.cross(a, np.array([1.0, 0.0, 0.0]))
        if np.linalg.norm(axis) < 1e-8:
            axis = np.cross(a, np.array([0.0, 1.0, 0.0]))
        axis = axis / (np.linalg.norm(axis) + EPS)
        K = np.array([[0, -axis[2], axis[1]], [axis[2], 0, -axis[0]], [-axis[1], axis[0], 0]], dtype=float)
        return np.eye(3) + 2.0 * (K @ K)
    v = np.cross(a, b)
    s = float(np.linalg.norm(v))
    K = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]], dtype=float)
    return np.eye(3) + K + K @ K * ((1.0 - c) / (s * s + EPS))


def parse_pair_faces(log: dict) -> Tuple[Face | None, Face | None]:
    try:
        fa = core.parse_face_string(str(log.get('face_a', '')))
        fb_txt = str(log.get('face_b', '')).split('perm=')[0].strip()
        fb = core.parse_face_string(fb_txt)
        if len(fa) == 3 and len(fb) == 3:
            return tuple(sorted(fa)), tuple(sorted(fb))
    except Exception:
        pass
    return None, None


def vector_field_projection(H: np.ndarray, W: np.ndarray) -> np.ndarray:
    if H.size == 0 or H.shape[1] == 0:
        return np.zeros_like(W)
    return H @ (H.T @ W)


def scalar_projection(H: np.ndarray, w: np.ndarray) -> np.ndarray:
    if H.size == 0 or H.shape[1] == 0:
        return np.zeros_like(w)
    return H @ (H.T @ w)


def support_coherence(W: np.ndarray, threshold_scale: float = 1e-8) -> Tuple[float, float, int]:
    if W.size == 0:
        return 0.0, 0.0, 0
    norms = np.linalg.norm(W, axis=1)
    if len(norms) == 0 or float(np.max(norms)) < EPS:
        return 0.0, 0.0, 0
    mask = norms > max(1e-10, threshold_scale * float(np.max(norms)))
    if not np.any(mask):
        return 0.0, 0.0, 0
    U = W[mask] / (norms[mask][:, None] + EPS)
    mean_vec = np.mean(U, axis=0)
    coherence = float(np.linalg.norm(mean_vec))
    # signless axis coherence detects undirected alignment without importing an orientation.
    Q = U.T @ U / max(1, U.shape[0])
    evals = np.linalg.eigvalsh(Q)
    axis_coherence = float(max(evals[-1], 0.0)) if len(evals) else 0.0
    return coherence, axis_coherence, int(np.sum(mask))


def incident_faces_by_edge(faces: List[Face]) -> Dict[Edge, List[int]]:
    out: Dict[Edge, List[int]] = {}
    for i, f in enumerate(faces):
        a, b, c = f
        for e in (tuple(sorted((a, b))), tuple(sorted((a, c))), tuple(sorted((b, c)))):
            out.setdefault(e, []).append(i)
    return out


def three_face_coherence(faces: List[Face], W: np.ndarray) -> Tuple[float, float, int, List[dict]]:
    norms = np.linalg.norm(W, axis=1) if len(faces) else np.array([])
    edge_map = incident_faces_by_edge(faces)
    vals: List[float] = []
    defects: List[float] = []
    rows: List[dict] = []
    for e, inds in edge_map.items():
        active = [i for i in inds if i < len(norms) and norms[i] > max(1e-10, 1e-8 * (float(np.max(norms)) if len(norms) else 1.0))]
        if len(active) < 3:
            continue
        # use all deterministic triples, but at L2 this remains tiny.
        for ia in range(len(active)):
            for ib in range(ia + 1, len(active)):
                for ic in range(ib + 1, len(active)):
                    tri = [active[ia], active[ib], active[ic]]
                    U = np.array([W[j] / (norms[j] + EPS) for j in tri], dtype=float)
                    coh = float(np.linalg.norm(np.mean(U, axis=0)))
                    pair_cos = [float(np.dot(U[p], U[q])) for p, q in ((0, 1), (0, 2), (1, 2))]
                    defect = float(1.0 - coh)
                    vals.append(coh)
                    defects.append(defect)
                    rows.append({
                        'shared_edge': str(list(e)),
                        'faces': str([list(faces[j]) for j in tri]),
                        'three_face_coherence': coh,
                        'three_face_phase_defect': defect,
                        'pair_cosines': str(pair_cos),
                        'mean_pair_cosine': float(np.mean(pair_cos)),
                    })
    rows.sort(key=lambda r: r['three_face_coherence'], reverse=True)
    if not vals:
        return 0.0, 0.0, 0, rows
    return float(np.mean(vals)), float(np.mean(defects)), len(vals), rows


def transported_pair_fields_directed(
    model: core.DynamicProvenanceGrowth,
    K: core.SimplicialComplex,
    pairing_log: List[dict],
    source: str,
    phase_sign: int,
    antisym_eta: float,
    erase_phase_for_strict_sym: bool,
) -> Tuple[np.ndarray, np.ndarray, List[dict]]:
    faces = K.faces()
    idx = {tuple(f): i for i, f in enumerate(faces)}
    W_pair = np.zeros((len(faces), 3), dtype=float)
    scalar_pair = np.zeros(len(faces), dtype=float)
    rows: List[dict] = []
    for k, log in enumerate(pairing_log):
        if not log.get('applied'):
            continue
        fa, fb = parse_pair_faces(log)
        if fa is None or fb is None or fa not in idx or fb not in idx:
            continue
        ka = axial(face_K_directed(model, fa, source, phase_sign, antisym_eta, erase_phase_for_strict_sym))
        kb = axial(face_K_directed(model, fb, source, phase_sign, antisym_eta, erase_phase_for_strict_sym))
        na = hk.face_normal(model, fa, 'outward')
        nb = hk.face_normal(model, fb, 'outward')
        R_b_to_a = rotation_from_to(nb, -na)
        kb_to_a_reversed = R_b_to_a @ kb
        pair_vec_a = ka + kb_to_a_reversed
        pair_vec_b = -(R_b_to_a.T @ pair_vec_a)
        ia, ib = idx[fa], idx[fb]
        W_pair[ia] += pair_vec_a
        W_pair[ib] += pair_vec_b
        scalar_strength = float(np.linalg.norm(pair_vec_a))
        scalar_pair[ia] += scalar_strength
        scalar_pair[ib] += scalar_strength
        rows.append({
            'pair_index': k,
            'event_t': log.get('event_t', ''),
            'cascade_index': log.get('cascade_index', ''),
            'move_class': log.get('move_class', ''),
            'face_a': str(list(fa)),
            'face_b': str(list(fb)),
            'A_rank_score': log.get('A_rank_score', ''),
            'A_invariant': log.get('A_invariant', ''),
            'directed_imbalance': log.get('directed_imbalance', ''),
            'transverse_complementarity': log.get('transverse_complementarity', ''),
            'measured_delta_beta1': log.get('measured_delta_beta1', ''),
            'measured_delta_beta2': log.get('measured_delta_beta2', ''),
            'decision_used_delta_beta': log.get('decision_used_delta_beta', ''),
            'K_axial_a_norm': float(np.linalg.norm(ka)),
            'K_axial_b_norm': float(np.linalg.norm(kb)),
            'transported_reversed_b_norm': float(np.linalg.norm(kb_to_a_reversed)),
            'pair_vec_a_norm': float(np.linalg.norm(pair_vec_a)),
            'pair_vec_b_norm': float(np.linalg.norm(pair_vec_b)),
            'a_b_transport_cosine': float(np.dot(ka, kb_to_a_reversed) / ((np.linalg.norm(ka) * np.linalg.norm(kb_to_a_reversed)) + EPS)),
        })
    return W_pair, scalar_pair, rows


def kappa_ratios(model: core.DynamicProvenanceGrowth, faces: List[Face], W: np.ndarray) -> dict:
    if not faces:
        return {
            'normal_flux_signed_ratio': 0.0, 'normal_flux_abs_ratio': 0.0, 'kappa_orientation_ratio': 0.0,
            'birth_normal_flux_signed_ratio': 0.0, 'birth_normal_flux_abs_ratio': 0.0, 'kappa_birth_orientation_ratio': 0.0,
        }
    norms = np.linalg.norm(W, axis=1)
    normals_out = np.array([hk.face_normal(model, f, 'outward') for f in faces], dtype=float)
    normals_birth = np.array([hk.face_normal(model, f, 'birth_order') for f in faces], dtype=float)
    areas = np.array([max(hk.face_area(model, f), EPS) for f in faces], dtype=float)
    denom = float(np.sum(norms * areas)) + EPS
    dot_out = np.einsum('ij,ij->i', W, normals_out)
    dot_birth = np.einsum('ij,ij->i', W, normals_birth)
    signed_out = float(np.sum(dot_out * areas)) / denom
    abs_out = float(np.sum(np.abs(dot_out) * areas)) / denom
    signed_birth = float(np.sum(dot_birth * areas)) / denom
    abs_birth = float(np.sum(np.abs(dot_birth) * areas)) / denom
    return {
        'normal_flux_signed_ratio': signed_out,
        'normal_flux_abs_ratio': abs_out,
        'kappa_orientation_ratio': abs(signed_out) / (abs_out + EPS),
        'birth_normal_flux_signed_ratio': signed_birth,
        'birth_normal_flux_abs_ratio': abs_birth,
        'kappa_birth_orientation_ratio': abs(signed_birth) / (abs_birth + EPS),
    }


def directed_metrics(
    model: core.DynamicProvenanceGrowth,
    K: core.SimplicialComplex,
    pairing_log: List[dict],
    args: argparse.Namespace,
) -> Tuple[dict, List[dict], List[dict], List[dict]]:
    faces = K.faces()
    topo = core.topology(K)
    H, eigs = hk.harmonic_basis_faces(K)
    W_raw = np.array([
        axial(face_K_directed(model, f, args.source, args.phase_sign, args.antisym_eta, args.erase_phase_for_strict_sym))
        for f in faces
    ], dtype=float) if faces else np.zeros((0, 3), dtype=float)
    W_raw_H = vector_field_projection(H, W_raw)
    raw_total = float(np.linalg.norm(W_raw)) + EPS
    raw_h = float(np.linalg.norm(W_raw_H))
    raw_coh, raw_axis_coh, raw_support = support_coherence(W_raw)
    raw_H_coh, raw_H_axis_coh, raw_H_support = support_coherence(W_raw_H)
    raw3_coh, raw3_defect, raw3_count, raw3_rows = three_face_coherence(faces, W_raw)

    W_pair, scalar_pair, pair_rows = transported_pair_fields_directed(
        model, K, pairing_log, args.source, args.phase_sign, args.antisym_eta, args.erase_phase_for_strict_sym
    )
    W_pair_H = vector_field_projection(H, W_pair)
    scalar_H = scalar_projection(H, scalar_pair)
    pair_total = float(np.linalg.norm(W_pair)) + EPS
    pair_H_total = float(np.linalg.norm(W_pair_H))
    scalar_total = float(np.linalg.norm(scalar_pair)) + EPS
    scalar_H_total = float(np.linalg.norm(scalar_H))
    pair_raw_coh, pair_raw_axis_coh, pair_raw_support = support_coherence(W_pair)
    pair_H_coh, pair_H_axis_coh, pair_H_support = support_coherence(W_pair_H)
    pair3_coh, pair3_defect, pair3_count, pair3_rows = three_face_coherence(faces, W_pair)

    raw_kappa = kappa_ratios(model, faces, W_raw)
    raw_H_kappa = kappa_ratios(model, faces, W_raw_H)
    pair_raw_kappa = kappa_ratios(model, faces, W_pair)
    pair_H_kappa = kappa_ratios(model, faces, W_pair_H)

    top_rows: List[dict] = []
    H_norms = np.linalg.norm(W_pair_H, axis=1) if len(faces) else np.array([])
    raw_norms = np.linalg.norm(W_raw, axis=1) if len(faces) else np.array([])
    pair_norms = np.linalg.norm(W_pair, axis=1) if len(faces) else np.array([])
    for i, f in enumerate(faces):
        if float(pair_norms[i]) <= 0 and float(H_norms[i]) <= 0 and float(raw_norms[i]) <= 0:
            continue
        top_rows.append({
            'face': str(list(f)),
            'birth_orders': str([model.nodes[v].birth_order for v in f]),
            'birth_times': str([model.nodes[v].birth_time for v in f]),
            'raw_K_norm': float(raw_norms[i]),
            'pair_transport_norm': float(pair_norms[i]),
            'pair_H_norm': float(H_norms[i]),
            'pair_scalar_strength': float(scalar_pair[i]) if len(scalar_pair) else 0.0,
            'pair_scalar_H_value': float(scalar_H[i]) if len(scalar_H) else 0.0,
        })
    top_rows.sort(key=lambda r: (r['pair_H_norm'], r['pair_transport_norm'], r['raw_K_norm']), reverse=True)

    metrics = {
        'beta0': topo['beta0'], 'beta1': topo['beta1'], 'beta2': topo['beta2'], 'beta3': topo['beta3'],
        'harmonic_dim_real': int(H.shape[1]) if H.ndim == 2 else 0,
        'applied_pair_count': len(pair_rows),
        'raw_K_total_norm': raw_total - EPS,
        'raw_K_harmonic_norm': raw_h,
        'raw_K_harmonic_ratio': raw_h / raw_total,
        'raw_local_orientation_coherence': raw_coh,
        'raw_local_axis_coherence': raw_axis_coh,
        'raw_support_count': raw_support,
        'raw_H_orientation_coherence': raw_H_coh,
        'raw_H_axis_coherence': raw_H_axis_coh,
        'raw_H_support_count': raw_H_support,
        'raw_3face_coherence': raw3_coh,
        'shared_edge_3face_phase_defect': raw3_defect,
        'raw_3face_count': raw3_count,
        'pair_transport_total_norm': pair_total - EPS,
        'pair_transport_harmonic_norm': pair_H_total,
        'pair_transport_harmonic_ratio': pair_H_total / pair_total,
        'pair_scalar_total_norm': scalar_total - EPS,
        'pair_scalar_harmonic_norm': scalar_H_total,
        'pair_scalar_harmonic_ratio': scalar_H_total / scalar_total,
        'pair_raw_orientation_coherence': pair_raw_coh,
        'pair_raw_axis_coherence': pair_raw_axis_coh,
        'pair_raw_support_count': pair_raw_support,
        'pair_orientation_coherence': pair_H_coh,
        'pair_axis_coherence': pair_H_axis_coh,
        'pair_H_support_count': pair_H_support,
        'pair_3face_coherence': pair3_coh,
        'pair_shared_edge_3face_phase_defect': pair3_defect,
        'pair_3face_count': pair3_count,
        'mean_pair_transport_cosine': float(np.mean([float(r['a_b_transport_cosine']) for r in pair_rows])) if pair_rows else 0.0,
        'interfan_phase_transport_residual': float(np.mean([1.0 - abs(float(r['a_b_transport_cosine'])) for r in pair_rows])) if pair_rows else 0.0,
        'decision_used_delta_beta_any': any(str(r.get('decision_used_delta_beta', '')).lower() == 'true' for r in pair_rows),
        'measured_delta_beta2_sum': sum(int(float(r.get('measured_delta_beta2', 0) or 0)) for r in pair_rows),
    }
    metrics.update({f'raw_{k}': v for k, v in raw_kappa.items()})
    metrics.update({f'raw_H_{k}': v for k, v in raw_H_kappa.items()})
    metrics.update({f'pair_raw_{k}': v for k, v in pair_raw_kappa.items()})
    metrics.update({f'pair_{k}': v for k, v in pair_H_kappa.items()})
    return metrics, pair_rows, top_rows, raw3_rows[:args.keep_top_faces] + pair3_rows[:args.keep_top_faces]


def build_variant(variant: str, args: argparse.Namespace, out: Path):
    model = nl.build_model(variant, args)
    model.grow(args.max_level)
    baseline_K = core.build_dynamic_outward_ngf_complex(model)
    baseline_metrics = core.full_metrics(model, baseline_K, args.source)
    auto_K, birth_log, pairing_log, candidate_sample, scans = nl.build_nonlinear_auto_complex(model, args, out / variant, variant)
    auto_metrics = core.full_metrics(model, auto_K, args.source)
    return model, auto_K, baseline_metrics, auto_metrics, birth_log, pairing_log, candidate_sample, scans


def run_variant(variant: str, args: argparse.Namespace, out: Path) -> dict:
    vout = out / variant
    vout.mkdir(parents=True, exist_ok=True)
    model, K, baseline, auto, birth_log, pairing_log, candidate_sample, scans = build_variant(variant, args, out)
    metrics, pair_rows, face_rows, three_rows = directed_metrics(model, K, pairing_log, args)
    write_csv(vout / 'birth_geometry_log.csv', birth_log)
    write_csv(vout / 'nonlinear_pairing_cascade_log.csv', pairing_log)
    write_csv(vout / 'directed_antisym_pair_transport_pairs.csv', pair_rows)
    write_csv(vout / 'directed_antisym_pair_transport_faces_top.csv', face_rows[:args.keep_top_faces])
    write_csv(vout / 'three_face_interfan_coherence_top.csv', three_rows[:args.keep_top_faces])
    summary = {
        'variant': variant,
        'max_level': args.max_level,
        'source': args.source,
        'antisym_eta': args.antisym_eta,
        'phase_sign': args.phase_sign,
        'erase_phase_for_strict_sym': args.erase_phase_for_strict_sym,
        'baseline_metrics_legacy_core': baseline,
        'auto_metrics_legacy_core': auto,
        'directed_antisym_metrics': metrics,
        'automatic_pairings_applied': sum(1 for x in pairing_log if x.get('applied')),
        'automatic_pairing_attempts_logged': len(pairing_log),
        'births_with_cascade_logs': scans,
        'interpretation_flags': {
            'beta2_opened': auto['beta2'] > baseline['beta2'],
            'raw_local_orientation_present': metrics['raw_local_orientation_coherence'] > args.coherence_threshold,
            'raw_3face_coherence_present': metrics['raw_3face_coherence'] > args.coherence_threshold,
            'pair_transport_harmonic_positive': metrics['pair_transport_harmonic_ratio'] > args.harmonic_positive_threshold,
            'pair_kappa_signed_bias_positive': metrics['pair_kappa_orientation_ratio'] > args.kappa_orientation_threshold,
            'pair_birth_kappa_signed_bias_positive': metrics['pair_kappa_birth_orientation_ratio'] > args.kappa_orientation_threshold,
            'pair_orientation_coherence_positive': metrics['pair_orientation_coherence'] > args.coherence_threshold,
            'strict_sym_should_kill': variant == 'strict_symmetrized_control' and metrics['pair_transport_total_norm'] <= args.zero_threshold,
            'decision_used_delta_beta_any': metrics['decision_used_delta_beta_any'],
        },
    }
    (vout / 'variant_directed_antisym_pair_transport_summary.json').write_text(json.dumps(summary, indent=2), encoding='utf-8')
    return summary


def write_comparative(out: Path, rows: List[dict]) -> None:
    flat = []
    for r in rows:
        a = r['auto_metrics_legacy_core']; m = r['directed_antisym_metrics']
        flat.append({
            'variant': r['variant'],
            'beta0': a['beta0'], 'beta1': a['beta1'], 'beta2': a['beta2'], 'beta3': a['beta3'],
            'pairings': r['automatic_pairings_applied'],
            'raw_K_harmonic_ratio': m['raw_K_harmonic_ratio'],
            'raw_local_orientation_coherence': m['raw_local_orientation_coherence'],
            'raw_3face_coherence': m['raw_3face_coherence'],
            'shared_edge_3face_phase_defect': m['shared_edge_3face_phase_defect'],
            'pair_transport_harmonic_ratio': m['pair_transport_harmonic_ratio'],
            'pair_scalar_harmonic_ratio': m['pair_scalar_harmonic_ratio'],
            'pair_kappa_orientation_ratio': m['pair_kappa_orientation_ratio'],
            'pair_kappa_birth_orientation_ratio': m['pair_kappa_birth_orientation_ratio'],
            'pair_raw_orientation_coherence': m['pair_raw_orientation_coherence'],
            'pair_orientation_coherence': m['pair_orientation_coherence'],
            'pair_3face_coherence': m['pair_3face_coherence'],
            'interfan_phase_transport_residual': m['interfan_phase_transport_residual'],
            'mean_pair_transport_cosine': m['mean_pair_transport_cosine'],
            'decision_used_delta_beta_any': m['decision_used_delta_beta_any'],
        })
    write_csv(out / 'comparative_directed_antisym_pair_transport_summary.csv', flat)


def make_docs(summary: dict) -> Tuple[str, str, str, str]:
    rows = summary['variant_rows']
    lines = [
        '| variant | beta auto | pairings | raw K harm | raw local coh | raw 3face coh | pair harm | pair kappa | pair birth kappa | pair H coh | interfan residual | used Δβ? |',
        '|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|',
    ]
    for r in rows:
        a = r['auto_metrics_legacy_core']; m = r['directed_antisym_metrics']
        lines.append(
            f"| {r['variant']} | ({a['beta0']},{a['beta1']},{a['beta2']},{a['beta3']}) | {r['automatic_pairings_applied']} | "
            f"{m['raw_K_harmonic_ratio']:.6g} | {m['raw_local_orientation_coherence']:.6g} | {m['raw_3face_coherence']:.6g} | "
            f"{m['pair_transport_harmonic_ratio']:.6g} | {m['pair_kappa_orientation_ratio']:.6g} | {m['pair_kappa_birth_orientation_ratio']:.6g} | "
            f"{m['pair_orientation_coherence']:.6g} | {m['interfan_phase_transport_residual']:.6g} | {m['decision_used_delta_beta_any']} |"
        )
    table = '\n'.join(lines)
    readme = """# Directed antisymmetric birth transport + pairing transport coherence gate

Run:

```bash
python3 test_pairing_transport_antisym_birth_coherence_gate.py
```

This package deliberately does **not** apply a final `sym(M)` in the tested vertex operator.  The strict control is the growth/control variant, not an operator symmetrization.
"""
    smd = f"""# SUMMARY — directed antisymmetric birth transport through pairing/inter-fan gate

This package combines the missing pieces from the recent sequence:

```text
51: asymmetry-gated beta2 carrier + pairing transport + harmonic projection
55: local birth-order-derived antisymmetric transport operator without final sym(M)
```

The tested operator is real and directed:

```text
M = metric_part + eta * strength * (q⊗h - h⊗q)
```

No final `sym(M)` is applied in the tested path.  The antisymmetric term is derived from ternary birth order; it is not an input `J`, orientation, Hodge star, positivity, norm, spin, or complex scalar.

{table}

Interpretation must remain conservative: a positive pair harmonic component is not a J-derivation.  The gate asks whether local directed birth orientation survives pairing transport and inter-face/3-face coherence checks.
"""
    rmd = f"""# RESULTS — directed antisymmetric birth transport through pairing/inter-fan gate

## Comparative table

{table}

## Gate reading

- `raw_K_harmonic_ratio`: harmonic projection of the local directed face-K field.
- `raw_local_orientation_coherence`: pre-H2 coherence of local axial K vectors.
- `raw_3face_coherence`: coherence of triples of incident faces around shared edges.
- `pair_transport_harmonic_ratio`: H2 projection of the actual asymmetry-gated pairing transport field.
- `pair_kappa_orientation_ratio` / `pair_kappa_birth_orientation_ratio`: signed normal/birth-normal bias of the harmonic transported field.
- `interfan_phase_transport_residual`: 1 - |cos| averaged over transported actual pairings.
- `decision_used_delta_beta_any` must remain false.

## Conservative status

This is not a complex-structure derivation and not a real `*`-structure.  It is a falsifiable diagnostic for the missing combination: local directed birth orientation + beta2 carrier + pairing/inter-fan transport.
"""
    audit = """# SOURCE AUDIT

Hard methodological constraints:

- no i/J/Hodge/star/positivity/norm/spin/orientation as input;
- no final sym(M) in the tested directed vertex operator;
- the antisymmetric part is not free but built from birth-order phase q and h;
- beta changes are measured after moves and not used in the pairing decision.

What this package combines:

- Package 51 style nonlinear asymmetry-gated complement pairing and pairing transport;
- Package 55 style antisymmetric birth transport operator;
- additional raw local and 3-face/inter-fan coherence diagnostics before H2 projection.

Known limitation:

The geometry and pairing rule are still deterministic test scaffolding, not a final CNNA theorem.  NGF/CQNM remains comparison language only.
"""
    return smd, rmd, audit, readme


def package(out: Path, zip_path: Path) -> None:
    files = [
        Path(__file__).name,
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
    ap.add_argument('--source', default='live', choices=['record', 'live', 'full', 'handoff', 'aging'])
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
    ap.add_argument('--keep-top-faces', type=int, default=80)
    ap.add_argument('--harmonic-positive-threshold', type=float, default=1e-4)
    ap.add_argument('--coherence-threshold', type=float, default=0.15)
    ap.add_argument('--kappa-orientation-threshold', type=float, default=0.15)
    ap.add_argument('--zero-threshold', type=float, default=1e-10)
    ap.add_argument('--antisym-eta', type=float, default=1.0)
    ap.add_argument('--phase-sign', type=int, default=1, choices=[-1, 1])
    ap.add_argument('--erase-phase-for-strict-sym', action='store_true', default=True)
    ap.add_argument('--variants', nargs='*', default=['real_growth', 'strict_symmetrized_control', 'no_backreaction'])
    ap.add_argument('--out', default='pairing_transport_antisym_birth_coherence_out_L2')
    ap.add_argument('--zip', default='cnna_pairing_transport_antisym_birth_coherence_gate_pkg_L2.zip')
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
                'auto_beta': [r['auto_metrics_legacy_core'][f'beta{i}'] for i in range(4)],
                'pairings': r['automatic_pairings_applied'],
                'raw_K_harmonic_ratio': r['directed_antisym_metrics']['raw_K_harmonic_ratio'],
                'raw_local_orientation_coherence': r['directed_antisym_metrics']['raw_local_orientation_coherence'],
                'raw_3face_coherence': r['directed_antisym_metrics']['raw_3face_coherence'],
                'pair_transport_harmonic_ratio': r['directed_antisym_metrics']['pair_transport_harmonic_ratio'],
                'pair_kappa_orientation_ratio': r['directed_antisym_metrics']['pair_kappa_orientation_ratio'],
                'pair_kappa_birth_orientation_ratio': r['directed_antisym_metrics']['pair_kappa_birth_orientation_ratio'],
                'pair_orientation_coherence': r['directed_antisym_metrics']['pair_orientation_coherence'],
                'interfan_phase_transport_residual': r['directed_antisym_metrics']['interfan_phase_transport_residual'],
                'decision_used_delta_beta_any': r['directed_antisym_metrics']['decision_used_delta_beta_any'],
            } for r in rows
        ]
    }, indent=2))


if __name__ == '__main__':
    main()
