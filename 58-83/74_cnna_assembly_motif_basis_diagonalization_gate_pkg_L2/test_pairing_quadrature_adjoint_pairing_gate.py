#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import shutil
import zipfile
from pathlib import Path
from typing import List, Tuple

import numpy as np

import cnna_non_shelling_core as core
import test_nonlinear_asymmetry_cascade_growth as nl
import test_harmonic_k_orientation_kappa_gate as hk
import test_pairing_transport_antisym_birth_coherence_gate as p56
import test_pairing_quadrature_split_symplectic_defect_gate as p58

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


def frob(A: np.ndarray) -> float:
    return float(np.linalg.norm(A))


def rel(A: np.ndarray, denom: float) -> float:
    return frob(A) / (denom + EPS)


def block_J(R_b_to_a: np.ndarray) -> np.ndarray:
    R = R_b_to_a
    Z = np.zeros((3, 3), dtype=float)
    return np.block([[Z, R], [-R.T, Z]])


def block_C(R_b_to_a: np.ndarray) -> np.ndarray:
    R = R_b_to_a
    Z = np.zeros((3, 3), dtype=float)
    return np.block([[Z, R], [R.T, Z]])


def block_vec(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    return np.concatenate([a, b])


def safe_ratio(num: float, den: float) -> float:
    return float(num) / (float(den) + EPS)


def signed_commutator_stats(Q: np.ndarray, P: np.ndarray, n_birth: np.ndarray, n_out: np.ndarray) -> dict:
    cross = np.cross(Q, P)
    abs_area = float(np.linalg.norm(cross))
    return {
        'comm_abs_area': abs_area,
        'comm_signed_birth': float(np.dot(cross, n_birth)),
        'comm_signed_outward': float(np.dot(cross, n_out)),
        'comm_signed_birth_over_abs': float(np.dot(cross, n_birth)) / (abs_area + EPS),
        'comm_signed_outward_over_abs': float(np.dot(cross, n_out)) / (abs_area + EPS),
    }


def adjoint_pairing_metrics(
    model: core.DynamicProvenanceGrowth,
    K: core.SimplicialComplex,
    pairing_log: List[dict],
    args: argparse.Namespace,
) -> tuple[dict, List[dict], List[dict]]:
    faces = K.faces()
    idx = {tuple(f): i for i, f in enumerate(faces)}
    topo = core.topology(K)
    H, _ = hk.harmonic_basis_faces(K)

    W_Q = np.zeros((len(faces), 3), dtype=float)
    W_P = np.zeros((len(faces), 3), dtype=float)
    W_comm = np.zeros((len(faces), 3), dtype=float)

    pair_rows: List[dict] = []
    j_square_resids: List[float] = []
    c_square_resids: List[float] = []
    cjc_resids: List[float] = []
    orth_resids: List[float] = []
    q_even_resids: List[float] = []
    p_odd_resids: List[float] = []
    jq_to_p_resids: List[float] = []
    jp_to_minus_q_resids: List[float] = []
    comm_abs: List[float] = []
    comm_signed_birth: List[float] = []
    comm_signed_out: List[float] = []
    block_comm_grade: List[float] = []

    I6 = np.eye(6)
    G = np.diag([1.0, 1.0, 1.0, -1.0, -1.0, -1.0])
    G_norm = frob(G)

    for k, log in enumerate(pairing_log):
        if not log.get('applied'):
            continue
        fa, fb = p56.parse_pair_faces(log)
        if fa is None or fb is None or fa not in idx or fb not in idx:
            continue

        ka = p56.axial(p56.face_K_directed(
            model, fa, args.source, args.phase_sign, args.antisym_eta, args.erase_phase_for_strict_sym
        ))
        kb = p56.axial(p56.face_K_directed(
            model, fb, args.source, args.phase_sign, args.antisym_eta, args.erase_phase_for_strict_sym
        ))
        na_out = hk.face_normal(model, fa, 'outward')
        nb_out = hk.face_normal(model, fb, 'outward')
        na_birth = hk.face_normal(model, fa, 'birth_order')
        R = p56.rotation_from_to(nb_out, -na_out)
        kb_r = R @ kb

        Q_a = ka + kb_r
        P_a = ka - kb_r
        Q_b = -(R.T @ Q_a)
        P_b = +(R.T @ P_a)

        ia, ib = idx[fa], idx[fb]
        W_Q[ia] += Q_a
        W_Q[ib] += Q_b
        W_P[ia] += P_a
        W_P[ib] += P_b
        cross_a = np.cross(Q_a, P_a)
        cross_b = -(R.T @ cross_a)
        W_comm[ia] += cross_a
        W_comm[ib] += cross_b

        J = block_J(R)
        C = block_C(R)
        q_block = block_vec(Q_a, Q_b)
        p_block = block_vec(P_a, P_b)
        qn = float(np.linalg.norm(q_block))
        pn = float(np.linalg.norm(p_block))
        jn = frob(J)
        cn = frob(C)

        j_square = rel(J @ J + I6, frob(I6))
        c_square = rel(C @ C - I6, frob(I6))
        cjc = rel(C @ J @ C + J, jn)
        orth = rel(J.T @ J - I6, frob(I6))

        # Derived parity/conjugation check on the actual Q/P fields.  The sign convention
        # is the cochain sign carried by the orientation-reversing face-pair gluing:
        # Q is even after cochain compensation, P is odd.
        q_even = float(np.linalg.norm(Q_a + R @ Q_b)) / (float(np.linalg.norm(Q_a)) + float(np.linalg.norm(Q_b)) + EPS)
        p_odd = float(np.linalg.norm(P_a - R @ P_b)) / (float(np.linalg.norm(P_a)) + float(np.linalg.norm(P_b)) + EPS)

        # Does the local pair-exchange J map the actual Q-channel to the actual P-channel?
        # A positive derived complex/quadrature lock would require this residual to be small,
        # or an independently derived replacement for J.  This test does not import such a map.
        jq = J @ q_block
        jp = J @ p_block
        jq_to_p = float(np.linalg.norm(jq - p_block)) / (qn + pn + EPS)
        jp_to_minus_q = float(np.linalg.norm(jp + q_block)) / (qn + pn + EPS)

        comm = signed_commutator_stats(Q_a, P_a, na_birth, na_out)
        comm_abs.append(comm['comm_abs_area'])
        comm_signed_birth.append(comm['comm_signed_birth'])
        comm_signed_out.append(comm['comm_signed_outward'])

        # Parity exchange algebra diagnostic.  E is the oriented cochain exchange (J-like),
        # O is the unoriented pair conjugation (C-like), scaled by the actual Q/P sizes.
        qscale = float(np.linalg.norm(Q_a))
        pscale = float(np.linalg.norm(P_a))
        E = qscale * J
        O = pscale * C
        B = E @ O - O @ E
        grade_signed = float(np.trace(B.T @ G)) / ((frob(B) * G_norm) + EPS) if frob(B) > EPS else 0.0
        block_comm_grade.append(grade_signed)

        j_square_resids.append(j_square)
        c_square_resids.append(c_square)
        cjc_resids.append(cjc)
        orth_resids.append(orth)
        q_even_resids.append(q_even)
        p_odd_resids.append(p_odd)
        jq_to_p_resids.append(jq_to_p)
        jp_to_minus_q_resids.append(jp_to_minus_q)

        pair_rows.append({
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
            'ka_norm': float(np.linalg.norm(ka)),
            'kb_reversed_norm': float(np.linalg.norm(kb_r)),
            'transport_cosine_ka_kb_reversed': float(np.dot(ka, kb_r) / ((np.linalg.norm(ka) * np.linalg.norm(kb_r)) + EPS)),
            'Q_even_norm': float(np.linalg.norm(Q_a)),
            'P_odd_norm': float(np.linalg.norm(P_a)),
            'P_over_Q_norm_ratio': safe_ratio(float(np.linalg.norm(P_a)), float(np.linalg.norm(Q_a))),
            'local_J_square_minus_neg_I_resid': j_square,
            'local_C_square_minus_I_resid': c_square,
            'local_CJC_plus_J_resid': cjc,
            'local_J_orthogonality_resid': orth,
            'Q_even_pair_conjugation_resid': q_even,
            'P_odd_pair_conjugation_resid': p_odd,
            'JQ_to_P_resid': jq_to_p,
            'JP_to_minusQ_resid': jp_to_minus_q,
            'block_exchange_comm_grade_signed': grade_signed,
            **comm,
        })

    Qm = p58.channel_metrics(model, faces, H, W_Q, 'Q_even')
    Pm = p58.channel_metrics(model, faces, H, W_P, 'P_odd')
    Cm = p58.channel_metrics(model, faces, H, W_comm, 'QP_comm')

    comm_abs_total = float(np.sum(comm_abs)) + EPS
    metrics = {
        'beta0': topo['beta0'], 'beta1': topo['beta1'], 'beta2': topo['beta2'], 'beta3': topo['beta3'],
        'harmonic_dim_real': int(H.shape[1]) if H.ndim == 2 else 0,
        'applied_pair_count': len(pair_rows),
        'decision_used_delta_beta_any': any(str(r.get('decision_used_delta_beta', '')).lower() == 'true' for r in pair_rows),
        'measured_delta_beta2_sum': sum(int(float(r.get('measured_delta_beta2', 0) or 0)) for r in pair_rows),
        'mean_transport_cosine_ka_kb_reversed': float(np.mean([r['transport_cosine_ka_kb_reversed'] for r in pair_rows])) if pair_rows else 0.0,
        'mean_P_over_Q_norm_ratio': float(np.mean([r['P_over_Q_norm_ratio'] for r in pair_rows])) if pair_rows else 0.0,
        'mean_local_J_square_resid': float(np.mean(j_square_resids)) if j_square_resids else 0.0,
        'mean_local_C_square_resid': float(np.mean(c_square_resids)) if c_square_resids else 0.0,
        'mean_local_CJC_plus_J_resid': float(np.mean(cjc_resids)) if cjc_resids else 0.0,
        'mean_local_J_orthogonality_resid': float(np.mean(orth_resids)) if orth_resids else 0.0,
        'mean_Q_even_pair_conjugation_resid': float(np.mean(q_even_resids)) if q_even_resids else 0.0,
        'mean_P_odd_pair_conjugation_resid': float(np.mean(p_odd_resids)) if p_odd_resids else 0.0,
        'mean_JQ_to_P_resid': float(np.mean(jq_to_p_resids)) if jq_to_p_resids else 0.0,
        'mean_JP_to_minusQ_resid': float(np.mean(jp_to_minus_q_resids)) if jp_to_minus_q_resids else 0.0,
        'comm_abs_total': comm_abs_total - EPS,
        'comm_signed_birth_over_abs_sum_ratio': float(np.sum(comm_signed_birth)) / comm_abs_total,
        'comm_signed_outward_over_abs_sum_ratio': float(np.sum(comm_signed_out)) / comm_abs_total,
        'mean_abs_comm_signed_birth_over_abs': float(np.mean([abs(r['comm_signed_birth_over_abs']) for r in pair_rows])) if pair_rows else 0.0,
        'mean_abs_comm_signed_outward_over_abs': float(np.mean([abs(r['comm_signed_outward_over_abs']) for r in pair_rows])) if pair_rows else 0.0,
        'mean_block_exchange_comm_grade_signed': float(np.mean(block_comm_grade)) if block_comm_grade else 0.0,
        'mean_abs_block_exchange_comm_grade_signed': float(np.mean(np.abs(block_comm_grade))) if block_comm_grade else 0.0,
    }
    metrics.update(Qm)
    metrics.update(Pm)
    metrics.update(Cm)

    face_rows: List[dict] = []
    qn = np.linalg.norm(W_Q, axis=1) if len(faces) else np.array([])
    pn = np.linalg.norm(W_P, axis=1) if len(faces) else np.array([])
    cn = np.linalg.norm(W_comm, axis=1) if len(faces) else np.array([])
    for i, f in enumerate(faces):
        if max(float(qn[i]), float(pn[i]), float(cn[i])) <= 0:
            continue
        face_rows.append({
            'face': str(list(f)),
            'birth_orders': str([model.nodes[v].birth_order for v in f]),
            'birth_times': str([model.nodes[v].birth_time for v in f]),
            'Q_even_norm': float(qn[i]),
            'P_odd_norm': float(pn[i]),
            'QP_comm_norm': float(cn[i]),
        })
    face_rows.sort(key=lambda r: (r['QP_comm_norm'], r['P_odd_norm']), reverse=True)
    return metrics, pair_rows, face_rows


def build_variant(variant: str, args: argparse.Namespace, out: Path):
    model = nl.build_model(variant, args)
    model.grow(args.max_level)
    baseline_K = core.build_dynamic_outward_ngf_complex(model)
    baseline_metrics = core.full_metrics(model, baseline_K, args.source)
    auto_K, birth_log, pairing_log, candidate_sample, scans = nl.build_nonlinear_auto_complex(model, args, out / variant, variant)
    auto_metrics = core.full_metrics(model, auto_K, args.source)
    return model, auto_K, baseline_metrics, auto_metrics, birth_log, pairing_log, candidate_sample, scans


def run_variant(variant: str, phase_sign: int, args: argparse.Namespace, out: Path) -> dict:
    vname = f'{variant}_phase{phase_sign:+d}'.replace('+', 'plus').replace('-', 'minus')
    vout = out / vname
    vout.mkdir(parents=True, exist_ok=True)
    local_args = argparse.Namespace(**vars(args))
    local_args.phase_sign = phase_sign
    model, K, baseline, auto, birth_log, pairing_log, candidate_sample, scans = build_variant(variant, local_args, out)
    metrics, pair_rows, face_rows = adjoint_pairing_metrics(model, K, pairing_log, local_args)
    write_csv(vout / 'birth_geometry_log.csv', birth_log)
    write_csv(vout / 'nonlinear_pairing_cascade_log.csv', pairing_log)
    write_csv(vout / 'adjoint_pairing_pair_rows.csv', pair_rows)
    write_csv(vout / 'adjoint_pairing_face_support_top.csv', face_rows[:args.keep_top_faces])
    summary = {
        'variant': variant,
        'variant_phase': vname,
        'phase_sign': phase_sign,
        'max_level': args.max_level,
        'source': args.source,
        'antisym_eta': args.antisym_eta,
        'erase_phase_for_strict_sym': args.erase_phase_for_strict_sym,
        'baseline_metrics_legacy_core': baseline,
        'auto_metrics_legacy_core': auto,
        'adjoint_pairing_metrics': metrics,
        'automatic_pairings_applied': sum(1 for x in pairing_log if x.get('applied')),
        'automatic_pairing_attempts_logged': len(pairing_log),
        'births_with_cascade_logs': scans,
        'interpretation_flags': {
            'beta2_opened': auto['beta2'] > baseline['beta2'],
            'Q_even_harmonic_positive': metrics['Q_even_harmonic_ratio'] > args.harmonic_positive_threshold,
            'P_odd_harmonic_positive': metrics['P_odd_harmonic_ratio'] > args.harmonic_positive_threshold,
            'local_pair_J_square_minus_I_positive': metrics['mean_local_J_square_resid'] < args.residual_threshold,
            'local_pair_C_square_positive': metrics['mean_local_C_square_resid'] < args.residual_threshold,
            'local_pair_C_anticommutes_with_J': metrics['mean_local_CJC_plus_J_resid'] < args.residual_threshold,
            'actual_QP_locked_by_local_J': metrics['mean_JQ_to_P_resid'] < args.lock_residual_threshold and metrics['mean_JP_to_minusQ_resid'] < args.lock_residual_threshold,
            'comm_signed_birth_nontrivial': abs(metrics['comm_signed_birth_over_abs_sum_ratio']) > args.signed_comm_threshold,
            'decision_used_delta_beta_any': metrics['decision_used_delta_beta_any'],
        },
    }
    (vout / 'variant_adjoint_pairing_summary.json').write_text(json.dumps(summary, indent=2), encoding='utf-8')
    return summary


def flip_comparisons(rows: List[dict]) -> List[dict]:
    by_variant = {}
    for r in rows:
        by_variant.setdefault(r['variant'], {})[r['phase_sign']] = r
    comps = []
    for variant, d in sorted(by_variant.items()):
        if 1 not in d or -1 not in d:
            continue
        p = d[1]['adjoint_pairing_metrics']
        m = d[-1]['adjoint_pairing_metrics']
        def flip_score(key: str) -> float:
            denom = abs(float(p[key])) + abs(float(m[key])) + EPS
            return (float(p[key]) + float(m[key])) / denom
        comps.append({
            'variant': variant,
            'comm_signed_birth_plus': p['comm_signed_birth_over_abs_sum_ratio'],
            'comm_signed_birth_minus': m['comm_signed_birth_over_abs_sum_ratio'],
            'comm_signed_birth_flip_score_zero_if_perfect_flip': flip_score('comm_signed_birth_over_abs_sum_ratio'),
            'block_comm_grade_plus': p['mean_block_exchange_comm_grade_signed'],
            'block_comm_grade_minus': m['mean_block_exchange_comm_grade_signed'],
            'block_comm_grade_flip_score_zero_if_perfect_flip': flip_score('mean_block_exchange_comm_grade_signed'),
            'JQ_to_P_plus': p['mean_JQ_to_P_resid'],
            'JQ_to_P_minus': m['mean_JQ_to_P_resid'],
            'P_harm_plus': p['P_odd_harmonic_ratio'],
            'P_harm_minus': m['P_odd_harmonic_ratio'],
        })
    return comps


def write_comparative(out: Path, rows: List[dict], comps: List[dict]) -> None:
    flat = []
    for r in rows:
        a = r['auto_metrics_legacy_core']; m = r['adjoint_pairing_metrics']
        flat.append({
            'variant_phase': r['variant_phase'],
            'variant': r['variant'],
            'phase_sign': r['phase_sign'],
            'beta0': a['beta0'], 'beta1': a['beta1'], 'beta2': a['beta2'], 'beta3': a['beta3'],
            'pairings': r['automatic_pairings_applied'],
            'Q_even_harmonic_ratio': m['Q_even_harmonic_ratio'],
            'P_odd_harmonic_ratio': m['P_odd_harmonic_ratio'],
            'QP_comm_harmonic_ratio': m['QP_comm_harmonic_ratio'],
            'mean_local_J_square_resid': m['mean_local_J_square_resid'],
            'mean_local_C_square_resid': m['mean_local_C_square_resid'],
            'mean_local_CJC_plus_J_resid': m['mean_local_CJC_plus_J_resid'],
            'mean_Q_even_pair_conjugation_resid': m['mean_Q_even_pair_conjugation_resid'],
            'mean_P_odd_pair_conjugation_resid': m['mean_P_odd_pair_conjugation_resid'],
            'mean_JQ_to_P_resid': m['mean_JQ_to_P_resid'],
            'mean_JP_to_minusQ_resid': m['mean_JP_to_minusQ_resid'],
            'comm_signed_birth_over_abs_sum_ratio': m['comm_signed_birth_over_abs_sum_ratio'],
            'comm_signed_outward_over_abs_sum_ratio': m['comm_signed_outward_over_abs_sum_ratio'],
            'mean_block_exchange_comm_grade_signed': m['mean_block_exchange_comm_grade_signed'],
            'mean_abs_block_exchange_comm_grade_signed': m['mean_abs_block_exchange_comm_grade_signed'],
            'decision_used_delta_beta_any': m['decision_used_delta_beta_any'],
        })
    write_csv(out / 'comparative_adjoint_pairing_summary.csv', flat)
    write_csv(out / 'phase_flip_comparison.csv', comps)


def make_docs(summary: dict, comps: List[dict]) -> tuple[str, str, str, str]:
    rows = summary['variant_rows']
    table_lines = [
        '| variant/phase | beta | pairs | Q harm | P harm | J^2+I resid | C^2-I resid | CJC+J resid | JQ->P resid | comm signed birth | block comm grade | used dBeta? |',
        '|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|',
    ]
    for r in rows:
        a = r['auto_metrics_legacy_core']; m = r['adjoint_pairing_metrics']
        table_lines.append(
            f"| {r['variant_phase']} | ({a['beta0']},{a['beta1']},{a['beta2']},{a['beta3']}) | "
            f"{r['automatic_pairings_applied']} | {m['Q_even_harmonic_ratio']:.6g} | {m['P_odd_harmonic_ratio']:.6g} | "
            f"{m['mean_local_J_square_resid']:.3g} | {m['mean_local_C_square_resid']:.3g} | {m['mean_local_CJC_plus_J_resid']:.3g} | "
            f"{m['mean_JQ_to_P_resid']:.6g} | {m['comm_signed_birth_over_abs_sum_ratio']:.6g} | "
            f"{m['mean_block_exchange_comm_grade_signed']:.6g} | {m['decision_used_delta_beta_any']} |"
        )
    table = '\n'.join(table_lines)
    comp_lines = [
        '| variant | comm birth + | comm birth - | comm flip-score | block grade + | block grade - | grade flip-score |',
        '|---|---:|---:|---:|---:|---:|---:|',
    ]
    for c in comps:
        comp_lines.append(
            f"| {c['variant']} | {c['comm_signed_birth_plus']:.6g} | {c['comm_signed_birth_minus']:.6g} | "
            f"{c['comm_signed_birth_flip_score_zero_if_perfect_flip']:.6g} | {c['block_comm_grade_plus']:.6g} | "
            f"{c['block_comm_grade_minus']:.6g} | {c['block_comm_grade_flip_score_zero_if_perfect_flip']:.6g} |"
        )
    comp_table = '\n'.join(comp_lines)
    smd = f"""# SUMMARY — pairing quadrature adjoint-pairing gate

Model label:
CNNA growing primal simplicial complex with deterministic sequential provenance growth,
nonlinear asymmetry-gated complement pairing, and directed antisymmetric birth-transport
operators.  This package tests a real pair-conjugation/quadrature structure.  It does not
claim complex scalars, a physical adjoint, positivity, norm, or C*-structure.

Core real pair operators for each actual glued face-pair:

```text
C_pair = [[0, R], [ R^T, 0]]       C_pair^2 = +I
J_pair = [[0, R], [-R^T, 0]]       J_pair^2 = -I
C_pair J_pair C_pair = -J_pair
```

Here R is the orientation-reversing transport determined by the actual face-pair gluing.
This is a local cochain-pair algebra, not a global i.

{table}

## Phase-sign flip comparison

{comp_table}

Conservative reading:
If the local pair algebra is exact but `JQ->P` remains large and signed commutators do not
flip, then the test found a local real conjugation/J scaffold but not a derived Q/P
complex lock.
"""
    rmd = f"""# RESULTS — pairing quadrature adjoint-pairing gate

## Comparative table

{table}

## Phase-sign flip comparison

{comp_table}

## Interpretation protocol

This package separates four claims:

```text
1. Pair conjugation exists:       C^2 = +I.
2. Oriented cochain exchange:     J^2 = -I locally on each glued face-pair.
3. Conjugation compatibility:     C J C = -J.
4. Actual Q/P lock:               J(Q) = P and J(P) = -Q.
```

Only (1)-(3) are structural local pair-algebra tests.  Claim (4) is the nontrivial
quadrature-lock gate.  A nonzero commutator magnitude is not enough; signed/kappa-like
behavior is logged separately.

## Anti-smuggling conditions

- no `i`, no imported `J`, no Hodge star, no physical adjoint, no positivity;
- no final `sym(M)` in the directed birth-transport operator;
- `J_pair` is the oriented cochain-exchange map induced by the actual paired faces;
- local `J_pair^2=-I` is not interpreted as a global complex structure unless Q/P lock
  and coherence gates also pass;
- `decision_used_delta_beta_any` must remain false.
"""
    audit = """# SOURCE AUDIT

The previous package corrected the magnitude-only quadrature area test by measuring
signed area.  The signed spatial area did not flip coherently.  This package therefore
moves from spatial area to real pair algebra:

- pair exchange/conjugation C;
- oriented cochain exchange J_pair;
- their algebraic residuals;
- whether the actual derived Q/P channels are locked by that local J_pair;
- signed real commutator diagnostics.

This tests the user's hypothesis that the missing structure is closer to an
erzeuger/vernichter or quadrature-adjunction split than to a spatial face orientation.
"""
    readme = """# Pairing quadrature adjoint-pairing gate

Run:

```bash
python3 test_pairing_quadrature_adjoint_pairing_gate.py
```

Outputs:

- comparative_summary.json
- comparative_adjoint_pairing_summary.csv
- phase_flip_comparison.csv
- RESULTS.md
- SUMMARY.md
- SOURCE_AUDIT.md
- per-variant pair/face logs

Positive result would require more than local `J_pair^2=-I`: the actual Q/P channels
must be locked by the pair algebra and the signed commutator should transform coherently
under the phase/kappa proxy while strict_sym remains killed.
"""
    return smd, rmd, audit, readme


def package(out: Path, zip_path: Path) -> None:
    files = [
        Path(__file__).name,
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
    ap.add_argument('--residual-threshold', type=float, default=1e-9)
    ap.add_argument('--lock-residual-threshold', type=float, default=0.1)
    ap.add_argument('--signed-comm-threshold', type=float, default=0.15)
    ap.add_argument('--antisym-eta', type=float, default=1.0)
    ap.add_argument('--erase-phase-for-strict-sym', action='store_true', default=True)
    ap.add_argument('--phase-signs', nargs='*', type=int, default=[1, -1], choices=[-1, 1])
    ap.add_argument('--variants', nargs='*', default=['real_growth', 'strict_symmetrized_control', 'no_backreaction'])
    ap.add_argument('--out', default='pairing_quadrature_adjoint_pairing_out_L2')
    ap.add_argument('--zip', default='cnna_pairing_quadrature_adjoint_pairing_gate_pkg_L2.zip')
    args = ap.parse_args()

    out = Path(args.out)
    if out.exists():
        shutil.rmtree(out)
    out.mkdir(parents=True)
    rows = []
    for v in args.variants:
        for s in args.phase_signs:
            rows.append(run_variant(v, s, args, out))
    comps = flip_comparisons(rows)
    summary = {'args': vars(args), 'variant_rows': rows, 'phase_flip_comparisons': comps}
    (out / 'comparative_summary.json').write_text(json.dumps(summary, indent=2), encoding='utf-8')
    write_comparative(out, rows, comps)
    smd, rmd, audit, readme = make_docs(summary, comps)
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
                'variant_phase': r['variant_phase'],
                'auto_beta': [r['auto_metrics_legacy_core'][f'beta{i}'] for i in range(4)],
                'pairings': r['automatic_pairings_applied'],
                'Q_harm': r['adjoint_pairing_metrics']['Q_even_harmonic_ratio'],
                'P_harm': r['adjoint_pairing_metrics']['P_odd_harmonic_ratio'],
                'J_square_resid': r['adjoint_pairing_metrics']['mean_local_J_square_resid'],
                'C_square_resid': r['adjoint_pairing_metrics']['mean_local_C_square_resid'],
                'CJC_plus_J_resid': r['adjoint_pairing_metrics']['mean_local_CJC_plus_J_resid'],
                'JQ_to_P_resid': r['adjoint_pairing_metrics']['mean_JQ_to_P_resid'],
                'comm_signed_birth': r['adjoint_pairing_metrics']['comm_signed_birth_over_abs_sum_ratio'],
                'block_comm_grade_signed': r['adjoint_pairing_metrics']['mean_block_exchange_comm_grade_signed'],
                'decision_used_delta_beta_any': r['adjoint_pairing_metrics']['decision_used_delta_beta_any'],
            } for r in rows
        ],
        'phase_flip_comparisons': comps,
    }, indent=2))


if __name__ == '__main__':
    main()
