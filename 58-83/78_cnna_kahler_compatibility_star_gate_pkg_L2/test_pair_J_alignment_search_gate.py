#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import shutil
import zipfile
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np

import cnna_non_shelling_core as core
import test_nonlinear_asymmetry_cascade_growth as nl
import test_harmonic_k_orientation_kappa_gate as hk
import test_pairing_transport_antisym_birth_coherence_gate as p56
import test_pairing_quadrature_split_symplectic_defect_gate as p58
import test_pairing_quadrature_adjoint_pairing_gate as p60

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


def norm(x: np.ndarray) -> float:
    return float(np.linalg.norm(x))


def block(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    return np.concatenate([np.asarray(a, dtype=float), np.asarray(b, dtype=float)])


def safe_unit(v: np.ndarray) -> np.ndarray:
    n = norm(v)
    if n < EPS:
        return np.zeros_like(v)
    return v / n


def J_lock_residual(J: np.ndarray, q: np.ndarray, p: np.ndarray) -> Dict[str, float]:
    qn = norm(q)
    pn = norm(p)
    den = qn + pn + EPS
    if qn < EPS or pn < EPS:
        return {
            'JQ_to_P_resid': 0.0 if qn < EPS and pn < EPS else 1.0,
            'JP_to_minusQ_resid': 0.0 if qn < EPS and pn < EPS else 1.0,
            'J_lock_mean_resid': 0.0 if qn < EPS and pn < EPS else 1.0,
            'J_lock_max_resid': 0.0 if qn < EPS and pn < EPS else 1.0,
            'Q_norm': qn,
            'P_norm': pn,
        }
    r1 = norm(J @ q - p) / den
    r2 = norm(J @ p + q) / den
    return {
        'JQ_to_P_resid': r1,
        'JP_to_minusQ_resid': r2,
        'J_lock_mean_resid': 0.5 * (r1 + r2),
        'J_lock_max_resid': max(r1, r2),
        'Q_norm': qn,
        'P_norm': pn,
    }


def scalar_matched_pair(J: np.ndarray, q: np.ndarray, p: np.ndarray) -> Tuple[np.ndarray, np.ndarray, dict]:
    # Not a fitted rotation: only a derived positive scale based on norms of Jq and p.
    # This checks whether the obstruction is only a Q/P amplitude mismatch.
    jq = J @ q
    jp = J @ p
    alpha = norm(jq) / (norm(p) + EPS)
    beta = norm(jp) / (norm(q) + EPS)
    p_scaled = alpha * p
    q_scaled = beta * q
    return q_scaled, p_scaled, {'alpha_P_scale_from_norms': alpha, 'beta_Q_scale_from_norms': beta}


def project_face_field(H: np.ndarray, W: np.ndarray) -> np.ndarray:
    if H.size == 0 or H.shape[1] == 0:
        return np.zeros_like(W)
    return H @ (H.T @ W)


def kappa_like_flip_score(plus: float, minus: float) -> float:
    return (plus + minus) / (abs(plus) + abs(minus) + EPS)


def build_variant(variant: str, args: argparse.Namespace, out: Path):
    model = nl.build_model(variant, args)
    model.grow(args.max_level)
    baseline_K = core.build_dynamic_outward_ngf_complex(model)
    baseline_metrics = core.full_metrics(model, baseline_K, args.source)
    auto_K, birth_log, pairing_log, candidate_sample, scans = nl.build_nonlinear_auto_complex(model, args, out / variant, variant)
    auto_metrics = core.full_metrics(model, auto_K, args.source)
    return model, auto_K, baseline_metrics, auto_metrics, birth_log, pairing_log, candidate_sample, scans


def local_pair_data(model, K, pairing_log, args):
    faces = K.faces()
    idx = {tuple(f): i for i, f in enumerate(faces)}
    H, _ = hk.harmonic_basis_faces(K)
    W_Q = np.zeros((len(faces), 3), dtype=float)
    W_P = np.zeros((len(faces), 3), dtype=float)
    W_C = np.zeros((len(faces), 3), dtype=float)
    pairs = []
    for k, log in enumerate(pairing_log):
        if not log.get('applied'):
            continue
        fa, fb = p56.parse_pair_faces(log)
        if fa is None or fb is None or fa not in idx or fb not in idx:
            continue
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
        ia, ib = idx[fa], idx[fb]
        W_Q[ia] += Q_a; W_Q[ib] += Q_b
        W_P[ia] += P_a; W_P[ib] += P_b
        W_C[ia] += C_a; W_C[ib] += C_b
        pairs.append({
            'pair_index': k,
            'log': log,
            'fa': fa,
            'fb': fb,
            'ia': ia,
            'ib': ib,
            'R': R,
            'J': p60.block_J(R),
            'C': p60.block_C(R),
            'ka': ka,
            'kb': kb,
            'kb_r': kb_r,
            'Q_a': Q_a,
            'P_a': P_a,
            'Q_b': Q_b,
            'P_b': P_b,
            'C_a': C_a,
            'C_b': C_b,
            'na_birth': na_birth,
            'na_out': na_out,
        })
    return faces, H, W_Q, W_P, W_C, pairs


def candidate_blocks(pair: dict, HQ: np.ndarray, HP: np.ndarray, HC: np.ndarray) -> List[Tuple[str, np.ndarray, np.ndarray, dict]]:
    J = pair['J']; C = pair['C']
    I = np.eye(6)
    ia, ib = pair['ia'], pair['ib']
    q = block(pair['Q_a'], pair['Q_b'])
    p = block(pair['P_a'], pair['P_b'])
    c = block(pair['C_a'], pair['C_b'])
    hq = block(HQ[ia], HQ[ib])
    hp = block(HP[ia], HP[ib])
    hc = block(HC[ia], HC[ib])
    raw = block(pair['ka'], pair['kb'])
    transported_raw = block(pair['ka'], pair['kb_r'])
    E = 0.5 * (I + C)
    O = 0.5 * (I - C)
    even_raw = E @ raw
    odd_raw = O @ raw
    even_trans = E @ transported_raw
    odd_trans = O @ transported_raw
    out: List[Tuple[str, np.ndarray, np.ndarray, dict]] = []
    out.append(('raw_QP', q, p, {}))
    out.append(('unit_QP', safe_unit(q), safe_unit(p), {}))
    qs, ps, extra = scalar_matched_pair(J, q, p)
    out.append(('norm_matched_QP', qs, ps, extra))
    out.append(('harmonic_QP', hq, hp, {}))
    out.append(('unit_harmonic_QP', safe_unit(hq), safe_unit(hp), {}))
    qhs, phs, extra_h = scalar_matched_pair(J, hq, hp)
    out.append(('norm_matched_harmonic_QP', qhs, phs, extra_h))
    out.append(('Q_comm', q, c, {}))
    out.append(('P_comm', p, c, {}))
    out.append(('harmonic_Q_comm', hq, hc, {}))
    out.append(('harmonic_P_comm', hp, hc, {}))
    out.append(('C_eigen_raw_even_odd', even_raw, odd_raw, {}))
    out.append(('C_eigen_transported_even_odd', even_trans, odd_trans, {}))
    out.append(('unit_C_eigen_raw_even_odd', safe_unit(even_raw), safe_unit(odd_raw), {}))
    out.append(('unit_C_eigen_transported_even_odd', safe_unit(even_trans), safe_unit(odd_trans), {}))
    # Non-circular check: use J only as a diagnostic, not to define the target P.  Therefore
    # no candidate of the form (v, Jv) is included.
    return out


def alignment_search_metrics(model, K, pairing_log, args):
    topo = core.topology(K)
    faces, H, W_Q, W_P, W_C, pairs = local_pair_data(model, K, pairing_log, args)
    HQ = project_face_field(H, W_Q)
    HP = project_face_field(H, W_P)
    HC = project_face_field(H, W_C)
    Qm = p58.channel_metrics(model, faces, H, W_Q, 'Q_even')
    Pm = p58.channel_metrics(model, faces, H, W_P, 'P_odd')
    Cm = p58.channel_metrics(model, faces, H, W_C, 'QP_comm')

    pair_rows: List[dict] = []
    cand_rows: List[dict] = []
    by_candidate: Dict[str, List[dict]] = {}
    signed_comm_birth = []
    abs_comm = []
    for pair in pairs:
        log = pair['log']
        J = pair['J']
        cross = np.cross(pair['Q_a'], pair['P_a'])
        comm_abs = norm(cross)
        comm_birth = float(np.dot(cross, pair['na_birth']))
        abs_comm.append(comm_abs)
        signed_comm_birth.append(comm_birth)
        candidate_results = []
        for name, qv, pv, extra in candidate_blocks(pair, HQ, HP, HC):
            res = J_lock_residual(J, qv, pv)
            row = {
                'pair_index': pair['pair_index'],
                'candidate': name,
                'JQ_to_P_resid': res['JQ_to_P_resid'],
                'JP_to_minusQ_resid': res['JP_to_minusQ_resid'],
                'J_lock_mean_resid': res['J_lock_mean_resid'],
                'J_lock_max_resid': res['J_lock_max_resid'],
                'Q_norm': res['Q_norm'],
                'P_norm': res['P_norm'],
                **extra,
            }
            cand_rows.append(row)
            by_candidate.setdefault(name, []).append(row)
            candidate_results.append(row)
        best = min(candidate_results, key=lambda r: (r['J_lock_mean_resid'], r['J_lock_max_resid'])) if candidate_results else None
        pair_rows.append({
            'pair_index': pair['pair_index'],
            'face_a': str(list(pair['fa'])),
            'face_b': str(list(pair['fb'])),
            'A_rank_score': log.get('A_rank_score', ''),
            'A_invariant': log.get('A_invariant', ''),
            'directed_imbalance': log.get('directed_imbalance', ''),
            'transverse_complementarity': log.get('transverse_complementarity', ''),
            'measured_delta_beta1': log.get('measured_delta_beta1', ''),
            'measured_delta_beta2': log.get('measured_delta_beta2', ''),
            'decision_used_delta_beta': log.get('decision_used_delta_beta', ''),
            'transport_cosine_ka_kb_reversed': float(np.dot(pair['ka'], pair['kb_r']) / ((norm(pair['ka']) * norm(pair['kb_r'])) + EPS)),
            'comm_abs_area': comm_abs,
            'comm_signed_birth': comm_birth,
            'comm_signed_birth_over_abs': comm_birth / (comm_abs + EPS),
            'best_candidate': best['candidate'] if best else '',
            'best_J_lock_mean_resid': best['J_lock_mean_resid'] if best else 0.0,
            'best_J_lock_max_resid': best['J_lock_max_resid'] if best else 0.0,
        })

    candidate_summary: List[dict] = []
    for name, rows in sorted(by_candidate.items()):
        if not rows:
            continue
        means = {k: float(np.mean([float(r[k]) for r in rows])) for k in ['JQ_to_P_resid','JP_to_minusQ_resid','J_lock_mean_resid','J_lock_max_resid','Q_norm','P_norm']}
        maxs = {f'max_{k}': float(np.max([float(r[k]) for r in rows])) for k in ['J_lock_mean_resid','J_lock_max_resid']}
        mins = {f'min_{k}': float(np.min([float(r[k]) for r in rows])) for k in ['J_lock_mean_resid','J_lock_max_resid']}
        candidate_summary.append({'candidate': name, 'pair_count': len(rows), **means, **maxs, **mins})
    candidate_summary.sort(key=lambda r: (r['J_lock_mean_resid'], r['J_lock_max_resid']))
    best_global = candidate_summary[0] if candidate_summary else {}
    comm_abs_total = float(np.sum(abs_comm)) + EPS
    best_pair_mean = float(np.mean([r['best_J_lock_mean_resid'] for r in pair_rows])) if pair_rows else 0.0
    best_pair_max = float(np.max([r['best_J_lock_max_resid'] for r in pair_rows])) if pair_rows else 0.0
    metrics = {
        'beta0': topo['beta0'], 'beta1': topo['beta1'], 'beta2': topo['beta2'], 'beta3': topo['beta3'],
        'harmonic_dim_real': int(H.shape[1]) if H.ndim == 2 else 0,
        'applied_pair_count': len(pairs),
        'candidate_count': len(candidate_summary),
        'decision_used_delta_beta_any': any(str(p['log'].get('decision_used_delta_beta', '')).lower() == 'true' for p in pairs),
        'measured_delta_beta2_sum': sum(int(float(p['log'].get('measured_delta_beta2', 0) or 0)) for p in pairs),
        'best_global_candidate': best_global.get('candidate', ''),
        'best_global_J_lock_mean_resid': best_global.get('J_lock_mean_resid', 0.0),
        'best_global_J_lock_max_resid': best_global.get('J_lock_max_resid', 0.0),
        'best_per_pair_mean_J_lock_resid': best_pair_mean,
        'best_per_pair_max_J_lock_resid': best_pair_max,
        'comm_signed_birth_over_abs_sum_ratio': float(np.sum(signed_comm_birth)) / comm_abs_total,
        'mean_abs_comm_signed_birth_over_abs': float(np.mean([abs(r['comm_signed_birth_over_abs']) for r in pair_rows])) if pair_rows else 0.0,
    }
    metrics.update(Qm); metrics.update(Pm); metrics.update(Cm)
    return metrics, pair_rows, cand_rows, candidate_summary


def run_variant(variant: str, phase_sign: int, args: argparse.Namespace, out: Path) -> dict:
    vname = f'{variant}_phase{phase_sign:+d}'.replace('+', 'plus').replace('-', 'minus')
    vout = out / vname
    vout.mkdir(parents=True, exist_ok=True)
    local_args = argparse.Namespace(**vars(args))
    local_args.phase_sign = phase_sign
    model, K, baseline, auto, birth_log, pairing_log, candidate_sample, scans = build_variant(variant, local_args, out)
    metrics, pair_rows, cand_rows, candidate_summary = alignment_search_metrics(model, K, pairing_log, local_args)
    write_csv(vout / 'birth_geometry_log.csv', birth_log)
    write_csv(vout / 'nonlinear_pairing_cascade_log.csv', pairing_log)
    write_csv(vout / 'alignment_pair_rows.csv', pair_rows)
    write_csv(vout / 'alignment_candidate_rows.csv', cand_rows)
    write_csv(vout / 'alignment_candidate_summary.csv', candidate_summary)
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
        'alignment_metrics': metrics,
        'candidate_summary': candidate_summary,
        'automatic_pairings_applied': sum(1 for x in pairing_log if x.get('applied')),
        'automatic_pairing_attempts_logged': len(pairing_log),
        'births_with_cascade_logs': scans,
        'interpretation_flags': {
            'beta2_opened': auto['beta2'] > baseline['beta2'],
            'best_alignment_passed': metrics['best_per_pair_mean_J_lock_resid'] < args.lock_residual_threshold and metrics['best_per_pair_max_J_lock_resid'] < args.lock_max_threshold,
            'global_candidate_alignment_passed': metrics['best_global_J_lock_mean_resid'] < args.lock_residual_threshold and metrics['best_global_J_lock_max_resid'] < args.lock_max_threshold,
            'signed_comm_birth_nontrivial': abs(metrics['comm_signed_birth_over_abs_sum_ratio']) > args.signed_comm_threshold,
            'decision_used_delta_beta_any': metrics['decision_used_delta_beta_any'],
        },
    }
    (vout / 'variant_alignment_search_summary.json').write_text(json.dumps(summary, indent=2), encoding='utf-8')
    return summary


def flip_comparisons(rows: List[dict]) -> List[dict]:
    by = {}
    for r in rows:
        by.setdefault(r['variant'], {})[r['phase_sign']] = r
    comps = []
    for variant, d in sorted(by.items()):
        if 1 not in d or -1 not in d:
            continue
        p = d[1]['alignment_metrics']; m = d[-1]['alignment_metrics']
        comps.append({
            'variant': variant,
            'best_candidate_plus': p['best_global_candidate'],
            'best_candidate_minus': m['best_global_candidate'],
            'best_mean_plus': p['best_global_J_lock_mean_resid'],
            'best_mean_minus': m['best_global_J_lock_mean_resid'],
            'best_pair_mean_plus': p['best_per_pair_mean_J_lock_resid'],
            'best_pair_mean_minus': m['best_per_pair_mean_J_lock_resid'],
            'comm_signed_birth_plus': p['comm_signed_birth_over_abs_sum_ratio'],
            'comm_signed_birth_minus': m['comm_signed_birth_over_abs_sum_ratio'],
            'comm_signed_birth_flip_score_zero_if_perfect_flip': kappa_like_flip_score(p['comm_signed_birth_over_abs_sum_ratio'], m['comm_signed_birth_over_abs_sum_ratio']),
        })
    return comps


def write_comparative(out: Path, rows: List[dict], comps: List[dict]) -> None:
    flat = []
    cand_flat = []
    for r in rows:
        a = r['auto_metrics_legacy_core']; m = r['alignment_metrics']
        flat.append({
            'variant_phase': r['variant_phase'],
            'variant': r['variant'],
            'phase_sign': r['phase_sign'],
            'beta0': a['beta0'], 'beta1': a['beta1'], 'beta2': a['beta2'], 'beta3': a['beta3'],
            'pairings': r['automatic_pairings_applied'],
            'Q_even_harmonic_ratio': m['Q_even_harmonic_ratio'],
            'P_odd_harmonic_ratio': m['P_odd_harmonic_ratio'],
            'QP_comm_harmonic_ratio': m['QP_comm_harmonic_ratio'],
            'best_global_candidate': m['best_global_candidate'],
            'best_global_J_lock_mean_resid': m['best_global_J_lock_mean_resid'],
            'best_global_J_lock_max_resid': m['best_global_J_lock_max_resid'],
            'best_per_pair_mean_J_lock_resid': m['best_per_pair_mean_J_lock_resid'],
            'best_per_pair_max_J_lock_resid': m['best_per_pair_max_J_lock_resid'],
            'comm_signed_birth_over_abs_sum_ratio': m['comm_signed_birth_over_abs_sum_ratio'],
            'decision_used_delta_beta_any': m['decision_used_delta_beta_any'],
        })
        for c in r['candidate_summary']:
            cand_flat.append({'variant_phase': r['variant_phase'], 'variant': r['variant'], 'phase_sign': r['phase_sign'], **c})
    write_csv(out / 'comparative_alignment_search_summary.csv', flat)
    write_csv(out / 'comparative_candidate_summary.csv', cand_flat)
    write_csv(out / 'phase_flip_comparison.csv', comps)


def make_docs(summary: dict, comps: List[dict]) -> tuple[str, str, str, str]:
    rows = summary['variant_rows']
    table_lines = [
        '| variant/phase | beta | pairs | Q harm | P harm | best candidate | best mean resid | best max resid | per-pair best mean | comm signed birth | used dBeta? |',
        '|---|---:|---:|---:|---:|---|---:|---:|---:|---:|---:|',
    ]
    for r in rows:
        a = r['auto_metrics_legacy_core']; m = r['alignment_metrics']
        table_lines.append(
            f"| {r['variant_phase']} | ({a['beta0']},{a['beta1']},{a['beta2']},{a['beta3']}) | {r['automatic_pairings_applied']} | "
            f"{m['Q_even_harmonic_ratio']:.6g} | {m['P_odd_harmonic_ratio']:.6g} | {m['best_global_candidate']} | "
            f"{m['best_global_J_lock_mean_resid']:.6g} | {m['best_global_J_lock_max_resid']:.6g} | "
            f"{m['best_per_pair_mean_J_lock_resid']:.6g} | {m['comm_signed_birth_over_abs_sum_ratio']:.6g} | {m['decision_used_delta_beta_any']} |"
        )
    table = '\n'.join(table_lines)
    comp_lines = [
        '| variant | best + | best - | mean + | mean - | pair mean + | pair mean - | comm + | comm - | flip-score |',
        '|---|---|---|---:|---:|---:|---:|---:|---:|---:|',
    ]
    for c in comps:
        comp_lines.append(
            f"| {c['variant']} | {c['best_candidate_plus']} | {c['best_candidate_minus']} | {c['best_mean_plus']:.6g} | "
            f"{c['best_mean_minus']:.6g} | {c['best_pair_mean_plus']:.6g} | {c['best_pair_mean_minus']:.6g} | "
            f"{c['comm_signed_birth_plus']:.6g} | {c['comm_signed_birth_minus']:.6g} | {c['comm_signed_birth_flip_score_zero_if_perfect_flip']:.6g} |"
        )
    comp_table = '\n'.join(comp_lines)
    smd = f"""# SUMMARY — pair J-alignment search gate

Model label:
CNNA growing primal simplicial complex with deterministic sequential provenance growth,
nonlinear asymmetry-gated complement pairing, directed antisymmetric birth-transport
operators, and local pair-exchange algebra.

This package does **not** fit an arbitrary rotation and does **not** define P as JQ.  It tests
only already-derived candidate pairs:

```text
raw Q/P, unit Q/P, norm-matched Q/P,
harmonic Q/P, harmonic unit/norm-matched Q/P,
Q vs commutator, P vs commutator,
C-pair eigen even/odd projections of raw pair data.
```

Gate:
`J_pair(Q') ≈ P'` and `J_pair(P') ≈ -Q'` must hold for an allowed candidate, while
strict_sym remains killed.

{table}

## Phase-sign comparison

{comp_table}

Conservative reading: if all allowed candidates keep large lock residuals, the local pair
algebra exists but does not dynamically select the Q/P or a/a† split.
"""
    rmd = f"""# RESULTS — pair J-alignment search gate

## Comparative table

{table}

## Phase-sign comparison

{comp_table}

## Interpretation protocol

A positive result would require an allowed, derived candidate channel with small residuals:

```text
J_pair(Q') -> P'
J_pair(P') -> -Q'
```

This package deliberately rejects the circular construction `P' := J_pair(Q')`.  The question
is whether Q/P-like channels already present in the data align with the local pair J.

## Anti-smuggling conditions

- no `i`, no complex scalars;
- no imported Hodge star, positivity, physical adjoint, or norm axiom;
- no final `sym(M)` in the directed birth-transport operator;
- no arbitrary fitted rotation;
- no topology/H2/kappa used in move decisions;
- `decision_used_delta_beta_any` must remain false.
"""
    audit = """# SOURCE AUDIT

Previous gate found a local pair algebra:

```text
C_pair^2 = +I
J_pair^2 = -I
C_pair J_pair C_pair = -J_pair
```

but the actual raw Q/P channels were not locked by this J.  This package searches for
alignment only among derived channels already present in the pair data, not among arbitrary
rotations or user-chosen complex structures.
"""
    readme = """# Pair J-alignment search gate

Run:

```bash
python3 test_pair_J_alignment_search_gate.py
```

Outputs:

- comparative_summary.json
- comparative_alignment_search_summary.csv
- comparative_candidate_summary.csv
- phase_flip_comparison.csv
- RESULTS.md
- SUMMARY.md
- SOURCE_AUDIT.md
- per-variant candidate logs

Positive result requires a derived candidate pair aligned by local J_pair and killed by
strict_sym.  A local algebra alone is not sufficient.
"""
    return smd, rmd, audit, readme


def package(out: Path, zip_path: Path) -> None:
    files = [
        Path(__file__).name,
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
    ap.add_argument('--keep-top-candidates', type=int, default=80)
    ap.add_argument('--keep-top-faces', type=int, default=80)
    ap.add_argument('--harmonic-positive-threshold', type=float, default=1e-4)
    ap.add_argument('--antisym-eta', type=float, default=1.0)
    ap.add_argument('--erase-phase-for-strict-sym', action='store_true', default=True)
    ap.add_argument('--phases', nargs='*', type=int, default=[1, -1])
    ap.add_argument('--lock-residual-threshold', type=float, default=0.20)
    ap.add_argument('--lock-max-threshold', type=float, default=0.30)
    ap.add_argument('--signed-comm-threshold', type=float, default=0.15)
    ap.add_argument('--variants', nargs='*', default=['real_growth', 'strict_symmetrized_control', 'no_backreaction'])
    ap.add_argument('--out', default='pair_J_alignment_search_out_L2')
    ap.add_argument('--zip', default='cnna_pair_J_alignment_search_gate_pkg_L2.zip')
    args = ap.parse_args()

    out = Path(args.out)
    if out.exists():
        shutil.rmtree(out)
    out.mkdir(parents=True)
    rows = []
    for v in args.variants:
        for ph in args.phases:
            rows.append(run_variant(v, ph, args, out))
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
                'beta': [r['auto_metrics_legacy_core'][f'beta{i}'] for i in range(4)],
                'pairings': r['automatic_pairings_applied'],
                'best_candidate': r['alignment_metrics']['best_global_candidate'],
                'best_global_J_lock_mean_resid': r['alignment_metrics']['best_global_J_lock_mean_resid'],
                'best_per_pair_mean_J_lock_resid': r['alignment_metrics']['best_per_pair_mean_J_lock_resid'],
                'Q_harm': r['alignment_metrics']['Q_even_harmonic_ratio'],
                'P_harm': r['alignment_metrics']['P_odd_harmonic_ratio'],
                'decision_used_delta_beta_any': r['alignment_metrics']['decision_used_delta_beta_any'],
            } for r in rows
        ]
    }, indent=2))


if __name__ == '__main__':
    main()
