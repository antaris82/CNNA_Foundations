#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import math
import shutil
import zipfile
from pathlib import Path
from typing import List, Tuple

import numpy as np

import cnna_non_shelling_core as core
import test_nonlinear_asymmetry_cascade_growth as nl
import test_harmonic_k_orientation_kappa_gate as hk
import test_pairing_transport_antisym_birth_coherence_gate as p56

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


def skew_from_axial(w: np.ndarray) -> np.ndarray:
    wx, wy, wz = [float(x) for x in w]
    return np.array([[0.0, -wz, wy], [wz, 0.0, -wx], [-wy, wx, 0.0]], dtype=float)


def vector_field_projection(H: np.ndarray, W: np.ndarray) -> np.ndarray:
    if H.size == 0 or H.ndim != 2 or H.shape[1] == 0:
        return np.zeros_like(W)
    return H @ (H.T @ W)


def scalar_projection(H: np.ndarray, w: np.ndarray) -> np.ndarray:
    if H.size == 0 or H.ndim != 2 or H.shape[1] == 0:
        return np.zeros_like(w)
    return H @ (H.T @ w)


def channel_metrics(model: core.DynamicProvenanceGrowth, faces: List[Face], H: np.ndarray, W: np.ndarray, prefix: str) -> dict:
    total = float(np.linalg.norm(W)) + EPS
    WH = vector_field_projection(H, W)
    htotal = float(np.linalg.norm(WH))
    coh, axis_coh, support = p56.support_coherence(W)
    H_coh, H_axis_coh, H_support = p56.support_coherence(WH)
    raw_k = p56.kappa_ratios(model, faces, W)
    H_k = p56.kappa_ratios(model, faces, WH)
    return {
        f'{prefix}_total_norm': total - EPS,
        f'{prefix}_harmonic_norm': htotal,
        f'{prefix}_harmonic_ratio': htotal / total,
        f'{prefix}_raw_orientation_coherence': coh,
        f'{prefix}_raw_axis_coherence': axis_coh,
        f'{prefix}_raw_support_count': support,
        f'{prefix}_H_orientation_coherence': H_coh,
        f'{prefix}_H_axis_coherence': H_axis_coh,
        f'{prefix}_H_support_count': H_support,
        **{f'{prefix}_raw_{k}': v for k, v in raw_k.items()},
        **{f'{prefix}_H_{k}': v for k, v in H_k.items()},
    }


def transported_quadrature_fields(
    model: core.DynamicProvenanceGrowth,
    K: core.SimplicialComplex,
    pairing_log: List[dict],
    args: argparse.Namespace,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, List[dict]]:
    faces = K.faces()
    idx = {tuple(f): i for i, f in enumerate(faces)}
    W_Q = np.zeros((len(faces), 3), dtype=float)
    W_P = np.zeros((len(faces), 3), dtype=float)
    W_cross = np.zeros((len(faces), 3), dtype=float)
    scalar_area = np.zeros(len(faces), dtype=float)
    rows: List[dict] = []

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
        na = hk.face_normal(model, fa, 'outward')
        nb = hk.face_normal(model, fb, 'outward')
        R_b_to_a = p56.rotation_from_to(nb, -na)
        kb_to_a_reversed = R_b_to_a @ kb

        # Real quadrature split in the a-face chart.
        # Q is the even / observable-like pairing channel already measured in package 56.
        # P is the odd / complementary channel that package 56 did not isolate.
        Q_a = ka + kb_to_a_reversed
        P_a = ka - kb_to_a_reversed

        # Put fields back on the actual face set with the cochain sign convention.
        # For Q, the orientation-reversing gluing gives an anti-oriented copy on b.
        # For P, the pair-odd exchange sign cancels one gluing sign, hence +R.T below.
        Q_b = -(R_b_to_a.T @ Q_a)
        P_b = +(R_b_to_a.T @ P_a)

        cross_a = np.cross(Q_a, P_a)
        cross_b = -(R_b_to_a.T @ cross_a)
        area_strength = float(np.linalg.norm(cross_a))

        ia, ib = idx[fa], idx[fb]
        W_Q[ia] += Q_a
        W_Q[ib] += Q_b
        W_P[ia] += P_a
        W_P[ib] += P_b
        W_cross[ia] += cross_a
        W_cross[ib] += cross_b
        scalar_area[ia] += area_strength
        scalar_area[ib] += area_strength

        qn = float(np.linalg.norm(Q_a))
        pn = float(np.linalg.norm(P_a))
        kan = float(np.linalg.norm(ka))
        kbn = float(np.linalg.norm(kb_to_a_reversed))
        qp_cos = float(np.dot(Q_a, P_a) / ((qn * pn) + EPS))
        transport_cos = float(np.dot(ka, kb_to_a_reversed) / ((kan * kbn) + EPS))
        area_ratio = float(area_strength / ((qn * pn) + EPS))
        comm_axial = p56.axial(skew_from_axial(Q_a) @ skew_from_axial(P_a) - skew_from_axial(P_a) @ skew_from_axial(Q_a))
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
            'ka_norm': kan,
            'transported_reversed_kb_norm': kbn,
            'transport_cosine_ka_kb_reversed': transport_cos,
            'Q_even_norm': qn,
            'P_odd_norm': pn,
            'P_over_Q_norm_ratio': pn / (qn + EPS),
            'Q_P_cosine': qp_cos,
            'symplectic_area_candidate': area_strength,
            'symplectic_area_ratio': area_ratio,
            'cross_commutator_axial_norm': float(np.linalg.norm(comm_axial)),
        })

    return W_Q, W_P, W_cross, scalar_area, rows


def quadrature_metrics(
    model: core.DynamicProvenanceGrowth,
    K: core.SimplicialComplex,
    pairing_log: List[dict],
    args: argparse.Namespace,
) -> Tuple[dict, List[dict], List[dict]]:
    faces = K.faces()
    topo = core.topology(K)
    H, eigs = hk.harmonic_basis_faces(K)

    W_Q, W_P, W_cross, scalar_area, rows = transported_quadrature_fields(model, K, pairing_log, args)
    Qm = channel_metrics(model, faces, H, W_Q, 'Q_even')
    Pm = channel_metrics(model, faces, H, W_P, 'P_odd')
    Cm = channel_metrics(model, faces, H, W_cross, 'QP_cross')

    QH = vector_field_projection(H, W_Q)
    PH = vector_field_projection(H, W_P)
    area_H = scalar_projection(H, scalar_area)
    Q_norm0 = Qm['Q_even_total_norm']
    P_norm0 = Pm['P_odd_total_norm']
    Q_total = Q_norm0 + EPS
    P_total = P_norm0 + EPS
    QH_total = Qm['Q_even_harmonic_norm']
    PH_total = Pm['P_odd_harmonic_norm']

    dot_QP = float(np.sum(W_Q * W_P))
    dot_QH_PH = float(np.sum(QH * PH))
    pair_area_total = float(np.sum([float(r['symplectic_area_candidate']) for r in rows]))
    pair_area_ratio_mean = float(np.mean([float(r['symplectic_area_ratio']) for r in rows])) if rows else 0.0
    mean_transport_cos = float(np.mean([float(r['transport_cosine_ka_kb_reversed']) for r in rows])) if rows else 0.0
    mean_qp_cos = float(np.mean([float(r['Q_P_cosine']) for r in rows])) if rows else 0.0
    mean_p_over_q = float(np.mean([float(r['P_over_Q_norm_ratio']) for r in rows])) if rows else 0.0

    # A purely real diagnostic: the complementary P-channel is significant if it is
    # nonzero and at least comparable to Q. This is not a positivity/norm axiom; it
    # is only a finite-dimensional numerical size comparison of two derived channels.
    metrics = {
        'beta0': topo['beta0'], 'beta1': topo['beta1'], 'beta2': topo['beta2'], 'beta3': topo['beta3'],
        'harmonic_dim_real': int(H.shape[1]) if H.ndim == 2 else 0,
        'applied_pair_count': len(rows),
        'decision_used_delta_beta_any': any(str(r.get('decision_used_delta_beta', '')).lower() == 'true' for r in rows),
        'measured_delta_beta2_sum': sum(int(float(r.get('measured_delta_beta2', 0) or 0)) for r in rows),
        'mean_transport_cosine_ka_kb_reversed': mean_transport_cos,
        'mean_Q_P_cosine': mean_qp_cos,
        'mean_P_over_Q_norm_ratio': mean_p_over_q,
        'Q_P_global_cosine': dot_QP / (Q_total * P_total),
        'QH_PH_global_cosine': dot_QH_PH / ((QH_total * PH_total) + EPS),
        'symplectic_area_candidate_total': pair_area_total,
        'symplectic_area_candidate_mean_ratio': pair_area_ratio_mean,
        'scalar_area_total_norm': float(np.linalg.norm(scalar_area)),
        'scalar_area_harmonic_norm': float(np.linalg.norm(area_H)),
        'scalar_area_harmonic_ratio': float(np.linalg.norm(area_H)) / (float(np.linalg.norm(scalar_area)) + EPS),
        'P_not_killed_relative_to_Q': (P_norm0 / (Q_norm0 + EPS)) if rows else 0.0,
        'P_H_not_killed_relative_to_Q_H': (PH_total / (QH_total + EPS)) if rows else 0.0,
    }
    metrics.update(Qm)
    metrics.update(Pm)
    metrics.update(Cm)

    face_rows: List[dict] = []
    qn = np.linalg.norm(W_Q, axis=1) if len(faces) else np.array([])
    pn = np.linalg.norm(W_P, axis=1) if len(faces) else np.array([])
    qhn = np.linalg.norm(QH, axis=1) if len(faces) else np.array([])
    phn = np.linalg.norm(PH, axis=1) if len(faces) else np.array([])
    cn = np.linalg.norm(W_cross, axis=1) if len(faces) else np.array([])
    for i, f in enumerate(faces):
        if max(float(qn[i]), float(pn[i]), float(qhn[i]), float(phn[i]), float(cn[i])) <= 0:
            continue
        face_rows.append({
            'face': str(list(f)),
            'birth_orders': str([model.nodes[v].birth_order for v in f]),
            'birth_times': str([model.nodes[v].birth_time for v in f]),
            'Q_even_norm': float(qn[i]),
            'P_odd_norm': float(pn[i]),
            'Q_even_H_norm': float(qhn[i]),
            'P_odd_H_norm': float(phn[i]),
            'QP_cross_norm': float(cn[i]),
            'scalar_area': float(scalar_area[i]) if len(scalar_area) else 0.0,
        })
    face_rows.sort(key=lambda r: (r['P_odd_H_norm'], r['Q_even_H_norm'], r['QP_cross_norm']), reverse=True)
    return metrics, rows, face_rows


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
    metrics, pair_rows, face_rows = quadrature_metrics(model, K, pairing_log, args)
    write_csv(vout / 'birth_geometry_log.csv', birth_log)
    write_csv(vout / 'nonlinear_pairing_cascade_log.csv', pairing_log)
    write_csv(vout / 'quadrature_pair_rows.csv', pair_rows)
    write_csv(vout / 'quadrature_face_support_top.csv', face_rows[:args.keep_top_faces])
    summary = {
        'variant': variant,
        'max_level': args.max_level,
        'source': args.source,
        'antisym_eta': args.antisym_eta,
        'phase_sign': args.phase_sign,
        'erase_phase_for_strict_sym': args.erase_phase_for_strict_sym,
        'baseline_metrics_legacy_core': baseline,
        'auto_metrics_legacy_core': auto,
        'quadrature_metrics': metrics,
        'automatic_pairings_applied': sum(1 for x in pairing_log if x.get('applied')),
        'automatic_pairing_attempts_logged': len(pairing_log),
        'births_with_cascade_logs': scans,
        'interpretation_flags': {
            'beta2_opened': auto['beta2'] > baseline['beta2'],
            'Q_even_harmonic_positive': metrics['Q_even_harmonic_ratio'] > args.harmonic_positive_threshold,
            'P_odd_harmonic_positive': metrics['P_odd_harmonic_ratio'] > args.harmonic_positive_threshold,
            'P_channel_comparable_to_Q': metrics['P_not_killed_relative_to_Q'] > args.comparable_threshold,
            'P_H_comparable_to_Q_H': metrics['P_H_not_killed_relative_to_Q_H'] > args.comparable_threshold,
            'symplectic_area_positive': metrics['symplectic_area_candidate_mean_ratio'] > args.symplectic_area_threshold,
            'decision_used_delta_beta_any': metrics['decision_used_delta_beta_any'],
        },
    }
    (vout / 'variant_quadrature_summary.json').write_text(json.dumps(summary, indent=2), encoding='utf-8')
    return summary


def write_comparative(out: Path, rows: List[dict]) -> None:
    flat = []
    for r in rows:
        a = r['auto_metrics_legacy_core']
        m = r['quadrature_metrics']
        flat.append({
            'variant': r['variant'],
            'beta0': a['beta0'], 'beta1': a['beta1'], 'beta2': a['beta2'], 'beta3': a['beta3'],
            'pairings': r['automatic_pairings_applied'],
            'Q_even_total_norm': m['Q_even_total_norm'],
            'P_odd_total_norm': m['P_odd_total_norm'],
            'P_over_Q': m['P_not_killed_relative_to_Q'],
            'Q_even_harmonic_ratio': m['Q_even_harmonic_ratio'],
            'P_odd_harmonic_ratio': m['P_odd_harmonic_ratio'],
            'P_H_over_Q_H': m['P_H_not_killed_relative_to_Q_H'],
            'Q_even_H_kappa_orientation_ratio': m['Q_even_H_kappa_orientation_ratio'],
            'P_odd_H_kappa_orientation_ratio': m['P_odd_H_kappa_orientation_ratio'],
            'Q_even_H_orientation_coherence': m['Q_even_H_orientation_coherence'],
            'P_odd_H_orientation_coherence': m['P_odd_H_orientation_coherence'],
            'QP_cross_total_norm': m['QP_cross_total_norm'],
            'QP_cross_harmonic_ratio': m['QP_cross_harmonic_ratio'],
            'symplectic_area_candidate_mean_ratio': m['symplectic_area_candidate_mean_ratio'],
            'scalar_area_harmonic_ratio': m['scalar_area_harmonic_ratio'],
            'mean_transport_cosine_ka_kb_reversed': m['mean_transport_cosine_ka_kb_reversed'],
            'mean_Q_P_cosine': m['mean_Q_P_cosine'],
            'decision_used_delta_beta_any': m['decision_used_delta_beta_any'],
        })
    write_csv(out / 'comparative_pairing_quadrature_split_summary.csv', flat)


def make_docs(summary: dict) -> tuple[str, str, str, str]:
    rows = summary['variant_rows']
    table_lines = [
        '| variant | beta | pairs | Q harm | P harm | P/Q | P_H/Q_H | sympl area | Q kappa | P kappa | used dBeta? |',
        '|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|',
    ]
    for r in rows:
        a = r['auto_metrics_legacy_core']; m = r['quadrature_metrics']
        table_lines.append(
            f"| {r['variant']} | ({a['beta0']},{a['beta1']},{a['beta2']},{a['beta3']}) | "
            f"{r['automatic_pairings_applied']} | {m['Q_even_harmonic_ratio']:.6g} | "
            f"{m['P_odd_harmonic_ratio']:.6g} | {m['P_not_killed_relative_to_Q']:.6g} | "
            f"{m['P_H_not_killed_relative_to_Q_H']:.6g} | "
            f"{m['symplectic_area_candidate_mean_ratio']:.6g} | "
            f"{m['Q_even_H_kappa_orientation_ratio']:.6g} | "
            f"{m['P_odd_H_kappa_orientation_ratio']:.6g} | "
            f"{m['decision_used_delta_beta_any']} |"
        )
    table = '\n'.join(table_lines)
    smd = f"""# SUMMARY — pairing quadrature split symplectic defect gate

Model label:
CNNA growing primal simplicial complex with deterministic sequential provenance growth,
nonlinear asymmetry-gated complement pairing, and a directed antisymmetric birth-transport
operator.  This is not SG/ST as a global geometry, not a finished NGF/CQNM model,
and not a complex/J/*/positivity derivation.

Purpose:

```text
Test whether the apparent cancellation in the paired channel is only the even
real observable-like projection Q = ka + transport(kb), while the odd channel
P = ka - transport(kb) carries a second real quadrature.
```

No `i`, no `J`, no Hodge star, no positive-frequency split, no imported adjoint,
no positivity, and no final sym(M) is used.  Norms are finite-dimensional diagnostics
only, not axioms.

{table}

Conservative reading:
If P is nonzero and harmonic while strict_sym remains zero, the pair cancellation is
not just a failure; it exposes a candidate real two-quadrature split.  This does not
derive a complex structure.  It only moves the question from spatial orientation-lock
to a possible real operator/quadrature split.
"""
    rmd = f"""# RESULTS — pairing quadrature split symplectic defect gate

## Comparative table

{table}

## Interpretation protocol

```text
Q_even = ka + R(kb)
P_odd  = ka - R(kb)
```

where R transports face b into the orientation-reversed gluing chart of face a.

The decisive question is whether P survives when Q is antikohärent / cancellation-prone.
A positive P channel is not a proof of creation-annihilation structure.  It is only a
derived real pre-structure which could later be tested for operator composition,
anti-automorphism, CCR-like or symplectic closure.

## Anti-smuggling checks

- `decision_used_delta_beta_any` must remain false.
- `strict_symmetrized_control` must not produce the same Q/P signal.
- The antisymmetric transport term is derived from birth_order via q and h.
- No final symmetrization is applied to the directed vertex operator.
"""
    audit = """# SOURCE AUDIT

Inherited result chain:

- Package 51: pair-transport fills H2 but is weakly oriented.
- Package 55: antisymmetric birth transport creates local orientation but does not
  project as raw face-K into H2.
- Package 56: antisymmetric birth transport plus pair transport still has weak
  oriented/kappa lock.
- Package 57: Z3 interfan propagation reduces residuals but does not orient H2.

New hypothesis:
The cancellation may be the even observable-like channel Q, while the missing
complementary information sits in the odd channel P.  This is motivated by the
real-quadrature pattern behind creation/annihilation splitting, but the test itself
does not import a, a†, i, J, *, positivity, frequency, or Hilbert space structure.
"""
    readme = """# Pairing quadrature split symplectic defect gate

Run:

```bash
python3 test_pairing_quadrature_split_symplectic_defect_gate.py
```

Outputs:

- comparative_summary.json
- comparative_pairing_quadrature_split_summary.csv
- RESULTS.md
- SUMMARY.md
- SOURCE_AUDIT.md
- per-variant pair and face support logs

Default L2 run uses variants:

```text
real_growth
strict_symmetrized_control
no_backreaction
```

Next test suggested by this package:
`test_real_quadrature_operator_closure_gate.py`, but only if P survives with a
nontrivial harmonic component and strict_sym remains killed.
"""
    return smd, rmd, audit, readme


def package(out: Path, zip_path: Path) -> None:
    files = [
        Path(__file__).name,
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
    ap.add_argument('--comparable-threshold', type=float, default=0.25)
    ap.add_argument('--symplectic-area-threshold', type=float, default=1e-4)
    ap.add_argument('--kappa-orientation-threshold', type=float, default=0.15)
    ap.add_argument('--antisym-eta', type=float, default=1.0)
    ap.add_argument('--phase-sign', type=int, default=1, choices=[-1, 1])
    ap.add_argument('--erase-phase-for-strict-sym', action='store_true', default=True)
    ap.add_argument('--variants', nargs='*', default=['real_growth', 'strict_symmetrized_control', 'no_backreaction'])
    ap.add_argument('--out', default='pairing_quadrature_split_symplectic_defect_out_L2')
    ap.add_argument('--zip', default='cnna_pairing_quadrature_split_symplectic_defect_gate_pkg_L2.zip')
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
                'Q_even_harmonic_ratio': r['quadrature_metrics']['Q_even_harmonic_ratio'],
                'P_odd_harmonic_ratio': r['quadrature_metrics']['P_odd_harmonic_ratio'],
                'P_over_Q': r['quadrature_metrics']['P_not_killed_relative_to_Q'],
                'P_H_over_Q_H': r['quadrature_metrics']['P_H_not_killed_relative_to_Q_H'],
                'symplectic_area_candidate_mean_ratio': r['quadrature_metrics']['symplectic_area_candidate_mean_ratio'],
                'decision_used_delta_beta_any': r['quadrature_metrics']['decision_used_delta_beta_any'],
            } for r in rows
        ]
    }, indent=2))


if __name__ == '__main__':
    main()
