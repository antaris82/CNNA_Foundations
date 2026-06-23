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


def signed_ratio(values: List[float], abs_values: List[float]) -> float:
    denom = float(np.sum(abs_values)) + EPS
    return float(np.sum(values)) / denom if values else 0.0


def mean_abs(values: List[float]) -> float:
    return float(np.mean(np.abs(values))) if values else 0.0


def mean_signed(values: List[float]) -> float:
    return float(np.mean(values)) if values else 0.0


def signed_quadrature_rows(
    model: core.DynamicProvenanceGrowth,
    K: core.SimplicialComplex,
    pairing_log: List[dict],
    args: argparse.Namespace,
) -> Tuple[dict, List[dict], List[dict]]:
    faces = K.faces()
    idx = {tuple(f): i for i, f in enumerate(faces)}
    topo = core.topology(K)
    H, eigs = hk.harmonic_basis_faces(K)

    W_Q = np.zeros((len(faces), 3), dtype=float)
    W_P = np.zeros((len(faces), 3), dtype=float)
    W_cross = np.zeros((len(faces), 3), dtype=float)
    scalar_abs_area = np.zeros(len(faces), dtype=float)
    scalar_signed_birth = np.zeros(len(faces), dtype=float)
    scalar_signed_out = np.zeros(len(faces), dtype=float)

    pair_rows: List[dict] = []
    signed_birth_vals: List[float] = []
    signed_out_vals: List[float] = []
    signed_pair_axis_vals: List[float] = []
    abs_vals: List[float] = []

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
        nb_birth = hk.face_normal(model, fb, 'birth_order')
        R_b_to_a = p56.rotation_from_to(nb_out, -na_out)
        kb_to_a_reversed = R_b_to_a @ kb

        Q_a = ka + kb_to_a_reversed
        P_a = ka - kb_to_a_reversed
        Q_b = -(R_b_to_a.T @ Q_a)
        P_b = +(R_b_to_a.T @ P_a)
        cross_a = np.cross(Q_a, P_a)
        cross_b = -(R_b_to_a.T @ cross_a)

        abs_area = float(np.linalg.norm(cross_a))
        birth_signed = float(np.dot(cross_a, na_birth))
        out_signed = float(np.dot(cross_a, na_out))

        # Pair-axis reference is derived from the actual gluing geometry: the axis
        # from face-a centroid to face-b centroid, expressed in the a-chart.  It is
        # not a global orientation choice; it is logged as an auxiliary diagnostic.
        ca = sum((model.nodes[v].pos for v in fa), np.zeros(3)) / 3.0
        cb = sum((model.nodes[v].pos for v in fb), np.zeros(3)) / 3.0
        pair_axis = cb - ca
        pan = float(np.linalg.norm(pair_axis))
        pair_axis = pair_axis / (pan + EPS) if pan > EPS else np.zeros(3)
        pair_axis_signed = float(np.dot(cross_a, pair_axis))

        ia, ib = idx[fa], idx[fb]
        W_Q[ia] += Q_a
        W_Q[ib] += Q_b
        W_P[ia] += P_a
        W_P[ib] += P_b
        W_cross[ia] += cross_a
        W_cross[ib] += cross_b
        scalar_abs_area[ia] += abs_area
        scalar_abs_area[ib] += abs_area
        scalar_signed_birth[ia] += birth_signed
        scalar_signed_birth[ib] -= birth_signed
        scalar_signed_out[ia] += out_signed
        scalar_signed_out[ib] -= out_signed

        qn = float(np.linalg.norm(Q_a))
        pn = float(np.linalg.norm(P_a))
        kan = float(np.linalg.norm(ka))
        kbn = float(np.linalg.norm(kb_to_a_reversed))
        transport_cos = float(np.dot(ka, kb_to_a_reversed) / ((kan * kbn) + EPS))
        qp_cos = float(np.dot(Q_a, P_a) / ((qn * pn) + EPS))
        area_ratio = abs_area / ((qn * pn) + EPS)
        signed_birth_ratio_local = birth_signed / (abs_area + EPS)
        signed_out_ratio_local = out_signed / (abs_area + EPS)
        signed_pair_axis_ratio_local = pair_axis_signed / (abs_area + EPS)
        signed_birth_vals.append(birth_signed)
        signed_out_vals.append(out_signed)
        signed_pair_axis_vals.append(pair_axis_signed)
        abs_vals.append(abs_area)
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
            'ka_norm': kan,
            'transported_reversed_kb_norm': kbn,
            'transport_cosine_ka_kb_reversed': transport_cos,
            'Q_even_norm': qn,
            'P_odd_norm': pn,
            'P_over_Q_norm_ratio': pn / (qn + EPS),
            'Q_P_cosine': qp_cos,
            'abs_area_norm_Q_cross_P': abs_area,
            'abs_area_ratio': area_ratio,
            'signed_birth_area': birth_signed,
            'signed_outward_area': out_signed,
            'signed_pair_axis_area': pair_axis_signed,
            'signed_birth_over_abs': signed_birth_ratio_local,
            'signed_outward_over_abs': signed_out_ratio_local,
            'signed_pair_axis_over_abs': signed_pair_axis_ratio_local,
        })

    Qm = p58.channel_metrics(model, faces, H, W_Q, 'Q_even')
    Pm = p58.channel_metrics(model, faces, H, W_P, 'P_odd')
    Cm = p58.channel_metrics(model, faces, H, W_cross, 'QP_cross')
    abs_area_H = p58.scalar_projection(H, scalar_abs_area)
    signed_birth_H = p58.scalar_projection(H, scalar_signed_birth)
    signed_out_H = p58.scalar_projection(H, scalar_signed_out)

    abs_total = float(np.sum(abs_vals)) + EPS
    signed_birth_sum = float(np.sum(signed_birth_vals)) if signed_birth_vals else 0.0
    signed_out_sum = float(np.sum(signed_out_vals)) if signed_out_vals else 0.0
    signed_pair_axis_sum = float(np.sum(signed_pair_axis_vals)) if signed_pair_axis_vals else 0.0

    metrics = {
        'beta0': topo['beta0'], 'beta1': topo['beta1'], 'beta2': topo['beta2'], 'beta3': topo['beta3'],
        'harmonic_dim_real': int(H.shape[1]) if H.ndim == 2 else 0,
        'applied_pair_count': len(pair_rows),
        'decision_used_delta_beta_any': any(str(r.get('decision_used_delta_beta', '')).lower() == 'true' for r in pair_rows),
        'measured_delta_beta2_sum': sum(int(float(r.get('measured_delta_beta2', 0) or 0)) for r in pair_rows),
        'mean_transport_cosine_ka_kb_reversed': float(np.mean([r['transport_cosine_ka_kb_reversed'] for r in pair_rows])) if pair_rows else 0.0,
        'mean_Q_P_cosine': float(np.mean([r['Q_P_cosine'] for r in pair_rows])) if pair_rows else 0.0,
        'mean_P_over_Q_norm_ratio': float(np.mean([r['P_over_Q_norm_ratio'] for r in pair_rows])) if pair_rows else 0.0,
        'abs_area_total': abs_total - EPS,
        'abs_area_mean_ratio': float(np.mean([r['abs_area_ratio'] for r in pair_rows])) if pair_rows else 0.0,
        'signed_birth_area_sum': signed_birth_sum,
        'signed_outward_area_sum': signed_out_sum,
        'signed_pair_axis_area_sum': signed_pair_axis_sum,
        'signed_birth_over_abs_sum_ratio': signed_birth_sum / abs_total,
        'signed_outward_over_abs_sum_ratio': signed_out_sum / abs_total,
        'signed_pair_axis_over_abs_sum_ratio': signed_pair_axis_sum / abs_total,
        'mean_signed_birth_over_abs': float(np.mean([r['signed_birth_over_abs'] for r in pair_rows])) if pair_rows else 0.0,
        'mean_signed_outward_over_abs': float(np.mean([r['signed_outward_over_abs'] for r in pair_rows])) if pair_rows else 0.0,
        'mean_signed_pair_axis_over_abs': float(np.mean([r['signed_pair_axis_over_abs'] for r in pair_rows])) if pair_rows else 0.0,
        'mean_abs_signed_birth_over_abs': float(np.mean([abs(r['signed_birth_over_abs']) for r in pair_rows])) if pair_rows else 0.0,
        'mean_abs_signed_outward_over_abs': float(np.mean([abs(r['signed_outward_over_abs']) for r in pair_rows])) if pair_rows else 0.0,
        'mean_abs_signed_pair_axis_over_abs': float(np.mean([abs(r['signed_pair_axis_over_abs']) for r in pair_rows])) if pair_rows else 0.0,
        'scalar_abs_area_harmonic_ratio': float(np.linalg.norm(abs_area_H)) / (float(np.linalg.norm(scalar_abs_area)) + EPS),
        'scalar_signed_birth_harmonic_ratio': float(np.linalg.norm(signed_birth_H)) / (float(np.linalg.norm(scalar_signed_birth)) + EPS),
        'scalar_signed_outward_harmonic_ratio': float(np.linalg.norm(signed_out_H)) / (float(np.linalg.norm(scalar_signed_out)) + EPS),
        'scalar_signed_birth_total_norm': float(np.linalg.norm(scalar_signed_birth)),
        'scalar_signed_outward_total_norm': float(np.linalg.norm(scalar_signed_out)),
    }
    metrics.update(Qm)
    metrics.update(Pm)
    metrics.update(Cm)

    face_rows: List[dict] = []
    qn = np.linalg.norm(W_Q, axis=1) if len(faces) else np.array([])
    pn = np.linalg.norm(W_P, axis=1) if len(faces) else np.array([])
    cn = np.linalg.norm(W_cross, axis=1) if len(faces) else np.array([])
    for i, f in enumerate(faces):
        if max(float(qn[i]), float(pn[i]), float(cn[i]), abs(float(scalar_signed_birth[i]))) <= 0:
            continue
        face_rows.append({
            'face': str(list(f)),
            'birth_orders': str([model.nodes[v].birth_order for v in f]),
            'birth_times': str([model.nodes[v].birth_time for v in f]),
            'Q_even_norm': float(qn[i]),
            'P_odd_norm': float(pn[i]),
            'QP_cross_norm': float(cn[i]),
            'abs_area_scalar': float(scalar_abs_area[i]),
            'signed_birth_area_scalar': float(scalar_signed_birth[i]),
            'signed_outward_area_scalar': float(scalar_signed_out[i]),
        })
    face_rows.sort(key=lambda r: (abs(r['signed_birth_area_scalar']), r['QP_cross_norm']), reverse=True)
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
    metrics, pair_rows, face_rows = signed_quadrature_rows(model, K, pairing_log, local_args)
    write_csv(vout / 'birth_geometry_log.csv', birth_log)
    write_csv(vout / 'nonlinear_pairing_cascade_log.csv', pairing_log)
    write_csv(vout / 'signed_quadrature_pair_rows.csv', pair_rows)
    write_csv(vout / 'signed_quadrature_face_support_top.csv', face_rows[:args.keep_top_faces])
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
        'signed_quadrature_metrics': metrics,
        'automatic_pairings_applied': sum(1 for x in pairing_log if x.get('applied')),
        'automatic_pairing_attempts_logged': len(pairing_log),
        'births_with_cascade_logs': scans,
        'interpretation_flags': {
            'beta2_opened': auto['beta2'] > baseline['beta2'],
            'Q_even_harmonic_positive': metrics['Q_even_harmonic_ratio'] > args.harmonic_positive_threshold,
            'P_odd_harmonic_positive': metrics['P_odd_harmonic_ratio'] > args.harmonic_positive_threshold,
            'abs_area_positive': metrics['abs_area_mean_ratio'] > args.area_positive_threshold,
            'signed_birth_area_nontrivial': abs(metrics['signed_birth_over_abs_sum_ratio']) > args.signed_area_threshold,
            'signed_outward_area_nontrivial': abs(metrics['signed_outward_over_abs_sum_ratio']) > args.signed_area_threshold,
            'decision_used_delta_beta_any': metrics['decision_used_delta_beta_any'],
        },
    }
    (vout / 'variant_signed_quadrature_summary.json').write_text(json.dumps(summary, indent=2), encoding='utf-8')
    return summary


def flip_comparisons(rows: List[dict]) -> List[dict]:
    by_variant = {}
    for r in rows:
        by_variant.setdefault(r['variant'], {})[r['phase_sign']] = r
    comps = []
    for variant, d in sorted(by_variant.items()):
        if 1 not in d or -1 not in d:
            continue
        p = d[1]['signed_quadrature_metrics']
        m = d[-1]['signed_quadrature_metrics']
        def flip_score(key: str) -> float:
            denom = abs(float(p[key])) + abs(float(m[key])) + EPS
            return (float(p[key]) + float(m[key])) / denom
        comps.append({
            'variant': variant,
            'signed_birth_plus': p['signed_birth_over_abs_sum_ratio'],
            'signed_birth_minus': m['signed_birth_over_abs_sum_ratio'],
            'signed_birth_flip_score_zero_if_perfect_flip': flip_score('signed_birth_over_abs_sum_ratio'),
            'signed_outward_plus': p['signed_outward_over_abs_sum_ratio'],
            'signed_outward_minus': m['signed_outward_over_abs_sum_ratio'],
            'signed_outward_flip_score_zero_if_perfect_flip': flip_score('signed_outward_over_abs_sum_ratio'),
            'signed_pair_axis_plus': p['signed_pair_axis_over_abs_sum_ratio'],
            'signed_pair_axis_minus': m['signed_pair_axis_over_abs_sum_ratio'],
            'signed_pair_axis_flip_score_zero_if_perfect_flip': flip_score('signed_pair_axis_over_abs_sum_ratio'),
            'abs_area_plus': p['abs_area_mean_ratio'],
            'abs_area_minus': m['abs_area_mean_ratio'],
            'Q_harm_plus': p['Q_even_harmonic_ratio'],
            'Q_harm_minus': m['Q_even_harmonic_ratio'],
            'P_harm_plus': p['P_odd_harmonic_ratio'],
            'P_harm_minus': m['P_odd_harmonic_ratio'],
        })
    return comps


def write_comparative(out: Path, rows: List[dict], comps: List[dict]) -> None:
    flat = []
    for r in rows:
        a = r['auto_metrics_legacy_core']; m = r['signed_quadrature_metrics']
        flat.append({
            'variant_phase': r['variant_phase'],
            'variant': r['variant'],
            'phase_sign': r['phase_sign'],
            'beta0': a['beta0'], 'beta1': a['beta1'], 'beta2': a['beta2'], 'beta3': a['beta3'],
            'pairings': r['automatic_pairings_applied'],
            'Q_even_harmonic_ratio': m['Q_even_harmonic_ratio'],
            'P_odd_harmonic_ratio': m['P_odd_harmonic_ratio'],
            'Q_even_H_kappa_orientation_ratio': m['Q_even_H_kappa_orientation_ratio'],
            'P_odd_H_kappa_orientation_ratio': m['P_odd_H_kappa_orientation_ratio'],
            'abs_area_mean_ratio': m['abs_area_mean_ratio'],
            'signed_birth_over_abs_sum_ratio': m['signed_birth_over_abs_sum_ratio'],
            'signed_outward_over_abs_sum_ratio': m['signed_outward_over_abs_sum_ratio'],
            'signed_pair_axis_over_abs_sum_ratio': m['signed_pair_axis_over_abs_sum_ratio'],
            'mean_abs_signed_birth_over_abs': m['mean_abs_signed_birth_over_abs'],
            'mean_abs_signed_outward_over_abs': m['mean_abs_signed_outward_over_abs'],
            'mean_abs_signed_pair_axis_over_abs': m['mean_abs_signed_pair_axis_over_abs'],
            'scalar_signed_birth_harmonic_ratio': m['scalar_signed_birth_harmonic_ratio'],
            'scalar_signed_outward_harmonic_ratio': m['scalar_signed_outward_harmonic_ratio'],
            'mean_transport_cosine_ka_kb_reversed': m['mean_transport_cosine_ka_kb_reversed'],
            'decision_used_delta_beta_any': m['decision_used_delta_beta_any'],
        })
    write_csv(out / 'comparative_signed_quadrature_area_summary.csv', flat)
    write_csv(out / 'phase_flip_comparison.csv', comps)


def make_docs(summary: dict, comps: List[dict]) -> tuple[str, str, str, str]:
    rows = summary['variant_rows']
    table_lines = [
        '| variant/phase | beta | pairs | Q harm | P harm | abs area | signed birth/abs | signed outward/abs | Q kappa | P kappa | used dBeta? |',
        '|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|',
    ]
    for r in rows:
        a = r['auto_metrics_legacy_core']; m = r['signed_quadrature_metrics']
        table_lines.append(
            f"| {r['variant_phase']} | ({a['beta0']},{a['beta1']},{a['beta2']},{a['beta3']}) | "
            f"{r['automatic_pairings_applied']} | {m['Q_even_harmonic_ratio']:.6g} | "
            f"{m['P_odd_harmonic_ratio']:.6g} | {m['abs_area_mean_ratio']:.6g} | "
            f"{m['signed_birth_over_abs_sum_ratio']:.6g} | {m['signed_outward_over_abs_sum_ratio']:.6g} | "
            f"{m['Q_even_H_kappa_orientation_ratio']:.6g} | {m['P_odd_H_kappa_orientation_ratio']:.6g} | "
            f"{m['decision_used_delta_beta_any']} |"
        )
    table = '\n'.join(table_lines)
    comp_lines = [
        '| variant | birth plus | birth minus | birth flip-score | outward plus | outward minus | outward flip-score |',
        '|---|---:|---:|---:|---:|---:|---:|',
    ]
    for c in comps:
        comp_lines.append(
            f"| {c['variant']} | {c['signed_birth_plus']:.6g} | {c['signed_birth_minus']:.6g} | "
            f"{c['signed_birth_flip_score_zero_if_perfect_flip']:.6g} | {c['signed_outward_plus']:.6g} | "
            f"{c['signed_outward_minus']:.6g} | {c['signed_outward_flip_score_zero_if_perfect_flip']:.6g} |"
        )
    comp_table = '\n'.join(comp_lines)
    smd = f"""# SUMMARY — signed quadrature area kappa gate

Model label:
CNNA growing primal simplicial complex with deterministic sequential provenance growth,
nonlinear asymmetry-gated complement pairing, and directed antisymmetric birth-transport
operators.  This is not a complex/J/*/positivity derivation.

Purpose:
Package 58 measured `||Q x P||`, a magnitude.  This package keeps that diagnostic
but adds signed tests:

```text
signed_birth  = <Q x P, n_birth>
signed_out    = <Q x P, n_outward>
signed_axis   = <Q x P, pair_axis>
```

It also runs phase_sign = +1 and -1 to test whether signed area flips.  No final
symmetrization is used in the directed vertex operator.

{table}

## Phase-sign flip comparison

{comp_table}

Conservative reading:
A true oriented/symplectic candidate would need a nontrivial signed ratio and a
controlled flip.  A large abs-area alone is only a magnitude.
"""
    rmd = f"""# RESULTS — signed quadrature area kappa gate

## Comparative table

{table}

## Phase-sign flip comparison

{comp_table}

## Interpretation protocol

This package explicitly distinguishes:

```text
abs area     = ||Q x P||
signed area  = <Q x P, n_ref>
```

The abs area is not allowed to count as oriented symplectic structure.  The gate is
only positive if the signed ratios are nontrivial and transform coherently under the
phase/kappa proxy while strict_sym remains killed.

## Anti-smuggling conditions

- no `i`, no `J`, no Hodge star, no imported adjoint, no positivity;
- no final `sym(M)` in the antisymmetric birth-transport operator;
- `decision_used_delta_beta_any` must remain false;
- abs-area is treated as magnitude only.
"""
    audit = """# SOURCE AUDIT

Package 58 found a robust Q/P split, but its `symplectic_area_candidate` was
`np.linalg.norm(Q x P)`.  This package tests Claude's objection directly by replacing
that magnitude-only quantity with signed projections against derived local references:
face birth-order normal, outward normal, and pair-axis.

A positive abs-area with near-zero or non-flipping signed area is not a J/i/symplectic
orientation derivation.  It is only a nonzero two-quadrature magnitude.
"""
    readme = """# Signed quadrature area kappa gate

Run:

```bash
python3 test_signed_quadrature_area_kappa_gate.py
```

Outputs:

- comparative_summary.json
- comparative_signed_quadrature_area_summary.csv
- phase_flip_comparison.csv
- RESULTS.md
- SUMMARY.md
- SOURCE_AUDIT.md
- per-variant signed pair/face logs

Next gate depends on this result:

- if signed area is small/non-flipping: treat Q/P as magnitude-only and search for a
  native operator involution/quadrature pairing instead of spatial orientation;
- if signed area is nontrivial and flips: build a strict derived real symplectic-form
  closure test.
"""
    return smd, rmd, audit, readme


def package(out: Path, zip_path: Path) -> None:
    files = [
        Path(__file__).name,
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
    ap.add_argument('--area-positive-threshold', type=float, default=1e-4)
    ap.add_argument('--signed-area-threshold', type=float, default=0.15)
    ap.add_argument('--kappa-orientation-threshold', type=float, default=0.15)
    ap.add_argument('--antisym-eta', type=float, default=1.0)
    ap.add_argument('--erase-phase-for-strict-sym', action='store_true', default=True)
    ap.add_argument('--phase-signs', nargs='*', type=int, default=[1, -1], choices=[-1, 1])
    ap.add_argument('--variants', nargs='*', default=['real_growth', 'strict_symmetrized_control', 'no_backreaction'])
    ap.add_argument('--out', default='signed_quadrature_area_kappa_out_L2')
    ap.add_argument('--zip', default='cnna_signed_quadrature_area_kappa_gate_pkg_L2.zip')
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
                'Q_harm': r['signed_quadrature_metrics']['Q_even_harmonic_ratio'],
                'P_harm': r['signed_quadrature_metrics']['P_odd_harmonic_ratio'],
                'abs_area': r['signed_quadrature_metrics']['abs_area_mean_ratio'],
                'signed_birth_over_abs': r['signed_quadrature_metrics']['signed_birth_over_abs_sum_ratio'],
                'signed_outward_over_abs': r['signed_quadrature_metrics']['signed_outward_over_abs_sum_ratio'],
                'decision_used_delta_beta_any': r['signed_quadrature_metrics']['decision_used_delta_beta_any'],
            } for r in rows
        ],
        'phase_flip_comparisons': comps,
    }, indent=2))


if __name__ == '__main__':
    main()
