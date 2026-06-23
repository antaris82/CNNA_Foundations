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
        # choose a deterministic orthogonal axis
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


def transported_pair_fields(model: core.DynamicProvenanceGrowth, K: core.SimplicialComplex, pairing_log: List[dict], source: str) -> Tuple[np.ndarray, np.ndarray, List[dict]]:
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
        ka = hk.skew_to_axial(core.face_K(model, fa, source))
        kb = hk.skew_to_axial(core.face_K(model, fb, source))
        na = hk.face_normal(model, fa, 'outward')
        nb = hk.face_normal(model, fb, 'outward')
        # Pairing is orientation reversing on the face.  Compare b to a by first
        # rotating b's outward normal into -a's outward normal; the cochain sign
        # reversal is then represented by adding the transported b-axial vector.
        R_b_to_a = rotation_from_to(nb, -na)
        kb_to_a_reversed = R_b_to_a @ kb
        pair_vec_a = ka + kb_to_a_reversed
        # Put an anti-oriented transported copy back on b so the field lives on
        # the actual final face set, not only in a local chart of face a.
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


def pair_transport_metrics(model: core.DynamicProvenanceGrowth, K: core.SimplicialComplex, pairing_log: List[dict], source: str) -> Tuple[dict, List[dict], List[dict]]:
    faces = K.faces()
    topo = core.topology(K)
    H, eigs = hk.harmonic_basis_faces(K)
    W_pair, scalar_pair, pair_rows = transported_pair_fields(model, K, pairing_log, source)
    Wh = vector_field_projection(H, W_pair)
    sh = scalar_projection(H, scalar_pair)
    total = float(np.linalg.norm(W_pair)) + EPS
    htotal = float(np.linalg.norm(Wh))
    stotal = float(np.linalg.norm(scalar_pair)) + EPS
    shtotal = float(np.linalg.norm(sh))
    hn = np.linalg.norm(Wh, axis=1) if len(faces) else np.array([])
    normals_out = np.array([hk.face_normal(model, f, 'outward') for f in faces], dtype=float) if faces else np.zeros((0,3))
    normals_birth = np.array([hk.face_normal(model, f, 'birth_order') for f in faces], dtype=float) if faces else np.zeros((0,3))
    areas = np.array([max(hk.face_area(model, f), EPS) for f in faces], dtype=float) if faces else np.array([])
    denom = float(np.sum(hn * areas)) + EPS
    dot_out = np.einsum('ij,ij->i', Wh, normals_out) if len(faces) else np.array([])
    dot_birth = np.einsum('ij,ij->i', Wh, normals_birth) if len(faces) else np.array([])
    signed_out = float(np.sum(dot_out * areas)) / denom if len(faces) else 0.0
    abs_out = float(np.sum(np.abs(dot_out) * areas)) / denom if len(faces) else 0.0
    signed_birth = float(np.sum(dot_birth * areas)) / denom if len(faces) else 0.0
    abs_birth = float(np.sum(np.abs(dot_birth) * areas)) / denom if len(faces) else 0.0
    mask = hn > max(1e-10, 1e-8 * (float(np.max(hn)) if len(hn) else 1.0))
    if np.any(mask):
        unit_vecs = Wh[mask] / (hn[mask][:, None] + EPS)
        coherence = float(np.linalg.norm(np.mean(unit_vecs, axis=0)))
    else:
        coherence = 0.0
    top_rows: List[dict] = []
    for i, f in enumerate(faces):
        if float(hn[i]) <= 0 and abs(float(sh[i])) <= 0:
            continue
        top_rows.append({
            'face': str(list(f)),
            'birth_times': str([model.nodes[v].birth_time for v in f]),
            'birth_orders': str([model.nodes[v].birth_order for v in f]),
            'pair_vector_norm': float(np.linalg.norm(W_pair[i])),
            'pair_harmonic_axial_norm': float(hn[i]),
            'pair_scalar_strength': float(scalar_pair[i]),
            'pair_scalar_harmonic_value': float(sh[i]),
            'H_dot_outward_normal': float(dot_out[i]) if len(faces) else 0.0,
            'H_dot_birth_normal': float(dot_birth[i]) if len(faces) else 0.0,
            'area': float(areas[i]) if len(faces) else 0.0,
        })
    top_rows.sort(key=lambda r: (r['pair_harmonic_axial_norm'], abs(r['pair_scalar_harmonic_value'])), reverse=True)
    metrics = {
        'beta0': topo['beta0'], 'beta1': topo['beta1'], 'beta2': topo['beta2'], 'beta3': topo['beta3'],
        'harmonic_dim_real': int(H.shape[1]) if H.ndim == 2 else 0,
        'applied_pair_count': len(pair_rows),
        'pair_transport_total_norm': total - EPS,
        'pair_transport_harmonic_norm': htotal,
        'pair_transport_harmonic_ratio': htotal / total,
        'pair_scalar_total_norm': stotal - EPS,
        'pair_scalar_harmonic_norm': shtotal,
        'pair_scalar_harmonic_ratio': shtotal / stotal,
        'pair_kappa_orientation_ratio': abs(signed_out) / (abs_out + EPS),
        'pair_kappa_birth_orientation_ratio': abs(signed_birth) / (abs_birth + EPS),
        'pair_normal_flux_signed_ratio': signed_out,
        'pair_normal_flux_abs_ratio': abs_out,
        'pair_birth_normal_flux_signed_ratio': signed_birth,
        'pair_birth_normal_flux_abs_ratio': abs_birth,
        'pair_orientation_coherence': coherence,
        'pair_support_fraction': float(np.mean(mask)) if len(mask) else 0.0,
        'decision_used_delta_beta_any': any(str(r.get('decision_used_delta_beta', '')).lower() == 'true' for r in pair_rows),
        'measured_delta_beta2_sum': sum(int(float(r.get('measured_delta_beta2', 0) or 0)) for r in pair_rows),
        'mean_pair_transport_cosine': float(np.mean([float(r['a_b_transport_cosine']) for r in pair_rows])) if pair_rows else 0.0,
    }
    return metrics, pair_rows, top_rows


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
    metrics, pair_rows, face_rows = pair_transport_metrics(model, K, pairing_log, args.source)
    write_csv(vout / 'birth_geometry_log.csv', birth_log)
    write_csv(vout / 'nonlinear_pairing_cascade_log.csv', pairing_log)
    write_csv(vout / 'pairing_transport_pairs.csv', pair_rows)
    write_csv(vout / 'pairing_transport_harmonic_faces_top.csv', face_rows[:args.keep_top_faces])
    summary = {
        'variant': variant,
        'max_level': args.max_level,
        'source': args.source,
        'baseline_metrics': baseline,
        'auto_metrics': auto,
        'pairing_transport_metrics': metrics,
        'automatic_pairings_applied': sum(1 for x in pairing_log if x.get('applied')),
        'automatic_pairing_attempts_logged': len(pairing_log),
        'births_with_cascade_logs': scans,
        'interpretation_flags': {
            'beta2_opened': auto['beta2'] > baseline['beta2'],
            'scalar_K_harmonic_positive': auto['harmonic_ratio'] > args.harmonic_positive_threshold,
            'pair_transport_axial_harmonic_positive': metrics['pair_transport_harmonic_ratio'] > args.harmonic_positive_threshold,
            'pair_transport_scalar_harmonic_positive': metrics['pair_scalar_harmonic_ratio'] > args.harmonic_positive_threshold,
            'pair_kappa_signed_bias_positive': metrics['pair_kappa_orientation_ratio'] > args.kappa_orientation_threshold,
            'pair_birth_kappa_signed_bias_positive': metrics['pair_kappa_birth_orientation_ratio'] > args.kappa_orientation_threshold,
            'decision_used_delta_beta_any': metrics['decision_used_delta_beta_any'],
        },
    }
    (vout / 'variant_pairing_transport_summary.json').write_text(json.dumps(summary, indent=2), encoding='utf-8')
    return summary


def write_comparative(out: Path, rows: List[dict]) -> None:
    flat = []
    for r in rows:
        a = r['auto_metrics']; m = r['pairing_transport_metrics']
        flat.append({
            'variant': r['variant'],
            'beta0': a['beta0'], 'beta1': a['beta1'], 'beta2': a['beta2'], 'beta3': a['beta3'],
            'automatic_pairings_applied': r['automatic_pairings_applied'],
            'scalar_K_harmonic_ratio': a['harmonic_ratio'],
            'pair_transport_harmonic_ratio': m['pair_transport_harmonic_ratio'],
            'pair_scalar_harmonic_ratio': m['pair_scalar_harmonic_ratio'],
            'pair_kappa_orientation_ratio': m['pair_kappa_orientation_ratio'],
            'pair_kappa_birth_orientation_ratio': m['pair_kappa_birth_orientation_ratio'],
            'pair_orientation_coherence': m['pair_orientation_coherence'],
            'harmonic_dim_real': m['harmonic_dim_real'],
            'decision_used_delta_beta_any': m['decision_used_delta_beta_any'],
            'measured_delta_beta2_sum': m['measured_delta_beta2_sum'],
            'mean_pair_transport_cosine': m['mean_pair_transport_cosine'],
        })
    write_csv(out / 'comparative_pairing_transport_kappa_summary.csv', flat)


def make_docs(summary: dict) -> Tuple[str, str, str, str]:
    rows = summary['variant_rows']
    table_lines = ['| variant | beta auto | pairings | scalar K harm | pair axial harm | pair scalar harm | Hdim | pair kappa | pair birth kappa | delta beta2 sum | used Δβ? |', '|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|']
    for r in rows:
        a = r['auto_metrics']; m = r['pairing_transport_metrics']
        table_lines.append(f"| {r['variant']} | ({a['beta0']},{a['beta1']},{a['beta2']},{a['beta3']}) | {r['automatic_pairings_applied']} | {a['harmonic_ratio']:.6g} | {m['pair_transport_harmonic_ratio']:.6g} | {m['pair_scalar_harmonic_ratio']:.6g} | {m['harmonic_dim_real']} | {m['pair_kappa_orientation_ratio']:.6g} | {m['pair_kappa_birth_orientation_ratio']:.6g} | {m['measured_delta_beta2_sum']} | {m['decision_used_delta_beta_any']} |")
    table = '\n'.join(table_lines)
    smd = f"""# SUMMARY — pairing-transport harmonic kappa gate

This package tests the next anti-smuggling step after nonlinear asymmetry-gated complement pairing.

Previous result: beta2 opens and the scalar |K| field has a harmonic component.  This package asks whether the **actual applied pairings** carry a transported axial K-flow:

```text
K_pair(face_a, face_b) = K_face_a - orientation_reversed_transport(K_face_b)
```

The pairing logs are used only after the nonlinear growth has already selected and applied moves.  The transport/harmonic/kappa projection is therefore diagnostic and not part of the move decision.

{table}

Read the result conservatively: positive beta2 is a carrier; positive pair-transport axial harmonic ratio is evidence that the pairing operation itself carries an oriented skew-flow into H2.
"""
    rmd = f"""# RESULTS — pairing-transport harmonic kappa gate

## Comparative table

{table}

## Interpretation

The test distinguishes three levels:

```text
1. beta2 carrier opens.
2. scalar |K| has a harmonic residual.
3. transported axial K_pair has a harmonic and kappa-biased component.
```

A strong candidate requires all three in real_growth and a strict kill in strict_symmetrized_control.

Important: `decision_used_delta_beta_any` must remain false.  If true, the topology would have entered the selection rule.  In this package, topology and harmonic projection are measured after the fact.

## Current status

See `comparative_pairing_transport_kappa_summary.csv` and each variant directory for pair-level logs and top harmonic-face support.
"""
    audit = """# SOURCE AUDIT 1–40

Relevant inherited threads:

- Script 1/2: sequential birth environment and older-sibling asymmetry.
- Script 35: real local `K_abc=[A_ab,A_bc]` skew sector.
- Script 40: local tetrahedral boundary closure as an obstruction.
- Nonlinear cascade package: asymmetry-gated complement pairing opens beta2 without using beta in the decision.

This package does not claim a final J derivation.  It tests whether the pairing operation transports the axial skew K-sector into the newly opened H2 carrier.
"""
    readme = f"""# Pairing transport harmonic kappa gate

Run:

```bash
python3 test_pairing_transport_harmonic_kappa_gate.py
```

Example L3:

```bash
python3 test_pairing_transport_harmonic_kappa_gate.py \\
  --max-level 3 \\
  --max-cascade-per-birth 3 \\
  --max-auto-pairings 2 \\
  --out pairing_transport_harmonic_kappa_out_L3 \\
  --zip cnna_pairing_transport_harmonic_kappa_gate_pkg_L3.zip
```

Outputs include JSON, CSV, RESULTS.md, SUMMARY.md and per-variant pair transport logs.
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
    ap.add_argument('--kappa-orientation-threshold', type=float, default=0.15)
    ap.add_argument('--variants', nargs='*', default=['real_growth', 'strict_symmetrized_control', 'no_backreaction'])
    ap.add_argument('--out', default='pairing_transport_harmonic_kappa_out_L2')
    ap.add_argument('--zip', default='cnna_pairing_transport_harmonic_kappa_gate_pkg_L2.zip')
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
    (out / 'SOURCE_AUDIT_1_40.md').write_text(audit, encoding='utf-8')
    (out / 'README.md').write_text(readme, encoding='utf-8')
    package(out, Path(args.zip))
    print(json.dumps({
        'zip': args.zip,
        'out': args.out,
        'summary': [
            {
                'variant': r['variant'],
                'auto_beta': [r['auto_metrics'][f'beta{i}'] for i in range(4)],
                'pairings': r['automatic_pairings_applied'],
                'scalar_K_harmonic_ratio': r['auto_metrics']['harmonic_ratio'],
                'pair_transport_harmonic_ratio': r['pairing_transport_metrics']['pair_transport_harmonic_ratio'],
                'pair_scalar_harmonic_ratio': r['pairing_transport_metrics']['pair_scalar_harmonic_ratio'],
                'pair_kappa_orientation_ratio': r['pairing_transport_metrics']['pair_kappa_orientation_ratio'],
                'decision_used_delta_beta_any': r['pairing_transport_metrics']['decision_used_delta_beta_any'],
            } for r in rows
        ]
    }, indent=2))


if __name__ == '__main__':
    main()
