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


def skew_to_axial(M: np.ndarray) -> np.ndarray:
    # For skew matrix [[0,-z,y],[z,0,-x],[-y,x,0]], return (x,y,z).
    return np.array([
        0.5 * (M[2, 1] - M[1, 2]),
        0.5 * (M[0, 2] - M[2, 0]),
        0.5 * (M[1, 0] - M[0, 1]),
    ], dtype=float)


def face_normal(model: core.DynamicProvenanceGrowth, f: Face, mode: str = 'outward') -> np.ndarray:
    pts = [model.nodes[v].pos for v in f]
    if mode == 'birth_order':
        order = sorted(f, key=lambda v: (model.nodes[v].birth_time, v))
        pts = [model.nodes[v].pos for v in order]
    a, b, c = pts
    n = np.cross(b - a, c - a)
    nn = float(np.linalg.norm(n))
    if nn < EPS:
        return np.zeros(3)
    n = n / nn
    if mode in {'outward', 'birth_order'}:
        centroid = sum((model.nodes[v].pos for v in f), np.zeros(3)) / 3.0
        root = model.nodes[model.root].pos
        if float(np.dot(n, centroid - root)) < 0.0:
            n = -n
    return n


def face_area(model: core.DynamicProvenanceGrowth, f: Face) -> float:
    a, b, c = [model.nodes[v].pos for v in f]
    return 0.5 * float(np.linalg.norm(np.cross(b - a, c - a)))


def harmonic_basis_faces(K: core.SimplicialComplex) -> Tuple[np.ndarray, np.ndarray]:
    faces = K.faces()
    edges = K.edges()
    tets = sorted(K.tets)
    if not faces:
        return np.zeros((0, 0)), np.zeros(0)
    B2 = core.boundary_matrix_real([tuple(f) for f in faces], [tuple(e) for e in edges]) if edges else np.zeros((0, len(faces)))
    B3 = core.boundary_matrix_real([tuple(t) for t in tets], [tuple(f) for f in faces]) if tets else np.zeros((len(faces), 0))
    if not (B2.size and B3.size):
        return np.zeros((len(faces), 0)), np.zeros(0)
    L2 = B2.T @ B2 + B3 @ B3.T
    vals, vecs = np.linalg.eigh(L2)
    mask = vals < 1e-9
    return vecs[:, mask], vals


def vector_harmonic_projection(H: np.ndarray, W: np.ndarray) -> np.ndarray:
    if H.size == 0 or H.shape[1] == 0:
        return np.zeros_like(W)
    return H @ (H.T @ W)


def orientation_metrics(model: core.DynamicProvenanceGrowth, K: core.SimplicialComplex, source: str) -> dict:
    faces = K.faces()
    topo = core.topology(K)
    if not faces:
        return {
            'harmonic_dim_real': 0,
            'harmonic_axial_ratio': 0.0,
            'orientation_coherence': 0.0,
            'normal_flux_signed_ratio': 0.0,
            'normal_flux_abs_ratio': 0.0,
            'birth_normal_flux_signed_ratio': 0.0,
            'birth_normal_flux_abs_ratio': 0.0,
            'kappa_orientation_ratio': 0.0,
            'kappa_birth_orientation_ratio': 0.0,
            'support_fraction': 0.0,
        }
    W = np.array([skew_to_axial(core.face_K(model, f, source)) for f in faces], dtype=float)
    total = float(np.linalg.norm(W)) + EPS
    H, vals = harmonic_basis_faces(K)
    Wh = vector_harmonic_projection(H, W)
    hn = np.linalg.norm(Wh, axis=1)
    htotal = float(np.linalg.norm(Wh))
    mask = hn > max(1e-10, 1e-8 * (float(np.max(hn)) if len(hn) else 1.0))
    if np.any(mask):
        unit_vecs = Wh[mask] / (hn[mask][:, None] + EPS)
        coherence = float(np.linalg.norm(np.mean(unit_vecs, axis=0)))
    else:
        coherence = 0.0
    normals_out = np.array([face_normal(model, f, 'outward') for f in faces], dtype=float)
    normals_birth = np.array([face_normal(model, f, 'birth_order') for f in faces], dtype=float)
    areas = np.array([max(face_area(model, f), EPS) for f in faces], dtype=float)
    denom = float(np.sum(hn * areas)) + EPS
    signed_out = float(np.sum(np.einsum('ij,ij->i', Wh, normals_out) * areas)) / denom
    abs_out = float(np.sum(np.abs(np.einsum('ij,ij->i', Wh, normals_out)) * areas)) / denom
    signed_birth = float(np.sum(np.einsum('ij,ij->i', Wh, normals_birth) * areas)) / denom
    abs_birth = float(np.sum(np.abs(np.einsum('ij,ij->i', Wh, normals_birth)) * areas)) / denom
    return {
        'beta0': topo['beta0'], 'beta1': topo['beta1'], 'beta2': topo['beta2'], 'beta3': topo['beta3'],
        'harmonic_dim_real': int(H.shape[1]) if H.ndim == 2 else 0,
        'harmonic_axial_ratio': htotal / total,
        'orientation_coherence': coherence,
        'normal_flux_signed_ratio': signed_out,
        'normal_flux_abs_ratio': abs_out,
        'birth_normal_flux_signed_ratio': signed_birth,
        'birth_normal_flux_abs_ratio': abs_birth,
        'kappa_orientation_ratio': abs(signed_out) / (abs_out + EPS),
        'kappa_birth_orientation_ratio': abs(signed_birth) / (abs_birth + EPS),
        'support_fraction': float(np.mean(mask)) if len(mask) else 0.0,
        'harmonic_axis_mean_norm': float(np.mean(hn)) if len(hn) else 0.0,
        'harmonic_axis_p95_norm': float(np.percentile(hn, 95)) if len(hn) else 0.0,
        'laplacian_zero_eigs': int(H.shape[1]) if H.ndim == 2 else 0,
    }


def face_projection_rows(model: core.DynamicProvenanceGrowth, K: core.SimplicialComplex, source: str, topn: int = 80) -> List[dict]:
    faces = K.faces()
    if not faces:
        return []
    W = np.array([skew_to_axial(core.face_K(model, f, source)) for f in faces], dtype=float)
    H, _ = harmonic_basis_faces(K)
    Wh = vector_harmonic_projection(H, W)
    hn = np.linalg.norm(Wh, axis=1)
    rows = []
    for i, f in enumerate(faces):
        n_out = face_normal(model, f, 'outward')
        n_birth = face_normal(model, f, 'birth_order')
        rows.append({
            'face': str(list(f)),
            'birth_times': str([model.nodes[v].birth_time for v in f]),
            'birth_orders': str([model.nodes[v].birth_order for v in f]),
            'K_axial_norm': float(np.linalg.norm(W[i])),
            'H_axial_norm': float(hn[i]),
            'H_dot_outward_normal': float(np.dot(Wh[i], n_out)),
            'H_dot_birth_normal': float(np.dot(Wh[i], n_birth)),
            'area': face_area(model, f),
        })
    rows.sort(key=lambda r: r['H_axial_norm'], reverse=True)
    return rows[:topn]


def build_auto_variant(variant: str, args: argparse.Namespace, out: Path) -> Tuple[core.DynamicProvenanceGrowth, core.SimplicialComplex, dict, List[dict], List[dict]]:
    model = nl.build_model(variant, args)
    model.grow(args.max_level)
    baseline_K = core.build_dynamic_outward_ngf_complex(model)
    baseline_metrics = core.full_metrics(model, baseline_K, args.source)
    K, birth_log, pairing_log, candidate_sample, scans = nl.build_nonlinear_auto_complex(model, args, out / variant, variant)
    auto_metrics = core.full_metrics(model, K, args.source)
    return model, K, {'baseline_metrics': baseline_metrics, 'auto_metrics': auto_metrics, 'pairings': sum(1 for x in pairing_log if x.get('applied')), 'cascade_log': pairing_log}, birth_log, candidate_sample


def run_variant(variant: str, args: argparse.Namespace, out: Path) -> dict:
    vout = out / variant
    vout.mkdir(parents=True, exist_ok=True)
    model, K, base_summary, birth_log, candidate_sample = build_auto_variant(variant, args, out)
    orient = orientation_metrics(model, K, args.source)
    rows = face_projection_rows(model, K, args.source, args.keep_top_faces)
    write_csv(vout / 'harmonic_face_projection_top.csv', rows)
    summary = {
        'variant': variant,
        'max_level': args.max_level,
        'source': args.source,
        'response_mode': args.response_mode,
        **base_summary,
        'orientation_metrics': orient,
        'interpretation_flags': {
            'beta2_opened': base_summary['auto_metrics']['beta2'] > base_summary['baseline_metrics']['beta2'],
            'scalar_harmonic_positive': base_summary['auto_metrics']['harmonic_ratio'] > args.harmonic_positive_threshold,
            'axial_harmonic_positive': orient['harmonic_axial_ratio'] > args.harmonic_positive_threshold,
            'kappa_signed_bias_positive': orient['kappa_orientation_ratio'] > args.kappa_orientation_threshold,
            'birth_kappa_signed_bias_positive': orient['kappa_birth_orientation_ratio'] > args.kappa_orientation_threshold,
        },
    }
    (vout / 'variant_kappa_summary.json').write_text(json.dumps(summary, indent=2), encoding='utf-8')
    return summary


def make_docs(summary: dict) -> Tuple[str, str, str, str]:
    rows = summary['variant_rows']
    table_lines = ['| variant | beta auto | pairings | scalar harmonic | axial harmonic | Hdim | kappa out | kappa birth | signed out | signed birth |', '|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|']
    for r in rows:
        a = r['auto_metrics']; o = r['orientation_metrics']
        table_lines.append(f"| {r['variant']} | ({a['beta0']},{a['beta1']},{a['beta2']},{a['beta3']}) | {r['pairings']} | {a['harmonic_ratio']:.6g} | {o['harmonic_axial_ratio']:.6g} | {o['harmonic_dim_real']} | {o['kappa_orientation_ratio']:.6g} | {o['kappa_birth_orientation_ratio']:.6g} | {o['normal_flux_signed_ratio']:.6g} | {o['birth_normal_flux_signed_ratio']:.6g} |")
    table = '\n'.join(table_lines)
    smd = f"""# SUMMARY — harmonic K orientation / kappa gate

This package takes the nonlinear asymmetry-cascade growth result and tests the next question:

```text
beta2 opened -> harmonic K exists;
but does the harmonic K-component carry an oriented kappa/J-like bias?
```

It computes a vector-valued harmonic projection of the skew operator field.  The scalar harmonic ratio from the previous package used Frobenius norms on faces; this package also projects the axial vectors of the skew matrices into the harmonic 2-cochain space.

{table}

A positive result requires beta2>0, scalar and axial harmonic components >0, strict symmetry negative, and preferably a nonzero signed orientation bias.  The signed kappa metrics are diagnostic only; the face orientation convention is still a model choice and not yet a theorem.
"""
    rmd = f"""# RESULTS — harmonic K orientation / kappa gate

{table}

## What is measured

For every face, the skew matrix `K_face` is converted to an axial vector.  The vector field over faces is projected componentwise onto the real harmonic 2-cochain basis of the final complex:

```text
W_harm = H (H^T W)
```

Then the script measures:

- `harmonic_axial_ratio`: vector-valued harmonic K fraction.
- `orientation_coherence`: whether harmonic axial vectors align with one another.
- `normal_flux_signed_ratio`: signed coupling of harmonic K to outward-oriented face normals.
- `birth_normal_flux_signed_ratio`: same, but face orientation is seeded by vertex birth order before outward correction.
- `kappa_*_ratio`: absolute signed bias divided by absolute normal coupling.  Values near 1 are sign-coherent; values near 0 indicate cancellation.

## Critical reading

This is not yet a proof of a canonical J-sign.  It is a gate test.  If beta2 opens and harmonic K is nonzero but signed kappa ratios cancel, then topology has opened a real harmonic carrier but not yet an oriented complex structure.  If strict symmetry remains zero while real growth has a stable signed bias, the next step is a sign/flip test.

## Compact JSON

```json
{json.dumps(rows, indent=2)[:18000]}
```
"""
    amd = """# SOURCE_AUDIT_1_40

Carried-forward audit:

- Script 1/2: local sequential birth asymmetry and older-sibling environment are the first source of asymmetry.
- Script 35: K must be the real skew commutator sector, not a synthetic scalar field.  This package uses `face_K` from the core, then converts the skew matrix to its axial vector.
- Script 40: local tetrahedral closure alone was a Korand trap; beta2 must open before harmonic K can exist.
- Package 49: nonlinear asymmetry-gated complement pairing opens beta2 in real/sequential growth and is killed by strict symmetry.
- Current package: tests whether the opened harmonic K-sector is oriented or merely a real harmonic flow.
"""
    readme = """# Harmonic K orientation / kappa gate

Default:

```bash
python3 test_harmonic_k_orientation_kappa_gate.py
```

L3 local run:

```bash
python3 test_harmonic_k_orientation_kappa_gate.py \\
  --max-level 3 \\
  --max-cascade-per-birth 3 \\
  --max-auto-pairings 2 \\
  --out harmonic_kappa_out_L3 \\
  --zip cnna_harmonic_k_orientation_kappa_gate_pkg_L3.zip
```

This test is diagnostic.  A kappa sign bias depends on the current face-orientation convention and must later be replaced by a Lean-level/provenance-level orientation definition.
"""
    return smd, rmd, amd, readme


def write_comparative(out: Path, rows: List[dict]) -> None:
    compact = []
    for r in rows:
        a = r['auto_metrics']; b = r['baseline_metrics']; o = r['orientation_metrics']
        compact.append({
            'variant': r['variant'],
            'baseline_beta': f"({b['beta0']},{b['beta1']},{b['beta2']},{b['beta3']})",
            'auto_beta': f"({a['beta0']},{a['beta1']},{a['beta2']},{a['beta3']})",
            'pairings': r['pairings'],
            'scalar_harmonic_ratio': a['harmonic_ratio'],
            'axial_harmonic_ratio': o['harmonic_axial_ratio'],
            'harmonic_dim_real': o['harmonic_dim_real'],
            'orientation_coherence': o['orientation_coherence'],
            'normal_flux_signed_ratio': o['normal_flux_signed_ratio'],
            'normal_flux_abs_ratio': o['normal_flux_abs_ratio'],
            'birth_normal_flux_signed_ratio': o['birth_normal_flux_signed_ratio'],
            'birth_normal_flux_abs_ratio': o['birth_normal_flux_abs_ratio'],
            'kappa_orientation_ratio': o['kappa_orientation_ratio'],
            'kappa_birth_orientation_ratio': o['kappa_birth_orientation_ratio'],
            'support_fraction': o['support_fraction'],
            **{f"flag_{k}": v for k, v in r['interpretation_flags'].items()},
        })
    write_csv(out / 'comparative_kappa_orientation_summary.csv', compact)


def package(out: Path, zip_path: Path) -> None:
    if zip_path.exists():
        zip_path.unlink()
    base = Path(__file__).parent
    with zipfile.ZipFile(zip_path, 'w', zipfile.ZIP_DEFLATED) as z:
        for name in [
            'test_harmonic_k_orientation_kappa_gate.py',
            'test_nonlinear_asymmetry_cascade_growth.py',
            'cnna_non_shelling_core.py',
            'test_interfan_transport_from_asymmetry_invariants.py',
            'test_growth_with_asymmetry_gated_complement_pairing.py',
        ]:
            z.write(base / name, name)
        for p in out.rglob('*'):
            if p.is_file():
                z.write(p, str(p.relative_to(out.parent)))


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
    ap.add_argument('--out', default='harmonic_k_orientation_kappa_out_L2')
    ap.add_argument('--zip', default='cnna_harmonic_k_orientation_kappa_gate_pkg_L2.zip')
    args = ap.parse_args()

    out = Path(args.out)
    if out.exists():
        shutil.rmtree(out)
    out.mkdir(parents=True)
    rows = [run_variant(v, args, out) for v in args.variants]
    summary = {'args': vars(args), 'variant_rows': rows}
    (out / 'comparative_summary.json').write_text(json.dumps(summary, indent=2), encoding='utf-8')
    write_comparative(out, rows)
    smd, rmd, amd, readme = make_docs(summary)
    (out / 'SUMMARY.md').write_text(smd, encoding='utf-8')
    (out / 'RESULTS.md').write_text(rmd, encoding='utf-8')
    (out / 'SOURCE_AUDIT_1_40.md').write_text(amd, encoding='utf-8')
    (out / 'README.md').write_text(readme, encoding='utf-8')
    package(out, Path(args.zip))
    print(json.dumps({
        'zip': args.zip,
        'out': args.out,
        'summary': [
            {
                'variant': r['variant'],
                'auto_beta': [r['auto_metrics'][f'beta{i}'] for i in range(4)],
                'pairings': r['pairings'],
                'scalar_harmonic_ratio': r['auto_metrics']['harmonic_ratio'],
                'axial_harmonic_ratio': r['orientation_metrics']['harmonic_axial_ratio'],
                'kappa_orientation_ratio': r['orientation_metrics']['kappa_orientation_ratio'],
                'kappa_birth_orientation_ratio': r['orientation_metrics']['kappa_birth_orientation_ratio'],
            } for r in rows
        ]
    }, indent=2))


if __name__ == '__main__':
    main()
