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
import test_dual_assembly_order_context_ablation_gate as p69
import test_signed_Jlock_role_coupling_gate as p70
import test_assembly_motif_basis_diagonalization_gate as p71
import test_edge_interface_motif_operator_gate as p72
import test_real_symplectic_before_star_gate as p74
import test_kahler_compatibility_star_gate as p75

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


def orth(cols: List[np.ndarray], tol: float=1e-10) -> np.ndarray:
    return p71.orth_basis(cols, tol=tol)


def projector(U: np.ndarray) -> np.ndarray:
    if U.size == 0 or U.shape[1] == 0:
        return np.zeros((U.shape[0], U.shape[0]), dtype=float)
    return U @ U.T


def subspace_image_residual(J: np.ndarray, source_cols: List[np.ndarray], target_cols: List[np.ndarray]) -> float:
    return p71.subspace_image_residual(J, source_cols, target_cols)


def signed_scaled_residual(M: np.ndarray, target_sign: int) -> tuple[float, float]:
    # Fit M ≈ target_sign * alpha * I.  For target_sign=-1 this is the spinor 2π sign test.
    # Zero operators are NOT accepted as projective closure.
    d = M.shape[0]
    I = np.eye(d)
    if d == 0 or norm(M) < EPS:
        return 1.0, 0.0
    alpha = float(target_sign * np.trace(M) / (d + EPS))
    denom = norm(M) + abs(alpha) * norm(I) + EPS
    resid = norm(M - target_sign * alpha * I) / denom
    return float(resid), float(alpha)


def vector_scaled_residual(M: np.ndarray, cols: List[np.ndarray], target_sign: int) -> tuple[float, float]:
    U = orth(cols)
    if U.shape[1] == 0:
        return 0.0, 0.0
    A = U.T @ M @ U
    return signed_scaled_residual(A, target_sign)


def matrix_power_metrics(T: np.ndarray, qproj: List[np.ndarray], pproj: List[np.ndarray], args: argparse.Namespace) -> dict:
    d = T.shape[0]
    I = np.eye(d)
    T2 = T @ T
    T4 = T2 @ T2
    T_norm = norm(T)
    T2_norm = norm(T2)
    T4_norm = norm(T4)
    minus_resid, alpha = signed_scaled_residual(T2, -1)
    plus4_alpha = alpha * alpha
    plus4_resid = norm(T4 - plus4_alpha * I) / (norm(T4) + abs(plus4_alpha) * norm(I) + EPS)
    raw_J2 = norm(T2 + I) / (norm(I) + EPS)
    raw_J4 = norm(T4 - I) / (norm(I) + EPS)
    q2_resid, q2_alpha = vector_scaled_residual(T2, qproj, -1)
    p2_resid, p2_alpha = vector_scaled_residual(T2, pproj, -1)
    q4_resid, q4_alpha = vector_scaled_residual(T4, qproj, +1)
    p4_resid, p4_alpha = vector_scaled_residual(T4, pproj, +1)
    q_to_p = subspace_image_residual(T, qproj, pproj)
    p_to_q = subspace_image_residual(T, pproj, qproj)
    eigvals = np.linalg.eigvals(T) if d else np.array([], dtype=complex)
    eig_imag_ratio = 0.0
    if len(eigvals):
        eig_imag_ratio = float(np.mean(np.abs(np.imag(eigvals)) / (np.abs(eigvals) + EPS)))
    nontrivial = bool(T_norm >= args.min_operator_norm and T2_norm >= args.min_operator_norm)
    spinor_gate = bool(
        nontrivial
        and alpha > args.min_spinor_alpha
        and minus_resid <= args.spinor_T2_threshold
        and plus4_resid <= args.spinor_T4_threshold
        and max(q2_resid, p2_resid) <= args.spinor_subspace_threshold
        and 0.5 * (q_to_p + p_to_q) <= args.spinor_QP_lock_threshold
    )
    # Strong spinor would be unscaled: T²≈-I and T⁴≈I.  Projective spinor only requires alpha scaling.
    strong_gate = bool(
        nontrivial
        and raw_J2 <= args.strong_J2_threshold
        and raw_J4 <= args.strong_J4_threshold
        and max(q2_resid, p2_resid) <= args.spinor_subspace_threshold
    )
    return {
        'basis_dim': d,
        'T_norm': float(T_norm),
        'T2_norm': float(T2_norm),
        'T4_norm': float(T4_norm),
        'nontrivial_operator': bool(nontrivial),
        'T_Q_to_P_resid': float(q_to_p),
        'T_P_to_Q_resid': float(p_to_q),
        'T_QP_mean_resid': float(0.5 * (q_to_p + p_to_q)),
        'T2_plus_I_raw_resid': float(raw_J2),
        'T4_minus_I_raw_resid': float(raw_J4),
        'T2_projective_minus_I_resid': float(minus_resid),
        'T2_projective_alpha': float(alpha),
        'T4_projective_plus_I_resid': float(plus4_resid),
        'T2_Q_minus_resid': float(q2_resid),
        'T2_Q_alpha': float(q2_alpha),
        'T2_P_minus_resid': float(p2_resid),
        'T2_P_alpha': float(p2_alpha),
        'T4_Q_plus_resid': float(q4_resid),
        'T4_Q_alpha': float(q4_alpha),
        'T4_P_plus_resid': float(p4_resid),
        'T4_P_alpha': float(p4_alpha),
        'eig_imag_ratio_mean': float(eig_imag_ratio),
        'eigvals_json': json.dumps([[float(np.real(z)), float(np.imag(z))] for z in eigvals]),
        'spinor_projective_gate_pass': spinor_gate,
        'spinor_strong_unscaled_gate_pass': strong_gate,
    }


def spinor_rows_for_assembly(model, K, assembly_row: dict, args: argparse.Namespace) -> List[dict]:
    parsed = p74.parse_assembly_pairs(model, assembly_row, args)
    if parsed is None:
        return []
    pA, pB = parsed
    faces = p72.union_faces_from_pairs(pA, pB)
    Jbase, Cbase = p71.union_JC(faces, [pA, pB])
    vecs, qcols, pcols, allcols, U = p75.projected_data(faces, pA, pB)
    n = len(faces) * 3
    rows: List[dict] = []
    base = {
        'context': assembly_row.get('context',''),
        'A_face_a': str(list(pA['fa'])), 'A_face_b': str(list(pA['fb'])),
        'B_face_a': str(list(pB['fa'])), 'B_face_b': str(list(pB['fb'])),
        'union_faces': str([list(f) for f in faces]),
    }
    omegas = p75.omega_candidates(model, K, faces, pA, pB, args)
    metrics = p75.metric_candidates(n, U, faces, pA, pB, Jbase, Cbase, qcols, pcols, args)
    for oname, Omega, oprim, oreason, oextra in omegas:
        for gname, G, gprim, greason in metrics:
            all_basis = p75.orth(allcols)
            if all_basis.shape[1] == 0:
                continue
            Op = all_basis.T @ Omega @ all_basis
            Gp = all_basis.T @ G @ all_basis
            try:
                T = np.linalg.pinv(Gp, rcond=1e-10) @ Op
            except Exception:
                continue
            qproj = [all_basis.T @ q for q in qcols]
            pproj = [all_basis.T @ p for p in pcols]
            compat = p75.compatibility_metrics('compat', Omega, G, Cbase, qcols, pcols, allcols, args)
            spins = matrix_power_metrics(T, qproj, pproj, args)
            rows.append({
                **base,
                'omega_candidate': oname,
                'omega_primary': bool(oprim),
                'omega_source_note': oreason,
                'metric_candidate': gname,
                'metric_primary': bool(gprim),
                'metric_source_note': greason,
                'primary_pair': bool(oprim and gprim),
                **oextra,
                **compat,
                **spins,
            })
    return rows


def summarize_spinor(rows: List[dict]) -> dict:
    if not rows:
        return {
            'spinor_row_count': 0,
            'primary_spinor_row_count': 0,
            'projective_spinor_pass_count': 0,
            'primary_projective_spinor_pass_count': 0,
            'strong_spinor_pass_count': 0,
            'primary_strong_spinor_pass_count': 0,
            'best_primary_T2_projective_resid': 0.0,
            'best_primary_T4_projective_resid': 0.0,
            'best_primary_alpha': 0.0,
            'best_primary_raw_J2': 0.0,
            'best_primary_candidate': '',
            'best_primary_metric': '',
        }
    primary = [r for r in rows if fbool(r.get('primary_pair'))]
    def count(rs, key): return sum(1 for r in rs if fbool(r.get(key)))
    def score(r):
        # We score actual double-cover behavior, not isolated Q->P map.
        if (not fbool(r.get('nontrivial_operator'))) or float(r.get('T2_projective_alpha',0.0)) <= 0:
            return (100.0, 100.0, 100.0, 100.0)
        return (
            float(r.get('T2_projective_minus_I_resid', 99.0)),
            float(r.get('T4_projective_plus_I_resid', 99.0)),
            float(r.get('T_QP_mean_resid', 99.0)),
            float(r.get('T2_plus_I_raw_resid', 99.0)),
        )
    bestp = min(primary, key=score, default=None)
    besta = min(rows, key=score, default=None)
    def pick(r,k,d=0.0): return float(r.get(k,d)) if r else d
    return {
        'spinor_row_count': len(rows),
        'primary_spinor_row_count': len(primary),
        'projective_spinor_pass_count': count(rows, 'spinor_projective_gate_pass'),
        'primary_projective_spinor_pass_count': count(primary, 'spinor_projective_gate_pass'),
        'strong_spinor_pass_count': count(rows, 'spinor_strong_unscaled_gate_pass'),
        'primary_strong_spinor_pass_count': count(primary, 'spinor_strong_unscaled_gate_pass'),
        'best_primary_T2_projective_resid': pick(bestp, 'T2_projective_minus_I_resid'),
        'best_primary_T4_projective_resid': pick(bestp, 'T4_projective_plus_I_resid'),
        'best_primary_alpha': pick(bestp, 'T2_projective_alpha'),
        'best_primary_QP_mean_resid': pick(bestp, 'T_QP_mean_resid'),
        'best_primary_raw_J2': pick(bestp, 'T2_plus_I_raw_resid'),
        'best_primary_raw_J4': pick(bestp, 'T4_minus_I_raw_resid'),
        'best_primary_T2_Q_minus_resid': pick(bestp, 'T2_Q_minus_resid'),
        'best_primary_T2_P_minus_resid': pick(bestp, 'T2_P_minus_resid'),
        'best_primary_eig_imag_ratio': pick(bestp, 'eig_imag_ratio_mean'),
        'best_primary_candidate': bestp.get('omega_candidate','') if bestp else '',
        'best_primary_metric': bestp.get('metric_candidate','') if bestp else '',
        'best_any_T2_projective_resid': pick(besta, 'T2_projective_minus_I_resid'),
        'best_any_T4_projective_resid': pick(besta, 'T4_projective_plus_I_resid'),
        'best_any_alpha': pick(besta, 'T2_projective_alpha'),
        'best_any_candidate': besta.get('omega_candidate','') if besta else '',
        'best_any_metric': besta.get('metric_candidate','') if besta else '',
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
    dm, pair_rows, top_rows, three_rows = p75.p56.directed_metrics(model, K, pairing_log, args)
    sm, signed_rows, signed_face_rows = p69.p59.signed_quadrature_rows(model, K, pairing_log, args)
    am, align_pair_rows, align_candidate_rows, align_candidate_summary = p75.p61.alignment_search_metrics(model, K, pairing_log, args)
    sel = p69.summarize_selection(pairing_log)
    motif_rows, motif_summary = p71.assembly_motif_rows(model, K, pairing_log, assembly_log, args)
    compat_rows: List[dict] = []
    for i, a in enumerate(assembly_log):
        if not fbool(a.get('assembly_applied')):
            continue
        for r in spinor_rows_for_assembly(model, K, a, args):
            r['assembly_index'] = i
            compat_rows.append(r)
    spin_summary = summarize_spinor(compat_rows)
    write_csv(vout / 'birth_geometry_log.csv', birth_log)
    write_csv(vout / 'assembly_pairing_log.csv', pairing_log)
    write_csv(vout / 'assembly_ablation_log.csv', assembly_log)
    write_csv(vout / 'candidate_eval_rows.csv', candidate_rows)
    write_csv(vout / 'spinor_double_cover_rows.csv', compat_rows)
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
        'spinor_double_cover_metrics': spin_summary,
        'automatic_pairings_applied': sum(1 for x in pairing_log if fbool(x.get('applied'))),
        'assemblies_applied': sum(1 for x in assembly_log if fbool(x.get('assembly_applied'))),
        'assemblies_attempted': len(assembly_log),
        'decision_used_delta_beta_any': sel['decision_used_delta_beta_any'],
    }
    (vout / 'variant_spinor_double_cover_summary.json').write_text(json.dumps(summary, indent=2), encoding='utf-8')
    return summary


def run_experiment(args: argparse.Namespace, out: Path) -> List[dict]:
    p70.patch_modules()
    rows = []
    opts = []
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
                opts.append(clone_args(args, **updates))
    for opt in opts:
        for variant in opt.variants:
            rows.append(run_variant(variant, opt, out))
    return rows


def slim(r: dict) -> dict:
    a = r['auto_metrics']; dm = r['directed_metrics']; sm = r['signed_quadrature_metrics']; am = r['alignment_metrics']; mm = r['motif_basis_metrics']; sp = r['spinor_double_cover_metrics']
    return {
        'option': r['option'], 'variant': r['variant'],
        'beta': [a['beta0'], a['beta1'], a['beta2'], a['beta3']],
        'pairings': r['automatic_pairings_applied'], 'assemblies': r['assemblies_applied'],
        'pair_harm': dm['pair_transport_harmonic_ratio'],
        'Q_harm': am['Q_even_harmonic_ratio'], 'P_harm': am['P_odd_harmonic_ratio'],
        'pair_local_J_lock': am['best_per_pair_mean_J_lock_resid'],
        'motif_union_J_lock': mm.get('union_sum_best_mean_resid',0.0),
        'signed_birth': sm['signed_birth_over_abs_sum_ratio'],
        'spinor_projective_pass': sp.get('primary_projective_spinor_pass_count',0),
        'spinor_strong_pass': sp.get('primary_strong_spinor_pass_count',0),
        'best_T2_proj_resid': sp.get('best_primary_T2_projective_resid',0.0),
        'best_T4_proj_resid': sp.get('best_primary_T4_projective_resid',0.0),
        'best_alpha': sp.get('best_primary_alpha',0.0),
        'best_QP_mean': sp.get('best_primary_QP_mean_resid',0.0),
        'best_raw_J2': sp.get('best_primary_raw_J2',0.0),
        'best_raw_J4': sp.get('best_primary_raw_J4',0.0),
        'best_T2_Q_minus': sp.get('best_primary_T2_Q_minus_resid',0.0),
        'best_T2_P_minus': sp.get('best_primary_T2_P_minus_resid',0.0),
        'best_eig_imag_ratio': sp.get('best_primary_eig_imag_ratio',0.0),
        'best_candidate': sp.get('best_primary_candidate',''),
        'best_metric': sp.get('best_primary_metric',''),
        'used_delta_beta': r['decision_used_delta_beta_any'],
    }


def write_comparative(out: Path, rows: List[dict]) -> None:
    flat = []
    for r in rows:
        s = slim(r)
        flat.append({**s, 'beta0': s['beta'][0], 'beta1': s['beta'][1], 'beta2': s['beta'][2], 'beta3': s['beta'][3]})
    write_csv(out / 'comparative_spinor_double_cover_summary.csv', flat)


def make_docs(summary: dict) -> tuple[str,str,str,str]:
    rows = summary['variant_rows']
    lines = ['| option | variant | beta | pairs | asm | Q harm | P harm | pair J-lock | motif lock | projective pass | strong pass | T²≈-αI | α | T⁴≈α²I | raw J²+I | QP mean | T²Q≈-Q | T²P≈-P | best Ω/g | used dβ? |',
             '|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|---:|']
    for r in rows:
        s = slim(r)
        lines.append(f"| {s['option']} | {s['variant']} | ({s['beta'][0]},{s['beta'][1]},{s['beta'][2]},{s['beta'][3]}) | {s['pairings']} | {s['assemblies']} | {s['Q_harm']:.6g} | {s['P_harm']:.6g} | {s['pair_local_J_lock']:.6g} | {s['motif_union_J_lock']:.6g} | {s['spinor_projective_pass']} | {s['spinor_strong_pass']} | {s['best_T2_proj_resid']:.6g} | {s['best_alpha']:.6g} | {s['best_T4_proj_resid']:.6g} | {s['best_raw_J2']:.6g} | {s['best_QP_mean']:.6g} | {s['best_T2_Q_minus']:.6g} | {s['best_T2_P_minus']:.6g} | {s['best_candidate']} / {s['best_metric']} | {s['used_delta_beta']} |")
    table = '\n'.join(lines)
    nonstrict = [r for r in rows if r['variant'] != 'strict_symmetrized_control']
    any_proj = any(r['spinor_double_cover_metrics'].get('primary_projective_spinor_pass_count',0) > 0 for r in nonstrict)
    any_strong = any(r['spinor_double_cover_metrics'].get('primary_strong_spinor_pass_count',0) > 0 for r in nonstrict)
    smd = f"""# SUMMARY — Spinor / double-cover closure gate

Model tag: `CQNM/s=-1 saturated geometry reference, provenance-growth L2 diagnostic`.

This package tests the user's spin-1/2 suspicion in a strict derived-only form.  It does not set spin.  For each existing Ω/g-derived operator `T = g^-1 Ω` on the actual Q/P motif span it audits whether the operator is projectively order-four:

```text
T(Q) -> P,
T² ≈ -α I,
T⁴ ≈ α² I,
```

with α > 0.  This is the real double-cover/spinor-like signature: a sign after the second application and closure after the fourth application.  The strong unscaled gate additionally requires `T²≈-I` and `T⁴≈I`.

{table}

Decision:

```json
{{
  "any_non_strict_primary_projective_spinor_pass": {str(any_proj).lower()},
  "any_non_strict_primary_strong_unscaled_spinor_pass": {str(any_strong).lower()}
}}
```
"""
    rmd = f"""# RESULTS — Spinor / double-cover closure gate

## Comparative table

{table}

## Gate definition

A primary projective spinor row passes only if:

```text
α > min_spinor_alpha,
T² is close to -α I,
T⁴ is close to α² I,
T² maps Q and P back to their own subspaces with negative sign,
T maps Q/P subspaces into each other sufficiently well.
```

This distinguishes three situations:

1. ordinary failed J: Q/P may exist, but T² is not scalar negative;
2. projective spinor-like closure: T²≈-αI and T⁴≈α²I, even if α≠1;
3. strong complex/spinor closure: T²≈-I and T⁴≈I.

No normalization or polar correction is used to force closure.
"""
    audit = """# SOURCE AUDIT

No i, global J, Hodge star, physical adjoint, positivity, C*-norm, final sym(M), or delta-beta/H2 decision is introduced.

This test reuses the Ω/g candidates from the Kähler compatibility gate and asks a different question: maybe the operator is not a good unscaled J but still closes projectively like a double-cover/spinor object.

Important limitation: the Ω candidates are not independent of earlier pair-exchange J-like operators.  The test therefore evaluates closure properties of those existing derived structural operators; it does not prove spin.

This is a Python L2 diagnostic, not a Lean theorem.
"""
    readme = """# Spinor / double-cover closure gate

Run:

```bash
python3 test_spinor_double_cover_closure_gate.py
```

The package tests whether the derived Q/P operator behaves projectively as order four: T²≈-αI and T⁴≈α²I.
"""
    return smd, rmd, audit, readme


def package(out: Path, zip_path: Path) -> None:
    files = [
        Path(__file__).name,
        'test_kahler_compatibility_star_gate.py',
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
    # Compatibility args reused by p75.
    ap.add_argument('--compat-J2-threshold', type=float, default=0.20)
    ap.add_argument('--compat-lock-mean-threshold', type=float, default=0.20)
    ap.add_argument('--compat-lock-max-threshold', type=float, default=0.30)
    ap.add_argument('--metric-orth-threshold', type=float, default=0.25)
    ap.add_argument('--hash-anti-threshold', type=float, default=0.25)
    ap.add_argument('--star-span-threshold', type=float, default=0.25)
    # Spinor/double-cover thresholds.
    ap.add_argument('--min-spinor-alpha', type=float, default=1e-6)
    ap.add_argument('--min-operator-norm', type=float, default=1e-8)
    ap.add_argument('--spinor-T2-threshold', type=float, default=0.20)
    ap.add_argument('--spinor-T4-threshold', type=float, default=0.20)
    ap.add_argument('--spinor-subspace-threshold', type=float, default=0.25)
    ap.add_argument('--spinor-QP-lock-threshold', type=float, default=0.25)
    ap.add_argument('--strong-J2-threshold', type=float, default=0.20)
    ap.add_argument('--strong-J4-threshold', type=float, default=0.20)
    # Args required by p74.symplectic_rows although this script does not use it directly.
    ap.add_argument('--symplectic-ratio-threshold', type=float, default=1e-3)
    ap.add_argument('--qp-ratio-threshold', type=float, default=1e-3)
    ap.add_argument('--isotropic-threshold', type=float, default=0.35)
    ap.add_argument('--out', default='spinor_double_cover_closure_out_L2')
    ap.add_argument('--zip', default='cnna_spinor_double_cover_closure_gate_pkg_L2.zip')
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
