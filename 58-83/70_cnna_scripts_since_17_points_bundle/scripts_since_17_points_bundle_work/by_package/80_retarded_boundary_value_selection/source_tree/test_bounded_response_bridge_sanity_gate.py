#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import math
import shutil
import zipfile
from pathlib import Path
from typing import List, Optional

import numpy as np

import cnna_non_shelling_core as core
import test_nonlinear_asymmetry_cascade_growth as nl
import test_dual_assembly_order_context_ablation_gate as p69
import test_signed_Jlock_role_coupling_gate as p70
import test_assembly_motif_basis_diagonalization_gate as p71
import test_real_symplectic_before_star_gate as p74
import test_kahler_compatibility_star_gate as p75
import test_pairing_transport_antisym_birth_coherence_gate as p56
import test_signed_quadrature_area_kappa_gate as p59
import test_pair_J_alignment_search_gate as p61
import test_response_monodromy_to_QP_transfer_gate as p78

EPS = 1e-12


def fbool(x) -> bool:
    if isinstance(x, str):
        return x.lower() in {'true','1','yes'}
    return bool(x)


def norm(x: np.ndarray) -> float:
    return float(np.linalg.norm(np.asarray(x, dtype=float)))


def opnorm(A: np.ndarray) -> float:
    if A.size == 0:
        return 0.0
    try:
        s = np.linalg.svd(A, compute_uv=False)
        return float(s[0]) if len(s) else 0.0
    except np.linalg.LinAlgError:
        return float('inf')


def spectral_radius(A: np.ndarray) -> float:
    if A.size == 0:
        return 0.0
    try:
        eig = np.linalg.eigvals(A)
        return float(max(abs(z) for z in eig)) if eig.size else 0.0
    except np.linalg.LinAlgError:
        return float('inf')


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


def residual_to_sign(A: np.ndarray, sign: int) -> float:
    if A.size == 0:
        return 1.0
    I = np.eye(A.shape[0])
    return norm(A - sign*I) / (norm(A) + norm(I) + EPS)


def matrix_diagnostics(A: np.ndarray, prefix: str) -> dict:
    if A.size == 0:
        return {
            f'{prefix}_opnorm': 0.0, f'{prefix}_fro_norm': 0.0, f'{prefix}_spectral_radius': 0.0,
            f'{prefix}_normality_resid': 1.0, f'{prefix}_skew_ratio': 0.0, f'{prefix}_sym_ratio': 0.0,
            f'{prefix}_max_real_eig': 0.0, f'{prefix}_min_real_eig': 0.0, f'{prefix}_max_abs_imag_eig': 0.0,
            f'{prefix}_eigvals_json': '[]', f'{prefix}_singular_values_json': '[]'
        }
    try:
        eig = np.linalg.eigvals(A)
        s = np.linalg.svd(A, compute_uv=False)
        normality = norm(A.T @ A - A @ A.T) / (norm(A.T @ A) + norm(A @ A.T) + EPS)
        skew_ratio = norm(0.5*(A-A.T)) / (norm(A)+EPS)
        sym_ratio = norm(0.5*(A+A.T)) / (norm(A)+EPS)
        return {
            f'{prefix}_opnorm': float(s[0]) if len(s) else 0.0,
            f'{prefix}_fro_norm': norm(A),
            f'{prefix}_spectral_radius': float(max(abs(z) for z in eig)) if eig.size else 0.0,
            f'{prefix}_normality_resid': normality,
            f'{prefix}_skew_ratio': skew_ratio,
            f'{prefix}_sym_ratio': sym_ratio,
            f'{prefix}_max_real_eig': float(max(np.real(z) for z in eig)) if eig.size else 0.0,
            f'{prefix}_min_real_eig': float(min(np.real(z) for z in eig)) if eig.size else 0.0,
            f'{prefix}_max_abs_imag_eig': float(max(abs(np.imag(z)) for z in eig)) if eig.size else 0.0,
            f'{prefix}_eigvals_json': json.dumps([[float(np.real(z)), float(np.imag(z))] for z in eig]),
            f'{prefix}_singular_values_json': json.dumps([float(x) for x in s]),
        }
    except np.linalg.LinAlgError:
        return {
            f'{prefix}_opnorm': float('inf'), f'{prefix}_fro_norm': norm(A), f'{prefix}_spectral_radius': float('inf'),
            f'{prefix}_normality_resid': float('inf'), f'{prefix}_skew_ratio': 0.0, f'{prefix}_sym_ratio': 0.0,
            f'{prefix}_max_real_eig': 0.0, f'{prefix}_min_real_eig': 0.0, f'{prefix}_max_abs_imag_eig': 0.0,
            f'{prefix}_eigvals_json': '[]', f'{prefix}_singular_values_json': '[]'
        }


def projected_metrics(T: np.ndarray, qcols: List[np.ndarray], pcols: List[np.ndarray], allcols: List[np.ndarray], args: argparse.Namespace) -> dict:
    qp_cols = qcols + pcols
    q_to_p = p78.subspace_image_residual(T, qcols, pcols)
    p_to_q = p78.subspace_image_residual(T, pcols, qcols)
    mean_lock = 0.5*(q_to_p+p_to_q)
    max_lock = max(q_to_p,p_to_q)
    j2 = p78.projected_J2_residual(T, allcols)
    spec = p78.projected_spectrum_metrics(T, qp_cols, 'projected_QP')
    U = p78.orth(qp_cols)
    A = U.T @ T @ U if U.shape[1] else np.zeros((0,0))
    diag = matrix_diagnostics(A, 'projected_QP')
    # bounded/semigroup style audits; no interpretation as i/J.
    contraction = bool(diag['projected_QP_opnorm'] <= args.bounded_opnorm_threshold)
    complex_like_stable = bool(
        spec['projected_QP_eig_imag_max'] >= args.complex_imag_threshold and
        diag['projected_QP_normality_resid'] <= args.normality_threshold and
        diag['projected_QP_opnorm'] <= args.bounded_opnorm_threshold
    )
    c3_like_stable = bool(
        spec['projected_QP_cube_identity_resid'] <= args.c3_identity_threshold and
        spec['projected_QP_identity_resid'] >= args.not_identity_threshold and
        diag['projected_QP_opnorm'] <= args.bounded_opnorm_threshold
    )
    minus_like_stable = bool(
        spec['projected_QP_minus_identity_resid'] <= args.minus_identity_threshold and
        diag['projected_QP_opnorm'] <= args.bounded_opnorm_threshold
    )
    qp_lock_stable = bool(mean_lock <= args.qp_lock_threshold and max_lock <= args.qp_lock_max_threshold and diag['projected_QP_opnorm'] <= args.bounded_opnorm_threshold)
    # crude half-plane / semigroup audit: does the projected transfer have a preferred dissipative real side?
    eigvals = []
    try:
        eigvals = list(np.linalg.eigvals(A)) if A.size else []
    except np.linalg.LinAlgError:
        eigvals = []
    pos_real = sum(1 for z in eigvals if np.real(z) > args.real_halfplane_tol)
    neg_real = sum(1 for z in eigvals if np.real(z) < -args.real_halfplane_tol)
    halfplane_bias = (pos_real - neg_real) / (len(eigvals) + EPS) if eigvals else 0.0
    return {
        'Q_to_P_transfer_resid': q_to_p,
        'P_to_Q_transfer_resid': p_to_q,
        'QP_transfer_mean_resid': mean_lock,
        'QP_transfer_max_resid': max_lock,
        'projected_J2_plus_I_resid': j2,
        'bounded_contraction_audit_pass': contraction,
        'QP_lock_stable_pass': qp_lock_stable,
        'complex_like_stable_pass': complex_like_stable,
        'C3_like_stable_pass': c3_like_stable,
        'minus_identity_stable_pass': minus_like_stable,
        'halfplane_pos_real_count': pos_real,
        'halfplane_neg_real_count': neg_real,
        'halfplane_bias_pos_minus_neg': halfplane_bias,
        **spec,
        **diag,
    }


def bounded_lift_rows(model, K, assembly_row: dict, args: argparse.Namespace) -> List[dict]:
    parsed = p74.parse_assembly_pairs(model, assembly_row, args)
    if parsed is None:
        return []
    pA, pB = parsed
    faces = p78.p72.union_faces_from_pairs(pA, pB)
    vecs, qcols, pcols, allcols, U = p75.projected_data(faces, pA, pB)
    parents = p78.motif_parent_fans(model, faces, args)
    base = {
        'event_t': assembly_row.get('event_t',''), 'scan_id': assembly_row.get('scan_id',''),
        'context': assembly_row.get('context',''), 'face_overlap': assembly_row.get('face_overlap',''),
        'edge_overlap': assembly_row.get('edge_overlap',''), 'vertex_overlap': assembly_row.get('vertex_overlap',''),
        'A_face_a': str(list(pA['fa'])), 'A_face_b': str(list(pA['fb'])),
        'B_face_a': str(list(pB['fa'])), 'B_face_b': str(list(pB['fb'])),
        'union_faces': str([list(f) for f in faces]),
        'motif_basis_rank': int(U.shape[1]), 'Q_rank': int(p78.orth(qcols).shape[1]), 'P_rank': int(p78.orth(pcols).shape[1]),
        'candidate_parent_fans': str(parents), 'parent_fan_count': len(parents),
    }
    if not parents or U.shape[1] == 0:
        return [{**base, 'transfer_status': 'no_parent_fan_or_basis'}]
    rows = []
    for bridge_mode in args.bridge_modes:
        B, brow = p78.build_bridge_matrix(model, faces, parents, bridge_mode)
        if B.size == 0 or p78.norm(B) < EPS:
            continue
        s = np.linalg.svd(B, compute_uv=False)
        max_sv = float(s[0]) if len(s) else 0.0
        min_sv = float(s[-1]) if len(s) else 0.0
        rank = int(np.sum(s > args.bridge_rank_tol * max(max_sv, 1.0)))
        cond = float(max_sv / (min_sv + EPS)) if len(s) else float('inf')
        well_conditioned = bool(cond <= args.max_bridge_cond and rank >= args.min_bridge_rank and min_sv >= args.min_bridge_sv)
        Binv = np.linalg.pinv(B, rcond=args.bridge_pinv_rcond)
        for op_mode in args.response_ops:
            R, rstats = p78.block_response_operator(model, parents, op_mode)
            if R.size == 0:
                continue
            R_opnorm = opnorm(R)
            R_radius = spectral_radius(R)
            response_active = bool(rstats['max_abs_response_circulation'] > args.min_response_circulation and R_opnorm > EPS)
            response_normalizations = []
            response_normalizations.append(('raw_unbounded_audit', R, False))
            response_normalizations.append(('response_opnorm_bounded', R / max(R_opnorm, 1.0), True))
            response_normalizations.append(('response_radius_bounded', R / max(R_radius, 1.0), True))
            for norm_mode, Rb, bounded_mode in response_normalizations:
                T = Binv @ Rb @ B
                T_raw_norm = opnorm(T)
                variants = [(norm_mode, T, bounded_mode)]
                if args.include_transfer_post_bounded:
                    variants.append((norm_mode + '__transfer_post_opnorm_bounded', T / max(T_raw_norm, 1.0), True))
                for lift_mode, Tb, is_bounded in variants:
                    pmet = projected_metrics(Tb, qcols, pcols, allcols, args)
                    active = bool(response_active and p78.norm(Tb) > args.min_transfer_norm)
                    artifact_suspect = bool((not well_conditioned) or cond > args.artifact_cond_threshold or T_raw_norm > args.artifact_raw_opnorm_threshold)
                    structural_pass = bool(
                        active and well_conditioned and is_bounded and not artifact_suspect and (
                            pmet['QP_lock_stable_pass'] or pmet['complex_like_stable_pass'] or pmet['C3_like_stable_pass'] or pmet['minus_identity_stable_pass']
                        )
                    )
                    rows.append({
                        **base,
                        'transfer_status': 'ok', 'bridge_mode': bridge_mode, 'response_op': op_mode,
                        'bridge_rows': str(brow),
                        'lift_mode': lift_mode,
                        'response_bounded_mode': is_bounded,
                        'transfer_active': active,
                        'bridge_rank': rank, 'bridge_min_sv': min_sv, 'bridge_max_sv': max_sv, 'bridge_cond': cond,
                        'bridge_well_conditioned_pass': well_conditioned,
                        'response_opnorm': R_opnorm, 'response_spectral_radius': R_radius,
                        'response_active': response_active,
                        'transfer_opnorm_raw_before_post_bound': T_raw_norm,
                        'artifact_suspect': artifact_suspect,
                        'bounded_structural_transfer_pass': structural_pass,
                        **rstats,
                        **matrix_diagnostics(Tb, 'lifted_transfer'),
                        **pmet,
                    })
    return rows


def slim_row(r: Optional[dict]) -> dict:
    if not r:
        return {}
    keys = ['bridge_mode','response_op','lift_mode','context','bridge_rank','bridge_cond','bridge_well_conditioned_pass','response_opnorm','response_spectral_radius','transfer_opnorm_raw_before_post_bound','lifted_transfer_opnorm','QP_transfer_mean_resid','projected_QP_eig_imag_max','projected_QP_normality_resid','projected_QP_opnorm','projected_QP_cube_identity_resid','projected_QP_minus_identity_resid','halfplane_bias_pos_minus_neg','artifact_suspect','bounded_structural_transfer_pass','complex_like_stable_pass','QP_lock_stable_pass','C3_like_stable_pass','minus_identity_stable_pass']
    return {k: r.get(k) for k in keys}


def summarize(rows: List[dict]) -> dict:
    ok = [r for r in rows if r.get('transfer_status') == 'ok']
    good = [r for r in ok if fbool(r.get('bridge_well_conditioned_pass')) and fbool(r.get('response_bounded_mode')) and not fbool(r.get('artifact_suspect'))]
    def count(rows_, k): return sum(1 for r in rows_ if fbool(r.get(k)))
    def finite_vals(rows_, k):
        vals=[]
        for r in rows_:
            try:
                v=float(r.get(k, math.nan))
                if np.isfinite(v): vals.append(v)
            except Exception:
                pass
        return vals
    def mn(rows_, k, default=0.0):
        vals=finite_vals(rows_, k); return float(min(vals)) if vals else default
    def mx(rows_, k, default=0.0):
        vals=finite_vals(rows_, k); return float(max(vals)) if vals else default
    best_lock = min(good, key=lambda r: float(r.get('QP_transfer_mean_resid', 99.0)), default=None)
    best_complex = max(good, key=lambda r: float(r.get('projected_QP_eig_imag_max', 0.0)) if float(r.get('projected_QP_normality_resid', 99.0)) <= 0.75 else -1, default=None)
    best_struct = next((r for r in good if fbool(r.get('bounded_structural_transfer_pass'))), None)
    return {
        'transfer_row_count': len(rows),
        'ok_transfer_row_count': len(ok),
        'well_conditioned_bounded_nonartifact_count': len(good),
        'ill_conditioned_or_artifact_count': len(ok)-len(good),
        'bounded_structural_transfer_pass_count': count(good, 'bounded_structural_transfer_pass'),
        'QP_lock_stable_pass_count': count(good, 'QP_lock_stable_pass'),
        'complex_like_stable_pass_count': count(good, 'complex_like_stable_pass'),
        'C3_like_stable_pass_count': count(good, 'C3_like_stable_pass'),
        'minus_identity_stable_pass_count': count(good, 'minus_identity_stable_pass'),
        'bounded_contraction_audit_pass_count': count(good, 'bounded_contraction_audit_pass'),
        'best_good_QP_transfer_mean_resid': mn(good, 'QP_transfer_mean_resid'),
        'best_good_QP_transfer_max_resid': mn(good, 'QP_transfer_max_resid'),
        'best_good_projected_J2_plus_I_resid': mn(good, 'projected_J2_plus_I_resid'),
        'max_good_projected_QP_eig_imag': mx(good, 'projected_QP_eig_imag_max'),
        'best_good_C3_cube_identity_resid': mn(good, 'projected_QP_cube_identity_resid'),
        'best_good_minus_identity_resid': mn(good, 'projected_QP_minus_identity_resid'),
        'max_good_halfplane_abs_bias': max([abs(v) for v in finite_vals(good,'halfplane_bias_pos_minus_neg')] or [0.0]),
        'best_lock_good_row': slim_row(best_lock),
        'best_complex_good_row': slim_row(best_complex),
        'first_structural_pass_row': slim_row(best_struct),
    }


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
    rows: List[dict] = []
    for a in assembly_log:
        if fbool(a.get('assembly_applied')):
            rows.extend(bounded_lift_rows(model, K, a, args))
    summ = summarize(rows)
    write_csv(vout / 'birth_geometry_log.csv', birth_log)
    write_csv(vout / 'assembly_pairing_log.csv', pairing_log)
    write_csv(vout / 'assembly_log.csv', assembly_log)
    write_csv(vout / 'candidate_eval_rows.csv', candidate_rows)
    write_csv(vout / 'bounded_response_bridge_rows.csv', rows)
    write_csv(vout / 'directed_pair_rows.csv', pair_rows)
    write_csv(vout / 'signed_quadrature_rows.csv', signed_rows)
    decision_used_delta = any(fbool(x.get('decision_used_delta_beta')) for x in pairing_log + assembly_log + candidate_rows)
    res = {
        'variant': variant,
        'baseline_metrics': baseline_metrics,
        'auto_metrics': auto_metrics,
        'directed_metrics': dm,
        'signed_metrics': sm,
        'alignment_metrics': am,
        'assemblies_applied': sum(1 for x in assembly_log if fbool(x.get('assembly_applied'))),
        'pairings_applied': sum(1 for x in pairing_log if fbool(x.get('applied'))),
        'bounded_bridge_summary': summ,
        'decision_used_delta_beta_any': decision_used_delta,
    }
    (vout / 'summary.json').write_text(json.dumps(res, indent=2), encoding='utf-8')
    return res


def table_rows(summaries: List[dict]) -> str:
    lines = [
        '| variant | beta | assemblies | pair_harm | Q_harm | P_harm | base J-lock | good bounded rows | structural pass | QP pass | complex pass | C3 pass | -I pass | best good QP resid | max good ImEig | best C3 resid | max halfplane bias | used Δβ? |',
        '|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|'
    ]
    for r in summaries:
        a = r['auto_metrics']; dm = r['directed_metrics']; am = r['alignment_metrics']; bs = r['bounded_bridge_summary']
        lines.append(f"| {r['variant']} | ({a['beta0']},{a['beta1']},{a['beta2']},{a['beta3']}) | {r['assemblies_applied']} | {dm['pair_transport_harmonic_ratio']:.6g} | {am['Q_even_harmonic_ratio']:.6g} | {am['P_odd_harmonic_ratio']:.6g} | {am['best_per_pair_mean_J_lock_resid']:.6g} | {bs['well_conditioned_bounded_nonartifact_count']} | {bs['bounded_structural_transfer_pass_count']} | {bs['QP_lock_stable_pass_count']} | {bs['complex_like_stable_pass_count']} | {bs['C3_like_stable_pass_count']} | {bs['minus_identity_stable_pass_count']} | {bs['best_good_QP_transfer_mean_resid']:.6g} | {bs['max_good_projected_QP_eig_imag']:.6g} | {bs['best_good_C3_cube_identity_resid']:.6g} | {bs['max_good_halfplane_abs_bias']:.6g} | {r['decision_used_delta_beta_any']} |")
    return '\n'.join(lines)


def docs(summaries: List[dict]) -> tuple[str,str,str,str]:
    table = table_rows(summaries)
    smd = f"""# SUMMARY — bounded response bridge sanity gate

This package follows the response-monodromy-to-Q/P transfer package, but it does **not** treat every complex-like projected eigenvalue as evidence.  It first asks whether the response→Q/P bridge is numerically and structurally sane:

- bridge singular values and condition number are audited;
- ill-conditioned or pseudo-inverse-sensitive rows are separated;
- response operators are bounded by operator norm / spectral radius without fitting a target structure;
- Q/P-lock, complex-like, C3-like, and minus-identity-like signatures are only counted on well-conditioned bounded non-artifact rows;
- contraction/leakage/half-plane bias are logged as possible semigroup/boundary-value traces, not as complex structure.

{table}
"""
    rmd = f"""# RESULTS — bounded response bridge sanity gate

{table}

## Interpretation rule

A positive result requires a bounded, well-conditioned, non-artifact bridge row that carries one of the audited structures:

- stable Q/P-lock transfer,
- stable complex-/rotation-like projected transfer,
- stable C3-like cyclic transfer,
- stable minus-identity-like transfer.

Rows dominated by poor bridge conditioning or large unbounded transfer norm are not counted as structural evidence.  This is directly motivated by the warning that half-plane / boundary-value / damping conventions can hide inside analytic continuation choices; here we therefore distinguish real semigroup/contraction traces from mere pseudo-inverse activity.

## Compact JSON

```json
{json.dumps(summaries, indent=2)[:28000]}
```
"""
    audit = """# SOURCE AUDIT

This package is derived from `test_response_monodromy_to_QP_transfer_gate.py` and intentionally tightens it.

It does not:

- set `i`, global `J`, Hodge, positivity, physical adjoint, or C*-norm;
- use delta-beta/H2 as a selection input;
- infer spin/double-cover structure from powers of one operator;
- count ill-conditioned pseudo-inverse rows as evidence.

It does:

- use the same dynamic growth / A-B assembly path as the prior packages;
- build response operators from directed sibling fan data;
- bridge to Q/P motif carriers by incidence and birth-frame projections;
- audit bridge singular values, condition numbers, boundedness, normality, contraction, leakage, Q/P-lock, complex eigen activity, C3 closure, and minus-identity proximity.

Model label: CQNM/s=-1 inspired growing real complement-network diagnostic, not an SG/static geometry proof.  The script remains deterministic and provenance-only in the sense used by the package series.
"""
    readme = """# Bounded response bridge sanity gate

Run:

```bash
python3 test_bounded_response_bridge_sanity_gate.py
```

Outputs are written to `bounded_response_bridge_sanity_out_L2/`.
"""
    return smd, rmd, audit, readme


def package(out: Path, zip_path: Path) -> None:
    files = [
        Path(__file__).name,
        'test_response_monodromy_to_QP_transfer_gate.py',
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
    if zip_path.exists(): zip_path.unlink()
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
    ap.add_argument('--min-bridge-rank', type=int, default=2)
    ap.add_argument('--min-bridge-sv', type=float, default=1e-10)
    ap.add_argument('--max-bridge-cond', type=float, default=1e4)
    ap.add_argument('--artifact-cond-threshold', type=float, default=1e5)
    ap.add_argument('--artifact-raw-opnorm-threshold', type=float, default=1e4)
    ap.add_argument('--include-transfer-post-bounded', action='store_true', default=True)
    ap.add_argument('--min-transfer-norm', type=float, default=1e-8)
    ap.add_argument('--min-response-circulation', type=float, default=1e-6)
    ap.add_argument('--bounded-opnorm-threshold', type=float, default=1.05)
    ap.add_argument('--normality-threshold', type=float, default=0.35)
    ap.add_argument('--qp-lock-threshold', type=float, default=0.25)
    ap.add_argument('--qp-lock-max-threshold', type=float, default=0.35)
    ap.add_argument('--complex-imag-threshold', type=float, default=0.05)
    ap.add_argument('--c3-identity-threshold', type=float, default=0.35)
    ap.add_argument('--not-identity-threshold', type=float, default=0.35)
    ap.add_argument('--minus-identity-threshold', type=float, default=0.35)
    ap.add_argument('--real-halfplane-tol', type=float, default=1e-9)
    ap.add_argument('--out', default='bounded_response_bridge_sanity_out_L2')
    ap.add_argument('--zip', default='cnna_bounded_response_bridge_sanity_gate_pkg_L2.zip')
    args = ap.parse_args()

    p70.patch_modules()
    out = Path(args.out)
    if out.exists(): shutil.rmtree(out)
    out.mkdir(parents=True)
    summaries = [run_variant(v, args, out) for v in args.variants]
    smd,rmd,audit,readme = docs(summaries)
    (out/'SUMMARY.md').write_text(smd, encoding='utf-8')
    (out/'RESULTS.md').write_text(rmd, encoding='utf-8')
    (out/'SOURCE_AUDIT.md').write_text(audit, encoding='utf-8')
    (out/'README.md').write_text(readme, encoding='utf-8')
    write_csv(out/'comparative_bounded_response_bridge_sanity_summary.csv', [
        {
            'variant': r['variant'],
            'beta0': r['auto_metrics']['beta0'], 'beta1': r['auto_metrics']['beta1'], 'beta2': r['auto_metrics']['beta2'], 'beta3': r['auto_metrics']['beta3'],
            'assemblies_applied': r['assemblies_applied'],
            'pair_harm': r['directed_metrics']['pair_transport_harmonic_ratio'],
            'Q_harm': r['alignment_metrics']['Q_even_harmonic_ratio'],
            'P_harm': r['alignment_metrics']['P_odd_harmonic_ratio'],
            **r['bounded_bridge_summary'],
            'decision_used_delta_beta_any': r['decision_used_delta_beta_any'],
        } for r in summaries
    ])
    (out/'summary.json').write_text(json.dumps(summaries, indent=2), encoding='utf-8')
    package(out, Path(args.zip))
    print(f'wrote {args.zip}')

if __name__ == '__main__':
    main()
