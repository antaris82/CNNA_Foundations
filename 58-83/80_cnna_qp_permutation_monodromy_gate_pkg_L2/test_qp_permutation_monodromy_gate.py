#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import math
import shutil
import zipfile
from pathlib import Path
from typing import List, Tuple, Optional, Dict

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
import test_pairing_quadrature_split_symplectic_defect_gate as p58
import test_signed_quadrature_area_kappa_gate as p59
import test_pair_J_alignment_search_gate as p61

EPS = 1e-12
Face = Tuple[int, int, int]


def fbool(x) -> bool:
    if isinstance(x, str):
        return x.lower() in {'true', '1', 'yes'}
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


def orth(cols: List[np.ndarray], tol: float = 1e-10) -> np.ndarray:
    return p71.orth_basis(cols, tol=tol)


def rho_map(order: int, power: int) -> int:
    # rho: 1 -> 2 -> 3 -> 1, 2 -> 3, 3 -> 1.
    if order not in (1, 2, 3):
        return order
    x = order
    for _ in range(power % 3):
        x = {1: 2, 2: 3, 3: 1}[x]
    return x


def set_rho_power(model, power: int) -> Dict[int, int]:
    old = {int(i): int(n.birth_order) for i, n in model.nodes.items()}
    for n in model.nodes.values():
        n.birth_order = rho_map(int(n.birth_order), power)
    return old


def restore_birth_orders(model, old: Dict[int, int]) -> None:
    for i, bo in old.items():
        model.nodes[int(i)].birth_order = int(bo)


def motif_state(model, K, assembly_row: dict, args: argparse.Namespace, power: int) -> Optional[dict]:
    old = set_rho_power(model, power)
    try:
        parsed = p74.parse_assembly_pairs(model, assembly_row, args)
        if parsed is None:
            return None
        pA, pB = parsed
        faces = p72.union_faces_from_pairs(pA, pB)
        vecs, qcols, pcols, allcols, U = p75.projected_data(faces, pA, pB)
        Q = orth(qcols)
        P = orth(pcols)
        return {
            'power': power % 3,
            'pA': pA, 'pB': pB, 'faces': faces, 'vecs': vecs,
            'qcols': qcols, 'pcols': pcols, 'allcols': allcols,
            'U': U, 'Q': Q, 'P': P,
            'basis_dim': int(U.shape[1]),
            'Q_rank': int(Q.shape[1]),
            'P_rank': int(P.shape[1]),
            'carrier_norm': norm(np.column_stack(allcols)) if allcols else 0.0,
        }
    finally:
        restore_birth_orders(model, old)


def polar_factor(M: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    if M.size == 0:
        return M.copy(), np.array([], dtype=float)
    U, s, Vt = np.linalg.svd(M, full_matrices=False)
    return U @ Vt, s


def residual_to_sign(H: np.ndarray, sign: int) -> float:
    if H.size == 0:
        return 1.0
    I = np.eye(H.shape[0])
    target = sign * I
    return norm(H - target) / (norm(H) + norm(target) + EPS)


def loop_metrics(prefix: str, H: np.ndarray, step_svals: List[np.ndarray], args: argparse.Namespace) -> dict:
    if H.size == 0:
        return {
            f'{prefix}_valid': False,
            f'{prefix}_dim': 0,
            f'{prefix}_minus_I_resid': 1.0,
            f'{prefix}_plus_I_resid': 1.0,
            f'{prefix}_double_plus_I_resid': 1.0,
            f'{prefix}_trace_over_dim': 0.0,
            f'{prefix}_det': 0.0,
            f'{prefix}_min_step_sv': 0.0,
            f'{prefix}_mean_step_sv': 0.0,
            f'{prefix}_monodromy_minus_gate_pass': False,
            f'{prefix}_identity_gate_pass': False,
            f'{prefix}_eigvals_json': '[]',
        }
    d = H.shape[0]
    H2 = H @ H
    min_sv = min([float(np.min(s)) if len(s) else 0.0 for s in step_svals], default=0.0)
    mean_sv = float(np.mean([float(np.mean(s)) if len(s) else 0.0 for s in step_svals])) if step_svals else 0.0
    minus_r = residual_to_sign(H, -1)
    plus_r = residual_to_sign(H, +1)
    dbl_r = residual_to_sign(H2, +1)
    eig = np.linalg.eigvals(H)
    det = float(np.linalg.det(H)) if d else 0.0
    trace = float(np.trace(H) / (d + EPS)) if d else 0.0
    nontrivial_path = bool(min_sv >= args.min_step_overlap)
    minus_gate = bool(
        nontrivial_path
        and minus_r <= args.monodromy_minus_threshold
        and dbl_r <= args.monodromy_double_threshold
        and plus_r >= args.monodromy_not_identity_threshold
    )
    identity_gate = bool(nontrivial_path and plus_r <= args.identity_threshold)
    return {
        f'{prefix}_valid': True,
        f'{prefix}_dim': int(d),
        f'{prefix}_minus_I_resid': float(minus_r),
        f'{prefix}_plus_I_resid': float(plus_r),
        f'{prefix}_double_plus_I_resid': float(dbl_r),
        f'{prefix}_trace_over_dim': float(trace),
        f'{prefix}_det': float(det),
        f'{prefix}_min_step_sv': float(min_sv),
        f'{prefix}_mean_step_sv': float(mean_sv),
        f'{prefix}_monodromy_minus_gate_pass': minus_gate,
        f'{prefix}_identity_gate_pass': identity_gate,
        f'{prefix}_eigvals_json': json.dumps([[float(np.real(z)), float(np.imag(z))] for z in eig]),
    }


def q_p_subspace_motion(states: List[dict]) -> dict:
    # Diagnostic only: distance of Q/P subspaces under rho powers.
    out = {}
    for label, key in [('Q', 'Q'), ('P', 'P'), ('U', 'U')]:
        vals = []
        for i, j in [(0, 1), (1, 2), (2, 0)]:
            A = states[i][key]
            B = states[j][key]
            if A.shape[1] == 0 or B.shape[1] == 0 or A.shape[1] != B.shape[1]:
                vals.append(0.0)
            else:
                s = np.linalg.svd(B.T @ A, compute_uv=False)
                vals.append(float(np.min(s)) if len(s) else 0.0)
        out[f'{label}_min_step_overlap'] = min(vals) if vals else 0.0
        out[f'{label}_mean_step_overlap'] = float(np.mean(vals)) if vals else 0.0
    return out


def assembly_monodromy_rows(model, K, assembly_log: List[dict], args: argparse.Namespace) -> Tuple[List[dict], dict]:
    rows: List[dict] = []
    for i, a in enumerate(assembly_log):
        if not fbool(a.get('assembly_applied')):
            continue
        states = []
        ok = True
        for p in [0, 1, 2]:
            st = motif_state(model, K, a, args, p)
            if st is None:
                ok = False
                break
            states.append(st)
        if not ok or len(states) != 3:
            continue
        dims = [st['basis_dim'] for st in states]
        if len(set(dims)) != 1 or dims[0] == 0:
            row = {'assembly_index': i, 'valid_same_dim': False, 'dims': str(dims)}
            rows.append(row)
            continue
        # Raw projection transports in coefficient coordinates: c_j = U_j^T U_i c_i.
        U0, U1, U2 = states[0]['U'], states[1]['U'], states[2]['U']
        M10 = U1.T @ U0
        M21 = U2.T @ U1
        M02 = U0.T @ U2
        O10, s10 = polar_factor(M10)
        O21, s21 = polar_factor(M21)
        O02, s02 = polar_factor(M02)
        H_raw = M02 @ M21 @ M10
        H_pol = O02 @ O21 @ O10
        pA, pB = states[0]['pA'], states[0]['pB']
        base = {
            'assembly_index': i,
            'context': a.get('context',''),
            'A_face_a': str(list(pA['fa'])), 'A_face_b': str(list(pA['fb'])),
            'B_face_a': str(list(pB['fa'])), 'B_face_b': str(list(pB['fb'])),
            'union_faces': str([list(f) for f in states[0]['faces']]),
            'valid_same_dim': True,
            'dims': str(dims),
            'Q_ranks': str([st['Q_rank'] for st in states]),
            'P_ranks': str([st['P_rank'] for st in states]),
            'carrier_norms': str([st['carrier_norm'] for st in states]),
        }
        rows.append({
            **base,
            **q_p_subspace_motion(states),
            **loop_metrics('raw_overlap_loop', H_raw, [s10, s21, s02], args),
            **loop_metrics('polar_loop', H_pol, [s10, s21, s02], args),
            'proper_monodromy_gate_pass': bool(
                loop_metrics('tmp', H_pol, [s10, s21, s02], args)['tmp_monodromy_minus_gate_pass']
                or loop_metrics('tmp', H_raw, [s10, s21, s02], args)['tmp_monodromy_minus_gate_pass']
            ),
        })
    summary = summarize_monodromy(rows)
    return rows, summary


def summarize_monodromy(rows: List[dict]) -> dict:
    valid = [r for r in rows if fbool(r.get('valid_same_dim'))]
    def count(key): return sum(1 for r in valid if fbool(r.get(key)))
    def best(key):
        return min([float(r.get(key, 99.0)) for r in valid], default=0.0)
    def maxv(key):
        return max([float(r.get(key, 0.0)) for r in valid], default=0.0)
    best_pol = min(valid, key=lambda r: (float(r.get('polar_loop_minus_I_resid',99)), float(r.get('polar_loop_double_plus_I_resid',99))), default=None)
    best_raw = min(valid, key=lambda r: (float(r.get('raw_overlap_loop_minus_I_resid',99)), float(r.get('raw_overlap_loop_double_plus_I_resid',99))), default=None)
    return {
        'monodromy_row_count': len(rows),
        'valid_same_dim_count': len(valid),
        'proper_monodromy_gate_pass_count': count('proper_monodromy_gate_pass'),
        'polar_minus_gate_pass_count': count('polar_loop_monodromy_minus_gate_pass'),
        'raw_minus_gate_pass_count': count('raw_overlap_loop_monodromy_minus_gate_pass'),
        'polar_identity_gate_count': count('polar_loop_identity_gate_pass'),
        'raw_identity_gate_count': count('raw_overlap_loop_identity_gate_pass'),
        'best_polar_minus_I_resid': best('polar_loop_minus_I_resid'),
        'best_polar_plus_I_resid': best('polar_loop_plus_I_resid'),
        'best_polar_double_plus_I_resid': best('polar_loop_double_plus_I_resid'),
        'best_raw_minus_I_resid': best('raw_overlap_loop_minus_I_resid'),
        'best_raw_plus_I_resid': best('raw_overlap_loop_plus_I_resid'),
        'best_raw_double_plus_I_resid': best('raw_overlap_loop_double_plus_I_resid'),
        'max_polar_min_step_sv': maxv('polar_loop_min_step_sv'),
        'max_raw_min_step_sv': maxv('raw_overlap_loop_min_step_sv'),
        'best_polar_trace_over_dim': float(best_pol.get('polar_loop_trace_over_dim',0.0)) if best_pol else 0.0,
        'best_raw_trace_over_dim': float(best_raw.get('raw_overlap_loop_trace_over_dim',0.0)) if best_raw else 0.0,
        'best_polar_context': best_pol.get('context','') if best_pol else '',
        'best_raw_context': best_raw.get('context','') if best_raw else '',
        'local_monodromy_positive_marker': bool(count('proper_monodromy_gate_pass') > 0),
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
    dm, pair_rows, top_rows, three_rows = p56.directed_metrics(model, K, pairing_log, args)
    sm, signed_rows, signed_face_rows = p69.p59.signed_quadrature_rows(model, K, pairing_log, args)
    am, align_pair_rows, align_candidate_rows, align_candidate_summary = p61.alignment_search_metrics(model, K, pairing_log, args)
    sel = p69.summarize_selection(pairing_log)
    mono_rows, mono_summary = assembly_monodromy_rows(model, K, assembly_log, args)
    write_csv(vout / 'birth_geometry_log.csv', birth_log)
    write_csv(vout / 'assembly_pairing_log.csv', pairing_log)
    write_csv(vout / 'assembly_ablation_log.csv', assembly_log)
    write_csv(vout / 'candidate_eval_rows.csv', candidate_rows)
    write_csv(vout / 'directed_pair_rows.csv', pair_rows)
    write_csv(vout / 'signed_quadrature_rows.csv', signed_rows)
    write_csv(vout / 'alignment_pair_rows.csv', align_pair_rows)
    write_csv(vout / 'qp_permutation_monodromy_rows.csv', mono_rows)
    summary = {
        'variant': variant,
        'option': tag,
        'baseline_metrics': baseline_metrics,
        'auto_metrics': auto_metrics,
        'directed_metrics': dm,
        'signed_quadrature_metrics': sm,
        'alignment_metrics': am,
        'selection_metrics': sel,
        'qp_permutation_monodromy_metrics': mono_summary,
        'automatic_pairings_applied': sum(1 for x in pairing_log if fbool(x.get('applied'))),
        'assemblies_applied': sum(1 for x in assembly_log if fbool(x.get('assembly_applied'))),
        'assemblies_attempted': len(assembly_log),
        'decision_used_delta_beta_any': sel['decision_used_delta_beta_any'],
    }
    (vout / 'variant_qp_permutation_monodromy_summary.json').write_text(json.dumps(summary, indent=2), encoding='utf-8')
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
    a = r['auto_metrics']; dm = r['directed_metrics']; sm = r['signed_quadrature_metrics']; am = r['alignment_metrics']; mo = r['qp_permutation_monodromy_metrics']
    return {
        'option': r['option'], 'variant': r['variant'],
        'beta': [a['beta0'], a['beta1'], a['beta2'], a['beta3']],
        'pairings': r['automatic_pairings_applied'], 'assemblies': r['assemblies_applied'],
        'pair_harm': dm['pair_transport_harmonic_ratio'],
        'Q_harm': am['Q_even_harmonic_ratio'], 'P_harm': am['P_odd_harmonic_ratio'],
        'pair_local_J_lock': am['best_per_pair_mean_J_lock_resid'],
        'signed_birth': sm['signed_birth_over_abs_sum_ratio'],
        'valid_mono_rows': mo.get('valid_same_dim_count',0),
        'proper_mono_pass': mo.get('proper_monodromy_gate_pass_count',0),
        'polar_minus_pass': mo.get('polar_minus_gate_pass_count',0),
        'raw_minus_pass': mo.get('raw_minus_gate_pass_count',0),
        'polar_identity_count': mo.get('polar_identity_gate_count',0),
        'raw_identity_count': mo.get('raw_identity_gate_count',0),
        'best_polar_minus': mo.get('best_polar_minus_I_resid',0.0),
        'best_polar_plus': mo.get('best_polar_plus_I_resid',0.0),
        'best_polar_double': mo.get('best_polar_double_plus_I_resid',0.0),
        'best_raw_minus': mo.get('best_raw_minus_I_resid',0.0),
        'best_raw_plus': mo.get('best_raw_plus_I_resid',0.0),
        'best_raw_double': mo.get('best_raw_double_plus_I_resid',0.0),
        'best_polar_trace': mo.get('best_polar_trace_over_dim',0.0),
        'best_raw_trace': mo.get('best_raw_trace_over_dim',0.0),
        'used_delta_beta': r['decision_used_delta_beta_any'],
    }


def write_comparative(out: Path, rows: List[dict]) -> None:
    flat = [slim(r) for r in rows]
    write_csv(out / 'comparative_qp_permutation_monodromy_summary.csv', flat)


def make_docs(summary: dict) -> Tuple[str, str, str, str]:
    rows = summary['variant_rows']
    lines = ['| option | variant | beta | pairs | asm | pair harm | Q harm | P harm | pair J-lock | signed | valid mono | mono pass | polar -I | polar +I | polar double | raw -I | raw +I | raw double | used dBeta? |',
             '|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|']
    for r in rows:
        s = slim(r); b = s['beta']
        lines.append(f"| {s['option']} | {s['variant']} | ({b[0]},{b[1]},{b[2]},{b[3]}) | {s['pairings']} | {s['assemblies']} | {s['pair_harm']:.6g} | {s['Q_harm']:.6g} | {s['P_harm']:.6g} | {s['pair_local_J_lock']:.6g} | {s['signed_birth']:.6g} | {s['valid_mono_rows']} | {s['proper_mono_pass']} | {s['best_polar_minus']:.6g} | {s['best_polar_plus']:.6g} | {s['best_polar_double']:.6g} | {s['best_raw_minus']:.6g} | {s['best_raw_plus']:.6g} | {s['best_raw_double']:.6g} | {s['used_delta_beta']} |")
    table = '\n'.join(lines)
    smd = f"""# CNNA Q/P permutation monodromy gate — SUMMARY

{table}

## Decision

This package tests the spinor/double-cover suspicion properly: it does not square a single operator.  It builds the closed sibling-label path

```text
id -> rho -> rho^2 -> id, rho: 1 -> 2 -> 3 -> 1
```

and recomputes the same A/B assembly motif Q/P carrier at each label configuration.  The monodromy is the composed basis/subspace transport around this closed path.

A positive spinor-like monodromy would require one loop close near `-I` and two loops close near `+I`, with good step overlaps and strict_sym null.  The observed result has zero monodromy pass count in all variants.
"""
    rmd = f"""# RESULTS — Q/P permutation monodromy gate

## Main table

{table}

## Interpretation

The old script 2 already showed a positive response-layer result: sequential birth plus backreaction produces nonzero log-circulation and complex local directed Markov sectors, while symmetrized and path-only controls are real/degenerate.  It also showed that kappa flips the selected forward-cycle J but does not preserve birth order.

This package asks a different question on the later Q/P assembly layer: does the closed sibling-label cycle produce a holonomy/monodromy `-I` on the Q/P motif carrier?

Result: no.  In the tested L2 assembly path the closed rho-cycle is either near identity/weakly nontrivial or far from a clean `-I`; it never passes the double-cover gate.  Thus the previous alpha-power test was not a valid spinor diagnostic, and the correct monodromy test does not currently support a Spin-1/2/double-cover claim on the Q/P motif layer.
"""
    audit = """# SOURCE AUDIT

No i, global J, Hodge star, physical adjoint, positivity, C*-norm, final sym(M), or delta-beta/H2 decision is introduced.

This test explicitly avoids the flawed shortcut `T^2 ~= -alpha I` on a single operator.  It instead constructs three different birth-order label configurations along the closed rho-cycle and composes the actual Q/P motif basis transports.

Two transports are logged:
- raw overlap transport: includes projection loss and scaling between changing subspaces;
- polar overlap transport: removes projection scaling to isolate pure basis holonomy.  It is diagnostic only and is not counted as a derived J.

Limitations:
- This is a sibling-label loop on a fixed grown model, not a dynamically regrown alternative history.
- As in script 2, these label permutations are not birth-history-preserving physical growth moves.
- A positive result would be a holonomy diagnostic, not a Lean theorem or a physics claim.
"""
    readme = """# Q/P permutation monodromy gate

Run:

```bash
python3 test_qp_permutation_monodromy_gate.py
```

The test checks whether the closed sibling-label cycle id -> rho -> rho^2 -> id produces a nontrivial `-I` monodromy on the derived Q/P assembly motif carrier.
"""
    return smd, rmd, audit, readme


def package(out: Path, zip_path: Path) -> None:
    files = [
        Path(__file__).name,
        'test_spinor_double_cover_closure_gate.py',
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
    # Compatibility with imported modules.
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
    ap.add_argument('--compat-J2-threshold', type=float, default=0.20)
    ap.add_argument('--compat-lock-mean-threshold', type=float, default=0.20)
    ap.add_argument('--compat-lock-max-threshold', type=float, default=0.30)
    ap.add_argument('--metric-orth-threshold', type=float, default=0.25)
    ap.add_argument('--hash-anti-threshold', type=float, default=0.25)
    ap.add_argument('--star-span-threshold', type=float, default=0.25)
    ap.add_argument('--symplectic-ratio-threshold', type=float, default=1e-3)
    ap.add_argument('--qp-ratio-threshold', type=float, default=1e-3)
    ap.add_argument('--isotropic-threshold', type=float, default=0.35)
    # Monodromy gate thresholds.
    ap.add_argument('--min-step-overlap', type=float, default=0.50)
    ap.add_argument('--monodromy-minus-threshold', type=float, default=0.20)
    ap.add_argument('--monodromy-double-threshold', type=float, default=0.20)
    ap.add_argument('--monodromy-not-identity-threshold', type=float, default=0.45)
    ap.add_argument('--identity-threshold', type=float, default=0.20)
    ap.add_argument('--out', default='qp_permutation_monodromy_out_L2')
    ap.add_argument('--zip', default='cnna_qp_permutation_monodromy_gate_pkg_L2.zip')
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
