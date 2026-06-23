#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import shutil
import zipfile
from pathlib import Path
from typing import Dict, List, Tuple, Optional

import numpy as np

import cnna_non_shelling_core as core
import test_nonlinear_asymmetry_cascade_growth as nl
from test_growth_with_asymmetry_gated_complement_pairing import ordinary_outward_step
import test_harmonic_k_orientation_kappa_gate as hk
import test_pairing_transport_antisym_birth_coherence_gate as p56
import test_pairing_quadrature_adjoint_pairing_gate as p60
import test_pair_J_alignment_search_gate as p61

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


def parse_candidate_faces(row: dict) -> Tuple[Optional[Face], Optional[Face]]:
    try:
        fa = tuple(sorted(core.parse_face_string(str(row.get('face_a', '')))))
        fb_txt = str(row.get('face_b', '')).split('perm=')[0].strip()
        fb = tuple(sorted(core.parse_face_string(fb_txt)))
        if len(fa) == 3 and len(fb) == 3:
            return fa, fb
    except Exception:
        pass
    return None, None


def candidate_alignment(model, K: core.SimplicialComplex, row: dict, args: argparse.Namespace) -> dict:
    fa, fb = parse_candidate_faces(row)
    if fa is None or fb is None:
        return {'candidate_eval_status': 'bad_face_parse'}
    faces = set(K.faces())
    if fa not in faces or fb not in faces:
        return {'candidate_eval_status': 'faces_not_in_current_K'}
    try:
        ka = p56.axial(p56.face_K_directed(model, fa, args.source, args.phase_sign, args.antisym_eta, args.erase_phase_for_strict_sym))
        kb = p56.axial(p56.face_K_directed(model, fb, args.source, args.phase_sign, args.antisym_eta, args.erase_phase_for_strict_sym))
        na_out = hk.face_normal(model, fa, 'outward')
        nb_out = hk.face_normal(model, fb, 'outward')
        na_birth = hk.face_normal(model, fa, 'birth_order')
        R = p56.rotation_from_to(nb_out, -na_out)
        J = p60.block_J(R)
        C = p60.block_C(R)
        I = np.eye(6)
        kb_r = R @ kb
        Q_a = ka + kb_r
        P_a = ka - kb_r
        Q_b = -(R.T @ Q_a)
        P_b = +(R.T @ P_a)
        q = p61.block(Q_a, Q_b)
        p = p61.block(P_a, P_b)
        raw = p61.block(ka, kb)
        transported_raw = p61.block(ka, kb_r)
        E = 0.5 * (I + C)
        O = 0.5 * (I - C)
        even_raw = E @ raw
        odd_raw = O @ raw
        even_trans = E @ transported_raw
        odd_trans = O @ transported_raw
        candidates = [
            ('raw_QP', q, p),
            ('unit_QP', p61.safe_unit(q), p61.safe_unit(p)),
            ('C_eigen_raw_even_odd', even_raw, odd_raw),
            ('C_eigen_transported_even_odd', even_trans, odd_trans),
            ('unit_C_eigen_raw_even_odd', p61.safe_unit(even_raw), p61.safe_unit(odd_raw)),
            ('unit_C_eigen_transported_even_odd', p61.safe_unit(even_trans), p61.safe_unit(odd_trans)),
        ]
        res_rows = []
        for name, qv, pv in candidates:
            res = p61.J_lock_residual(J, qv, pv)
            res_rows.append((name, res))
        best_name, best = min(res_rows, key=lambda z: (z[1]['J_lock_mean_resid'], z[1]['J_lock_max_resid']))
        # Focus metric requested by the gate: best C-eigen only, separate from raw_QP.
        ceig = [z for z in res_rows if z[0].startswith('C_eigen') or z[0].startswith('unit_C_eigen')]
        best_c_name, best_c = min(ceig, key=lambda z: (z[1]['J_lock_mean_resid'], z[1]['J_lock_max_resid']))
        cross = np.cross(Q_a, P_a)
        comm_abs = norm(cross)
        comm_signed_birth = float(np.dot(cross, na_birth))
        transport_cos = float(np.dot(ka, kb_r) / ((norm(ka) * norm(kb_r)) + EPS))
        return {
            'candidate_eval_status': 'ok',
            'face_a_tuple': str(list(fa)),
            'face_b_tuple': str(list(fb)),
            'K_axial_a_norm': norm(ka),
            'K_axial_b_norm': norm(kb),
            'transport_cosine_ka_kb_reversed': transport_cos,
            'Q_norm': norm(q),
            'P_norm': norm(p),
            'comm_abs_area': comm_abs,
            'comm_signed_birth': comm_signed_birth,
            'comm_signed_birth_over_abs': comm_signed_birth / (comm_abs + EPS),
            'best_any_candidate': best_name,
            'best_any_J_lock_mean_resid': best['J_lock_mean_resid'],
            'best_any_J_lock_max_resid': best['J_lock_max_resid'],
            'best_C_eigen_candidate': best_c_name,
            'best_C_eigen_J_lock_mean_resid': best_c['J_lock_mean_resid'],
            'best_C_eigen_J_lock_max_resid': best_c['J_lock_max_resid'],
            'raw_QP_J_lock_mean_resid': dict(res_rows)['raw_QP']['J_lock_mean_resid'],
            'unit_QP_J_lock_mean_resid': dict(res_rows)['unit_QP']['J_lock_mean_resid'],
        }
    except Exception as e:
        return {'candidate_eval_status': f'error:{type(e).__name__}:{e}'}


def candidate_base_row(row: dict) -> dict:
    keys = [
        'candidate_id','move_class','status','face_a','face_b','A_gate','A_invariant','A_rank_score',
        'directed_imbalance','transverse_complementarity','response_score','response_rank_legal',
        'delta_beta1','delta_beta2','delta_beta3','delta_boundary_faces','new_beta1','new_beta2',
        'A_nonreciprocal_norm','A_live_record_gap','A_shell_gap','A_fan_directed_mean',
    ]
    return {k: row.get(k, '') for k in keys}


def should_eval_candidate(row: dict, args: argparse.Namespace) -> bool:
    if row.get('status') != 'ok':
        return False
    if row.get('move_class') not in {'handle_candidate', 'quotient_candidate'}:
        return False
    if args.only_A_gate and not bool(row.get('A_gate')):
        return False
    return True


def scan_evaluate_and_apply(model, args: argparse.Namespace, variant: str, out: Path):
    K = core.SimplicialComplex(f'{variant}_C_eigen_refinement')
    root_seeded = False
    birth_log: List[dict] = []
    pairing_log: List[dict] = []
    candidate_rows: List[dict] = []
    scan_summaries: List[dict] = []
    used_faces: set[Face] = set()
    global_pair_count = 0
    scan_id = 0

    for ev in sorted(model.birth_events, key=lambda x: int(x['t'])):
        root_seeded, added, encoded = ordinary_outward_step(model, K, ev, root_seeded)
        event_t = int(ev['t'])
        birth_log.append({
            't': event_t, 'parent': int(ev['parent']), 'child': int(ev['child']), 'level': int(ev['level']),
            'ordinary_added': added, 'ordinary_encoded': encoded, 'pair_count_before': global_pair_count,
            'tet_count_after_ordinary': len(K.tets),
        })
        if not added or event_t < args.min_birth_time_before_pairing or global_pair_count >= args.max_auto_pairings:
            continue
        cascade_index = 1
        while cascade_index <= args.max_cascade_per_birth and global_pair_count < args.max_auto_pairings:
            if len(K.tets) < args.min_tets_before_pairing:
                break
            rows = nl.scan_rows(model, K, args)
            eval_rows = [r for r in rows if should_eval_candidate(r, args)]
            eval_rows.sort(key=lambda r: float(r.get('A_rank_score', 0.0) or 0.0), reverse=True)
            if args.max_eval_candidates > 0:
                eval_rows = eval_rows[:args.max_eval_candidates]
            chosen = nl.pick_nonlinear_pair(rows, args, cascade_index, used_faces)
            chosen_id = chosen.get('candidate_id') if chosen else None
            chosen_key = None
            if chosen is not None:
                fa_ch, fb_ch = parse_candidate_faces(chosen)
                chosen_key = (fa_ch, fb_ch, chosen.get('move_class'))
            local_eval = []
            for rank, r in enumerate(eval_rows, start=1):
                fa, fb = parse_candidate_faces(r)
                key = (fa, fb, r.get('move_class'))
                row = {
                    'variant': variant,
                    'phase_sign': args.phase_sign,
                    'event_t': event_t,
                    'scan_id': scan_id,
                    'cascade_index': cascade_index,
                    'eval_rank_A_score': rank,
                    'selected_by_current_rule': bool(r.get('candidate_id') == chosen_id or key == chosen_key),
                    'decision_used_delta_beta': False,
                    **candidate_base_row(r),
                    **candidate_alignment(model, K, r, args),
                }
                # Derived classification for the gate; not used by selection.
                d2 = int(float(row.get('delta_beta2', 0) or 0)) if str(row.get('delta_beta2','')).strip() != '' else 0
                row['beta2_opening_candidate'] = d2 > 0
                row['C_eigen_lock_lt_threshold'] = float(row.get('best_C_eigen_J_lock_mean_resid', 1.0) or 1.0) < args.lock_residual_threshold
                row['C_eigen_lock_max_lt_threshold'] = float(row.get('best_C_eigen_J_lock_max_resid', 1.0) or 1.0) < args.lock_max_threshold
                row['beta2_and_C_eigen_lock'] = bool(row['beta2_opening_candidate'] and row['C_eigen_lock_lt_threshold'] and row['C_eigen_lock_max_lt_threshold'])
                candidate_rows.append(row)
                local_eval.append(row)
            beta2_rows = [r for r in local_eval if r.get('beta2_opening_candidate')]
            good_rows = [r for r in local_eval if r.get('beta2_and_C_eigen_lock')]
            best_all = min([r for r in local_eval if r.get('candidate_eval_status') == 'ok'], key=lambda r: (float(r.get('best_C_eigen_J_lock_mean_resid', 1.0)), float(r.get('best_C_eigen_J_lock_max_resid', 1.0))), default=None)
            best_beta2 = min([r for r in beta2_rows if r.get('candidate_eval_status') == 'ok'], key=lambda r: (float(r.get('best_C_eigen_J_lock_mean_resid', 1.0)), float(r.get('best_C_eigen_J_lock_max_resid', 1.0))), default=None)
            scan_summaries.append({
                'variant': variant,
                'phase_sign': args.phase_sign,
                'event_t': event_t,
                'scan_id': scan_id,
                'cascade_index': cascade_index,
                'pair_count_before': global_pair_count,
                'candidate_count_evaluated': len(local_eval),
                'beta2_opening_candidate_count': len(beta2_rows),
                'beta2_and_C_eigen_lock_count': len(good_rows),
                'selected_candidate_id': chosen_id if chosen is not None else '',
                'selected_face_a': chosen.get('face_a','') if chosen else '',
                'selected_face_b': chosen.get('face_b','') if chosen else '',
                'best_all_C_eigen_lock_mean': best_all.get('best_C_eigen_J_lock_mean_resid','') if best_all else '',
                'best_all_candidate_id': best_all.get('candidate_id','') if best_all else '',
                'best_all_delta_beta2': best_all.get('delta_beta2','') if best_all else '',
                'best_beta2_C_eigen_lock_mean': best_beta2.get('best_C_eigen_J_lock_mean_resid','') if best_beta2 else '',
                'best_beta2_C_eigen_lock_max': best_beta2.get('best_C_eigen_J_lock_max_resid','') if best_beta2 else '',
                'best_beta2_candidate_id': best_beta2.get('candidate_id','') if best_beta2 else '',
                'best_beta2_transport_cosine': best_beta2.get('transport_cosine_ka_kb_reversed','') if best_beta2 else '',
                'best_beta2_comm_signed_birth_over_abs': best_beta2.get('comm_signed_birth_over_abs','') if best_beta2 else '',
            })
            scan_id += 1
            if chosen is None:
                break
            if float(chosen.get('nonlinear_cascade_score', 0.0) or 0.0) < args.min_nonlinear_score:
                break
            K2, log, applied = nl.apply_pair(model, K, chosen, args, event_t, cascade_index)
            pairing_log.append(log)
            if not applied:
                break
            fa, fb = parse_candidate_faces(chosen)
            if fa is not None and fb is not None:
                used_faces.add(fa); used_faces.add(fb)
            K = K2
            global_pair_count += 1
            cascade_index += 1
            if not args.cascade_rescan:
                break
    return K, birth_log, pairing_log, candidate_rows, scan_summaries


def summarize_candidates(candidate_rows: List[dict], args: argparse.Namespace) -> dict:
    ok = [r for r in candidate_rows if r.get('candidate_eval_status') == 'ok']
    beta2 = [r for r in ok if r.get('beta2_opening_candidate')]
    selected = [r for r in ok if r.get('selected_by_current_rule')]
    good = [r for r in ok if r.get('beta2_and_C_eigen_lock')]
    def best(rows: List[dict]) -> dict:
        if not rows:
            return {}
        return min(rows, key=lambda r: (float(r.get('best_C_eigen_J_lock_mean_resid', 1.0)), float(r.get('best_C_eigen_J_lock_max_resid', 1.0))))
    def mean_val(rows: List[dict], k: str) -> float:
        vals = [float(r.get(k, 0.0) or 0.0) for r in rows if str(r.get(k, '')).strip() != '']
        return float(np.mean(vals)) if vals else 0.0
    b_all = best(ok); b_beta2 = best(beta2); b_sel = best(selected)
    return {
        'candidate_count_ok': len(ok),
        'candidate_count_beta2_opening': len(beta2),
        'candidate_count_selected_by_rule': len(selected),
        'candidate_count_beta2_and_C_eigen_lock': len(good),
        'exists_beta2_and_C_eigen_lock': len(good) > 0,
        'best_all_C_eigen_lock_mean': b_all.get('best_C_eigen_J_lock_mean_resid', 0.0),
        'best_all_C_eigen_lock_max': b_all.get('best_C_eigen_J_lock_max_resid', 0.0),
        'best_all_delta_beta2': b_all.get('delta_beta2', ''),
        'best_all_candidate_id': b_all.get('candidate_id', ''),
        'best_beta2_C_eigen_lock_mean': b_beta2.get('best_C_eigen_J_lock_mean_resid', 0.0),
        'best_beta2_C_eigen_lock_max': b_beta2.get('best_C_eigen_J_lock_max_resid', 0.0),
        'best_beta2_candidate_id': b_beta2.get('candidate_id', ''),
        'best_beta2_selected_by_current_rule': b_beta2.get('selected_by_current_rule', False),
        'best_beta2_transport_cosine': b_beta2.get('transport_cosine_ka_kb_reversed', 0.0),
        'best_beta2_directed_imbalance': b_beta2.get('directed_imbalance', 0.0),
        'best_beta2_transverse_complementarity': b_beta2.get('transverse_complementarity', 0.0),
        'best_beta2_comm_signed_birth_over_abs': b_beta2.get('comm_signed_birth_over_abs', 0.0),
        'selected_C_eigen_lock_mean': b_sel.get('best_C_eigen_J_lock_mean_resid', 0.0),
        'selected_C_eigen_lock_max': b_sel.get('best_C_eigen_J_lock_max_resid', 0.0),
        'selected_transport_cosine_mean': mean_val(selected, 'transport_cosine_ka_kb_reversed'),
        'beta2_transport_cosine_mean': mean_val(beta2, 'transport_cosine_ka_kb_reversed'),
        'beta2_C_eigen_lock_mean_mean': mean_val(beta2, 'best_C_eigen_J_lock_mean_resid'),
        'beta2_comm_signed_birth_over_abs_mean': mean_val(beta2, 'comm_signed_birth_over_abs'),
    }


def run_variant(variant: str, phase_sign: int, args: argparse.Namespace, out: Path) -> dict:
    vname = f'{variant}_phase{phase_sign:+d}'.replace('+', 'plus').replace('-', 'minus')
    vout = out / vname
    vout.mkdir(parents=True, exist_ok=True)
    local_args = argparse.Namespace(**vars(args))
    local_args.phase_sign = phase_sign
    model = nl.build_model(variant, local_args)
    model.grow(local_args.max_level)
    baseline_K = core.build_dynamic_outward_ngf_complex(model)
    baseline_metrics = core.full_metrics(model, baseline_K, local_args.source)
    K, birth_log, pairing_log, candidate_rows, scan_summaries = scan_evaluate_and_apply(model, local_args, variant, out)
    auto_metrics = core.full_metrics(model, K, local_args.source)
    cand_summary = summarize_candidates(candidate_rows, local_args)
    write_csv(vout / 'birth_geometry_log.csv', birth_log)
    write_csv(vout / 'nonlinear_pairing_cascade_log.csv', pairing_log)
    write_csv(vout / 'all_pair_candidate_C_eigen_refinement_rows.csv', candidate_rows)
    write_csv(vout / 'scan_summaries.csv', scan_summaries)
    # Focus rows: beta2 opening candidates sorted by C-eigen residual.
    focus = [r for r in candidate_rows if r.get('candidate_eval_status') == 'ok' and r.get('beta2_opening_candidate')]
    focus.sort(key=lambda r: (float(r.get('best_C_eigen_J_lock_mean_resid', 1.0)), float(r.get('best_C_eigen_J_lock_max_resid', 1.0))))
    write_csv(vout / 'beta2_candidates_sorted_by_C_eigen_lock.csv', focus[:args.keep_top_candidates])
    summary = {
        'variant': variant,
        'variant_phase': vname,
        'phase_sign': phase_sign,
        'max_level': args.max_level,
        'source': args.source,
        'baseline_metrics': baseline_metrics,
        'auto_metrics': auto_metrics,
        'automatic_pairings_applied': sum(1 for x in pairing_log if x.get('applied')),
        'automatic_pairing_attempts_logged': len(pairing_log),
        'candidate_refinement_metrics': cand_summary,
        'interpretation_flags': {
            'beta2_opened_by_current_rule': auto_metrics['beta2'] > baseline_metrics['beta2'],
            'exists_beta2_and_C_eigen_lock': cand_summary['exists_beta2_and_C_eigen_lock'],
            'current_rule_missed_better_beta2_candidate': bool(cand_summary.get('best_beta2_selected_by_current_rule') is False and float(cand_summary.get('best_beta2_C_eigen_lock_mean', 1.0) or 1.0) < float(cand_summary.get('selected_C_eigen_lock_mean', 1.0) or 1.0)),
            'decision_used_delta_beta_any': any(str(r.get('decision_used_delta_beta','')).lower() == 'true' for r in pairing_log),
        },
    }
    (vout / 'variant_C_eigen_refinement_summary.json').write_text(json.dumps(summary, indent=2), encoding='utf-8')
    return summary


def write_comparative(out: Path, rows: List[dict]) -> None:
    flat = []
    for r in rows:
        m = r['candidate_refinement_metrics']; a = r['auto_metrics']
        flat.append({
            'variant_phase': r['variant_phase'],
            'variant': r['variant'],
            'phase_sign': r['phase_sign'],
            'beta0': a['beta0'], 'beta1': a['beta1'], 'beta2': a['beta2'], 'beta3': a['beta3'],
            'pairings': r['automatic_pairings_applied'],
            'candidate_count_ok': m['candidate_count_ok'],
            'candidate_count_beta2_opening': m['candidate_count_beta2_opening'],
            'exists_beta2_and_C_eigen_lock': m['exists_beta2_and_C_eigen_lock'],
            'candidate_count_beta2_and_C_eigen_lock': m['candidate_count_beta2_and_C_eigen_lock'],
            'best_beta2_C_eigen_lock_mean': m['best_beta2_C_eigen_lock_mean'],
            'best_beta2_C_eigen_lock_max': m['best_beta2_C_eigen_lock_max'],
            'best_beta2_candidate_id': m['best_beta2_candidate_id'],
            'best_beta2_selected_by_current_rule': m['best_beta2_selected_by_current_rule'],
            'selected_C_eigen_lock_mean': m['selected_C_eigen_lock_mean'],
            'selected_C_eigen_lock_max': m['selected_C_eigen_lock_max'],
            'best_beta2_transport_cosine': m['best_beta2_transport_cosine'],
            'best_beta2_directed_imbalance': m['best_beta2_directed_imbalance'],
            'best_beta2_transverse_complementarity': m['best_beta2_transverse_complementarity'],
            'best_beta2_comm_signed_birth_over_abs': m['best_beta2_comm_signed_birth_over_abs'],
            'beta2_transport_cosine_mean': m['beta2_transport_cosine_mean'],
        })
    write_csv(out / 'comparative_C_eigen_refinement_summary.csv', flat)


def make_docs(summary: dict) -> Tuple[str, str, str, str]:
    rows = summary['variant_rows']
    lines = [
        '| variant phase | beta auto | pairings | ok candidates | beta2 candidates | beta2 + C-lock | best beta2 lock | selected lock | best beta2 selected? | best beta2 cos |',
        '|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|',
    ]
    for r in rows:
        a = r['auto_metrics']; m = r['candidate_refinement_metrics']
        lines.append(f"| {r['variant_phase']} | ({a['beta0']},{a['beta1']},{a['beta2']},{a['beta3']}) | {r['automatic_pairings_applied']} | {m['candidate_count_ok']} | {m['candidate_count_beta2_opening']} | {m['candidate_count_beta2_and_C_eigen_lock']} | {float(m['best_beta2_C_eigen_lock_mean']):.6g} / {float(m['best_beta2_C_eigen_lock_max']):.6g} | {float(m['selected_C_eigen_lock_mean']):.6g} | {m['best_beta2_selected_by_current_rule']} | {float(m['best_beta2_transport_cosine']):.6g} |")
    table = '\n'.join(lines)
    summary_md = f"""# SUMMARY — C-eigen quadrature refinement gate

This package audits whether the current nonlinear asymmetry-gated pairing rule misses beta2-opening candidates that have a better native C-pair even/odd J-lock.

It evaluates legal pair candidates at each scan state before the current selection rule applies.  The gate does not fit a new rotation and does not introduce i, Hodge, positivity, a physical adjoint, or a final symmetrization.

{table}

Main question:

```text
Do there exist candidates with delta_beta2 > 0 and C-eigen J-lock residual below the threshold?
```
"""
    results_md = f"""# RESULTS — C-eigen quadrature refinement gate

## Comparative table

{table}

## Conservative interpretation

The test separates two possibilities:

```text
1. The current pairing rule selects the wrong pairs.
2. beta2-opening and C-eigen J-lock are structurally in tension in the available L2 candidate space.
```

A positive hit would require at least one beta2-opening candidate with both mean and max C-eigen J-lock residual below the configured thresholds.
"""
    audit_md = """# SOURCE AUDIT

Inherited structure:

- nonlinear asymmetry-gated complement pairing supplies the beta2-opening candidate space.
- the antisymmetric birth-transport vertex operator is used without final sym(M).
- C-pair/J-pair even/odd tests are inherited from the local pair algebra gate.

This package does not claim a complex structure.  It audits whether the existing candidate space already contains a beta2-opening pair with native C-eigen J alignment.
"""
    readme_md = """# C-eigen quadrature refinement gate

Run:

```bash
python3 test_C_eigen_quadrature_refinement_gate.py
```

Outputs include per-variant candidate scans, beta2 candidates sorted by C-eigen J-lock residual, comparative CSV/JSON, RESULTS.md, SUMMARY.md and a ZIP package.
"""
    return summary_md, results_md, audit_md, readme_md


def package(out: Path, zip_path: Path) -> None:
    files = [
        Path(__file__).name,
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
    ap.add_argument('--keep-top-candidates', type=int, default=160)
    ap.add_argument('--max-eval-candidates', type=int, default=0, help='0 means evaluate all legal pair candidates in each scan')
    ap.add_argument('--only-A-gate', action='store_true', default=True, help='default: evaluate only provenance A-gated pair candidates')
    ap.add_argument('--include-non-A-gate', action='store_false', dest='only_A_gate', help='audit all legal pair candidates, including candidates rejected by provenance A_gate')
    ap.add_argument('--antisym-eta', type=float, default=1.0)
    ap.add_argument('--erase-phase-for-strict-sym', action='store_true', default=True)
    ap.add_argument('--lock-residual-threshold', type=float, default=0.20)
    ap.add_argument('--lock-max-threshold', type=float, default=0.30)
    ap.add_argument('--variants', nargs='*', default=['real_growth', 'strict_symmetrized_control', 'no_backreaction'])
    ap.add_argument('--phase-signs', nargs='*', type=int, default=[1, -1])
    ap.add_argument('--out', default='C_eigen_quadrature_refinement_out_L2')
    ap.add_argument('--zip', default='cnna_C_eigen_quadrature_refinement_gate_pkg_L2.zip')
    args = ap.parse_args()

    out = Path(args.out)
    if out.exists():
        shutil.rmtree(out)
    out.mkdir(parents=True)
    rows = []
    for v in args.variants:
        for s in args.phase_signs:
            rows.append(run_variant(v, int(s), args, out))
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
                'variant_phase': r['variant_phase'],
                'beta': [r['auto_metrics'][f'beta{i}'] for i in range(4)],
                'pairings': r['automatic_pairings_applied'],
                'ok_candidates': r['candidate_refinement_metrics']['candidate_count_ok'],
                'beta2_candidates': r['candidate_refinement_metrics']['candidate_count_beta2_opening'],
                'beta2_and_C_eigen_lock': r['candidate_refinement_metrics']['candidate_count_beta2_and_C_eigen_lock'],
                'best_beta2_lock_mean': r['candidate_refinement_metrics']['best_beta2_C_eigen_lock_mean'],
                'best_beta2_lock_max': r['candidate_refinement_metrics']['best_beta2_C_eigen_lock_max'],
                'best_beta2_selected': r['candidate_refinement_metrics']['best_beta2_selected_by_current_rule'],
            } for r in rows
        ]
    }, indent=2))


if __name__ == '__main__':
    main()
