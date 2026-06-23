#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import math
import shutil
import zipfile
from pathlib import Path
from typing import Dict, List, Tuple, Optional

import numpy as np

import cnna_non_shelling_core as core
import test_nonlinear_asymmetry_cascade_growth as nl
from test_growth_with_asymmetry_gated_complement_pairing import ordinary_outward_step
import test_C_eigen_quadrature_refinement_gate as p62

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


def fval(row: dict, key: str, default: float = 0.0) -> float:
    try:
        v = row.get(key, default)
        if v is None or v == '':
            return default
        return float(v)
    except Exception:
        return default


def bval(row: dict, key: str) -> bool:
    return str(row.get(key, '')).lower() in {'true','1','yes'}


def parse_faces(row: dict) -> Tuple[Optional[Face], Optional[Face]]:
    try:
        fa = tuple(sorted(core.parse_face_string(str(row.get('face_a', '')))))
        fb_txt = str(row.get('face_b', '')).split('perm=')[0].strip()
        fb = tuple(sorted(core.parse_face_string(fb_txt)))
        if len(fa) == 3 and len(fb) == 3:
            return fa, fb
    except Exception:
        pass
    return None, None


def canonical_pair_key(row: dict) -> Tuple[str, str, str]:
    fa, fb = parse_faces(row)
    if fa is None or fb is None:
        return (str(row.get('move_class','')), str(row.get('face_a','')), str(row.get('face_b','')))
    a, b = tuple(fa), tuple(fb)
    # Pair candidates are unordered at the candidate-space audit level.
    if b < a:
        a, b = b, a
    return (str(row.get('move_class','')), ' '.join(map(str, a)), ' '.join(map(str, b)))


def apply_kappa_birth_order_reflection(model: core.DynamicProvenanceGrowth) -> None:
    """Apply the concrete sibling-birth-order reflection 1<->3, 2->2 in the model.

    This is stronger than a phase_sign flip: the Node.birth_order fields used by
    vertex_operator_directed(), birth_order normals, and child_ids_ordered() are
    actually changed.  Geometry, conductances, directed_edges and birth event ids
    are intentionally left fixed so exact face-pair matching is possible.
    """
    perm = {1: 3, 2: 2, 3: 1, 0: 0}
    for n in model.nodes.values():
        if n.parent is not None:
            n.birth_order = perm.get(int(n.birth_order), int(n.birth_order))
    for ev in model.birth_events:
        if 'order' in ev:
            ev['order'] = perm.get(int(ev['order']), int(ev['order']))


def candidate_base_row(row: dict) -> dict:
    keys = [
        'candidate_id','move_class','status','face_a','face_b','A_gate','A_invariant','A_rank_score',
        'directed_imbalance','transverse_complementarity','response_score','response_rank_legal',
        'A_rank_legal','A_rank_gated','delta_beta1','delta_beta2','delta_beta3','new_beta1','new_beta2',
        'A_nonreciprocal_norm','A_live_record_gap','A_shell_gap','A_fan_directed_mean',
        'address_similarity','centroid_distance','nonlinear_cascade_score'
    ]
    return {k: row.get(k, '') for k in keys}


def should_eval(row: dict, args: argparse.Namespace, used_faces: set[Face]) -> bool:
    if row.get('status') != 'ok':
        return False
    if row.get('move_class') not in {'handle_candidate','quotient_candidate'}:
        return False
    if args.only_A_gate and not bool(row.get('A_gate')):
        return False
    if not args.allow_reuse_faces:
        fa, fb = parse_faces(row)
        if fa in used_faces or fb in used_faces:
            return False
    return True


def build_scan_rows(variant: str, kappa_mode: str, args: argparse.Namespace, out: Path) -> Tuple[List[dict], List[dict]]:
    model = nl.build_model(variant, args)
    model.grow(args.max_level)
    if kappa_mode == 'kappa_reflect_birth_order':
        apply_kappa_birth_order_reflection(model)
    elif kappa_mode != 'identity':
        raise ValueError(kappa_mode)

    K = core.SimplicialComplex(f'{variant}_{kappa_mode}_candidate_space')
    root_seeded = False
    used_faces: set[Face] = set()
    scan_id = 0
    all_rows: List[dict] = []
    scan_summaries: List[dict] = []

    for ev in sorted(model.birth_events, key=lambda x: int(x['t'])):
        root_seeded, added, encoded = ordinary_outward_step(model, K, ev, root_seeded)
        event_t = int(ev['t'])
        if not added or event_t < args.min_birth_time_before_pairing:
            continue
        if len(K.tets) < args.min_tets_before_pairing:
            continue
        # Candidate-space audit only.  No candidate is applied here.
        rows = nl.scan_rows(model, K, args)
        legal = [r for r in rows if should_eval(r, args, used_faces)]
        legal.sort(key=lambda r: fval(r, 'A_rank_score', 0.0), reverse=True)
        if args.max_eval_candidates > 0:
            legal = legal[:args.max_eval_candidates]
        chosen_A = nl.pick_nonlinear_pair(rows, args, 1, used_faces)
        chosen_A_key = canonical_pair_key(chosen_A) if chosen_A is not None else None
        local: List[dict] = []
        for rank, r in enumerate(legal, start=1):
            rr = dict(r)
            rr['nonlinear_cascade_score'] = nl.nonlinear_score(rr, args, 1)
            base = candidate_base_row(rr)
            align = p62.candidate_alignment(model, K, rr, args)
            row = {
                'variant': variant,
                'kappa_mode': kappa_mode,
                'event_t': event_t,
                'scan_id': scan_id,
                'cascade_index': 1,
                'eval_rank_A_score': rank,
                'candidate_pair_key': '|'.join(canonical_pair_key(rr)),
                'selected_by_A_rank_rule': canonical_pair_key(rr) == chosen_A_key,
                'decision_used_delta_beta': False,
                **base,
                **align,
            }
            d2 = int(fval(row, 'delta_beta2', 0.0))
            row['beta2_opening_audit_only'] = d2 > 0
            row['C_eigen_lock_lt_threshold'] = fval(row, 'best_C_eigen_J_lock_mean_resid', 1.0) < args.lock_residual_threshold
            row['C_eigen_lock_max_lt_threshold'] = fval(row, 'best_C_eigen_J_lock_max_resid', 1.0) < args.lock_max_threshold
            row['beta2_and_C_eigen_lock'] = bool(row['beta2_opening_audit_only'] and row['C_eigen_lock_lt_threshold'] and row['C_eigen_lock_max_lt_threshold'])
            all_rows.append(row)
            local.append(row)
        ok = [r for r in local if r.get('candidate_eval_status') == 'ok']
        beta2 = [r for r in ok if r.get('beta2_opening_audit_only')]
        good = [r for r in ok if r.get('beta2_and_C_eigen_lock')]
        best = min(ok, key=lambda r: (fval(r, 'best_C_eigen_J_lock_mean_resid', 1.0), fval(r, 'best_C_eigen_J_lock_max_resid', 1.0)), default=None)
        scan_summaries.append({
            'variant': variant,
            'kappa_mode': kappa_mode,
            'event_t': event_t,
            'scan_id': scan_id,
            'tet_count': len(K.tets),
            'candidate_count_ok': len(ok),
            'candidate_count_beta2_opening_audit_only': len(beta2),
            'candidate_count_beta2_and_C_lock': len(good),
            'best_candidate_id': best.get('candidate_id','') if best else '',
            'best_pair_key': best.get('candidate_pair_key','') if best else '',
            'best_C_lock_mean': best.get('best_C_eigen_J_lock_mean_resid','') if best else '',
            'best_delta_beta2_audit_only': best.get('delta_beta2','') if best else '',
            'selected_A_candidate_id': chosen_A.get('candidate_id','') if chosen_A else '',
            'selected_A_pair_key': '|'.join(chosen_A_key) if chosen_A_key else '',
        })
        scan_id += 1
    vout = out / f'{variant}_{kappa_mode}'
    write_csv(vout / 'candidate_eval_rows.csv', all_rows)
    write_csv(vout / 'scan_summaries.csv', scan_summaries)
    return all_rows, scan_summaries


def match_identity_kappa(identity_rows: List[dict], kappa_rows: List[dict], args: argparse.Namespace) -> Tuple[List[dict], List[dict]]:
    groups: Dict[Tuple[str,str,str,str,str], Dict[str,dict]] = {}
    for r in identity_rows:
        key = (str(r.get('variant','')), str(r.get('event_t','')), str(r.get('scan_id','')), str(r.get('cascade_index','')), str(r.get('candidate_pair_key','')))
        groups.setdefault(key, {})['identity'] = r
    for r in kappa_rows:
        key = (str(r.get('variant','')), str(r.get('event_t','')), str(r.get('scan_id','')), str(r.get('cascade_index','')), str(r.get('candidate_pair_key','')))
        groups.setdefault(key, {})['kappa'] = r
    matched: List[dict] = []
    unmatched: List[dict] = []
    for key, d in groups.items():
        if 'identity' not in d or 'kappa' not in d:
            row = {'variant': key[0], 'event_t': key[1], 'scan_id': key[2], 'cascade_index': key[3], 'candidate_pair_key': key[4], 'has_identity': 'identity' in d, 'has_kappa': 'kappa' in d}
            unmatched.append(row)
            continue
        a, b = d['identity'], d['kappa']
        s_a = fval(a, 'comm_signed_birth_over_abs', 0.0)
        s_b = fval(b, 'comm_signed_birth_over_abs', 0.0)
        abs_a, abs_b = abs(s_a), abs(s_b)
        flip_score = (s_a + s_b) / (abs_a + abs_b + EPS)
        amp_balance = abs(abs_a - abs_b) / (abs_a + abs_b + EPS)
        c_a = fval(a, 'best_C_eigen_J_lock_mean_resid', 1.0)
        c_b = fval(b, 'best_C_eigen_J_lock_mean_resid', 1.0)
        cmax_a = fval(a, 'best_C_eigen_J_lock_max_resid', 1.0)
        cmax_b = fval(b, 'best_C_eigen_J_lock_max_resid', 1.0)
        q_a, p_a = fval(a, 'Q_norm', 0.0), fval(a, 'P_norm', 0.0)
        q_b, p_b = fval(b, 'Q_norm', 0.0), fval(b, 'P_norm', 0.0)
        qp_a = min(q_a, p_a) / (max(q_a, p_a) + EPS) if max(q_a,p_a) > EPS else 0.0
        qp_b = min(q_b, p_b) / (max(q_b, p_b) + EPS) if max(q_b,p_b) > EPS else 0.0
        row = {
            'variant': key[0],
            'event_t': key[1],
            'scan_id': key[2],
            'cascade_index': key[3],
            'candidate_pair_key': key[4],
            'candidate_id_identity': a.get('candidate_id',''),
            'candidate_id_kappa': b.get('candidate_id',''),
            'move_class_identity': a.get('move_class',''),
            'move_class_kappa': b.get('move_class',''),
            'A_gate_identity': a.get('A_gate',''),
            'A_gate_kappa': b.get('A_gate',''),
            'beta2_opening_identity': a.get('beta2_opening_audit_only',''),
            'beta2_opening_kappa': b.get('beta2_opening_audit_only',''),
            'delta_beta2_identity': a.get('delta_beta2',''),
            'delta_beta2_kappa': b.get('delta_beta2',''),
            'C_lock_identity': c_a,
            'C_lock_kappa': c_b,
            'C_lock_worst': max(c_a, c_b),
            'C_lock_max_worst': max(cmax_a, cmax_b),
            'signed_birth_identity': s_a,
            'signed_birth_kappa': s_b,
            'signed_flip_score_zero_if_perfect': flip_score,
            'signed_flip_abs': abs(flip_score),
            'signed_amplitude_min': min(abs_a, abs_b),
            'signed_amplitude_avg': 0.5 * (abs_a + abs_b),
            'signed_amplitude_balance_zero_if_equal': amp_balance,
            'Q_norm_identity': q_a,
            'P_norm_identity': p_a,
            'Q_norm_kappa': q_b,
            'P_norm_kappa': p_b,
            'QP_balance_identity': qp_a,
            'QP_balance_kappa': qp_b,
            'QP_balance_min': min(qp_a, qp_b),
            'comm_abs_area_identity': fval(a, 'comm_abs_area', 0.0),
            'comm_abs_area_kappa': fval(b, 'comm_abs_area', 0.0),
            'transport_cosine_identity': fval(a, 'transport_cosine_ka_kb_reversed', 0.0),
            'transport_cosine_kappa': fval(b, 'transport_cosine_ka_kb_reversed', 0.0),
            'directed_imbalance_identity': fval(a, 'directed_imbalance', 0.0),
            'directed_imbalance_kappa': fval(b, 'directed_imbalance', 0.0),
            'directed_imbalance_avg': 0.5*(fval(a,'directed_imbalance',0.0)+fval(b,'directed_imbalance',0.0)),
            'transverse_complementarity_identity': fval(a, 'transverse_complementarity', 0.0),
            'transverse_complementarity_kappa': fval(b, 'transverse_complementarity', 0.0),
            'transverse_complementarity_avg': 0.5*(fval(a,'transverse_complementarity',0.0)+fval(b,'transverse_complementarity',0.0)),
            'selected_A_identity': a.get('selected_by_A_rank_rule',''),
            'selected_A_kappa': b.get('selected_by_A_rank_rule',''),
            'decision_used_delta_beta': False,
        }
        row['passes_A_gate_both'] = bval(row, 'A_gate_identity') and bval(row, 'A_gate_kappa')
        row['passes_beta2_audit_both'] = bval(row, 'beta2_opening_identity') and bval(row, 'beta2_opening_kappa')
        row['passes_C_lock_worst'] = row['C_lock_worst'] < args.lock_residual_threshold and row['C_lock_max_worst'] < args.lock_max_threshold
        row['passes_QP_proxy'] = row['QP_balance_min'] > args.qp_balance_min and min(row['comm_abs_area_identity'], row['comm_abs_area_kappa']) > args.comm_abs_min
        row['passes_kappa_signed_flip'] = row['signed_flip_abs'] < args.flip_abs_threshold and row['signed_amplitude_min'] > args.signed_amp_min
        row['passes_all_kappa_pareto_gate'] = bool(row['passes_A_gate_both'] and row['passes_C_lock_worst'] and row['passes_QP_proxy'] and row['passes_kappa_signed_flip'])
        row['passes_all_including_beta2_audit'] = bool(row['passes_all_kappa_pareto_gate'] and row['passes_beta2_audit_both'])
        row['kappa_pareto_score'] = (
            (1.0 / (row['C_lock_worst'] + 0.05))
            * (1.0 - min(1.0, row['signed_flip_abs']))
            * row['signed_amplitude_min']
            * row['QP_balance_min']
            * (0.25 + row['directed_imbalance_avg'])
            * (0.25 + row['transverse_complementarity_avg'])
        )
        matched.append(row)
    return matched, unmatched


def pareto_front(rows: List[dict]) -> List[dict]:
    def dominates(a: dict, b: dict) -> bool:
        dims = [
            ('C_lock_worst', -1),
            ('signed_flip_abs', -1),
            ('signed_amplitude_min', 1),
            ('QP_balance_min', 1),
            ('directed_imbalance_avg', 1),
            ('transverse_complementarity_avg', 1),
        ]
        ge = True; gt = False
        for k, sign in dims:
            av, bv = fval(a,k), fval(b,k)
            if sign == 1:
                if av + 1e-12 < bv: return False
                if av > bv + 1e-12: gt = True
            else:
                if av > bv + 1e-12: return False
                if av + 1e-12 < bv: gt = True
        return gt
    out=[]
    for r in rows:
        if not any(dominates(o, r) for o in rows if o is not r):
            out.append(r)
    out.sort(key=lambda r: (not r.get('passes_all_kappa_pareto_gate'), fval(r,'signed_flip_abs',1), fval(r,'C_lock_worst',1), -fval(r,'kappa_pareto_score',0)))
    return out


def summarize_variant(variant: str, matched: List[dict], unmatched: List[dict], args: argparse.Namespace) -> dict:
    both_A = [r for r in matched if r.get('passes_A_gate_both')]
    both_beta = [r for r in both_A if r.get('passes_beta2_audit_both')]
    Cpass = [r for r in both_A if r.get('passes_C_lock_worst')]
    flip = [r for r in both_A if r.get('passes_kappa_signed_flip')]
    allp = [r for r in both_A if r.get('passes_all_kappa_pareto_gate')]
    allbeta = [r for r in both_A if r.get('passes_all_including_beta2_audit')]
    best = min(both_A, key=lambda r: (fval(r,'signed_flip_abs',1), fval(r,'C_lock_worst',1), -fval(r,'signed_amplitude_min',0)), default=None)
    bestC = min(both_A, key=lambda r: (fval(r,'C_lock_worst',1), fval(r,'signed_flip_abs',1)), default=None)
    bestAllScore = max(both_A, key=lambda r: fval(r,'kappa_pareto_score',0), default=None)
    return {
        'variant': variant,
        'matched_candidate_count': len(matched),
        'unmatched_candidate_count': len(unmatched),
        'A_gated_both_count': len(both_A),
        'beta2_opening_audit_both_count': len(both_beta),
        'C_lock_pass_count': len(Cpass),
        'kappa_signed_flip_pass_count': len(flip),
        'all_kappa_pareto_gate_pass_count': len(allp),
        'all_kappa_pareto_plus_beta2_audit_pass_count': len(allbeta),
        'decision_used_delta_beta_any': False,
        'best_flip_candidate_key': best.get('candidate_pair_key','') if best else '',
        'best_flip_signed_flip_abs': best.get('signed_flip_abs','') if best else '',
        'best_flip_signed_identity': best.get('signed_birth_identity','') if best else '',
        'best_flip_signed_kappa': best.get('signed_birth_kappa','') if best else '',
        'best_flip_C_lock_worst': best.get('C_lock_worst','') if best else '',
        'best_C_candidate_key': bestC.get('candidate_pair_key','') if bestC else '',
        'best_C_lock_worst': bestC.get('C_lock_worst','') if bestC else '',
        'best_C_signed_flip_abs': bestC.get('signed_flip_abs','') if bestC else '',
        'best_score_candidate_key': bestAllScore.get('candidate_pair_key','') if bestAllScore else '',
        'best_score': bestAllScore.get('kappa_pareto_score','') if bestAllScore else '',
        'best_score_C_lock_worst': bestAllScore.get('C_lock_worst','') if bestAllScore else '',
        'best_score_signed_flip_abs': bestAllScore.get('signed_flip_abs','') if bestAllScore else '',
        'best_score_beta2_both': bestAllScore.get('passes_beta2_audit_both','') if bestAllScore else '',
    }


def make_docs(summary: dict, rows: List[dict]) -> Tuple[str,str,str,str]:
    table_lines = ['| variant | matched | A both | beta2 audit both | C-lock pass | kappa flip pass | all pass | all+beta2 | best C lock | best C flip | used dBeta? |', '|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|']
    for r in rows:
        table_lines.append(f"| {r['variant']} | {r['matched_candidate_count']} | {r['A_gated_both_count']} | {r['beta2_opening_audit_both_count']} | {r['C_lock_pass_count']} | {r['kappa_signed_flip_pass_count']} | {r['all_kappa_pareto_gate_pass_count']} | {r['all_kappa_pareto_plus_beta2_audit_pass_count']} | {float(r['best_C_lock_worst'] or 0):.6g} | {float(r['best_C_signed_flip_abs'] or 0):.6g} | {r['decision_used_delta_beta_any']} |")
    table = '\n'.join(table_lines)
    smd = f"""# SUMMARY — kappa-permuted candidate Pareto gate

Model label:
CNNA growing primal simplicial complex with deterministic sequential provenance growth, A-gated complement-pair candidate space, directed antisymmetric birth-transport operators, local C/J pair algebra, and explicit sibling birth-order reflection `1 <-> 3` as a model-level kappa audit.

This package replaces the earlier internal `phase_sign +/-1` audit with a concrete model-label transformation:

```text
identity model:       original birth_order labels
kappa-reflected model: Node.birth_order is reflected inside every sibling fan, 1<->3, 2->2
```

The reflection is applied to the model's birth_order fields, not to a standalone phase-sign parameter.  Geometry, conductances, directed_edges and vertex IDs are kept fixed so the same face-pair candidates can be matched exactly.

{table}

Read conservatively: this is a label-kappa audit on the same grown model, not a full re-growth with reversed birth sequence.
"""
    rmd = f"""# RESULTS — kappa-permuted candidate Pareto gate

## Comparative table

{table}

## Gate logic

A successful candidate would satisfy:

```text
A_gate in identity and kappa-reflected model,
C-eigen J-lock residual below threshold in both,
Q/P support in both,
signed_birth flips under kappa reflection,
optional beta2-opening audit in both.
```

`delta_beta2` is audit-only and is never used to form the candidate space or decide a move.
"""
    audit = """# SOURCE AUDIT

This test is a derived-only candidate audit.

Anti-smuggling constraints:

- no complex scalar and no `i`;
- no Hodge star, no positivity, no physical adjoint;
- no final `sym(M)` in the directed birth-transport vertex operator;
- no arbitrary fitted rotation;
- no delta_beta/H2/kappa used as a move decision input;
- kappa is implemented as the concrete sibling birth-order reflection 1<->3 on `Node.birth_order`, not as a mere phase-sign flag.

Limitation:
The kappa audit keeps geometry, conductances, directed_edges and node IDs fixed to allow exact candidate matching.  It is stronger than a phase-sign flip but weaker than a full re-growth under reversed birth sequence.  A full address-permuted re-growth would require a separate address-level candidate matcher.
"""
    readme = """# Kappa-permuted candidate Pareto gate

Run:

```bash
python3 test_kappa_permuted_candidate_pareto_gate.py
```

Outputs include per-variant identity/kappa candidate rows, matched candidate rows, Pareto fronts, RESULTS.md, SUMMARY.md, SOURCE_AUDIT.md, and a ZIP package.
"""
    return smd, rmd, audit, readme


def package(out: Path, zip_path: Path) -> None:
    files = [
        Path(__file__).name,
        'test_C_eigen_guided_phase_robustness_pareto_gate.py',
        'test_C_eigen_guided_pairing_rule_gate.py',
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
    ap.add_argument('--keep-top-candidates', type=int, default=120)
    ap.add_argument('--keep-top-faces', type=int, default=80)
    ap.add_argument('--max-eval-candidates', type=int, default=0, help='0 means evaluate all legal candidates after filters')
    ap.add_argument('--only-A-gate', action='store_true', default=True)
    ap.add_argument('--harmonic-positive-threshold', type=float, default=1e-4)
    ap.add_argument('--phase-sign', type=int, default=1, help='Kept fixed; kappa is tested by birth_order reflection, not phase sign.')
    ap.add_argument('--antisym-eta', type=float, default=1.0)
    ap.add_argument('--erase-phase-for-strict-sym', action='store_true', default=True)
    ap.add_argument('--lock-residual-threshold', type=float, default=0.20)
    ap.add_argument('--lock-max-threshold', type=float, default=0.30)
    ap.add_argument('--flip-abs-threshold', type=float, default=0.25)
    ap.add_argument('--signed-amp-min', type=float, default=0.08)
    ap.add_argument('--qp-balance-min', type=float, default=0.25)
    ap.add_argument('--comm-abs-min', type=float, default=1e-8)
    ap.add_argument('--variants', nargs='*', default=['real_growth','strict_symmetrized_control','no_backreaction'])
    ap.add_argument('--out', default='kappa_permuted_candidate_pareto_out_L2')
    ap.add_argument('--zip', default='cnna_kappa_permuted_candidate_pareto_gate_pkg_L2.zip')
    args = ap.parse_args()

    out = Path(args.out)
    if out.exists():
        shutil.rmtree(out)
    out.mkdir(parents=True)

    all_summary_rows: List[dict] = []
    all_matched: List[dict] = []
    all_unmatched: List[dict] = []
    all_pareto: List[dict] = []
    for variant in args.variants:
        ident_rows, ident_scans = build_scan_rows(variant, 'identity', args, out)
        kappa_rows, kappa_scans = build_scan_rows(variant, 'kappa_reflect_birth_order', args, out)
        matched, unmatched = match_identity_kappa(ident_rows, kappa_rows, args)
        for r in unmatched:
            r['variant'] = variant
        front = pareto_front([r for r in matched if r.get('passes_A_gate_both')])
        vout = out / variant
        write_csv(vout / 'matched_identity_kappa_candidate_rows.csv', matched)
        write_csv(vout / 'unmatched_identity_kappa_candidate_rows.csv', unmatched)
        write_csv(vout / 'pareto_front_identity_kappa_rows.csv', front)
        write_csv(vout / 'top_kappa_candidate_rows.csv', sorted(matched, key=lambda r: (not r.get('passes_all_kappa_pareto_gate'), fval(r,'signed_flip_abs',1), fval(r,'C_lock_worst',1), -fval(r,'kappa_pareto_score',0)))[:80])
        s = summarize_variant(variant, matched, unmatched, args)
        all_summary_rows.append(s)
        all_matched.extend(matched)
        all_unmatched.extend(unmatched)
        all_pareto.extend(front)
    write_csv(out / 'comparative_kappa_pareto_summary.csv', all_summary_rows)
    write_csv(out / 'matched_identity_kappa_candidate_rows.csv', all_matched)
    write_csv(out / 'unmatched_identity_kappa_candidate_rows.csv', all_unmatched)
    write_csv(out / 'pareto_front_identity_kappa_rows.csv', all_pareto)
    summary = {'args': vars(args), 'summary_rows': all_summary_rows}
    (out / 'comparative_summary.json').write_text(json.dumps(summary, indent=2), encoding='utf-8')
    smd, rmd, audit, readme = make_docs(summary, all_summary_rows)
    (out / 'SUMMARY.md').write_text(smd, encoding='utf-8')
    (out / 'RESULTS.md').write_text(rmd, encoding='utf-8')
    (out / 'SOURCE_AUDIT.md').write_text(audit, encoding='utf-8')
    (out / 'README.md').write_text(readme, encoding='utf-8')
    package(out, Path(args.zip))
    print(json.dumps({'zip': args.zip, 'out': args.out, 'summary': all_summary_rows}, indent=2))


if __name__ == '__main__':
    main()
