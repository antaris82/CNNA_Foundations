#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import math
import shutil
import zipfile
from pathlib import Path
from typing import List, Optional, Dict, Tuple

import numpy as np

import cnna_non_shelling_core as core
import test_nonlinear_asymmetry_cascade_growth as nl
import test_dual_assembly_order_context_ablation_gate as p69
import test_signed_Jlock_role_coupling_gate as p70
import test_assembly_motif_basis_diagonalization_gate as p71
import test_pairing_transport_antisym_birth_coherence_gate as p56
import test_signed_quadrature_area_kappa_gate as p59
import test_pair_J_alignment_search_gate as p61
import test_response_monodromy_to_QP_transfer_gate as p78
import test_bounded_response_bridge_sanity_gate as p79

EPS = 1e-12

RETARDED_OPS = {'markov','forward_cycle','cycle_difference','signed_C3_skew','skew_markov'}
ADVANCED_OPS = {'reverse_cycle'}
SYMMETRIC_OPS = {'cycle_sum'}
LONGITUDINAL_BRIDGES = {'radial','birth_h'}
TRANSVERSE_BRIDGES = {'birth_q','conductance_birth_q','record_birth_q'}


def fbool(x) -> bool:
    if isinstance(x, str):
        return x.lower() in {'true','1','yes'}
    return bool(x)


def finite_float(x, default=math.nan) -> float:
    try:
        v = float(x)
        return v if np.isfinite(v) else default
    except Exception:
        return default


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


def parse_fan_circs(row: dict) -> List[float]:
    try:
        fans = json.loads(row.get('fans_json','[]'))
        return [float(f.get('circulation_log_ratio', 0.0)) for f in fans]
    except Exception:
        return []


def causalize_row(row: dict, args: argparse.Namespace) -> dict:
    op = str(row.get('response_op',''))
    bridge = str(row.get('bridge_mode',''))
    good = bool(
        row.get('transfer_status') == 'ok' and
        fbool(row.get('bridge_well_conditioned_pass')) and
        fbool(row.get('response_bounded_mode')) and
        not fbool(row.get('artifact_suspect')) and
        fbool(row.get('transfer_active'))
    )
    retarded = op in RETARDED_OPS
    advanced = op in ADVANCED_OPS
    symmetric = op in SYMMETRIC_OPS
    longitudinal = bridge in LONGITUDINAL_BRIDGES
    transverse = bridge in TRANSVERSE_BRIDGES
    circs = parse_fan_circs(row)
    signed_mean_circ = float(np.mean(circs)) if circs else 0.0
    signed_sum_circ = float(np.sum(circs)) if circs else 0.0
    same_sign = bool(circs and (all(c >= -args.circ_sign_tol for c in circs) or all(c <= args.circ_sign_tol for c in circs)))
    circ_bias = abs(signed_sum_circ) / (sum(abs(c) for c in circs) + EPS) if circs else 0.0
    halfplane_bias = finite_float(row.get('halfplane_bias_pos_minus_neg'), 0.0)
    opnorm = finite_float(row.get('projected_QP_opnorm'), math.inf)
    normality = finite_float(row.get('projected_QP_normality_resid'), math.inf)
    qp_resid = finite_float(row.get('QP_transfer_mean_resid'), math.inf)
    leak_like = finite_float(row.get('projected_J2_plus_I_resid'), math.inf)  # not used as J gate; retained diagnostic.
    bounded_contraction = bool(opnorm <= args.retarded_opnorm_threshold)
    nonzero_circ = bool(finite_float(row.get('max_abs_response_circulation'), 0.0) >= args.min_response_circulation)
    stable_halfplane = bool(abs(halfplane_bias) >= args.halfplane_bias_threshold)
    causal_order_signal = bool(good and retarded and nonzero_circ and same_sign and circ_bias >= args.circ_bias_threshold)
    longitudinal_retarded_signal = bool(causal_order_signal and longitudinal and bounded_contraction)
    boundary_value_signal = bool(longitudinal_retarded_signal and stable_halfplane)
    # Retarded selection is deliberately weaker than Q/P-lock: it asks for bounded forward pre-causal half-plane trace.
    retarded_boundary_pass = bool(
        boundary_value_signal and
        normality <= args.retarded_normality_threshold and
        qp_resid <= args.loose_qp_leakage_threshold
    )
    # Mark advanced/symmetric leakage separately; not failure by itself, only evidence against clean retarded selection.
    advanced_boundary_signal = bool(good and advanced and nonzero_circ and bounded_contraction and stable_halfplane)
    symmetric_boundary_signal = bool(good and symmetric and bounded_contraction and stable_halfplane)
    return {
        **row,
        'causal_class': 'retarded' if retarded else ('advanced' if advanced else ('symmetric' if symmetric else 'other')),
        'bridge_axis_class': 'longitudinal' if longitudinal else ('transverse' if transverse else 'other'),
        'good_bounded_nonartifact_active': good,
        'retarded_op_audit': retarded,
        'advanced_op_audit': advanced,
        'symmetric_op_audit': symmetric,
        'longitudinal_bridge_audit': longitudinal,
        'transverse_bridge_audit': transverse,
        'signed_mean_response_circulation': signed_mean_circ,
        'signed_sum_response_circulation': signed_sum_circ,
        'response_circulation_same_sign_across_fans': same_sign,
        'response_circulation_directional_bias': circ_bias,
        'bounded_contraction_precausal_pass': bounded_contraction,
        'stable_halfplane_bias_pass': stable_halfplane,
        'causal_order_signal_pass': causal_order_signal,
        'longitudinal_retarded_signal_pass': longitudinal_retarded_signal,
        'boundary_value_signal_pass': boundary_value_signal,
        'advanced_boundary_signal_pass': advanced_boundary_signal,
        'symmetric_boundary_signal_pass': symmetric_boundary_signal,
        'retarded_boundary_value_gate_pass': retarded_boundary_pass,
        'loose_QP_leakage_resid_diagnostic': qp_resid,
        'J2_resid_not_primary_diagnostic': leak_like,
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

    raw_rows: List[dict] = []
    for a in assembly_log:
        if fbool(a.get('assembly_applied')):
            raw_rows.extend(p79.bounded_lift_rows(model, K, a, args))
    rows = [causalize_row(r, args) for r in raw_rows]
    summary = summarize(rows)

    write_csv(vout / 'birth_geometry_log.csv', birth_log)
    write_csv(vout / 'assembly_pairing_log.csv', pairing_log)
    write_csv(vout / 'assembly_log.csv', assembly_log)
    write_csv(vout / 'candidate_eval_rows.csv', candidate_rows)
    write_csv(vout / 'retarded_boundary_bridge_rows.csv', rows)
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
        'retarded_boundary_summary': summary,
        'decision_used_delta_beta_any': decision_used_delta,
    }
    (vout / 'summary.json').write_text(json.dumps(res, indent=2), encoding='utf-8')
    return res


def vals(rows: List[dict], key: str) -> List[float]:
    out=[]
    for r in rows:
        v=finite_float(r.get(key), math.nan)
        if np.isfinite(v): out.append(v)
    return out


def summarize(rows: List[dict]) -> dict:
    ok = [r for r in rows if r.get('transfer_status') == 'ok']
    good = [r for r in ok if fbool(r.get('good_bounded_nonartifact_active'))]
    ret_good = [r for r in good if fbool(r.get('retarded_op_audit'))]
    adv_good = [r for r in good if fbool(r.get('advanced_op_audit'))]
    sym_good = [r for r in good if fbool(r.get('symmetric_op_audit'))]
    long_ret = [r for r in ret_good if fbool(r.get('longitudinal_bridge_audit'))]
    trans_ret = [r for r in ret_good if fbool(r.get('transverse_bridge_audit'))]
    def count(rs,k): return sum(1 for r in rs if fbool(r.get(k)))
    def mn(rs,k,default=0.0):
        v=vals(rs,k); return float(min(v)) if v else default
    def mx(rs,k,default=0.0):
        v=vals(rs,k); return float(max(v)) if v else default
    def mean(rs,k,default=0.0):
        v=vals(rs,k); return float(np.mean(v)) if v else default
    best_ret = min(long_ret, key=lambda r: finite_float(r.get('loose_QP_leakage_resid_diagnostic'), 99.0) - 0.1*abs(finite_float(r.get('halfplane_bias_pos_minus_neg'), 0.0)), default=None)
    best_boundary = next((r for r in long_ret if fbool(r.get('retarded_boundary_value_gate_pass'))), None)
    # Paired comparison: same bridge/lift/motif but forward vs reverse_cycle.
    groups: Dict[Tuple, Dict[str, List[dict]]] = {}
    for r in good:
        key = (r.get('event_t'), r.get('scan_id'), r.get('bridge_mode'), r.get('lift_mode'), r.get('union_faces'))
        groups.setdefault(key, {'forward': [], 'reverse': [], 'sym': []})
        if r.get('response_op') == 'forward_cycle': groups[key]['forward'].append(r)
        if r.get('response_op') == 'reverse_cycle': groups[key]['reverse'].append(r)
        if r.get('response_op') == 'cycle_sum': groups[key]['sym'].append(r)
    pair_rows=[]
    for key,g in groups.items():
        if g['forward'] and g['reverse']:
            f = min(g['forward'], key=lambda x: finite_float(x.get('QP_transfer_mean_resid'),99.0))
            rv = min(g['reverse'], key=lambda x: finite_float(x.get('QP_transfer_mean_resid'),99.0))
            pair_rows.append({
                'key': str(key),
                'forward_qp': finite_float(f.get('QP_transfer_mean_resid'), math.nan),
                'reverse_qp': finite_float(rv.get('QP_transfer_mean_resid'), math.nan),
                'forward_halfplane': finite_float(f.get('halfplane_bias_pos_minus_neg'), 0.0),
                'reverse_halfplane': finite_float(rv.get('halfplane_bias_pos_minus_neg'), 0.0),
                'forward_opnorm': finite_float(f.get('projected_QP_opnorm'), math.inf),
                'reverse_opnorm': finite_float(rv.get('projected_QP_opnorm'), math.inf),
            })
    forward_better_count=sum(1 for p in pair_rows if p['forward_qp'] < p['reverse_qp'])
    halfplane_opposed_count=sum(1 for p in pair_rows if abs(p['forward_halfplane'] + p['reverse_halfplane']) <= 0.25 and abs(p['forward_halfplane']) >= 0.5 and abs(p['reverse_halfplane']) >= 0.5)
    return {
        'transfer_row_count': len(rows),
        'ok_transfer_row_count': len(ok),
        'good_bounded_nonartifact_active_count': len(good),
        'good_retarded_count': len(ret_good),
        'good_advanced_count': len(adv_good),
        'good_symmetric_count': len(sym_good),
        'good_longitudinal_retarded_count': len(long_ret),
        'good_transverse_retarded_count': len(trans_ret),
        'causal_order_signal_count': count(good, 'causal_order_signal_pass'),
        'longitudinal_retarded_signal_count': count(good, 'longitudinal_retarded_signal_pass'),
        'boundary_value_signal_count': count(good, 'boundary_value_signal_pass'),
        'retarded_boundary_value_gate_pass_count': count(good, 'retarded_boundary_value_gate_pass'),
        'advanced_boundary_signal_count': count(good, 'advanced_boundary_signal_pass'),
        'symmetric_boundary_signal_count': count(good, 'symmetric_boundary_signal_pass'),
        'bounded_contraction_good_count': count(good, 'bounded_contraction_precausal_pass'),
        'stable_halfplane_good_count': count(good, 'stable_halfplane_bias_pass'),
        'best_longitudinal_retarded_QP_resid': mn(long_ret, 'loose_QP_leakage_resid_diagnostic'),
        'best_longitudinal_retarded_J2_diag_resid': mn(long_ret, 'J2_resid_not_primary_diagnostic'),
        'max_longitudinal_retarded_halfplane_abs_bias': max([abs(x) for x in vals(long_ret,'halfplane_bias_pos_minus_neg')] or [0.0]),
        'mean_longitudinal_retarded_circ_bias': mean(long_ret, 'response_circulation_directional_bias'),
        'max_longitudinal_retarded_circ_bias': mx(long_ret, 'response_circulation_directional_bias'),
        'forward_reverse_pair_count': len(pair_rows),
        'forward_qp_better_than_reverse_count': forward_better_count,
        'forward_reverse_halfplane_opposed_count': halfplane_opposed_count,
        'best_retarded_row': slim(best_ret),
        'first_boundary_pass_row': slim(best_boundary),
    }


def slim(r: Optional[dict]) -> dict:
    if not r:
        return {}
    keys = ['bridge_mode','response_op','lift_mode','causal_class','bridge_axis_class','QP_transfer_mean_resid','projected_J2_plus_I_resid','projected_QP_opnorm','projected_QP_normality_resid','halfplane_bias_pos_minus_neg','response_circulation_directional_bias','signed_mean_response_circulation','retarded_boundary_value_gate_pass','boundary_value_signal_pass','longitudinal_retarded_signal_pass','context']
    return {k:r.get(k) for k in keys}


def table_rows(summaries: List[dict]) -> str:
    lines=[
        '| variant | beta | assemblies | pair_harm | Q_harm | P_harm | good rows | long-ret rows | causal signal | boundary signal | retarded gate | advanced boundary | best long-ret QP | max halfplane bias | fwd/rev pairs | fwd better | used Δβ? |',
        '|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|'
    ]
    for r in summaries:
        a=r['auto_metrics']; dm=r['directed_metrics']; am=r['alignment_metrics']; s=r['retarded_boundary_summary']
        lines.append(f"| {r['variant']} | ({a['beta0']},{a['beta1']},{a['beta2']},{a['beta3']}) | {r['assemblies_applied']} | {dm['pair_transport_harmonic_ratio']:.6g} | {am['Q_even_harmonic_ratio']:.6g} | {am['P_odd_harmonic_ratio']:.6g} | {s['good_bounded_nonartifact_active_count']} | {s['good_longitudinal_retarded_count']} | {s['causal_order_signal_count']} | {s['boundary_value_signal_count']} | {s['retarded_boundary_value_gate_pass_count']} | {s['advanced_boundary_signal_count']} | {s['best_longitudinal_retarded_QP_resid']:.6g} | {s['max_longitudinal_retarded_halfplane_abs_bias']:.6g} | {s['forward_reverse_pair_count']} | {s['forward_qp_better_than_reverse_count']} | {r['decision_used_delta_beta_any']} |")
    return '\n'.join(lines)


def docs(summaries: List[dict]) -> Tuple[str,str,str,str]:
    table=table_rows(summaries)
    short = {r['variant']: r['retarded_boundary_summary'] for r in summaries}
    smd=f"""# SUMMARY — retarded boundary-value selection gate

This package shifts the question from geometric `J`/spin signatures to a pre-causal boundary-value audit.  The hypothesis is that the sequential birth order and the longitudinal/radial growth axis may already behave like a retarded/advanced selection layer.

Primary gate: bounded, well-conditioned, non-artifact **forward/retarded** response on a longitudinal bridge with stable half-plane bias and no reverse-cycle symmetrization.  Q/P-lock and `J²=-I` are diagnostics, not primary requirements.

{table}
"""
    rmd=f"""# RESULTS — retarded boundary-value selection gate

{table}

## Interpretation rule

The test does **not** search for spin or `J`.  It asks whether the growth semi-order itself creates a robust boundary-value/retarded selection analogous to choosing an upper/lower half-plane prescription.

A positive row must be:

- well-conditioned, bounded, non-artifact, and active;
- forward/retarded rather than reverse/advanced or symmetrized;
- longitudinal/radial rather than purely transverse;
- directionally biased by sibling response circulation;
- half-plane biased on the Q/P projection;
- not dependent on `delta_beta` as a decision input.

## Compact JSON

```json
{json.dumps(short, indent=2)}
```

## Conservative reading

If the retarded/boundary signal is nonzero but does not order Q/P, then the result is not `i`; it is only a candidate pre-causal boundary-value layer.  This directly tracks the distinction raised by the logarithm/branch-cut discussion: the hidden convention may live in a retarded/advanced prescription rather than in an explicit geometric orientation.
"""
    audit="""# SOURCE AUDIT

This package is derived from `test_bounded_response_bridge_sanity_gate.py` and reuses its bounded bridge rows, but changes the primary question.

It does not:

- set `i`, global `J`, Hodge star, physical adjoint, positivity, C*-norm, or a target spin signature;
- use `delta_beta`/H2 as a selection input;
- count Q/P-lock or `J²=-I` as primary success.

It does:

- distinguish retarded/forward response operators from advanced/reverse and symmetrized cycle operators;
- distinguish longitudinal/radial bridge modes from transverse/birth-q bridge modes;
- audit bounded contraction, bridge conditioning, half-plane bias, and response-circulation directional bias;
- keep strict-symmetrized control as a null test.

Model label: CQNM/s=-1 inspired growing real complement-network diagnostic.  The growth order / longitudinal axis is treated as a pre-causal provenance structure, not as physical time.
"""
    readme="""# Retarded boundary-value selection gate

Run:

```bash
python3 test_retarded_boundary_value_selection_gate.py
```

Outputs are written to `retarded_boundary_value_selection_out_L2/`.
"""
    return smd,rmd,audit,readme


def package(out: Path, zip_path: Path) -> None:
    files=[
        Path(__file__).name,
        'test_bounded_response_bridge_sanity_gate.py',
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
    with zipfile.ZipFile(zip_path,'w',zipfile.ZIP_DEFLATED) as z:
        for f in files:
            if Path(f).exists(): z.write(f,f)
        for p in sorted(out.rglob('*')):
            if p.is_file(): z.write(p,p.resolve().relative_to(Path.cwd()))


def main() -> None:
    ap=argparse.ArgumentParser()
    # Core/growth arguments compatible with p69/p79 stack.
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
    # Bridge/response args inherited from p79.
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
    # Retarded-boundary thresholds.
    ap.add_argument('--circ-sign-tol', type=float, default=1e-10)
    ap.add_argument('--circ-bias-threshold', type=float, default=0.80)
    ap.add_argument('--halfplane-bias-threshold', type=float, default=0.50)
    ap.add_argument('--retarded-opnorm-threshold', type=float, default=1.05)
    ap.add_argument('--retarded-normality-threshold', type=float, default=0.75)
    ap.add_argument('--loose-qp-leakage-threshold', type=float, default=0.95)
    ap.add_argument('--out', default='retarded_boundary_value_selection_out_L2')
    ap.add_argument('--zip', default='cnna_retarded_boundary_value_selection_gate_pkg_L2.zip')
    args=ap.parse_args()
    p70.patch_modules()
    out=Path(args.out)
    if out.exists(): shutil.rmtree(out)
    out.mkdir(parents=True)
    summaries=[run_variant(v,args,out) for v in args.variants]
    smd,rmd,audit,readme=docs(summaries)
    (out/'SUMMARY.md').write_text(smd,encoding='utf-8')
    (out/'RESULTS.md').write_text(rmd,encoding='utf-8')
    (out/'SOURCE_AUDIT.md').write_text(audit,encoding='utf-8')
    (out/'README.md').write_text(readme,encoding='utf-8')
    write_csv(out/'comparative_retarded_boundary_value_selection_summary.csv', [
        {
            'variant': r['variant'],
            'beta0': r['auto_metrics']['beta0'], 'beta1': r['auto_metrics']['beta1'], 'beta2': r['auto_metrics']['beta2'], 'beta3': r['auto_metrics']['beta3'],
            'assemblies_applied': r['assemblies_applied'],
            'pair_harm': r['directed_metrics']['pair_transport_harmonic_ratio'],
            'Q_harm': r['alignment_metrics']['Q_even_harmonic_ratio'],
            'P_harm': r['alignment_metrics']['P_odd_harmonic_ratio'],
            **r['retarded_boundary_summary'],
            'decision_used_delta_beta_any': r['decision_used_delta_beta_any'],
        } for r in summaries
    ])
    (out/'summary.json').write_text(json.dumps(summaries,indent=2),encoding='utf-8')
    package(out, Path(args.zip))
    print(f'wrote {args.zip}')


if __name__ == '__main__':
    main()
