#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import math
import shutil
import zipfile
import csv
from pathlib import Path
from typing import List, Tuple

import numpy as np

import test_dual_pairing_assembly_growth_rule_gate as p68
import test_dual_assembly_order_context_ablation_gate as p69

EPS = 1e-12
ORIG_P69_SUMMARIZE_SELECTION = p69.summarize_selection


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
    return p68.fval(row, key, default)


def fbool(x) -> bool:
    return p68.fbool(x)


def qpb(row: dict) -> float:
    return p68.qp_balance(row)


def bounded01(x: float) -> float:
    return max(0.0, min(1.0, float(x)))


def signed_jlock_components(row: dict) -> dict:
    c_mean = fval(row, 'best_C_eigen_J_lock_mean_resid', 9.0)
    c_max = fval(row, 'best_C_eigen_J_lock_max_resid', 9.0)
    flip = fval(row, 'kappa_signed_flip_abs', 9.0)
    kappa_amp = fval(row, 'kappa_signed_amp_min', 0.0)
    signed_id = abs(fval(row, 'comm_signed_birth_over_abs', 0.0))
    signed_k = abs(fval(row, 'kappa_comm_signed_birth_over_abs', 0.0))
    signed_min = min(signed_id, signed_k) if (signed_id > 0 or signed_k > 0) else 0.0
    signed_max = max(signed_id, signed_k)
    c_quality = bounded01((0.65 - c_mean) / 0.65)
    cmax_quality = bounded01((0.95 - c_max) / 0.95)
    flip_quality = bounded01((0.55 - flip) / 0.55)
    amp_quality = bounded01(max(kappa_amp, signed_min) / 0.25)
    signed_quality = bounded01(signed_max / 0.35)
    qp_quality = qpb(row)
    # Joint term is deliberately multiplicative: B-role only scores strongly when C-lock,
    # signed amplitude, and kappa flip are jointly present.
    joint = (c_quality * (0.60*cmax_quality + 0.40) * (0.55*flip_quality + 0.45) * (0.60*amp_quality + 0.40))
    return {
        'c_mean': c_mean,
        'c_max': c_max,
        'flip_abs': flip,
        'kappa_amp_min': kappa_amp,
        'signed_identity_abs': signed_id,
        'signed_kappa_abs': signed_k,
        'signed_min_abs': signed_min,
        'signed_max_abs': signed_max,
        'c_quality': c_quality,
        'cmax_quality': cmax_quality,
        'flip_quality': flip_quality,
        'amp_quality': amp_quality,
        'signed_quality': signed_quality,
        'qp_quality': qp_quality,
        'joint_quality': joint,
    }


def signed_jlock_proxy_B_score(row: dict) -> float:
    c = signed_jlock_components(row)
    directed = bounded01(fval(row, 'directed_imbalance', 0.0))
    transverse = bounded01(fval(row, 'transverse_complementarity', 0.0))
    # C-lock + signed amplitude are primary; kappa flip is important but not allowed to
    # erase the amplitude term.  No beta/H2/harmonic quantity is used here.
    return (
        3.50 * c['joint_quality']
        + 1.00 * c['c_quality']
        + 0.85 * c['amp_quality']
        + 0.55 * c['flip_quality']
        + 0.35 * c['qp_quality']
        + 0.15 * directed
        + 0.10 * transverse
        - 0.25 * min(2.0, c['c_mean'])
    )


def summarize_selection_with_signed(pairing_log: List[dict]) -> dict:
    base = ORIG_P69_SUMMARIZE_SELECTION(pairing_log)
    b_rows = [r for r in pairing_log if 'B_' in str(r.get('assembly_role','')) and fbool(r.get('applied'))]
    def avg(key: str) -> float:
        vals = [fval(r, key, math.nan) for r in b_rows if str(r.get(key,'')).strip() != '']
        vals = [v for v in vals if math.isfinite(v)]
        return float(np.mean(vals)) if vals else 0.0
    base.update({
        'B_signed_identity_abs_avg': avg('B_signed_identity_abs_eval'),
        'B_signed_kappa_abs_avg': avg('B_signed_kappa_abs_eval'),
        'B_signed_min_abs_avg': avg('B_signed_min_abs_eval'),
        'B_joint_quality_avg': avg('B_joint_quality_eval'),
        'B_c_quality_avg': avg('B_c_quality_eval'),
        'B_flip_quality_avg': avg('B_flip_quality_eval'),
        'B_amp_quality_avg': avg('B_amp_quality_eval'),
    })
    return base


def annotate_pair_logs(out: Path) -> None:
    # p69 logs selected eval metrics for B, but not the new components.  Annotate logs
    # after-the-fact by matching candidate_eval_rows for readability.
    for cand_path in out.rglob('candidate_eval_rows.csv'):
        cand_rows = []
        with cand_path.open(newline='', encoding='utf-8') as f:
            for r in csv.DictReader(f):
                key = (r.get('candidate_id',''), r.get('face_a',''), r.get('face_b',''), r.get('move_class',''), r.get('stage',''))
                cc = signed_jlock_components(r)
                r.update({f'newB_{k}': v for k,v in cc.items()})
                cand_rows.append(r)
        write_csv(cand_path, cand_rows)


def patch_modules() -> None:
    # Monkey-patch the B-role scoring used by p68.choose_pair_B and by p69 wrappers.
    p68.proxy_B_score = signed_jlock_proxy_B_score
    p69.p68.proxy_B_score = signed_jlock_proxy_B_score

    old_apply = p69.apply_pair
    def apply_pair_annotated(model, K, orig, args, event_t: int, cascade_index: int, role: str, rule: str):
        K2, log, applied = old_apply(model, K, orig, args, event_t, cascade_index, role, rule)
        # eval components are written below by choose wrappers; keep hook simple.
        return K2, log, applied
    p69.apply_pair = apply_pair_annotated
    p69.summarize_selection = summarize_selection_with_signed


def enrich_variant_summaries(out: Path) -> None:
    for vjson in out.rglob('variant_option_summary.json'):
        data = json.loads(vjson.read_text(encoding='utf-8'))
        # Add selected B component averages from the pairing log by matching B rows.
        log_path = vjson.parent / 'assembly_pairing_log.csv'
        if log_path.exists():
            rows = list(csv.DictReader(log_path.open(newline='', encoding='utf-8')))
            # If component columns are absent in older log, approximate from selected C/kappa only.
            b = [r for r in rows if 'B_' in str(r.get('assembly_role','')) and fbool(r.get('applied'))]
            comps = []
            for r in b:
                c = {
                    'c_mean': fval(r,'C_eigen_J_lock_mean_resid_eval',9.0),
                    'flip_abs': fval(r,'kappa_signed_flip_abs_eval',9.0),
                }
                comps.append(c)
            if comps:
                data.setdefault('signed_Jlock_role_coupling_metrics', {})
                data['signed_Jlock_role_coupling_metrics'].update({
                    'selected_B_count': len(comps),
                    'selected_B_C_lock_mean_avg': float(np.mean([x['c_mean'] for x in comps])),
                    'selected_B_kappa_flip_abs_avg': float(np.mean([x['flip_abs'] for x in comps])),
                })
        vjson.write_text(json.dumps(data, indent=2), encoding='utf-8')


def collect_comparative(out: Path, rows: List[dict]) -> None:
    flat = []
    for r in rows:
        a = r['auto_metrics']; dm = r['directed_metrics']; sm = r['signed_quadrature_metrics']; am = r['alignment_metrics']; sel = r['selection_metrics']
        flat.append({
            'option': r['option'],
            'variant': r['variant'],
            'beta0': a['beta0'], 'beta1': a['beta1'], 'beta2': a['beta2'], 'beta3': a['beta3'],
            'pairings': r['automatic_pairings_applied'],
            'assemblies': r['assemblies_applied'],
            'pair_harm': dm['pair_transport_harmonic_ratio'],
            'Q_harm': am['Q_even_harmonic_ratio'],
            'P_harm': am['P_odd_harmonic_ratio'],
            'J_lock': am['best_per_pair_mean_J_lock_resid'],
            'signed_birth': sm['signed_birth_over_abs_sum_ratio'],
            'selected_C_lock_avg': sel.get('selected_C_eigen_J_lock_mean_avg',0.0),
            'selected_kappa_flip_abs_avg': sel.get('selected_kappa_flip_abs_avg',0.0),
            'B_joint_quality_avg': sel.get('B_joint_quality_avg',0.0),
            'B_signed_min_abs_avg': sel.get('B_signed_min_abs_avg',0.0),
            'used_delta_beta': r['decision_used_delta_beta_any'],
        })
    write_csv(out / 'comparative_signed_Jlock_role_coupling_summary.csv', flat)


def slim(r: dict) -> dict:
    return {
        'option': r['option'], 'variant': r['variant'],
        'beta': [r['auto_metrics'][f'beta{i}'] for i in range(4)],
        'pairings': r['automatic_pairings_applied'],
        'assemblies': r['assemblies_applied'],
        'pair_harm': r['directed_metrics']['pair_transport_harmonic_ratio'],
        'Q_harm': r['alignment_metrics']['Q_even_harmonic_ratio'],
        'P_harm': r['alignment_metrics']['P_odd_harmonic_ratio'],
        'J_lock': r['alignment_metrics']['best_per_pair_mean_J_lock_resid'],
        'signed_birth': r['signed_quadrature_metrics']['signed_birth_over_abs_sum_ratio'],
        'selected_C_lock_avg': r['selection_metrics'].get('selected_C_eigen_J_lock_mean_avg',0.0),
        'selected_kappa_flip_avg': r['selection_metrics'].get('selected_kappa_flip_abs_avg',0.0),
        'used_delta_beta': r['decision_used_delta_beta_any'],
    }


def make_docs(summary: dict) -> tuple[str,str,str,str]:
    rows = summary['variant_rows']
    lines = ['| option | variant | beta | pairs | asm | pair harm | Q harm | P harm | J-lock | signed birth | selected C | selected kappa | used dBeta? |',
             '|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|']
    for r in rows:
        s = slim(r)
        lines.append(f"| {s['option']} | {s['variant']} | ({s['beta'][0]},{s['beta'][1]},{s['beta'][2]},{s['beta'][3]}) | {s['pairings']} | {s['assemblies']} | {s['pair_harm']:.6g} | {s['Q_harm']:.6g} | {s['P_harm']:.6g} | {s['J_lock']:.6g} | {s['signed_birth']:.6g} | {s['selected_C_lock_avg']:.6g} | {s['selected_kappa_flip_avg']:.6g} | {s['used_delta_beta']} |")
    table = '\n'.join(lines)
    best_signed = max([r for r in rows if r['variant']!='strict_symmetrized_control'], key=lambda r: abs(r['signed_quadrature_metrics']['signed_birth_over_abs_sum_ratio']), default=None)
    best_j = min([r for r in rows if r['variant']!='strict_symmetrized_control'], key=lambda r: r['alignment_metrics']['best_per_pair_mean_J_lock_resid'], default=None)
    best_joint = None
    nonstrict = [r for r in rows if r['variant']!='strict_symmetrized_control']
    if nonstrict:
        best_joint = max(nonstrict, key=lambda r: abs(r['signed_quadrature_metrics']['signed_birth_over_abs_sum_ratio']) / (r['alignment_metrics']['best_per_pair_mean_J_lock_resid'] + EPS))
    best = {'best_signed': slim(best_signed) if best_signed else {}, 'best_J_lock': slim(best_j) if best_j else {}, 'best_signed_over_J': slim(best_joint) if best_joint else {}}
    smd = f"""# SUMMARY — signed J-lock role-coupling gate

Model label:
CNNA growing primal simplicial complex with deterministic sequential provenance growth, directed antisymmetric birth-transport operators, and a dynamic two-pair assembly rule whose B-role couples C-lock with signed kappa/birth amplitude.

This test keeps the A-role as provenance/QP-carrier proxy and changes only the B-role scoring:

```text
B role = low C-eigen J-lock residual + signed birth/kappa amplitude + kappa flip,
         with shared face/edge context to A.
```

No decision uses delta_beta, H2, harmonic projections, complex scalars, Hodge, positivity, physical adjoint, or final sym(M).

{table}

## Best rows

```json
{json.dumps(best, indent=2)}
```
"""
    rmd = f"""# RESULTS — signed J-lock role-coupling gate

## Comparative table

{table}

## Gate criterion

A constructive result would require the same non-strict row to show:

```text
beta2 open,
Q/P harmonic positive,
signed_birth strong,
J-lock low,
strict_sym killed,
decision_used_delta_beta_any = false.
```

The previous order/context ablation showed that signed orientation and J-lock separate.  This test asks whether explicitly coupling those two quantities in the B-role can make them coincide dynamically.
"""
    audit = """# SOURCE AUDIT

Carried forward:

- Pair-property tradeoff gate: beta2/QP, C-lock, and kappa-flip split across single-pair candidates.
- Dual-pair audit/growth/order-context gates: two-pair assemblies can open beta2 and Q/P strongly, but signed orientation and J-lock do not stabilize together.

This package changes only the B-role score.  It does not set i, J, Hodge, *, positivity, C*-norm, or final sym(M).  Delta-beta/H2/harmonic values are measured after the fact and are not used in the B-role score.

Caveat: the B-role uses signed-birth/kappa diagnostics as a selector.  This is an intentional diagnostic gate, not yet a claimed canonical CNNA rule.  If it succeeds, the selected signed/J-lock coupling still requires a derived-only formal justification.
"""
    readme = """# Signed J-lock role-coupling gate

Run:

```bash
python3 test_signed_Jlock_role_coupling_gate.py
```

The script monkey-patches only the B-role score of the existing dynamic assembly engine and writes a full ZIP package with summaries and logs.
"""
    return smd, rmd, audit, readme


def package(out: Path, zip_path: Path) -> None:
    files = [
        Path(__file__).name,
        'test_dual_assembly_order_context_ablation_gate.py',
        'test_dual_pairing_assembly_growth_rule_gate.py',
        'test_dual_pairing_two_edge_assembly_gate.py',
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
    # Use p69 defaults, but reduced to the relevant dynamic A->B gate unless explicitly changed.
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
    ap.add_argument('--context-modes', nargs='*', default=['connected','strong'])
    ap.add_argument('--reuse-modes', nargs='*', default=['reuseB','noReuseB'])
    ap.add_argument('--variants', nargs='*', default=['real_growth','strict_symmetrized_control','no_backreaction'])
    ap.add_argument('--phase-sign', type=int, default=1)
    ap.add_argument('--out', default='signed_Jlock_role_coupling_out_L2')
    ap.add_argument('--zip', default='cnna_signed_Jlock_role_coupling_gate_pkg_L2.zip')
    args = ap.parse_args()

    patch_modules()
    out = Path(args.out)
    if out.exists():
        shutil.rmtree(out)
    out.mkdir(parents=True)
    rows = p69.run_experiment(args, out)
    # Keep p69 comparative too, then write our comparison/docs.
    try:
        p69.write_comparative(out, rows)
    except Exception:
        pass
    summary = {'args': vars(args), 'variant_rows': rows}
    summary['best_summary'] = p69.summarize_best(rows)
    collect_comparative(out, rows)
    (out / 'comparative_summary.json').write_text(json.dumps(summary, indent=2), encoding='utf-8')
    smd, rmd, audit, readme = make_docs(summary)
    (out / 'SUMMARY.md').write_text(smd, encoding='utf-8')
    (out / 'RESULTS.md').write_text(rmd, encoding='utf-8')
    (out / 'SOURCE_AUDIT.md').write_text(audit, encoding='utf-8')
    (out / 'README.md').write_text(readme, encoding='utf-8')
    package(out, Path(args.zip))
    print(json.dumps({
        'zip': args.zip,
        'out': args.out,
        'rows': [slim(r) for r in rows],
        'best_summary': summary['best_summary'],
    }, indent=2))


if __name__ == '__main__':
    main()
