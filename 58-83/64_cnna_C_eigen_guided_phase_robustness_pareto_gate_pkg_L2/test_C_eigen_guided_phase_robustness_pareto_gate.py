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

EPS = 1e-12


def read_csv(path: Path) -> List[dict]:
    if not path.exists() or path.stat().st_size == 0:
        return []
    with path.open(newline='', encoding='utf-8') as f:
        return list(csv.DictReader(f))


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


def candidate_key(r: dict) -> Tuple[str, str, str, str, str, str, str, str]:
    # A candidate must match across phase signs without using beta/H2 as a selector.
    return (
        str(r.get('variant','')),
        str(r.get('rule','')),
        str(r.get('event_t','')),
        str(r.get('scan_id','')),
        str(r.get('cascade_index','')),
        str(r.get('candidate_id','')),
        str(r.get('face_a','')),
        str(r.get('face_b','')),
    )


def load_candidate_rows(input_dir: Path) -> List[dict]:
    rows: List[dict] = []
    for p in sorted(input_dir.glob('*_phaseplus1/candidate_eval_rows.csv')) + sorted(input_dir.glob('*_phaseminus1/candidate_eval_rows.csv')):
        for r in read_csv(p):
            r = dict(r)
            r['source_candidate_file'] = str(p.relative_to(input_dir))
            rows.append(r)
    return rows


def match_phase_rows(rows: List[dict], args: argparse.Namespace) -> Tuple[List[dict], List[dict]]:
    groups: Dict[Tuple[str,...], Dict[int, dict]] = {}
    unmatched: List[dict] = []
    for r in rows:
        try:
            ph = int(float(r.get('phase_sign', 0) or 0))
        except Exception:
            ph = 0
        if ph not in {1, -1}:
            unmatched.append(r)
            continue
        groups.setdefault(candidate_key(r), {})[ph] = r
    out: List[dict] = []
    for key, d in groups.items():
        if 1 not in d or -1 not in d:
            continue
        p, m = d[1], d[-1]
        c_plus = fval(p, 'best_C_eigen_J_lock_mean_resid', 1.0)
        c_minus = fval(m, 'best_C_eigen_J_lock_mean_resid', 1.0)
        cmax_plus = fval(p, 'best_C_eigen_J_lock_max_resid', 1.0)
        cmax_minus = fval(m, 'best_C_eigen_J_lock_max_resid', 1.0)
        s_plus = fval(p, 'comm_signed_birth_over_abs', 0.0)
        s_minus = fval(m, 'comm_signed_birth_over_abs', 0.0)
        abs_plus = abs(s_plus)
        abs_minus = abs(s_minus)
        # zero is a perfect sign flip of comparable amplitude; one means same sign.
        flip_score = (s_plus + s_minus) / (abs_plus + abs_minus + EPS)
        # amplitude balance zero is perfect equal magnitude; one is very unbalanced.
        amp_balance = abs(abs_plus - abs_minus) / (abs_plus + abs_minus + EPS)
        # bounded quality score, not used by any growth rule; audit only.
        phase_flip_quality = (1.0 - min(1.0, abs(flip_score))) * (1.0 - min(1.0, amp_balance))
        q_plus, p_plus = fval(p, 'Q_norm', 0.0), fval(p, 'P_norm', 0.0)
        q_minus, p_minus = fval(m, 'Q_norm', 0.0), fval(m, 'P_norm', 0.0)
        # Candidate-level harmonic cannot be known without applying the move.
        # This proxy only asks whether both quadrature channels are present before applying the candidate.
        qp_proxy_plus = min(q_plus, p_plus) / (max(q_plus, p_plus) + EPS)
        qp_proxy_minus = min(q_minus, p_minus) / (max(q_minus, p_minus) + EPS)
        qp_proxy = min(qp_proxy_plus, qp_proxy_minus)
        comm_abs_plus = fval(p, 'comm_abs_area', 0.0)
        comm_abs_minus = fval(m, 'comm_abs_area', 0.0)
        comm_abs_proxy = min(comm_abs_plus, comm_abs_minus) / (max(comm_abs_plus, comm_abs_minus) + EPS) if max(comm_abs_plus, comm_abs_minus) > EPS else 0.0
        # Pareto objective signs: minimize C-lock and flip_abs, maximize signed amplitude, Q/P balance, provenance scores.
        row = {
            'variant': p.get('variant',''),
            'rule': p.get('rule',''),
            'event_t': p.get('event_t',''),
            'scan_id': p.get('scan_id',''),
            'cascade_index': p.get('cascade_index',''),
            'candidate_id': p.get('candidate_id',''),
            'move_class': p.get('move_class',''),
            'face_a': p.get('face_a',''),
            'face_b': p.get('face_b',''),
            'A_gate_plus': p.get('A_gate',''),
            'A_gate_minus': m.get('A_gate',''),
            'beta2_opening_audit_plus': p.get('beta2_opening_audit_only',''),
            'beta2_opening_audit_minus': m.get('beta2_opening_audit_only',''),
            'delta_beta2_audit_plus': p.get('delta_beta2_audit_only',''),
            'delta_beta2_audit_minus': m.get('delta_beta2_audit_only',''),
            'C_lock_mean_plus': c_plus,
            'C_lock_mean_minus': c_minus,
            'C_lock_mean_worst': max(c_plus, c_minus),
            'C_lock_mean_avg': 0.5 * (c_plus + c_minus),
            'C_lock_max_plus': cmax_plus,
            'C_lock_max_minus': cmax_minus,
            'C_lock_max_worst': max(cmax_plus, cmax_minus),
            'signed_birth_plus': s_plus,
            'signed_birth_minus': s_minus,
            'signed_flip_score_zero_if_perfect': flip_score,
            'signed_flip_abs': abs(flip_score),
            'signed_amplitude_min': min(abs_plus, abs_minus),
            'signed_amplitude_avg': 0.5 * (abs_plus + abs_minus),
            'signed_amplitude_balance_zero_if_equal': amp_balance,
            'phase_flip_quality_0_1': phase_flip_quality,
            'Q_norm_plus': q_plus,
            'P_norm_plus': p_plus,
            'Q_norm_minus': q_minus,
            'P_norm_minus': p_minus,
            'QP_balance_proxy_min': qp_proxy,
            'comm_abs_area_plus': comm_abs_plus,
            'comm_abs_area_minus': comm_abs_minus,
            'comm_abs_area_balance_proxy': comm_abs_proxy,
            'transport_cosine_plus': fval(p, 'transport_cosine_ka_kb_reversed', 0.0),
            'transport_cosine_minus': fval(m, 'transport_cosine_ka_kb_reversed', 0.0),
            'A_rank_score_plus': fval(p, 'A_rank_score', 0.0),
            'A_rank_score_minus': fval(m, 'A_rank_score', 0.0),
            'A_rank_score_avg': 0.5 * (fval(p, 'A_rank_score', 0.0) + fval(m, 'A_rank_score', 0.0)),
            'directed_imbalance_plus': fval(p, 'directed_imbalance', 0.0),
            'directed_imbalance_minus': fval(m, 'directed_imbalance', 0.0),
            'directed_imbalance_avg': 0.5 * (fval(p, 'directed_imbalance', 0.0) + fval(m, 'directed_imbalance', 0.0)),
            'transverse_complementarity_plus': fval(p, 'transverse_complementarity', 0.0),
            'transverse_complementarity_minus': fval(m, 'transverse_complementarity', 0.0),
            'transverse_complementarity_avg': 0.5 * (fval(p, 'transverse_complementarity', 0.0) + fval(m, 'transverse_complementarity', 0.0)),
            'selected_plus': any(str(k).startswith('selected_by_') and str(v).lower() == 'true' for k,v in p.items()),
            'selected_minus': any(str(k).startswith('selected_by_') and str(v).lower() == 'true' for k,v in m.items()),
            'decision_used_delta_beta': False,
        }
        row['passes_C_lock_worst'] = row['C_lock_mean_worst'] < args.lock_residual_threshold and row['C_lock_max_worst'] < args.lock_max_threshold
        row['passes_phase_flip'] = row['signed_flip_abs'] < args.flip_abs_threshold and row['signed_amplitude_min'] > args.signed_amp_min
        row['passes_QP_proxy'] = row['QP_balance_proxy_min'] > args.qp_balance_min and min(comm_abs_plus, comm_abs_minus) > args.comm_abs_min
        row['passes_A_gate_both'] = bval(row, 'A_gate_plus') and bval(row, 'A_gate_minus')
        row['passes_beta2_audit_both'] = bval(row, 'beta2_opening_audit_plus') and bval(row, 'beta2_opening_audit_minus')
        row['passes_all_pareto_gate'] = bool(row['passes_C_lock_worst'] and row['passes_phase_flip'] and row['passes_QP_proxy'] and row['passes_A_gate_both'])
        # Compact scalar for ranking; not a selection rule for growth.
        row['pareto_audit_score'] = (
            (1.0 / (row['C_lock_mean_worst'] + 0.05))
            * (1.0 - min(1.0, row['signed_flip_abs']))
            * row['signed_amplitude_min']
            * row['QP_balance_proxy_min']
            * (0.25 + row['directed_imbalance_avg'])
            * (0.25 + row['transverse_complementarity_avg'])
        )
        out.append(row)
    return out, unmatched


def dominates(a: dict, b: dict) -> bool:
    # Pareto dimensions: lower C-lock-worst, lower flip_abs, higher signed amplitude, higher Q/P proxy, higher directed/transverse.
    dims = [
        ('C_lock_mean_worst', -1),
        ('signed_flip_abs', -1),
        ('signed_amplitude_min', 1),
        ('QP_balance_proxy_min', 1),
        ('directed_imbalance_avg', 1),
        ('transverse_complementarity_avg', 1),
    ]
    better_or_equal = True
    strictly_better = False
    for k, sign in dims:
        av, bv = fval(a,k), fval(b,k)
        if sign == 1:
            if av + 1e-12 < bv:
                better_or_equal = False
                break
            if av > bv + 1e-12:
                strictly_better = True
        else:
            if av > bv + 1e-12:
                better_or_equal = False
                break
            if av + 1e-12 < bv:
                strictly_better = True
    return better_or_equal and strictly_better


def pareto_front(rows: List[dict]) -> List[dict]:
    front = []
    for i, r in enumerate(rows):
        if not bval(r, 'passes_A_gate_both'):
            continue
        if any(dominates(o, r) for j, o in enumerate(rows) if j != i and bval(o, 'passes_A_gate_both')):
            continue
        rr = dict(r)
        rr['pareto_front'] = True
        front.append(rr)
    front.sort(key=lambda r: fval(r, 'pareto_audit_score'), reverse=True)
    return front


def summarize(rows: List[dict], front: List[dict]) -> Tuple[List[dict], List[dict]]:
    grouped: Dict[Tuple[str,str], List[dict]] = {}
    for r in rows:
        grouped.setdefault((r['variant'], r['rule']), []).append(r)
    summaries = []
    top = []
    for (variant, rule), rs in sorted(grouped.items()):
        ag = [r for r in rs if bval(r, 'passes_A_gate_both')]
        beta = [r for r in ag if bval(r, 'passes_beta2_audit_both')]
        cpass = [r for r in ag if bval(r, 'passes_C_lock_worst')]
        fpass = [r for r in ag if bval(r, 'passes_phase_flip')]
        allp = [r for r in ag if bval(r, 'passes_all_pareto_gate')]
        fr = [r for r in front if r['variant'] == variant and r['rule'] == rule]
        best_score = max(ag, key=lambda r: fval(r,'pareto_audit_score'), default=None)
        best_c = min(ag, key=lambda r: fval(r,'C_lock_mean_worst',1.0), default=None)
        best_flip = min(ag, key=lambda r: fval(r,'signed_flip_abs',1.0), default=None)
        best_all = min(allp, key=lambda r: fval(r,'C_lock_mean_worst',1.0), default=None)
        def g(r,k,d=0.0): return fval(r,k,d) if r else d
        summaries.append({
            'variant': variant,
            'rule': rule,
            'matched_candidate_count': len(rs),
            'A_gated_both_count': len(ag),
            'beta2_opening_both_audit_count': len(beta),
            'C_lock_pass_count': len(cpass),
            'phase_flip_pass_count': len(fpass),
            'all_pareto_gate_pass_count': len(allp),
            'pareto_front_count': len(fr),
            'best_score_candidate_id': best_score.get('candidate_id','') if best_score else '',
            'best_score': g(best_score,'pareto_audit_score'),
            'best_score_C_lock_worst': g(best_score,'C_lock_mean_worst'),
            'best_score_flip_abs': g(best_score,'signed_flip_abs'),
            'best_score_signed_amp_min': g(best_score,'signed_amplitude_min'),
            'best_score_QP_proxy': g(best_score,'QP_balance_proxy_min'),
            'best_C_candidate_id': best_c.get('candidate_id','') if best_c else '',
            'best_C_lock_worst': g(best_c,'C_lock_mean_worst'),
            'best_C_flip_abs': g(best_c,'signed_flip_abs'),
            'best_flip_candidate_id': best_flip.get('candidate_id','') if best_flip else '',
            'best_flip_abs': g(best_flip,'signed_flip_abs'),
            'best_flip_C_lock_worst': g(best_flip,'C_lock_mean_worst'),
            'best_all_gate_candidate_id': best_all.get('candidate_id','') if best_all else '',
            'best_all_gate_C_lock_worst': g(best_all,'C_lock_mean_worst'),
            'best_all_gate_flip_abs': g(best_all,'signed_flip_abs'),
            'decision_used_delta_beta_any': False,
        })
        for r in sorted(ag, key=lambda x: fval(x,'pareto_audit_score'), reverse=True)[:10]:
            tr = dict(r)
            tr['top_rank_by_pareto_audit_score'] = len([x for x in top if x.get('variant') == variant and x.get('rule') == rule]) + 1
            top.append(tr)
    return summaries, top


def make_docs(summary_rows: List[dict], top_rows: List[dict], args: argparse.Namespace) -> Tuple[str,str,str,str]:
    lines = ['| variant/rule | matched | A-gated | beta2 audit | C-pass | flip-pass | all-pass | best score | best C-lock | best flip |', '|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|']
    for r in summary_rows:
        lines.append(f"| {r['variant']} / {r['rule']} | {r['matched_candidate_count']} | {r['A_gated_both_count']} | {r['beta2_opening_both_audit_count']} | {r['C_lock_pass_count']} | {r['phase_flip_pass_count']} | {r['all_pareto_gate_pass_count']} | {r['best_score']:.6g} | {r['best_C_lock_worst']:.6g} | {r['best_flip_abs']:.6g} |")
    table = '\n'.join(lines)
    tops = ['| variant/rule | cand | C-lock worst | flip abs | signed amp min | Q/P proxy | A-rank avg | directed avg | transv avg | beta2 audit | selected? |','|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|']
    for r in top_rows[:20]:
        tops.append(f"| {r['variant']} / {r['rule']} | {r['candidate_id']} | {float(r['C_lock_mean_worst']):.6g} | {float(r['signed_flip_abs']):.6g} | {float(r['signed_amplitude_min']):.6g} | {float(r['QP_balance_proxy_min']):.6g} | {float(r['A_rank_score_avg']):.6g} | {float(r['directed_imbalance_avg']):.6g} | {float(r['transverse_complementarity_avg']):.6g} | {r['passes_beta2_audit_both']} | {r['selected_plus'] or r['selected_minus']} |")
    top_table = '\n'.join(tops)
    smd = f"""# SUMMARY — C-eigen guided phase robustness Pareto gate

Model label:
CNNA growing primal simplicial complex with deterministic sequential provenance growth, A-gated nonlinear complement-pair candidate space, directed antisymmetric birth-transport operators, and local C/J pair algebra.

This package does **not** apply a new growth rule. It audits the already enumerated A-gated candidate space from the C-eigen-guided pairing package and matches the same candidates across phase_sign +1/-1.

Anti-smuggling rule:
`delta_beta2` is included only as an audit column. It is not used in the Pareto score and not used as a decision input.

## Comparative summary

{table}

## Top candidates by Pareto audit score

{top_table}
"""
    rmd = f"""# RESULTS — C-eigen guided phase robustness Pareto gate

## Gate being tested

The gate asks whether there are A-gated candidates that simultaneously have:

```text
low C-eigen J-lock residual,
nontrivial Q/P quadrature support,
phase-robust signed-birth flip,
reasonable directed-imbalance / transverse-complementarity scores,
without using delta_beta/H2/kappa in the selection.
```

The candidate-level pair-harmonic projection cannot be known without applying the candidate. Therefore this audit uses Q/P balance and commutator-area support as candidate-level proxies, while `delta_beta2` remains audit-only.

## Comparative summary

{table}

## Top candidates

{top_table}

## Interpretation

A robust positive result would require nonzero `all_pareto_gate_pass_count`, especially for `real_growth`, with strict-sym remaining empty in the upstream package. A negative result means that the candidate space contains good C-lock candidates and good phase-flip candidates, but not the same candidates.
"""
    audit = """# SOURCE AUDIT

Derived-only constraints for this package:

- no complex scalars, no `i`;
- no Hodge star, positivity, norm axiom, physical adjoint, or C*-claim;
- no final sym(M) is introduced here; this is a downstream audit of the antisymmetric birth-transport candidate data;
- no arbitrary fitted rotation;
- no delta_beta/H2/kappa used in the Pareto audit score.

This package audits candidate records generated by `test_C_eigen_guided_pairing_rule_gate.py`. If the expected input output directory is missing, rerun that package first or pass `--input` to the directory containing `candidate_eval_rows.csv` subdirectories.
"""
    readme = """# C-eigen guided phase robustness Pareto gate

Run after generating `C_eigen_guided_pairing_rule_out_L2`:

```bash
python3 test_C_eigen_guided_phase_robustness_pareto_gate.py \
  --input C_eigen_guided_pairing_rule_out_L2
```

Outputs:

- `matched_phase_candidate_rows.csv`
- `pareto_front_rows.csv`
- `top_pareto_candidate_rows.csv`
- `comparative_pareto_summary.csv`
- `RESULTS.md`
- `SUMMARY.md`
- `SOURCE_AUDIT.md`
"""
    return smd, rmd, audit, readme


def package(out: Path, zip_path: Path, script_path: Path, include_paths: List[Path]) -> None:
    if zip_path.exists():
        zip_path.unlink()
    with zipfile.ZipFile(zip_path, 'w', zipfile.ZIP_DEFLATED) as z:
        z.write(script_path, script_path.name)
        for p in include_paths:
            if p.exists() and p.is_file():
                z.write(p, p.name)
        for p in sorted(out.rglob('*')):
            if p.is_file():
                z.write(p, p.resolve().relative_to(Path.cwd()))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument('--input', default='C_eigen_guided_pairing_rule_out_L2')
    ap.add_argument('--out', default='C_eigen_guided_phase_robustness_pareto_out_L2')
    ap.add_argument('--zip', default='cnna_C_eigen_guided_phase_robustness_pareto_gate_pkg_L2.zip')
    ap.add_argument('--lock-residual-threshold', type=float, default=0.20)
    ap.add_argument('--lock-max-threshold', type=float, default=0.30)
    ap.add_argument('--flip-abs-threshold', type=float, default=0.25)
    ap.add_argument('--signed-amp-min', type=float, default=0.10)
    ap.add_argument('--qp-balance-min', type=float, default=0.50)
    ap.add_argument('--comm-abs-min', type=float, default=1e-6)
    args = ap.parse_args()
    input_dir = Path(args.input)
    out = Path(args.out)
    if out.exists():
        shutil.rmtree(out)
    out.mkdir(parents=True)
    rows = load_candidate_rows(input_dir)
    matched, unmatched = match_phase_rows(rows, args)
    front = pareto_front(matched)
    summaries, top_rows = summarize(matched, front)
    write_csv(out / 'matched_phase_candidate_rows.csv', matched)
    write_csv(out / 'pareto_front_rows.csv', front)
    write_csv(out / 'top_pareto_candidate_rows.csv', top_rows)
    write_csv(out / 'comparative_pareto_summary.csv', summaries)
    write_csv(out / 'unmatched_candidate_rows.csv', unmatched)
    summary = {
        'args': vars(args),
        'input_dir': str(input_dir),
        'candidate_rows_loaded': len(rows),
        'matched_phase_candidate_rows': len(matched),
        'unmatched_rows': len(unmatched),
        'summary_rows': summaries,
        'top_pareto_rows': top_rows[:50],
    }
    (out / 'comparative_summary.json').write_text(json.dumps(summary, indent=2), encoding='utf-8')
    smd, rmd, audit, readme = make_docs(summaries, top_rows, args)
    (out / 'SUMMARY.md').write_text(smd, encoding='utf-8')
    (out / 'RESULTS.md').write_text(rmd, encoding='utf-8')
    (out / 'SOURCE_AUDIT.md').write_text(audit, encoding='utf-8')
    (out / 'README.md').write_text(readme, encoding='utf-8')
    include_paths = [
        Path('test_C_eigen_guided_pairing_rule_gate.py'),
        Path('test_C_eigen_quadrature_refinement_gate.py'),
        Path('test_pair_J_alignment_search_gate.py'),
        Path('test_pairing_quadrature_adjoint_pairing_gate.py'),
        Path('test_signed_quadrature_area_kappa_gate.py'),
        Path('test_pairing_transport_antisym_birth_coherence_gate.py'),
        Path('test_nonlinear_asymmetry_cascade_growth.py'),
        Path('cnna_non_shelling_core.py'),
    ]
    package(out, Path(args.zip), Path(__file__), include_paths)
    print(json.dumps({
        'zip': args.zip,
        'out': args.out,
        'candidate_rows_loaded': len(rows),
        'matched_phase_candidate_rows': len(matched),
        'summary': summaries,
        'top_candidates': top_rows[:8],
    }, indent=2))


if __name__ == '__main__':
    main()
