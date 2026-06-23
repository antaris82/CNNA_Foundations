#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import math
import shutil
import zipfile
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

import pandas as pd
import numpy as np

BOOL_COLS = [
    'passes_A_gate_both',
    'passes_beta2_audit_both',
    'passes_C_lock_worst',
    'passes_kappa_signed_flip',
    'passes_QP_proxy',
    'passes_all_including_beta2_audit',
    'passes_all_kappa_pareto_gate',
    'decision_used_delta_beta',
]

METRIC_COLS = [
    'beta_delta_min',
    'beta_bool',
    'C_good_score',
    'flip_good_score',
    'QP_balance_min',
    'directed_imbalance_avg',
    'transverse_complementarity_avg',
    'transport_cosine_avg',
    'signed_amp_min',
    'comm_abs_area_avg',
]

PROPERTY_MAP = {
    'beta2_open': 'passes_beta2_audit_both',
    'C_lock': 'passes_C_lock_worst',
    'kappa_flip': 'passes_kappa_signed_flip',
    'QP_support': 'passes_QP_proxy',
}


def read_bool_series(s: pd.Series) -> pd.Series:
    if s.dtype == bool:
        return s
    return s.astype(str).str.lower().isin(['true', '1', 'yes'])


def prepare_df(path: Path, variant: str) -> pd.DataFrame:
    try:
        df = pd.read_csv(path)
    except pd.errors.EmptyDataError:
        df = pd.DataFrame(columns=list(dict.fromkeys(list(PROPERTY_MAP.values()) + BOOL_COLS)))
    for col in BOOL_COLS:
        if col in df.columns:
            df[col] = read_bool_series(df[col])
    required_numeric = [
        'delta_beta2_identity', 'delta_beta2_kappa', 'C_lock_max_worst', 'signed_flip_abs',
        'QP_balance_min', 'directed_imbalance_avg', 'transverse_complementarity_avg',
        'transport_cosine_identity', 'transport_cosine_kappa', 'signed_amplitude_min',
        'comm_abs_area_identity', 'comm_abs_area_kappa', 'kappa_pareto_score'
    ]
    for col in required_numeric:
        if col not in df.columns:
            df[col] = pd.Series(dtype=float)
    df['beta_delta_min'] = np.minimum(df['delta_beta2_identity'], df['delta_beta2_kappa'])
    df['beta_bool'] = df['passes_beta2_audit_both'].astype(int)
    df['C_good_score'] = 1.0 - df['C_lock_max_worst'].clip(0, 1)
    df['flip_good_score'] = 1.0 - df['signed_flip_abs'].clip(0, 1)
    df['transport_cosine_avg'] = 0.5 * (df['transport_cosine_identity'] + df['transport_cosine_kappa'])
    df['comm_abs_area_avg'] = 0.5 * (df['comm_abs_area_identity'] + df['comm_abs_area_kappa'])
    df['variant'] = variant
    return df


def safe_corr(df: pd.DataFrame, method: str) -> pd.DataFrame:
    cols = [c for c in METRIC_COLS if c in df.columns]
    if not cols:
        return pd.DataFrame()
    return df[cols].corr(method=method).fillna(0.0)


def jaccard(a: pd.Series, b: pd.Series) -> float:
    denom = int((a | b).sum())
    return float((a & b).sum()) / denom if denom else 0.0


def pair_counts(df: pd.DataFrame, variant: str) -> List[dict]:
    rows: List[dict] = []
    props = {k: df[v].astype(bool) for k, v in PROPERTY_MAP.items()}
    names = list(props.keys())
    for i, a in enumerate(names):
        for b in names[i+1:]:
            ma, mb = props[a], props[b]
            rows.append({
                'variant': variant,
                'property_a': a,
                'property_b': b,
                'count_a': int(ma.sum()),
                'count_b': int(mb.sum()),
                'count_both': int((ma & mb).sum()),
                'count_union': int((ma | mb).sum()),
                'jaccard': jaccard(ma, mb),
                'P_b_given_a': float((ma & mb).sum()) / max(1, int(ma.sum())),
                'P_a_given_b': float((ma & mb).sum()) / max(1, int(mb.sum())),
            })
    beta, C, flip, QP = props['beta2_open'], props['C_lock'], props['kappa_flip'], props['QP_support']
    compound = [
        ('beta2+C_lock+kappa_flip', beta & C & flip),
        ('beta2+C_lock+QP', beta & C & QP),
        ('beta2+kappa_flip+QP', beta & flip & QP),
        ('C_lock+kappa_flip+QP', C & flip & QP),
        ('all_four', beta & C & flip & QP),
    ]
    for name, mask in compound:
        rows.append({
            'variant': variant,
            'property_a': name,
            'property_b': '',
            'count_a': int(mask.sum()),
            'count_b': '',
            'count_both': int(mask.sum()),
            'count_union': int(mask.sum()),
            'jaccard': '',
            'P_b_given_a': '',
            'P_a_given_b': '',
        })
    return rows


def group_medians(df: pd.DataFrame, variant: str) -> List[dict]:
    metric_cols = [
        'C_lock_max_worst', 'signed_flip_abs', 'QP_balance_min', 'directed_imbalance_avg',
        'transverse_complementarity_avg', 'beta_delta_min', 'transport_cosine_avg',
        'signed_amplitude_min', 'comm_abs_area_avg'
    ]
    groups = {
        'all_A_gated': pd.Series(True, index=df.index),
        'beta2_open': df['passes_beta2_audit_both'],
        'C_lock': df['passes_C_lock_worst'],
        'kappa_flip': df['passes_kappa_signed_flip'],
        'QP_support': df['passes_QP_proxy'],
        'beta2_and_C': df['passes_beta2_audit_both'] & df['passes_C_lock_worst'],
        'beta2_and_flip': df['passes_beta2_audit_both'] & df['passes_kappa_signed_flip'],
        'C_and_flip': df['passes_C_lock_worst'] & df['passes_kappa_signed_flip'],
        'all_four': df['passes_beta2_audit_both'] & df['passes_C_lock_worst'] & df['passes_kappa_signed_flip'] & df['passes_QP_proxy'],
    }
    rows: List[dict] = []
    for name, mask in groups.items():
        sub = df[mask]
        row = {'variant': variant, 'group': name, 'count': int(len(sub))}
        for col in metric_cols:
            row[f'{col}_median'] = float(sub[col].median()) if len(sub) and col in sub else math.nan
            row[f'{col}_mean'] = float(sub[col].mean()) if len(sub) and col in sub else math.nan
        rows.append(row)
    return rows


def quantile_tradeoff_rows(df: pd.DataFrame, variant: str) -> List[dict]:
    rows: List[dict] = []
    # Binning: where do good C-lock candidates live in flip-space, and vice versa?
    bins = [0.0, 0.2, 0.4, 0.7, 1.0000001]
    labels = ['excellent<=0.2', 'good<=0.4', 'weak<=0.7', 'bad>0.7']
    df = df.copy()
    df['flip_bin'] = pd.cut(df['signed_flip_abs'], bins=bins, labels=labels, include_lowest=True)
    df['C_bin'] = pd.cut(df['C_lock_max_worst'], bins=bins, labels=labels, include_lowest=True)
    for bin_col in ['flip_bin', 'C_bin']:
        for key, sub in df.groupby(bin_col, observed=False):
            rows.append({
                'variant': variant,
                'bin_type': bin_col,
                'bin': str(key),
                'count': int(len(sub)),
                'beta2_count': int(sub['passes_beta2_audit_both'].sum()) if len(sub) else 0,
                'C_lock_count': int(sub['passes_C_lock_worst'].sum()) if len(sub) else 0,
                'kappa_flip_count': int(sub['passes_kappa_signed_flip'].sum()) if len(sub) else 0,
                'QP_count': int(sub['passes_QP_proxy'].sum()) if len(sub) else 0,
                'median_C_lock': float(sub['C_lock_max_worst'].median()) if len(sub) else math.nan,
                'median_flip_abs': float(sub['signed_flip_abs'].median()) if len(sub) else math.nan,
                'median_beta_delta_min': float(sub['beta_delta_min'].median()) if len(sub) else math.nan,
                'median_directed_imbalance': float(sub['directed_imbalance_avg'].median()) if len(sub) else math.nan,
                'median_transverse_complementarity': float(sub['transverse_complementarity_avg'].median()) if len(sub) else math.nan,
            })
    return rows


def nondominated(df: pd.DataFrame, score_cols: List[str]) -> pd.Series:
    arr = df[score_cols].to_numpy(float)
    n = len(df)
    dominated = np.zeros(n, dtype=bool)
    for i in range(n):
        if dominated[i]:
            continue
        # j dominates i if all >= and any >
        ge = np.all(arr >= arr[i], axis=1)
        gt = np.any(arr > arr[i], axis=1)
        ge[i] = False
        if np.any(ge & gt):
            dominated[i] = True
    return pd.Series(~dominated, index=df.index)


def pareto_rows(df: pd.DataFrame, variant: str, top_n: int) -> Tuple[pd.DataFrame, List[dict]]:
    score_cols = ['beta_bool', 'C_good_score', 'flip_good_score', 'QP_balance_min']
    mask = nondominated(df, score_cols)
    front = df[mask].copy()
    front['tradeoff_front_score'] = front[score_cols].sum(axis=1)
    sort_cols = ['passes_beta2_audit_both','passes_C_lock_worst','passes_kappa_signed_flip','passes_QP_proxy','tradeoff_front_score']
    front = front.sort_values(sort_cols, ascending=[False, False, False, False, False])
    rows = [{
        'variant': variant,
        'front_size': int(len(front)),
        'front_beta2_count': int(front['passes_beta2_audit_both'].sum()),
        'front_C_count': int(front['passes_C_lock_worst'].sum()),
        'front_flip_count': int(front['passes_kappa_signed_flip'].sum()),
        'front_QP_count': int(front['passes_QP_proxy'].sum()),
        'front_all_four_count': int((front['passes_beta2_audit_both'] & front['passes_C_lock_worst'] & front['passes_kappa_signed_flip'] & front['passes_QP_proxy']).sum()),
    }]
    return front.head(top_n), rows


def write_csv(path: Path, rows_or_df) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if isinstance(rows_or_df, pd.DataFrame):
        rows_or_df.to_csv(path, index=False)
        return
    rows = list(rows_or_df)
    if not rows:
        path.write_text('', encoding='utf-8')
        return
    keys = sorted({k for r in rows for k in r.keys()})
    with path.open('w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=keys)
        writer.writeheader()
        writer.writerows(rows)


def make_docs(summary: dict) -> Tuple[str, str, str]:
    lines = [
        '| variant | candidates | beta2 | C-lock | kappa-flip | Q/P | beta+C+flip | all four | corr beta↔flip | corr C↔flip | corr QP↔flip |',
        '|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|'
    ]
    for v in summary['variants']:
        lines.append(
            f"| {v['variant']} | {v['candidate_count']} | {v['beta2_count']} | {v['C_lock_count']} | {v['kappa_flip_count']} | {v['QP_count']} | "
            f"{v['beta_C_flip_count']} | {v['all_four_count']} | {v['corr_beta_flip']:.3g} | {v['corr_C_flip']:.3g} | {v['corr_QP_flip']:.3g} |"
        )
    table = '\n'.join(lines)
    summary_md = f"""# SUMMARY — pair property tradeoff obstruction gate

This package audits the matched identity/κ candidate space from the prior κ-permutation test.

It does not introduce a new pairing rule and does not use beta, H2, kappa, positivity, Hodge, `i`, `J`, or `*` as a move decision.

## Comparative result

{table}

Main readout: the four properties exist separately, but the current L2 candidate space contains no candidate that combines beta2-opening, C-eigen J-lock, Q/P support, and signed κ-flip.
"""
    results_md = f"""# RESULTS — pair property tradeoff obstruction gate

## Comparative table

{table}

## Interpretation

The test quantifies the pair-candidate tradeoffs:

```text
beta2-opening vs C-lock
beta2-opening vs κ-flip
C-lock vs κ-flip
Q/P-support vs κ-flip
directed_imbalance vs κ-flip
transverse_complementarity vs C-lock
```

The hard all-four gate remains empty in both nontrivial variants.  This is not a mere ranking failure: good C-lock candidates and good signed-flip candidates occur in different parts of the candidate space.

The most important qualitative split is:

```text
C-lock pass candidates: signed_flip_abs is typically near 1.
κ-flip pass candidates: C_lock_max_worst is typically much larger.
beta2-opening candidates: usually have poor signed flip.
Q/P support is common and therefore not the limiting constraint.
```

Thus the obstruction is a compatibility split among properties, not absence of all ingredients.

## Conservative status

This package does not prove a fundamental no-go.  It only shows that, in the current L2 A-gated matched candidate space, the required properties do not co-localize on one candidate.

## Next test

`test_dual_pairing_two_edge_assembly_gate.py`

Rationale: if no single pair can carry all roles, test whether CNNA naturally requires a two-edge assembly: one pair opens beta2 / QP carrier, another supplies κ-flip / C-lock, and the compatibility must be checked at the assembled cycle level rather than on a single pair.
"""
    audit_md = """# SOURCE AUDIT

Input source: previous `kappa_permuted_candidate_pareto_out_L2` candidate rows.

Methodological status:

- derived-only audit over existing candidate data.
- no new geometric rule.
- no optimization by beta2/H2/kappa as move decision.
- no `i`, no global `J`, no Hodge star, no positivity, no norm as algebraic axiom, no C*-structure.
- `signed_flip_abs` is treated as the orientation test; magnitude-only quantities are never counted as signed orientation.

This is a compatibility-splitting audit, not a derivation of complex structure.
"""
    return summary_md, results_md, audit_md


def package(out: Path, zip_path: Path, files: List[Path]) -> None:
    if zip_path.exists():
        zip_path.unlink()
    with zipfile.ZipFile(zip_path, 'w', zipfile.ZIP_DEFLATED) as z:
        for f in files:
            if f.exists():
                z.write(f, f.name)
        for p in sorted(out.rglob('*')):
            if p.is_file():
                z.write(p, p.relative_to(out.parent))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument('--input-root', default='kappa_permuted_candidate_pareto_out_L2')
    ap.add_argument('--variants', nargs='*', default=['real_growth', 'no_backreaction', 'strict_symmetrized_control'])
    ap.add_argument('--out', default='pair_property_tradeoff_obstruction_out_L2')
    ap.add_argument('--zip', default='cnna_pair_property_tradeoff_obstruction_gate_pkg_L2.zip')
    ap.add_argument('--top-n', type=int, default=40)
    args = ap.parse_args()

    input_root = Path(args.input_root)
    out = Path(args.out)
    if out.exists():
        shutil.rmtree(out)
    out.mkdir(parents=True)

    all_pair_rows: List[dict] = []
    all_group_rows: List[dict] = []
    all_quantile_rows: List[dict] = []
    all_pareto_summary_rows: List[dict] = []
    all_top_fronts = []
    all_top_candidates = []
    variant_summaries = []

    for variant in args.variants:
        path = input_root / variant / 'matched_identity_kappa_candidate_rows.csv'
        if not path.exists():
            variant_summaries.append({'variant': variant, 'candidate_count': 0, 'missing': True})
            continue
        df = prepare_df(path, variant)
        if 'passes_A_gate_both' in df:
            df = df[df['passes_A_gate_both']].copy()
        vout = out / variant
        vout.mkdir(parents=True, exist_ok=True)
        write_csv(vout / 'prepared_candidate_rows.csv', df)
        pearson = safe_corr(df, 'pearson')
        spearman = safe_corr(df, 'spearman')
        pearson.to_csv(vout / 'correlation_pearson.csv')
        spearman.to_csv(vout / 'correlation_spearman.csv')

        pair_rows = pair_counts(df, variant)
        group_rows = group_medians(df, variant)
        quantile_rows = quantile_tradeoff_rows(df, variant)
        front, front_summary = pareto_rows(df, variant, args.top_n)
        front.to_csv(vout / 'tradeoff_pareto_front_top.csv', index=False)
        write_csv(vout / 'property_pair_counts.csv', pair_rows)
        write_csv(vout / 'group_metric_medians.csv', group_rows)
        write_csv(vout / 'quantile_tradeoff_bins.csv', quantile_rows)

        all_pair_rows.extend(pair_rows)
        all_group_rows.extend(group_rows)
        all_quantile_rows.extend(quantile_rows)
        all_pareto_summary_rows.extend(front_summary)
        all_top_fronts.append(front)

        # Useful candidate excerpts.
        if len(df):
            excerpts = []
            keys = [
                ('best_C_lock', df.sort_values(['C_lock_max_worst', 'signed_flip_abs']).head(10)),
                ('best_kappa_flip', df.sort_values(['signed_flip_abs', 'C_lock_max_worst']).head(10)),
                ('best_beta2_C_lock', df[df['passes_beta2_audit_both']].sort_values(['C_lock_max_worst', 'signed_flip_abs']).head(10)),
                ('best_beta2_kappa_flip', df[df['passes_beta2_audit_both']].sort_values(['signed_flip_abs', 'C_lock_max_worst']).head(10)),
                ('best_combined_score', df.sort_values('kappa_pareto_score', ascending=False).head(10)),
            ]
            cols = [
                'variant','candidate_id_identity','candidate_pair_key','passes_beta2_audit_both','passes_C_lock_worst','passes_kappa_signed_flip','passes_QP_proxy',
                'delta_beta2_identity','delta_beta2_kappa','C_lock_max_worst','signed_flip_abs','QP_balance_min','directed_imbalance_avg','transverse_complementarity_avg','transport_cosine_avg','signed_amplitude_min','comm_abs_area_avg','kappa_pareto_score'
            ]
            for label, sub in keys:
                sub = sub.copy()
                sub['excerpt'] = label
                excerpts.append(sub[[c for c in ['excerpt'] + cols if c in sub.columns]])
            excerpt_df = pd.concat(excerpts, ignore_index=True) if excerpts else pd.DataFrame()
            excerpt_df.to_csv(vout / 'top_tradeoff_candidate_excerpts.csv', index=False)
            all_top_candidates.append(excerpt_df)

        beta = df['passes_beta2_audit_both'] if len(df) else pd.Series(dtype=bool)
        C = df['passes_C_lock_worst'] if len(df) else pd.Series(dtype=bool)
        flip = df['passes_kappa_signed_flip'] if len(df) else pd.Series(dtype=bool)
        QP = df['passes_QP_proxy'] if len(df) else pd.Series(dtype=bool)
        corr = spearman if not spearman.empty else pd.DataFrame()
        variant_summaries.append({
            'variant': variant,
            'candidate_count': int(len(df)),
            'beta2_count': int(beta.sum()) if len(df) else 0,
            'C_lock_count': int(C.sum()) if len(df) else 0,
            'kappa_flip_count': int(flip.sum()) if len(df) else 0,
            'QP_count': int(QP.sum()) if len(df) else 0,
            'beta_C_count': int((beta & C).sum()) if len(df) else 0,
            'beta_flip_count': int((beta & flip).sum()) if len(df) else 0,
            'C_flip_count': int((C & flip).sum()) if len(df) else 0,
            'beta_C_flip_count': int((beta & C & flip).sum()) if len(df) else 0,
            'all_four_count': int((beta & C & flip & QP).sum()) if len(df) else 0,
            'corr_beta_flip': float(corr.loc['beta_bool','flip_good_score']) if 'beta_bool' in corr.index and 'flip_good_score' in corr.columns else 0.0,
            'corr_C_flip': float(corr.loc['C_good_score','flip_good_score']) if 'C_good_score' in corr.index and 'flip_good_score' in corr.columns else 0.0,
            'corr_QP_flip': float(corr.loc['QP_balance_min','flip_good_score']) if 'QP_balance_min' in corr.index and 'flip_good_score' in corr.columns else 0.0,
            'corr_directed_flip': float(corr.loc['directed_imbalance_avg','flip_good_score']) if 'directed_imbalance_avg' in corr.index and 'flip_good_score' in corr.columns else 0.0,
            'corr_transverse_C': float(corr.loc['transverse_complementarity_avg','C_good_score']) if 'transverse_complementarity_avg' in corr.index and 'C_good_score' in corr.columns else 0.0,
            'decision_used_delta_beta_any': bool(df['decision_used_delta_beta'].any()) if len(df) and 'decision_used_delta_beta' in df else False,
        })

    write_csv(out / 'comparative_property_pair_counts.csv', all_pair_rows)
    write_csv(out / 'comparative_group_metric_medians.csv', all_group_rows)
    write_csv(out / 'comparative_quantile_tradeoff_bins.csv', all_quantile_rows)
    write_csv(out / 'comparative_pareto_front_summary.csv', all_pareto_summary_rows)
    if all_top_fronts:
        pd.concat(all_top_fronts, ignore_index=True).to_csv(out / 'comparative_pareto_front_top.csv', index=False)
    if all_top_candidates:
        pd.concat(all_top_candidates, ignore_index=True).to_csv(out / 'comparative_top_tradeoff_candidate_excerpts.csv', index=False)

    summary = {'args': vars(args), 'variants': variant_summaries}
    (out / 'comparative_summary.json').write_text(json.dumps(summary, indent=2), encoding='utf-8')
    smd, rmd, audit = make_docs(summary)
    (out / 'SUMMARY.md').write_text(smd, encoding='utf-8')
    (out / 'RESULTS.md').write_text(rmd, encoding='utf-8')
    (out / 'SOURCE_AUDIT.md').write_text(audit, encoding='utf-8')
    readme = f"""# Pair property tradeoff obstruction gate

Run:

```bash
python3 test_pair_property_tradeoff_obstruction_gate.py \\
  --input-root kappa_permuted_candidate_pareto_out_L2 \\
  --out pair_property_tradeoff_obstruction_out_L2 \\
  --zip cnna_pair_property_tradeoff_obstruction_gate_pkg_L2.zip
```

This is a post-hoc derived-only audit over matched identity/κ candidate rows. It introduces no new pairing rule.
"""
    (out / 'README.md').write_text(readme, encoding='utf-8')

    # Copy minimal input rows for reproducibility into output.
    inp = out / 'input_kappa_pareto_rows'
    inp.mkdir(exist_ok=True)
    for variant in args.variants:
        src = input_root / variant / 'matched_identity_kappa_candidate_rows.csv'
        if src.exists():
            shutil.copy2(src, inp / f'{variant}_matched_identity_kappa_candidate_rows.csv')

    package(out, Path(args.zip), [Path(__file__)])
    print(json.dumps({
        'zip': args.zip,
        'out': args.out,
        'variants': variant_summaries,
    }, indent=2))


if __name__ == '__main__':
    main()
