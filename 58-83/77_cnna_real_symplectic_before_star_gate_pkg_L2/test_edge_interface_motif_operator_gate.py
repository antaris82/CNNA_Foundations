#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import math
import shutil
import zipfile
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np

import cnna_non_shelling_core as core
import test_nonlinear_asymmetry_cascade_growth as nl
import test_harmonic_k_orientation_kappa_gate as hk
import test_pairing_transport_antisym_birth_coherence_gate as p56
import test_pairing_quadrature_adjoint_pairing_gate as p60
import test_pair_J_alignment_search_gate as p61
import test_dual_pairing_assembly_growth_rule_gate as p68
import test_dual_assembly_order_context_ablation_gate as p69
import test_signed_Jlock_role_coupling_gate as p70
import test_assembly_motif_basis_diagonalization_gate as p71

EPS = 1e-12
Face = Tuple[int, int, int]
Edge = Tuple[int, int]


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


def fbool(x) -> bool:
    return p68.fbool(x)


def norm(x: np.ndarray) -> float:
    return float(np.linalg.norm(np.asarray(x, dtype=float)))


def sign_nonzero(x: float) -> float:
    if x > 0:
        return 1.0
    if x < 0:
        return -1.0
    return 0.0


def face_edges(f: Face) -> List[Edge]:
    a, b, c = tuple(sorted(f))
    return [tuple(sorted((a, b))), tuple(sorted((a, c))), tuple(sorted((b, c)))]


def shared_edge(f: Face, g: Face) -> Optional[Edge]:
    s = sorted(set(f) & set(g))
    if len(s) == 2:
        return (s[0], s[1])
    return None


def edge_incidence_sign(face: Face, edge: Edge) -> float:
    """Boundary sign of oriented edge in oriented sorted face [v0,v1,v2].

    For [v0,v1,v2], d = [v1,v2] - [v0,v2] + [v0,v1].  Since the
    stored edge orientation is sorted, this gives a pure incidence sign.
    """
    v = tuple(sorted(face))
    e = tuple(sorted(edge))
    if e == (v[1], v[2]):
        return +1.0
    if e == (v[0], v[2]):
        return -1.0
    if e == (v[0], v[1]):
        return +1.0
    return 0.0


def face_birth_signature(model, face: Face) -> float:
    """Derived provenance scalar from sibling order labels only.

    birth_order 1,2,3 is centered to -1,0,+1.  No topology or harmonic data is
    used here.  The value is used only as a tie-breaking/signature factor for
    the edge-interface sign.
    """
    vals = []
    for v in face:
        bo = getattr(model.nodes[int(v)], 'birth_order', 0)
        vals.append(float(bo - 2))
    return float(sum(vals))


def provenance_pair_sign(model, f: Face, g: Face) -> float:
    sf = face_birth_signature(model, f)
    sg = face_birth_signature(model, g)
    raw = sf - sg
    if abs(raw) < EPS:
        # deterministic fallback from birth-time orientation, still provenance-only
        raw = sum(model.nodes[int(v)].birth_time for v in f) - sum(model.nodes[int(v)].birth_time for v in g)
    s = sign_nonzero(raw)
    return s if s != 0.0 else 1.0


def edge_unit(model, edge: Edge) -> np.ndarray:
    a, b = edge
    pa = model.nodes[int(a)].pos
    pb = model.nodes[int(b)].pos
    return core.unit(pb - pa)


def parse_face_from_field(x) -> Optional[Face]:
    return p71.parse_face_from_field(x)


def parse_faces_from_assembly(row: dict, prefix: str) -> Tuple[Optional[Face], Optional[Face]]:
    return p71.parse_faces_from_assembly(row, prefix)


def pair_fields(model, fa: Face, fb: Face, args: argparse.Namespace) -> Optional[dict]:
    return p71.pair_fields(model, fa, fb, args)


def interface_links_for_assembly(model, pA: dict, pB: dict) -> List[dict]:
    Afaces = [pA['fa'], pA['fb']]
    Bfaces = [pB['fa'], pB['fb']]
    links: List[dict] = []
    seen = set()
    for f in Afaces:
        for g in Bfaces:
            e = shared_edge(f, g)
            if e is None:
                continue
            key = (f, g, e)
            if key in seen:
                continue
            seen.add(key)
            si = edge_incidence_sign(f, e) * edge_incidence_sign(g, e)
            sp = provenance_pair_sign(model, f, g)
            s = si * sp
            links.append({
                'from_face': f,
                'to_face': g,
                'edge': e,
                'incidence_sign': si,
                'provenance_sign': sp,
                'interface_sign': s,
            })
    return links


def union_faces_from_pairs(pA: dict, pB: dict) -> List[Face]:
    return sorted({pA['fa'], pA['fb'], pB['fa'], pB['fb']})


def sl_for(idx: Dict[Face, int], f: Face) -> slice:
    i = idx[f]
    return slice(3*i, 3*i + 3)


def build_interface_operator(model, faces: List[Face], links: List[dict], mode: str) -> Tuple[np.ndarray, np.ndarray, dict]:
    n = 3 * len(faces)
    Jint = np.zeros((n, n), dtype=float)
    Cint = np.zeros((n, n), dtype=float)
    idx = {f: i for i, f in enumerate(faces)}
    if not links:
        return Jint, Cint, {'link_count': 0, 'interface_norm': 0.0}
    weight = 1.0 / math.sqrt(len(links))
    for L in links:
        f, g, e = L['from_face'], L['to_face'], L['edge']
        if f not in idx or g not in idx:
            continue
        s = float(L['interface_sign']) * weight
        if mode == 'incidence_identity':
            M = np.eye(3)
        elif mode == 'edge_projector':
            u = edge_unit(model, e)
            M = np.outer(u, u)
        elif mode == 'edge_complement_projector':
            u = edge_unit(model, e)
            M = np.eye(3) - np.outer(u, u)
        else:
            raise ValueError(mode)
        sf, sg = sl_for(idx, f), sl_for(idx, g)
        # Jint is a skew interface handoff: Face -> shared edge -> Face and back with opposite sign.
        Jint[sf, sg] += s * M
        Jint[sg, sf] += -s * M.T
        # Cint is the corresponding even exchange diagnostic.
        Cint[sf, sg] += s * M
        Cint[sg, sf] += s * M.T
    return Jint, Cint, {'link_count': len(links), 'interface_norm': norm(Jint)}


def scaled_add_base_interface(Jbase: np.ndarray, Cbase: np.ndarray, Jint: np.ndarray, Cint: np.ndarray, scale: str) -> Tuple[np.ndarray, np.ndarray, float]:
    nb, ni = norm(Jbase), norm(Jint)
    if ni < EPS:
        lam = 0.0
    elif scale == 'unit':
        lam = 1.0
    elif scale == 'base_norm':
        lam = nb / (ni + EPS)
    elif scale == 'half_base_norm':
        lam = 0.5 * nb / (ni + EPS)
    else:
        raise ValueError(scale)
    return Jbase + lam * Jint, Cbase + lam * Cint, lam


def interface_metrics_for_motif(model, pA: dict, pB: dict, args: argparse.Namespace) -> List[dict]:
    faces = union_faces_from_pairs(pA, pB)
    Jbase, Cbase = p71.union_JC(faces, [pA, pB])
    vecs = {
        'A_Q': p71.put_pair_vec_union(faces, pA, 'Q'),
        'A_P': p71.put_pair_vec_union(faces, pA, 'P'),
        'B_Q': p71.put_pair_vec_union(faces, pB, 'Q'),
        'B_P': p71.put_pair_vec_union(faces, pB, 'P'),
    }
    links = interface_links_for_assembly(model, pA, pB)
    base = p71.motif_metrics('base_union', Jbase, Cbase, vecs)
    out: List[dict] = []
    for mode in args.interface_modes:
        Jint, Cint, istat = build_interface_operator(model, faces, links, mode)
        int_only = p71.motif_metrics('interface_only', Jint, Cint, vecs)
        for scale in args.interface_scales:
            Jeff, Ceff, lam = scaled_add_base_interface(Jbase, Cbase, Jint, Cint, scale)
            eff = p71.motif_metrics('edge_interface', Jeff, Ceff, vecs)
            row = {
                'interface_mode': mode,
                'interface_scale': scale,
                'lambda_used': lam,
                'union_face_count': len(faces),
                'shared_edge_link_count': istat['link_count'],
                'interface_operator_norm': istat['interface_norm'],
                'edge_interface_improves_mean': bool(eff['edge_interface_J_QP_subspace_mean_resid'] < base['base_union_J_QP_subspace_mean_resid']),
                'edge_interface_improves_J2': bool(eff['edge_interface_projected_J2_plus_I_resid'] < base['base_union_projected_J2_plus_I_resid']),
                **{k:v for k,v in base.items()},
                **{k:v for k,v in int_only.items()},
                **{k:v for k,v in eff.items()},
            }
            out.append(row)
    return out


def edge_interface_rows(model, K, pairing_log: List[dict], assembly_log: List[dict], args: argparse.Namespace) -> Tuple[List[dict], dict]:
    rows: List[dict] = []
    for i, a in enumerate(assembly_log):
        if not fbool(a.get('assembly_applied')):
            continue
        A1, A2 = parse_faces_from_assembly(a, 'A')
        B1, B2 = parse_faces_from_assembly(a, 'B')
        if A1 is None or A2 is None or B1 is None or B2 is None:
            continue
        pA = pair_fields(model, A1, A2, args)
        pB = pair_fields(model, B1, B2, args)
        if pA is None or pB is None:
            continue
        links = interface_links_for_assembly(model, pA, pB)
        signed = p71.signed_motif_stats(pA, pB)
        singA = p61.J_lock_residual(pA['J'], np.concatenate([pA['Q_a'], pA['Q_b']]), np.concatenate([pA['P_a'], pA['P_b']]))
        singB = p61.J_lock_residual(pB['J'], np.concatenate([pB['Q_a'], pB['Q_b']]), np.concatenate([pB['P_a'], pB['P_b']]))
        for r in interface_metrics_for_motif(model, pA, pB, args):
            r.update({
                'assembly_index': i,
                'context': a.get('context',''),
                'A_face_a': str(list(A1)), 'A_face_b': str(list(A2)),
                'B_face_a': str(list(B1)), 'B_face_b': str(list(B2)),
                'pair_local_mean_J_lock_raw_QP': 0.5*(singA['J_lock_mean_resid'] + singB['J_lock_mean_resid']),
                'pair_local_max_J_lock_raw_QP': max(singA['J_lock_max_resid'], singB['J_lock_max_resid']),
                'interface_links': ';'.join([f"{list(L['from_face'])}->{list(L['to_face'])}@{list(L['edge'])}:s={L['interface_sign']:+.0f}" for L in links]),
                **signed,
            })
            rows.append(r)
    if not rows:
        return rows, {
            'assembly_count': 0,
            'interface_row_count': 0,
            'edge_interface_best_mean_resid': 0.0,
            'edge_interface_gate_pass_count': 0,
        }
    def vals(k): return [float(r[k]) for r in rows if k in r and np.isfinite(float(r[k]))]
    def avg(k):
        xs = vals(k); return float(np.mean(xs)) if xs else 0.0
    def mn(k):
        xs = vals(k); return float(np.min(xs)) if xs else 0.0
    best_row = min(rows, key=lambda r: float(r['edge_interface_J_QP_subspace_mean_resid']))
    summary = {
        'assembly_count': len({r['assembly_index'] for r in rows}),
        'interface_row_count': len(rows),
        'shared_edge_operator_rows': sum(1 for r in rows if int(r['shared_edge_link_count']) > 0),
        'base_union_best_mean_resid': mn('base_union_J_QP_subspace_mean_resid'),
        'base_union_avg_mean_resid': avg('base_union_J_QP_subspace_mean_resid'),
        'base_union_best_J2_resid': mn('base_union_projected_J2_plus_I_resid'),
        'interface_only_best_mean_resid': mn('interface_only_J_QP_subspace_mean_resid'),
        'interface_only_avg_mean_resid': avg('interface_only_J_QP_subspace_mean_resid'),
        'interface_only_best_J2_resid': mn('interface_only_projected_J2_plus_I_resid'),
        'edge_interface_best_mean_resid': mn('edge_interface_J_QP_subspace_mean_resid'),
        'edge_interface_avg_mean_resid': avg('edge_interface_J_QP_subspace_mean_resid'),
        'edge_interface_best_max_resid': mn('edge_interface_J_QP_subspace_max_resid'),
        'edge_interface_best_span_leakage': mn('edge_interface_J_span_leakage'),
        'edge_interface_best_J2_resid': mn('edge_interface_projected_J2_plus_I_resid'),
        'edge_interface_gate_pass_count': sum(1 for r in rows if fbool(r['edge_interface_gate_pass'])),
        'edge_interface_improves_mean_count': sum(1 for r in rows if fbool(r['edge_interface_improves_mean'])),
        'edge_interface_improves_J2_count': sum(1 for r in rows if fbool(r['edge_interface_improves_J2'])),
        'best_interface_mode': best_row.get('interface_mode',''),
        'best_interface_scale': best_row.get('interface_scale',''),
        'best_lambda_used': float(best_row.get('lambda_used',0.0)),
        'best_interface_links': best_row.get('interface_links',''),
        'signed_birth_weighted_avg': avg('motif_signed_birth_over_abs_weighted'),
    }
    return rows, summary


def run_variant_option_with_edge_interface(variant: str, args: argparse.Namespace, out: Path) -> dict:
    model = nl.build_model(variant, args)
    model.grow(args.max_level)
    baseline_K = core.build_dynamic_outward_ngf_complex(model)
    baseline_metrics = core.full_metrics(model, baseline_K, args.source)
    tag = p69.option_tag(args)
    vout = out / tag / variant
    vout.mkdir(parents=True, exist_ok=True)
    K, birth_log, pairing_log, assembly_log, candidate_rows = p69.build_ablation_complex(model, args, variant, vout)
    auto_metrics = core.full_metrics(model, K, args.source)
    dm, pair_rows, top_rows, three_rows = p56.directed_metrics(model, K, pairing_log, args)
    sm, signed_rows, signed_face_rows = p69.p59.signed_quadrature_rows(model, K, pairing_log, args)
    am, align_pair_rows, align_candidate_rows, align_candidate_summary = p61.alignment_search_metrics(model, K, pairing_log, args)
    sel = p69.summarize_selection(pairing_log)
    motif_rows, motif_summary = p71.assembly_motif_rows(model, K, pairing_log, assembly_log, args)
    edge_rows, edge_summary = edge_interface_rows(model, K, pairing_log, assembly_log, args)
    write_csv(vout / 'birth_geometry_log.csv', birth_log)
    write_csv(vout / 'assembly_pairing_log.csv', pairing_log)
    write_csv(vout / 'assembly_ablation_log.csv', assembly_log)
    write_csv(vout / 'candidate_eval_rows.csv', candidate_rows)
    write_csv(vout / 'directed_pair_rows.csv', pair_rows)
    write_csv(vout / 'signed_quadrature_rows.csv', signed_rows)
    write_csv(vout / 'alignment_pair_rows.csv', align_pair_rows)
    write_csv(vout / 'alignment_candidate_rows.csv', align_candidate_rows)
    write_csv(vout / 'alignment_candidate_summary.csv', align_candidate_summary)
    write_csv(vout / 'assembly_motif_basis_rows.csv', motif_rows)
    write_csv(vout / 'edge_interface_motif_operator_rows.csv', edge_rows)
    summary = {
        'variant': variant,
        'option': tag,
        'assembly_order': args.assembly_order,
        'require_connected_assembly': args.require_connected_assembly,
        'require_strong_assembly_context': args.require_strong_assembly_context,
        'allow_B_reuse_A_faces': args.allow_B_reuse_A_faces,
        'baseline_metrics': baseline_metrics,
        'auto_metrics': auto_metrics,
        'directed_metrics': dm,
        'signed_quadrature_metrics': sm,
        'alignment_metrics': am,
        'selection_metrics': sel,
        'motif_basis_metrics': motif_summary,
        'edge_interface_metrics': edge_summary,
        'automatic_pairings_applied': sum(1 for x in pairing_log if fbool(x.get('applied'))),
        'assemblies_applied': sum(1 for x in assembly_log if fbool(x.get('assembly_applied'))),
        'assemblies_attempted': len(assembly_log),
        'decision_used_delta_beta_any': sel['decision_used_delta_beta_any'],
    }
    (vout / 'variant_option_edge_interface_summary.json').write_text(json.dumps(summary, indent=2), encoding='utf-8')
    return summary


def clone_args(args: argparse.Namespace, **updates) -> argparse.Namespace:
    d = vars(args).copy(); d.update(updates); return argparse.Namespace(**d)


def run_experiment(args: argparse.Namespace, out: Path) -> List[dict]:
    p70.patch_modules()
    options = []
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
                options.append(clone_args(args, **updates))
    rows = []
    for opt in options:
        for variant in args.variants:
            rows.append(run_variant_option_with_edge_interface(variant, opt, out))
    return rows


def slim(r: dict) -> dict:
    a = r['auto_metrics']; dm = r['directed_metrics']; sm = r['signed_quadrature_metrics']; am = r['alignment_metrics']; mm = r['motif_basis_metrics']; em = r['edge_interface_metrics']
    return {
        'option': r['option'], 'variant': r['variant'],
        'beta': [a['beta0'], a['beta1'], a['beta2'], a['beta3']],
        'pairings': r['automatic_pairings_applied'],
        'assemblies': r['assemblies_applied'],
        'pair_harm': dm['pair_transport_harmonic_ratio'],
        'Q_harm': am['Q_even_harmonic_ratio'],
        'P_harm': am['P_odd_harmonic_ratio'],
        'pair_local_J_lock': am['best_per_pair_mean_J_lock_resid'],
        'signed_birth': sm['signed_birth_over_abs_sum_ratio'],
        'motif_count': mm.get('assembly_count',0),
        'union_motif_lock': mm.get('union_sum_best_mean_resid',0.0),
        'base_union_lock': em.get('base_union_best_mean_resid',0.0),
        'edge_interface_lock': em.get('edge_interface_best_mean_resid',0.0),
        'edge_interface_J2': em.get('edge_interface_best_J2_resid',0.0),
        'edge_interface_pass': em.get('edge_interface_gate_pass_count',0),
        'edge_interface_improve_count': em.get('edge_interface_improves_mean_count',0),
        'best_mode': em.get('best_interface_mode',''),
        'best_scale': em.get('best_interface_scale',''),
        'used_delta_beta': r['decision_used_delta_beta_any'],
    }


def write_comparative(out: Path, rows: List[dict]) -> None:
    flat = []
    for r in rows:
        s = slim(r); em = r['edge_interface_metrics']
        flat.append({
            **s,
            'beta0': s['beta'][0], 'beta1': s['beta'][1], 'beta2': s['beta'][2], 'beta3': s['beta'][3],
            'interface_row_count': em.get('interface_row_count',0),
            'shared_edge_operator_rows': em.get('shared_edge_operator_rows',0),
            'base_union_avg_mean_resid': em.get('base_union_avg_mean_resid',0.0),
            'interface_only_best_mean_resid': em.get('interface_only_best_mean_resid',0.0),
            'interface_only_best_J2_resid': em.get('interface_only_best_J2_resid',0.0),
            'edge_interface_avg_mean_resid': em.get('edge_interface_avg_mean_resid',0.0),
            'edge_interface_best_max_resid': em.get('edge_interface_best_max_resid',0.0),
            'edge_interface_best_span_leakage': em.get('edge_interface_best_span_leakage',0.0),
            'edge_interface_best_J2_resid': em.get('edge_interface_best_J2_resid',0.0),
            'edge_interface_improves_J2_count': em.get('edge_interface_improves_J2_count',0),
            'best_lambda_used': em.get('best_lambda_used',0.0),
            'best_interface_links': em.get('best_interface_links',''),
        })
    write_csv(out / 'comparative_edge_interface_motif_operator_summary.csv', flat)


def make_docs(summary: dict) -> Tuple[str, str, str, str]:
    rows = summary['variant_rows']
    lines = ['| option | variant | beta | pairs | asm | pair harm | Q harm | P harm | pair J-lock | signed | base union | edge-if lock | edge-if J2 | pass | best mode | scale | used dBeta? |',
             '|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|---|---:|']
    for r in rows:
        s = slim(r)
        lines.append(f"| {s['option']} | {s['variant']} | ({s['beta'][0]},{s['beta'][1]},{s['beta'][2]},{s['beta'][3]}) | {s['pairings']} | {s['assemblies']} | {s['pair_harm']:.6g} | {s['Q_harm']:.6g} | {s['P_harm']:.6g} | {s['pair_local_J_lock']:.6g} | {s['signed_birth']:.6g} | {s['base_union_lock']:.6g} | {s['edge_interface_lock']:.6g} | {s['edge_interface_J2']:.6g} | {s['edge_interface_pass']} | {s['best_mode']} | {s['best_scale']} | {s['used_delta_beta']} |")
    table = '\n'.join(lines)
    nonstrict = [r for r in rows if r['variant'] != 'strict_symmetrized_control']
    best = min(nonstrict, key=lambda r: r['edge_interface_metrics'].get('edge_interface_best_mean_resid',9.0), default=None)
    best_payload = slim(best) if best else {}
    smd = f"""# SUMMARY — edge interface motif operator gate

Model label:
CNNA growing primal simplicial complex with deterministic sequential provenance growth, directed antisymmetric birth-transport operators, signed-Jlock two-pair assembly motifs, and a derived shared-edge interface diagnostic.

This test checks whether the previous motif-basis obstruction is caused by missing A/B interface coupling.  The new operator is built on the actual union of faces of each complete A/B assembly:

```text
Face -> shared edge -> Face
```

Primary interface mode:

```text
incidence_identity:
  skew handoff from face to face through the common edge,
  sign = boundary-incidence sign * birth/provenance signature sign.
```

Secondary diagnostic modes:

```text
edge_projector
edge_complement_projector
```

No i, global J, Hodge star, *, positivity, C*-norm, final sym(M), or delta-beta/H2 decision is introduced.  Delta-beta/H2/harmonic quantities are measured after the fact only.

{table}

## Best non-strict row

```json
{json.dumps(best_payload, indent=2)}
```
"""
    rmd = f"""# RESULTS — edge interface motif operator gate

## Comparative table

{table}

## Gate criterion

A positive edge-interface result requires the effective operator

```text
J_eff = J_pair_union + lambda * J_edge_interface
```

to improve the motif Q/P subspace lock and give acceptable projected J^2 behavior, while strict_sym stays null and used_delta_beta remains false.

The tested interface signs come only from oriented boundary incidence of faces on their shared edge plus birth/provenance signatures.  The primary `incidence_identity` mode avoids Hodge and uses no metric edge rotation.  Edge-direction projector modes are included only as secondary geometric diagnostics.

## Interpretation rules

- If `edge_interface_lock` drops strongly below `base_union_lock` and projected J2 improves, the obstruction was likely a missing derived edge-interface handoff.
- If Q/P and beta2 remain positive but `edge_interface_lock` stays high, then merely adding incidence-level edge coupling is not enough.
- If strict_sym is zero, the path remains tied to nonsymmetric provenance growth.
"""
    audit = """# SOURCE AUDIT

Carried forward:

- Single-pair tests found Q/P channels and local C/J pair algebra but no dynamic J_pair(Q)=P lock.
- Kappa and tradeoff tests split beta2, C-lock, Q/P, and signed kappa flip across different single-pair candidates.
- Dual-pair assembly tests showed those roles can coexist in a two-pair motif.
- Motif-basis diagonalization improved over raw pair-local lock but did not close the gate.

This package changes only the motif operator by adding a derived shared-edge interface.  The primary interface uses only simplicial boundary incidence and birth/provenance signatures.  No i, global J, Hodge, *, positivity, C*-norm, final sym(M), or delta-beta decision is introduced.

Caveat: this is still a Python diagnostic, not a formal theorem.  A positive result would still require Lean formalization and a derived-only proof that the interface operator is forced by the CNNA generator/provenance chain.
"""
    readme = """# Edge interface motif operator gate

Run:

```bash
python3 test_edge_interface_motif_operator_gate.py
```

The script evaluates complete dynamic two-pair A/B assemblies and tests whether a derived shared-edge interface operator improves the motif-level Q/P lock.
"""
    return smd, rmd, audit, readme


def package(out: Path, zip_path: Path) -> None:
    files = [
        Path(__file__).name,
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
    ap.add_argument('--out', default='edge_interface_motif_operator_out_L2')
    ap.add_argument('--zip', default='cnna_edge_interface_motif_operator_gate_pkg_L2.zip')
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
