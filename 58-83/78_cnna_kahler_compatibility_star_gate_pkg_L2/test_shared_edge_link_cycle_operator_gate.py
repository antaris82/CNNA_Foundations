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
import test_pair_J_alignment_search_gate as p61
import test_dual_assembly_order_context_ablation_gate as p69
import test_signed_Jlock_role_coupling_gate as p70
import test_assembly_motif_basis_diagonalization_gate as p71
import test_edge_interface_motif_operator_gate as p72

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
    return p69.fbool(x)


def norm(x: np.ndarray) -> float:
    return float(np.linalg.norm(np.asarray(x, dtype=float)))


def sign_nonzero(x: float) -> float:
    if x > 0:
        return 1.0
    if x < 0:
        return -1.0
    return 0.0


def edge_vertices(face: Face, edge: Edge) -> Optional[int]:
    xs = [v for v in face if v not in set(edge)]
    return xs[0] if len(xs) == 1 else None


def all_faces_containing_edge(K, edge: Edge) -> List[Face]:
    e = set(edge)
    return sorted([tuple(sorted(f)) for f in K.faces() if e.issubset(set(f))])


def all_tets_containing_edge(K, edge: Edge) -> List[Tuple[int, int, int, int]]:
    e = set(edge)
    return sorted([tuple(sorted(t)) for t in K.tets if e.issubset(set(t))])


def link_graph_data(K, edge: Edge) -> dict:
    faces = all_faces_containing_edge(K, edge)
    tets = all_tets_containing_edge(K, edge)
    verts = sorted({edge_vertices(f, edge) for f in faces if edge_vertices(f, edge) is not None})
    link_edges = []
    for t in tets:
        rest = sorted([v for v in t if v not in set(edge)])
        if len(rest) == 2:
            link_edges.append(tuple(rest))
    link_edges = sorted(set(link_edges))
    deg = {v: 0 for v in verts}
    for a, b in link_edges:
        deg[a] = deg.get(a, 0) + 1
        deg[b] = deg.get(b, 0) + 1
    closed_cycle = bool(len(verts) >= 3 and len(link_edges) >= len(verts) and all(deg.get(v, 0) >= 2 for v in verts))
    return {'faces': faces, 'tets': tets, 'link_vertices': verts, 'link_edges': link_edges, 'degrees': deg, 'closed_cycle': closed_cycle}


def face_for_side_vertex(edge: Edge, x: int) -> Face:
    return tuple(sorted((edge[0], edge[1], x)))


def edge_unit(model, edge: Edge) -> np.ndarray:
    a, b = edge
    return core.unit(model.nodes[int(b)].pos - model.nodes[int(a)].pos)


def link_order(model, K, edge: Edge, order_mode: str) -> List[Face]:
    data = link_graph_data(K, edge)
    verts = list(data['link_vertices'])
    if len(verts) < 2:
        return [face_for_side_vertex(edge, x) for x in verts]
    if order_mode == 'birth_order':
        verts.sort(key=lambda v: (model.nodes[int(v)].birth_time, model.nodes[int(v)].birth_order, model.address_tuple(int(v)), int(v)))
    elif order_mode == 'address_order':
        verts.sort(key=lambda v: (model.address_tuple(int(v)), model.nodes[int(v)].birth_time, int(v)))
    elif order_mode == 'geometric_angle':
        a, b = edge
        pa = model.nodes[int(a)].pos
        u = edge_unit(model, edge)
        _, e1, e2 = core.frame_from_radial(u)
        vals = []
        for v in verts:
            r = model.nodes[int(v)].pos - pa
            r = r - np.dot(r, u) * u
            vals.append((math.atan2(float(np.dot(r, e2)), float(np.dot(r, e1))), v))
        vals.sort()
        verts = [v for _, v in vals]
    else:
        raise ValueError(order_mode)
    return [face_for_side_vertex(edge, x) for x in verts]


def set_kappa_mirror(model) -> dict:
    old = {i: n.birth_order for i, n in model.nodes.items()}
    for n in model.nodes.values():
        if n.birth_order in (1, 2, 3):
            n.birth_order = 4 - n.birth_order
    return old


def restore_birth_orders(model, old: dict) -> None:
    for i, bo in old.items():
        model.nodes[i].birth_order = bo


def face_axial_on_edge(model, face: Face, edge: Edge, args: argparse.Namespace) -> float:
    try:
        v = p56.axial(p56.face_K_directed(model, face, args.source, args.phase_sign, args.antisym_eta, args.erase_phase_for_strict_sym))
        return float(np.dot(v, edge_unit(model, edge)))
    except Exception:
        return 0.0


def link_circulation_stats(model, K, edge: Edge, order_mode: str, args: argparse.Namespace) -> dict:
    faces = link_order(model, K, edge, order_mode)
    if len(faces) < 2:
        return {
            'link_face_count': len(faces), 'link_circulation_birth': 0.0,
            'link_circulation_K_edge': 0.0, 'link_circulation_abs': 0.0,
        }
    cb = 0.0
    ck = 0.0
    cabs = 0.0
    for i, f in enumerate(faces):
        inc = p72.edge_incidence_sign(f, edge)
        bs = p72.face_birth_signature(model, f)
        kv = face_axial_on_edge(model, f, edge, args)
        # Pure cyclic side sum: no target J, no Hodge, no harmonic projection.
        cb += inc * bs
        ck += inc * kv
        cabs += abs(kv) + abs(bs)
    return {
        'link_face_count': len(faces),
        'link_circulation_birth': float(cb),
        'link_circulation_K_edge': float(ck),
        'link_circulation_abs': float(cabs),
        'link_birth_over_abs': float(cb / (cabs + EPS)),
        'link_K_edge_over_abs': float(ck / (cabs + EPS)),
    }


def kappa_flip_stats(model, K, edge: Edge, order_mode: str, args: argparse.Namespace) -> dict:
    ident = link_circulation_stats(model, K, edge, order_mode, args)
    old = set_kappa_mirror(model)
    try:
        kap = link_circulation_stats(model, K, edge, order_mode, args)
    finally:
        restore_birth_orders(model, old)
    s0 = float(ident['link_circulation_K_edge'])
    s1 = float(kap['link_circulation_K_edge'])
    b0 = float(ident['link_circulation_birth'])
    b1 = float(kap['link_circulation_birth'])
    return {
        **{f'ident_{k}': v for k, v in ident.items()},
        **{f'kappa_{k}': v for k, v in kap.items()},
        'K_edge_kappa_flip_abs': abs(s0 + s1) / (abs(s0) + abs(s1) + EPS),
        'birth_kappa_flip_abs': abs(b0 + b1) / (abs(b0) + abs(b1) + EPS),
        'K_edge_sign_identity': sign_nonzero(s0),
        'K_edge_sign_kappa': sign_nonzero(s1),
        'birth_sign_identity': sign_nonzero(b0),
        'birth_sign_kappa': sign_nonzero(b1),
        'K_edge_has_nonzero_circulation': bool(abs(s0) > args.link_circulation_threshold),
        'K_edge_kappa_flips': bool(abs(s0) > args.link_circulation_threshold and abs(s1) > args.link_circulation_threshold and (abs(s0 + s1) / (abs(s0) + abs(s1) + EPS)) < args.kappa_flip_threshold),
        'birth_kappa_flips': bool(abs(b0) > args.link_circulation_threshold and abs(b1) > args.link_circulation_threshold and (abs(b0 + b1) / (abs(b0) + abs(b1) + EPS)) < args.kappa_flip_threshold),
    }


def sl_for(idx: Dict[Face, int], f: Face) -> slice:
    i = idx[f]
    return slice(3*i, 3*i+3)


def build_link_cycle_operator(model, K, carrier_faces: List[Face], edge: Edge, order_mode: str, block_mode: str) -> Tuple[np.ndarray, np.ndarray, dict]:
    n = 3 * len(carrier_faces)
    J = np.zeros((n, n), dtype=float)
    C = np.zeros((n, n), dtype=float)
    idx = {f: i for i, f in enumerate(carrier_faces)}
    ordered = link_order(model, K, edge, order_mode)
    if len(ordered) < 2:
        return J, C, {'cycle_edge_count': 0, 'cycle_operator_norm': 0.0, 'cycle_order': ''}
    closed = len(ordered) >= 3
    pairs = list(zip(ordered, ordered[1:]))
    if closed:
        pairs.append((ordered[-1], ordered[0]))
    count = 0
    for f, g in pairs:
        if f not in idx or g not in idx:
            continue
        s = p72.edge_incidence_sign(f, edge) * p72.edge_incidence_sign(g, edge) * p72.provenance_pair_sign(model, f, g)
        if block_mode == 'identity':
            M = np.eye(3)
        elif block_mode == 'edge_projector':
            u = edge_unit(model, edge)
            M = np.outer(u, u)
        elif block_mode == 'edge_complement':
            u = edge_unit(model, edge)
            M = np.eye(3) - np.outer(u, u)
        else:
            raise ValueError(block_mode)
        sf, sg = sl_for(idx, f), sl_for(idx, g)
        J[sf, sg] += s * M
        J[sg, sf] += -s * M.T
        C[sf, sg] += s * M
        C[sg, sf] += s * M.T
        count += 1
    return J, C, {
        'cycle_edge_count': count,
        'cycle_operator_norm': norm(J),
        'cycle_order': ' -> '.join(str(list(f)) for f in ordered),
        'cycle_closed_by_order': bool(closed),
        **{f'link_graph_{k}': (str(v) if isinstance(v, dict) else v) for k, v in link_graph_data(K, edge).items() if k in ('link_vertices','link_edges','closed_cycle')},
    }


def scaled_add(Jbase: np.ndarray, Cbase: np.ndarray, Jcyc: np.ndarray, Ccyc: np.ndarray, scale: str) -> Tuple[np.ndarray, np.ndarray, float]:
    nb, nc = norm(Jbase), norm(Jcyc)
    if nc < EPS:
        lam = 0.0
    elif scale == 'unit':
        lam = 1.0
    elif scale == 'half_base_norm':
        lam = 0.5 * nb / (nc + EPS)
    elif scale == 'base_norm':
        lam = nb / (nc + EPS)
    elif scale == 'quarter_base_norm':
        lam = 0.25 * nb / (nc + EPS)
    else:
        raise ValueError(scale)
    return Jbase + lam * Jcyc, Cbase + lam * Ccyc, lam


def carrier_faces_for_mode(K, union_faces: List[Face], edge: Edge, carrier_mode: str) -> List[Face]:
    if carrier_mode == 'motif_only':
        return sorted(set(union_faces))
    if carrier_mode == 'with_link_faces':
        return sorted(set(union_faces) | set(all_faces_containing_edge(K, edge)))
    raise ValueError(carrier_mode)


def motif_vecs(faces: List[Face], pA: dict, pB: dict) -> Dict[str, np.ndarray]:
    return {
        'A_Q': p71.put_pair_vec_union(faces, pA, 'Q'),
        'A_P': p71.put_pair_vec_union(faces, pA, 'P'),
        'B_Q': p71.put_pair_vec_union(faces, pB, 'Q'),
        'B_P': p71.put_pair_vec_union(faces, pB, 'P'),
    }


def shared_edges_between_pairs(pA: dict, pB: dict) -> List[Edge]:
    out = []
    for f in [pA['fa'], pA['fb']]:
        for g in [pB['fa'], pB['fb']]:
            e = p72.shared_edge(f, g)
            if e is not None:
                out.append(e)
    return sorted(set(out))


def link_cycle_metrics_for_motif(model, K, pA: dict, pB: dict, args: argparse.Namespace) -> List[dict]:
    union = p72.union_faces_from_pairs(pA, pB)
    rows: List[dict] = []
    for edge in shared_edges_between_pairs(pA, pB):
        for carrier_mode in args.carrier_modes:
            faces = carrier_faces_for_mode(K, union, edge, carrier_mode)
            Jbase, Cbase = p71.union_JC(faces, [pA, pB])
            vecs = motif_vecs(faces, pA, pB)
            base = p71.motif_metrics('base_union', Jbase, Cbase, vecs)
            for order_mode in args.link_order_modes:
                flip = kappa_flip_stats(model, K, edge, order_mode, args)
                for block_mode in args.link_block_modes:
                    Jcyc, Ccyc, cstat = build_link_cycle_operator(model, K, faces, edge, order_mode, block_mode)
                    only = p71.motif_metrics('link_cycle_only', Jcyc, Ccyc, vecs)
                    for scale in args.link_scales:
                        Jeff, Ceff, lam = scaled_add(Jbase, Cbase, Jcyc, Ccyc, scale)
                        eff = p71.motif_metrics('link_cycle_eff', Jeff, Ceff, vecs)
                        decision_pass = bool(
                            flip['K_edge_has_nonzero_circulation']
                            and flip['K_edge_kappa_flips']
                            and eff['link_cycle_eff_J_QP_subspace_mean_resid'] < args.decision_lock_threshold
                            and eff['link_cycle_eff_projected_J2_plus_I_resid'] < args.decision_J2_threshold
                            and eff['link_cycle_eff_J_span_leakage'] < args.decision_leak_threshold
                        )
                        rows.append({
                            'shared_edge': str(list(edge)),
                            'carrier_mode': carrier_mode,
                            'carrier_face_count': len(faces),
                            'order_mode': order_mode,
                            'block_mode': block_mode,
                            'scale_mode': scale,
                            'lambda_used': lam,
                            'improves_lock_over_base': bool(eff['link_cycle_eff_J_QP_subspace_mean_resid'] < base['base_union_J_QP_subspace_mean_resid']),
                            'improves_J2_over_base': bool(eff['link_cycle_eff_projected_J2_plus_I_resid'] < base['base_union_projected_J2_plus_I_resid']),
                            'decision_gate_pass': decision_pass,
                            'local_negative_abort_marker': bool(not decision_pass and eff['link_cycle_eff_J_QP_subspace_mean_resid'] >= args.local_negative_residual_floor),
                            **base,
                            **only,
                            **eff,
                            **cstat,
                            **flip,
                        })
    return rows


def link_cycle_rows(model, K, pairing_log: List[dict], assembly_log: List[dict], args: argparse.Namespace) -> Tuple[List[dict], dict]:
    rows: List[dict] = []
    for i, a in enumerate(assembly_log):
        if not fbool(a.get('assembly_applied')):
            continue
        A1, A2 = p71.parse_faces_from_assembly(a, 'A')
        B1, B2 = p71.parse_faces_from_assembly(a, 'B')
        if A1 is None or A2 is None or B1 is None or B2 is None:
            continue
        pA = p71.pair_fields(model, A1, A2, args)
        pB = p71.pair_fields(model, B1, B2, args)
        if pA is None or pB is None:
            continue
        signed = p71.signed_motif_stats(pA, pB)
        singA = p61.J_lock_residual(pA['J'], np.concatenate([pA['Q_a'], pA['Q_b']]), np.concatenate([pA['P_a'], pA['P_b']]))
        singB = p61.J_lock_residual(pB['J'], np.concatenate([pB['Q_a'], pB['Q_b']]), np.concatenate([pB['P_a'], pB['P_b']]))
        for r in link_cycle_metrics_for_motif(model, K, pA, pB, args):
            r.update({
                'assembly_index': i,
                'context': a.get('context',''),
                'A_face_a': str(list(A1)), 'A_face_b': str(list(A2)),
                'B_face_a': str(list(B1)), 'B_face_b': str(list(B2)),
                'pair_local_mean_J_lock_raw_QP': 0.5*(singA['J_lock_mean_resid'] + singB['J_lock_mean_resid']),
                'pair_local_max_J_lock_raw_QP': max(singA['J_lock_max_resid'], singB['J_lock_max_resid']),
                **signed,
            })
            rows.append(r)
    if not rows:
        return rows, {
            'assembly_count': 0,
            'link_cycle_row_count': 0,
            'decision_gate_pass_count': 0,
            'local_negative_abort_marker': True,
        }
    def vals(k): return [float(r[k]) for r in rows if k in r and np.isfinite(float(r[k]))]
    def mn(k):
        xs = vals(k); return float(np.min(xs)) if xs else 0.0
    def avg(k):
        xs = vals(k); return float(np.mean(xs)) if xs else 0.0
    best = min(rows, key=lambda r: float(r['link_cycle_eff_J_QP_subspace_mean_resid']))
    active_rows = [r for r in rows if int(r.get('cycle_edge_count', 0)) > 0]
    best_active = min(active_rows, key=lambda r: float(r['link_cycle_eff_J_QP_subspace_mean_resid'])) if active_rows else None
    best_decision = min(rows, key=lambda r: (not fbool(r['decision_gate_pass']), float(r['link_cycle_eff_J_QP_subspace_mean_resid'])))
    summary = {
        'assembly_count': len({r['assembly_index'] for r in rows}),
        'link_cycle_row_count': len(rows),
        'base_union_best_mean_resid': mn('base_union_J_QP_subspace_mean_resid'),
        'base_union_best_J2_resid': mn('base_union_projected_J2_plus_I_resid'),
        'link_cycle_eff_best_mean_resid': mn('link_cycle_eff_J_QP_subspace_mean_resid'),
        'link_cycle_eff_best_active_mean_resid': float(best_active['link_cycle_eff_J_QP_subspace_mean_resid']) if best_active else 0.0,
        'link_cycle_eff_best_active_J2_resid': float(best_active['link_cycle_eff_projected_J2_plus_I_resid']) if best_active else 0.0,
        'link_cycle_eff_best_active_cycle_edge_count': int(best_active['cycle_edge_count']) if best_active else 0,
        'link_cycle_eff_avg_mean_resid': avg('link_cycle_eff_J_QP_subspace_mean_resid'),
        'link_cycle_eff_best_max_resid': mn('link_cycle_eff_J_QP_subspace_max_resid'),
        'link_cycle_eff_best_span_leakage': mn('link_cycle_eff_J_span_leakage'),
        'link_cycle_eff_best_J2_resid': mn('link_cycle_eff_projected_J2_plus_I_resid'),
        'link_cycle_only_best_mean_resid': mn('link_cycle_only_J_QP_subspace_mean_resid'),
        'link_cycle_only_best_J2_resid': mn('link_cycle_only_projected_J2_plus_I_resid'),
        'improves_lock_over_base_count': sum(1 for r in rows if fbool(r['improves_lock_over_base'])),
        'improves_J2_over_base_count': sum(1 for r in rows if fbool(r['improves_J2_over_base'])),
        'nonzero_K_circulation_count': sum(1 for r in rows if fbool(r['K_edge_has_nonzero_circulation'])),
        'kappa_flip_count': sum(1 for r in rows if fbool(r['K_edge_kappa_flips'])),
        'birth_kappa_flip_count': sum(1 for r in rows if fbool(r['birth_kappa_flips'])),
        'decision_gate_pass_count': sum(1 for r in rows if fbool(r['decision_gate_pass'])),
        'best_carrier_mode': best.get('carrier_mode',''),
        'best_order_mode': best.get('order_mode',''),
        'best_block_mode': best.get('block_mode',''),
        'best_scale_mode': best.get('scale_mode',''),
        'best_active_carrier_mode': best_active.get('carrier_mode','') if best_active else '',
        'best_active_order_mode': best_active.get('order_mode','') if best_active else '',
        'best_active_block_mode': best_active.get('block_mode','') if best_active else '',
        'best_active_scale_mode': best_active.get('scale_mode','') if best_active else '',
        'best_shared_edge': best.get('shared_edge',''),
        'best_cycle_order': best.get('cycle_order',''),
        'best_identity_K_circ': float(best.get('ident_link_circulation_K_edge',0.0)),
        'best_kappa_K_circ': float(best.get('kappa_link_circulation_K_edge',0.0)),
        'best_K_flip_abs': float(best.get('K_edge_kappa_flip_abs',0.0)),
        'best_decision_like_row_lock': float(best_decision.get('link_cycle_eff_J_QP_subspace_mean_resid',0.0)),
        'local_negative_abort_marker': bool(sum(1 for r in rows if fbool(r['decision_gate_pass'])) == 0 and mn('link_cycle_eff_J_QP_subspace_mean_resid') >= 0.35),
    }
    return rows, summary


def option_tag(args: argparse.Namespace) -> str:
    return p69.option_tag(args)


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
    motif_rows, motif_summary = p71.assembly_motif_rows(model, K, pairing_log, assembly_log, args)
    edge_rows, edge_summary = p72.edge_interface_rows(model, K, pairing_log, assembly_log, args)
    link_rows, link_summary = link_cycle_rows(model, K, pairing_log, assembly_log, args)
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
    write_csv(vout / 'shared_edge_link_cycle_operator_rows.csv', link_rows)
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
        'edge_interface_metrics': edge_summary,
        'shared_edge_link_cycle_metrics': link_summary,
        'automatic_pairings_applied': sum(1 for x in pairing_log if fbool(x.get('applied'))),
        'assemblies_applied': sum(1 for x in assembly_log if fbool(x.get('assembly_applied'))),
        'assemblies_attempted': len(assembly_log),
        'decision_used_delta_beta_any': sel['decision_used_delta_beta_any'],
    }
    (vout / 'variant_shared_edge_link_cycle_summary.json').write_text(json.dumps(summary, indent=2), encoding='utf-8')
    return summary


def clone_args(args: argparse.Namespace, **updates) -> argparse.Namespace:
    d = vars(args).copy(); d.update(updates); return argparse.Namespace(**d)


def run_experiment(args: argparse.Namespace, out: Path) -> List[dict]:
    p70.patch_modules()
    rows = []
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
    for opt in options:
        for variant in opt.variants:
            rows.append(run_variant(variant, opt, out))
    return rows


def slim(r: dict) -> dict:
    a = r['auto_metrics']; dm = r['directed_metrics']; sm = r['signed_quadrature_metrics']; am = r['alignment_metrics']; mm = r['motif_basis_metrics']; em = r['edge_interface_metrics']; lm = r['shared_edge_link_cycle_metrics']
    return {
        'option': r['option'], 'variant': r['variant'],
        'beta': [a['beta0'], a['beta1'], a['beta2'], a['beta3']],
        'pairings': r['automatic_pairings_applied'], 'assemblies': r['assemblies_applied'],
        'pair_harm': dm['pair_transport_harmonic_ratio'],
        'Q_harm': am['Q_even_harmonic_ratio'], 'P_harm': am['P_odd_harmonic_ratio'],
        'pair_local_J_lock': am['best_per_pair_mean_J_lock_resid'],
        'signed_birth': sm['signed_birth_over_abs_sum_ratio'],
        'union_motif_lock': mm.get('union_sum_best_mean_resid',0.0),
        'edge_interface_lock': em.get('edge_interface_best_mean_resid',0.0),
        'link_cycle_lock': lm.get('link_cycle_eff_best_mean_resid',0.0),
        'link_cycle_active_lock': lm.get('link_cycle_eff_best_active_mean_resid',0.0),
        'link_cycle_J2': lm.get('link_cycle_eff_best_J2_resid',0.0),
        'link_cycle_active_J2': lm.get('link_cycle_eff_best_active_J2_resid',0.0),
        'nonzero_circ': lm.get('nonzero_K_circulation_count',0),
        'kappa_flip_count': lm.get('kappa_flip_count',0),
        'decision_pass': lm.get('decision_gate_pass_count',0),
        'local_negative_abort_marker': lm.get('local_negative_abort_marker',False),
        'best_order_mode': lm.get('best_order_mode',''),
        'best_block_mode': lm.get('best_block_mode',''),
        'best_carrier_mode': lm.get('best_carrier_mode',''),
        'used_delta_beta': r['decision_used_delta_beta_any'],
    }


def write_comparative(out: Path, rows: List[dict]) -> None:
    flat = []
    for r in rows:
        s = slim(r); lm = r['shared_edge_link_cycle_metrics']
        flat.append({
            **s,
            'beta0': s['beta'][0], 'beta1': s['beta'][1], 'beta2': s['beta'][2], 'beta3': s['beta'][3],
            'link_cycle_row_count': lm.get('link_cycle_row_count',0),
            'base_union_best_mean_resid': lm.get('base_union_best_mean_resid',0.0),
            'base_union_best_J2_resid': lm.get('base_union_best_J2_resid',0.0),
            'link_cycle_eff_avg_mean_resid': lm.get('link_cycle_eff_avg_mean_resid',0.0),
            'link_cycle_eff_best_active_mean_resid': lm.get('link_cycle_eff_best_active_mean_resid',0.0),
            'link_cycle_eff_best_active_J2_resid': lm.get('link_cycle_eff_best_active_J2_resid',0.0),
            'link_cycle_eff_best_active_cycle_edge_count': lm.get('link_cycle_eff_best_active_cycle_edge_count',0),
            'link_cycle_eff_best_max_resid': lm.get('link_cycle_eff_best_max_resid',0.0),
            'link_cycle_eff_best_span_leakage': lm.get('link_cycle_eff_best_span_leakage',0.0),
            'link_cycle_only_best_mean_resid': lm.get('link_cycle_only_best_mean_resid',0.0),
            'improves_lock_over_base_count': lm.get('improves_lock_over_base_count',0),
            'improves_J2_over_base_count': lm.get('improves_J2_over_base_count',0),
            'birth_kappa_flip_count': lm.get('birth_kappa_flip_count',0),
            'best_shared_edge': lm.get('best_shared_edge',''),
            'best_cycle_order': lm.get('best_cycle_order',''),
            'best_identity_K_circ': lm.get('best_identity_K_circ',0.0),
            'best_kappa_K_circ': lm.get('best_kappa_K_circ',0.0),
            'best_K_flip_abs': lm.get('best_K_flip_abs',0.0),
        })
    write_csv(out / 'comparative_shared_edge_link_cycle_operator_summary.csv', flat)


def make_docs(summary: dict) -> Tuple[str, str, str, str]:
    rows = summary['variant_rows']
    lines = ['| option | variant | beta | pairs | asm | pair harm | Q harm | P harm | pair J-lock | signed | union lock | edge-if lock | link-cycle best | active link lock | link J2 | active J2 | nonzero circ | kappa flips | decision pass | abort marker | used dBeta? |',
             '|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|']
    for r in rows:
        s = slim(r)
        lines.append(f"| {s['option']} | {s['variant']} | ({s['beta'][0]},{s['beta'][1]},{s['beta'][2]},{s['beta'][3]}) | {s['pairings']} | {s['assemblies']} | {s['pair_harm']:.6g} | {s['Q_harm']:.6g} | {s['P_harm']:.6g} | {s['pair_local_J_lock']:.6g} | {s['signed_birth']:.6g} | {s['union_motif_lock']:.6g} | {s['edge_interface_lock']:.6g} | {s['link_cycle_lock']:.6g} | {s['link_cycle_active_lock']:.6g} | {s['link_cycle_J2']:.6g} | {s['link_cycle_active_J2']:.6g} | {s['nonzero_circ']} | {s['kappa_flip_count']} | {s['decision_pass']} | {s['local_negative_abort_marker']} | {s['used_delta_beta']} |")
    table = '\n'.join(lines)
    nonstrict = [r for r in rows if r['variant'] != 'strict_symmetrized_control']
    best = min(nonstrict, key=lambda r: r['shared_edge_link_cycle_metrics'].get('link_cycle_eff_best_mean_resid',9.0), default=None)
    best_payload = slim(best) if best else {}
    all_nonstrict_abort = all(r['shared_edge_link_cycle_metrics'].get('local_negative_abort_marker', False) for r in nonstrict) if nonstrict else False
    smd = f"""# SUMMARY — shared edge link-cycle operator decision gate

Model label:
CNNA growing primal simplicial complex with deterministic sequential provenance growth, directed antisymmetric birth-transport operators, signed-Jlock two-pair assembly motifs, and a shared-edge link-cycle diagnostic.

This package is intentionally a **decision test**, not another open-ended move to a larger local environment.  It tests whether the shared-edge link around an A/B assembly edge provides the missing local orientation/alignment operator.

Positive decision gate requires all of the following:

```text
1. nonzero signed circulation on the shared-edge link,
2. circulation flips under sibling birth-order kappa mirror,
3. strict_sym remains null,
4. used_delta_beta remains false,
5. effective motif Q/P J-lock residual < {summary['args']['decision_lock_threshold']},
6. projected J^2 + I residual < {summary['args']['decision_J2_threshold']}.
```

If the link-cycle still gives residuals around 0.4--0.6, the documented interpretation is not "try the next bigger local structure", but: the current local line supports Q/P carrier structure without deriving a good local J-orientation on the actual Q/P motif space.

{table}

## Best non-strict row

```json
{json.dumps(best_payload, indent=2)}
```

## Local negative decision marker

```json
{{"all_non_strict_variants_abort_marked": {str(all_nonstrict_abort).lower()}}}
```
"""
    rmd = f"""# RESULTS — shared edge link-cycle operator decision gate

## Comparative table

{table}

## Interpretation

`link-cycle lock` is a residual.  Smaller is better.  Values around 0.4--0.6 mean the tested operator does not act as a good derived complex structure on the actual Q/P motif space.

This gate is stricter than the previous edge-interface package.  It requires signed circulation and kappa flip in addition to a low Q/P-J residual and projected J^2 behavior.  Merely finding nonzero magnitude or beta2 is not counted as success.

## Stop/continue rule

- If `decision pass > 0`, the shared-edge link is a serious candidate for the missing local alignment operator.
- If `decision pass = 0` and the best non-strict residual remains >= 0.35, this package marks the local link-cycle path as negative for J derivation on the actual Q/P motif space.  The next scientific step should be interpretation/formal obstruction documentation, not automatic escalation to a larger fan.
"""
    audit = """# SOURCE AUDIT

Carried forward:

- Q/P channels and beta2 carrier survive in the real-growth path.
- Local pair C/J algebra exists, but pair-local J_pair(Q)=P does not lock.
- Two-pair assemblies carry beta2/QP/signed magnitude better than single-pair candidates.
- Motif-basis and simple edge-interface operators did not reduce the residual below the decision range.

This package adds only the shared-edge link-cycle diagnostic.  The cycle order is derived either from birth/address order or from the generated geometric embedding angle around the shared edge.  The operator uses face-to-face cyclic handoff blocks with signs from boundary incidence and birth/provenance signatures.  No i, global J, Hodge star, physical adjoint, positivity, C*-norm, final sym(M), or delta-beta/H2 decision is introduced.

Caveat: this is still a Python diagnostic on the L2 model, not a Lean theorem.  A positive result would require formalization; a negative result localizes the obstruction only for this tested derived operator class and model regime.
"""
    readme = """# Shared edge link-cycle operator decision gate

Run:

```bash
python3 test_shared_edge_link_cycle_operator_gate.py
```

The script evaluates dynamic A/B assemblies and tests whether the shared-edge link cycle supplies a derived local orientation/alignment operator for the Q/P motif space.
"""
    return smd, rmd, audit, readme


def package(out: Path, zip_path: Path) -> None:
    files = [
        Path(__file__).name,
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
    ap.add_argument('--carrier-modes', nargs='*', default=['motif_only','with_link_faces'])
    ap.add_argument('--link-order-modes', nargs='*', default=['birth_order','address_order','geometric_angle'])
    ap.add_argument('--link-block-modes', nargs='*', default=['identity','edge_projector','edge_complement'])
    ap.add_argument('--link-scales', nargs='*', default=['unit','quarter_base_norm','half_base_norm','base_norm'])
    ap.add_argument('--link-circulation-threshold', type=float, default=1e-6)
    ap.add_argument('--kappa-flip-threshold', type=float, default=0.20)
    ap.add_argument('--decision-lock-threshold', type=float, default=0.20)
    ap.add_argument('--decision-J2-threshold', type=float, default=0.25)
    ap.add_argument('--decision-leak-threshold', type=float, default=0.35)
    ap.add_argument('--local-negative-residual-floor', type=float, default=0.35)
    ap.add_argument('--out', default='shared_edge_link_cycle_operator_out_L2')
    ap.add_argument('--zip', default='cnna_shared_edge_link_cycle_operator_gate_pkg_L2.zip')
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
