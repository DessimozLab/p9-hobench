"""Per-family hierarchical-orthology metrics (DESIGN §6, review round 1).

Merge levels are species-tree node indices; ``INF = -1`` = never in one HOG; reference merge
levels are **intervals** ``[lo, hi]`` (see :mod:`hobench.reference`). A predicted merge level
``p`` (lowest level at which both genes share a predicted HOG, ``INF`` if none) is *correct*
iff ``lo ⊆ p ⊆ hi`` along the pair's lineage (``p = INF`` is correct only if ``hi = INF``).

Signed level error (edges) = distance from ``p`` to the nearest interval bound:
``depth(bound) − depth(p)``; **positive = predicted too old = over-splitting** (genes kept
apart too long; ``+inf`` = never merged although the reference merges them), **negative =
too young = over-merging** (e.g. pseudo-orthologs called orthologs; ``-inf`` = merged although
the reference says never).

Public API (all documented below): :func:`reference_matrices`, :func:`representative_matrix`,
:func:`hog_codes`, :func:`predicted_merge_matrix`, :func:`interval_correct`,
:func:`interval_error`, :func:`reference_level_codes`, :func:`evaluate_family`,
:func:`family_integrity`, :func:`panel_context`, :func:`family_scalars`, ``SCALARS``,
``PRIMARY``.

Headline (primary) per-family endpoints (SYNTHESIS B8):
1. ``level_bcubed_f1`` — per-level B-cubed F1 (gene-balanced), mean over levels;
2. ``merge_accuracy_node`` — node-balanced merge-level accuracy;
3. ``event_recovery_0.5`` (and ``_0.8``) — two-sided duplication-event recovery;
4. ``membership_f1`` — family-membership F1 at the family level.
"""
from __future__ import annotations

import math
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from typing import Iterable, Optional, Sequence

import numpy as np

from .predictions import Prediction
from .reference import (ABOVE, INF, NONHOMOLOGOUS_TOKENS, NOT_ASSERTED, PER_PAIR_S0,
                        ReconciledFamily, ReferenceFamily, reconcile)
from .speciestree import SpeciesTree

_BLOCK_FANOUT = 16  # above this many children, write a node's whole block at once
THETAS = (0.5, 0.8)


class ScoringError(ValueError):
    """Raised instead of silently repairing inconsistent inputs during scoring."""


# ---------------------------------------------------------------------- reference side
@dataclass
class RefMatrices:
    lo: np.ndarray     # n×n int32 merge-interval lower bound (NOT_ASSERTED for U pairs)
    hi: np.ndarray     # n×n int32 upper bound (INF allowed)
    node: np.ndarray   # n×n int32 index into rf.events of the pair's LCA node
    xeno: np.ndarray   # n×n bool: pair involves a transferred (recipient) gene vs a gene outside its clade


def _fill_blocks(M: np.ndarray, ch: list[np.ndarray], value) -> None:
    """Write ``value(A, B)`` for all cross-child pairs of a node (both orientations)."""
    if len(ch) > _BLOCK_FANOUT:
        allv = np.concatenate(ch)
        M[np.ix_(allv, allv)] = value(allv, allv)
        return
    for a in range(len(ch)):
        for b in range(a + 1, len(ch)):
            A, B = ch[a], ch[b]
            v = value(A, B)
            M[np.ix_(A, B)] = v
            M[np.ix_(B, A)] = v.T if isinstance(v, np.ndarray) else v


def reference_matrices(rf: ReconciledFamily) -> RefMatrices:
    """Reference merge intervals for every gene pair, written per reference-tree node over
    products of child leaf sets (preorder; descendants overwrite whole-block writes)."""
    st = rf.st
    n = len(rf.genes)
    lo = np.full((n, n), NOT_ASSERTED, dtype=np.int32)
    hi = np.full((n, n), NOT_ASSERTED, dtype=np.int32)
    node = np.full((n, n), -1, dtype=np.int32)
    xeno = np.zeros((n, n), dtype=bool)
    sp = rf.species_idx
    for k, ev in enumerate(rf.events):
        if ev.ev == "S":
            val_lo = val_hi = lambda A, B: st.lca_idx[np.ix_(sp[A], sp[B])]
        elif ev.ev in ("D", "H"):
            val_lo = (lambda v: (lambda A, B: v))(ev.merge_lo)
            val_hi = (lambda v: (lambda A, B: v))(ev.merge_hi)
        else:
            val_lo = val_hi = lambda A, B: NOT_ASSERTED
        _fill_blocks(lo, ev.child_leaves, val_lo)
        _fill_blocks(hi, ev.child_leaves, val_hi)
        _fill_blocks(node, ev.child_leaves, lambda A, B, k=k: k)
    # xenolog pairs: transferred genes vs every gene outside their recipient clade
    for ev in rf.events:
        if ev.ev == "H" and ev.recipient_child is not None:
            rec = ev.child_leaves[ev.recipient_child]
            mask = np.ones(n, dtype=bool)
            mask[rec] = False
            out = np.flatnonzero(mask)
            xeno[np.ix_(rec, out)] = True
            xeno[np.ix_(out, rec)] = True
    return RefMatrices(lo, hi, node, xeno)


def _older(st: SpeciesTree, a: int, b: int) -> int:
    """The older (closer to the root) of two comparable levels; INF is oldest."""
    if a == INF or b == INF:
        return INF
    return a if st.depth[a] <= st.depth[b] else b


def representative_matrix(rf: ReconciledFamily) -> np.ndarray:
    """One consistent merge level per pair, inside every reference interval, used to build
    reference-derived HOGs (positive control, RefOGs). Duplications are placed at the youngest
    level compatible with their descendants (bottom-up), U pairs are grouped at the oldest
    level implied inside the U node (INF if any INF inside); U pairs stay unasserted for
    scoring."""
    st = rf.st
    n = len(rf.genes)
    chosen: dict[int, int] = {}   # event index -> chosen duplication level (D) or merge (H)
    idx_of = {id(e.node): k for k, e in enumerate(rf.events)}
    for k in reversed(range(len(rf.events))):
        e = rf.events[k]
        if e.ev != "D":
            continue
        c = e.lvl_lo
        stack = list(e.node.children)
        while stack:
            nd = stack.pop()
            if nd.is_leaf:
                continue
            ce = rf.events[idx_of[id(nd)]]
            if ce.ev == "D":
                c = _older(st, c, chosen[idx_of[id(nd)]])
            elif ce.ev == "U":
                stack.extend(nd.children)
        if e.lvl_hi != ABOVE and c != INF and not st.anc[e.lvl_hi, c]:
            raise ScoringError(f"[{rf.id}] duplication {e.id}: descendants force a level above its lvl_hi")
        chosen[k] = c
    M = np.full((n, n), NOT_ASSERTED, dtype=np.int32)
    sp = rf.species_idx
    for k, ev in enumerate(rf.events):
        if ev.ev == "S":
            _fill_blocks(M, ev.child_leaves, lambda A, B: st.lca_idx[np.ix_(sp[A], sp[B])])
        elif ev.ev == "D":
            p = chosen[k]
            v = INF if p == INF else (int(st.parent_idx[p]) if st.parent_idx[p] >= 0 else INF)
            _fill_blocks(M, ev.child_leaves, lambda A, B, v=v: v)
        elif ev.ev == "H":
            _fill_blocks(M, ev.child_leaves, lambda A, B, v=ev.merge_lo: v)
        else:
            _fill_blocks(M, ev.child_leaves, lambda A, B: NOT_ASSERTED)
    for ev in reversed(rf.events):  # inner U nodes first
        if ev.ev != "U":
            continue
        allv = np.concatenate(ev.child_leaves)
        sub = M[np.ix_(allv, allv)]
        off = ~np.eye(len(allv), dtype=bool)
        vals = np.unique(sub[off & (sub != NOT_ASSERTED)])
        if (vals == INF).any():
            v = INF
        else:
            cand = np.append(vals, ev.mapping).astype(np.int64)
            v = int(cand[np.argmin(st.depth[cand])])
        ch = ev.child_leaves
        for a in range(len(ch)):
            for b in range(a + 1, len(ch)):
                M[np.ix_(ch[a], ch[b])] = v
                M[np.ix_(ch[b], ch[a])] = v
    return M


def level_components(rf: ReconciledFamily, Mrep: np.ndarray, L: int) -> np.ndarray:
    """Reference-derived HOG labels at level ``L`` (components of "merge ⊆ L" among genes
    below ``L``; -1 for genes not below ``L``)."""
    from scipy.sparse import csr_matrix
    from scipy.sparse.csgraph import connected_components

    st = rf.st
    below = np.flatnonzero(st.anc[L, rf.species_idx])
    lab = np.full(len(rf.genes), -1, dtype=np.int64)
    if len(below) == 0:
        return lab
    sub = Mrep[np.ix_(below, below)]
    adj = (sub >= 0) & st.anc[L, np.where(sub >= 0, sub, 0)]
    np.fill_diagonal(adj, False)
    _, comp = connected_components(csr_matrix(adj), directed=False)
    lab[below] = comp
    return lab


def ancestor_levels(rf: ReconciledFamily) -> set[int]:
    """Internal levels that are ancestors of at least one member species."""
    st = rf.st
    out: set[int] = set()
    for s in set(rf.species_idx.tolist()):
        p = int(st.parent_idx[s])
        while p >= 0:
            out.add(p)
            p = int(st.parent_idx[p])
    return out


def reference_level_codes(rf: ReconciledFamily, levels: Optional[Iterable[int]] = None) -> dict[int, np.ndarray]:
    """Reference-derived HOG labels per internal level (default: all ancestor levels)."""
    Mrep = representative_matrix(rf)
    return {L: level_components(rf, Mrep, L) for L in (ancestor_levels(rf) if levels is None else levels)}


# ---------------------------------------------------------------------- prediction side
def hog_codes(rf: ReconciledFamily, pred: Prediction) -> dict[int, np.ndarray]:
    """Integer HOG code per tree gene for every ancestor level (-1 = no HOG there).
    Raises :class:`ScoringError` if a gene sits in a HOG at a level its species is not below."""
    st = rf.st
    out: dict[int, np.ndarray] = {}
    for L in ancestor_levels(rf):
        m = pred.levels.get(st.names[L], {})
        below = st.anc[L, rf.species_idx]
        codes = np.full(len(rf.genes), -1, dtype=np.int64)
        table: dict[str, int] = {}
        for k, g in enumerate(rf.genes):
            h = m.get(g)
            if h is None:
                continue
            if not below[k]:
                raise ScoringError(f"[{rf.id}] predicted HOG {h!r} puts {g!r} at level {st.names[L]!r}, "
                                   f"which its species is not below")
            codes[k] = table.setdefault(h, len(table))
        out[L] = codes
    # a gene in a HOG at a level that is not an ancestor of its species at all
    anc_names = {st.names[L] for L in out}
    for lvname, m in pred.levels.items():
        if lvname in anc_names:
            continue
        for g in rf.genes:
            if g in m:
                raise ScoringError(f"[{rf.id}] predicted HOG puts {g!r} at level {lvname!r}, "
                                   f"which its species is not below")
    return out


def predicted_merge_matrix(rf: ReconciledFamily, codes: dict[int, np.ndarray]) -> np.ndarray:
    """Lowest level at which each pair shares a predicted HOG (INF if never)."""
    st = rf.st
    n = len(rf.genes)
    M = np.full((n, n), INF, dtype=np.int32)
    for L in sorted(codes, key=lambda x: -st.depth[x]):  # deepest first => lowest level wins
        c = codes[L]
        ok = np.flatnonzero(c >= 0)
        if len(ok) < 2:
            continue
        order = ok[np.argsort(c[ok], kind="stable")]
        bounds = np.flatnonzero(np.diff(c[order])) + 1
        for grp in np.split(order, bounds):
            if len(grp) < 2:
                continue
            sub = M[np.ix_(grp, grp)]
            sub[sub == INF] = L
            M[np.ix_(grp, grp)] = sub
    return M


# ------------------------------------------------------------------- interval scoring
def _safe(x: np.ndarray) -> np.ndarray:
    return np.where(x >= 0, x, 0)


def interval_correct(st: SpeciesTree, p: np.ndarray, lo: np.ndarray, hi: np.ndarray) -> np.ndarray:
    """``lo ⊆ p ⊆ hi`` (vectorised; all arrays of level indices with INF = -1)."""
    pf, lf, hf = p >= 0, lo >= 0, hi >= 0
    below_lo = pf & lf & ~st.anc[_safe(p), _safe(lo)]
    above_hi = hf & (~pf | ~st.anc[_safe(hi), _safe(p)])
    return np.where(lf, ~below_lo & ~above_hi, ~pf)


def interval_error(st: SpeciesTree, p: np.ndarray, lo: np.ndarray, hi: np.ndarray) -> np.ndarray:
    """Signed error in edges to the nearest bound (float, ±inf; 0 when correct)."""
    pf, lf, hf = p >= 0, lo >= 0, hi >= 0
    err = np.zeros(len(p), dtype=float)
    d = st.depth
    below_lo = pf & lf & ~st.anc[_safe(p), _safe(lo)]
    err[below_lo] = (d[_safe(lo)] - d[_safe(p)])[below_lo]
    above_fin = hf & pf & ~st.anc[_safe(hi), _safe(p)]
    err[above_fin] = (d[_safe(hi)] - d[_safe(p)])[above_fin]
    err[hf & ~pf] = np.inf
    err[~lf & pf] = -np.inf
    return err


def _f1(tp: int, fp: int, fn: int) -> Optional[float]:
    """B1: F1 = 2tp/(2tp+fp+fn) whenever anything is counted; None only if tp=fp=fn=0."""
    d = 2 * tp + fp + fn
    return 2 * tp / d if d else None


def _prf(tp: int, fp: int, fn: int) -> dict:
    return {"tp": int(tp), "fp": int(fp), "fn": int(fn),
            "precision": tp / (tp + fp) if tp + fp else None,
            "recall": tp / (tp + fn) if tp + fn else None,
            "f1": _f1(tp, fp, fn)}


def _level_name(st: SpeciesTree, x: int) -> str:
    return "inf" if x == INF else st.names[x]


def _histogram(err: np.ndarray) -> dict[str, int]:
    h: dict[str, int] = {}
    vals, cnts = np.unique(err, return_counts=True)
    for v, c in zip(vals.tolist(), cnts.tolist()):
        if math.isinf(v):
            key = "+inf" if v > 0 else "-inf"
        else:
            key = f"{int(v):+d}" if v else "0"
        h[key] = int(c)
    return h


def _mean(x) -> Optional[float]:
    x = list(x)
    return float(np.mean(x)) if x else None


# ------------------------------------------------------------------- panel context
@dataclass
class PanelContext:
    """Panel-wide information for integrity metrics.

    ``gene_family``: every reference member → family id. ``decoys``: family → declared
    non-homologous gene IDs. ``distinct``: family → {other family: at_level index or INF (=all
    levels, non-homologous)}; symmetric. ``count_species``: species whose per-HOG counts the
    parsers must record (union of ``complete_species``)."""
    gene_family: dict[str, str]
    decoys: dict[str, set[str]] = field(default_factory=dict)
    distinct: dict[str, dict[str, int]] = field(default_factory=dict)
    count_species: set[str] = field(default_factory=set)


def panel_context(families: Sequence[ReferenceFamily], st: Optional[SpeciesTree] = None) -> PanelContext:
    gf: dict[str, str] = {}
    dec: dict[str, set[str]] = {}
    distinct: dict[str, dict[str, int]] = defaultdict(dict)
    cs: set[str] = set()
    for f in families:
        for g in f.member_ids:
            gf.setdefault(g, f.id)
        dec[f.id] = set(f.decoys)
        cs.update(f.complete_species)
        for d in f.distinct_from:
            lv = str(d["at_level"])
            if lv.upper() in NONHOMOLOGOUS_TOKENS:
                v = INF
            else:
                if st is None:
                    raise ScoringError("panel_context: a species tree is needed to resolve distinct_from levels")
                v = st.idx(lv)
            other = str(d["family"])
            for a, b in ((f.id, other), (other, f.id)):
                prev = distinct[a].get(b)
                # keep the most permissive-to-penalise (oldest) level; INF = all levels
                distinct[a][b] = v if prev is None else (INF if INF in (prev, v) else
                                                         (v if st.depth[v] < st.depth[prev] else prev))
    return PanelContext(gf, dec, dict(distinct), cs)


# --------------------------------------------------------------------- family scoring
def evaluate_family(
    fam: ReferenceFamily | ReconciledFamily,
    st: SpeciesTree,
    pred: Prediction,
    ctx: Optional[PanelContext] = None,
    thetas: Sequence[float] = THETAS,
    per_gene: bool = True,
    compact: bool = False,
) -> dict:
    """Score one reference family. ``compact`` drops per-gene/per-species/per-node detail
    (used for units of genome-scale collections)."""
    rf = fam if isinstance(fam, ReconciledFamily) else reconcile(fam, st)
    F = rf.family
    ctx = ctx or panel_context([F], st)
    n = len(rf.genes)
    R = reference_matrices(rf)
    codes = hog_codes(rf, pred)
    Mpred = predicted_merge_matrix(rf, codes)
    gsp = F.gene_species
    mism = [g for g in rf.genes if pred.gene_species.get(g) and pred.gene_species[g] != gsp[g]]
    if mism:
        raise ScoringError(f"[{F.id}] {len(mism)} gene(s) have a different species in the prediction, "
                           f"e.g. {mism[:3]} ({pred.gene_species[mism[0]]} vs {gsp[mism[0]]})")

    iu, ju = np.triu_indices(n, 1)
    sel = R.lo[iu, ju] != NOT_ASSERTED
    i, j = iu[sel], ju[sel]
    lo = R.lo[i, j].astype(np.int64)
    hi = R.hi[i, j].astype(np.int64)
    p = Mpred[i, j].astype(np.int64)
    node = R.node[i, j]
    s0 = st.lca_idx[rf.species_idx[i], rf.species_idx[j]]
    correct = interval_correct(st, p, lo, hi)
    err = interval_error(st, p, lo, hi)
    exact = lo == hi
    npairs = len(i)
    fam_depth = max(1, int(max(st.depth[s] for s in set(rf.species_idx.tolist())) - st.depth[rf.family_level])) \
        if n else 1

    res: dict = {
        "id": F.id, "name": F.name, "panel": F.panel, "collection": F.collection, "phenomena": F.phenomena,
        "family_level": st.names[rf.family_level], "family_depth": fam_depth,
        "n_members": len(F.member_ids), "n_tree_genes": n, "n_pairs_total": n * (n - 1) // 2,
        "n_asserted_pairs": int(npairs), "n_exact_pairs": int(exact.sum()), "n_interval_pairs": int((~exact).sum()),
    }
    absent = [g for g in F.member_ids if not any(g in m for m in pred.levels.values())]
    res["n_members_absent_from_prediction"] = len(absent)

    # ---- 1. merge-level accuracy: pair-pooled, node-balanced (B2), exact vs interval, tolerance (B12)
    tol1 = correct | (np.isfinite(err) & (np.abs(err) <= 1))
    # pairs grouped by their reference LCA node (one sort instead of a scan per node)
    order = np.argsort(node, kind="stable")
    node_sorted = node[order]
    node_ids, starts, counts = np.unique(node_sorted, return_index=True, return_counts=True)
    pairs_of = {int(k): order[s_:s_ + c_] for k, s_, c_ in zip(node_ids, starts, counts)}
    nb = len(rf.events)
    n_by = np.bincount(node, minlength=nb)
    acc_by = np.bincount(node, weights=correct, minlength=nb)
    tol_by = np.bincount(node, weights=tol1, minlength=nb)
    node_acc = [float(acc_by[k] / n_by[k]) for k in node_ids]
    node_tol = [float(tol_by[k] / n_by[k]) for k in node_ids]
    nodes_detail = [{"node": rf.events[k].id, "ev": rf.events[k].ev, "n_pairs": int(n_by[k]),
                     "accuracy": float(acc_by[k] / n_by[k]), "interval": bool(rf.events[k].is_interval)}
                    for k in node_ids]
    fin = np.isfinite(err)
    res["pair_merge"] = {
        "n_pairs": npairs,
        "accuracy_node_balanced": _mean(node_acc),
        "accuracy_pairs": float(correct.mean()) if npairs else None,
        "accuracy_exact_pairs": float(correct[exact].mean()) if exact.any() else None,
        "accuracy_interval_pairs": float(correct[~exact].mean()) if (~exact).any() else None,
        "n_nodes": len(node_acc),
        "error_histogram": _histogram(err),
        "frac_over_split": float((err > 0).mean()) if npairs else None,
        "frac_over_merged": float((err < 0).mean()) if npairs else None,
        "mean_abs_error_finite": float(np.abs(err[fin]).mean()) if fin.any() else None,
        "mean_abs_depth_frac_error": float((np.abs(err[fin]) / fam_depth).mean()) if fin.any() else None,
        "accuracy_tol1_pairs": float(tol1.mean()) if npairs else None,
        "accuracy_tol1_node_balanced": _mean(node_tol),
    }
    if not compact:
        res["nodes"] = nodes_detail

    # ---- 2. per-level: pair P/R/F1 (B1) and B-cubed (B2); levels inside an interval skipped per pair
    per_level = []
    tot = np.zeros(3, dtype=np.int64)
    for L in sorted(codes, key=lambda x: st.depth[x]):
        if not st.anc[rf.family_level, L]:
            continue
        inL = st.anc[L, s0]
        if not inL.any():
            continue
        ref_same = (hi >= 0) & st.anc[L, _safe(hi)]
        ref_apart = (lo < 0) | ~st.anc[L, _safe(lo)]
        scored = inL & (ref_same | ref_apart)
        c = codes[L]
        ii, jj = i[scored], j[scored]
        rs = ref_same[scored]
        ps = (c[ii] == c[jj]) & (c[ii] >= 0)
        tp = int((rs & ps).sum()); fp = int((~rs & ps).sum()); fn = int((rs & ~ps).sum())
        tot += (tp, fp, fn)
        d = {"level": st.names[L], "depth": int(st.depth[L]), "n_pairs": int(scored.sum()),
             "n_skipped_interval": int((inL & ~scored).sum()), "n_ref_same": int(rs.sum()), "n_pred_same": int(ps.sum())}
        d.update(_prf(tp, fp, fn))
        # B-cubed over genes with ≥1 scored pair at L (self counted in both clusters)
        if len(ii):
            tpg = np.bincount(ii, weights=rs & ps, minlength=n) + np.bincount(jj, weights=rs & ps, minlength=n)
            ppg = np.bincount(ii, weights=ps, minlength=n) + np.bincount(jj, weights=ps, minlength=n)
            rpg = np.bincount(ii, weights=rs, minlength=n) + np.bincount(jj, weights=rs, minlength=n)
            has = (np.bincount(ii, minlength=n) + np.bincount(jj, minlength=n)) > 0
            bp = float(((tpg + 1) / (ppg + 1))[has].mean())
            br = float(((tpg + 1) / (rpg + 1))[has].mean())
            d.update({"bcubed_precision": bp, "bcubed_recall": br, "bcubed_f1": 2 * bp * br / (bp + br),
                      "n_genes": int(has.sum())})
        else:
            d.update({"bcubed_precision": None, "bcubed_recall": None, "bcubed_f1": None, "n_genes": 0})
        per_level.append(d)
    res["per_level"] = per_level
    res["per_level_summary"] = {
        "mean_bcubed_f1": _mean(d["bcubed_f1"] for d in per_level if d["bcubed_f1"] is not None),
        "mean_pair_f1": _mean(d["f1"] for d in per_level if d["f1"] is not None),
        "n_levels": len(per_level), "pooled": _prf(*tot)}

    # ---- 3. within-family orthology P/R (B13): cross-species, non-xenolog pairs
    xen = R.xeno[i, j]
    cross = (rf.species_idx[i] != rf.species_idx[j]) & ~xen
    ro = ((lo == s0) & (hi == s0))[cross]
    po = (p == s0)[cross]
    res["within_family_orthology"] = {
        "note": "negatives are within-family paralogs only; not QfO-comparable",
        "n_pairs": int(cross.sum()), "n_xenolog_pairs_excluded": int(xen.sum()),
        "n_ref_orthologs": int(ro.sum()), "n_ref_paralogs": int((~ro).sum()),
        "ortholog": _prf(int((ro & po).sum()), int((~ro & po).sum()), int((ro & ~po).sum())),
        "paralog": _prf(int((~ro & ~po).sum()), int((ro & ~po).sum()), int((~ro & po).sum())),
        "accuracy": float((ro == po).mean()) if cross.any() else None,
    }

    # ---- 4. event recovery, two-sided (B3): recovered@θ iff ≥θ of cross-copy pairs fall in the
    #      reference merge interval (separation below + merge side above) and every copy keeps
    #      ≥θ of its reference-together pairs together at the duplication level (cohesion)
    events = []
    for k, ev in enumerate(rf.events):
        if ev.ev not in ("D", "H"):
            continue
        m = pairs_of.get(k, np.array([], dtype=np.int64))
        if len(m) == 0:
            continue
        pv, lov, hiv, ev_err = p[m], lo[m], hi[m], err[m]
        ev_corr = correct[m]
        D = ev.lvl_lo
        # separation: copies apart at the duplication level (younger bound)
        sep = float(((pv < 0) | ~st.anc[D, _safe(pv)]).mean())
        # cohesion: within each copy, pairs the reference puts together at D are together at D
        cohesions = []
        for leaves in ev.child_leaves:
            if len(leaves) < 2 or st.is_leaf[D] or D not in codes:
                cohesions.append(None)
                continue
            a, b = np.triu_indices(len(leaves), 1)
            ga, gb = leaves[a], leaves[b]
            hv = R.hi[ga, gb].astype(np.int64)
            together = (R.lo[ga, gb] != NOT_ASSERTED) & (hv >= 0) & st.anc[D, _safe(hv)]
            if not together.any():
                cohesions.append(None)
                continue
            c = codes[D]
            cohesions.append(float(((c[ga] == c[gb]) & (c[ga] >= 0))[together].mean()))
        min_coh = min([x for x in cohesions if x is not None], default=1.0)
        # modal predicted merge level; ties -> smallest |error|, then smallest level index (B13)
        cnt = Counter(pv.tolist())
        top = max(cnt.values())
        tied = [v for v, c in cnt.items() if c == top]
        e_of = {v: float(interval_error(st, np.array([v]), lov[:1], hiv[:1])[0]) for v in tied}
        modal = sorted(tied, key=lambda v: (abs(e_of[v]), v))[0]
        mis = e_of[modal]
        d = {
            "event": ev.id, "ev": ev.ev, "level": st.names[ev.lvl_lo],
            "level_hi": "ABOVE_RANGE" if ev.lvl_hi == ABOVE else st.names[ev.lvl_hi],
            "level_source": ev.level_source, "ref_merge_lo": _level_name(st, ev.merge_lo),
            "ref_merge_hi": _level_name(st, ev.merge_hi), "tags": ev.tags, "src": ev.src,
            "n_copies": len(ev.child_leaves), "n_pairs": int(len(m)),
            "frac_correct": float(ev_corr.mean()), "separation": sep,
            "cohesion_min": min_coh, "cohesion": cohesions,
            "frac_over_split": float((ev_err > 0).mean()), "frac_over_merged": float((ev_err < 0).mean()),
            "modal_pred_merge": _level_name(st, modal),
            "misplacement": ("+inf" if mis > 0 else "-inf") if math.isinf(mis) else int(mis),
        }
        d["merge_side"] = float((((hiv < 0) | ((pv >= 0) & st.anc[_safe(hiv), _safe(pv)]))).mean())
        # two-sided: cross-copy pairs inside the merge interval (apart below it AND merged by its
        # upper bound) and every copy cohesive at the duplication level
        for t in thetas:
            d[f"recovered_{t}"] = bool(d["frac_correct"] >= t and min_coh >= t)
        events.append(d)
    res["events"] = events
    res["events_summary"] = {"n_events": len(events),
                             "mean_frac_correct": _mean(e["frac_correct"] for e in events)}
    for t in thetas:
        res["events_summary"][f"recovery_rate_{t}"] = _mean(e[f"recovered_{t}"] for e in events)

    # ---- spurious duplications inside pure-speciation reference clades (B3 / ILS signal)
    res["spurious_duplications"] = _spurious_duplications(rf, codes)

    # ---- 5. family integrity (B5–B7)
    res["integrity"] = family_integrity(rf, st, pred, ctx)

    # ---- 6. per-gene / per-species correctness
    n_pairs_g = np.bincount(i, minlength=n) + np.bincount(j, minlength=n)
    n_corr_g = np.bincount(i, weights=correct, minlength=n) + np.bincount(j, weights=correct, minlength=n)
    with np.errstate(invalid="ignore", divide="ignore"):
        frac_g = np.where(n_pairs_g > 0, n_corr_g / np.maximum(n_pairs_g, 1), np.nan)
    if per_gene and not compact:
        res["genes"] = [{"gene": g, "species": gsp[g], "n_pairs": int(n_pairs_g[k]),
                         "frac_correct": None if np.isnan(frac_g[k]) else float(frac_g[k])}
                        for k, g in enumerate(rf.genes)]
    if not compact:
        by_sp: dict[str, list[float]] = defaultdict(list)
        for k, g in enumerate(rf.genes):
            if not np.isnan(frac_g[k]):
                by_sp[gsp[g]].append(float(frac_g[k]))
        res["species"] = {s: {"n_genes": len(v), "mean_frac_correct": float(np.mean(v)),
                              "n_perfect": int(sum(x == 1.0 for x in v)),
                              "n_majority": int(sum(x >= 0.5 for x in v))} for s, v in sorted(by_sp.items())}
    valid = frac_g[~np.isnan(frac_g)]
    res["genes_summary"] = {"n_genes_with_pairs": int(len(valid)),
                            "mean_frac_correct": float(valid.mean()) if len(valid) else None,
                            "frac_perfect": float((valid == 1.0).mean()) if len(valid) else None}
    return res


def _spurious_duplications(rf: ReconciledFamily, codes: dict[int, np.ndarray]) -> dict:
    """For every maximal reference clade made only of speciation nodes, count predicted
    splits: at each internal level L inside the clade, (#predicted HOGs among the clade's genes
    below L) − 1, genes without a HOG counting as their own HOG."""
    st = rf.st
    pure: dict[int, bool] = {}
    idx_of = {id(e.node): e for e in rf.events}
    for e in reversed(rf.events):
        pure[id(e.node)] = e.ev == "S" and all(c.is_leaf or pure.get(id(c), False) for c in e.node.children)
    count = opp = clades = 0
    for e in rf.events:
        if not pure[id(e.node)]:
            continue
        par = e.node.parent
        if par is not None and pure.get(id(par), False):
            continue
        clades += 1
        genes = np.concatenate(e.child_leaves)
        for L, c in codes.items():
            if not st.anc[e.mapping, L]:
                continue
            g = genes[st.anc[L, rf.species_idx[genes]]]
            if len(g) < 2:
                continue
            lab = c[g]
            k = len(set(lab[lab >= 0].tolist())) + int((lab < 0).sum())
            count += k - 1
            opp += len(g) - 1
    return {"n_clades": clades, "count": count, "opportunities": opp,
            "rate": count / opp if opp else None}


# ------------------------------------------------------------------------- integrity
def family_integrity(rf: ReconciledFamily, st: SpeciesTree, pred: Prediction, ctx: PanelContext) -> dict:
    """Root-level measures at the family level ``L_F`` (B5–B7).

    *RefOGs* = reference-derived HOGs at ``L_F`` (members not placed in the tree join the
    single RefOG if there is one, else they are excluded and counted). *PredOGs* = each
    member's predicted HOG at the highest level ≤ ``L_F`` containing it; a member in no HOG of
    size ≥ 2 is *missing* (as in OrthoBench). Judgeable non-own genes in a PredOG (false
    positives) are: members of another RefOG of the same family, contaminants (declared decoys;
    members of a ``distinct_from`` family when the HOG's level ⊆ its ``at_level``) and
    non-members from declared-complete species. Other genes are *unlabelled* and reported only.

    Returns OF3/OrthoBench-style measures: membership precision/recall/F1 (overlap-weighted per
    RefOG, size-weighted over RefOGs), OrthoBench pairwise P/R/F ((n−1)-normalised), % fission,
    % fusion, % exact, % missing genes, size-weighted normalised entropy, contamination =
    contaminants / (own + contaminants), unsigned split/merge fractions.
    """
    F = rf.family
    L_F = rf.family_level
    Lname = st.names[L_F]
    gsp = F.gene_species
    members = F.member_ids
    member_set = set(members)
    # RefOGs
    lab = reference_level_codes(rf, [L_F])[L_F] if rf.genes else np.array([], dtype=np.int64)
    refog_of: dict[str, int] = {g: int(lab[k]) for k, g in enumerate(rf.genes) if lab[k] >= 0}
    unplaced = [g for g in members if g not in rf.gene_index]
    n_refogs = len(set(refog_of.values()))
    excluded = []
    for g in unplaced:
        if n_refogs == 1:
            refog_of[g] = next(iter(refog_of.values()))
        elif n_refogs == 0:
            refog_of[g] = 0
        else:
            excluded.append(g)
    refogs: dict[int, list[str]] = defaultdict(list)
    for g, r_ in refog_of.items():
        refogs[r_].append(g)

    # PredOG of each member
    def top_hog(g: str) -> Optional[tuple[str, str]]:
        top = None
        for lv in st.ancestors(gsp[g], include_self=False):
            if not st.anc[L_F, st.index[lv]]:
                break
            h = pred.hog(lv, g)
            if h is not None:
                top = (lv, h)
        if top is not None and pred.hog_size(*top) < 2:
            return None
        return top
    # a member of a single-gene RefOG that is alone (no HOG / singleton HOG) is exactly right:
    # it is its own PredOG; for larger RefOGs such members are "missing" (as in OrthoBench)

    key_of = {g: top_hog(g) for g in refog_of}
    decoys = ctx.decoys.get(F.id, set()) | set(F.decoys)
    distinct = ctx.distinct.get(F.id, {})
    complete = set(F.complete_species)

    pog: dict[tuple, dict] = {}
    for key in {k for k in key_of.values() if k is not None}:
        lv, h = key
        L = st.index[lv]
        own_by_refog: Counter = Counter()
        other_refog = contam = 0
        contam_species: Counter = Counter()
        n_member_any = 0
        n_other_family = 0
        fam_members_by_species: Counter = Counter()
        for g in pred.hog_members(lv, h):
            fam = ctx.gene_family.get(g)
            sp = pred.gene_species.get(g)
            if g in member_set:
                fam_members_by_species[gsp[g]] += 1
                if g in refog_of:
                    own_by_refog[refog_of[g]] += 1
                n_member_any += 1
            elif g in decoys or (fam is not None and fam in distinct and
                                 (distinct[fam] == INF or st.anc[distinct[fam], L])):
                contam += 1
                if sp:
                    contam_species[sp] += 1
            elif fam is not None:
                n_other_family += 1
        complete_nonmembers = 0
        counts = pred.hog_species.get(lv, {}).get(h, {})
        for s in complete:
            if s not in pred.count_species:
                raise ScoringError(f"[{F.id}] complete species {s!r} was not counted by the parser "
                                   f"(pass count_species)")
            complete_nonmembers += counts.get(s, 0) - fam_members_by_species.get(s, 0) - contam_species.get(s, 0)
        if complete_nonmembers < 0:
            raise ScoringError(f"[{F.id}] negative non-member count in HOG {h!r} (inconsistent species counts)")
        size = pred.hog_size(lv, h)
        pog[key] = {"own": own_by_refog, "n_members": n_member_any, "contam": contam,
                    "complete_nonmembers": complete_nonmembers, "other_family": n_other_family, "size": size}

    def fp_for(key, r_) -> int:
        d = pog[key]
        return (d["n_members"] - d["own"][r_]) + d["contam"] + d["complete_nonmembers"]

    # per RefOG measures
    tot_w = 0
    PS = RC = H = 0.0
    w_ent = 0
    n_fission = n_fusion = n_exact = 0
    n_missing = 0
    ob_tp = ob_fp = ob_fn = 0.0
    split_off = assigned = 0
    for r_, genes in refogs.items():
        nR = len(genes)
        keys = Counter(key_of[g] for g in genes)
        missing = keys.pop(None, 0)
        if nR == 1 and missing:
            k1 = ("<single>", genes[0])
            pog[k1] = {"own": Counter({r_: 1}), "n_members": 1, "contam": 0, "complete_nonmembers": 0,
                       "other_family": 0, "size": 1}
            keys, missing = Counter({k1: 1}), 0
        n_missing += missing
        overl = {k: c for k, c in keys.items()}
        tot_w += nR
        if overl:
            s_t = sum(overl.values())
            PS_i = sum(t * (t / (t + fp_for(k, r_))) for k, t in overl.items()) / s_t
            RC_i = sum(t * t for t in overl.values()) / (nR * s_t)
            PS += nR * PS_i
            RC += nR * RC_i
            split_off += s_t - max(overl.values())
            assigned += s_t
        if len(overl) > 1:
            n_fission += 1
        if any((pog[k]["n_members"] - pog[k]["own"][r_]) + pog[k]["contam"] > 0 for k in overl):
            n_fusion += 1
        if missing == 0 and len(overl) == 1 and fp_for(next(iter(overl)), r_) == 0:
            n_exact += 1
        if nR >= 2:
            probs = [c / nR for c in overl.values()] + ([missing / nR] if missing else [])
            Hi = -sum(q * math.log2(q) for q in probs if q > 0)
            H += nR * Hi / math.log2(nR)
            w_ent += nR
            # OrthoBench benchmark.py (q_even): per-RefOG pair counts normalised by n-1
            tp = sum(t * (t - 1) / 2 for t in overl.values())
            fp = sum(t * fp_for(k, r_) for k, t in overl.items())
            fn = (sum((nR - t) * t for t in overl.values()) + missing * (nR - 1)) / 2
            ob_tp += tp / (nR - 1); ob_fp += fp / (nR - 1); ob_fn += fn / (nR - 1)
    S = len(refogs)
    mp = PS / tot_w if tot_w and assigned else None
    mr = RC / tot_w if tot_w else None
    if mp is None:
        mf1 = 0.0 if tot_w else None
    else:
        mf1 = 2 * mp * mr / (mp + mr) if (mp + mr) else 0.0
    own_total = sum(d["n_members"] for d in pog.values())
    contam_total = sum(d["contam"] for d in pog.values())
    judge_fp = sum(d["contam"] + d["complete_nonmembers"] for d in pog.values())
    ob_p = ob_tp / (ob_tp + ob_fp) if ob_tp + ob_fp else None
    ob_r = ob_tp / (ob_tp + ob_fn) if ob_tp + ob_fn else None
    return {
        "level": Lname, "n_members": len(members), "n_refogs": S, "n_unplaced_excluded": len(excluded),
        "n_predogs": len(pog),
        "membership_precision": mp, "membership_recall": mr, "membership_f1": mf1,
        "orthobench_precision": ob_p, "orthobench_recall": ob_r,
        "orthobench_f": (2 * ob_p * ob_r / (ob_p + ob_r) if (ob_p and ob_r) else (0.0 if ob_p is not None or ob_r is not None else None)),
        "fission_pct": 100 * n_fission / S if S else None,
        "fusion_pct": 100 * n_fusion / S if S else None,
        "exact_pct": 100 * n_exact / S if S else None,
        "missing_genes_pct": 100 * n_missing / tot_w if tot_w else None,
        "entropy": H / w_ent if w_ent else None,
        "contamination": contam_total / (own_total + contam_total) if own_total + contam_total else None,
        "n_contaminants": contam_total,
        "n_complete_species_nonmembers": sum(d["complete_nonmembers"] for d in pog.values()),
        "n_other_reference_family_genes": sum(d["other_family"] for d in pog.values()),
        "n_unlabelled": sum(d["size"] - d["n_members"] - d["contam"] - d["complete_nonmembers"] for d in pog.values()),
        "frac_members_split_off": split_off / assigned if assigned else None,
        "frac_over_merged": judge_fp / (own_total + judge_fp) if own_total + judge_fp else None,
    }


# --------------------------------------------------------------------------- scalars
PRIMARY = ["level_bcubed_f1", "merge_accuracy_node", "event_recovery_0.5", "membership_f1"]

SCALARS = {
    # primary
    "level_bcubed_f1": ("per_level_summary", "mean_bcubed_f1"),
    "merge_accuracy_node": ("pair_merge", "accuracy_node_balanced"),
    "event_recovery_0.5": ("events_summary", "recovery_rate_0.5"),
    "membership_f1": ("integrity", "membership_f1"),
    # secondary / sensitivity
    "event_recovery_0.8": ("events_summary", "recovery_rate_0.8"),
    "merge_accuracy_pairs": ("pair_merge", "accuracy_pairs"),
    "merge_accuracy_exact": ("pair_merge", "accuracy_exact_pairs"),
    "merge_accuracy_interval": ("pair_merge", "accuracy_interval_pairs"),
    "merge_accuracy_tol1": ("pair_merge", "accuracy_tol1_node_balanced"),
    "mean_abs_depth_frac_error": ("pair_merge", "mean_abs_depth_frac_error"),
    "frac_over_split": ("pair_merge", "frac_over_split"),
    "frac_over_merged": ("pair_merge", "frac_over_merged"),
    "level_pair_f1_mean": ("per_level_summary", "mean_pair_f1"),
    "level_pair_f1_pooled": ("per_level_summary", "pooled", "f1"),
    "wf_ortholog_precision": ("within_family_orthology", "ortholog", "precision"),
    "wf_ortholog_recall": ("within_family_orthology", "ortholog", "recall"),
    "wf_ortholog_f1": ("within_family_orthology", "ortholog", "f1"),
    "wf_paralog_f1": ("within_family_orthology", "paralog", "f1"),
    "event_mean_frac_correct": ("events_summary", "mean_frac_correct"),
    "spurious_dup_rate": ("spurious_duplications", "rate"),
    "membership_precision": ("integrity", "membership_precision"),
    "membership_recall": ("integrity", "membership_recall"),
    "orthobench_f": ("integrity", "orthobench_f"),
    "fission_pct": ("integrity", "fission_pct"),
    "fusion_pct": ("integrity", "fusion_pct"),
    "exact_pct": ("integrity", "exact_pct"),
    "missing_genes_pct": ("integrity", "missing_genes_pct"),
    "entropy": ("integrity", "entropy"),
    "contamination": ("integrity", "contamination"),
    "frac_members_split_off": ("integrity", "frac_members_split_off"),
    "integrity_frac_over_merged": ("integrity", "frac_over_merged"),
    "gene_mean_frac_correct": ("genes_summary", "mean_frac_correct"),
    "gene_frac_perfect": ("genes_summary", "frac_perfect"),
}


def family_scalars(res: dict) -> dict[str, Optional[float]]:
    """Flat dict of every scalar in ``SCALARS`` (None when undefined) plus skill scores."""
    out: dict[str, Optional[float]] = {}
    for name, path in SCALARS.items():
        v = res
        for k in path:
            v = v.get(k) if isinstance(v, dict) else None
        out[name] = None if v is None else float(v)
    for k, v in (res.get("skill") or {}).items():
        out[f"skill_{k}"] = v
    return out
