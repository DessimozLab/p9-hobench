"""Hierarchical-orthology metrics (core_spec.md, metrics 1–6) and aggregation.

Everything is driven by *merge levels* (species-tree node indices; ``INF = -1`` = never in
one HOG; ``NOT_ASSERTED = -2`` = pair below an unresolved ``U`` node):

* reference: from the reconciled reference tree, computed per tree node over products of
  child leaf sets (``reference_merge_matrix``) — O(n²) writes, no Python loop over pairs;
* prediction: lowest level ``L ⊇ lca(sp(a), sp(b))`` at which both genes share a HOG,
  computed level by level (deepest first) over predicted HOG blocks
  (``predicted_merge_matrix``).

Sign convention for level errors: ``error = depth(ref merge) − depth(pred merge)`` in
species-tree edges, so **positive = predicted merge level too old** (genes kept apart for
too long: *over-splitting*, e.g. missed orthology) and **negative = too young** (genes merged
too early: *over-merging*, e.g. pseudo-orthologs called orthologs). ``+inf``: reference merges
but prediction never does; ``-inf``: reference never merges (duplication at the family root)
but prediction does. (core_spec.md's parenthetical "(under-split)/(over-split)" labels are
swapped relative to this; the sign convention itself is the spec's — see ENGINEERING_NOTES.)
"""
from __future__ import annotations

import math
from collections import Counter, defaultdict
from dataclasses import dataclass
from typing import Iterable, Optional, Sequence

import numpy as np

from .predictions import Prediction
from .reference import (INF, NOT_ASSERTED, PER_PAIR_S0, ReconciledFamily, ReferenceFamily,
                        reconcile)
from .seeds import SEEDS
from .speciestree import SpeciesTree

_BLOCK_FANOUT = 16  # above this many children, write a node's whole block at once


# ---------------------------------------------------------------------- merge matrices
def reference_merge_matrix(rf: ReconciledFamily) -> np.ndarray:
    """n×n matrix of reference merge levels (diagonal meaningless)."""
    st = rf.st
    n = len(rf.genes)
    M = np.full((n, n), NOT_ASSERTED, dtype=np.int64)
    sp = rf.species_idx
    for ev in rf.events:  # preorder: descendants overwrite blocks written in "whole block" mode
        ch = ev.child_leaves

        def value(A, B):
            if ev.ev == "S":
                return st.lca_idx[np.ix_(sp[A], sp[B])]
            if ev.ev == "D":
                return ev.merge
            return NOT_ASSERTED

        if len(ch) > _BLOCK_FANOUT:
            allv = np.concatenate(ch)
            M[np.ix_(allv, allv)] = value(allv, allv)
        else:
            for a in range(len(ch)):
                for b in range(a + 1, len(ch)):
                    A, B = ch[a], ch[b]
                    v = value(A, B)
                    M[np.ix_(A, B)] = v
                    M[np.ix_(B, A)] = v.T if isinstance(v, np.ndarray) else v
    return M


def hog_codes(rf: ReconciledFamily, pred: Prediction) -> dict[int, np.ndarray]:
    """For every internal level that is an ancestor of some member species: integer HOG code
    per tree gene (-1 = no HOG at that level, or species not below the level)."""
    st = rf.st
    out: dict[int, np.ndarray] = {}
    anc_levels = set()
    for s in set(rf.species_idx.tolist()):
        p = int(st.parent_idx[s])
        while p >= 0:
            anc_levels.add(p)
            p = int(st.parent_idx[p])
    for L in anc_levels:
        lvname = st.names[L]
        m = pred.levels.get(lvname, {})
        below = st.anc[L, rf.species_idx]
        codes = np.full(len(rf.genes), -1, dtype=np.int64)
        table: dict[str, int] = {}
        for k, g in enumerate(rf.genes):
            if not below[k]:
                continue
            h = m.get(g)
            if h is not None:
                codes[k] = table.setdefault(h, len(table))
        out[L] = codes
    return out


def predicted_merge_matrix(rf: ReconciledFamily, codes: dict[int, np.ndarray]) -> np.ndarray:
    st = rf.st
    n = len(rf.genes)
    M = np.full((n, n), INF, dtype=np.int64)
    for L in sorted(codes, key=lambda x: -st.depth[x]):  # deepest first => lowest level wins
        c = codes[L]
        ok = np.flatnonzero(c >= 0)
        if len(ok) < 2:
            continue
        order = ok[np.argsort(c[ok], kind="stable")]
        cs = c[order]
        bounds = np.flatnonzero(np.diff(cs)) + 1
        for grp in np.split(order, bounds):
            if len(grp) < 2:
                continue
            sub = M[np.ix_(grp, grp)]
            sub[sub == INF] = L
            M[np.ix_(grp, grp)] = sub
    return M


def fill_unresolved(rf: ReconciledFamily, M: np.ndarray) -> np.ndarray:
    """Copy of the reference merge matrix where pairs under a ``U`` node get a *grouping*
    level: the oldest of the U node's LCA mapping and every asserted merge level inside its
    subtree (``INF`` if any inside is ``INF``). Used only to build reference-derived HOGs
    (positive control, expected integrity); these pairs remain unasserted for scoring."""
    st = rf.st
    M = M.copy()
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
            cand = np.append(vals, ev.mapping)
            v = int(cand[np.argmin(st.depth[cand])])
        ch = ev.child_leaves
        for a in range(len(ch)):
            for b in range(a + 1, len(ch)):
                M[np.ix_(ch[a], ch[b])] = v
                M[np.ix_(ch[b], ch[a])] = v
    return M


def level_components(rf: ReconciledFamily, Mfilled: np.ndarray, L: int) -> np.ndarray:
    """Reference-derived HOG labels at level ``L``: connected components of "together at L"
    (merge level ⊆ L) among genes whose species is below ``L``; -1 for other genes."""
    from scipy.sparse import csr_matrix
    from scipy.sparse.csgraph import connected_components

    st = rf.st
    below = np.flatnonzero(st.anc[L, rf.species_idx])
    lab = np.full(len(rf.genes), -1, dtype=np.int64)
    if len(below) == 0:
        return lab
    sub = Mfilled[np.ix_(below, below)]
    adj = (sub >= 0) & st.anc[L, np.where(sub >= 0, sub, 0)]
    np.fill_diagonal(adj, False)
    _, comp = connected_components(csr_matrix(adj), directed=False)
    lab[below] = comp
    return lab


def reference_level_codes(rf: ReconciledFamily, levels: Optional[Iterable[int]] = None) -> dict[int, np.ndarray]:
    """Reference-derived HOG labels for every internal ancestor level of the member species."""
    st = rf.st
    if levels is None:
        levels = set()
        for s in set(rf.species_idx.tolist()):
            p = int(st.parent_idx[s])
            while p >= 0:
                levels.add(p)
                p = int(st.parent_idx[p])
    Mf = fill_unresolved(rf, reference_merge_matrix(rf))
    return {L: level_components(rf, Mf, L) for L in levels}


# --------------------------------------------------------------------------- helpers
def _prf(tp: int, fp: int, fn: int) -> dict:
    p = tp / (tp + fp) if tp + fp else None
    r = tp / (tp + fn) if tp + fn else None
    f = (2 * p * r / (p + r) if (p + r) else 0.0) if (p is not None and r is not None) else None
    return {"tp": int(tp), "fp": int(fp), "fn": int(fn), "precision": p, "recall": r, "f1": f}


def _level_name(st: SpeciesTree, x: int) -> str:
    return "inf" if x == INF else st.names[x]


def _signed_errors(st: SpeciesTree, ref: np.ndarray, pred: np.ndarray):
    """Return (finite error array for both-finite pairs, n_plus_inf, n_minus_inf)."""
    fin = (ref >= 0) & (pred >= 0)
    err = st.depth[ref[fin]] - st.depth[pred[fin]]
    plus_inf = int(((ref >= 0) & (pred == INF)).sum())
    minus_inf = int(((ref == INF) & (pred >= 0)).sum())
    return err, plus_inf, minus_inf


def _histogram(err: np.ndarray, plus_inf: int, minus_inf: int, both_inf: int = 0) -> dict[str, int]:
    """Signed level-error histogram; "0" counts every exact agreement (including ∞ = ∞)."""
    h: dict[str, int] = {}
    if minus_inf:
        h["-inf"] = minus_inf
    cnt = Counter(err.tolist())
    if both_inf:
        cnt[0] += both_inf
    for v, c in sorted(cnt.items()):
        h[f"{v:+d}" if v else "0"] = int(c)
    if plus_inf:
        h["+inf"] = plus_inf
    return h


# ------------------------------------------------------------------- family evaluation
@dataclass
class PanelContext:
    """Panel-wide information needed for the contamination metric."""
    gene_family: dict[str, str]      # every reference member -> family id
    decoys: set[str]                 # declared non-homologous decoys (any family)


def panel_context(families: Sequence[ReferenceFamily]) -> PanelContext:
    gf: dict[str, str] = {}
    dec: set[str] = set()
    for f in families:
        for g in f.member_ids:
            gf.setdefault(g, f.id)
        dec.update(f.decoys)
    return PanelContext(gf, dec)


def evaluate_family(
    fam: ReferenceFamily | ReconciledFamily,
    st: SpeciesTree,
    pred: Prediction,
    ctx: Optional[PanelContext] = None,
    cover_fraction: float = 0.9,
    recovery_threshold: float = 0.5,
    per_gene: bool = True,
) -> dict:
    rf = fam if isinstance(fam, ReconciledFamily) else reconcile(fam, st)
    F = rf.family
    ctx = ctx or panel_context([F])
    n = len(rf.genes)
    Mref = reference_merge_matrix(rf)
    codes = hog_codes(rf, pred)
    Mpred = predicted_merge_matrix(rf, codes)

    iu, ju = np.triu_indices(n, 1)
    sel = Mref[iu, ju] != NOT_ASSERTED
    i, j = iu[sel], ju[sel]
    ref, prd = Mref[i, j], Mpred[i, j]
    s0 = st.lca_idx[rf.species_idx[i], rf.species_idx[j]]
    correct = ref == prd
    res: dict = {
        "id": F.id, "name": F.name, "panel": F.panel, "phenomena": F.phenomena,
        "family_level": st.names[rf.family_level], "n_members": len(F.member_ids),
        "n_tree_genes": n, "n_pairs_total": n * (n - 1) // 2, "n_asserted_pairs": int(len(i)),
    }
    warn = []
    gsp = F.gene_species
    mism = [g for g in rf.genes if g in pred.gene_species and pred.gene_species[g]
            and pred.gene_species[g] != gsp[g]]
    if mism:
        warn.append(f"{len(mism)} gene(s) have a different species in the prediction, e.g. {mism[:3]}")
    absent = [g for g in F.member_ids if g not in pred.gene_species and
              not any(g in m for m in pred.levels.values())]
    res["n_members_absent_from_prediction"] = len(absent)

    # 1. pair merge-level accuracy
    err, pinf, minf = _signed_errors(st, ref, prd)
    npairs = len(i)
    res["pair_merge"] = {
        "n_pairs": npairs,
        "accuracy": float(correct.mean()) if npairs else None,
        "error_histogram": _histogram(err, pinf, minf, int(((ref == INF) & (prd == INF)).sum())),
        "frac_over_split": float(((err > 0).sum() + pinf) / npairs) if npairs else None,
        "frac_over_merged": float(((err < 0).sum() + minf) / npairs) if npairs else None,
        "mean_abs_error_finite": float(np.abs(err).mean()) if len(err) else None,
    }

    # 2. per-level HOG pair precision/recall/F1 (levels within the family's root level)
    per_level = []
    tot = np.zeros(3, dtype=np.int64)
    for L in sorted(codes, key=lambda x: st.depth[x]):
        if not st.anc[rf.family_level, L]:
            continue
        mask = st.anc[L, s0]
        if not mask.any():
            continue
        r_l = ref[mask]
        ref_same = (r_l >= 0) & st.anc[L, np.where(r_l >= 0, r_l, 0)]
        c = codes[L]
        ci, cj = c[i[mask]], c[j[mask]]
        pred_same = (ci == cj) & (ci >= 0)
        tp = int((ref_same & pred_same).sum())
        fp = int((~ref_same & pred_same).sum())
        fn = int((ref_same & ~pred_same).sum())
        tot += (tp, fp, fn)
        d = {"level": st.names[L], "depth": int(st.depth[L]), "n_pairs": int(mask.sum()),
             "n_ref_same": int(ref_same.sum()), "n_pred_same": int(pred_same.sum())}
        d.update(_prf(tp, fp, fn))
        per_level.append(d)
    f1s = [d["f1"] for d in per_level if d["f1"] is not None]
    res["per_level"] = per_level
    res["per_level_summary"] = {"mean_f1": float(np.mean(f1s)) if f1s else None,
                                "n_levels": len(per_level), "pooled": _prf(*tot)}

    # 3. orthology (cross-species pairs only; within a family "not ortholog" = paralog)
    cross = rf.species_idx[i] != rf.species_idx[j]
    ro, po = (ref == s0)[cross], (prd == s0)[cross]
    res["orthology"] = {
        "n_pairs": int(cross.sum()), "n_ref_orthologs": int(ro.sum()), "n_ref_paralogs": int((~ro).sum()),
        "ortholog": _prf(int((ro & po).sum()), int((~ro & po).sum()), int((ro & ~po).sum())),
        "paralog": _prf(int((~ro & ~po).sum()), int((ro & ~po).sum()), int((~ro & po).sum())),
        "accuracy": float((ro == po).mean()) if cross.any() else None,
    }

    # 4. event recovery for every asserted duplication
    events = []
    for ev in rf.events:
        if ev.ev != "D":
            continue
        vals = []
        ch = ev.child_leaves
        for a in range(len(ch)):
            for b in range(a + 1, len(ch)):
                vals.append(Mpred[np.ix_(ch[a], ch[b])].ravel())
        pv = np.concatenate(vals)
        rv = np.full(len(pv), ev.merge, dtype=np.int64)
        e_err, e_pinf, e_minf = _signed_errors(st, rv, pv)
        frac = float((pv == ev.merge).mean())
        modal = Counter(pv.tolist()).most_common(1)[0][0]
        if ev.merge >= 0 and modal >= 0:
            mis: Optional[int | str] = int(st.depth[ev.merge] - st.depth[modal])
        elif ev.merge == modal:
            mis = 0
        else:
            mis = "+inf" if modal == INF else "-inf"
        events.append({
            "event": ev.id, "level": st.names[ev.level], "level_explicit": ev.level_explicit,
            "ref_merge": _level_name(st, ev.merge), "tags": ev.tags, "src": ev.src,
            "n_copies": len(ch), "n_pairs": int(len(pv)), "frac_correct": frac,
            "recovered": frac >= recovery_threshold,
            "frac_over_split": float(((e_err > 0).sum() + e_pinf) / len(pv)),
            "frac_over_merged": float(((e_err < 0).sum() + e_minf) / len(pv)),
            "median_signed_error_finite": float(np.median(e_err)) if len(e_err) else None,
            "modal_pred_merge": _level_name(st, modal), "misplacement": mis,
        })
    res["events"] = events
    res["events_summary"] = {
        "n_events": len(events),
        "recovery_rate": float(np.mean([e["recovered"] for e in events])) if events else None,
        "mean_frac_correct": float(np.mean([e["frac_correct"] for e in events])) if events else None,
    }

    # 5. family integrity at the family's root level, with the reference's own expectation
    integ = family_integrity(F, rf.family_level, st, pred, ctx, cover_fraction)
    ref_lab = reference_level_codes(rf, [rf.family_level])[rf.family_level] if n else np.array([])
    ref_pred = Prediction()
    for k, g in enumerate(rf.genes):
        if ref_lab[k] >= 0:
            ref_pred.add(st.names[rf.family_level], f"R{ref_lab[k]}", g)
    exp = family_integrity(F, rf.family_level, st, ref_pred, PanelContext({}, set()), cover_fraction)
    integ["ref_n_hogs_cover"] = exp["n_hogs_cover"]
    integ["ref_n_hogs_all"] = exp["n_hogs_all"]
    integ["ref_coverage"] = exp["coverage"]
    integ["split_excess"] = integ["n_hogs_cover"] - exp["n_hogs_cover"]
    res["integrity"] = integ

    # 6. gene-level / species-level correctness
    n_pairs_g = np.bincount(i, minlength=n) + np.bincount(j, minlength=n)
    n_corr_g = np.bincount(i, weights=correct, minlength=n) + np.bincount(j, weights=correct, minlength=n)
    with np.errstate(invalid="ignore", divide="ignore"):
        frac_g = np.where(n_pairs_g > 0, n_corr_g / np.maximum(n_pairs_g, 1), np.nan)
    if per_gene:
        res["genes"] = [{"gene": g, "species": gsp[g], "n_pairs": int(n_pairs_g[k]),
                         "frac_correct": None if np.isnan(frac_g[k]) else float(frac_g[k])}
                        for k, g in enumerate(rf.genes)]
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
    res["warnings"] = warn
    return res


def family_integrity(F: ReferenceFamily, family_level: int, st: SpeciesTree, pred: Prediction,
                     ctx: PanelContext, cover_fraction: float = 0.9) -> dict:
    """Split / contamination / coverage of a family's members among predicted HOGs.

    Each member is assigned its *top* predicted HOG at or below the family level (the HOG at
    the highest level ≤ family level on its species' path that contains it); members with
    no HOG count as singletons. ``n_hogs_cover`` = fewest such HOGs covering ≥ ``cover_fraction``
    of members; ``contamination`` = fraction of all genes in those HOGs that are members of
    *other* reference families or declared decoys.
    """
    members = F.member_ids
    gsp = F.gene_species
    keys: dict[str, tuple] = {}
    for g in members:
        top = None
        for lv in st.ancestors(gsp[g], include_self=False):
            if not st.anc[family_level, st.index[lv]]:
                break  # above the family level (or member outside the family clade)
            h = pred.hog(lv, g)
            if h is not None:
                top = (lv, h)
        keys[g] = top if top is not None else ("<singleton>", g)
    cnt = Counter(keys.values())
    ranked = sorted(cnt.items(), key=lambda kv: (-kv[1], str(kv[0])))
    n = len(members)
    need = cover_fraction * n - 1e-9
    cum, k90 = 0, 0
    for _, c in ranked:
        if cum >= need:
            break
        cum += c
        k90 += 1
    cover_keys = [k for k, _ in ranked[:k90]]
    total = contam = own = 0
    decoys = ctx.decoys | set(F.decoys)
    for lv, h in cover_keys:
        if lv == "<singleton>":
            total += 1
            own += 1
            continue
        total += pred.hog_size(lv, h)
        for g in pred.hog_members(lv, h):
            fam = ctx.gene_family.get(g)
            if fam == F.id:
                own += 1
            elif fam is not None or g in decoys:
                contam += 1
    covered = sum(1 for g in members if keys[g][0] != "<singleton>"
                  and pred.hog_size(*keys[g]) >= 2)
    return {
        "level": st.names[family_level], "n_members": n,
        "n_hogs_cover": k90, "cover_fraction": cover_fraction, "n_hogs_all": len(cnt),
        "largest_hog_frac": ranked[0][1] / n if n else None,
        "n_genes_in_cover_hogs": total, "n_contaminants": contam,
        "contamination": contam / total if total else None,
        "frac_unlabelled": (total - own - contam) / total if total else None,
        "coverage": covered / n if n else None,
    }


# --------------------------------------------------------------------------- aggregation
SCALARS = {
    "merge_accuracy": ("pair_merge", "accuracy"),
    "frac_over_split": ("pair_merge", "frac_over_split"),
    "frac_over_merged": ("pair_merge", "frac_over_merged"),
    "level_mean_f1": ("per_level_summary", "mean_f1"),
    "level_pooled_f1": ("per_level_summary", "pooled", "f1"),
    "ortholog_precision": ("orthology", "ortholog", "precision"),
    "ortholog_recall": ("orthology", "ortholog", "recall"),
    "ortholog_f1": ("orthology", "ortholog", "f1"),
    "paralog_f1": ("orthology", "paralog", "f1"),
    "event_recovery_rate": ("events_summary", "recovery_rate"),
    "event_mean_frac_correct": ("events_summary", "mean_frac_correct"),
    "split_n_hogs_cover": ("integrity", "n_hogs_cover"),
    "split_excess": ("integrity", "split_excess"),
    "split_n_hogs_all": ("integrity", "n_hogs_all"),
    "contamination": ("integrity", "contamination"),
    "coverage": ("integrity", "coverage"),
    "gene_mean_frac_correct": ("genes_summary", "mean_frac_correct"),
    "gene_frac_perfect": ("genes_summary", "frac_perfect"),
}


def family_scalars(res: dict) -> dict[str, Optional[float]]:
    out: dict[str, Optional[float]] = {}
    for name, path in SCALARS.items():
        v = res
        for k in path:
            v = v.get(k) if isinstance(v, dict) else None
        out[name] = None if v is None else float(v)
    return out


def bootstrap_ci(values: Sequence[float], n_boot: int = 2000, seed: int = SEEDS[0],
                 alpha: float = 0.05) -> tuple[Optional[float], Optional[float]]:
    """Percentile bootstrap CI of the mean, resampling families with replacement."""
    v = np.asarray([x for x in values if x is not None and not (isinstance(x, float) and math.isnan(x))],
                   dtype=float)
    if len(v) == 0:
        return None, None
    if len(v) == 1:
        return float(v[0]), float(v[0])
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, len(v), size=(n_boot, len(v)))
    means = v[idx].mean(axis=1)
    return float(np.quantile(means, alpha / 2)), float(np.quantile(means, 1 - alpha / 2))


def macro_average(results: Sequence[dict], n_boot: int = 2000, seed: int = SEEDS[0]) -> dict:
    """Mean over families (each family weight 1) of every scalar, with bootstrap CI."""
    rows = [family_scalars(r) for r in results]
    out = {}
    for name in SCALARS:
        vals = [r[name] for r in rows if r[name] is not None]
        lo, hi = bootstrap_ci(vals, n_boot=n_boot, seed=seed)
        out[name] = {"mean": float(np.mean(vals)) if vals else None, "ci_low": lo, "ci_high": hi,
                     "n_families": len(vals)}
    return out


def aggregate(results: Sequence[dict], n_boot: int = 2000, seed: int = SEEDS[0]) -> dict:
    """Macro averages overall and per phenomenon tag (family tags), plus event-level
    recovery per event tag (pooled over events of all families)."""
    out = {"n_families": len(results), "macro": macro_average(results, n_boot, seed), "per_tag": {},
           "events_by_tag": {}}
    tags = sorted({t for r in results for t in r.get("phenomena", [])})
    for t in tags:
        sub = [r for r in results if t in r.get("phenomena", [])]
        out["per_tag"][t] = macro_average(sub, n_boot, seed)
    ev_by_tag: dict[str, list[dict]] = defaultdict(list)
    for r in results:
        for e in r.get("events", []):
            for t in (e["tags"] or ["<untagged>"]):
                ev_by_tag[t].append(e)
    for t, evs in sorted(ev_by_tag.items()):
        out["events_by_tag"][t] = {
            "n_events": len(evs),
            "recovery_rate": float(np.mean([e["recovered"] for e in evs])),
            "mean_frac_correct": float(np.mean([e["frac_correct"] for e in evs])),
        }
    return out
