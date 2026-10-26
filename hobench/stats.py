"""Aggregation and inference over families (SYNTHESIS B4, B9, B10).

* **Items.** A family without ``collection`` is one item. All families sharing a
  ``collection`` id are *units* of one item (two-stage aggregation, B10): the collection's
  value is the mean over its units (CI by resampling units), and it counts as ONE item in the
  combined macro mean. Ordinary families form the main ``families`` track; collections are
  reported in their own track.
* **Bootstrap (B9).** Percentile bootstrap of the mean over item IDs, stratified by panel.
  Indices are a deterministic function of (sorted item IDs, strata, seed), so every method
  evaluated on the same items gets the same resamples (paired). No CI when n < ``MIN_N_CI``
  (= 5): the range is reported instead. n is always reported.
* **Paired comparison (B9).** Mean difference with paired bootstrap CI, Wilcoxon signed-rank
  p-value, Holm correction over a set of comparisons.
* **Agreement (B4).** Adjusted mutual information (AMI, arithmetic normalisation as in
  scikit-learn) of two partitions of the reference members per level.
"""
from __future__ import annotations

import math
from collections import defaultdict
from typing import Iterable, Optional, Sequence

import numpy as np
from scipy.special import gammaln

from .metrics import PRIMARY, SCALARS, THETAS, family_scalars
from .seeds import SEEDS

MIN_N_CI = 5


# ------------------------------------------------------------------------- bootstrap
def bootstrap_indices(ids: Sequence[str], n_boot: int = 2000, seed: int = SEEDS[0],
                      strata: Optional[dict[str, str]] = None) -> tuple[list[str], np.ndarray]:
    """Return (sorted ids, n_boot × n index matrix): resampling with replacement within each
    stratum (keeping stratum sizes). Deterministic in (ids, strata, seed)."""
    ids = sorted(ids)
    rng = np.random.default_rng(seed)
    groups: dict[str, list[int]] = defaultdict(list)
    for k, x in enumerate(ids):
        groups[(strata or {}).get(x, "")].append(k)
    idx = np.empty((n_boot, len(ids)), dtype=np.int64)
    col = 0
    for g in sorted(groups):
        members = np.array(groups[g])
        draw = rng.integers(0, len(members), size=(n_boot, len(members)))
        idx[:, col:col + len(members)] = members[draw]
        col += len(members)
    return ids, idx


def mean_ci(values: dict[str, Optional[float]], n_boot: int = 2000, seed: int = SEEDS[0],
            strata: Optional[dict[str, str]] = None, alpha: float = 0.05, min_n: int = MIN_N_CI) -> dict:
    """Mean over items with a stratified percentile bootstrap CI (None when n < ``min_n``)."""
    v = {k: x for k, x in values.items() if x is not None and not (isinstance(x, float) and math.isnan(x))}
    n = len(v)
    out = {"mean": None, "ci_low": None, "ci_high": None, "n": n, "min": None, "max": None}
    if not n:
        return out
    arr = np.array(list(v.values()), dtype=float)
    out.update(mean=float(arr.mean()), min=float(arr.min()), max=float(arr.max()))
    if n >= min_n:
        ids, idx = bootstrap_indices(list(v), n_boot, seed, strata)
        a = np.array([v[k] for k in ids])
        means = a[idx].mean(axis=1)
        out["ci_low"], out["ci_high"] = float(np.quantile(means, alpha / 2)), float(np.quantile(means, 1 - alpha / 2))
    return out


def bootstrap_ci(values: Sequence[float], n_boot: int = 2000, seed: int = SEEDS[0], alpha: float = 0.05,
                 min_n: int = MIN_N_CI) -> tuple[Optional[float], Optional[float]]:
    """Convenience wrapper: CI of the mean of a plain list (None, None when n < ``min_n``)."""
    d = mean_ci({str(k): x for k, x in enumerate(values)}, n_boot, seed, None, alpha, min_n)
    return d["ci_low"], d["ci_high"]


# ------------------------------------------------------------------------ items
def items(results: Sequence[dict]) -> dict[str, dict]:
    """Collapse family results into items: ``{item_id: {"panel", "collection", "tags",
    "scalars", "n_units", "events"}}``. Collection items average their units' scalars."""
    out: dict[str, dict] = {}
    by_coll: dict[str, list[dict]] = defaultdict(list)
    for r in results:
        if r.get("collection"):
            by_coll[r["collection"]].append(r)
        else:
            out[r["id"]] = {"panel": r.get("panel"), "collection": None, "tags": list(r.get("phenomena", [])),
                            "scalars": family_scalars(r), "n_units": 1, "events": r.get("events", [])}
    for c, units in by_coll.items():
        sc = [family_scalars(u) for u in units]
        keys = sc[0].keys() if sc else []
        mean = {}
        for k in keys:
            xs = [s[k] for s in sc if s.get(k) is not None]
            mean[k] = float(np.mean(xs)) if xs else None
        panels = {u.get("panel") for u in units}
        out[f"collection:{c}"] = {"panel": ",".join(sorted(map(str, panels))), "collection": c,
                                  "tags": sorted({t for u in units for t in u.get("phenomena", [])}),
                                  "scalars": mean, "n_units": len(units),
                                  "events": [e for u in units for e in u.get("events", [])]}
    return out


def _track(its: dict[str, dict], n_boot: int, seed: int, metrics: Iterable[str]) -> dict:
    strata = {k: str(v["panel"]) for k, v in its.items()}
    return {m: mean_ci({k: v["scalars"].get(m) for k, v in its.items()}, n_boot, seed, strata)
            for m in metrics}


def aggregate(results: Sequence[dict], n_boot: int = 2000, seed: int = SEEDS[0]) -> dict:
    """Two-stage aggregation (B10) with stratified bootstrap CIs (B9).

    ``families``: ordinary families only. ``collections``: per collection, mean over units
    with a unit-bootstrap CI. ``combined``: every collection counts as one item.
    ``per_tag``: combined items carrying the tag. ``events_by_tag``: each item's mean event
    recovery for events with that tag, then averaged over items (family/collection weight 1).
    """
    metrics = list(SCALARS) + sorted({k for r in results for k in family_scalars(r) if k.startswith("skill_")})
    its = items(results)
    fam_items = {k: v for k, v in its.items() if v["collection"] is None}
    out = {"n_families": len(results), "n_items": len(its), "primary": PRIMARY,
           "families": _track(fam_items, n_boot, seed, metrics), "collections": {}, "combined": _track(its, n_boot, seed, metrics),
           "per_tag": {}, "events_by_tag": {}}
    by_coll: dict[str, list[dict]] = defaultdict(list)
    for r in results:
        if r.get("collection"):
            by_coll[r["collection"]].append(r)
    for c, units in sorted(by_coll.items()):
        uits = {u["id"]: {"panel": u.get("panel"), "scalars": family_scalars(u)} for u in units}
        out["collections"][c] = {"n_units": len(units), **_track(uits, n_boot, seed, metrics)}
    for t in sorted({t for v in its.values() for t in v["tags"]}):
        sub = {k: v for k, v in its.items() if t in v["tags"]}
        out["per_tag"][t] = _track(sub, n_boot, seed, metrics)
    ev_items: dict[str, dict[str, list]] = defaultdict(lambda: defaultdict(list))
    for k, v in its.items():
        for e in v["events"]:
            for t in (e["tags"] or ["<untagged>"]):
                ev_items[t][k].append(e)
    for t, per_item in sorted(ev_items.items()):
        d = {"n_events": sum(len(x) for x in per_item.values()), "n_items": len(per_item)}
        for th in THETAS:
            d[f"recovery_rate_{th}"] = float(np.mean([np.mean([e[f"recovered_{th}"] for e in evs])
                                                      for evs in per_item.values()]))
        d["mean_frac_correct"] = float(np.mean([np.mean([e["frac_correct"] for e in evs]) for evs in per_item.values()]))
        out["events_by_tag"][t] = d
    return out


# ----------------------------------------------------------------- paired comparisons
def compare_methods(results_a: Sequence[dict], results_b: Sequence[dict], metric: str,
                    n_boot: int = 2000, seed: int = SEEDS[0], alpha: float = 0.05, min_n: int = MIN_N_CI) -> dict:
    """Paired comparison of two methods on the items both have a value for (B9): mean
    difference a − b, stratified paired bootstrap CI (same resampled item IDs for both), and a
    two-sided Wilcoxon signed-rank p-value."""
    from scipy.stats import wilcoxon

    ia, ib = items(results_a), items(results_b)
    common = sorted(k for k in set(ia) & set(ib)
                    if ia[k]["scalars"].get(metric) is not None and ib[k]["scalars"].get(metric) is not None)
    diffs = {k: ia[k]["scalars"][metric] - ib[k]["scalars"][metric] for k in common}
    strata = {k: str(ia[k]["panel"]) for k in common}
    d = mean_ci(diffs, n_boot, seed, strata, alpha, min_n)
    out = {"metric": metric, "n": len(common), "mean_diff": d["mean"], "ci_low": d["ci_low"],
           "ci_high": d["ci_high"], "min_diff": d["min"], "max_diff": d["max"],
           "mean_a": float(np.mean([ia[k]["scalars"][metric] for k in common])) if common else None,
           "mean_b": float(np.mean([ib[k]["scalars"][metric] for k in common])) if common else None,
           "wilcoxon_p": None}
    nz = [x for x in diffs.values() if x != 0]
    if common:
        out["wilcoxon_p"] = 1.0 if not nz else float(wilcoxon(list(diffs.values()), zero_method="wilcox").pvalue)
    return out


def holm(pvalues: dict[str, Optional[float]]) -> dict[str, Optional[float]]:
    """Holm–Bonferroni adjusted p-values (None passes through)."""
    items_ = sorted(((p, k) for k, p in pvalues.items() if p is not None))
    m = len(items_)
    adj: dict[str, Optional[float]] = {k: None for k, p in pvalues.items() if p is None}
    running = 0.0
    for r, (p, k) in enumerate(items_):
        running = max(running, min(1.0, (m - r) * p))
        adj[k] = running
    return adj


# ---------------------------------------------------------------------------- AMI
def _entropy(counts: np.ndarray) -> float:
    n = counts.sum()
    p = counts[counts > 0] / n
    return float(-(p * np.log(p)).sum())


def ami(labels_a: Sequence, labels_b: Sequence) -> float:
    """Adjusted mutual information (Vinh et al. 2010; arithmetic-mean normalisation, as
    scikit-learn's default). Returns 1.0 for identical partitions, ~0 for independent ones."""
    a = np.unique(np.asarray(labels_a, dtype=object).astype(str), return_inverse=True)[1]
    b = np.unique(np.asarray(labels_b, dtype=object).astype(str), return_inverse=True)[1]
    n = len(a)
    if n != len(b):
        raise ValueError("label vectors differ in length")
    R, C = a.max() + 1 if n else 0, b.max() + 1 if n else 0
    if n == 0 or (R == C == 1) or (R == C == n):
        return 1.0
    cont = np.zeros((R, C), dtype=np.int64)
    np.add.at(cont, (a, b), 1)
    ai, bj = cont.sum(1), cont.sum(0)
    nz = cont > 0
    mi = float((cont[nz] / n * np.log(n * cont[nz] / np.outer(ai, bj)[nz])).sum())
    # expected MI
    emi = 0.0
    lg_n = gammaln(n + 1)
    for x in ai:
        for y in bj:
            lo, hi = max(1, x + y - n), min(x, y)
            if lo > hi:
                continue
            nij = np.arange(lo, hi + 1)
            term = nij / n * np.log(n * nij / (x * y))
            lp = (gammaln(x + 1) + gammaln(y + 1) + gammaln(n - x + 1) + gammaln(n - y + 1) - lg_n
                  - gammaln(nij + 1) - gammaln(x - nij + 1) - gammaln(y - nij + 1) - gammaln(n - x - y + nij + 1))
            emi += float((term * np.exp(lp)).sum())
    h = 0.5 * (_entropy(ai) + _entropy(bj))
    denom = h - emi
    if abs(denom) < 1e-15:
        return 1.0 if abs(mi - emi) < 1e-15 else 0.0
    return float((mi - emi) / denom)


def agreement(fam, st, pred_a, pred_b) -> dict:
    """Per-level AMI of two predictions on a family's members (genes with no HOG at a level
    are singletons); levels at or below the family level with ≥2 members below them."""
    from .reference import reconcile
    rf = reconcile(fam, st, check=False)
    out = {}
    for L in sorted({p for s in set(rf.species_idx.tolist()) for p in st.ancestors(st.names[s], include_self=False)},
                    key=lambda x: st.depth_of(x)):
        Li = st.idx(L)
        if not st.anc[rf.family_level, Li]:
            continue
        genes = [g for k, g in enumerate(rf.genes) if st.anc[Li, rf.species_idx[k]]]
        if len(genes) < 2:
            continue
        la = [pred_a.hog(L, g) or f"_single_{g}" for g in genes]
        lb = [pred_b.hog(L, g) or f"_single_{g}" for g in genes]
        out[L] = ami(la, lb)
    return {"per_level": out, "mean": float(np.mean(list(out.values()))) if out else None}
