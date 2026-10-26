"""Control predictions with known behaviour.

* :func:`reference_prediction` — **positive control**: HOGs derived from the reference trees
  (at each level, connected components of "merge level ⊆ L"); all pair metrics are perfect.
* :func:`all_singletons` — **negative**: no HOGs at all (every pair merges at ∞).
* :func:`one_big_hog` — **negative**: at every level, one HOG with every gene below it
  (every pair merges at its species LCA; same-species pairs at the species' parent).
* :func:`random_relabel` — **negative**: a prediction with gene labels permuted uniformly
  *within species* (keeps HOG sizes and species composition, destroys gene identity).
* :func:`perturbed_reference_prediction` — **dose-response**: positive control built from
  reference trees after ``k`` random leaf swaps between genes of the same species (anywhere
  in the panel, so swaps can move genes between families).

Every random control takes an explicit seed (use values from :data:`hobench.seeds.SEEDS`).
"""
from __future__ import annotations

from collections import defaultdict
from typing import Optional, Sequence

import numpy as np
import pandas as pd

from .metrics import NOT_ASSERTED, reference_level_codes, reference_matrices
from .newick import copy_tree
from .predictions import Prediction
from .reference import ReferenceFamily, reconcile
from .speciestree import SpeciesTree


def _all_genes(families: Sequence[ReferenceFamily]) -> dict[str, str]:
    out: dict[str, str] = {}
    for f in families:
        out.update(f.gene_species)
    return out


def _count_species(families: Sequence[ReferenceFamily]) -> set[str]:
    return {s for f in families for s in f.complete_species}


def reference_prediction(families: Sequence[ReferenceFamily], st: SpeciesTree,
                         check: bool = True, label: str = "positive") -> Prediction:
    """Positive control: reference-derived HOGs at every level (``metrics.representative_matrix``:
    one consistent merge level inside every reference interval; U pairs grouped at the oldest
    level implied inside the U node). Members not placed in the tree stay singletons.
    ``meta['inconsistent_pairs']`` counts asserted pairs whose reference status at a level
    (surely together / surely apart) disagrees with the components (0 for valid references)."""
    pred = Prediction(meta={"control": label}, count_species=_count_species(families))
    bad = 0
    for fam in families:
        rf = reconcile(fam, st, check=check)
        codes = reference_level_codes(rf)
        R = reference_matrices(rf)
        gsp = fam.gene_species
        for L, lab in codes.items():
            lv = st.names[L]
            for k, g in enumerate(rf.genes):
                if lab[k] >= 0:
                    pred.add(lv, f"{fam.id}|{lv}|{lab[k]}", g, species=gsp[g])
            idx = np.flatnonzero(lab >= 0)
            if len(idx) > 1:
                lo = R.lo[np.ix_(idx, idx)].astype(np.int64)
                hi = R.hi[np.ix_(idx, idx)].astype(np.int64)
                same = lab[idx][:, None] == lab[idx][None, :]
                asserted = lo != NOT_ASSERTED
                surely_same = asserted & (hi >= 0) & st.anc[L, np.where(hi >= 0, hi, 0)]
                surely_apart = asserted & ((lo < 0) | ~st.anc[L, np.where(lo >= 0, lo, 0)])
                np.fill_diagonal(surely_apart, False)
                bad += int(((surely_apart & same) | (surely_same & ~same)).sum() // 2)
        pred.gene_species.update({g: gsp[g] for g in rf.genes})
    pred.meta["inconsistent_pairs"] = bad
    return pred


def all_singletons(families: Sequence[ReferenceFamily] = ()) -> Prediction:
    pred = Prediction(meta={"control": "all_singletons"}, count_species=_count_species(families))
    pred.gene_species.update(_all_genes(families))
    return pred


def one_big_hog(families: Sequence[ReferenceFamily], st: SpeciesTree) -> Prediction:
    pred = Prediction(meta={"control": "one_big_hog"}, count_species=_count_species(families))
    genes = _all_genes(families)
    for lv in st.levels:
        clade = st.clade_species(lv)
        for g, s in genes.items():
            if s in clade:
                pred.add(lv, f"ALL|{lv}", g, species=s)
    pred.gene_species.update(genes)
    return pred


def random_relabel(pred: Prediction, seed: int, gene_species: Optional[dict[str, str]] = None) -> Prediction:
    """Permute gene labels uniformly within each species (genes of ``gene_species``, default
    the prediction's own kept genes). A gene absent from every HOG can swap with one present."""
    gsp = dict(gene_species or pred.gene_species)
    rng = np.random.default_rng(seed)
    by_sp: dict[str, list[str]] = defaultdict(list)
    for g, s in sorted(gsp.items()):
        by_sp[s].append(g)
    perm: dict[str, str] = {}
    for s in sorted(by_sp):
        gs = by_sp[s]
        p = rng.permutation(len(gs))
        for a, b in zip(gs, p):
            perm[a] = gs[b]   # gene `a` takes the HOG memberships of gene `gs[b]`
    inv = {v: k for k, v in perm.items()}
    out = Prediction(meta={"control": "random_relabel", "seed": seed}, count_species=set(pred.count_species))
    for lv, m in pred.levels.items():
        for g, h in m.items():
            out.add(lv, h, inv.get(g, g))
    # sizes and per-species counts are unchanged (relabelling stays within species)
    out.sizes = {lv: dict(d) for lv, d in pred.sizes.items()}
    out.hog_species = {lv: {h: dict(c) for h, c in d.items()} for lv, d in pred.hog_species.items()}
    out.gene_species = dict(gsp)
    out.gene_species.update({g: s for g, s in pred.gene_species.items() if g not in gsp})
    return out


def swap_leaves(families: Sequence[ReferenceFamily], k: int, seed: int) -> list[ReferenceFamily]:
    """Return copies of ``families`` after ``k`` random swaps of two same-species genes
    (the gene's species is unchanged, so trees stay species-consistent)."""
    rng = np.random.default_rng(seed)
    fams = [ReferenceFamily(id=f.id, tree=copy_tree(f.tree), members=f.members.copy(), meta=dict(f.meta),
                            path=f.path) for f in families]
    leaf_of: dict[str, tuple[int, object]] = {}
    for fi, f in enumerate(fams):
        for lf in f.tree.leaves():
            leaf_of[lf.name] = (fi, lf)
    gsp = {}
    for f in fams:
        gsp.update(f.gene_species)
    by_sp: dict[str, list[str]] = defaultdict(list)
    for g in sorted(leaf_of):
        by_sp[gsp[g]].append(g)
    pool = [g for s in sorted(by_sp) if len(by_sp[s]) > 1 for g in by_sp[s]]
    if not pool:
        return fams
    for _ in range(k):
        g = pool[rng.integers(len(pool))]
        others = [x for x in by_sp[gsp[g]] if x != g]
        h = others[rng.integers(len(others))]
        (fg, lg), (fh, lh) = leaf_of[g], leaf_of[h]
        lg.name, lh.name = h, g
        leaf_of[g], leaf_of[h] = (fh, lh), (fg, lg)
        if fg != fh:  # move membership rows between families
            mg, mh = fams[fg].members, fams[fh].members
            rg = mg[mg["gene_id"] == g]
            rh = mh[mh["gene_id"] == h]
            fams[fg].members = pd.concat([mg[mg["gene_id"] != g], rh], ignore_index=True)
            fams[fh].members = pd.concat([mh[mh["gene_id"] != h], rg], ignore_index=True)
    return fams


def perturbed_reference_prediction(families: Sequence[ReferenceFamily], st: SpeciesTree,
                                   k: int, seed: int) -> Prediction:
    """Dose-response control: positive control computed on reference trees with ``k`` random
    same-species leaf swaps."""
    pert = swap_leaves(families, k, seed)
    pred = reference_prediction(pert, st, check=False, label="perturbed_reference")
    pred.meta.update({"k_swaps": k, "seed": seed})
    return pred
