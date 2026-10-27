#!/usr/bin/env python
"""Pair-level extraction (re-uses the frozen scorer's functions; no metric is changed).
For every tool run of panel P/mode (same loading as scripts/evaluate_panel.py; FastOMA 'declared' placement only):
  nodes_<P>_<mode>.tsv.gz : per (method, family, reference LCA node): ev, tags, node level/depth, merge interval,
                            n asserted pairs, n correct (merge level in interval), n over-split, n over-merged.
  ami_<P>_<mode>.tsv      : AMI between OrthoFinder-default and FastOMA-declared HOG partitions of the reference tree genes,
                            per family and level (genes without a HOG at that level = singletons).
Output dir: results/analysis/cache/.  Usage: python scripts/analysis/extract_nodes.py --panel Y --mode full"""
import argparse, csv, gzip, sys
from pathlib import Path
import numpy as np
from sklearn.metrics import adjusted_mutual_info_score
from common import ROOT, OUT
import evaluate_panel as ep
from hobench.ids import make_normalizer
from hobench.metrics import (reference_matrices, hog_codes, predicted_merge_matrix, interval_correct,
                             interval_error, NOT_ASSERTED)
from hobench.reference import load_families, reconcile
from hobench.speciestree import SpeciesTree


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--panel", required=True); ap.add_argument("--mode", default="full")
    ap.add_argument("--cache", default="/tmp/hobench_eval_cache")
    a = ap.parse_args(); P, mode = a.panel, a.mode
    st = SpeciesTree.from_newick(str(ROOT / "panels" / P / "species_tree.nwk"))
    fams = load_families(str(ROOT / "families"), panel=P)
    if mode == "subset":
        sub = {r["id"] for r in csv.DictReader(open(ROOT / "panels" / P / "subset" / "members.tsv"), delimiter="\t")}
        fams = [f for f in fams if set(f.member_ids) <= sub]
    keep, cs = set(), set()
    for f in fams:
        keep |= set(f.member_ids) | set(f.decoys); cs |= set(f.complete_species)
    kw = {"normalize": make_normalizer("auto"), "keep_genes": keep, "count_species": cs}
    rfs = [reconcile(f, st, check=False) for f in fams]
    R = {rf.family.id: reference_matrices(rf) for rf in rfs}
    out = OUT / "cache"; out.mkdir(exist_ok=True)
    codes_by = {}
    with gzip.open(out / f"nodes_{P}_{mode}.tsv.gz", "wt") as fh:
        fh.write("method\tvariant\tfamily\tcollection\tnode\tev\ttags\tnode_level\tnode_depth\tmerge_lo\tmerge_hi\t"
                 "n_pairs\tn_correct\tn_over_split\tn_over_merged\n")
        rundir = ROOT / "runs" / P / mode
        for run in sorted(p for p in rundir.iterdir() if (p / "manifest.json").exists()):
            if not run.name.startswith(("orthofinder", "fastoma")):
                continue
            for tool, variant, pred, *_ in ep.load_run(run, st, kw, Path(a.cache), ["declared"]):
                print("scoring", tool, variant, file=sys.stderr)
                for rf in rfs:
                    Rm = R[rf.family.id]; n = len(rf.genes)
                    codes = hog_codes(rf, pred)
                    if (tool, variant) in (("orthofinder", "default"), ("fastoma", "default_declared")):
                        codes_by[(tool, rf.family.id)] = codes
                    M = predicted_merge_matrix(rf, codes)
                    iu, ju = np.triu_indices(n, 1)
                    sel = Rm.lo[iu, ju] != NOT_ASSERTED
                    i, j = iu[sel], ju[sel]
                    lo, hi, p, node = Rm.lo[i, j].astype(np.int64), Rm.hi[i, j].astype(np.int64), M[i, j].astype(np.int64), Rm.node[i, j]
                    cor = interval_correct(st, p, lo, hi); err = interval_error(st, p, lo, hi)
                    for k in np.unique(node):
                        m = node == k; ev = rf.events[k]
                        lvl = ev.mapping
                        fh.write(f"{tool}\t{variant}\t{rf.family.id}\t{rf.family.collection or ''}\t{ev.id}\t{ev.ev}\t"
                                 f"{','.join(ev.tags or [])}\t{st.names[lvl]}\t{int(st.depth[lvl])}\t"
                                 f"{'INF' if ev.merge_lo < 0 else st.names[ev.merge_lo]}\t{'INF' if ev.merge_hi < 0 else st.names[ev.merge_hi]}\t"
                                 f"{int(m.sum())}\t{int(cor[m].sum())}\t{int((err[m] > 0).sum())}\t{int((err[m] < 0).sum())}\n")
    with open(out / f"ami_{P}_{mode}.tsv", "w") as fh:
        fh.write("family\tcollection\tlevel\tdepth\tis_family_level\tn_genes\tami\n")
        for rf in rfs:
            a_, b_ = codes_by.get(("orthofinder", rf.family.id)), codes_by.get(("fastoma", rf.family.id))
            if a_ is None or b_ is None:
                continue
            for L in sorted(a_, key=lambda x: st.depth[x]):
                below = np.flatnonzero(st.anc[L, rf.species_idx])
                if len(below) < 2:
                    continue
                x, y = a_[L][below].copy(), b_[L][below].copy()
                x[x < 0] = -1 - np.arange((x < 0).sum()); y[y < 0] = -1 - np.arange((y < 0).sum())
                ami = adjusted_mutual_info_score(x, y) if (len(set(x)) > 1 or len(set(y)) > 1) else 1.0
                fh.write(f"{rf.family.id}\t{rf.family.collection or ''}\t{st.names[L]}\t{int(st.depth[L])}\t"
                         f"{int(L == rf.family_level)}\t{len(below)}\t{ami:.4f}\n")
    print("done", P, mode, file=sys.stderr)


if __name__ == "__main__":
    main()
