#!/usr/bin/env python
"""Pre-registered subset-vs-full equivalence (DESIGN sec. 7), panels Y and K, per tool run.
(a) TOST on paired per-family deltas (subset - full) of the 4 primary endpoints: equivalent iff the 90% paired
    bootstrap CI of the mean delta lies inside +-0.05 (units = reference families/units with a defined value in both modes;
    10,000 resamples, seed hobench.seeds.SEEDS[0]).
(b) per-pair merge-level identity: over all member pairs of the scored families, the fraction whose predicted merge
    level (deepest level at which both genes share a HOG; none = unmerged) is identical in the full and subset runs.
Writes results/equivalence.md.  Usage: python scripts/equivalence.py [--panels Y K]"""
import argparse, csv, sys
from collections import defaultdict
from itertools import combinations
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / "scripts"))
from hobench.metrics import PRIMARY                      # noqa: E402
from hobench.seeds import SEEDS                          # noqa: E402
from hobench.ids import make_normalizer                  # noqa: E402
from hobench.reference import load_families             # noqa: E402
from hobench.speciestree import SpeciesTree              # noqa: E402
import evaluate_panel as ep                              # noqa: E402

MARGIN, NB = 0.05, 10000
NOTE = '\nNotes: K event recovery has n=1 event (K02 TU_G_dup), so a delta of +1.0 is a single flip, not an estimate; treat as inconclusive. Y membership F1 is NA (no complete cells/decoys). Pair identity is computed once per tool (shown on the level_bcubed_f1 row). FastOMA variants: declared/late; dendroblast and OF variants other than default have no subset runs.\n'


def read_pf(P, mode):
    d = defaultdict(dict)
    with open(ROOT / "results" / f"{P}_{mode}" / "per_family.tsv") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["metric"] in PRIMARY:
                d[(r["method"], r["variant"], r["metric"])][r["family"]] = float(r["value"])
    return d


def merge_levels(pred, st, genes_by_fam):
    order = sorted(pred.levels, key=lambda L: st.depth_of(L))          # root first
    out = {}
    for fid, genes in genes_by_fam.items():
        lab = {L: pred.levels[L] for L in order}
        lv = np.full((len(genes), len(genes)), -1, dtype=np.int16)
        for k, L in enumerate(order):
            h = [lab[L].get(g) for g in genes]
            for i, j in combinations(range(len(genes)), 2):
                if h[i] is not None and h[i] == h[j]:
                    lv[i, j] = st.depth_of(L)
        out[fid] = lv[np.triu_indices(len(genes), 1)]
    return out, order


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--panels", nargs="+", default=["Y", "K"])
    a = ap.parse_args()
    rng = np.random.default_rng(SEEDS[0])
    md = ["# Subset vs full equivalence (pre-registered; TOST margin +-0.05, 90% paired bootstrap CI)", "",
          "Delta = subset - full, paired per reference family/unit (all tracks; collections are bootstrapped over their units). "
          "Equivalent = 90% CI inside +-0.05. Pair identity = fraction of member pairs with identical predicted merge level.", "",
          "| panel | method | variant | metric | n | mean delta | 90% CI | equivalent | pair identity (n pairs) |", "|---|---|---|---|---|---|---|---|---|"]
    for P in a.panels:
        if not (ROOT / "results" / f"{P}_subset" / "per_family.tsv").exists():
            md.append(f"| {P} | subset results missing |  |  |  |  |  |  |  |"); continue
        full, sub = read_pf(P, "full"), read_pf(P, "subset")
        st = SpeciesTree.from_newick(str(ROOT / "panels" / P / "species_tree.nwk"))
        fams = [f for f in load_families(str(ROOT / "families"), panel=P)]
        genes = {f.id: sorted(f.member_ids) for f in fams}
        keep, cs = set(), set()
        for f in fams:
            keep |= set(f.member_ids) | set(f.decoys); cs |= set(f.complete_species)
        kw = {"normalize": make_normalizer("auto"), "keep_genes": keep, "count_species": cs}
        ident = {}
        for mode in ("full", "subset"):
            for run in sorted(p for p in (ROOT / "runs" / P / mode).iterdir() if (p / "manifest.json").exists()):
                for tool, var, pred, _, _ in ep.load_run(run, st, kw, Path("/tmp/hobench_eval_cache"), ["declared"]):
                    ident[(mode, tool, var)] = merge_levels(pred, st, genes)
        for (meth, var, metric), fv in sorted(full.items()):
            sv = sub.get((meth, var, metric))
            if sv is None or meth in ("positive_control", "all_singletons", "one_big_hog", "random_relabel"):
                continue
            ids = sorted(set(fv) & set(sv))
            if not ids:
                continue
            d = np.array([sv[i] - fv[i] for i in ids])
            idx = rng.integers(0, len(d), size=(NB, len(d)))
            bm = d[idx].mean(axis=1)
            lo, hi = np.percentile(bm, [5, 95])
            eq = lo > -MARGIN and hi < MARGIN
            pid = ""
            if metric == PRIMARY[0] and ("full", meth, var) in ident and ("subset", meth, var) in ident:
                a_, b_ = ident[("full", meth, var)][0], ident[("subset", meth, var)][0]
                same = tot = 0
                for fid in a_:
                    same += int((a_[fid] == b_[fid]).sum()); tot += len(a_[fid])
                pid = f"{same / tot:.4f} ({tot})" if tot else "NA"
            md.append(f"| {P} | {meth} | {var} | {metric} | {len(d)} | {d.mean():+.4f} | [{lo:+.4f}, {hi:+.4f}] | "
                      f"{'YES' if eq else 'no'} | {pid} |")
    (ROOT / "results" / "equivalence.md").write_text("\n".join(md) + "\n" + NOTE)
    print("\n".join(md))


if __name__ == "__main__":
    main()
