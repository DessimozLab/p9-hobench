#!/usr/bin/env python
"""Error anatomy case studies -> results/analysis/error_anatomy.md (+ error_anatomy_*.tsv).
(a) P02_bHLH membership (P subset): members' predicted family-level HOGs (same rule as hobench.metrics.family_integrity:
    the member's HOG at the highest level <= family level; size < 2 = missing), their sizes and the non-member genes
    in them, classified by panels/P/subset/members.tsv category (neighbourhood = DIAMOND hit of a family member,
    domain_only, background) or membership of another reference family. Predictions read WITHOUT gene filtering.
(b) E panel merge accuracy (E subset): over-split vs over-merged pairs per reference node level (cache/nodes_E_subset).
(c) K02 EF-Tu/EF-G (K full): per copy, predicted HOGs at the family level (LUCA) and their composition.
(d) YGOB testable WGD events (Y full events.tsv): failure side among non-recovered events.
Usage: python scripts/analysis/error_anatomy.py"""
import csv, sys
from collections import Counter, defaultdict
from pathlib import Path
import numpy as np
import pandas as pd
from common import ROOT, OUT, read_tsv, fmt, write_md
import evaluate_panel as ep
from hobench.ids import make_normalizer
from hobench.metrics import reference_level_codes, panel_context
from hobench.reference import load_families, reconcile
from hobench.speciestree import SpeciesTree

TOOLS = [("orthofinder", "default"), ("fastoma", "default_declared")]
CACHE = Path("/tmp/hobench_eval_cache")


def load(P, mode, filtered=True):
    st = SpeciesTree.from_newick(str(ROOT / "panels" / P / "species_tree.nwk"))
    fams = load_families(str(ROOT / "families"), panel=P)
    if mode == "subset":
        sub = {r["id"] for r in csv.DictReader(open(ROOT / "panels" / P / "subset" / "members.tsv"), delimiter="\t")}
        fams = [f for f in fams if set(f.member_ids) <= sub]
    keep, cs = set(), set()
    for f in fams:
        keep |= set(f.member_ids) | set(f.decoys); cs |= set(f.complete_species)
    kw = {"normalize": make_normalizer("auto"), "keep_genes": keep if filtered else None, "count_species": cs}
    preds = {}
    for run in sorted(p for p in (ROOT / "runs" / P / mode).iterdir() if (p / "manifest.json").exists()):
        if run.name in ("orthofinder_default", "fastoma_default"):
            for tool, variant, pred, *_ in ep.load_run(run, st, kw, CACHE, ["declared"]):
                preds[(tool, variant)] = pred
    return st, fams, preds


def top_hog(st, pred, g, sp, LF):
    top = None
    for lv in st.ancestors(sp, include_self=False):
        if not st.anc[LF, st.index[lv]]:
            break
        h = pred.hog(lv, g)
        if h is not None:
            top = (lv, h)
    return top if top is not None and pred.hog_size(*top) >= 2 else None


def case_a():
    st, fams, preds = load("P", "subset", filtered=False)
    F = next(f for f in fams if f.id == "P02_bHLH"); rf = reconcile(F, st, check=False)
    ctx = panel_context(fams, st)
    cat = {r["id"]: r["category"] for r in csv.DictReader(open(ROOT / "panels/P/subset/members.tsv"), delimiter="\t")}
    LF = rf.family_level; lab = reference_level_codes(rf, [LF])[LF]
    refog = {g: int(lab[k]) for k, g in enumerate(rf.genes) if lab[k] >= 0}
    gsp = F.gene_species; members = set(F.member_ids); complete = set(F.complete_species)
    L = ["## (a) P02_bHLH: FastOMA membership F1 ~0.02 (P subset)", "",
         f"Family level {st.names[LF]}; {len(F.member_ids)} members, {len(rf.genes)} in the reference tree, "
         f"{len(set(refog.values()))} reference HOGs (RefOGs) at the family level; complete species {sorted(complete)}.", "",
         "| tool | members missing (no HOG >=2) | distinct predicted HOGs holding members | members per predicted HOG: median / max | "
         "largest RefOG covered by one HOG (frac) | non-members in those HOGs: total / complete-species / other ref. family / "
         "neighbourhood / domain_only / background |", "|---|---|---|---|---|---|"]
    rows, sub_lines = [], []
    for (m, v), pred in preds.items():
        key = {g: top_hog(st, pred, g, gsp[g], LF) for g in refog}
        miss = sum(k is None for k in key.values())
        hogs = Counter(k for k in key.values() if k is not None)
        nm = Counter()
        for k in hogs:
            for g in pred.hog_members(*k):
                if g in members:
                    continue
                sp = pred.gene_species.get(g)
                nm["total"] += 1
                nm["complete"] += sp in complete
                fam = ctx.gene_family.get(g)
                nm["other_family" if fam else cat.get(g, "not_in_subset_table")] += 1
                rows.append({"tool": f"{m}/{v}", "hog_level": k[0], "hog": k[1], "gene": g, "species": sp,
                             "class": "other_family:" + fam if fam else cat.get(g, "NA")})
        big = Counter(refog.values()).most_common(1)[0]
        bigcov = max(Counter(key[g] for g in refog if refog[g] == big[0] and key[g]).values() or [0]) / big[1]
        sz = list(hogs.values())
        L.append(f"| {m}/{v} | {miss}/{len(refog)} | {len(hogs)} | {np.median(sz) if sz else 0:.0f} / {max(sz) if sz else 0} | "
                 f"{bigcov:.3f} (RefOG n={big[1]}) | {nm['total']} / {nm['complete']} / {nm['other_family']} / {nm['neighbourhood']} / "
                 f"{nm['domain_only']} / {nm['background']} |")
        root = rf.events[0]
        fr = []
        for leaves in root.child_leaves:
            if len(leaves) < 2:
                continue
            c = Counter(key[rf.genes[k]] for k in leaves)
            fr.append(max((n_ for k_, n_ in c.items() if k_ is not None), default=0) / len(leaves))
        sub_lines.append(f"| {m}/{v} | {len(fr)} | {np.mean(fr):.3f} | {np.median(fr):.3f} | {sum(x == 1 for x in fr)} |")
    L += ["", f"Reference root is ev=U ({rf.events[0].id}: unresolved polytomy of subfamilies). Membership scores the whole "
          "superfamily as ONE family-level RefOG. Clade view (children of the U root with >=2 genes; subfamilies are nested below them): largest fraction "
          "of a subfamily's members inside one predicted HOG:", "",
          "| tool | n root-child clades | mean | median | n fully in one HOG |", "|---|---|---|---|---|"] + sub_lines
    pd.DataFrame(rows).to_csv(OUT / "error_anatomy_P02_nonmembers.tsv", sep="\t", index=False)
    return L


def case_b(P="E", tag="(b) E panel"):
    nd = pd.read_csv(OUT / f"cache/nodes_{P}_subset.tsv.gz", sep="\t", keep_default_na=False, dtype={"variant": str})
    hl = read_tsv(P, "subset", "per_family", usecols=["family", "track", "headline"]).drop_duplicates("family")
    keep = set(hl[(hl.track != "control") & (hl.headline == 1)].family)
    nd = nd[nd.family.isin(keep) & nd.method.isin(["orthofinder", "fastoma"])]
    nd["acc"] = nd.n_correct / nd.n_pairs
    L = [f"## {tag}: FastOMA vs OrthoFinder merge-level accuracy ({P} subset, headline families/units: " + (", ".join(sorted(keep)) if len(keep) < 12 else f"{len(keep)}") + ")", "",
         "Node-balanced accuracy and, over incorrect pairs, the share that is over-split (merged too late or never) vs "
         "over-merged (merged below the reference interval). Pair counts pooled over families.", "",
         "| tool | node ev | n nodes | node-bal. acc | pairs | frac over-split | frac over-merged |", "|---|---|---|---|---|---|---|"]
    for (m, ev), g in nd.groupby(["method", "ev"]):
        L.append(f"| {m} | {ev} | {len(g)} | {fmt(g.acc.mean())} | {g.n_pairs.sum()} | {fmt(g.n_over_split.sum() / g.n_pairs.sum())} | "
                 f"{fmt(g.n_over_merged.sum() / g.n_pairs.sum())} |")
    nd["grp"] = np.where(nd.collection.isin(["XLAEVIS", "YGOB"]), "units:" + nd.collection, nd.family)
    L += ["", "By item (node-balanced accuracy, mean over units for a collection; pairs pooled):", "",
          "| item | tool | node ev | n nodes | acc | pairs | frac over-split | frac over-merged |", "|---|---|---|---|---|---|---|---|"]
    for (gname, m, ev), g in nd.groupby(["grp", "method", "ev"]):
        L.append(f"| {gname} | {m} | {ev} | {len(g)} | {fmt(g.groupby('family').acc.mean().mean())} | {g.n_pairs.sum()} | "
                 f"{fmt(g.n_over_split.sum() / g.n_pairs.sum())} | {fmt(g.n_over_merged.sum() / g.n_pairs.sum())} |")
    t = nd.groupby(["method", "node_level", "node_depth"]).agg(n_nodes=("acc", "size"), acc=("acc", "mean"), pairs=("n_pairs", "sum"),
                                                            split=("n_over_split", "sum"), merged=("n_over_merged", "sum")).reset_index()
    t["frac_split"] = t.split / t.pairs; t["frac_merged"] = t.merged / t.pairs
    t.to_csv(OUT / f"error_anatomy_{P}_by_level.tsv", sep="\t", index=False, float_format="%.4f")
    w = t.pivot_table(index=["node_depth", "node_level"], columns="method", values=["acc", "frac_split", "frac_merged", "n_nodes"])
    L += ["", "By reference node level (S and D nodes together; levels with >=3 nodes for both tools):", "",
          "| depth | level | n nodes | acc OF / FO | over-split OF / FO | over-merged OF / FO |", "|---|---|---|---|---|---|"]
    for (d, lv), r in w.iterrows():
        if min(r[("n_nodes", "orthofinder")], r[("n_nodes", "fastoma")]) < 3:
            continue
        L.append(f"| {d} | {lv} | {int(r[('n_nodes', 'fastoma')])} | {fmt(r[('acc', 'orthofinder')])} / {fmt(r[('acc', 'fastoma')])} | "
                 f"{fmt(r[('frac_split', 'orthofinder')])} / {fmt(r[('frac_split', 'fastoma')])} | "
                 f"{fmt(r[('frac_merged', 'orthofinder')])} / {fmt(r[('frac_merged', 'fastoma')])} |")
    return L


def case_c():
    st, fams, preds = load("K", "full")
    F = next(f for f in fams if f.id == "K02_trGTPases"); rf = reconcile(F, st, check=False)
    ev = next(e for e in rf.events if e.id == "TU_G_dup")
    LF = rf.family_level; gsp = F.gene_species
    L = ["## (c) K02 EF-Tu / EF-G duplication (K full; control track)", "",
         f"Event TU_G_dup (D at {st.names[ev.mapping]}); family level {st.names[LF]}. Per copy: predicted HOGs at the family "
         "level, i.e. each gene's HOG at the highest level <= family level (as in membership scoring), holding the copy's genes (HOG id: n copy genes [domains]; other-copy genes in the same HOG).", ""]
    dom = lambda s: "A" if st.anc[st.index["Archaea"], st.index[s]] else ("B" if st.anc[st.index["Bacteria"], st.index[s]] else "?")
    for (m, v), pred in preds.items():
        L.append(f"**{m}/{v}**")
        copy_of = {}
        for c, leaves in enumerate(ev.child_leaves):
            for k in leaves:
                copy_of[rf.genes[k]] = c
        for c, leaves in enumerate(ev.child_leaves):
            hs = defaultdict(list)
            for k in leaves:
                g = rf.genes[k]; h = top_hog(st, pred, g, gsp[g], LF)
                hs[h].append(g)
            parts = []
            for h, gs in sorted(hs.items(), key=lambda x: -len(x[1])):
                other = 0 if h is None else sum(copy_of.get(x) not in (None, c) for x in pred.hog_members(*h))
                doms = "".join(sorted(Counter(dom(gsp[g]) for g in gs).elements()))
                d2 = Counter(doms)
                parts.append(f"{'no HOG >=2' if h is None else h[0] + ':' + h[1]}: {len(gs)} [A{d2['A']}/B{d2['B']}]" + (f"; +{other} other-copy" if other else ""))
            name = "EF-Tu" if c == 0 else "EF-G"
            L.append(f"- copy {c} ({name}, {len(leaves)} genes): " + "; ".join(parts[:6]) + (" ..." if len(parts) > 6 else ""))
        L.append("")
    return L


def case_d():
    ev = read_tsv("Y", "full", "events")
    ev = ev[(ev.collection == "YGOB") & (ev.testable == 1) & (ev.event_branch == "internal")]
    L = ["## (d) YGOB testable WGD events (Y full)", "",
         "Recovery (theta=0.5) needs copies apart at the event level (separation = 1) AND each copy's reference-together "
         "pairs >= 50% together (cohesion). Failure side among non-recovered events:", "",
         "| method | testable events | recovered | not recovered | apart fails only | cohesion fails only | both fail | "
         "modal predicted merge level of cross-copy pairs (non-recovered, top 3) |", "|---|---|---|---|---|---|---|---|"]
    for (m, v) in [("orthofinder", "default"), ("orthofinder", "dendroblast"), ("fastoma", "default_declared"), ("of_flat", "default")]:
        s = ev[(ev.method == m) & (ev.variant.fillna("") == v)]
        if s.empty:
            continue
        nr = s[s["recovered_0.5"] == 0]
        ap = nr.separation < 1 - 1e-9; co = nr.cohesion_min < 0.5
        top = ", ".join(f"{k} {c}" for k, c in nr.modal_pred_merge.value_counts().head(3).items())
        L.append(f"| {m}/{v} | {len(s)} | {int(s['recovered_0.5'].sum())} | {len(nr)} | {int((ap & ~co).sum())} | "
                 f"{int((~ap & co).sum())} | {int((ap & co).sum())} | {top} |")
    return L


def main():
    L = ["# Error anatomy (case studies)", "", "Generated by scripts/analysis/error_anatomy.py; gene IDs are UniProt accessions.", ""]
    case_v = lambda: case_b("V", "(e) V panel")
    for f in (case_b, case_v, case_d, case_c, case_a) if len(sys.argv) < 2 else [globals()[x] for x in sys.argv[1:]]:
        L += f() + [""]
        write_md(OUT / "error_anatomy.md", L)


if __name__ == "__main__":
    main()
