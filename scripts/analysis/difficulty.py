#!/usr/bin/env python
"""Per-family difficulty vs tool agreement -> results/analysis/difficulty.{md,tsv} + difficulty_vs_agreement.png.
difficulty = 1 - mean skill over {OrthoFinder default, FastOMA declared} x primary endpoints defined for the item
(skill = (s-b)/(1-b) from per_family.tsv; unbounded below, so difficulty can exceed 1).
agreement = AMI of the two tools' HOG partitions of the reference-tree genes (results/analysis/cache/ami_*.tsv,
scripts/analysis/extract_nodes.py): at the family level, and mean over all levels at/below it.
All scored families (any track); a collection is one item (unit means). Scatter: headline items only."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from common import OUT, OKABE, PRIMARY, available, items, read_tsv, fmt, write_md


TOOLS3 = [("orthofinder", "default"), ("fastoma", "default_declared"), ("orthofinder", "dendroblast")]


def main():
    rows, noskill = [], []
    for P, mode in available():
        it = items(P, mode, headline_only=False)
        hi_ = items(P, mode)
        if hi_.empty:
            noskill.append(f"(panel {P} {mode}: no headline non-control item with any scored metric)")
        noskill += sorted(set(hi_["item"]) - set(hi_[hi_.metric.str.startswith("skill_")]["item"]))
        pf = read_tsv(P, mode, "per_family", usecols=["family", "collection", "track", "headline"]).drop_duplicates("family")
        pf["item"] = np.where((pf.track == "collection") & (pf.collection.fillna("") != ""), "collection:" + pf.collection.fillna(""), pf.family)
        hl = pf.groupby("item").agg(track=("track", "first"), headline=("headline", "max"))
        sk = it[it.metric.str.startswith("skill_") & it.metric.str[6:].isin(PRIMARY)
                & it.set_index(["method", "variant"]).index.isin(TOOLS3)]
        d = sk.groupby("item").agg(mean_skill=("value", "mean"), n_skill_values=("value", "size"),
                                   n_endpoints=("metric", "nunique"), n_methods=("method", lambda x: len(set(zip(x, sk.loc[x.index, "variant"])))))
        for e in PRIMARY:
            for (m, v) in TOOLS3[:2]:
                x = it[(it.metric == e) & (it.method == m) & (it.variant == v)].set_index("item")["value"]
                d[f"{m}_{e}"] = x
        f = OUT / "cache" / f"ami_{P}_{mode}.tsv"
        if f.exists():
            am = pd.read_csv(f, sep="\t", keep_default_na=False, dtype={"collection": str})
            am["item"] = am.family.map(dict(zip(pf.family, pf["item"])))
            fd = am[am.is_family_level == 1].set_index("family").depth
            am = am[am.depth >= am.family.map(fd)]          # family level and below only
            u = am.groupby(["item", "family"]).apply(lambda g: pd.Series({
                "ami_family_level": g.loc[g.is_family_level == 1, "ami"].mean(), "ami_mean_levels": g.ami.mean(),
                "ami3_family_level": g.loc[g.is_family_level == 1, "ami3_mean"].mean(), "ami3_mean_levels": g.ami3_mean.mean()}),
                include_groups=False)
            d = d.join(u.groupby("item").mean(), how="left")
        d = d.join(hl, how="left")
        d.insert(0, "mode", mode); d.insert(0, "panel", P)
        d["difficulty"] = 1 - d.mean_skill
        rows.append(d.reset_index())
    D = pd.concat(rows).sort_values("difficulty", ascending=False)
    D.to_csv(OUT / "difficulty.tsv", sep="\t", index=False, float_format="%.4f")
    HH = D[(D.track != "control") & (D.headline == 1)]
    H = HH.dropna(subset=["difficulty", "ami3_mean_levels"])
    dropped = HH[~HH.index.isin(H.index)]
    fig, ax = plt.subplots(figsize=(5.2, 3.8))
    for k, (P, g) in enumerate(H.groupby("panel")):
        ax.scatter(g.ami3_mean_levels, g.difficulty, s=28, color=OKABE[k], label=f"{P} (n={len(g)})",
                   marker="osD^v"[k % 5], edgecolor="black", lw=.4)
        for _, r in g.iterrows():
            ax.annotate(r["item"].replace("collection:", ""), (r.ami3_mean_levels, r.difficulty), fontsize=6,
                        xytext=(3, 2), textcoords="offset points")
    ax.set_xlabel("3-way agreement: mean pairwise AMI (OF, OF-dendroblast, FastOMA), mean over levels")
    ax.set_ylabel("difficulty = 1 - mean skill (3 methods)")
    ax.legend(fontsize=7, frameon=False); ax.grid(alpha=.3)
    fig.tight_layout(); fig.savefig(OUT / "difficulty_vs_agreement.png", dpi=110); plt.close(fig)
    def rho(col):
        x = H.dropna(subset=[col])
        if len(x) < 5:
            return f"NA (n={len(x)})"
        r = spearmanr(x[col], x.difficulty)
        return f"{r.statistic:.2f} (p={r.pvalue:.2g}, n={len(x)})"
    L = ["# Family difficulty and tool agreement", "",
         "difficulty = 1 - mean skill over the (method, primary endpoint) values DEFINED for the item, methods = OF default, "
         "FastOMA declared, OF dendroblast (skill unbounded below, so difficulty > 1 occurs when tools are worse than the best "
         "trivial baseline). AMI = adjusted mutual information of two methods' HOG partitions of reference-tree genes (genes without a "
         "HOG at a level = singletons), at the family level and averaged over levels at/below it; 3-way agreement = mean of the three "
         "pairwise AMIs. Collections: unit means. Generated by scripts/analysis/difficulty.py.", "",
         f"Headline non-control items: {len(HH)}; with difficulty and 3-way AMI: {len(H)}"
         + (f" (dropped: {', '.join(dropped['item'])} - no AMI)" if len(dropped) else "") + ". Headline items with no defined "
         f"skill on any primary endpoint (excluded; e.g. membership NA after E5 or best baseline = 1): {', '.join(noskill) or 'none'}.", "",
         "Spearman rho(difficulty, agreement): 3-way AMI mean over levels " + rho("ami3_mean_levels")
         + "; 3-way AMI family level " + rho("ami3_family_level") + "; OF-FastOMA AMI mean over levels " + rho("ami_mean_levels")
         + "; OF-FastOMA AMI family level " + rho("ami_family_level") + ".", "",
         "![difficulty vs agreement](difficulty_vs_agreement.png)", "",
         "Ranking below: non-control items (headline flag shown); control-track items (K01 RPROT units, K02, E05) only in difficulty.tsv.", "",
         "| rank | item | panel | track | headline | difficulty | mean skill (n values; n endpoints) | AMI OF-FO family / mean levels | 3-way AMI family / mean levels | "
         "OF/FO bcubed | OF/FO merge acc | OF/FO event rec | OF/FO membership |", "|" + "---|" * 13]
    L += ["", "Ranking: non-control items (headline flag shown); control-track items (K01 units, K02, E05) are in difficulty.tsv.", ""] if False else []
    for k, (_, r) in enumerate(D[D.track != "control"].iterrows(), 1):
        pair = lambda e: f"{fmt(r.get('orthofinder_' + e))}/{fmt(r.get('fastoma_' + e))}"
        L.append(f"| {k} | {r['item']} | {r.panel} | {r.track} | {int(r.headline) if r.headline == r.headline else ''} | "
                 f"{fmt(r.difficulty)} | {fmt(r.mean_skill)} ({int(r.n_skill_values) if r.n_skill_values == r.n_skill_values else 0}; "
                 f"{int(r.n_endpoints) if r.n_endpoints == r.n_endpoints else 0}) | "
                 f"{fmt(r.get('ami_family_level'))} / {fmt(r.get('ami_mean_levels'))} | {fmt(r.get('ami3_family_level'))} / "
                 f"{fmt(r.get('ami3_mean_levels'))} | {pair('level_bcubed_f1')} | "
                 f"{pair('merge_accuracy_node')} | {pair('event_recovery_0.5')} | {pair('membership_f1')} |")
    write_md(OUT / "difficulty.md", L)


if __name__ == "__main__":
    main()
