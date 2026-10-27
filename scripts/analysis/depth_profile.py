#!/usr/bin/env python
"""Depth profile of per-level B-cubed F1 (headline track) per panel and method.
Depth = edges from the panel species-tree root (column `depth` of per_level.tsv). Units of a collection are averaged
into one item per (level); then items are averaged per depth. Writes results/analysis/depth_profile.tsv and
results/analysis/depth_profile_<P>.png."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from common import OUT, OKABE, METHOD_ORDER, available, read_tsv, item_of, mlabel

SHOW = [("orthofinder", "default"), ("orthofinder", "dendroblast"), ("fastoma", "default_declared"),
        ("of_flat", "default"), ("all_singletons", ""), ("one_big_hog", "")]
STY = {k: (OKABE[i], ls) for i, (k, ls) in enumerate(zip(SHOW, ["-", "--", "-", ":", "-.", "-."]))}


def main():
    rows = []
    for P, mode in available():
        pf = read_tsv(P, mode, "per_family", usecols=["family", "track", "headline"]).drop_duplicates("family")
        keep = set(pf[(pf.track != "control") & (pf.headline == 1)].family)
        trk = dict(zip(pf.family, pf.track))
        d = read_tsv(P, mode, "per_level", usecols=["family", "collection", "method", "variant", "level", "depth",
                                                     "n_genes", "bcubed_f1"])
        d["variant"] = d["variant"].fillna("")
        d = d[d.family.isin(keep)].dropna(subset=["bcubed_f1"]).copy()
        d["track"] = d.family.map(trk)
        d["item"] = item_of(d)
        it = d.groupby(["method", "variant", "item", "level", "depth"]).agg(f1=("bcubed_f1", "mean"),
                                                                           genes=("n_genes", "sum")).reset_index()
        g = it.groupby(["method", "variant", "depth"]).agg(bcubed_f1=("f1", "mean"), n_items=("item", "nunique"),
                                                           n_item_levels=("f1", "size"), n_genes=("genes", "sum"),
                                                           levels=("level", lambda x: ",".join(sorted(set(x))[:6]))).reset_index()
        g.insert(0, "mode", mode); g.insert(0, "panel", P)
        rows.append(g)
        if g.empty:   # e.g. K: no headline item has asserted pairs
            (OUT / f"depth_profile_{P}.png").unlink(missing_ok=True); continue
        fig, ax = plt.subplots(figsize=(5.5, 3.6))
        for k in SHOW:
            s = g[(g.method == k[0]) & (g.variant == k[1])].sort_values("depth")
            if s.empty:
                continue
            c, ls = STY[k]
            ax.plot(s.depth, s.bcubed_f1, ls, color=c, marker="o", ms=4, lw=1.6, label=mlabel(*k))
        ax.set_xlabel("level depth (edges from panel species-tree root)")
        ax.set_ylabel("per-level B-cubed F1 (mean over items)")
        nit = g.groupby("depth").n_items.max()
        ax.set_title(f"Panel {P} ({mode}): headline items per depth " + ",".join(f"{int(x)}:{int(n)}" for x, n in nit.items()),
                     fontsize=8)
        ax.set_ylim(0, 1.02); ax.grid(alpha=.3)
        ax.legend(fontsize=7, frameon=False, loc="lower left")
        fig.tight_layout(); fig.savefig(OUT / f"depth_profile_{P}.png", dpi=110); plt.close(fig)
    pd.concat(rows).to_csv(OUT / "depth_profile.tsv", sep="\t", index=False, float_format="%.4f")
    print("wrote depth_profile.tsv + PNGs")


if __name__ == "__main__":
    main()
