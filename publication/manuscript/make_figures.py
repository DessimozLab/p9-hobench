#!/usr/bin/env python
"""Figures for the HOBench manuscript. Reads committed result tables only (results/analysis, results/robustness).

Run from the project root:  python publication/manuscript/make_figures.py
Writes publication/manuscript/figures/Fig{1..4}.{pdf,png}.
Palette: validated categorical slots (blue/orange/aqua; dataviz validator, light mode, all checks pass; contrast warning
for aqua handled with direct labels and distinct marker shapes). Baselines in neutral grey.
"""
from __future__ import annotations

import ast
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

ROOT = Path(__file__).resolve().parents[2] if (Path(__file__).resolve().parents[2] / "results").exists() else Path(__file__).resolve().parents[1]
A = ROOT / "results" / "analysis"
OUT = Path(__file__).resolve().parent / "figures"
OUT.mkdir(exist_ok=True)

INK, INK2, GRID, SURF = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"
M = {  # method -> (label, colour, marker)
    "orthofinder/default": ("OrthoFinder", "#2a78d6", "o"),
    "orthofinder/dendroblast": ("OrthoFinder-DB", "#1baf7a", "s"),
    "fastoma/default_declared": ("FastOMA", "#eb6834", "D"),
}
BASE = "#9a9994"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8, "axes.edgecolor": INK2, "axes.labelcolor": INK,
                     "xtick.color": INK2, "ytick.color": INK2, "axes.spines.top": False,
                     "axes.spines.right": False, "figure.facecolor": "white", "savefig.dpi": 300})


def save(fig, name):
    fig.savefig(OUT / f"{name}.pdf", bbox_inches="tight")
    fig.savefig(OUT / f"{name}.png", bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------- Fig 1: workflow + merge-level concept
def fig1():
    fig = plt.figure(figsize=(7.5, 3.6))
    ax = fig.add_axes([0.0, 0.0, 0.56, 1.0]); ax.set_axis_off(); ax.set_xlim(0, 10); ax.set_ylim(0, 10)
    steps = [
        ("Literature survey", "verified citations;\n~85 candidate families"),
        ("Paper retrieval", "Europe PMC full text,\nsupplements, human requests"),
        ("Extraction by parsers", "IDs, labels, trees;\nquotes for every event"),
        ("Mapping to panels", "same species; ≥98% id,\n≥80% coverage, RBH"),
        ("Reference tree", "S / D / H / U nodes,\nlevel intervals"),
        ("Verification & freeze", "blind re-derivation,\naudit, checksums"),
    ]
    y = 9.3
    for i, (t, d) in enumerate(steps):
        ax.add_patch(FancyBboxPatch((0.3, y - 1.05), 9.2, 1.15, boxstyle="round,pad=0.02,rounding_size=0.15",
                                    fc="#eef4fc" if i % 2 == 0 else "#f6f5f2", ec=INK2, lw=0.6))
        ax.text(0.6, y - 0.47, t, fontsize=8, fontweight="bold", va="center", color=INK)
        ax.text(5.0, y - 0.47, d, fontsize=7, va="center", color=INK2)
        if i < len(steps) - 1:
            ax.add_patch(FancyArrowPatch((2.0, y - 1.08), (2.0, y - 1.42), arrowstyle="-|>", mutation_scale=7,
                                         color=INK2, lw=0.8))
        y -= 1.5
    ax.text(0.3, 9.75, "A", fontsize=11, fontweight="bold")

    # B: species tree with a WGD and merge-level intervals
    bx = fig.add_axes([0.6, 0.08, 0.4, 0.86]); bx.set_axis_off(); bx.set_xlim(0, 10); bx.set_ylim(0, 10)
    bx.text(0, 9.6, "B", fontsize=11, fontweight="bold")
    lw = 1.2
    # species tree: ((s1,s2)L1,(s3,s4)L2)L0 ; leaves at y=1..4
    leaves = {"sp1": 8.0, "sp2": 6.5, "sp3": 5.0, "sp4": 3.5}
    for s, yy in leaves.items():
        bx.plot([6.0, 8.0], [yy, yy], color=INK2, lw=lw); bx.text(8.2, yy, s, va="center", fontsize=7, color=INK)
    bx.plot([6.0, 6.0], [8.0, 6.5], color=INK2, lw=lw); bx.plot([6.0, 6.0], [5.0, 3.5], color=INK2, lw=lw)
    bx.plot([4.0, 6.0], [7.25, 7.25], color=INK2, lw=lw); bx.plot([4.0, 6.0], [4.25, 4.25], color=INK2, lw=lw)
    bx.plot([4.0, 4.0], [7.25, 4.25], color=INK2, lw=lw); bx.plot([2.8, 4.0], [5.75, 5.75], color=INK2, lw=lw)
    bx.text(6.15, 7.45, "L1", fontsize=7, color=INK2); bx.text(6.15, 4.45, "L2", fontsize=7, color=INK2)
    bx.text(3.6, 6.1, "L0", ha="right", fontsize=7, color=INK2)
    bx.plot([5.0], [7.25], marker="*", ms=9, color="#eb6834")
    bx.text(5.0, 7.75, "WGD", ha="center", fontsize=7, color=INK)
    rows = [("a1–a2 (sp1–sp2 orthologs)", "merge = [L1, L1]"),
            ("a1–b1 (ohnologs, WGD on L1 branch)", "merge = [L0, L0]"),
            ("a1–c (sp1–sp3, unresolved)", "not scored"),
            ("timing contested", "merge = [L1, L0] interval")]
    yy = 2.0
    for t, v in rows:
        bx.text(0.0, yy, t, fontsize=6.5, color=INK); bx.text(10.0, yy, v, fontsize=6.5, color=INK2, ha="right")
        yy -= 0.62
    bx.text(0.0, 2.6, "Reference merge levels of gene pairs:", fontsize=7, fontweight="bold", color=INK)
    save(fig, "Fig1")


# ---------------------------------------------------------------- Fig 2: primary endpoints by panel
def fig2():
    t = pd.read_csv(A / "main_table.tsv", sep="\t")
    t["key"] = t["method"] + "/" + t["variant"].fillna("")
    ends = [("level_bcubed_f1", "Per-level B-cubed F1"), ("merge_accuracy_node", "Merge-level accuracy"),
            ("event_recovery_0.5", "Event recovery"), ("membership_f1", "Membership F1")]
    panels = ["Y", "E", "V", "P"]
    fig, axs = plt.subplots(1, 4, figsize=(7.5, 2.6), sharey=True)
    for ax, (m, title) in zip(axs, ends):
        for i, P in enumerate(panels):
            sub = t[t.panel == P]
            bl = sub[sub.method.isin(["all_singletons", "one_big_hog", "random_relabel"])][m].max()
            if pd.notna(bl):
                ax.plot([bl, bl], [i - 0.32, i + 0.32], color=BASE, lw=2, solid_capstyle="butt")
            for j, (k, (lab, col, mk)) in enumerate(M.items()):
                r = sub[sub.key == k]
                if r.empty or pd.isna(r[m].iloc[0]):
                    continue
                n = int(r[f"n_{m}"].iloc[0])
                ax.plot(r[m].iloc[0], i + (j - 1) * 0.18, marker=mk, ms=4.5, color=col, mec=SURF, mew=0.5,
                        ls="none")
        ax.set_title(title, fontsize=8, color=INK)
        ax.set_xlim(-0.02, 1.02); ax.grid(axis="x", color=GRID, lw=0.5); ax.set_axisbelow(True)
        ax.set_yticks(range(len(panels))); ax.set_yticklabels([f"{p} ({'full' if p == 'Y' else 'subset'})"
                                                              for p in panels])
    axs[0].invert_yaxis()
    h = [plt.Line2D([], [], marker=mk, color=col, ls="none", ms=5, label=lab) for lab, col, mk in M.values()]
    h.append(plt.Line2D([], [], color=BASE, lw=2, label="best trivial baseline"))
    fig.legend(handles=h, loc="lower center", ncol=4, frameon=False, bbox_to_anchor=(0.5, -0.08), fontsize=7)
    fig.tight_layout(w_pad=1.6)
    save(fig, "Fig2")


# ---------------------------------------------------------------- Fig 3: phenomena
def ci(s):
    lo, hi = ast.literal_eval(s)
    return lo, hi


def fig3():
    fig, axs = plt.subplots(1, 3, figsize=(7.5, 2.5), gridspec_kw={"width_ratios": [1.2, 1.2, 1.0]})
    # A: hidden paralogy trap
    tr = pd.read_csv(A / "phenomena_ygob_trap.tsv", sep="\t").set_index("method")
    ax = axs[0]
    for j, (k, (lab, col, mk)) in enumerate(M.items()):
        r = tr.loc[k]
        for x, kind in enumerate(["same_track", "pseudo"]):
            v = r[f"frac_correct_{kind}"]; lo, hi = ci(r[f"ci_{kind}"])
            xx = x + (j - 1) * 0.2
            ax.errorbar(xx, v, yerr=[[v - lo], [hi - v]], marker=mk, color=col, ms=4.5, lw=0.8, capsize=0,
                        mec=SURF, mew=0.5, label=lab if x == 0 else None)
    ax.set_xticks([0, 1]); ax.set_xticklabels(["true orthologs\n(n=1,189)", "pseudo-orthologs\n(n=49)"])
    ax.set_ylim(0, 1); ax.set_ylabel("Pairs classified correctly"); ax.grid(axis="y", color=GRID, lw=0.5)
    ax.set_title("A  Yeast WGD hidden paralogy", fontsize=8, loc="left")
    # B: WGD events, failure side (Y, testable events)
    ev = pd.read_csv(A / "phenomena_event_tags.tsv", sep="\t")
    ev = ev[(ev.panel == "Y") & (ev.event_tags == "WGD")].set_index("method")
    ax = axs[1]
    for j, (k, (lab, col, mk)) in enumerate(M.items()):
        r = ev.loc[k]
        for x, c in enumerate(["separation", "cohesion_min", "recovery_0.5"]):
            ax.plot(x + (j - 1) * 0.2, r[c], marker=mk, color=col, ms=4.5, ls="none", mec=SURF, mew=0.5)
    ax.set_xticks([0, 1, 2]); ax.set_xticklabels(["copies\napart", "each copy\ntogether", "event\nrecovered"])
    ax.set_ylim(0, 1); ax.grid(axis="y", color=GRID, lw=0.5)
    ax.set_title(f"B  YGOB WGD events (n={int(ev['n_events'].iloc[0])})", fontsize=8, loc="left")
    # C: Hox 2R / 3R recovery (V03)
    hx = pd.read_csv(A / "phenomena_V03_hox_events.tsv", sep="\t")
    ax = axs[2]
    for j, (k, (lab, col, mk)) in enumerate(M.items()):
        for x, lev in enumerate(["Gnathostomata", "Teleostei"]):
            r = hx[(hx.method == k) & (hx.level == lev)]
            if not r.empty:
                ax.plot(x + (j - 1) * 0.2, r["recovery_0.5"].iloc[0], marker=mk, color=col, ms=4.5, ls="none",
                        mec=SURF, mew=0.5)
    ax.set_xticks([0, 1]); ax.set_xticklabels(["2R\n(n=13)", "3R\n(n=11)"]); ax.set_ylim(0, 1)
    ax.grid(axis="y", color=GRID, lw=0.5); ax.set_title("C  Hox cluster duplications", fontsize=8, loc="left")
    axs[0].legend(frameon=False, fontsize=6.5, loc="center left")
    fig.tight_layout()
    save(fig, "Fig3")


# ---------------------------------------------------------------- Fig 4: robustness
def fig4():
    rb = pd.read_csv(A / "robustness.tsv", sep="\t")
    rb = rb[(rb.metric == "merge_accuracy_node")]
    c200 = ROOT / "results" / "robustness"
    rows = [("frag", "both", "Fragmentation\n(both pieces)"), ("frag", "either", "Fragmentation\n(either piece)"),
            ("contam", "", "Contamination\n(20 proteins)"), ("alttree", "", "Alternative\nspecies tree"),
            ("taxrm", "", "Taxon removal\n(2 species)")]
    fig, ax = plt.subplots(figsize=(4.6, 2.6))
    tools = {"orthofinder": M["orthofinder/default"], "fastoma": M["fastoma/default_declared"]}
    for i, (p, cr, lab) in enumerate(rows):
        for j, (tool, (tl, col, mk)) in enumerate(tools.items()):
            r = rb[(rb.pert == p) & (rb.tool == tool) & (rb.credit.fillna("") == cr)]
            if r.empty:
                continue
            d, lo, hi = r["delta"].iloc[0], r["ci_lo"].iloc[0], r["ci_hi"].iloc[0]
            yy = i + (j - 0.5) * 0.3
            ax.errorbar(d, yy, xerr=[[d - lo], [hi - d]], marker=mk, color=col, ms=4.5, lw=0.8, capsize=0,
                        mec=SURF, mew=0.5, label=tl if i == 0 else None)
    ax.axvline(0, color=INK2, lw=0.6)
    ax.set_yticks(range(len(rows))); ax.set_yticklabels([r[2] for r in rows]); ax.invert_yaxis()
    ax.set_xlabel("Δ merge-level accuracy (perturbed − clean), yeast panel")
    ax.grid(axis="x", color=GRID, lw=0.5); ax.legend(frameon=False, fontsize=7, loc="lower left")
    fig.tight_layout()
    save(fig, "Fig4")


if __name__ == "__main__":
    fig1(); fig2(); fig3(); fig4()
    print("figures written to", OUT)
