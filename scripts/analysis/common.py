"""Shared helpers for scripts/analysis/*.py (read-only use of results/<P>_<mode>/ and hobench)."""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / "scripts"))
from hobench.metrics import PRIMARY   # noqa: E402,F401
from hobench.seeds import SEEDS       # noqa: E402,F401

OUT = ROOT / "results" / "analysis"
OUT.mkdir(parents=True, exist_ok=True)
# main analysis mode per panel (DESIGN sec. 7: Y,K full; E,P,V subset; modes never pooled within a panel)
SPEC = [("Y", "full"), ("K", "full"), ("E", "subset"), ("P", "subset"), ("V", "subset")]
METHOD_ORDER = [("orthofinder", "default"), ("orthofinder", "dendroblast"), ("fastoma", "default_declared"),
                ("fastoma", "default_late"), ("of_flat", "default"), ("of_flat", "dendroblast"),
                ("all_singletons", ""), ("one_big_hog", ""), ("random_relabel", "mean10"),
                ("positive_control", "reference")]
TOOLS2 = [("orthofinder", "default"), ("fastoma", "default_declared")]
# Okabe-Ito colour-blind-safe palette
OKABE = ["#0072B2", "#56B4E9", "#D55E00", "#E69F00", "#009E73", "#CC79A7", "#000000", "#999999", "#F0E442"]


def available(spec=SPEC):
    return [(P, m) for P, m in spec if (ROOT / "results" / f"{P}_{m}" / "per_family.tsv").exists()]


def read_tsv(P, mode, name, **kw):
    return pd.read_csv(ROOT / "results" / f"{P}_{mode}" / f"{name}.tsv", sep="\t", keep_default_na=False,
                       na_values=["", "NA", "None", "nan"], dtype={"variant": str, "collection": str}, **kw)


def item_of(df):
    """Item id: a collection counts as one item (its units are averaged), else the family."""
    coll = df["collection"].fillna("")
    istrue = (df["track"] == "collection") if "track" in df else (coll != "")
    return np.where(istrue & (coll != ""), "collection:" + coll, df["family"])


def items(P, mode, headline_only=True):
    """per-item means of every per_family metric: columns panel, mode, item, track, method, variant, metric, value, n_units."""
    d = read_tsv(P, mode, "per_family")
    d["variant"] = d["variant"].fillna("")
    if headline_only:
        d = d[(d["track"] != "control") & (d["headline"] == 1)]
    d = d.dropna(subset=["value"]).copy()
    d["item"] = item_of(d)
    g = d.groupby(["item", "track", "method", "variant", "metric"])["value"].agg(["mean", "size"]).reset_index()
    g = g.rename(columns={"mean": "value", "size": "n_units"})
    g.insert(0, "mode", mode); g.insert(0, "panel", P)
    return g


def fmt(x, nd=3):
    return "NA" if x is None or (isinstance(x, float) and not np.isfinite(x)) else f"{x:.{nd}f}"


def mlabel(m, v):
    return f"{m}" + (f"/{v}" if v else "")


def write_md(path, lines):
    Path(path).write_text("\n".join(lines) + "\n")
    print("wrote", path)
