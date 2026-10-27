#!/usr/bin/env python
"""Main table: per panel (analysis mode) x method, primary endpoints (headline track, items; a collection = one item)
with n and baseline-normalised skill; plus pooled headline means over all items of all available panels.
Writes results/analysis/main_table.{md,tsv}."""
import pandas as pd
from common import OUT, PRIMARY, METHOD_ORDER, available, items, fmt, mlabel, write_md

SK = {m: f"skill_{m}" for m in PRIMARY}


def table(df):
    rows = []
    for (m, v) in METHOD_ORDER:
        s = df[(df.method == m) & (df.variant == v)]
        if s.empty:
            continue
        r = {"method": m, "variant": v, "panels": "".join(sorted(s.panel.unique()))}
        for e in PRIMARY:
            x = s[s.metric == e]["value"]; k = s[s.metric == SK[e]]["value"]
            r[e], r[f"n_{e}"] = (x.mean() if len(x) else None), len(x)
            r[SK[e]] = k.mean() if len(k) else None
        rows.append(r)
    return rows


def md(rows):
    h = "| method | panels | " + " | ".join(f"{e} (n) | skill" for e in PRIMARY) + " |"
    out = [h, "|" + "---|" * (2 + 2 * len(PRIMARY))]
    for r in rows:
        out.append(f"| {mlabel(r['method'], r['variant'])} | {r['panels']} | " + " | ".join(
            f"{fmt(r[e])} ({r[f'n_{e}']}) | {fmt(r[SK[e]])}" for e in PRIMARY) + " |")
    return out


def main():
    av = available()
    alld, tsv = [], []
    L = ["# HOBench main table (headline track)", "",
         "Source: results/<P>_<mode>/per_family.tsv (headline = track family+collection, headline-eligible; controls K01/K02/E05 excluded).",
         "Value = mean over items (a collection = one item = mean of its units); n = items with a defined value; "
         "skill = (s-b)/(1-b), b = per-family max of trivial baselines (mean over items). NA = no item defines the endpoint. "
         "Descriptive only; no CIs (n<5 almost everywhere). Panels: " + ", ".join(f"{P} {m}" for P, m in av), ""]
    for P, mode in av:
        d = items(P, mode); alld.append(d)
        rows = table(d)
        L += [f"## {P} ({mode})", ""] + md(rows) + [""]
        tsv += [dict(panel=P, mode=mode, **r) for r in rows]
    pooled = pd.concat(alld)
    rows = table(pooled)
    L += ["## Pooled headline means (all panels above, each in its analysis mode; items pooled, unweighted)", ""] + md(rows) + [""]
    L += ["Note: modes are never pooled within a panel; pooling across panels mixes full (Y,K) and subset (E,P,V) modes "
          "and is descriptive only (confirmatory test: scripts/confirmatory.py)."]
    tsv += [dict(panel="pooled", mode="mixed", **r) for r in rows]
    pd.DataFrame(tsv).to_csv(OUT / "main_table.tsv", sep="\t", index=False, float_format="%.4f")
    write_md(OUT / "main_table.md", L)


if __name__ == "__main__":
    main()
