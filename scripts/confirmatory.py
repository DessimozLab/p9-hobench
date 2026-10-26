#!/usr/bin/env python
"""Confirmatory test (R2-7): OrthoFinder-default vs FastOMA-default(declared), one paired sign-flip
permutation test (10,000 flips, seed hobench.seeds.SEEDS[0]) per primary endpoint over all pooled
headline items across the panels' results/<P>_<mode>/per_family.tsv; Holm over the 4 endpoints.
Items: tracks family + collection only, headline-eligible; a collection is one item (mean of its
units); items with an undefined value for either method are dropped.
Usage: python scripts/confirmatory.py [--mode full] [--panels Y K]"""
import argparse, csv, sys
from collections import defaultdict
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from hobench.metrics import PRIMARY      # noqa: E402
from hobench.seeds import SEEDS          # noqa: E402
from hobench.stats import holm           # noqa: E402

A, B = ("orthofinder", "default"), ("fastoma", "default_declared")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", default="full")
    ap.add_argument("--panels", nargs="+", default=["Y", "K"])
    ap.add_argument("--n-perm", type=int, default=10000)
    ap.add_argument("--out", default=str(ROOT / "results" / "confirmatory.tsv"))
    a = ap.parse_args()
    vals = {m: defaultdict(lambda: defaultdict(list)) for m in PRIMARY}   # metric -> item -> method -> [unit values]
    for P in a.panels:
        with open(ROOT / "results" / f"{P}_{a.mode}" / "per_family.tsv") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                if r["metric"] not in PRIMARY or r["track"] == "control" or r["headline"] != "1":
                    continue
                key = (r["method"], r["variant"])
                if key in (A, B):
                    item = f"{P}:collection:{r['collection']}" if r["track"] == "collection" else f"{P}:{r['family']}"
                    vals[r["metric"]][item][key].append(float(r["value"]))
    rng = np.random.default_rng(SEEDS[0])
    res = {}
    for m in PRIMARY:
        d = np.array([np.mean(v[A]) - np.mean(v[B]) for v in vals[m].values() if v[A] and v[B]])
        n = len(d)
        if n == 0:
            res[m] = (0, None, None); continue
        flips = rng.choice([-1.0, 1.0], size=(a.n_perm, n))
        obs = abs(d.mean())
        p = (1 + np.sum(np.abs((flips * d).mean(axis=1)) >= obs - 1e-12)) / (a.n_perm + 1)
        res[m] = (n, float(d.mean()), float(p))
    hp = holm({m: r[2] for m, r in res.items()})
    with open(a.out, "w") as fh:
        fh.write("metric\tn_items\tmean_diff_OF_minus_FastOMA\tp_perm\tp_holm\n")
        for m in PRIMARY:
            n, md, p = res[m]
            fh.write(f"{m}\t{n}\t{md}\t{p}\t{hp[m]}\n")
            print(f"{m:22s} n={n} diff={md} p={p} holm={hp[m]}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
