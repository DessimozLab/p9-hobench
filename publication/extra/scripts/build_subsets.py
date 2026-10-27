#!/usr/bin/env python
"""Rebuild the frozen subset proteomes (panels/<P>/subset/<CODE>.fa) from the committed ID lists.

The subsets used in the paper (E, V, P; Y and K for validation) were built once with wrapper.subset (DIAMOND
neighbourhood of literature members + capped domain-only hits + seeded background) and frozen as ID lists in
panels/<P>/subset/ids/<CODE>.txt. This script only extracts those IDs from the verified full proteomes, so the subsets
are identical to the paper's without re-running DIAMOND.

Usage: python scripts/build_subsets.py [--panels E V P]
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--panels", nargs="+", default=["Y", "K", "E", "V", "P"])
    a = ap.parse_args()
    bad = 0
    for P in a.panels:
        sdir = ROOT / "panels" / P / "subset"
        for r in csv.DictReader(open(ROOT / "panels" / P / "species.tsv"), delimiter="\t"):
            code = r["code"]
            idf = sdir / "ids" / f"{code}.txt"
            if not idf.exists():
                continue
            want = set(idf.read_text().split())
            src = ROOT / "panels" / P / "proteomes" / f"{code}.fa"
            n, keep = 0, False
            with open(src) as fin, open(sdir / f"{code}.fa", "w") as fout:
                for line in fin:
                    if line.startswith(">"):
                        keep = line[1:].split()[0] in want
                        n += keep
                    if keep:
                        fout.write(line)
            if n != len(want):
                print(f"[{P}] {code}: {n}/{len(want)} IDs found"); bad += 1
        print(f"[{P}] subset written to {sdir}")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
