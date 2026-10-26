#!/usr/bin/env python
"""Audit reference families for denylisted (tool-derived) provenance and for leakage of tool outputs.

Scans families/**/{family.yaml,collection.yaml,assertions.tsv,build.py,curation_log.md} for
(1) denylisted sources used as truth (CURATION_PROTOCOL §0.3) and (2) references to runs/ or results/ paths.
A denylisted name that appears only in a sentence that rejects/excludes it is reported as a NOTE, not a violation;
every hit is printed with file:line so a human can judge. Exit code 1 if any VIOLATION.

Usage: python scripts/audit_sources.py [--families families] [--out families/AUDIT.md]
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

DENY = {
    "OMA": r"\bOMA\b(?!mer)|omabrowser",
    "OrthoFinder": r"OrthoFinder",
    "Ensembl Compara/HCOP": r"Compara|\bHCOP\b",
    "VGNC": r"\bVGNC\b",
    "OHNOLOGS": r"OHNOLOGS|ohnologs\.curie",
    "PANTHER": r"\bPANTHER\b",
    "TreeFam": r"TreeFam",
    "eggNOG": r"eggNOG",
    "OrthoDB": r"OrthoDB",
    "PLAZA": r"\bPLAZA\b",
    "LOC names": r"\bLOC\d{5,}\b",
    "KRAB age_MA": r"age_MA",
}
LEAK = r"(?<![\w/])(runs|results)/"
NEGATION = re.compile(r"\b(not|never|no|exclud\w*|avoid\w*|denylist\w*|reject\w*|instead of|rather than|"
                      r"cannot|can't|unusable|lack\w*|is not used|not used)\b", re.I)
FILES = ("family.yaml", "collection.yaml", "assertions.tsv", "build.py", "curation_log.md")


def scan(root: Path):
    hits = []
    for f in sorted(root.rglob("*")):
        if f.name not in FILES or "verification" in f.parts:
            continue
        for i, line in enumerate(f.read_text(errors="replace").splitlines(), 1):
            for name, pat in DENY.items():
                if re.search(pat, line):
                    kind = "NOTE" if NEGATION.search(line) or f.name == "curation_log.md" else "VIOLATION"
                    hits.append((kind, name, f, i, line.strip()[:160]))
            if re.search(LEAK, line) and f.name != "curation_log.md":
                hits.append(("VIOLATION", "tool-output path", f, i, line.strip()[:160]))
            elif re.search(LEAK, line):
                hits.append(("NOTE", "tool-output path (log)", f, i, line.strip()[:160]))
    return hits


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--families", default="families")
    ap.add_argument("--out", default="families/AUDIT.md")
    a = ap.parse_args(argv)
    hits = scan(Path(a.families))
    nv = sum(h[0] == "VIOLATION" for h in hits)
    lines = ["# Source audit (scripts/audit_sources.py)", "",
             f"VIOLATIONS: {nv}; NOTES: {len(hits) - nv}. Every hit listed for human judgement.", "",
             "| kind | rule | file:line | text |", "|---|---|---|---|"]
    lines += [f"| {k} | {n} | {f}:{i} | {t.replace('|', '/')} |" for k, n, f, i, t in hits]
    Path(a.out).write_text("\n".join(lines) + "\n")
    print(f"{nv} violations, {len(hits) - nv} notes -> {a.out}")
    return 1 if nv else 0


if __name__ == "__main__":
    sys.exit(main())
