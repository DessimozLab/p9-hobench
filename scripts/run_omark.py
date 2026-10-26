#!/usr/bin/env python3
"""OMArk completeness/contamination QC of every panel proteome from the precomputed OMAmer hogmaps (SYNTHESIS D7;
covariates only, never used to filter the benchmark).

The FastOMA-compatible hogmaps are `omamer search -n 10` (up to 10 families per query); OMArk expects one placement
per protein, so the best family (the first row of each query, rows are ordered by decreasing family_p) is
extracted to data/omark/<CODE>/<CODE>.top1.hogmap first. `--check-top1 CODE` verifies on one species that this
equals a fresh `omamer search -n 1` (needs ~14 GB RAM).

OMArk 0.5.0 (pip wheel installed with --no-deps into conda_envs/omark_site; deps from conda_envs/ortho):
    omark -f <top1.hogmap> -d LUCA.h5 -t <taxid> -e data/raw/ncbi/taxa.sqlite -of <CODE>.fa -o data/omark/<CODE>
Summary -> panels/omark_summary.tsv (committed). Runs one OMArk at a time (~1.6 GB RSS, ~2.5 min each).
"""
from __future__ import annotations

import argparse
import csv
import os
import re
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
HOG = ROOT / "data" / "hogmaps"
OUT = ROOT / "data" / "omark"
SITE = ROOT / "conda_envs" / "omark_site"
PY = ROOT / "conda_envs" / "ortho" / "bin" / "python"
TAXDB = ROOT / "data" / "raw" / "ncbi" / "taxa.sqlite"


def species() -> dict:
    sp = {}
    for st in sorted((ROOT / "panels").glob("*/species.tsv")):
        for r in csv.DictReader(open(st), delimiter="\t"):
            sp.setdefault(r["code"], (r["taxid"], st.parent / "proteomes" / f"{r['code']}.fa", st.parent.name))
    return sp


def top1(src: Path, dst: Path) -> None:
    seen = set()
    with open(src) as f, open(dst, "w") as o:
        for ln in f:
            if ln.startswith("!") or ln.startswith("qseqid"):
                o.write(ln)
                continue
            q = ln.split("\t", 1)[0]
            if q not in seen:
                seen.add(q)
                o.write(ln)


def parse_sum(p: Path) -> dict:
    t = p.read_text()
    d = {}
    m = re.search(r"#The selected clade was (.+)", t)
    d["selected_clade"] = m.group(1).strip() if m else ""
    m = re.search(r"#Number of conserved HOGs is: (\d+)", t)
    d["n_conserved_hogs"] = m.group(1) if m else ""
    m = re.search(r"S:([\d.]+)%,D:([\d.]+)%\[U:([\d.]+)%,E:([\d.]+)%\],M:([\d.]+)%", t)
    if m:
        d.update(zip(["single_pct", "dup_pct", "dup_unexpected_pct", "dup_expected_pct", "missing_pct"], m.groups()))
    m = re.search(r"#On the whole proteome, there are (\d+) proteins", t)
    d["n_proteins_omark"] = m.group(1) if m else ""
    m = re.search(r"A:([\d.]+)%\[P:([\d.]+)%,F:([\d.]+)%\],I:([\d.]+)%\[P:[\d.]+%,F:[\d.]+%\],C:([\d.]+)%"
                  r"\[P:[\d.]+%,F:[\d.]+%\],U:([\d.]+)%", t)
    if m:
        d.update(zip(["consistent_pct", "partial_pct", "fragment_pct", "inconsistent_pct", "contamination_pct",
                      "unknown_pct"], m.groups()))
    return d


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--db", default="/tmp/hobench_omamer/LUCA.h5")
    ap.add_argument("--codes", nargs="*")
    a = ap.parse_args(argv)
    sp = species()
    env = dict(os.environ, PYTHONPATH=str(SITE), PATH=f"{PY.parent}:{os.environ['PATH']}")
    rows = []
    for code, (taxid, fa, panel) in sp.items():
        if a.codes and code not in a.codes:
            continue
        hm = HOG / f"{code}.fa.hogmap"
        od = OUT / code
        summ = od / f"{code}.top1.sum"
        if hm.exists() and not summ.exists():
            od.mkdir(parents=True, exist_ok=True)
            t1 = od / f"{code}.top1.hogmap"
            top1(hm, t1)
            t0 = time.time()
            r = subprocess.run([str(PY), str(SITE / "bin" / "omark"), "-f", str(t1), "-d", a.db, "-t", str(taxid),
                                "-e", str(TAXDB), "-of", str(fa), "-o", str(od)], env=env,
                               capture_output=True, text=True)
            (od / "omark.log").write_text(r.stdout + r.stderr)
            print(f"{code}: exit {r.returncode} {time.time() - t0:.0f}s", flush=True)
        row = {"code": code, "first_panel": panel, "taxid": taxid}
        row.update(parse_sum(summ) if summ.exists() else {"selected_clade": "NOT_RUN"})
        rows.append(row)
    if not a.codes:
        cols = ["code", "first_panel", "taxid", "selected_clade", "n_conserved_hogs", "single_pct", "dup_pct",
                "dup_unexpected_pct", "dup_expected_pct", "missing_pct", "n_proteins_omark", "consistent_pct",
                "partial_pct", "fragment_pct", "inconsistent_pct", "contamination_pct", "unknown_pct"]
        with open(ROOT / "panels" / "omark_summary.tsv", "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=cols, delimiter="\t", lineterminator="\n", restval="")
            w.writeheader()
            w.writerows(rows)
    return 0


if __name__ == "__main__":
    sys.exit(main())
