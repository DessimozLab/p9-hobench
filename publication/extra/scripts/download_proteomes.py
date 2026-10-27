#!/usr/bin/env python
"""Install the pinned panel proteomes (UniProt Reference Proteomes, release 2026_03) and link them into panels/.

Two sources, tried in this order:
  1. --archive PATH : the release archive hobench_proteomes_uniprot2026_03.tar.gz (exact files used in the paper).
  2. UniProt FTP via curation_tools.fetch_proteome (only reproduces the paper if UniProt still serves 2026_03).

Every FASTA is checked against the sha256 recorded in panels/<P>/species.tsv. Scoring refuses mismatching files,
because a different UniProt release changes protein sets and accessions.

Usage:
  python scripts/download_proteomes.py --archive hobench_proteomes_uniprot2026_03.tar.gz [--panels Y K]
  python scripts/download_proteomes.py --panels K            # try UniProt
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import subprocess
import sys
import tarfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PANELS = ["Y", "K", "E", "V", "P"]


def sha256(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def rows(panel: str):
    with open(ROOT / "panels" / panel / "species.tsv") as fh:
        yield from csv.DictReader(fh, delimiter="\t")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--panels", nargs="+", default=PANELS)
    ap.add_argument("--archive", default=None, help="release archive with data/proteomes/... (recommended)")
    a = ap.parse_args()
    if a.archive:
        with tarfile.open(a.archive) as t:
            t.extractall(ROOT, filter="data")
        print(f"unpacked {a.archive}")
    bad = 0
    for P in a.panels:
        pdir = ROOT / "panels" / P / "proteomes"
        pdir.mkdir(parents=True, exist_ok=True)
        for r in rows(P):
            fa = ROOT / r["fasta_path"]
            if not fa.exists():
                print(f"[{P}] {r['code']}: fetching {r['upid']} from UniProt ...", flush=True)
                subprocess.run([sys.executable, "-m", "curation_tools.fetch_proteome", r["upid"],
                                "--out", str(fa.parent.parent)], cwd=ROOT, check=False)
            if not fa.exists():
                print(f"[{P}] {r['code']}: MISSING {fa}"); bad += 1; continue
            if sha256(fa) != r["sha256"]:
                print(f"[{P}] {r['code']}: sha256 MISMATCH (not UniProt 2026_03?) {fa}"); bad += 1
            link = pdir / f"{r['code']}.fa"
            if link.is_symlink() or link.exists():
                link.unlink()
            link.symlink_to(Path("../../..") / r["fasta_path"])
    print("all proteomes verified" if not bad else f"{bad} problem(s): use --archive with the release archive")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
