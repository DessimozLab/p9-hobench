#!/usr/bin/env python3
"""Precompute OMAmer placements (hogmaps) once per species for all panels (SYNTHESIS E4), reusable by FastOMA.

Command replicates FastOMA's process omamer_run (tools/FastOMA/FastOMA.nf):
    omamer search -n 10 --db <omamer_db> --query <CODE>.fa --out <CODE>.fa.hogmap
plus `-t THREADS` (threads only parallelise the search; placements are unchanged). The query is the panel
symlink panels/<P>/proteomes/<CODE>.fa so the hogmap name is exactly what FastOMA looks for with --hogmap_in.

* species are taken from panels/<P>/species.tsv in the order given by --panels (Y K V E P by default),
  de-duplicated by code (the same code must point to the same FASTA sha256 in every panel);
* existing data/hogmaps/<CODE>.fa.hogmap whose recorded FASTA sha256 matches are skipped;
* at most --jobs omamer processes at a time (each holds the LUCA DB in RAM, ~14 GB RSS);
* the DB is read from a local-NVMe copy (--local-db) after checking its sha256 against the original;
* outputs are written to a temp name and renamed when omamer exits 0; runtime/RSS logged to
  data/hogmaps/hogmap_runs.tsv. Panels whose species.tsv does not exist yet are re-scanned between jobs.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DB = ROOT / "data" / "raw" / "omamer" / "LUCA.h5"
OUT = ROOT / "data" / "hogmaps"


def sha256_file(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 22), b""):
            h.update(b)
    return h.hexdigest()


def todo(panels: list[str]) -> list[tuple[str, Path, str]]:
    seen: dict[str, str] = {}
    out = []
    for p in panels:
        st = ROOT / "panels" / p / "species.tsv"
        if not st.exists():
            continue
        for r in csv.DictReader(open(st), delimiter="\t"):
            c = r["code"]
            if c in seen:
                if seen[c] != r["sha256"]:
                    raise SystemExit(f"ERROR: code {c} has different FASTAs in different panels")
                continue
            seen[c] = r["sha256"]
            out.append((c, ROOT / "panels" / p / "proteomes" / f"{c}.fa", r["sha256"]))
    return out


def done_ok(code: str, sha: str) -> bool:
    hm, meta = OUT / f"{code}.fa.hogmap", OUT / f"{code}.fa.hogmap.sha256"
    return hm.exists() and meta.exists() and meta.read_text().split()[0] == sha


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--panels", nargs="+", default=["Y", "K", "V", "E", "P"])
    ap.add_argument("--jobs", type=int, default=2)
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--local-db", default="/tmp/hobench_omamer/LUCA.h5")
    ap.add_argument("--db-sha256", default="", help="expected sha256 of LUCA.h5 (computed if empty)")
    a = ap.parse_args(argv)
    OUT.mkdir(parents=True, exist_ok=True)
    db = Path(a.local_db)
    ref_sha = a.db_sha256 or sha256_file(DB)
    (OUT / "LUCA.h5.sha256").write_text(f"{ref_sha}  {DB.relative_to(ROOT)}\n")
    if not db.exists() or db.stat().st_size != DB.stat().st_size:
        db.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(DB, db)
    if sha256_file(db) != ref_sha:
        raise SystemExit("ERROR: local DB copy sha256 differs from data/raw/omamer/LUCA.h5")
    log = OUT / "hogmap_runs.tsv"
    if not log.exists():
        log.write_text("code\tquery\tfasta_sha256\tcommand\tstart\twall_s\tmax_rss_kb\texit\n")
    running: dict[str, tuple[subprocess.Popen, float, Path, str, str]] = {}
    failed: set[str] = set()
    while True:
        pending = [(c, f, s) for c, f, s in todo(a.panels)
                   if c not in running and c not in failed and not done_ok(c, s)]
        while pending and len(running) < a.jobs:
            c, f, s = pending.pop(0)
            tmp = OUT / f".{c}.fa.hogmap.tmp"
            cmd = ["omamer", "search", "-n", "10", "--db", str(db), "--query", str(f), "--out", str(tmp),
                   "-t", str(a.threads)]
            timef = OUT / f".{c}.time"
            p = subprocess.Popen(["/usr/bin/time", "-f", "%M", "-o", str(timef)] + cmd,
                                 stdout=open(OUT / f"{c}.omamer.log", "w"), stderr=subprocess.STDOUT)
            running[c] = (p, time.time(), tmp, s, " ".join(cmd).replace(str(db), "LUCA.h5"))
            print(f"[{time.strftime('%H:%M:%S')}] start {c} ({len(running)} running)", flush=True)
        if not running:
            break
        time.sleep(10)
        for c in list(running):
            p, t0, tmp, s, cmd = running[c]
            if p.poll() is None:
                continue
            wall = time.time() - t0
            timef = OUT / f".{c}.time"
            rss = timef.read_text().strip().splitlines()[-1] if timef.exists() else ""
            if p.returncode == 0 and tmp.exists():
                os.replace(tmp, OUT / f"{c}.fa.hogmap")
                (OUT / f"{c}.fa.hogmap.sha256").write_text(f"{s}  {c}.fa\n")
            else:
                failed.add(c)
            with open(log, "a") as fh:
                fh.write(f"{c}\t{c}.fa\t{s}\t{cmd}\t{time.strftime('%Y-%m-%dT%H:%M:%S', time.localtime(t0))}"
                         f"\t{wall:.0f}\t{rss}\t{p.returncode}\n")
            timef.unlink(missing_ok=True)
            print(f"[{time.strftime('%H:%M:%S')}] done {c} exit={p.returncode} {wall:.0f}s rss={rss}kB", flush=True)
            del running[c]
    print(f"finished; failed: {sorted(failed)}")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
