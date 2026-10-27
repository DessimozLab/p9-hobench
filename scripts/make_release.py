#!/usr/bin/env python
"""Build release/MANIFEST.tsv (sha256 + size of every file in the HOBench-core data release) and
release/file_list.txt (input to the tarball command in release/README.md).

    python scripts/make_release.py            # run scripts/export_references.py first

Covers: full + subset input FASTAs of every panel (incl. robustness variants Y_*), subset member tables,
species trees (main + alternative), levels.tsv, species.tsv, perturbation tables, and the reference export.
FASTAs are NOT committed to git (see .gitignore); MANIFEST.tsv and file_list.txt are.
"""
from __future__ import annotations
import hashlib, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
REL = ROOT / "release"


def sha256(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def kind(p: Path) -> str:
    r = p.relative_to(ROOT).as_posix()
    n = p.name
    if r.startswith("release/references/"):
        return "reference_export"
    if "/subset/" in r:
        return "fasta_subset" if n.endswith((".fa", ".fasta")) else "subset_table"
    if "/proteomes/" in r:
        return "fasta_full"
    if n.startswith("species_tree"):
        return "species_tree"
    return {"levels.tsv": "levels", "species.tsv": "species_table"}.get(n, "panel_table")


def main():
    files = []
    for pd in sorted(d for d in (ROOT / "panels").iterdir() if d.is_dir()):
        for p in sorted(pd.rglob("*")):
            if not p.is_file() or "/raw/" in p.as_posix() or p.name in ("manifest.json", "header_map.tsv"):
                continue
            if p.suffix in (".fa", ".fasta", ".nwk") or p.name in ("levels.tsv", "species.tsv", "members.tsv",
                                                                  "perturbation.json", "subset_seeds.txt"):
                files.append(p)
    files += sorted((REL / "references").rglob("*.gz"))
    rows = []
    for p in files:
        rel = p.relative_to(ROOT).as_posix()
        panel = rel.split("/")[1] if rel.startswith(("panels/", "release/references/")) else ""
        rows.append((kind(p), panel, rel, p.stat().st_size, sha256(p)))
        print(rows[-1][0], rel, file=sys.stderr) if len(rows) % 50 == 0 else None
    with open(REL / "MANIFEST.tsv", "w") as fh:
        fh.write("kind\tpanel\tpath\tbytes\tsha256\n")
        for r in rows:
            fh.write("\t".join(map(str, r)) + "\n")
    (REL / "file_list.txt").write_text("".join(r[2] + "\n" for r in rows) + "release/MANIFEST.tsv\nDATA.md\nLICENSE\n")
    tot = {}
    for r in rows:
        tot[r[0]] = tot.get(r[0], 0) + r[3]
    print(len(rows), "files;", {k: f"{v/1e6:.0f} MB" for k, v in tot.items()})


if __name__ == "__main__":
    main()
