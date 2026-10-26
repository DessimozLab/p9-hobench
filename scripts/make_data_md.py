#!/usr/bin/env python3
"""Generate DATA.md (one row per dataset: source, version, retrieval date, checksum, licence, redistributed in git,
generating script) from the manifests on disk. Nothing is typed by hand except the constant facts in KNOWN below.

Inputs
  data/proteomes/*/manifest.json            (curation_tools.fetch_proteome)
  panels/*/species.tsv                       (scripts/build_panels.py) -> panel membership
  data/raw/omamer/LUCA.h5(.provenance.json)  (written here with --refresh-omamer: URL redirect, Zenodo record
                                              metadata, sha256 + md5 of the local file vs Zenodo md5)
  taxdump.tar.gz                             (NCBI taxonomy dump; sha256 + file mtime)
  literature/raw/C/*.newick.json             (per-file provenance sidecars)
  literature/papers/*/meta.json (+ MANIFEST.sha256), literature/raw/B2/MANIFEST.tsv, literature/raw/*/*

Usage: python scripts/make_data_md.py [--refresh-omamer]
"""
from __future__ import annotations

import argparse
import csv
import glob
import hashlib
import json
import os
import subprocess
import sys
import time
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OMAMER = ROOT / "data" / "raw" / "omamer" / "LUCA.h5"
OMAMER_PROV = OMAMER.with_name("LUCA.h5.provenance.json")
OMAMER_URL = "https://omabrowser.org/All/LUCA.h5"
TAXDUMP = ROOT / "taxdump.tar.gz"
TAXDUMP_PROV = ROOT / "data" / "raw" / "ncbi" / "taxdump.provenance.json"


# Source URLs of raw files not produced by a fetch script (verified by re-download + sha256 match).
RAW_SOURCES = {
    "literature/raw/P_shared/all.locus_brief_info.7.0":
        "http://rice.uga.edu/pub/data/Eukaryotic_Projects/o_sativa/annotation_dbs/pseudomolecules/version_7.0/all.dir/"
        "all.locus_brief_info.7.0 (MSU Rice Genome Annotation Project r7; sha256 c8c546993cf583e8...; re-downloaded "
        "2026-10-27 and matched; public research data, not redistributed)",
}


def hashes(p: Path, algos=("sha256",)) -> dict:
    hs = {a: hashlib.new(a) for a in algos}
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 22), b""):
            for h in hs.values():
                h.update(b)
    return {a: h.hexdigest() for a, h in hs.items()}


def git_tracked(path: Path) -> int:
    out = subprocess.run(["git", "ls-files", str(path.relative_to(ROOT))], cwd=ROOT, capture_output=True, text=True)
    return len([x for x in out.stdout.splitlines() if x])


def git_ignored(path: Path) -> bool:
    return subprocess.run(["git", "check-ignore", "-q", str(path.relative_to(ROOT))], cwd=ROOT).returncode == 0


def refresh_omamer() -> None:
    import requests
    location = requests.head(OMAMER_URL, allow_redirects=False, timeout=60).headers.get("Location") or OMAMER_URL
    rec_id = location.split("/records/")[1].split("/")[0] if "/records/" in (location or "") else None
    rec = json.load(urllib.request.urlopen(f"https://zenodo.org/api/records/{rec_id}", timeout=60)) if rec_id else {}
    zf = next((f for f in rec.get("files", []) if f["key"] == "LUCA.h5"), {})
    h = hashes(OMAMER, ("sha256", "md5"))
    st = OMAMER.stat()
    prov = {"dataset": "OMAmer precomputed database LUCA.h5", "url": OMAMER_URL, "redirects_to": location,
            "zenodo_record": rec_id, "zenodo_doi": rec.get("doi"), "zenodo_title": rec.get("metadata", {}).get("title"),
            "zenodo_publication_date": rec.get("metadata", {}).get("publication_date"),
            "license": (rec.get("metadata", {}).get("license") or {}).get("id"),
            "zenodo_size": zf.get("size"), "zenodo_checksum": zf.get("checksum"),
            "local_size": st.st_size, "local_mtime": time.strftime("%Y-%m-%dT%H:%M:%S", time.localtime(st.st_mtime)),
            "sha256": h["sha256"], "md5": h["md5"],
            "md5_matches_zenodo": zf.get("checksum") == f"md5:{h['md5']}" if zf else None,
            "checked": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}
    OMAMER_PROV.write_text(json.dumps(prov, indent=1) + "\n")
    if TAXDUMP.exists():
        st = TAXDUMP.stat()
        TAXDUMP_PROV.parent.mkdir(parents=True, exist_ok=True)
        TAXDUMP_PROV.write_text(json.dumps({
            "dataset": "NCBI Taxonomy dump", "url": "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz",
            "local_mtime": time.strftime("%Y-%m-%dT%H:%M:%S", time.localtime(st.st_mtime)), "bytes": st.st_size,
            "sha256": hashes(TAXDUMP)["sha256"], "license": "public domain (NCBI; US government work)",
            "checked": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}, indent=1) + "\n")


def panel_membership() -> dict:
    mem = {}
    for st in sorted(glob.glob(str(ROOT / "panels" / "*" / "species.tsv"))):
        p = Path(st).parent.name
        for r in csv.DictReader(open(st), delimiter="\t"):
            mem.setdefault(r["upid"], []).append(f"{p}:{r['code']}")
    return mem


def dir_digest(d: Path) -> tuple[int, int, str]:
    """(n files, bytes, sha256 over sorted 'sha256  relpath' lines)."""
    files = sorted(x for x in d.rglob("*") if x.is_file()) if d.is_dir() else [d]
    lines, nb = [], 0
    for f in files:
        nb += f.stat().st_size
        lines.append(f"{hashes(f)['sha256']}  {f.relative_to(d) if d.is_dir() else f.name}")
    return len(files), nb, hashlib.sha256("\n".join(lines).encode()).hexdigest()


def md_row(cells) -> str:
    return "| " + " | ".join(str(c).replace("|", "/").replace("\n", " ") for c in cells) + " |"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--refresh-omamer", action="store_true")
    a = ap.parse_args(argv)
    if a.refresh_omamer or not OMAMER_PROV.exists():
        refresh_omamer()
    out = ["# DATA.md - datasets used by the benchmark", "",
           f"Generated by `scripts/make_data_md.py` on {time.strftime('%Y-%m-%d')} from the manifests on disk; "
           "do not edit by hand. *Redistributed* = the file itself is committed to git (y), or only its checksum/"
           "manifest is (n). Licences: UniProt CC BY 4.0; OMA/OMAmer CC BY 4.0 (Zenodo record); NCBI Taxonomy "
           "public domain; literature files carry the publisher's terms (open-access licence where recorded).", ""]
    # --- proteomes
    mem = panel_membership()
    out += ["## 1. Proteomes (UniProt, pinned release 2026_03)", "",
            "Script: `python -m curation_tools.fetch_proteome <UPID> --out data/proteomes` (fails on a release "
            "mismatch). Files under `data/proteomes/` are git-ignored (redistributed: n); `panels/<P>/species.tsv` "
            "records the sha256.", "",
            md_row(["UPID", "organism", "taxid", "type", "panels", "source", "release", "retrieved", "n_proteins",
                    "proteome.fasta sha256", "licence", "redistributed"]),
            md_row(["---"] * 12)]
    n = 0
    for mf in sorted(glob.glob(str(ROOT / "data" / "proteomes" / "*" / "manifest.json"))):
        m = json.load(open(mf))
        src = (m["downloads"][0]["url"] if m["downloads"] else "")
        src_txt = {"ftp": "UniProt FTP reference_proteomes", "rest": "UniProtKB REST stream",
                   "uniparc": "UniParc REST stream"}.get(m["source"], m["source"])
        out.append(md_row([m["proteome_id"], m["organism"] or "", m["taxid"], m["proteome_type"],
                           ",".join(mem.get(m["proteome_id"], [])) or "(not in a panel)",
                           f"[{src_txt}]({src})", m["uniprot_release"], m["retrieval_date"][:10], m["n_proteins"],
                           f"`{m['outputs']['proteome_fasta']['sha256']}`", m["license"], "n"]))
        n += 1
    out += ["", f"{n} proteomes on disk.", ""]
    # --- reference databases
    o = json.load(open(OMAMER_PROV))
    out += ["## 2. Reference databases and trees", "",
            md_row(["dataset", "source", "version", "retrieved/checked", "size", "sha256", "licence", "redistributed",
                    "script"]), md_row(["---"] * 9),
            md_row(["OMAmer DB `data/raw/omamer/LUCA.h5`", f"{o['url']} -> {o['redirects_to']} (Zenodo doi:"
                    f"{o['zenodo_doi']})", f"{o['zenodo_title']} ({o['zenodo_publication_date']}); md5 matches "
                    f"Zenodo: {o['md5_matches_zenodo']}", f"file mtime {o['local_mtime']}; checked {o['checked'][:10]}",
                    f"{o['local_size']:,} B", f"`{o['sha256']}`", o["license"], "n",
                    "manual download; provenance by scripts/make_data_md.py --refresh-omamer"])]
    if TAXDUMP_PROV.exists():
        t = json.load(open(TAXDUMP_PROV))
        out.append(md_row(["NCBI Taxonomy `taxdump.tar.gz` (-> data/raw/ncbi/taxa.sqlite via ete3)", t["url"],
                           "snapshot of file mtime", t["local_mtime"][:10], f"{t['bytes']:,} B", f"`{t['sha256']}`",
                           t["license"], "n", "scripts/build_panels.py (NCBI level annotation)"]))
    for sc in sorted(glob.glob(str(ROOT / "literature" / "raw" / "C" / "*.json"))):
        j = json.load(open(sc))
        if "sha256" not in j or "license" not in j:
            continue
        f = Path(sc[:-5])
        out.append(md_row([f"{j['dataset']} `{f.relative_to(ROOT)}`", f"{j['source_url']} (doi:{j.get('figshare_doi')})",
                           f"paper doi:{j.get('paper_doi')}", j["retrieval_date"], f"{f.stat().st_size:,} B",
                           f"`{j['sha256']}`", j["license"], "y" if git_tracked(f) else "n", j.get("used_by", "")]))
    out += ["", "Derived species panels (`panels/<P>/`) are generated by `scripts/build_panels.py`; citation DOIs "
            "verified in Europe PMC are in `panels/citations.tsv`. OMAmer placements `data/hogmaps/<CODE>.fa.hogmap` "
            "(git-ignored) are generated by `scripts/precompute_hogmaps.py` from LUCA.h5 above; OMArk 0.5.0 QC "
            "(`scripts/run_omark.py`, outputs in git-ignored `data/omark/`, summary `panels/omark_summary.tsv`).", ""]
    # --- literature papers
    out += ["## 3. Literature: per-paper folders (`literature/papers/*`)", "",
            "Raw files (`literature/papers/*/raw/`) are git-ignored; metadata, derived tables and MANIFEST.sha256 are "
            "committed.", "",
            md_row(["folder", "id", "title", "year", "licence (EPMC)", "redistributable", "raw files", "raw sha256 "
                    "(digest)", "redistributed raw", "script"]), md_row(["---"] * 10)]
    for d in sorted(p for p in (ROOT / "literature" / "papers").iterdir() if p.is_dir()):
        mj = d / "meta.json"
        m = json.load(open(mj)) if mj.exists() else {}
        raw = d / "raw"
        nf, nb, dg = dir_digest(raw) if raw.exists() else (0, 0, "")
        pid = m.get("doi") or m.get("pmid") or ""
        out.append(md_row([d.name, f"doi:{pid}" if m.get("doi") else pid, (m.get("title") or "")[:90],
                           m.get("year", ""), m.get("license") or "not recorded", m.get("redistributable", ""),
                           f"{nf} ({nb:,} B)", f"`{dg[:16]}...`" if dg else "",
                           "n" if raw.exists() and git_ignored(raw) else ("y" if raw.exists() and git_tracked(raw)
                                                                            else "-"),
                           "curation_tools.fetch_paper" + (" (+MANIFEST.sha256)" if (d / "MANIFEST.sha256").exists()
                                                           else "")]))
    # --- literature raw
    b2 = {}
    mf = ROOT / "literature" / "raw" / "B2" / "MANIFEST.tsv"
    if mf.exists():
        for ln in open(mf):
            if not ln.startswith("#") and ln.strip():
                md5, nbytes, path, src = ln.rstrip("\n").split("\t")[:4]
                b2[path.split("/")[0]] = src.rsplit("/", 2)[0] + "/..."
    out += ["", "## 4. Literature: survey raw files (`literature/raw/<part>/*`)", "",
            "One row per top-level entry; checksum = sha256 over the sorted per-file sha256 list. Licence is that of "
            "the source (Europe PMC OA subset, publisher, or database terms) and was not verified per file unless "
            "stated; entries that are git-ignored are not redistributed.", "",
            md_row(["entry", "files", "bytes", "digest sha256", "in git", "redistributed", "source"]),
            md_row(["---"] * 7)]
    for part in sorted((ROOT / "literature" / "raw").iterdir()):
        for e in sorted(part.iterdir()):
            nf, nb, dg = dir_digest(e)
            tr = git_tracked(e)
            ign = git_ignored(e)
            src = b2.get(e.name, "") if part.name == "B2" else RAW_SOURCES.get(str(e.relative_to(ROOT)), "")
            out.append(md_row([str(e.relative_to(ROOT)), nf, f"{nb:,}", f"`{dg[:16]}...`",
                               f"{tr}/{nf} tracked" if tr else ("ignored" if ign else "untracked"),
                               "y" if tr else "n", src]))
    (ROOT / "DATA.md").write_text("\n".join(out) + "\n")
    print(f"wrote DATA.md ({len(out)} lines)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
