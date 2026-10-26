"""Download a UniProt proteome as a canonical, one-protein-per-gene FASTA with a reproducibility manifest.

CLI:
  python -m curation_tools.fetch_proteome UP000002311 [--out data/proteomes] [--source auto|ftp|rest]
                                          [--include-isoforms] [--markdown]

Choice of source (documented in the manifest as "source"):
  ftp  (default for Reference/Representative proteomes)
       .../knowledgebase/reference_proteomes/<Kingdom>/<UP>/<UP>_<taxid>.fasta.gz
       This is UniProt's official canonical set: exactly ONE sequence per gene (geneCount), versioned by the
       UniProt release that created it, and it ships a gene2acc table. Example: UP000002311 has 6066 sequences
       in the FTP canonical file = geneCount, while the REST query (proteome:UP000002311) returns 6067 entries
       (one additional entry, P9WEJ6/YGR016C-A, that is not a canonical gene representative).
  rest https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=(proteome:UPxxxx)
       used for non-reference proteomes (not in Reference_Proteomes) or when the FTP file is missing; returns
       one canonical sequence per UniProtKB *entry* (a gene may have >1 entry), release in header X-UniProt-Release.
  --include-isoforms adds the *_additional.fasta file (ftp) / includeIsoform=true (rest).

Outputs in <out>/<UP>_<taxid>/ :
  proteome.fasta   headers rewritten to '>ACCESSION' (isoforms: '>ACCESSION-N')
  header_map.tsv   accession, original header, entry name, protein name, organism, taxid, gene name, PE, SV, length,
                   sequence sha256, gene symbol/group from gene2acc (ftp)
  manifest.json    source URLs, UniProt release + date, retrieval date, sha256 of raw download and outputs,
                   protein count, proteome type (reference?), assembly, license  -> use for DATA.md
  raw/             files exactly as downloaded
"""
from __future__ import annotations

import argparse
import csv
import gzip
import json
import re
import sys
from pathlib import Path

from . import common
from .common import HttpError, get, log, now_iso, sha256_file, seq_sha256

REST = "https://rest.uniprot.org"
FTP = "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes"
HDR = re.compile(r"^(?P<db>sp|tr)\|(?P<acc>[^|]+)\|(?P<entry>\S+)\s+(?P<name>.*?)\s+OS=(?P<os>.*?)\s+OX=(?P<ox>\d+)"
                 r"(?:\s+GN=(?P<gn>\S+))?(?:\s+PE=(?P<pe>\d))?(?:\s+SV=(?P<sv>\d+))?\s*$")
REFERENCE_TYPES = ("Reference proteome", "Reference and representative proteome", "Representative proteome")


def proteome_info(up: str) -> dict:
    r = get(f"{REST}/proteomes/{up}", ok_status=(200, 404))
    if r.status_code == 404:
        raise LookupError(f"proteome {up} not found in UniProt")
    d = r.json()
    t = d.get("taxonomy", {})
    return {"proteome_id": d["id"], "taxid": t.get("taxonId"), "organism": t.get("scientificName"),
            "common_name": t.get("commonName"), "proteome_type": d.get("proteomeType"),
            "is_reference_or_representative": d.get("proteomeType") in REFERENCE_TYPES,
            "superkingdom": d.get("superkingdom"), "protein_count_reported": d.get("proteinCount"),
            "gene_count_reported": d.get("geneCount"), "proteome_modified": d.get("modified"),
            "genome_assembly": (d.get("genomeAssembly") or {}).get("assemblyId"),
            "busco": d.get("proteomeCompletenessReport", {}).get("buscoReport"),
            "release_header": r.headers.get("X-UniProt-Release")}


def _download(url: str, dest: Path, optional: bool = False) -> dict | None:
    try:
        r = get(url, stream=True, timeout=600, ok_status=(200, 404) if optional else (200,))
    except HttpError as e:
        raise
    if r.status_code == 404:
        r.close()
        return None
    dest.parent.mkdir(parents=True, exist_ok=True)
    with open(dest, "wb") as f:
        for chunk in r.iter_content(1 << 20):
            f.write(chunk)
    meta = {"url": url, "file": dest.name, "bytes": dest.stat().st_size, "sha256": sha256_file(dest),
            "x_uniprot_release": r.headers.get("X-UniProt-Release"),
            "x_uniprot_release_date": r.headers.get("X-UniProt-Release-Date"),
            "x_total_results": r.headers.get("X-Total-Results")}
    r.close()
    return meta


def ftp_release() -> dict:
    """UniProt release of the FTP 'current_release' (RELEASE.metalink version + Reference Proteomes README date)."""
    out = {}
    try:
        txt = get("https://ftp.uniprot.org/pub/databases/uniprot/current_release/RELEASE.metalink").text
        m = re.search(r"<version>([^<]+)</version>", txt)
        out["release"] = m.group(1) if m else None
        rd = get(f"{FTP}/README").text
        m = re.search(r"Release\s+(\d{4}_\d{2}),\s*([0-9A-Za-z-]+)", rd)
        if m:
            out["release_readme"], out["release_date"] = m.group(1), m.group(2)
    except HttpError as e:
        log.warning("could not read FTP release information: %s", e)
    return out


def parse_header(h: str) -> dict:
    m = HDR.match(h)
    if not m:
        acc = h.split()[0].split("|")[1] if "|" in h.split()[0] else h.split()[0]
        return dict(db="", acc=acc, entry="", name=h, os="", ox="", gn="", pe="", sv="")
    d = m.groupdict()
    return {k: (v or "") for k, v in d.items()}


def rewrite(fasta: Path, out_fa: Path, out_map: Path, gene2acc: dict | None, label: str) -> tuple[int, int]:
    seen, n, n_res = set(), 0, 0
    recs = []
    for h, s in common.read_fasta(fasta):
        p = parse_header(h)
        acc = p["acc"]
        if acc in seen:
            raise ValueError(f"duplicate accession {acc} in {fasta}")
        seen.add(acc)
        recs.append((acc, h, s, p))
    with open(out_fa, "w") as fa, open(out_map, "w", newline="") as fm:
        w = csv.writer(fm, delimiter="\t", lineterminator="\n")
        w.writerow(["accession", "original_header", "db", "entry_name", "protein_name", "organism", "taxid",
                    "gene_name", "pe", "sv", "length", "seq_sha256", "gene_symbol_gene2acc", "canonical_of_gene_group",
                    "set"])
        for acc, h, s, p in recs:
            fa.write(f">{acc}\n")
            for i in range(0, len(s), 60):
                fa.write(s[i:i + 60] + "\n")
            g = (gene2acc or {}).get(acc, ("", ""))
            w.writerow([acc, h, p["db"], p["entry"], p["name"], p["os"], p["ox"], p["gn"], p["pe"], p["sv"], len(s),
                        seq_sha256(s), g[0], g[1], label])
            n += 1
            n_res += len(s)
    return n, n_res


def load_gene2acc(path: Path) -> dict:
    d = {}
    with gzip.open(path, "rt") as f:
        for ln in f:
            p = ln.rstrip("\n").split("\t")
            if len(p) >= 3:
                d.setdefault(p[1], (p[0], p[2]))
    return d


def fetch(up: str, out_root: str | Path = "data/proteomes", source: str = "auto", include_isoforms: bool = False) -> dict:
    up = up.upper()
    info = proteome_info(up)
    taxid = info["taxid"]
    outdir = Path(out_root) / f"{up}_{taxid}"
    raw = outdir / "raw"
    raw.mkdir(parents=True, exist_ok=True)
    downloads, notes = [], []
    chosen = source
    if source == "auto":
        chosen = "ftp" if info["is_reference_or_representative"] else "rest"
        notes.append(f"source=auto -> {chosen} (proteome type: {info['proteome_type']})")
    gene2acc = None
    release = {}
    parts: list[tuple[Path, str]] = []
    if chosen == "ftp":
        kd = (info["superkingdom"] or "").capitalize()
        base = f"{FTP}/{kd}/{up}"
        m = _download(f"{base}/{up}_{taxid}.fasta.gz", raw / f"{up}_{taxid}.fasta.gz", optional=True)
        if m is None:
            notes.append(f"FTP file not found under {base}; falling back to REST")
            log.warning(notes[-1])
            chosen = "rest"
        else:
            downloads.append(m)
            parts.append((raw / m["file"], "canonical"))
            g = _download(f"{base}/{up}_{taxid}.gene2acc.gz", raw / f"{up}_{taxid}.gene2acc.gz", optional=True)
            if g:
                downloads.append(g)
                gene2acc = load_gene2acc(raw / g["file"])
            else:
                notes.append("gene2acc file not available")
            if include_isoforms:
                a = _download(f"{base}/{up}_{taxid}_additional.fasta.gz", raw / f"{up}_{taxid}_additional.fasta.gz",
                              optional=True)
                if a:
                    downloads.append(a)
                    parts.append((raw / a["file"], "additional_isoforms"))
                else:
                    notes.append("no *_additional.fasta.gz for this proteome (no isoforms listed)")
            release = ftp_release()
    if chosen == "rest":
        url = f"{REST}/uniprotkb/stream"
        q = f"(proteome:{up})"
        dest = raw / f"{up}_rest.fasta"
        full = url + f"?format=fasta&query={q}" + ("&includeIsoform=true" if include_isoforms else "")
        m = _download(full, dest)
        downloads.append(m)
        parts.append((dest, "canonical" if not include_isoforms else "canonical+isoforms"))
        release = {"release": m["x_uniprot_release"], "release_date": m["x_uniprot_release_date"]}
        if info["protein_count_reported"] and m["x_total_results"] and int(m["x_total_results"]) != info["protein_count_reported"]:
            notes.append(f"REST returned {m['x_total_results']} entries but proteome reports {info['protein_count_reported']}")
    # rewrite
    fa_out, map_out = outdir / "proteome.fasta", outdir / "header_map.tsv"
    if len(parts) == 1:
        n, n_res = rewrite(parts[0][0], fa_out, map_out, gene2acc, parts[0][1])
    else:  # canonical + additional: concatenate headers into one map
        tmp = outdir / ".combined.fasta"
        with open(tmp, "w") as f:
            for p, _ in parts:
                for h, s in common.read_fasta(p):
                    f.write(f">{h}\n{s}\n")
        n, n_res = rewrite(tmp, fa_out, map_out, gene2acc, "canonical+additional")
        tmp.unlink()
    n_can = sum(1 for _ in common.read_fasta(parts[0][0]))
    if info["gene_count_reported"] and chosen == "ftp" and n_can != info["gene_count_reported"]:
        notes.append(f"canonical count {n_can} != reported geneCount {info['gene_count_reported']}")
    if (rel_h := info.get("release_header")) and release.get("release") and rel_h != release["release"]:
        notes.append(f"release mismatch: REST header {rel_h} vs FTP {release['release']} (release switchover?)")
    for nt in notes:
        log.info("note: %s", nt)
    gn = set()
    with open(map_out) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if r["gene_name"]:
                gn.add(r["gene_name"])
    manifest = {
        "dataset": "UniProt proteome", **{k: info[k] for k in (
            "proteome_id", "taxid", "organism", "common_name", "proteome_type", "is_reference_or_representative",
            "genome_assembly", "protein_count_reported", "gene_count_reported", "proteome_modified")},
        "source": chosen, "include_isoforms": include_isoforms,
        "uniprot_release": release.get("release"), "uniprot_release_date": release.get("release_date"),
        "retrieval_date": now_iso(),
        "downloads": downloads,
        "outputs": {"proteome_fasta": {"file": "proteome.fasta", "sha256": sha256_file(fa_out)},
                    "header_map": {"file": "header_map.tsv", "sha256": sha256_file(map_out)}},
        "n_proteins": n, "n_canonical": n_can, "n_residues": n_res, "n_distinct_gene_names": len(gn),
        "header_format": ">ACCESSION (isoforms >ACCESSION-N); original headers in header_map.tsv",
        "license": "CC BY 4.0 (UniProt)", "citation": "The UniProt Consortium, Nucleic Acids Res. (UniProt release "
                                                       f"{release.get('release')})",
        "notes": notes,
    }
    (outdir / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    return manifest


def markdown_row(m: dict) -> str:
    return (f"| {m['organism']} ({m['proteome_id']}, taxid {m['taxid']}) | UniProt {m['proteome_type']} via {m['source']} "
            f"| release {m['uniprot_release']} ({m['uniprot_release_date']}) | {m['n_proteins']} proteins | "
            f"sha256 {m['outputs']['proteome_fasta']['sha256'][:16]}... | {m['license']} |")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("proteome", nargs="+", help="UniProt proteome ids, e.g. UP000002311")
    ap.add_argument("--out", default="data/proteomes")
    ap.add_argument("--source", choices=["auto", "ftp", "rest"], default="auto")
    ap.add_argument("--include-isoforms", action="store_true")
    ap.add_argument("--markdown", action="store_true", help="print a DATA.md table row")
    ap.add_argument("-v", "--verbose", action="store_true")
    a = ap.parse_args(argv)
    common.setup_logging(a.verbose)
    rc = 0
    for up in a.proteome:
        try:
            m = fetch(up, a.out, a.source, a.include_isoforms)
        except (LookupError, HttpError, ValueError) as e:
            log.error("%s: %s", up, e)
            rc = 1
            continue
        print(json.dumps({k: m[k] for k in ("proteome_id", "organism", "proteome_type", "source", "uniprot_release",
                                             "n_proteins", "n_distinct_gene_names")} |
                         {"sha256": m["outputs"]["proteome_fasta"]["sha256"]}))
        if a.markdown:
            print(markdown_row(m))
    return rc


if __name__ == "__main__":
    sys.exit(main())
