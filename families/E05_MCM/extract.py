#!/usr/bin/env python
"""E05_MCM - extract (paralogue label, species, accession) from Liu, Richards & Aves 2009 Additional files 1-4.

Additional files 1-4 are vector PDF trees (one page each), captions (fulltext.txt [SUPPLEMENT S1-S4]):
  S1 = Mcm4, Mcm7; S2 = Mcm5, Mcm6; S3 = Mcm3, Mcm2; S4 = Mcm8, Mcm9, archaeal MCMs.
Each leaf is printed as "<Species> [<accession>]". Leaves are read with their page coordinates (pypdf), sorted
top-down, and cut into paralogue blocks:
  1. a vertical gap > GAP points between consecutive leaves (blocks are separated by collapsed-clade labels), and
  2. archaeal leaves (the 8 archaea named in Methods) form their own block;
  3. in S2 (no gap between Mcm5 and Mcm6) a block is cut where the first species of the block re-occurs
     (Mus musculus opens both blocks).
Blocks are labelled in caption order (top to bottom); the paralogue labels printed on the figure (e.g. 'Mcm4' at
y=636, 'Mcm7' at y=300) must fall inside the y-range of their block -- checked in check_blocks().
"""
import re
from pathlib import Path

import pypdf

ROOT = Path(__file__).resolve().parents[2]
PAPER = "literature/papers/PMID19292915"
SUPP = ROOT / PAPER / "raw/supp"
CAPTION_ORDER = {1: ["Mcm4", "Mcm7"], 2: ["Mcm5", "Mcm6"], 3: ["Mcm3", "Mcm2"], 4: ["Mcm8", "Mcm9", "Archaea"]}
GAP = 12.0
ARCHAEA = {"Thermofilum pendens", "Sulfolobus acidocaldarius", "Pyrococcus horikoshii", "Halobacterium sp.",
           "Archaeoglobus fulgidus", "Thermoplasma volcanium", "Nanoarchaeum equitans", "Methanococcus jannaschii"}
# paper species name (as printed, typos included) -> panel E code (only species present in panel E)
SPECIES2CODE = {
    "Homo sapiens": "HUMAN", "Drosophila melanogaster": "DROME", "Nematostella vectensis": "NEMVE",
    "Trichoplax adhaerens": "TRIAD", "Monosiga brevicollis": "MONBE", "Saccharomyces cerevisiae": "YEAST",
    "Neurospora crassa": "NEUCR", "Batrachochytrium dendrobatidis": "BATDJ", "Encephalitozoon cuniculi": "ENCCU",
    "Dictyostelium discoideum": "DICDI", "Entamoeba histolytica": "ENTH1", "Entamoeba Histolytica": "ENTH1",
    "Arabidopsis thaliana": "ARATH", "Chlamydomonas reinhardtii": "CHLRE", "Chlaydomonas reinhardtti": "CHLRE",
    "Cyanidioschyzon merolae": "CYAME", "Thalassiosira pseudonana": "THAPS", "Phytophthora ramorum": "PHYRM",
    "Trypanosoma brucei": "TRYB2", "Trypansoma brucei": "TRYB2", "Naegleria gruberi": "NAEGR",
    "Giardia lamblia": "GIAIC", "Halobacterium sp.": "HALSA"}
LEAF = re.compile(r"^(?P<sp>[A-Z][a-z]+ (?:[a-z]+|sp\.)|[A-Z][a-z]+ [A-Z][a-z]+)\s*\[(?P<acc>[^\]]+)\]")


def read_items(i):
    items = []

    def v(text, cm, tm, fd, fs):
        if text.strip():
            items.append((cm[4] + tm[4] * cm[0], cm[5] + tm[5] * cm[3], text.strip()))
    pypdf.PdfReader(SUPP / f"1471-2148-9-60-S{i}.pdf").pages[0].extract_text(visitor_text=v)
    return items


def leaves(i):
    """Group text items on one baseline (|dy|<1) and parse '<species> [acc]'."""
    items = sorted(read_items(i), key=lambda t: (-round(t[1], 0), t[0]))
    lines = {}
    for x, y, t in items:
        key = next((k for k in lines if abs(k - y) < 1.0), y)
        lines.setdefault(key, []).append((x, t))
    out = []
    for y, parts in lines.items():
        txt = " ".join(t for _, t in sorted(parts))
        for m in re.finditer(r"([A-Z][a-z]+ (?:[a-z]+|[A-Z][a-z]+|sp\.))\s*\[([^\]]+)\]", txt):
            out.append(dict(y=round(y, 1), species=m.group(1), acc=m.group(2).strip()))
    return sorted(out, key=lambda d: -d["y"])


def labels(i):
    return [(round(y, 1), t) for x, y, t in read_items(i) if re.match(r"^\+? ?(Mcm|Archaea)", t)]


def blocks(i):
    lv = leaves(i)
    euk = [d for d in lv if d["species"] not in ARCHAEA]
    arc = [d for d in lv if d["species"] in ARCHAEA]
    bl = [[euk[0]]]
    for d in euk[1:]:
        if bl[-1][-1]["y"] - d["y"] > GAP:
            bl.append([])
        elif i == 2 and d["species"] == bl[-1][0]["species"]:
            bl.append([])
        bl[-1].append(d)
    if arc:
        bl.append(arc)
    names = CAPTION_ORDER[i]
    assert len(bl) == len(names), (i, [len(b) for b in bl])
    for b, n in zip(bl, names):
        for d in b:
            d["paralog"] = n
            d["figure"] = f"Additional file {i}"
            d["file"] = f"{PAPER}/raw/supp/1471-2148-9-60-S{i}.pdf"
    return bl


def check_blocks(i, bl):
    """The single-paralogue label printed on the figure must lie within its block's y-range."""
    probs = []
    for y, t in labels(i):
        n = t.replace(" ", "")
        for b in bl:
            if b[0]["paralog"] == n or (n == "ArchaeaMCMs" and b[0]["paralog"] == "Archaea" and i == 4):
                lo, hi = min(d["y"] for d in b), max(d["y"] for d in b)
                if not (lo - 5 <= y <= hi + 5):
                    probs.append((i, t, y, lo, hi))
    return probs


def all_leaves():
    rows = []
    for i in (1, 2, 3, 4):
        bl = blocks(i)
        assert not check_blocks(i, bl), check_blocks(i, bl)
        for b in bl:
            rows.extend(b)
    for d in rows:
        d["code"] = SPECIES2CODE.get(d["species"], "")
    return rows


if __name__ == "__main__":
    for d in all_leaves():
        print(d["figure"], d["paralog"], d["species"], d["acc"], d["code"], d["y"], sep="\t")


# ---- identifier lookup for genome-database ids (offline; panel header_map.tsv) ----
def id_pattern(code, acc):
    """Paper genome-database id -> (UniProt ORF/gene-name token, source db) for ids without NCBI records."""
    m = re.match(r"^(Naegr1|Triad1|Nemve1|Monbr1|Thaps3|Phyra1_1):(\d+)", acc)
    if m:
        db, n = m.groups()
        tok = {"Naegr1": f"NAEGRDRAFT_{n}", "Triad1": f"TRIADDRAFT_{n}", "Nemve1": f"NEMVEDRAFT_v1g{n}",
               "Monbr1": f"MONBRDRAFT_{n}", "Thaps3": f"THAPSDRAFT_{n}", "Phyra1_1": f"Phyra{n}"}[db]
        return tok, f"JGI {db} protein id"
    m = re.match(r"^(CM[A-Z]\d{3}[CR])", acc)
    if m:
        return f"CYME_{m.group(1)}", "CyanidioschyzonDB locus"
    m = re.match(r"^GL50803_(\d+)", acc)
    if m:
        return f"GL50803_#{int(m.group(1))}", "GiardiaDB locus"
    return None, None


def header_index(code):
    import csv
    sp = {r["code"]: r for r in csv.DictReader(open(ROOT / "panels/E/species.tsv"), delimiter="\t")}
    hm = Path(ROOT / sp[code]["fasta_path"]).parent / "header_map.tsv"
    idx = {}
    for r in csv.DictReader(open(hm), delimiter="\t"):
        toks = set((r.get("gene_name") or "").split()) | set((r.get("gene_symbol_gene2acc") or "").split())
        toks |= set(re.findall(r"GN=(\S+)", r["original_header"]))
        for t in toks:
            m = re.match(r"^GL50803_(\d+)$", t)  # GiardiaDB loci are zero-padded to 6 or 7 digits in UniProt
            idx.setdefault(f"GL50803_#{int(m.group(1))}" if m else t, set()).add(r["accession"])
    return idx, hm
