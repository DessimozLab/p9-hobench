#!/usr/bin/env python3
"""Build the benchmark species panels: panels/<P>/{species.tsv, proteomes/<CODE>.fa, species_tree.nwk,
species_tree_alt.nwk, levels.tsv} for P in Y K V E P (panel D dropped, reviews/round1/SYNTHESIS.md).

Sources (no silent assumptions):
  * species lists: literature/raw/C/panels_verified.tsv rows with tier == core (survey C, section 3), plus the
    explicit EXTRA species below (SYNTHESIS A5/A6, PI decision 2026-10-26);
  * proteomes: data/proteomes/<UPID>_<taxid>/ written by `python -m curation_tools.fetch_proteome` (manifest.json
    gives release, counts, sha256; the sha256 is re-computed here and must match);
  * species trees: hand-written Newick strings below, each internal node named by a level id
    (^[A-Za-z0-9_-]{3,}$), topology and contested nodes cited in LEVELS (DOIs checked in Europe PMC,
    cache panels/citations.tsv; refresh with --verify-dois);
  * NCBI taxid/name per level: computed with ete3 from data/raw/ncbi/taxa.sqlite (built from taxdump.tar.gz)
    as the NCBI LCA of the level's species; `ncbi_exact` = y iff that NCBI taxon contains exactly the level's
    panel species (i.e. the level corresponds to an NCBI taxon restricted to the panel).

Usage: python scripts/build_panels.py [--panels Y K] [--verify-dois] [--no-sha]
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import re
import sys
import time
import urllib.parse
import urllib.request
from collections import Counter
from pathlib import Path

from ete3 import Tree

ROOT = Path(__file__).resolve().parents[1]
PROTEOMES = ROOT / "data" / "proteomes"
PANELS_DIR = ROOT / "panels"
VERIFIED = ROOT / "literature" / "raw" / "C" / "panels_verified.tsv"
TAXDB = ROOT / "data" / "raw" / "ncbi" / "taxa.sqlite"
CITATIONS = PANELS_DIR / "citations.tsv"
PINNED_RELEASE = "2026_03"
LEVEL_RE = re.compile(r"^[A-Za-z0-9_-]{3,}$")

# ---------------------------------------------------------------------------------------------------------------
# Citations (key -> DOI). Every key used in LEVELS must be here; --verify-dois checks each DOI in Europe PMC.
CIT = {
    "Shen2018": "10.1016/j.cell.2018.10.023",
    "Shen2016": "10.1534/g3.116.034744",
    "MarcetHouben2015": "10.1371/journal.pbio.1002220",
    "Wolfe2015": "10.1371/journal.pbio.1002221",
    "WolfeShields1997": "10.1038/42711",
    "Coleman2021": "10.1126/science.abe0511",
    "Hug2016": "10.1038/nmicrobiol.2016.48",
    "Parks2018": "10.1038/nbt.4229",
    "Parks2022GTDB": "10.1093/nar/gkab776",
    "Williams2020": "10.1038/s41559-019-1040-x",
    "Eme2023": "10.1038/s41586-023-06186-2",
}

# ---------------------------------------------------------------------------------------------------------------
# Extra species beyond the survey-C core lists: code -> (upid, why)
EXTRA = {
    "V": {},
    "P": {},
}

# Species-specific notes added to the `why` column
NOTES = {
    "ZYGRC": "FTP canonical file not found because REST omits superkingdom; fetched via REST stream (4,987 = geneCount)",
}

# ---------------------------------------------------------------------------------------------------------------
# Species trees. Leaves = panel codes; internal labels = level ids. Branch lengths deliberately omitted
# (topology-only reference; tools that need lengths estimate them).
#
# Panel Y: pruned from Shen et al. 2018 Cell (doi 10.1016/j.cell.2018.10.023), figshare 10.6084/m9.figshare.5854692
# file data_in_Figure2/343_2408OGs_unpartitioned_ML_species-names_updated.newick (CC BY 4.0), rooted on
# S. pombe, pruned with ete3 to the 16 panel species (Lachancea_thermotolerans, not Phaffomyces_thermotolerans).
# The 332-taxon time-calibrated tree in the same archive gives the identical topology for the 15 species it has.
# Bootstrap support after pruning: 100 at every node except Saccharomyces+Nakaseomyces vs
# Naumovozyma+Kazachstania (52) -> PostWGD_crown is flagged contested.
TREES = {
    "Y": "(SCHPO,(YARLI,(CANAL,((((ERECY,ASHGO)Eremothecium,KLULA)Kluyveromyces_Eremothecium,LACTH)KLE_clade,"
         "((TORDC,ZYGRC)ZT_clade,((((NAUDC,NAUCC)Naumovozyma,KAZAF)Naumovozyma_Kazachstania,"
         "(YEAST,CANGA)Saccharomyces_Nakaseomyces)PostWGD_crown,(VANPO,TETPH)Vanderwaltozyma_Tetrapisispora)PostWGD)"
         "ZT_PostWGD)Saccharomycetaceae)CUG_Saccharomycetaceae)Saccharomycotina)Ascomycota;",
    # Panel K: two-domain rooting between Bacteria and Archaea (Williams 2020). Bacteria rooted between
    # Terrabacteria and Gracilicutes (Coleman 2021). Phylum-level relationships inside the two bacterial
    # superphyla and inside archaea follow GTDB-based phylogenomics (Parks 2018/2022) and Hug 2016 and are
    # LOW CONFIDENCE (HGT-rich, short internal branches): every non-trivial internal bacterial node is flagged
    # contested in levels.tsv. Placement of Thermotogota (with Deinococcota as "DST" in Terrabacteria) and of
    # Aquificota (with Campylobacterota in Gracilicutes) is the most uncertain choice.
    "K": "((((DEIRA,THEMA)DST_clade,(SYNY3,(BACSU,MYCTU)Bacillota_Actinomycetota)Terrabacteria_core)Terrabacteria,"
         "((LEPIN,(BACTN,CHLTR)FCB_PVC)Spirochaetota_FCB_PVC,(((ECOLI,PSEAE)Gammaproteobacteria,BRADU)Pseudomonadota,"
         "(GEOSL,(HELPY,AQUAE)Campylobacterota_Aquificota)Desulfobacterota_Campylobacterota)Proteobacteria_sl)"
         "Gracilicutes)Bacteria,"
         "((THEKO,(METJA,(METAC,HALSA)Methanomicrobia_Halobacteria)Euryarchaeota_core)Euryarchaeota,"
         "(NITMS,(KORCO,SACS2)Crenarchaeota_Korarchaeota)TACK)Archaea)LUCA;",
}

# One plausible alternative topology per panel at a contested node (species-tree-error robustness experiment).
ALT_TREES = {
    # Alternative at ZT_PostWGD: post-WGD lineage sister to KLE instead of ZT, i.e. the signal of the KLE-like
    # parent of the hybrid WGD ancestor (Marcet-Houben & Gabaldon 2015; Wolfe 2015).
    "Y": "(SCHPO,(YARLI,(CANAL,((TORDC,ZYGRC)ZT_clade,((((ERECY,ASHGO)Eremothecium,KLULA)Kluyveromyces_Eremothecium,"
         "LACTH)KLE_clade,((((NAUDC,NAUCC)Naumovozyma,KAZAF)Naumovozyma_Kazachstania,"
         "(YEAST,CANGA)Saccharomyces_Nakaseomyces)PostWGD_crown,(VANPO,TETPH)Vanderwaltozyma_Tetrapisispora)PostWGD)"
         "KLE_PostWGD)Saccharomycetaceae)CUG_Saccharomycetaceae)Saccharomycotina)Ascomycota;",
    # Alternative: Aquificota sister to Thermotogota (classic hyperthermophile grouping of rRNA trees) inside
    # Terrabacteria instead of with Campylobacterota in Gracilicutes.
    "K": "((((DEIRA,(THEMA,AQUAE)Thermotogota_Aquificota)DST_clade,(SYNY3,(BACSU,MYCTU)Bacillota_Actinomycetota)"
         "Terrabacteria_core)Terrabacteria,"
         "((LEPIN,(BACTN,CHLTR)FCB_PVC)Spirochaetota_FCB_PVC,(((ECOLI,PSEAE)Gammaproteobacteria,BRADU)Pseudomonadota,"
         "(GEOSL,HELPY)Desulfobacterota_Campylobacterota)Proteobacteria_sl)Gracilicutes)Bacteria,"
         "((THEKO,(METJA,(METAC,HALSA)Methanomicrobia_Halobacteria)Euryarchaeota_core)Euryarchaeota,"
         "(NITMS,(KORCO,SACS2)Crenarchaeota_Korarchaeota)TACK)Archaea)LUCA;",
}

# Level annotations: level_id -> (literature_name, contested y/n, alternative topology, citation keys, note)
_ = None
LEVELS = {
    "Y": {
        "Ascomycota": ("Ascomycota (Taphrinomycotina + Saccharomycotina)", "n", "", ["Shen2018"], "root; SCHPO outgroup"),
        "Saccharomycotina": ("Saccharomycotina", "n", "", ["Shen2018", "Shen2016"], ""),
        "CUG_Saccharomycetaceae": ("CUG-Ser1 clade + Saccharomycetaceae (Shen 2018 major clades)", "n", "",
                                   ["Shen2018"], "Yarrowia (Dipodascaceae/Trigonopsidaceae) branches outside"),
        "Saccharomycetaceae": ("Saccharomycetaceae", "n", "", ["Shen2018"], ""),
        "KLE_clade": ("KLE clade (Kluyveromyces, Lachancea, Eremothecium)", "y",
                      "Lachancea outside KLE / KLE paraphyletic in some analyses; KLE-like parent of the WGD hybrid",
                      ["Shen2018", "MarcetHouben2015"], "BS 100 in Shen 2018; treat as soft for WGD-gene scoring"),
        "Kluyveromyces_Eremothecium": ("Kluyveromyces + Eremothecium", "y",
                                       "(Lachancea,Kluyveromyces) sister to Eremothecium",
                                       ["Shen2018"], "Lachancea is the unstable member (survey C 4.2)"),
        "Eremothecium": ("Eremothecium", "n", "", ["Shen2018"], ""),
        "ZT_PostWGD": ("ZT clade + post-WGD clade", "y",
                       "post-WGD sister to KLE (species_tree_alt.nwk); ~30% of post-WGD genes cluster with KLE "
                       "because the WGD followed an interspecies hybridization",
                       ["Shen2018", "MarcetHouben2015", "Wolfe2015"], "species tree = ZT parent signal"),
        "ZT_clade": ("ZT clade (Zygosaccharomyces, Torulaspora)", "n", "", ["Shen2018"], ""),
        "PostWGD": ("post-WGD clade", "n", "", ["Shen2018", "WolfeShields1997"], "WGD on the branch into this level"),
        "PostWGD_crown": ("Saccharomyces/Nakaseomyces + Naumovozyma/Kazachstania", "y",
                          "(Vanderwaltozyma,Tetrapisispora) sister to Naumovozyma/Kazachstania",
                          ["Shen2018"], "bootstrap 52 in the pruned Shen 2018 ML tree"),
        "Naumovozyma_Kazachstania": ("Naumovozyma + Kazachstania", "n", "", ["Shen2018"], ""),
        "Naumovozyma": ("Naumovozyma", "n", "", ["Shen2018"], ""),
        "Saccharomyces_Nakaseomyces": ("Saccharomyces + Nakaseomyces", "n", "", ["Shen2018"], ""),
        "Vanderwaltozyma_Tetrapisispora": ("Vanderwaltozyma + Tetrapisispora", "n", "", ["Shen2018"], ""),
        "KLE_PostWGD": ("KLE + post-WGD (alternative tree only)", "y", "alt only", ["MarcetHouben2015"], ""),
    },
    "K": {
        "LUCA": ("Last universal common ancestor (Bacteria/Archaea split)", "y",
                 "root position inside Bacteria (not assessed)", ["Williams2020", "Coleman2021"],
                 "two-domain rooting between Bacteria and Archaea"),
        "Bacteria": ("Bacteria", "y", "root of Bacteria elsewhere (e.g. CPR- or Thermotogota-first)",
                     ["Coleman2021"], "rooted between Terrabacteria and Gracilicutes"),
        "Terrabacteria": ("Terrabacteria", "y", "membership of Thermotogota disputed", ["Coleman2021", "Hug2016"], ""),
        "DST_clade": ("Deinococcota + Thermotogota (DST: Deinococcota, Synergistota, Thermotogota)", "y",
                      "Thermotogota near the bacterial root or with Aquificota (species_tree_alt.nwk)",
                      ["Coleman2021", "Hug2016"], "low confidence"),
        "Terrabacteria_core": ("Cyanobacteriota + Bacillota + Actinomycetota", "y",
                               "Cyanobacteriota sister to Actinomycetota or to Chloroflexota-group",
                               ["Coleman2021", "Parks2022GTDB"], "low confidence"),
        "Bacillota_Actinomycetota": ("Bacillota (Firmicutes) + Actinomycetota", "y",
                                     "Bacillota sister to the rest of Terrabacteria",
                                     ["Coleman2021", "Parks2022GTDB"], "low confidence"),
        "Gracilicutes": ("Gracilicutes", "n", "", ["Coleman2021"], ""),
        "Spirochaetota_FCB_PVC": ("Spirochaetota + FCB + PVC", "y", "Spirochaetota sister to Proteobacteria s.l.",
                                  ["Coleman2021", "Hug2016"], "low confidence"),
        "FCB_PVC": ("FCB (Bacteroidota) + PVC (Chlamydiota) superphyla", "y", "FCB and PVC not sisters",
                    ["Hug2016", "Coleman2021"], "low confidence"),
        "Proteobacteria_sl": ("Proteobacteria sensu lato (incl. Desulfobacterota, Campylobacterota)", "y",
                              "Proteobacteria s.l. paraphyletic (GTDB splits it)", ["Parks2018", "Parks2022GTDB"], ""),
        "Pseudomonadota": ("Pseudomonadota (Alpha + Gammaproteobacteria)", "n", "", ["Parks2018", "Parks2022GTDB"], ""),
        "Gammaproteobacteria": ("Gammaproteobacteria", "n", "", ["Parks2022GTDB"], ""),
        "Desulfobacterota_Campylobacterota": ("Desulfobacterota + Campylobacterota (+Aquificota)", "y",
                                              "Campylobacterota sister to Pseudomonadota",
                                              ["Parks2018", "Coleman2021"], "low confidence"),
        "Campylobacterota_Aquificota": ("Campylobacterota + Aquificota", "y",
                                        "Aquificota with Thermotogota (species_tree_alt.nwk) or near the root",
                                        ["Coleman2021", "Hug2016"], "low confidence; extensive HGT in Aquificota"),
        "Thermotogota_Aquificota": ("Thermotogota + Aquificota (alternative tree only)", "y", "alt only",
                                    ["Hug2016"], ""),
        "Archaea": ("Archaea", "n", "", ["Williams2020", "Eme2023"], ""),
        "Euryarchaeota": ("Euryarchaeota (Thermococci, Methanococci, Methanomicrobia, Halobacteria)", "y",
                          "Euryarchaeota paraphyletic w.r.t. archaeal root", ["Williams2020", "Parks2022GTDB"], ""),
        "Euryarchaeota_core": ("Methanococci + Methanomicrobia/Halobacteria", "y",
                               "Thermococci sister to Methanococci", ["Parks2022GTDB"], "low confidence"),
        "Methanomicrobia_Halobacteria": ("Methanomicrobia + Halobacteria (Stenosarchaea)", "n", "",
                                         ["Parks2022GTDB"], ""),
        "TACK": ("TACK superphylum (Thaumarchaeota, Crenarchaeota, Korarchaeota)", "n", "",
                 ["Williams2020", "Eme2023"], ""),
        "Crenarchaeota_Korarchaeota": ("Crenarchaeota + Korarchaeota", "y", "Korarchaeota sister to Thaum+Cren",
                                       ["Williams2020"], "low confidence"),
    },
}


# ---------------------------------------------------------------------------------------------------------------
def sha256_file(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def proteome_dir(upid: str) -> Path:
    c = sorted(PROTEOMES.glob(f"{upid}_*"))
    if len(c) != 1:
        raise SystemExit(f"ERROR: expected exactly one proteome dir for {upid}, found {c}")
    return c[0]


def uniprot_mnemonic(hmap: Path) -> str:
    cnt = Counter()
    with open(hmap) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if "_" in r["entry_name"]:
                cnt[r["entry_name"].rsplit("_", 1)[1]] += 1
    return cnt.most_common(1)[0][0] if cnt else ""


def panel_species(panel: str) -> list[dict]:
    rows = [r for r in csv.DictReader(open(VERIFIED), delimiter="\t") if r["panel"] == panel and r["tier"] == "core"]
    out = [dict(code=r["mnemonic"], upid=r["upid"], why=r["why"],
                source_set="QfO" if r["QfO2026_02"] == "yes" else "RP") for r in rows]
    for code, (upid, why) in EXTRA.get(panel, {}).items():
        out.append(dict(code=code, upid=upid, why=why, source_set="RP"))
    return out


def build_species(panel: str, check_sha: bool) -> list[dict]:
    pdir = PANELS_DIR / panel
    (pdir / "proteomes").mkdir(parents=True, exist_ok=True)
    rows = []
    for s in panel_species(panel):
        d = proteome_dir(s["upid"])
        m = json.loads((d / "manifest.json").read_text())
        if m["uniprot_release"] != PINNED_RELEASE:
            raise SystemExit(f"ERROR: {s['code']} {s['upid']} release {m['uniprot_release']} != {PINNED_RELEASE}")
        fa = d / "proteome.fasta"
        sha = m["outputs"]["proteome_fasta"]["sha256"]
        if check_sha and sha256_file(fa) != sha:
            raise SystemExit(f"ERROR: sha256 of {fa} differs from manifest")
        n = sum(1 for ln in open(fa) if ln.startswith(">"))
        if n != m["n_proteins"] or n == 0:
            raise SystemExit(f"ERROR: {fa}: {n} sequences vs manifest {m['n_proteins']}")
        link = pdir / "proteomes" / f"{s['code']}.fa"
        target = os.path.relpath(fa, link.parent)
        if link.is_symlink() or link.exists():
            link.unlink()
        link.symlink_to(target)
        why = s["why"] + (f"; {NOTES[s['code']]}" if s["code"] in NOTES else "")
        rows.append(dict(code=s["code"], species=m["organism"] or "", taxid=m["taxid"], upid=s["upid"],
                         mnemonic=uniprot_mnemonic(d / "header_map.tsv"), uniprot_release=m["uniprot_release"],
                         proteome_type=m["proteome_type"], fetch_source=m["source"],
                         fasta_path=os.path.relpath(fa, ROOT), n_proteins=n, sha256=sha,
                         source_set=s["source_set"], why=why))
    # remove stale links
    codes = {r["code"] for r in rows}
    for f in (pdir / "proteomes").glob("*.fa"):
        if f.stem not in codes:
            f.unlink()
    if len(codes) != len(rows):
        raise SystemExit(f"ERROR: duplicate codes in panel {panel}")
    with open(pdir / "species.tsv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    return rows


def validate_tree(nwk: str, codes: set[str], what: str) -> Tree:
    t = Tree(nwk, format=8)
    errs = []
    if len(t.children) != 2:
        errs.append(f"root has {len(t.children)} children (must be rooted binary)")
    for n in t.traverse():
        if not n.is_leaf() and len(n.children) != 2:
            errs.append(f"node {n.name} has {len(n.children)} children")
        if not n.is_leaf() and not LEVEL_RE.match(n.name or ""):
            errs.append(f"bad level id {n.name!r}")
    names = [n.name for n in t.traverse()]
    dup = [k for k, v in Counter(names).items() if v > 1]
    if dup:
        errs.append(f"duplicate node names {dup}")
    leaves = set(t.get_leaf_names())
    if leaves != codes:
        errs.append(f"leaves != species codes: missing {sorted(codes - leaves)}, extra {sorted(leaves - codes)}")
    if errs:
        raise SystemExit(f"ERROR in {what}:\n  " + "\n  ".join(errs))
    return t


def load_citations() -> dict:
    if not CITATIONS.exists():
        return {}
    return {r["key"]: r for r in csv.DictReader(open(CITATIONS), delimiter="\t")}


def verify_dois(keys) -> None:
    have = load_citations()
    for k in sorted(keys):
        doi = CIT[k]
        if k in have and have[k]["doi"] == doi and have[k]["status"] == "VERIFIED":
            continue
        url = ("https://www.ebi.ac.uk/europepmc/webservices/rest/search?" +
               urllib.parse.urlencode({"query": f'DOI:"{doi}"', "format": "json", "resultType": "lite"}))
        res = json.load(urllib.request.urlopen(url, timeout=60))["resultList"]["result"]
        r = res[0] if res else {}
        ok = bool(r) and (r.get("doi") or "").lower() == doi.lower()
        have[k] = dict(key=k, doi=doi, status="VERIFIED" if ok else "NOT_FOUND", title=r.get("title", ""),
                       journal=r.get("journalTitle", ""), year=r.get("pubYear", ""),
                       checked=time.strftime("%Y-%m-%d"), source="Europe PMC REST search DOI:")
        time.sleep(0.2)
    with open(CITATIONS, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["key", "doi", "status", "title", "journal", "year", "checked", "source"],
                           delimiter="\t", lineterminator="\n")
        w.writeheader()
        for k in sorted(have):
            w.writerow(have[k])


def ncbi_annot(t: Tree, taxid_of: dict, ncbi) -> dict:
    """level -> (ncbi_taxid, ncbi_name, exact y/n)."""
    lineages = {c: ncbi.get_lineage(int(x)) for c, x in taxid_of.items()}
    out = {}
    for n in t.traverse():
        if n.is_leaf():
            continue
        sp = n.get_leaf_names()
        common = set(lineages[sp[0]])
        for s in sp[1:]:
            common &= set(lineages[s])
        lca = [x for x in lineages[sp[0]] if x in common][-1]
        inside = {c for c, lin in lineages.items() if lca in lin}
        out[n.name] = (lca, ncbi.get_taxid_translator([lca])[lca], "y" if inside == set(sp) else "n")
    return out


def harmonise_alt(t: Tree, alt: Tree) -> None:
    """A level id shared by main and alt trees must denote the same species set; otherwise the alt node is
    renamed <id>_alt (so scoring never confuses two different clades under one name)."""
    main = {n.name: frozenset(n.get_leaf_names()) for n in t.traverse() if not n.is_leaf()}
    for n in alt.traverse():
        if not n.is_leaf() and n.name in main and main[n.name] != frozenset(n.get_leaf_names()):
            n.name = n.name + "_alt"


def write_levels(panel: str, t: Tree, alt: Tree, taxid_of: dict, ncbi, cits: dict) -> None:
    meta = LEVELS[panel]
    ann = ncbi_annot(t, taxid_of, ncbi)
    alt_ann = ncbi_annot(alt, taxid_of, ncbi)
    rows = []
    for tree, tag, a in ((t, "main", ann), (alt, "alt_only", alt_ann)):
        for n in tree.traverse("preorder"):
            if n.is_leaf() or (tag == "alt_only" and n.name in ann):
                continue
            base = n.name[:-4] if n.name.endswith("_alt") and n.name not in meta else n.name
            if base not in meta:
                raise SystemExit(f"ERROR: level {n.name} of panel {panel} has no LEVELS annotation")
            lit, contested, alt_txt, keys, note = meta[base]
            if base != n.name:
                contested, alt_txt, note = "y", "alt only", f"alternative-tree clade differing from main-tree {base}"
            for k in keys:
                if k not in CIT:
                    raise SystemExit(f"ERROR: citation key {k} not in CIT")
                if cits.get(k, {}).get("status") != "VERIFIED":
                    raise SystemExit(f"ERROR: citation {k} not verified (run --verify-dois)")
            rows.append(dict(level_id=n.name, tree=tag, ncbi_taxid=a[n.name][0], ncbi_name=a[n.name][1],
                             ncbi_exact=a[n.name][2], literature_name=lit,
                             children=",".join(c.name for c in n.children), n_species=len(n.get_leaf_names()),
                             contested=contested, alternative=alt_txt,
                             citation="; ".join(f"{k} doi:{CIT[k]}" for k in keys), note=note))
    unused = set(meta) - {r["level_id"] for r in rows} - {r["level_id"][:-4] for r in rows}
    if unused:
        raise SystemExit(f"ERROR: LEVELS[{panel}] has annotations for absent levels {sorted(unused)}")
    with open(PANELS_DIR / panel / "levels.tsv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--panels", nargs="+", default=list(TREES))
    ap.add_argument("--verify-dois", action="store_true")
    ap.add_argument("--no-sha", action="store_true", help="skip re-hashing FASTAs (still uses manifest sha256)")
    a = ap.parse_args(argv)
    from ete3 import NCBITaxa
    if not TAXDB.exists():
        raise SystemExit(f"ERROR: {TAXDB} missing; build with NCBITaxa(dbfile=..., taxdump_file='taxdump.tar.gz')")
    ncbi = NCBITaxa(dbfile=str(TAXDB))
    PANELS_DIR.mkdir(exist_ok=True)
    if a.verify_dois:
        verify_dois({k for p in a.panels for v in LEVELS[p].values() for k in v[3]})
    cits = load_citations()
    sys.path.insert(0, str(ROOT))
    from hobench.speciestree import SpeciesTree
    for p in a.panels:
        rows = build_species(p, not a.no_sha)
        codes = {r["code"] for r in rows}
        fa_codes = {f.stem for f in (PANELS_DIR / p / "proteomes").glob("*.fa")}
        assert fa_codes == codes, (p, fa_codes ^ codes)
        t = validate_tree(TREES[p], codes, f"{p} species_tree")
        alt = validate_tree(ALT_TREES[p], codes, f"{p} species_tree_alt")
        harmonise_alt(t, alt)
        if t.robinson_foulds(alt, unrooted_trees=False)[0] == 0:
            raise SystemExit(f"ERROR: {p} alt tree identical to main tree")
        for tr, fn in ((t, "species_tree.nwk"), (alt, "species_tree_alt.nwk")):
            (PANELS_DIR / p / fn).write_text(tr.write(format=8, format_root_node=True) + "\n")
            SpeciesTree.from_newick(PANELS_DIR / p / fn)  # project parser must accept it
        write_levels(p, t, alt, {r["code"]: r["taxid"] for r in rows}, ncbi, cits)
        print(f"panel {p}: {len(rows)} species, {sum(r['n_proteins'] for r in rows):,} proteins, "
              f"{len(t) - 1} internal levels; RF(main, alt) = {t.robinson_foulds(alt, unrooted_trees=False)[0]}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
