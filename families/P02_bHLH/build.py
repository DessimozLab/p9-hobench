#!/usr/bin/env python
"""P02_bHLH: plant basic-helix-loop-helix transcription factors (Pires & Dolan 2010, MBE 27:862).

Sources: File009 (species, protein name, locus), File011 (FASTA of all 566 proteins), Supplementary Figure 2
(msp288_2.pdf: table 'Subfamily / Proteins' assigning the 544 plant proteins to 26 subfamilies + orphans; text
extracted with pypdf).

Mapping onto panel P proteomes (same species only): sequence cross-reference (exact or >=98 % id, >=80 % cov,
RBH, >=2-point margin), else exact identifier mapping (AGI, LOC_Os, PHYPADRAFT_, CHLREDRAFT_), single panel hit.

Tree (subfamily dating quoted in assertions.tsv, T1):
  root U[ D_land(lo=Embryophyta)[20 subfamilies with moss members],
          D_vasc(lo=Tracheophyta)[4 subfamilies with lycophyte+angiosperm members only],
          D_angio(lo=Mesangiospermae)[2 angiosperm-only subfamilies], orphans ]
  (the parent subfamilies of the younger subfamilies are not identified -> pairs across these groups are U)
  every interval has lvl_hi=ABOVE_RANGE (the paper does not resolve inter-subfamily relationships);
  inside a subfamily: U (several paralogs per species).
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))  # plib.py: copy of families/_shared_P/plib.py (self-contained rebuild)
sys.path.insert(0, str(ROOT))
import plib  # noqa: E402
from plib import Node, mk, nhx, sha256  # noqa: E402

FAM = "P02_bHLH"
SRC = "PiresDolan2010"
SUPP = ROOT / "literature/raw/B2/supp/PMC2839125"
F009 = SUPP / "msp288_mbe-09-0678-File009.xls"
F011 = SUPP / "msp288_mbe-09-0678-File011.txt"
FIGS2 = SUPP / "msp288_2.pdf"
FULLTEXT = ROOT / "literature/papers/PMID19942615/fulltext.txt"
CACHE = HERE / "cache"
SPECIES = {"Arabidopsis thaliana": "ARATH", "Oryza sativa": "ORYSJ", "Physcomitrella patens": "PHYPA",
           "Selaginella moellendorffii": "SELMO", "Chlamydomonas reinhardtii": "CHLRE"}
PREFIX = {"At": "ARATH", "Os": "ORYSJ", "Pp": "PHYPA", "Sm": "SELMO", "Cr": "CHLRE"}


ASSERTIONS = [
    ("Abstract", "event", "D_subfamilies_Embryophyta", "T1",
     "Twenty of these subfamilies existed in the common ancestors of extant mosses and vascular plants, whereas six "
     "further subfamilies evolved among the vascular plants."),
    ("Results, Twenty bHLH Subfamilies Found in Flowering Plants Were also Present in Early Land Plants", "event",
     "D_subfamilies_Embryophyta", "T1",
     "20 of these subfamilies were already present in early land plants 443 Ma, by which time the mosses had "
     "diverged from the vascular plants."),
    ("Results, Twenty bHLH Subfamilies Found in Flowering Plants Were also Present in Early Land Plants", "event",
     "D_subfamilies_Tracheophyta", "T1",
     "Of the 26 plant bHLH subfamilies, 3 include only angiosperm proteins and 23 include angiosperm and lycophyte "
     "proteins (fig. 3)."),
    ("Results, Twenty bHLH Subfamilies Found in Flowering Plants Were also Present in Early Land Plants", "event",
     "D_subfamilies_Tracheophyta", "T1",
     "this implies that these 23 bHLH subfamilies are at least 415 million years (My) old"),
    ("Results, Twenty bHLH Subfamilies Found in Flowering Plants Were also Present in Early Land Plants", "event",
     "D_subfamilies_angiosperm", "T1",
     "Of the 26 plant bHLH subfamilies, 3 include only angiosperm proteins and 23 include angiosperm and lycophyte "
     "proteins (fig. 3)."),
    ("Abstract", "completeness", "ARATH,ORYSJ,SELMO", "T1",
     "In this study, we use whole-genome sequences from nine species of land plants and algae to define the "
     "relationships between these proteins in plants."),
]


def subfamilies() -> tuple[dict[str, str], dict[str, str]]:
    """protein -> subfamily, and subfamily -> clade-support string, from Supplementary Fig. S2 text."""
    from pypdf import PdfReader
    txt = PdfReader(str(FIGS2)).pages[0].extract_text()
    support = dict(re.findall(r"^([IVX]+[a-z]?(?:\([^)]*\))?(?:\([0-9]\))?) \(([^)]*;[^)]*)\)\s*$", txt, re.M))
    tail = txt.split("Subfamily Proteins", 1)[1].split("Supplementary Figure 2", 1)[0]
    tail = re.sub(r"bH\s+LH", "bHLH", tail)          # pypdf line-wrap artefacts ("OsbH LH045")
    tail = re.sub(r"\s+", " ", tail)
    parts = re.split(r"(?:^| )((?:[IVX]+[a-z]?(?:\([^)]*\))?)|Orphans) (?=[A-Z][a-z]bHLH)", " " + tail)
    out: dict[str, str] = {}
    for i in range(1, len(parts), 2):
        sf = parts[i]
        for p in re.findall(r"[A-Z][a-z]bHLH\d{3}", parts[i + 1]):
            if p in out:
                raise SystemExit(f"{p} assigned twice ({out[p]}, {sf})")
            out[p] = sf
    return out, support


def fasta(path) -> dict[str, str]:
    seqs, cur = {}, None
    for line in open(path, encoding="latin-1"):
        line = line.strip()
        if line.startswith(">"):
            cur = line[1:].split()[0]
            seqs[cur] = ""
        elif cur:
            seqs[cur] += line
    return seqs


def main() -> None:
    st = plib.species_tree()
    t = pd.read_excel(F009)
    t["Species"] = t["Species"].ffill()
    t["xls_row"] = range(2, len(t) + 2)
    sf_of, support = subfamilies()
    seqs = fasta(F011)
    report, members, info, queries = [], [], {}, {}
    for _, r in t.iterrows():
        prot = str(r["Protein"]).strip()
        code = SPECIES.get(r["Species"])
        locus = str(r["Locus name"]).strip() if pd.notna(r["Locus name"]) else ""
        sf = sf_of.get(prot, "")
        loc = f"File009 row {r['xls_row']}; subfamily: Suppl. Fig. S2 table" if sf else f"File009 row {r['xls_row']}"
        base = dict(original_id=prot, species=code or r["Species"], label=sf or "not in Fig. S2 table",
                    gene_id="", locator=loc, source=SRC)
        if code is None:
            report.append({**base, "status": "not_in_panel", "reason": "species not in panel P"})
            continue
        if not sf:
            report.append({**base, "status": "excluded", "reason": "not among the 544 plant proteins classified in Fig. S2"})
            continue
        info[prot] = dict(code=code, locus=locus, sf=sf, base=base)
        if prot in seqs:
            queries.setdefault(code, {})[prot] = seqs[prot]
    missing = sorted(set(sf_of) - set(t["Protein"].astype(str).str.strip()))
    if missing:
        print("WARNING: Fig. S2 proteins absent from File009:", missing)
    seqres = {}
    for code, q in sorted(queries.items()):
        seqres.update(plib.seqmap(code, q, CACHE / f"diamond_{code}.tsv"))
    cache = plib.JsonCache(CACHE / "uniprot_idmap.json")
    os_ids = sorted({(d["locus"], "gene_name") for d in info.values() if d["code"] == "ORYSJ" and d["locus"]})
    upmap = plib.uniprot_map(cache, os_ids, taxon="39947") if os_ids else {}
    ens = plib.JsonCache(CACHE / "ensembl_overlap.json")
    by_acc: dict[str, list[str]] = {}
    for prot, d in info.items():
        code, s = d["code"], seqres.get(prot)
        if s and s["status"] in ("exact", "unique"):
            d.update(acc=s["acc"], method=f"sequence_{s['status']}", ident=s["ident"], cov=f"{s['qcov']}/{s['tcov']}",
                     note=s["note"] + (f"; margin {s['margin']}" if s.get("margin") != "" else ""))
        else:
            hits = set(plib.id_index(code).get(d["locus"].upper(), set())) if d["locus"] else set()
            rapnote = ""
            if code == "ORYSJ" and d["locus"]:
                hits |= {x["uniprot_acc"] for x in upmap.get(d["locus"], []) if x.get("uniprot_acc") in plib.proteome(code)}
                if not hits:
                    rap, rapnote = plib.msu_to_rap(ens, d["locus"])
                    if rap:
                        hits = set(plib.id_index(code).get(rap.upper(), set()))
                        d["locus"] = f"{d['locus']} -> {rap}"
            seqnote = f"sequence: {s['status']} ({s.get('note', '')})" if s else "sequence: none"
            if len(hits) == 1:
                acc = next(iter(hits))
                d.update(acc=acc, method="exact_id", ident="", cov="",
                         note=f"exact identifier mapping {d['locus']}->{acc} (single panel hit); {seqnote}")
            elif hits:
                d.update(status="ambiguous", reason=f"{d['locus']} -> several panel proteins {sorted(hits)[:4]}; {seqnote}")
            else:
                d.update(status="unmapped", reason=f"{seqnote}; locus '{d['locus'] or 'none'}' gives no panel hit"
                         + (f" ({rapnote})" if rapnote else ""))
        if d.get("acc"):
            by_acc.setdefault(d["acc"], []).append(prot)
    for acc, ps in by_acc.items():
        if len(ps) > 1 and all(info[p]["method"] == "sequence_exact" for p in ps) and \
                len({plib.norm_seq(seqs[p]) for p in ps}) == 1:
            keep = sorted(ps)[0]
            for p in sorted(ps)[1:]:
                info[p].pop("acc")
                info[p].update(status="excluded", reason=(
                    f"identical sequence to {keep} (both exactly = panel protein {acc}); redundant entry for one panel "
                    f"gene (Selaginella v1 haplotype copy; protocol 0.6 identical copies) - {keep} kept"))
            continue
        if len(ps) > 1:
            for p in ps:
                info[p].pop("acc")
                info[p].update(status="ambiguous", reason=f"{len(ps)} paper proteins ({', '.join(ps)}) map to one panel protein {acc}")
    for prot, d in info.items():
        b = d["base"]
        if d.get("acc"):
            report.append({**b, "gene_id": d["acc"], "status": "mapped", "reason": f"{d['method']}: {d['note']}"})
            members.append(dict(gene_id=d["acc"], species=d["code"], original_id_in_paper=prot, label=d["sf"],
                                mapping_method=d["method"], mapping_identity=d["ident"], mapping_coverage=d["cov"],
                                source=SRC, locator=b["locator"]))
        else:
            report.append({**b, "status": d["status"], "reason": d["reason"]})
    mem = pd.DataFrame(members)

    # ---------------------------------------------------------------- subfamily ages from the paper's own table
    sp_by_sf: dict[str, set] = {}
    for p, sf in sf_of.items():
        sp_by_sf.setdefault(sf, set()).add(p[:2])
    land, vasc, angio = [], [], []
    for sf, sps in sorted(sp_by_sf.items()):
        if sf == "Orphans":
            continue
        if "Pp" in sps:
            land.append(sf)
        elif "Sm" in sps:
            vasc.append(sf)
        else:
            angio.append(sf)
    assert len(land) == 20 and len(vasc) + len(angio) == 6, (land, vasc, angio)

    T = dict(tier="T2", src=SRC)
    by_sf: dict[str, list] = {}
    for m in members:
        by_sf.setdefault(m["label"], []).append(m["gene_id"])

    def sfnode(sf):
        g = by_sf.get(sf, [])
        if not g:
            return None
        if len(g) == 1:
            return Node(g[0])
        n = Node(features=dict(ev="U", id="SF_" + re.sub(r"[^A-Za-z0-9]", "", sf), **T))
        for x in g:
            n.add_child(Node(x))
        return n

    d_land = mk("D", [sfnode(s) for s in land], id="D_subfamilies_Embryophyta", lvl_lo="Embryophyta",
                lvl_hi="ABOVE_RANGE", tier="T1", src=SRC)
    d_vasc = mk("D", [sfnode(s) for s in vasc], id="D_subfamilies_Tracheophyta", lvl_lo="Tracheophyta",
                lvl_hi="ABOVE_RANGE", tier="T1", src=SRC)
    d_ang = mk("D", [sfnode(s) for s in angio], id="D_subfamilies_angiosperm", lvl_lo="Mesangiospermae",
               lvl_hi="ABOVE_RANGE", tier="T1", src=SRC)
    root = mk("U", [d_land, d_vasc, d_ang] + [Node(x) for x in by_sf.get("Orphans", [])], id="bHLH_root", **T)
    (HERE / "reference_tree.nhx").write_text(nhx(root) + "\n")
    plib.write_assertions(HERE / "assertions.tsv", [
        dict(assertion_id=f"P02_A{i + 1}", type=a[1], target=a[2], tier=a[3], source_id=SRC, locator=a[0], quote=a[4])
        for i, a in enumerate(ASSERTIONS)], {SRC: FULLTEXT.read_text()})
    mem.to_csv(HERE / "members.tsv", sep="\t", index=False)
    rep = pd.DataFrame(report)[["original_id", "species", "label", "gene_id", "status", "reason", "locator", "source"]]
    rep.to_csv(HERE / "mapping_report.tsv", sep="\t", index=False)
    inp = rep[rep.status.isin(["mapped", "unmapped", "ambiguous"])]
    frac = inp.groupby("species").apply(lambda x: (x.status == "mapped").mean(), include_groups=False)
    complete = sorted(c for c in PREFIX.values() if frac.get(c, 0) >= 0.9)
    n_in, n_map = len(inp), int((inp.status == "mapped").sum())
    import yaml
    meta = {
        "id": FAM, "name": "Plant bHLH transcription factors (26 subfamilies)", "panel": "P",
        "phenomena": ["DOM", "DEEP"],
        "descriptors": {"n_members": int(len(mem)), "n_species": int(mem.species.nunique()),
                        "n_in_panel_paper_members": n_in, "mapped_fraction": round(n_map / n_in, 3),
                        "mapped_fraction_by_species": {k: round(float(v), 3) for k, v in frac.items()},
                        "subfamilies_land": land, "subfamilies_vascular_only": vasc, "subfamilies_angiosperm_only": angio,
                        "subfamily_support_ML_NJ_MP": support},
        "sources": [{
            "id": SRC, "doi": "10.1093/molbev/msp288", "pmid": "19942615", "pmcid": "PMC2839125",
            "title": "Origin and diversification of basic-helix-loop-helix proteins in plants", "year": 2010,
            "source_method": "phylogeny",
            "locator": "File009 (species/protein/locus), File011 (FASTA), Suppl. Fig. S2 subfamily table; Results",
            "quote": "Phylogenetic analyses suggest that plant bHLH proteins are monophyletic and constitute 26 subfamilies.",
            "files": [{"path": str(p.relative_to(ROOT)), "sha256": sha256(p)} for p in (F009, F011, FIGS2, FULLTEXT, plib.MSU_INFO)]}],
        "root_level": "Viridiplantae", "root_single_copy": False, "complete_species": complete,
        "decoys": [], "distinct_from": [],
        "tier_summary": {}, "protocol_version": 1,
        "notes": ("Subfamily membership from the Suppl. Fig. S2 table (pypdf text). Subfamily ages follow the paper: the 20 "
                  "subfamilies with moss members existed in the moss/vascular ancestor (D lo=Embryophyta), subfamilies "
                  "with lycophyte+angiosperm members only (lo=Tracheophyta) and angiosperm-only subfamilies "
                  "(lo=Mesangiospermae); all upper bounds ABOVE_RANGE because inter-subfamily relationships are "
                  "unresolved (paper: basal nodes have low support). Orphans and the chlorophyte proteins: U. Inside "
                  "subfamilies: U (paralogs; WGD lineages; no synteny curation). Pairs between the three age groups are U (the paper "
                  "does not identify the parent subfamily of the younger subfamilies)."),
    }
    txt = (HERE / "reference_tree.nhx").read_text()
    for ev, tier in re.findall(r"ev=([SDUH])[^\]]*?tier=(T\w+)", txt):
        meta["tier_summary"][f"{ev}_{tier}"] = meta["tier_summary"].get(f"{ev}_{tier}", 0) + 1
    (HERE / "family.yaml").write_text(yaml.safe_dump(meta, sort_keys=False, allow_unicode=True, width=120))
    print(f"{FAM}: {n_map}/{n_in} mapped ({n_map / n_in:.1%}); {frac.round(2).to_dict()}; complete={complete}; "
          f"land={len(land)} vasc={vasc} angio={angio}")


if __name__ == "__main__":
    main()
