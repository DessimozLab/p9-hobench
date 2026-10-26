#!/usr/bin/env python
"""P01_homeobox: plant homeobox superfamily (Mukherjee, Brocchieri & Buerglin 2009, MBE 26:2775).

Source data: supplementary table 1 (msp201_5.html: gene name, accession, alternative names, species code,
gene class, protein sequence). Members of panel-P species are mapped onto the pinned UniProt proteomes:
(1) same-species sequence cross-reference (exact, or >=98 % identity, >=80 % coverage of both, RBH, >=2-point
margin), else (2) exact identifier mapping (SELMODRAFT_/OSTLU_ JGI protein ids, Arabidopsis AGI codes, rice
MSU LOC ids, GenBank protein accessions via UniProt ID mapping), single hit in the panel proteome only.

Tree (T1 = class-origin statements quoted in assertions.tsv; T2 = class labels):
  root  D  TALE vs non-TALE (pre-Viridiplantae)
  TALE  U[ D(KNOX, BEL) (pre-Viridiplantae), unclassified TALE genes ]
        KNOX: U[ D(KNOX I, KNOX II) (Embryophyta..ABOVE_RANGE), other KNOX ]
  non-TALE U[ D(DDT, WOX, PINTOX) (pre-Viridiplantae),
              D(HD-ZIP I, II, III, IV) (Embryophyta..Streptophyta), PLINC, PHD, NDX, LD, SAWADEE ]
  DDT / WOX / PINTOX: S at Viridiplantae (chlorophyte gene vs land-plant genes) where both exist
  Inside a class: several genes per species and WGD-descended lineages -> U (protocol 3, rule 2b/2c).
"""
from __future__ import annotations

import csv
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

FAM = "P01_homeobox"
SRC = "Mukherjee2009"
TABLE = ROOT / "literature/raw/B2/supp/PMC2775110/msp201_5.html"
FULLTEXT = ROOT / "literature/papers/PMID19734295/fulltext.txt"
CACHE = HERE / "cache"

SPECIES = {"At": "ARATH", "Os": "ORYSJ", "Os_j": "ORYSJ", "Zm": "MAIZE", "Pp": "PHYPA", "Sm": "SELMO",
           "Cr": "CHLRE", "Ol": "OSTLU", "Vv": "VITVI", "Gm": "SOYBN", "Br": "BRARP", "Pv": "PHAVU"}
# species whose genomes the paper searched exhaustively (Table 1 / Abstract)
GENOME_WIDE = {"ARATH", "ORYSJ", "MAIZE", "PHYPA", "SELMO", "CHLRE", "OSTLU"}

EXCLUDE_NO_HD = {"PLINC ONLY", "SAWADEE ONLY", "SAWDEE ONLY", "START ONLY", "MEKHLA ONLY", "PH - START",
                 "RhoGAP - START", "TM - START", "START - HD-SAD", "START WITHOUT HOMEODOMAIN",
                 "PH + START WITHOUT HD", "no HD, other transcription factor", "SAWADEE - like",
                 "SAWADEE - like (SNF2?)", "DDT_WSD motif"}


def norm_class(lbl: str):
    """-> (class, subclass) or (None, reason)."""
    l = str(lbl).strip()
    if "pseudo" in l.lower():
        return None, "pseudogene (protocol 0.6)"
    if l in EXCLUDE_NO_HD:
        return None, f"not a homeobox gene in the paper's classification (label '{l}': codomain without HD / not counted in Table 2)"
    u = l.upper().replace("_", "-")
    if u.startswith("HD-ZIP III"):
        return "HD-ZIP III", ""
    if u.startswith("HD-ZIP IV"):
        return "HD-ZIP IV", ""
    if u.startswith("HD-ZIP 2"):
        return "HD-ZIP II", ""
    if u.startswith("HD-ZIP 1"):
        return "HD-ZIP I", ""
    if u in ("KNOX I", "KNOX II"):
        return "KNOX", u
    if u in ("KNOX", "TALE / KNOX"):
        return "KNOX", ""
    if u in ("BEL", "TALE / BEL"):
        return "BEL", ""
    if u == "TALE":
        return "TALE_unclassified", ""
    if u in ("DDT", "D-TOX A, NO DDT", "DDT (D-TOX A, NO DDT)"):
        return "DDT", ""
    if u in ("PHD ZINC FINGER", "PHD"):
        return "PHD", ""
    if u in ("PLINC", "WOX", "NDX", "LD", "PINTOX", "SAWADEE"):
        return u, ""
    return None, f"unrecognised class label '{l}'"


def load_table() -> pd.DataFrame:
    t = pd.read_html(TABLE)[0]
    t["row"] = range(1, len(t) + 1)
    return t


ID_GENE = re.compile(r"\b(At[1-5CM]g\d{5})\b", re.I)
ID_MSU = re.compile(r"\bOs(\d\d)g(\d{5})\b")
ID_RAP = re.compile(r"\b(Os\d\dg\d{7})\b")
ID_GB = re.compile(r"^(?:[A-Z]{3}\d{5}|[NX]P_\d{6,9})(?:\.\d+)?$")


def candidate_ids(r, sp: str) -> list[tuple[str, str]]:
    """(identifier, kind) usable for exact identifier mapping, from Gene Name / Accession / Alternative Names."""
    gname, acc, alt = (str(r[c]) if pd.notna(r[c]) else "" for c in ("Gene Name", "Accession No.", "Alternative Names"))
    blob = " ".join([gname, acc, alt])
    out: list[tuple[str, str]] = []
    if sp == "SELMO":
        m = re.search(r"jgi\|Selmo1\|(\d+)\|", alt)
        n = m.group(1) if m else (gname if re.fullmatch(r"\d{4,7}", gname) else "")
        if n:
            out.append((f"SELMODRAFT_{n}", "jgi_gene_name"))
    elif sp == "OSTLU":
        m = re.search(r"OSTLU_(\d+)|jgi_Ost9901_3_(\d+)_", alt)
        n = (m.group(1) or m.group(2)) if m else (gname if re.fullmatch(r"\d{3,6}", gname) else "")
        if n:
            out.append((f"OSTLU_{n}", "jgi_gene_name"))
    elif sp == "ARATH":
        for g in sorted({x.upper() for x in ID_GENE.findall(blob)}):
            out.append((g, "agi"))
    elif sp == "ORYSJ":
        for a, b in sorted(set(ID_MSU.findall(blob))):
            out.append((f"LOC_Os{a}g{b}", "msu_loc"))
        for g in sorted(set(ID_RAP.findall(blob))):
            out.append((g, "rap_locus"))
    if ID_GB.match(acc.strip()):
        out.append((acc.strip(), "refseq_protein" if "_" in acc else "genbank_protein"))
    return out


def id_lookup(code: str):
    """gene-name / locus index of a panel proteome: NAME(upper) -> set(acc)."""
    idx: dict[str, set] = {}
    for acc, h in plib.header_map(code).items():
        names = set()
        for k in ("gene_name", "gene_symbol_gene2acc", "canonical_of_gene_group"):
            for v in re.split(r"[;,\s]+", h.get(k) or ""):
                if v:
                    names.add(v.split(":", 1)[-1].upper())
        for m in re.finditer(r"(?:GN|OrderedLocusNames|ORFNames)=(\S+)", h.get("original_header", "")):
            names.add(m.group(1).upper())
        for v in names:
            idx.setdefault(v, set()).add(acc)
    return idx


# (locator, type, target node id, tier, quote) -- hand-entered T1 statements (protocol 0.4), checked against fulltext.txt
ASSERTIONS = [
    ("Discussion, HD Classes and Domains Common to All Examined Plant Groups", "event", "D_TALE_nonTALE", "T1",
     "suggesting that TALE and non-TALE HD proteins were already present in the common ancestor of plants and animals"),
    ("Abstract", "event", "D_KNOX_BEL", "T1",
     "Members of the BEL and KNOX classes were found in Chlorobionta (green plants) and in Rhodophyta."),
    ("Discussion, HD Classes and Domains Common to All Examined Plant Groups", "event", "D_KNOX_BEL", "T1",
     "In either case, KNOX and BEL HD proteins must have been already present in the last common ancestor of the plant clade."),
    ("Discussion, Appearance of Three New HD Classes in the Common Ancestor Lineage of Green Plants", "event",
     "D_DDT_WOX_PINTOX", "T1",
     "identifying the presence of the homeodomain protein DDT class in the last common ancestor of green plants (Chlorobionta)"),
    ("Discussion, Appearance of Three New HD Classes in the Common Ancestor Lineage of Green Plants", "event",
     "D_DDT_WOX_PINTOX", "T1",
     "We also infer that at this time HD was also associated with the green plant\u2013specific PINTOX and WUS motifs, "
     "giving rise to the corresponding PINTOX and WOX classes."),
    ("Results, The WOX Class", "event", "D_DDT_WOX_PINTOX", "T1",
     "indicating that the WOX class of proteins was present in the last common ancestors of green plants (Chlorobionta)"),
    ("Discussion, Multiple HD Genes and HD Classes in Embryophyta (Land Plants)", "event", "D_HDZIP_I_IV", "T1",
     "Embryophyta are characterized by the appearance, at least in the lineage leading to moss and vascular plants, of "
     "nine classes of HD proteins not present in Rhodophyta or Chlorophyta, namely HD-ZIP I to IV, LD, SAWADEE, NDX, PHD, and PLINC."),
    ("Discussion, Multiple HD Genes and HD Classes in Embryophyta (Land Plants)", "event", "D_HDZIP_I_IV", "T1",
     "suggesting that they originated early in land plant evolution (at least before separation of moss and vascular "
     "plants and possibly before transition to land)"),
    ("Results, The TALE Superclass (KNOX and BEL)", "event", "D_KNOXI_KNOXII", "T1",
     "suggesting that both KNOX I and KNOX II genes were present in the common ancestor of moss and vascular plants"),
    ("Discussion, Appearance of Three New HD Classes in the Common Ancestor Lineage of Green Plants", "event",
     "DDT_Viridiplantae_S", "T2",
     "Presumably each of these classes was originally represented by one gene, as observed in modern-day unicellular green algae."),
    ("Discussion, Appearance of Three New HD Classes in the Common Ancestor Lineage of Green Plants", "event",
     "WOX_Viridiplantae_S", "T2",
     "Presumably each of these classes was originally represented by one gene, as observed in modern-day unicellular green algae."),
    ("Results, The WOX Class", "tag", "WOX_Viridiplantae_S", "T1",
     "was successively lost in the Chlamydomonas lineage"),
    ("Abstract", "completeness", "", "T1",
     "The full complement of homeobox transcription factor sequences, including genes and pseudogenes, was determined "
     "from the analysis of 10 complete genomes"),
]


def main() -> None:
    t = load_table()
    st = plib.species_tree()
    cache = plib.JsonCache(CACHE / "uniprot_idmap.json")
    report, members = [], []
    queries: dict[str, dict[str, str]] = {}
    info: dict[str, dict] = {}
    for _, r in t.iterrows():
        code = SPECIES.get(r["Species Code"])
        rid = f"row{r['row']}"
        orig = str(r["Gene Name"]) if pd.notna(r["Gene Name"]) else rid
        cls, sub = norm_class(r["Gene Classes"])
        loc = f"msp201_5.html row {r['row']} (Gene Number {r['Gene Number']})"
        base = dict(original_id=orig, species=code or str(r["Species Code"]), label=str(r["Gene Classes"]),
                    gene_id="", locator=loc, source=SRC)
        if code is None:
            report.append({**base, "status": "not_in_panel", "reason": "species not in panel P"})
            continue
        if cls is None:
            report.append({**base, "status": "excluded", "reason": sub})
            continue
        info[rid] = dict(row=r, code=code, cls=cls, sub=sub, base=base)
        seq = r["Protein Sequence"]
        if isinstance(seq, str) and plib.norm_seq(seq):
            queries.setdefault(code, {})[rid] = seq
    seqres: dict[str, dict] = {}
    for code, q in sorted(queries.items()):
        seqres.update(plib.seqmap(code, q, CACHE / f"diamond_{code}.tsv"))

    # identifier route for rows not mapped by sequence
    need_uniprot = []
    for rid, d in info.items():
        d["ids"] = candidate_ids(d["row"], d["code"])
        for i, kind in d["ids"]:
            if kind in ("genbank_protein", "refseq_protein", "msu_loc"):
                need_uniprot.append((i, kind if kind != "msu_loc" else "gene_name"))
    taxa = {c: plib.species_table()[c]["taxid"] for c in SPECIES.values()}
    upmap = {}
    for code in sorted(set(SPECIES.values())):
        ids = sorted({(i, k) for rid, d in info.items() if d["code"] == code
                      for i, k0 in d["ids"] for k in [k0 if k0 != "msu_loc" else "gene_name"]
                      if k0 in ("genbank_protein", "refseq_protein", "msu_loc")})
        if ids:
            upmap.update({(code, i): v for i, v in plib.uniprot_map(cache, ids, taxon=taxa[code]).items()})
    idx_cache: dict[str, dict] = {}
    rows_by_acc: dict[str, list[str]] = {}
    for rid, d in sorted(info.items(), key=lambda x: int(x[0][3:])):
        code = d["code"]
        s = seqres.get(rid)
        acc, method, ident, cov, note = "", "", "", "", ""
        if s and s["status"] in ("exact", "unique"):
            acc, method, ident = s["acc"], f"sequence_{s['status']}", s["ident"]
            cov = f"{s['qcov']}/{s['tcov']}"
            note = s["note"] + (f"; margin {s['margin']}" if s.get("margin") != "" else "")
        else:
            hits = set()
            used = []
            if code not in idx_cache:
                idx_cache[code] = id_lookup(code)
            for i, kind in d["ids"]:
                if kind in ("jgi_gene_name", "agi", "rap_locus"):
                    h = idx_cache[code].get(i.upper(), set())
                elif kind == "msu_loc":
                    h = idx_cache[code].get(i.upper(), set()) | {x["uniprot_acc"] for x in upmap.get((code, i), [])
                                                                  if x.get("uniprot_acc") in plib.proteome(code)}
                else:
                    h = {x["uniprot_acc"] for x in upmap.get((code, i), []) if x.get("uniprot_acc") in plib.proteome(code)}
                if h:
                    used.append(f"{i}->{','.join(sorted(h))}")
                hits |= h
            seqnote = f"sequence: {s['status']} ({s.get('note', '')})" if s else "sequence: none in table"
            if len(hits) == 1:
                acc, method = hits.pop(), "exact_id"
                note = f"exact identifier mapping {'; '.join(used)} (single hit in panel proteome); {seqnote}"
                cand = (s or {}).get("acc_candidate") or (s or {}).get("acc")
                if cand and cand != acc:
                    note += f" [sequence best hit {cand} differs]"
            elif len(hits) > 1:
                d["status"], d["reason"] = "ambiguous", f"identifier maps to several panel proteins: {'; '.join(used)}; {seqnote}"
            else:
                d["status"] = "unmapped"
                d["reason"] = f"{seqnote}; identifiers {[i for i, _ in d['ids']] or 'none usable'} give no panel hit"
        if acc:
            d.update(acc=acc, method=method, ident=ident, cov=cov, note=note)
            rows_by_acc.setdefault(acc, []).append(rid)
    # one gene per accession: several paper rows on one UniProt entry -> ambiguous (all of them)
    for acc, rids in rows_by_acc.items():
        if len(rids) > 1:
            for rid in rids:
                d = info[rid]
                d["status"], d["reason"] = "ambiguous", f"{len(rids)} paper rows ({', '.join(rids)}) map to the same panel protein {acc}"
                d.pop("acc")
    for rid, d in sorted(info.items(), key=lambda x: int(x[0][3:])):
        b = d["base"]
        lab = d["cls"] + (f"|{d['sub']}" if d["sub"] else "")
        if d.get("acc"):
            report.append({**b, "gene_id": d["acc"], "status": "mapped", "reason": f"{d['method']}: {d['note']}"})
            members.append(dict(gene_id=d["acc"], species=d["code"], original_id_in_paper=b["original_id"], label=lab,
                                mapping_method=d["method"], mapping_identity=d["ident"], mapping_coverage=d["cov"],
                                source=SRC, locator=b["locator"]))
        else:
            report.append({**b, "status": d["status"], "reason": d["reason"]})
    mem = pd.DataFrame(members)

    # ------------------------------------------------------------------ tree
    T = dict(tier="T2", src=SRC)
    by_cls: dict[str, list[tuple[str, str, str]]] = {}
    for m in members:
        c, _, sub = m["label"].partition("|")
        by_cls.setdefault(c, []).append((m["gene_id"], m["species"], sub))
    CHLORO = {"CHLRE", "OSTLU"}

    def u(genes, nid=None):
        if not genes:
            return None
        if len(genes) == 1:
            return Node(genes[0][0])
        feats = dict(ev="U", **T)
        if nid:
            feats["id"] = nid
        n = Node(features=feats)
        for g in genes:
            n.add_child(Node(g[0]))
        return n

    def green_class(c):
        """S at Viridiplantae between chlorophyte and streptophyte genes (classes present as one gene in the
        green-plant ancestor); each side U unless single."""
        gs = by_cls.get(c, [])
        ch = [g for g in gs if g[1] in CHLORO]
        la = [g for g in gs if g[1] not in CHLORO]
        chn = (plib.species_s_tree(st, [(g[0], g[1]) for g in ch], feats=dict(T))
               if ch and len({g[1] for g in ch}) == len(ch) else u(ch))
        if chn is not None and not chn.is_leaf:
            chn.features["id"] = f"{c}_Chlorophyta"
        lan = u(la, f"{c}_Streptophyta")
        if chn is not None and lan is not None:
            return mk("S", [chn, lan], id=f"{c}_Viridiplantae_S", tier="T2", src=SRC)
        return chn if chn is not None else lan

    knox = by_cls.get("KNOX", [])
    k1 = [g for g in knox if g[2] == "KNOX I"]
    k2 = [g for g in knox if g[2] == "KNOX II"]
    k0 = [g for g in knox if not g[2]]
    knox_d = mk("D", [u(k1, "KNOX_I"), u(k2, "KNOX_II")], id="D_KNOXI_KNOXII", lvl_lo="Embryophyta",
                lvl_hi="ABOVE_RANGE", tier="T1", src=SRC) if k1 and k2 else None
    knox_n = mk("U", [knox_d] + [Node(g[0]) for g in (k0 + ([] if knox_d else k1 + k2))], id="KNOX", **T)
    knox_bel = mk("D", [knox_n, u(by_cls.get("BEL", []), "BEL")], id="D_KNOX_BEL", lvl_lo="Viridiplantae",
                  lvl_hi="ABOVE_RANGE", tier="T1", src=SRC)
    tale = mk("U", [knox_bel] + [Node(g[0]) for g in by_cls.get("TALE_unclassified", [])], id="TALE", **T)
    green = mk("D", [green_class("DDT"), green_class("WOX"), green_class("PINTOX")], id="D_DDT_WOX_PINTOX",
               lvl_lo="Viridiplantae", lvl_hi="ABOVE_RANGE", tier="T1", src=SRC)
    hdzip = mk("D", [u(by_cls.get(c, []), c.replace(" ", "_")) for c in ("HD-ZIP I", "HD-ZIP II", "HD-ZIP III", "HD-ZIP IV")],
               id="D_HDZIP_I_IV", lvl_lo="Embryophyta", lvl_hi="Streptophyta", tier="T1", src=SRC)
    land = [u(by_cls.get(c, []), c) for c in ("PLINC", "PHD", "NDX", "LD", "SAWADEE")]
    nontale = mk("U", [green, hdzip] + land, id="nonTALE", **T)
    root = mk("D", [tale, nontale], id="D_TALE_nonTALE", lvl_lo="Viridiplantae", lvl_hi="ABOVE_RANGE",
              tier="T1", src=SRC)
    (HERE / "reference_tree.nhx").write_text(nhx(root) + "\n")
    node_ids = set(re.findall(r"id=([^:\]]+)", nhx(root)))
    rows = [a for a in ASSERTIONS if a[2] in node_ids or a[1] == "completeness"]
    plib.write_assertions(HERE / "assertions.tsv", [
        dict(assertion_id=f"P01_A{i + 1}", type=a[1], target=a[2] if a[1] != "completeness" else ",".join(
            sorted(GENOME_WIDE)), tier=a[3], source_id=SRC, locator=a[0], quote=a[4]) for i, a in enumerate(rows)],
        {SRC: FULLTEXT.read_text()})

    # ------------------------------------------------------------------ tables
    mem.to_csv(HERE / "members.tsv", sep="\t", index=False)
    rep = pd.DataFrame(report)[["original_id", "species", "label", "gene_id", "status", "reason", "locator", "source"]]
    rep.to_csv(HERE / "mapping_report.tsv", sep="\t", index=False)

    # complete species: genome-wide search AND >= 90 % of in-panel (non-excluded) members uniquely mapped
    inp = rep[rep.status.isin(["mapped", "unmapped", "ambiguous"])]
    frac = inp.groupby("species").apply(lambda x: (x.status == "mapped").mean(), include_groups=False)
    complete = sorted(c for c in GENOME_WIDE if frac.get(c, 0) >= 0.9)
    n_in, n_map = len(inp), int((inp.status == "mapped").sum())
    write_yaml(mem, complete, frac.round(3).to_dict(), n_in, n_map, st)
    print(f"{FAM}: {n_map}/{n_in} in-panel members mapped ({n_map / n_in:.1%}); per species {frac.round(2).to_dict()}; "
          f"complete={complete}")


def write_yaml(mem, complete, frac, n_in, n_map, st):
    import yaml
    files = [{"path": str(p.relative_to(ROOT)), "sha256": sha256(p)} for p in (TABLE, FULLTEXT)]
    meta = {
        "id": FAM,
        "name": "Plant homeobox superfamily (14 classes: KNOX, BEL, HD-ZIP I-IV, WOX, PLINC, PHD, DDT, NDX, LD, SAWADEE, PINTOX)",
        "panel": "P",
        "phenomena": ["DOM", "LOSS", "DEEP"],
        "descriptors": {"n_members": int(len(mem)), "n_species": int(mem.species.nunique()),
                        "n_in_panel_paper_members": n_in, "mapped_fraction": round(n_map / n_in, 3),
                        "mapped_fraction_by_species": {k: float(v) for k, v in frac.items()},
                        "classes": sorted(mem.label.str.split("|").str[0].unique().tolist())},
        "sources": [{
            "id": SRC, "doi": "10.1093/molbev/msp201", "pmid": "19734295", "pmcid": "PMC2775110",
            "title": "A comprehensive classification and evolutionary analysis of plant homeobox genes",
            "year": 2009, "source_method": "phylogeny",
            "locator": "Supplementary table 1 (msp201_5.html); Results; Discussion 'Evolutionary History of Plant Homeobox Genes'",
            "quote": "All homeobox genes can be unambiguously classified by sequence evolutionary analysis into 14 distinct classes",
            "files": files}],
        "root_level": "Viridiplantae",
        "root_single_copy": False,
        "complete_species": complete,
        "decoys": [],
        "distinct_from": [],
        "tier_summary": tier_summary(),
        "protocol_version": 1,
        "notes": ("Classes are the paper's labels (Supplementary table 1). Inter-class relations are only asserted where "
                  "the paper dates the class origins (TALE vs non-TALE and KNOX vs BEL before plants; DDT, WOX, PINTOX in "
                  "the green-plant ancestor; HD-ZIP I-IV after the Chlorophyta/Streptophyta split and before the "
                  "moss/vascular split; KNOX I vs II before the moss/vascular split). Land-plant classes PLINC/PHD/NDX/LD/"
                  "SAWADEE relative to each other and to DDT/WOX/PINTOX: U (parent class not identified; the paper says "
                  "inter-class topology is poorly supported). Within classes: U (many paralogs; WGD-descended "
                  "lineages, no synteny curation), except chlorophyte-vs-land S nodes for DDT/WOX/PINTOX. Rows labelled "
                  "with codomains but no homeodomain (START ONLY, PLINC ONLY, MEKHLA ONLY, ...) and pseudogenes are "
                  "excluded. ORYSJ includes rows coded 'Os' and 'Os_j' (O. sativa; the searched genome is japonica). "
                  "Non-genome-wide species (VITVI, SOYBN, BRARP, PHAVU) contribute EST/cDNA-based rows only."),
    }
    (HERE / "family.yaml").write_text(yaml.safe_dump(meta, sort_keys=False, allow_unicode=True, width=120))


def tier_summary():
    txt = (HERE / "reference_tree.nhx").read_text()
    out = {}
    for ev, tier in re.findall(r"ev=([SDUH])[^\]]*?tier=(T\w+)", txt):
        out[f"{ev}_{tier}"] = out.get(f"{ev}_{tier}", 0) + 1
    return out


if __name__ == "__main__":
    main()
