#!/usr/bin/env python
"""V02_opsins: vertebrate visual opsins (LWS, SWS1, SWS2, RH2, RH1) from Lagman et al. 2013
(BMC Evol Biol 13:238, PMID 24180662, PMC3826523), Additional file 1 (S1.xlsx, sheet 'OPN':
species x subtype with Ensembl gene / NCBI / UniProt accessions) + text (2R scenario, local duplicates,
RH1 retrotransposition).

Tree (all quotes verified against fulltext.txt):
  root   D  local duplication before 2R giving ancestral SWS and LWS (lvl_lo=Vertebrata, lvl_hi=ABOVE_RANGE)
  SWS    D  2R -> SWS1, SWS2, RH1, RH2 (interval lvl_lo=Gnathostomata, lvl_hi=Vertebrata; Lagman: before
            lamprey divergence)
  RH1    S(sarcopterygian RH1, D[Neopterygii, retrotransposition](exorh clade, rho clade))
  subtypes otherwise S over species; within-species local duplicates: gar SWS1-1/-2 D(LEPOC,TAND);
           primate LWS/MWS D(HUMAN, TAND); teleost local duplicates U (sharing across teleosts not stated);
           zebrafish rho/rhol + medaka RH1-2 U ("cannot say for certain ... 3R").
Pinopsin / V/A opsin / Ci-opsin rows (sheet section 'ADDITIONAL PUTATIVE ...') and OPN3 (outgroup) are not
members (Lagman: synteny does not place them in the visual-opsin duplications).
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

import pandas as pd
import yaml

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))  # vlib.py: shared Panel V helpers (identical copy in each V0x family)
import vlib  # noqa: E402
from vlib import mk, species_s_tree  # noqa: E402

FID = "V02_opsins"
SRC = "Lagman2013"
PAPER = vlib.ROOT / "literature/papers/PMID24180662"
S1 = vlib.ROOT / "literature/raw/B1/supp/opsin_Lagman2013__PMC3826523/1471-2148-13-238-S1.xlsx"
TEXT = (PAPER / "fulltext.txt").read_text()
RAMIREZ = vlib.ROOT / "literature/raw/B1/supp/opsin_Ramirez2016__PMC5521729/evw248_Supp.zip"

SPECIES = {"Human": "HUMAN", "Mouse": "MOUSE", "Opossum": None, "Chicken": "CHICK", "Anole lizard": None,
           "Western clawed frog": "XENTR", "Coelacanth": "LATCH", "Spotted gar": "LEPOC", "Zebrafish": "DANRE",
           "Three-spined stickleback": None, "Medaka": "ORYLA", "Green spoted pufferfish": None,
           "Pouched lamprey": None}
TELEOSTS = {"DANRE", "ORYLA"}

Q = {
    "local_pre2R": "a local duplication occurred before 2R giving rise to two ancestral genes, SWS and LWS",
    "2R": "These genes later duplicated in 2R as part of a larger chromosomal block so that the ancestral SWS "
          "gene gave rise to the types SWS1, SWS2, RH1 and RH2 and the ancestral LWS gave rise to four copies "
          "out of which only one has been retained, namely the LWS type (Figure 5).",
    "lamprey": "In pouched lamprey all five visual opsin genes have previously been identified, suggesting that "
               "lampreys diverged from the jawed vertebrates after 2R.",
    "primate_lws": "Primates have a duplicate of LWS that has mutated to have its absorption maximum in green, "
                   "namely OPN1MW.",
    "gar_sws1": "The two SWS1 genes are located adjacent to each other on the same linkage group (LG8) "
                "approximately 8 Kb apart (see Table 1), and are thus most likely the result of a local duplication.",
    "retro": "Our analyses indicate that the retrotransposition event occurred before the divergence of holostean "
             "fishes (including gars) and teleosts (Figure 6).",
    "rh1_teleost": "We could identify these duplicate RH1 sequences in all investigated teleost genomes (RH1-1 and "
                   "RH1-2 in Table 1 and Figure 1).",
    "rho_gar": "the rhodopsin gene of the spotted gar clusters together with the teleost rhodopsin genes with high "
               "support, forming a well-defined actinopterygian cluster (marked rho in Figure 1).",
    "rhol_unsure": "While it is a likely scenario, we cannot say for certain that the duplicate intron-less RH1 "
                   "genes originated in 3R.",
    "teleost_local": "Additionally, LWS, SWS2, RH1 and RH2 sequences often occur as multiple local duplicates in "
                     "teleosts (Table 1).",
    "subtype_clusters": "In the phylogenetic analyses these sequences form five well-supported clusters: The LWS "
                        "cluster forms a basal separate branch while the SWS1, SWS2, RH2 and RH1 cluster together "
                        "with high support (Figure 1).",
    "va_pin": "the present synteny analysis does not indicate that pinopsins and V/A opsins arose in the same "
              "chromosome duplications that gave rise to the visual opsin gene family.",
    "mammal_loss": "the ancestral repertoire in jawed vertebrates consisted of five visual opsin genes with losses "
                   "of SWS2 and RH2 in mammals.",
}
for k, q in Q.items():
    assert vlib.quote_in(TEXT, q), f"quote {k} not found"

# chromosome name for gar loci (LepOcu1 = current Ensembl gar assembly)
GAR_LEN = 2000  # probe window from the given gene start


def parse_s1():
    v = pd.read_excel(S1, sheet_name="OPN", header=None, dtype=str).fillna("")
    recs, sp, section = [], None, "main"
    for i in range(2, len(v)):
        r = [str(c).strip() for c in v.iloc[i]]
        if r[0].startswith("ADDITIONAL"):
            section = "additional"
            continue
        if r[1] == "Species":
            continue
        if r[1]:
            lab = re.sub(r"\s*\(.*$", "", r[1]).strip()
            if lab in SPECIES:
                sp = lab
        if not r[3]:
            continue
        recs.append(dict(row_xls=i + 1, section=section, species_label=sp, symbol=r[2], subtype=r[3],
                         seqname=r[4], chrom=r[5], start=r[6], gene_id_paper=r[7], transcript=r[8], comment=r[9]))
    return recs


def main():
    st = vlib.species_tree()
    cache = vlib.JsonCache(HERE / "cache" / "lookups.json")
    recs = parse_s1()
    for x in recs:
        code = SPECIES.get(x["species_label"])
        x["species"] = code
        x.update(gene_id="", status="", reason="", mapping_method="", mapping_identity="", mapping_coverage="")
        oid = re.sub(r"^GenBank:\s*", "", x["gene_id_paper"]).strip()
        x["original_id"] = oid if oid not in ("-", "") else (x["transcript"] or x["seqname"])
        if x["seqname"].startswith("Cin "):
            x["species_label"], code = "Ciona intestinalis", "CIOIN"
            x["species"] = code
        if code is None:
            x.update(status="not_in_panel", reason=f"{x['species_label']} not in panel V")
            continue
        if x["section"] != "main" or x["subtype"] == "OPN3":
            x.update(status="excluded", reason="not a visual opsin member (OPN3 outgroup / pinopsin / V/A / "
                                                "Ci-opsin rows of the 'additional' section)")
            continue
        if "Possible pseudogene" in x["comment"]:
            x.update(status="excluded", reason="pseudogene per S1 comment")
            continue
        region = None
        if code == "LEPOC" and x["chrom"].startswith("LG") and x["start"].isdigit():
            s0 = int(x["start"])
            region = f"{x['chrom']}:{s0}-{s0 + 1}"
            oid = ""  # gar IDs in S1 are stickleback/zebrafish/coelacanth query proteins, not gar genes
            x["original_id"] = f"{x['seqname']} ({x['chrom']}:{s0})"
        res = vlib.map_one(cache, code, oid, region)
        alts = re.findall(r"NCBI sequence ([A-Z]{1,2}_?\d+\.?\d*)", x["comment"])
        if x["transcript"].startswith("ENS"):
            alts.append(x["transcript"])
        for alt in alts:  # S1 comment: paper used this NCBI sequence; else the S1 Ensembl transcript ID
            if res["status"] == "mapped":
                break
            r2 = vlib.map_one(cache, code, alt)
            if r2["status"] == "mapped" or r2.get("candidates"):
                r2["reason"] = f"via S1 alternative ID {alt}; " + r2["reason"]
                r2["candidates"] = res.get("candidates", []) + r2.get("candidates", [])
                res = r2
        x.update(res)
    vlib.resolve_by_sequence  # noqa
    for code in sorted({x["species"] for x in recs if x["species"]}):
        vlib.resolve_by_sequence(cache, code, [x for x in recs if x["species"] == code and x["status"] != "excluded"],
                                 strict_family=True)  # TAND-tagged family -> 99 % identity
    # identical copies / two S1 rows on one accession -> ambiguous (protocol: exclude identical copies)
    by: dict[str, list] = {}
    for x in recs:
        if x["status"] == "mapped":
            by.setdefault(x["gene_id"], []).append(x)
    for x in recs:
        if x["status"] in ("ambiguous", "unmapped") and "panel accessions" in x["reason"]:
            accs = re.findall(r"'([A-Z0-9]+)'", x["reason"])
            seqs = {vlib.proteome(x["species"]).get(a) for a in accs}
            if len(accs) > 1 and len(seqs) == 1:
                x.update(status="excluded", reason="identical copies in panel proteome (protocol 6): " + x["reason"])
    for acc, xs in by.items():
        if len(xs) > 1:
            for x in xs:
                x.update(status="ambiguous", gene_id="",
                         reason=f"{len(xs)} S1 rows map to {acc} (same gene ID in S1 or identical copies)")
    mapped = [x for x in recs if x["status"] == "mapped"]
    for x in mapped:
        x["label"] = x["subtype"] + (f":{x['symbol']}" if x["symbol"] else "") + f" [{x['seqname']}]"

    # ---------------------------------------------------------------- tree
    assertions = []
    nid = iter(range(1, 1000))

    def A(node_id, typ, target, qkey, loc, tier="T1"):
        assertions.append(dict(assertion_id=f"{FID}.A{len(assertions) + 1:03d}", type=typ, node_id=node_id,
                               target=node_id, note=target, tier=tier, source_id=SRC, locator=loc, quote=Q[qkey]))

    def genes(sub, pred=lambda x: True):
        return [(x["gene_id"], x["species"]) for x in mapped if x["subtype"] == sub and pred(x)]

    T2 = {"tier": "T2", "src": SRC}

    def subtype_tree(sub):
        """S over non-teleosts; teleost block: S over species, several copies per teleost species -> U.
        Species with 2 copies outside teleosts handled by caller."""
        non = genes(sub, lambda x: x["species"] not in TELEOSTS)
        tel = genes(sub, lambda x: x["species"] in TELEOSTS)
        tel_counts = pd.Series([s for _, s in tel]).value_counts() if tel else pd.Series(dtype=int)
        if tel and (tel_counts > 1).any():
            i = f"{FID}.U{next(nid)}"
            teleo = mk("U", [Node(g) for g, _ in tel], id=i, tier="T1", src=SRC)  # local duplicates; sharing unknown
            A(i, "event", f"{sub}: teleost local duplicates left unresolved (U)", "teleost_local", "Results p.1",
              tier="T1")
        else:
            teleo = species_s_tree(st, tel, T2) if tel else None
        return attach(non, teleo)

    def attach(non, teleo, inparalog=None):
        if teleo is None:
            return species_s_tree(st, non, T2) if non else None
        if not non:
            return teleo
        outer = species_s_tree(st, non + [("__T__", "DANRE")], T2)
        for nd in outer.preorder():
            for k, c in enumerate(nd.children):
                if c.name == "__T__":
                    nd.children[k] = teleo
        return teleo if outer.name == "__T__" else outer

    from hobench.newick import Node
    globals()["Node"] = Node

    # LWS: human OPN1LW/OPN1MW local duplicate
    lws_h = genes("LWS", lambda x: x["species"] == "HUMAN")
    lws_non = genes("LWS", lambda x: x["species"] not in TELEOSTS and x["species"] != "HUMAN")
    if len(lws_h) >= 2:
        i = f"{FID}.D{next(nid)}"
        hnode = mk("D", [Node(g) for g, _ in lws_h], id=i, tag="TAND", tier="T1", src=SRC)
        A(i, "event", "OPN1LW|OPN1MW primate duplicate (lvl=HUMAN; panel has one primate)", "primate_lws", "Background p.2")
        lws_non = lws_non + [("__H__", "HUMAN")]
    else:
        hnode = None
        lws_non = lws_non + lws_h
    tel = genes("LWS", lambda x: x["species"] in TELEOSTS)
    if len(tel) > 1:
        i = f"{FID}.U{next(nid)}"
        teleo = mk("U", [Node(g) for g, _ in tel], id=i, tier="T1", src=SRC)
        A(i, "event", "LWS: teleost local duplicates left unresolved (U)", "teleost_local", "Results p.1")
    else:
        teleo = Node(tel[0][0]) if tel else None
    lws = attach(lws_non, teleo)
    if hnode is not None:
        for nd in lws.preorder():
            for k, c in enumerate(nd.children):
                if c.name == "__H__":
                    nd.children[k] = hnode
    lws.name = "LWS"

    # SWS1: gar local duplicate
    gar = genes("SWS1", lambda x: x["species"] == "LEPOC")
    non = genes("SWS1", lambda x: x["species"] not in TELEOSTS | {"LEPOC"})
    tel = genes("SWS1", lambda x: x["species"] in TELEOSTS)
    if len(gar) >= 2:
        i = f"{FID}.D{next(nid)}"
        gnode = mk("D", [Node(g) for g, _ in gar], id=i, tag="TAND", tier="T1", src=SRC)
        A(i, "event", "gar SWS1-1|SWS1-2 local duplicate (lvl=LEPOC)", "gar_sws1", "Results p.1")
        # gar block + teleosts: S(gar D-node, teleosts)
        tb = species_s_tree(st, tel, T2) if tel else None
        act = mk("S", [gnode, tb], tier="T2", src=SRC)
        sws1 = attach(non, None)
        sws1 = mk("S", [species_s_tree(st, non, T2), act], tier="T2", src=SRC)
    else:
        sws1 = attach(non + gar, species_s_tree(st, tel, T2) if tel else None)
    sws1.name = "SWS1"
    sws2 = subtype_tree("SWS2")
    sws2.name = "SWS2"
    rh2 = subtype_tree("RH2")
    rh2.name = "RH2"

    # RH1
    def is_exo(x):
        return x["species"] in TELEOSTS | {"LEPOC"} and (re.search(r"RH1[ -]1$|exorh", x["seqname"]) is not None)

    def is_rho(x):
        return x["species"] in TELEOSTS | {"LEPOC"} and not is_exo(x)
    exo = genes("RH1", is_exo)
    rho_gar = genes("RH1", lambda x: is_rho(x) and x["species"] == "LEPOC")
    rho_tel = genes("RH1", lambda x: is_rho(x) and x["species"] in TELEOSTS)
    sarc = genes("RH1", lambda x: x["species"] not in TELEOSTS | {"LEPOC"})
    exo_node = species_s_tree(st, exo, T2) if exo else None
    if exo_node is not None and not exo_node.is_leaf:
        exo_node.name = "exorh"
    cnt = pd.Series([s for _, s in rho_tel]).value_counts()
    iu = f"{FID}.U{next(nid)}"
    rho_t = (mk("U", [Node(g) for g, _ in rho_tel], id=iu, tier="T1", src=SRC) if (cnt > 1).any()
             else species_s_tree(st, rho_tel, T2)) if rho_tel else None
    if (cnt > 1).any():
        A(iu, "event", "zebrafish rho|rhol (and medaka RH1-2) unresolved: 3R origin not certain -> U",
          "rhol_unsure", "Results p.1")
    rho_node = mk("S", [species_s_tree(st, rho_gar, T2) if rho_gar else None, rho_t], tier="T1", src=SRC)
    if rho_node is not None and not rho_node.is_leaf:
        rho_node.name = "rho"
    i = f"{FID}.D{next(nid)}"
    act = mk("D", [exo_node, rho_node], id=i, lvl="Neopterygii", tag="RETRO", tier="T1", src=SRC)
    A(i, "event", "exorh|rho split by retrotransposition before holostean-teleost divergence (lvl=Neopterygii)",
      "retro", "Results p.1-2")
    A(i, "event", "RH1-1/RH1-2 present in all investigated teleosts", "rh1_teleost", "Results p.1")
    rh1 = mk("S", [species_s_tree(st, sarc, T2) if sarc else None, act], tier="T2", src=SRC)
    rh1.name = "RH1"

    i = f"{FID}.D{next(nid)}"
    sws = mk("D", [sws1, sws2, rh1, rh2], id=i, lvl_lo="Gnathostomata", lvl_hi="Vertebrata", tag="WGD",
             tier="T1", src=SRC)
    sws.name = "SWS_2R"
    A(i, "event", "SWS -> SWS1, SWS2, RH1, RH2 in 2R (interval Gnathostomata..Vertebrata)", "2R", "Discussion")
    A(i, "event", "2R before lamprey divergence per Lagman (upper end of interval)", "lamprey", "Abstract")
    i = f"{FID}.D{next(nid)}"
    root = mk("D", [lws, sws], id=i, lvl_lo="Vertebrata", lvl_hi="ABOVE_RANGE", tag="TAND", tier="T1", src=SRC)
    A(i, "event", "SWS|LWS local duplication before 2R", "local_pre2R", "Discussion")

    # ---------------------------------------------------------------- outputs
    members = pd.DataFrame([dict(gene_id=x["gene_id"], species=x["species"], original_id_in_paper=x["original_id"],
                                 label=x["label"], mapping_method=x["mapping_method"],
                                 mapping_detail=x.get("mapping_detail", ""),
                                 mapping_identity=x["mapping_identity"], mapping_coverage=x["mapping_coverage"],
                                 source=SRC, locator=f"S1 sheet OPN row {x['row_xls']}") for x in mapped])
    members.to_csv(HERE / "members.tsv", sep="\t", index=False)
    rep = pd.DataFrame([dict(original_id=x["original_id"], species_label=x["species_label"], species=x["species"] or "",
                             label=f"{x['subtype']} {x['seqname']}", status=x["status"], gene_id=x["gene_id"],
                             reason=x["reason"], locator=f"S1 sheet OPN row {x['row_xls']}") for x in recs])
    rep.to_csv(HERE / "mapping_report.tsv", sep="\t", index=False)
    (HERE / "reference_tree.nhx").write_text(vlib.nhx(root) + "\n")
    pd.DataFrame(assertions).to_csv(HERE / "assertions.tsv", sep="\t", index=False)
    inpanel = rep[(rep.species != "") & (rep.status != "excluded")]
    n_in, n_map = len(inpanel), int((inpanel.status == "mapped").sum())
    per_sp = {k: f"{int((d.status == 'mapped').sum())}/{len(d)}" for k, d in inpanel.groupby("species")}
    meta = dict(
        id=FID, name="Vertebrate visual opsins (LWS, SWS1, SWS2, RH2, RH1)", panel="V",
        phenomena=["WGD", "TAND", "LOSS", "RETRO"],
        descriptors=dict(n_members=len(mapped), n_species=int(members.species.nunique()),
                         mapping=f"{n_map}/{n_in} in-panel S1 rows mapped ({n_map / max(n_in, 1):.0%})",
                         per_species_mapped=per_sp),
        sources=[dict(id=SRC, doi="10.1186/1471-2148-13-238", pmid="24180662",
                      title="The vertebrate ancestral repertoire of visual opsins, transducin alpha subunits and "
                            "oxytocin/vasopressin receptors was established by duplication of their shared genomic "
                            "region in the two rounds of early vertebrate genome duplications",
                      year=2013, source_method="synteny",
                      table="Additional file 1 (S1.xlsx), sheet OPN; Table 1", locator="S1 OPN; Results; Discussion",
                      quote=Q["2R"],
                      files={str(S1.relative_to(vlib.ROOT)): vlib.sha256(S1),
                             str((PAPER / 'fulltext.txt').relative_to(vlib.ROOT)): vlib.sha256(PAPER / 'fulltext.txt')}),
                 dict(id="Ramirez2016", doi="10.1093/gbe/evw248", pmid="28007952",
                      title="The last common ancestor of most bilaterian animals possessed at least 9 opsins",
                      year=2016, source_method="phylogeny", locator="Supplemental Table T1 (not used for members)",
                      quote="", files={str(RAMIREZ.relative_to(vlib.ROOT)): vlib.sha256(RAMIREZ)})],
        root_level="Vertebrata", root_single_copy=False, complete_species=[],
        decoys=[], distinct_from=[],
        tier_summary=dict(T1_events=sum(a["type"] == "event" and a["node_id"] != "" for a in assertions),
                          T2="subtype S polytomies (S1 subtype column)"),
        protocol_version=1,
        notes=("Scope = visual opsins only (Lagman S1 'OPN' main section). complete_species empty: Lagman searched "
               "Ensembl assemblies of 2013 (chicken lacked 3 of 5 genes in WASHUC2; frog RH2 possibly assembly gap). "
               "Gar members mapped by Ensembl locus overlap at the S1 gene-start coordinate (LepOcu1, still the "
               "Ensembl gar assembly); gar IDs printed in S1 are the query proteins of other species. Ramirez 2016 "
               "Table T1 lists database protein names only (no per-sequence classification) and was not used for "
               "membership; non-visual opsins therefore not covered. Root level Vertebrata with lvl_hi=ABOVE_RANGE: "
               "the SWS/LWS split predates 2R; Ciona/amphioxus visual-opsin relatives not classified by the source."),
    )
    (HERE / "family.yaml").write_text(yaml.safe_dump(meta, sort_keys=False, allow_unicode=True, width=120))
    print(f"{FID}: {n_map}/{n_in} mapped; {per_sp}")


if __name__ == "__main__":
    main()
