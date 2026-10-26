#!/usr/bin/env python
"""V03_Hox: vertebrate Hox genes (clusters A-D; paralog groups 1-13) from expert nomenclature with a synteny
basis, compiled into an S/D/U tree with events quoted from OA literature.

Members (nomenclature authorities; files stored under literature/papers/NOMEN_*/raw, sha256 in family.yaml):
  HUMAN  HGNC complete set (symbols HOXA1..HOXD13; gene_group 'HOXL subclass homeoboxes')
  MOUSE  MGI MRK_SwissProt_TrEMBL.rpt (Hoxa1..Hoxd13)
  DANRE  ZFIN uniprot.txt / gene.txt (hoxa1a..hoxd13a; a/b = teleost 3R cluster copies, Amores 1998 nomenclature)
  XENTR  Xenbase GeneExternalRef.txt (hoxa1..; taxon 8364)
  XENLA  Xenbase GeneExternalRef.txt (hoxa1.L / hoxa1.S; taxon 8355; L/S = subgenomes)
Human HOX gene-family assignment (Hox1..Hox5, Hox6-8, Hox9-13) cross-checked against Holland et al. 2007
Table 2 (PMID 17963489).

Tree:
  root            D  Holland gene families split before Bilateria (lvl=Chordata, lvl_hi=ABOVE_RANGE)
  Hox6-8, Hox9-13 D  paralog groups within the family, pre-2R tandem duplications (lvl_lo=Vertebrata,
                     lvl_hi=ABOVE_RANGE, TAND, T2)
  PGn             D  2R -> clusters A-D (lvl_lo=Gnathostomata, lvl_hi=Vertebrata, WGD)
  cluster gene    S over species; DANRE a/b copies under D(Teleostei, WGD); XENLA .L/.S under D(XENLA, ALLO)
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
from hobench.newick import Node  # noqa: E402

FID = "V03_Hox"
R = vlib.ROOT
P = R / "literature/papers"
FILES = {
    "HGNC": P / "NOMEN_HGNC/raw/hgnc_complete_set.txt",
    "MGI": P / "NOMEN_MGI/raw/MRK_SwissProt_TrEMBL.rpt",
    "ZFIN_uniprot": P / "NOMEN_ZFIN/raw/uniprot.txt",
    "ZFIN_gene": P / "NOMEN_ZFIN/raw/gene.txt",
    "Xenbase": P / "NOMEN_Xenbase/raw/GeneExternalRef.txt",
}
TXT = {
    "Holland2007": (P / "PMID17963489/fulltext.txt").read_text(),
    "Kumar2024": (P / "PMID39079985/fulltext.txt").read_text(),
    "PascualAnaya2012": (P / "PMID23198682/fulltext.txt").read_text(),
    "Session2016": (P / "PMID27762356/fulltext.txt").read_text(),
    "Kikuchi2025": (P / "PMID41268905/fulltext.txt").read_text(),
}
Q = {
    "families": ("Holland2007", "The human genome has 39 HOX genes, arranged into four Hox clusters. Here we divide the "
                                "HOX genes into seven gene families: Hox1, Hox2, Hox3, Hox4, Hox5, Hox6-8 and Hox9-13."),
    "family_def": ("Holland2007", "We define a gene family as a set of genes derived from a single gene in the most "
                                  "recent common ancestor of bilaterian animals"),
    "2R": ("Kumar2024", "In early vertebrates, two rounds of genome duplication events (2R) produced the four "
                        "paralogous Hox clusters (HoxA, HoxB, HoxC, and HoxD)"),
    "3R": ("Kumar2024", "Teleosts, however, have seven or eight Hox clusters due to the TS-3R WGD and subsequent gene "
                        "loss in some species"),
    "posterior": ("PascualAnaya2012", "Posterior Hox genes of cephalochordates, urochordates, echinoderms, and "
                                      "hemichordates do not have clear orthologous relationships to the posterior "
                                      "paralogy groups of vertebrates"),
    "amphioxus15": ("PascualAnaya2012", "possesses a single Hox cluster of 15 genes"),
    "human39": ("Holland2007", "The human genome has 39 HOX genes, arranged into four Hox clusters."),
    "zf49": ("Kikuchi2025", "Schematic representation of the 49 hox genes organized into seven hox clusters in zebrafish."),
    "LS": ("Session2016", "The L and S chromosome sets therefore represent the descendants of two distinct diploid "
                          "progenitors, confirming the allotetraploid hypothesis even in the absence of extant "
                          "progenitor species."),
    "LS_time": ("Session2016", "the L and S subgenomes diverged from each other ~34 Mya (T_2) and from X. tropicalis "
                               "~48 Mya (T_1)"),
}
for k, (src, q) in Q.items():
    assert vlib.quote_in(TXT[src], q), f"quote {k} not found in {src}"

FAMILY_OF_PG = {1: "Hox1", 2: "Hox2", 3: "Hox3", 4: "Hox4", 5: "Hox5", 6: "Hox6-8", 7: "Hox6-8", 8: "Hox6-8",
                9: "Hox9-13", 10: "Hox9-13", 11: "Hox9-13", 12: "Hox9-13", 13: "Hox9-13"}


def load_members():
    recs = []
    hg = pd.read_csv(FILES["HGNC"], sep="\t", dtype=str).fillna("")
    for i, r in hg.iterrows():
        if re.match(r"^HOX[ABCD]\d+$", r["symbol"]) and r["locus_type"] == "gene with protein product":
            recs.append(dict(species="HUMAN", symbol=r["symbol"], authority="HGNC", auth_id=r["hgnc_id"],
                             cands=[a for a in r["uniprot_ids"].split("|") if a], geneid=r["entrez_id"],
                             locator=f"hgnc_complete_set.txt row {i + 2}", quote=f"{r['hgnc_id']}\t{r['symbol']}\t{r['name']}"))
    with open(FILES["MGI"]) as fh:
        for i, line in enumerate(fh):
            f = line.rstrip("\n").split("\t")
            if len(f) > 6 and re.match(r"^Hox[abcd]\d+$", f[1]) and f[2] == "O":
                recs.append(dict(species="MOUSE", symbol=f[1], authority="MGI", auth_id=f[0], cands=f[6].split(),
                                 geneid="", locator=f"MRK_SwissProt_TrEMBL.rpt line {i + 1}",
                                 quote="\t".join(f[:4])))
    zu: dict[str, list] = {}
    zline: dict[str, int] = {}
    with open(FILES["ZFIN_uniprot"]) as fh:
        for i, line in enumerate(fh):
            f = line.rstrip("\n").split("\t")
            if len(f) >= 4 and re.match(r"^hox[abcd]\d+[ab]$", f[2]):
                zu.setdefault((f[0], f[2]), []).append(f[3])
                zline.setdefault((f[0], f[2]), i + 1)
    zg = {}
    with open(FILES["ZFIN_gene"]) as fh:
        for i, line in enumerate(fh):
            f = line.rstrip("\n").split("\t")
            if len(f) >= 4 and re.match(r"^hox[abcd]\d+[ab]$", f[2]):
                zg[(f[0], f[2])] = (f[3], i + 1)
    for key in sorted(set(zu) | set(zg)):
        zid, sym = key
        recs.append(dict(species="DANRE", symbol=sym, authority="ZFIN", auth_id=zid, cands=zu.get(key, []),
                         geneid=zg.get(key, ("", 0))[0],
                         locator=f"ZFIN gene.txt line {zg.get(key, ('', 0))[1]}; uniprot.txt line {zline.get(key, 0)}",
                         quote=f"{zid}\tSO:0001217\t{sym}"))
    with open(FILES["Xenbase"]) as fh:
        for i, line in enumerate(fh):
            f = line.rstrip("\n").split("\t")
            if len(f) < 5:
                continue
            if f[3] == "8364" and re.match(r"^hox[abcd]\d+$", f[1]):
                sp = "XENTR"
            elif f[3] == "8355" and re.match(r"^hox[abcd]\d+\.[LS]$", f[1]):
                sp = "XENLA"
            else:
                continue
            recs.append(dict(species=sp, symbol=f[1], authority="Xenbase", auth_id=f[0],
                             cands=[a for a in (f[6] if len(f) > 6 else "").split("|") if re.match(r"^[A-Z0-9]{6,10}$", a)],
                             geneid=f[4], locator=f"Xenbase GeneExternalRef.txt line {i + 1}",
                             quote="\t".join(f[:4])))
    for x in recs:
        m = re.match(r"^hox([abcd])(\d+)([ab]?)(?:\.([LS]))?$", x["symbol"], re.I)
        x["cluster"], x["pg"], x["tcopy"], x["sub"] = m.group(1).upper(), int(m.group(2)), m.group(3).lower(), m.group(4) or ""
    return recs


def map_members(recs, cache):
    for x in recs:
        x.update(status="unmapped", gene_id="", mapping_method="", mapping_identity="", mapping_coverage="", reason="")
        prot = vlib.proteome(x["species"])
        inp = sorted({a for a in x["cands"] if a in prot})
        if len(inp) == 1:
            x.update(status="mapped", gene_id=inp[0], mapping_method="exact_id",
                     mapping_detail=f"{x['authority']} UniProt cross-reference",
                     mapping_identity="100", mapping_coverage="1.0", reason=f"{len(x['cands'])} xref(s), 1 in proteome")
            continue
        if len(inp) > 1:
            gh = sorted({a for a in (cache.get("geneid|" + x["geneid"], lambda: _geneid_map(x["geneid"]))
                                     if x["geneid"] else []) if a in prot})
            both = sorted(set(gh) & set(inp))
            if len(both) == 1:
                x.update(status="mapped", gene_id=both[0], mapping_method="exact_id",
                         mapping_detail=f"{x['authority']} xrefs {inp} disambiguated by NCBI GeneID {x['geneid']}",
                         mapping_identity="100", mapping_coverage="1.0", reason="GeneID disambiguation")
            else:
                x.update(status="ambiguous", reason=f"{len(inp)} authority UniProt xrefs in panel proteome: {inp}; "
                                                    f"GeneID -> {gh}")
            continue
        if x["geneid"]:
            rows = cache.get("geneid|" + x["geneid"], lambda: _geneid_map(x["geneid"]))
            hit = sorted({a for a in rows if a in prot})
            if len(hit) == 1:
                x.update(status="mapped", gene_id=hit[0], mapping_method="exact_id",
                         mapping_detail=f"{x['authority']} -> NCBI GeneID {x['geneid']} -> UniProt ID mapping",
                         mapping_identity="100", mapping_coverage="1.0",
                         reason="authority UniProt xrefs not in proteome; NCBI GeneID mapping")
                continue
            x["reason"] = f"GeneID {x['geneid']} -> {len(hit)} panel accessions; xrefs {x['cands'][:4]} not in proteome"
            if x["cands"]:
                x["candidates"] = x["cands"]
        else:
            x["reason"] = f"xrefs {x['cands'][:4]} not in proteome; no GeneID"
            x["candidates"] = x["cands"]
    for sp in sorted({x["species"] for x in recs}):
        vlib.resolve_by_sequence(cache, sp, [x for x in recs if x["species"] == sp], strict_family=True)


def _geneid_map(gid):
    from curation_tools.map_ids import run_job
    try:
        rows, _ = run_job("GeneID", "UniProtKB", [gid])
    except Exception:
        return []
    return [r.get("Entry", "") for r in rows]


def main():
    st = vlib.species_tree()
    cache = vlib.JsonCache(HERE / "cache" / "lookups.json")
    recs = load_members()
    map_members(recs, cache)
    # identical accession for 2 symbols -> ambiguous
    by: dict[str, list] = {}
    for x in recs:
        if x["status"] == "mapped":
            by.setdefault(x["gene_id"], []).append(x)
    for acc, xs in by.items():
        if len(xs) > 1:
            for x in xs:
                x.update(status="ambiguous", gene_id="", reason=f"{acc} shared by {[y['symbol'] for y in xs]}")
    mapped = [x for x in recs if x["status"] == "mapped"]
    # cross-check Holland families vs PG rule for human genes
    hol_ff, cur = {}, None
    for line in TXT["Holland2007"].splitlines():
        f = line.split("\t")
        if len(f) >= 2 and re.match(r"^HOX[ABCD]\d+$", f[1]):
            cur = f[0] or cur
            hol_ff[f[1]] = cur
    mism = [s for s, fam in hol_ff.items() if FAMILY_OF_PG[int(re.sub(r"\D", "", s[4:]))] != fam]
    assert not mism, f"PG->Holland family rule disagrees with Table 2 for {mism}"

    assertions = []
    nid = iter(range(1, 10000))
    srcmap = {"HUMAN": "HGNC", "MOUSE": "MGI", "DANRE": "ZFIN", "XENTR": "Xenbase", "XENLA": "Xenbase"}

    def A(node_id, typ, target, qkey, loc, tier="T1"):
        src, q = Q[qkey]
        assertions.append(dict(assertion_id=f"{FID}.A{len(assertions) + 1:03d}", type=typ, node_id=node_id,
                               target=node_id if typ == "event" else target, note=target, tier=tier, source_id=src, locator=loc, quote=q))

    def gene_node(cl, pg):
        g = [x for x in mapped if x["cluster"] == cl and x["pg"] == pg]
        if not g:
            return None
        parts = []
        # zebrafish a/b
        za = [(x["gene_id"], x["species"]) for x in g if x["species"] == "DANRE" and x["tcopy"] == "a"]
        zb = [(x["gene_id"], x["species"]) for x in g if x["species"] == "DANRE" and x["tcopy"] == "b"]
        z = None
        if za and zb:
            i = f"{FID}.D{next(nid)}"
            z = mk("D", [Node(za[0][0]), Node(zb[0][0])], id=i, lvl="Teleostei", tag="WGD", tier="T2", src="ZFIN")
            A(i, "event", f"hox{cl.lower()}{pg}a|hox{cl.lower()}{pg}b = teleost 3R cluster copies", "3R",
              "Kumar2024 'gene clusters' paragraph + ZFIN symbols", tier="T2")
        elif za or zb:
            z = Node((za or zb)[0][0])
        xl = [x for x in g if x["species"] == "XENLA"]
        xln = None
        if len(xl) == 2 and {x["sub"] for x in xl} == {"L", "S"}:
            i = f"{FID}.D{next(nid)}"
            xln = mk("D", [Node(x["gene_id"]) for x in xl], id=i, tag="ALLO", tier="T2", src="Xenbase")
            A(i, "event", f"hox{cl.lower()}{pg}.L|.S homeologs (X. laevis allotetraploidy; lvl=XENLA)", "LS",
              "Session2016 results; Xenbase .L/.S symbols", tier="T2")
        elif len(xl) == 1:
            xln = Node(xl[0]["gene_id"])
        others = [(x["gene_id"], x["species"]) for x in g if x["species"] not in ("DANRE", "XENLA")]
        ph = list(others)
        if z is not None:
            ph.append(("__Z__", "DANRE"))
        if xln is not None:
            ph.append(("__X__", "XENLA"))
        tree = species_s_tree(st, ph, {"tier": "T2", "src": "nomenclature"})
        for nd in [tree] + list(tree.preorder()):
            for k, c in enumerate(nd.children):
                if c.name == "__Z__":
                    nd.children[k] = z
                elif c.name == "__X__":
                    nd.children[k] = xln
        if tree.name == "__Z__":
            tree = z
        elif tree.name == "__X__":
            tree = xln
        if not tree.is_leaf:
            tree.name = f"HOX{cl}{pg}"
        return tree

    fam_nodes = {}
    for pg in range(1, 15):
        cl_nodes = [gene_node(cl, pg) for cl in "ABCD"]
        cl_nodes = [c for c in cl_nodes if c is not None]
        if not cl_nodes:
            continue
        if len(cl_nodes) > 1:
            i = f"{FID}.D{next(nid)}"
            pgn = mk("D", cl_nodes, id=i, lvl_lo="Gnathostomata", lvl_hi="Vertebrata", tag="WGD", tier="T1",
                     src="Kumar2024")
            A(i, "event", f"PG{pg}: clusters A-D from 2R (interval Gnathostomata..Vertebrata)", "2R",
              "Kumar2024 'gene clusters' paragraph")
        else:
            pgn = cl_nodes[0]
        if not pgn.is_leaf:
            pgn.name = pgn.name or f"PG{pg}"
            if pgn.features.get("ev") == "D":
                pgn.name = f"PG{pg}"
        fam_nodes.setdefault(FAMILY_OF_PG[pg], []).append(pgn)
    top = []
    for fam in ["Hox1", "Hox2", "Hox3", "Hox4", "Hox5", "Hox6-8", "Hox9-13"]:
        kids = fam_nodes.get(fam, [])
        if len(kids) > 1:
            i = f"{FID}.D{next(nid)}"
            n = mk("D", kids, id=i, lvl_lo="Vertebrata", lvl_hi="ABOVE_RANGE", tag="TAND", tier="T2", src="Holland2007")
            n.name = fam.replace("-", "_")
            A(i, "event", f"{fam}: paralog groups predate 2R (clusters arose by 2R duplication of a cluster that "
                          f"already contained them)", "2R", "Kumar2024 + Holland 2007 Table 2", tier="T2")
            A(i, "event", f"{fam}: posterior PG orthology to invertebrate genes unclear -> no lower bound beyond Vertebrata",
              "posterior", "Pascual-Anaya 2012 Discussion", tier="T2")
        elif kids:
            n = kids[0]
        else:
            continue
        top.append(n)
    i = f"{FID}.D{next(nid)}"
    root = mk("D", top, id=i, lvl_lo="Chordata", lvl_hi="ABOVE_RANGE", tag="TAND", tier="T1", src="Holland2007")
    A(i, "event", "Holland gene families Hox1..Hox9-13 = distinct genes in the bilaterian ancestor", "families",
      "Holland 2007 Results/Table 2")
    A(i, "event", "gene family definition (single gene in bilaterian ancestor)", "family_def", "Holland 2007 Methods")

    # ---------------------------------------------------------------- outputs
    members = pd.DataFrame([dict(gene_id=x["gene_id"], species=x["species"], original_id_in_paper=f"{x['auth_id']}",
                                 label=x["symbol"], mapping_method=x["mapping_method"],
                                 mapping_detail=x.get("mapping_detail", ""),
                                 mapping_identity=x["mapping_identity"], mapping_coverage=x["mapping_coverage"],
                                 source=srcmap[x["species"]], locator=x["locator"], quote=x["quote"]) for x in mapped])
    members.to_csv(HERE / "members.tsv", sep="\t", index=False)
    rep = pd.DataFrame([dict(original_id=x["auth_id"], species=x["species"], label=x["symbol"], status=x["status"],
                             gene_id=x["gene_id"], reason=x["reason"], locator=x["locator"]) for x in recs])
    rep.to_csv(HERE / "mapping_report.tsv", sep="\t", index=False)
    (HERE / "reference_tree.nhx").write_text(vlib.nhx(root) + "\n")
    pd.DataFrame(assertions).to_csv(HERE / "assertions.tsv", sep="\t", index=False)
    n_in, n_map = len(rep), int((rep.status == "mapped").sum())
    per_sp = {k: f"{int((d.status == 'mapped').sum())}/{len(d)}" for k, d in rep.groupby("species")}
    complete = [k for k, v in per_sp.items() if int(v.split("/")[0]) / int(v.split("/")[1]) >= 0.9
                and k in ("HUMAN", "DANRE")]  # only species with a quoted complete count
    A("", "completeness", "HUMAN", "human39", "Holland 2007 Results (ANTP class)")
    A("", "completeness", "DANRE", "zf49", "Kikuchi 2025 Figure 5 legend")
    pd.DataFrame(assertions).to_csv(HERE / "assertions.tsv", sep="\t", index=False)
    paper = lambda pid: P / f"PMID{pid}/fulltext.txt"
    meta = dict(
        id=FID, name="Hox genes (clusters A-D, paralog groups 1-13)", panel="V",
        phenomena=["WGD", "TAND", "LOSS", "ALLO", "DEEP"],
        descriptors=dict(n_members=len(mapped), n_species=int(members.species.nunique()),
                         mapping=f"{n_map}/{n_in} nomenclature entries mapped ({n_map / n_in:.0%})",
                         per_species_mapped=per_sp),
        sources=[
            dict(id="HGNC", doi="10.1093/nar/gkac888", pmid="36243972", title="HGNC complete set (downloaded 2026-10-26)",
                 year=2026, source_method="nomenclature", table="hgnc_complete_set.txt",
                 locator="symbol HOXA1..HOXD13", quote="HOXA1\thomeobox A1",
                 files={str(FILES['HGNC'].relative_to(R)): vlib.sha256(FILES['HGNC'])}),
            dict(id="MGI", doi="10.1093/genetics/iyad031", pmid="36869748", title="MGI MRK_SwissProt_TrEMBL.rpt (2026-10-26)",
                 year=2026, source_method="nomenclature", table="MRK_SwissProt_TrEMBL.rpt", locator="Hoxa1..Hoxd13",
                 quote="Hoxa1\tO\thomeobox A1", files={str(FILES['MGI'].relative_to(R)): vlib.sha256(FILES['MGI'])}),
            dict(id="ZFIN", doi="10.1093/genetics/iyad032", pmid="36864651",
                 title="ZFIN gene.txt / uniprot.txt (2026-10-26); hox a/b cluster nomenclature (Amores et al. 1998)",
                 year=2026, source_method="nomenclature", table="gene.txt, uniprot.txt",
                 locator="hoxa1a..hoxd13a", quote="hoxa9b",
                 files={str(FILES[k].relative_to(R)): vlib.sha256(FILES[k]) for k in ("ZFIN_gene", "ZFIN_uniprot")}),
            dict(id="Xenbase", doi="10.1093/genetics/iyad018", pmid="36755307",
                 title="Xenbase GeneExternalRef.txt (2026-10-26); .L/.S subgenome nomenclature", year=2026,
                 source_method="nomenclature", table="GeneExternalRef.txt",
                 locator="hoxa1, hoxa1.L, hoxa1.S ...", quote="hoxa1.L\thomeobox A1",
                 files={str(FILES['Xenbase'].relative_to(R)): vlib.sha256(FILES['Xenbase'])}),
            dict(id="Holland2007", doi="10.1186/1741-7007-5-47", pmid="17963489",
                 title="Classification and nomenclature of all human homeobox genes", year=2007,
                 source_method="expert_classification", table="Table 2", locator="Results, ANTP class",
                 quote=Q["families"][1], files={str(paper(17963489).relative_to(R)): vlib.sha256(paper(17963489))}),
            dict(id="Kumar2024", doi="10.1007/s00427-024-00720-6", pmid="39079985",
                 title="Zebrafish: unraveling genetic complexity through duplicated genes", year=2024,
                 source_method="synteny", locator="section on gene clusters", quote=Q["2R"][1],
                 files={str(paper(39079985).relative_to(R)): vlib.sha256(paper(39079985))}),
            dict(id="PascualAnaya2012", doi="10.1186/2041-9139-3-28", pmid="23198682",
                 title="Broken colinearity of the amphioxus Hox cluster", year=2012, source_method="phylogeny",
                 locator="Discussion", quote=Q["posterior"][1],
                 files={str(paper(23198682).relative_to(R)): vlib.sha256(paper(23198682))}),
            dict(id="Session2016", doi="10.1038/nature19840", pmid="27762356",
                 title="Genome evolution in the allotetraploid frog Xenopus laevis", year=2016,
                 source_method="synteny", locator="Results 'Subgenome identity'", quote=Q["LS"][1],
                 files={str(paper(27762356).relative_to(R)): vlib.sha256(paper(27762356))}),
            dict(id="nomenclature", doi="10.1186/1741-7007-5-47", pmid="17963489",
                 title="(alias) S polytomies compiled from the nomenclature sources above", year=2026,
                 source_method="nomenclature", locator="members.tsv quote column", quote=Q["families"][1],
                 files={str(paper(17963489).relative_to(R)): vlib.sha256(paper(17963489))}),
            dict(id="Kikuchi2025", doi="10.7554/elife.105889", pmid="41268905",
                 title="HoxB-derived hoxba and hoxbb clusters are essential for the anterior-posterior positioning of "
                       "zebrafish pectoral fins", year=2025, source_method="nomenclature", locator="Figure 5 legend",
                 quote=Q["zf49"][1], files={str(paper(41268905).relative_to(R)): vlib.sha256(paper(41268905))}),
        ],
        root_level="Chordata", root_single_copy=False, complete_species=complete, decoys=[], distinct_from=[],
        tier_summary=dict(T1_events=sum(a["tier"] == "T1" and a["node_id"] != "" for a in assertions),
                          T2_events=sum(a["tier"] == "T2" and a["node_id"] != "" for a in assertions)),
        protocol_version=1,
        notes=("Scope v1: species with an expert nomenclature authority (HGNC, MGI, ZFIN, Xenbase); gar, medaka, "
               "salmon, coelacanth, shark, lamprey and amphioxus Hox lack machine-readable OA expert tables "
               "(Braasch 2016 gar-teleost tables use Ensembl Compara GeneTrees -> denylisted; Amores 1998, Kuraku & "
               "Meyer 2009, Mungpakdee 2008 not OA/not retrieved). The 2R interval and the pre-2R paralog-group splits "
               "lie above the members' LCA (Euteleostomi), so cross-cluster and cross-PG pairs have merge level "
               "'never within range'. complete_species = species with >=90% of authority Hox symbols mapped "
               "(authority lists are genome-wide nomenclature sets)."),
    )
    (HERE / "family.yaml").write_text(yaml.safe_dump(meta, sort_keys=False, allow_unicode=True, width=120))
    print(f"{FID}: {n_map}/{n_in} mapped; {per_sp}")


if __name__ == "__main__":
    main()
