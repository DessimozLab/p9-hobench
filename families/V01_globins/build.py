#!/usr/bin/env python
"""V01_globins: chordate globins (vertebrate Ngb, GbX, Cygb, Mb, GbE, GbY, alpha/beta-Hb; lamprey Hbs; Ciona and
amphioxus globins).

Members: Hoffmann et al. 2012b MBE 29:1735 (PMID 22319164) Supplementary Table S1 (labels + accessions; received
as requested item #04 and stored under literature/papers/PMID22319164/raw/received/). Events (all quoted):
  * Hoffmann et al. 2012b text: four deuterostome globin clades (clade 4 = Ngb + amphioxus Gb4 [+ POMT2 synteny];
    clade 2 = GbX; clade 3 = amphioxus Gb3, 6, 12-14; clade 1 = vertebrate-specific + urochordate + other amphioxus
    globins); clade 4 sister to the rest, clade 2 sister to clade 3; each clade contains ambulacrarian globins
    -> splits above the panel root.
  * Hoffmann et al. 2012a MBE 29:303 (PMID 21965344): 1R -> (Mb/GbE/Cygb) vs (alpha/beta-Hb + GbY); 2R -> Cygb vs
    Mb/GbE; tandem duplication -> proto alpha/beta before the gnathostome ancestor.
  * Hoffmann et al. 2010 PNAS (PMID 20660759): cyclostome Hbs are not orthologous to gnathostome Hbs (cooption of a
    Cygb-lineage gene); Mb/GbE duplication predates the gnathostome radiation; four lineages present in the
    vertebrate ancestor.
Contested 2R timing vs the cyclostome split (Simakov 2020 / hagfish genomes 2024 vs Hoffmann's 2R-before-split):
2R = interval Gnathostomata..Vertebrata and lamprey Hbs sit under a U with the (Cygb, Mb/GbE) subtree, so only the
uncontested claim (lamprey Hb vs gnathostome alpha/beta-Hb/GbY = paralogs, split at 1R) is scored.
Teleost members of one label are WGD (3R)-descended and not synteny-curated -> U blocks (protocol 3b).
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

FID = "V01_globins"
R = vlib.ROOT
P12B = R / "literature/papers/PMID22319164"
P12A = R / "literature/papers/PMID21965344"
P10 = R / "literature/papers/PMID20660759"
SUPP_PDF = P12B / "raw/received/supp/Deuterostome_globins_Supp_Tables.pdf"
SUPP_TXT = P12B / "raw/received/supp/Deuterostome_globins_Supp_Tables.txt"
ALN = P12B / "raw/received/supp/Sequence_alignment.txt"
TXT = {"Hoffmann2012b": (P12B / "raw/fulltext.txt").read_text(),
       "Hoffmann2012a": (P12A / "raw/fulltext.txt").read_text(),
       "Hoffmann2010": (P10 / "raw/fulltext.txt").read_text()}

Q = {
    "clades": ("Hoffmann2012b", "Relationships among these four clades are well resolved, with the Ngb-containing "
                                "clade 4 sister to the remaining three clades."),
    "clade_amb": ("Hoffmann2012b", "into four well-supported clades (fig. 2), each of them defined by the presence of "
                                   "both chordate and ambulacrarian globins."),
    "clade23": ("Hoffmann2012b", "The second group, clade 2, contains acorn worm Gbs 7–10, 16, acoel Gb2, and "
                                 "vertebrate GbX and is sister to the third group, clade 3, which contains acorn worm "
                                 "Gb6 and Gb15, and amphioxus Gbs 3, 6, 12–14."),
    "clade4": ("Hoffmann2012b", "Finally, the fourth group, clade 4, includes the sea urchin subdomains, plus "
                                "vertebrate Ngb and amphioxus Gb4."),
    "pomt2": ("Hoffmann2012b", "vertebrate Ngb, amphioxus Gb4, and the multidomain globin gene of the purple sea "
                               "urchin (pGb) are each found adjacent to an ortholog of the human gene POMT2"),
    "clade1": ("Hoffmann2012b", "The first of these groups, clade 1, contains all vertebrate-specific globins (i. e., "
                                "Cygb [including cyclostome Hbs], Mb, α- and β-Hb, GbE, and GbY), plus two separate "
                                "groups of amphioxus globins, all urochordate globins"),
    "gb8": ("Hoffmann2012b", "Bayesian phylogenies arranged all deuterostome globins other than amphioxus Gb8 and the "
                             "acorn worm intron-less globins"),
    "1R": ("Hoffmann2012a", "The first round of WGD produced a proto “Mb/GbE/Cygb” paralogon and a proto “α/β-Hb + "
                            "GbY” paralogon"),
    "2R": ("Hoffmann2012a", "Reduplication of the former paralogon in the subsequent round of WGD gave rise to the "
                            "progenitors of the Cygb and Mb/GbE gene lineages"),
    "ab": ("Hoffmann2012a", "The subsequent tandem duplication that produced the proto α- and β-globin genes then "
                            "promoted a further specialization of function with respect to blood–gas transport."),
    "ab_gn": ("Hoffmann2012a", "In the common ancestor of gnathostomes, functional divergence of the proto α- and "
                               "β-globin genes permitted the formation of multimeric Hbs"),
    "four": ("Hoffmann2010", "The common ancestor of extant vertebrates possessed a globin gene repertoire that "
                             "included progenitors of each of these four distinct gene lineages."),
    "not_orth": ("Hoffmann2010", "Results of this analysis revealed that the erythroid Hbs of cyclostomes and "
                                 "gnathostomes are not orthologous proteins."),
    "cygb_orth": ("Hoffmann2010", "Instead, the cyclostome Hbs are orthologous to the hexacoordinate Cygbs of "
                                  "gnathostome vertebrates."),
    "mbgbe": ("Hoffmann2010", "It thus appears that Mb and GbE derive from a duplication event that predated the "
                              "gnathostome radiation"),
    "lamprey_lse": ("Hoffmann2010", "the globin gene repertoires of lampreys and hagfish were then further expanded "
                                    "by subsequent rounds of lineage-specific gene duplication and divergence."),
    "redundant": ("Hoffmann2012b", "to identify and remove redundant sequences"),
}
for k, (s, q) in Q.items():
    assert vlib.quote_in(TXT[s], q), f"quote {k} not found in {s}"

SPECIES = {"Homo sapiens": "HUMAN", "Gallus gallus": "CHICK", "Xenopus tropicalis": "XENTR", "Danio rerio": "DANRE",
           "Oryzias latipes": "ORYLA", "Petromyzon marinus": "PETMA", "Ciona intestinalis": "CIOIN",
           "Branchiostoma floridae": "BRAFL"}
TELEOSTS = {"DANRE", "ORYLA"}
ROW = re.compile(r"^(?P<label>.+?)\s(?P<plus>\+\s)?(?P<sp>" + "|".join(SPECIES) + r")\s+\S+\s+\S+\s+(?P<acc>\S+?)"
                 r"(?P<star>\*)?\s+(?P<src>GenBank|Ensembl)\s*$")


def gene_type(label: str) -> str:
    t = re.split(r"[\s_]", label.strip())[-1]
    t = {"HbA": "Hba", "HbB": "Hbb", "GlbX": "GbX"}.get(t, t)
    if re.match(r"^Cygb\d$", t):
        return "Cygb"
    return t


def parse_s1():
    recs = []
    for i, line in enumerate(SUPP_TXT.read_text().splitlines()):
        m = ROW.match(line.strip())
        if not m:
            continue
        recs.append(dict(line=i + 1, label=m.group("label").strip(), plus=bool(m.group("plus")),
                         species=SPECIES[m.group("sp")], acc=m.group("acc"), star=bool(m.group("star")),
                         src=m.group("src"), raw=line.strip()))
    return recs


def aln_seqs():
    out, name = {}, None
    for line in ALN.read_text().splitlines():
        if line.startswith(">"):
            name = line[1:].strip()
            out[name] = ""
        elif name:
            out[name] += re.sub(r"[^A-Za-z]", "", line)
    return out


def main():
    st = vlib.species_tree()
    cache = vlib.JsonCache(HERE / "cache" / "lookups.json")
    recs = parse_s1()
    aln = aln_seqs()
    for x in recs:
        x["gtype"] = gene_type(x["label"])
        x.update(status="", gene_id="", reason="", mapping_method="", mapping_identity="", mapping_coverage="")
        if x["species"] == "BRAFL" and not x["plus"]:
            x.update(status="excluded", reason="amphioxus gene model not in the authors' pruned set (removed as "
                                                "redundant, Methods) - the '+' Gb1-15 set represents the repertoire")
            continue
        oid = x["acc"]
        t = vlib.id_type(oid)
        if t:
            x.update(vlib.map_one(cache, x["species"], oid))
        elif re.match(r"^[A-Z]{1,2}\d{5,6}$", oid):
            rows = cache.get("embl|" + oid, lambda: _embl(oid))
            acc, note = vlib.pick_in_proteome(x["species"], [{"uniprot_acc": a} for a in rows])
            if acc:
                x.update(status="mapped", gene_id=acc, mapping_method="exact_id",
                         mapping_detail=f"EMBL/GenBank nucleotide {oid} -> UniProt ID mapping",
                         mapping_identity="100", mapping_coverage="1.0", reason=note)
            else:
                x.update(status="unmapped", reason=f"EMBL {oid}: {note}", candidates=rows)
        else:
            x.update(status="unmapped", reason=f"accession {oid!r} is not a resolvable sequence identifier")
        # paper's own aligned sequence as an additional same-species candidate (ungapped)
        key = x["label"].replace(" ", "_").replace("Sea_squirt_(C._intestinalis)", "Ciona_intestinalis")
        if x["status"] != "mapped" and key in aln and aln[key]:
            cache.data["seq|aln:" + key] = aln[key]
            x["candidates"] = list(x.get("candidates", [])) + ["aln:" + key]
    for sp in sorted({x["species"] for x in recs}):
        _resolve(cache, sp, [x for x in recs if x["species"] == sp])
    by: dict[str, list] = {}
    for x in recs:
        if x["status"] == "mapped":
            by.setdefault(x["gene_id"], []).append(x)
    for acc, xs in by.items():
        if len(xs) > 1:
            for x in xs:
                x.update(status="ambiguous", gene_id="", reason=f"{acc} assigned to {[y['label'] for y in xs]}")
    mapped = [x for x in recs if x["status"] == "mapped"]

    # ---------------------------------------------------------------- tree
    assertions, nid = [], iter(range(1, 1000))

    def A(node_id, typ, note, qkey, loc, tier="T1", target=None):
        s, q = Q[qkey]
        assertions.append(dict(assertion_id=f"{FID}.A{len(assertions) + 1:03d}", type=typ,
                               target=node_id if typ == "event" else (target or ""), note=note, tier=tier,
                               source_id=s, locator=loc, quote=q))

    def G(pred):
        return [(x["gene_id"], x["species"]) for x in mapped if pred(x)]

    T2 = {"tier": "T2", "src": "Hoffmann2012b"}

    def label_tree(gt):
        non = G(lambda x: x["gtype"] == gt and x["species"] not in TELEOSTS)
        tel = G(lambda x: x["gtype"] == gt and x["species"] in TELEOSTS)
        if len(tel) > 1:
            i = f"{FID}.U{next(nid)}"
            tb = mk("U", [Node(g) for g, _ in tel], id=i, **T2)  # 3R-descended, not synteny-curated (protocol 3b)
        else:
            tb = Node(tel[0][0]) if tel else None
        if tb is None:
            n = species_s_tree(st, non, T2) if non else None
        elif not non:
            n = tb
        else:
            n = species_s_tree(st, non + [("__T__", "DANRE")], T2)
            for nd in [n] + list(n.preorder()):
                for k, c in enumerate(nd.children):
                    if c.name == "__T__":
                        nd.children[k] = tb
        if n is not None and not n.is_leaf:
            n.name = gt
        return n

    def Dn(kids, qkeys, note, loc, **lv):
        kids = [k for k in kids if k is not None]
        if len(kids) < 2:
            return kids[0] if kids else None
        i = f"{FID}.D{next(nid)}"
        n = mk("D", kids, id=i, tier="T1", src=Q[qkeys[0]][0], **lv)
        for q in qkeys:
            A(i, "event", note, q, loc)
        return n

    def Un(kids, name=""):
        kids = [k for k in kids if k is not None]
        if len(kids) < 2:
            return kids[0] if kids else None
        n = mk("U", kids, id=f"{FID}.U{next(nid)}", **T2)
        n.name = name
        return n

    # vertebrate-specific globins
    lam = G(lambda x: x["species"] == "PETMA")
    lam_node = Un([Node(g) for g, _ in lam], "cyclostome_Hb")
    mbgbe = Dn([label_tree("Mb"), label_tree("GbE")], ["mbgbe"], "Mb|GbE duplication predates gnathostome radiation",
               "Hoffmann 2010 Discussion", lvl_lo="Gnathostomata", lvl_hi="Vertebrata")
    r2 = Dn([label_tree("Cygb"), mbgbe], ["2R"], "Cygb|Mb/GbE = 2R (interval Gnathostomata..Vertebrata)",
            "Hoffmann 2012a Discussion", lvl_lo="Gnathostomata", lvl_hi="Vertebrata", tag="WGD")
    side_a = Un([lam_node, r2], "Cygb_Mb_GbE_paralogon")  # lamprey Hb: Cygb ortholog (Hoffmann) vs contested 2R
    ab = Dn([label_tree("Hba"), label_tree("Hbb")], ["ab", "ab_gn"], "alpha|beta-Hb tandem duplication before the "
            "gnathostome ancestor", "Hoffmann 2012a Discussion", lvl_lo="Gnathostomata", lvl_hi="Vertebrata",
            tag="TAND")
    side_b = Un([label_tree("GbY"), ab], "Hb_GbY_paralogon")
    vsg = Dn([side_a, side_b], ["1R", "four", "not_orth"], "1R: (Mb/GbE/Cygb) vs (alpha/beta-Hb + GbY); lamprey Hb not "
             "orthologous to gnathostome Hb", "Hoffmann 2012a Discussion; Hoffmann 2010 Results", lvl="Vertebrata",
             tag="WGD")
    vsg.name = "vertebrate_specific"
    amph = lambda nums: G(lambda x: x["species"] == "BRAFL" and re.sub(r"\D", "", x["gtype"]) in nums)
    c3 = {"3", "6", "12", "13", "14"}
    c1_amph = [Node(g) for g, _ in G(lambda x: x["species"] == "BRAFL" and x["gtype"].startswith("Gb")
                                     and re.sub(r"\D", "", x["gtype"]) not in c3 | {"4", "8"})]
    ciona = [Node(g) for g, _ in G(lambda x: x["species"] == "CIOIN")]
    clade1 = Un(c1_amph + ciona + [vsg], "clade1")
    A(clade1.features["id"], "event", "clade 1 = vertebrate-specific + urochordate + two amphioxus groups (internal "
      "relationships U)", "clade1", "Hoffmann 2012b Results", tier="T2")
    clade2 = label_tree("GbX")
    clade3 = Un([Node(g) for g, _ in amph(c3)], "clade3_amphioxus")
    ngb = label_tree("Ngb")
    gb4 = amph({"4"})
    clade4 = mk("S", [Node(gb4[0][0]), ngb], tier="T1", src="Hoffmann2012b") if gb4 and ngb is not None else (ngb or None)
    if clade4 is not None:
        clade4.name = "clade4_Ngb"
    d23 = Dn([clade2, clade3], ["clade23", "clade_amb"], "clade 2 (GbX) | clade 3 (amphioxus Gb3,6,12-14)",
             "Hoffmann 2012b Results", lvl_lo="Chordata", lvl_hi="ABOVE_RANGE")
    d123 = Dn([clade1, d23], ["clades", "clade_amb"], "clade 1 | clades 2+3", "Hoffmann 2012b Results",
              lvl_lo="Chordata", lvl_hi="ABOVE_RANGE")
    root = Dn([clade4, d123], ["clades", "clade4", "pomt2"], "clade 4 (Ngb + amphioxus Gb4) | clades 1-3",
              "Hoffmann 2012b Results", lvl_lo="Chordata", lvl_hi="ABOVE_RANGE")
    gb8 = amph({"8"})
    if gb8:
        root = Un([root, Node(gb8[0][0])], "")
        A(root.features["id"], "event", "amphioxus Gb8 not placed in the four clades -> U", "gb8",
          "Hoffmann 2012b Results", tier="T2")

    # ---------------------------------------------------------------- outputs
    members = pd.DataFrame([dict(gene_id=x["gene_id"], species=x["species"], original_id_in_paper=x["acc"],
                                 label=f"{x['label']} [{x['gtype']}]", mapping_method=x["mapping_method"],
                                 mapping_detail=x.get("mapping_detail", ""), mapping_identity=x["mapping_identity"],
                                 mapping_coverage=x["mapping_coverage"], source="Hoffmann2012b",
                                 locator=f"Supplementary Table S1 (PDF text line {x['line']})") for x in mapped])
    members.to_csv(HERE / "members.tsv", sep="\t", index=False)
    rep = pd.DataFrame([dict(original_id=x["acc"], species=x["species"], label=x["label"], status=x["status"],
                             gene_id=x["gene_id"], reason=x["reason"], locator=f"Table S1 line {x['line']}")
                        for x in recs])
    rep.to_csv(HERE / "mapping_report.tsv", sep="\t", index=False)
    (HERE / "reference_tree.nhx").write_text(vlib.nhx(root) + "\n")
    pd.DataFrame(assertions).to_csv(HERE / "assertions.tsv", sep="\t", index=False)
    inpanel = rep[rep.status != "excluded"]
    n_in, n_map = len(inpanel), int((inpanel.status == "mapped").sum())
    per_sp = {k: f"{int((d.status == 'mapped').sum())}/{len(d)}" for k, d in inpanel.groupby("species")}
    f = lambda p: {str(p.relative_to(R)): vlib.sha256(p)}
    meta = dict(
        id=FID, name="Chordate globins (Ngb, GbX, Cygb, Mb, GbE, GbY, alpha/beta-Hb, cyclostome Hb)", panel="V",
        phenomena=["WGD", "TAND", "LOSS", "CONV", "DEEP"],
        descriptors=dict(n_members=len(mapped), n_species=int(members.species.nunique()),
                         mapping=f"{n_map}/{n_in} in-panel Table S1 rows mapped ({n_map / max(n_in, 1):.0%})",
                         per_species_mapped=per_sp),
        sources=[
            dict(id="Hoffmann2012b", doi="10.1093/molbev/mss018", pmid="22319164",
                 title="Evolution of the globin gene family in deuterostomes: lineage-specific patterns of "
                       "diversification and attrition", year=2012, source_method="phylogeny",
                 table="Supplementary Table S1", locator="Table S1; Results (four globin clades)", quote=Q["clades"][1],
                 files={**f(SUPP_PDF), **f(SUPP_TXT), **f(ALN), **f(P12B / "raw/fulltext.txt"),
                        **f(P12B / "raw/received/04_supplement.zip")}),
            dict(id="Hoffmann2012a", doi="10.1093/molbev/msr207", pmid="21965344",
                 title="Whole-genome duplications spurred the functional diversification of the globin gene "
                       "superfamily in vertebrates", year=2012, source_method="synteny", locator="Discussion",
                 quote=Q["1R"][1], files=f(P12A / "raw/fulltext.txt")),
            dict(id="Hoffmann2010", doi="10.1073/pnas.1006756107", pmid="20660759",
                 title="Gene cooption and convergent evolution of oxygen transport hemoglobins in jawed and jawless "
                       "vertebrates", year=2010, source_method="phylogeny", locator="Results/Discussion",
                 quote=Q["not_orth"][1], files=f(P10 / "raw/fulltext.txt")),
        ],
        root_level="Chordata", root_single_copy=False, complete_species=[], decoys=[], distinct_from=[],
        tier_summary=dict(T1_events=sum(1 for a in assertions if a["tier"] == "T1")),
        protocol_version=1,
        notes=("Members are the representative sequences of Hoffmann 2012b Table S1 (one alpha- and one beta-Hb per "
               "gnathostome species, not the full clusters) -> no complete_species. Hb cluster tandem copies, Opazo "
               "2008 embryonic genes and teleost 3R Hb clusters (Opazo 2013) are not encoded. Lamprey Hbs (Table S1 "
               "sea-lamprey entries) placed under U with the Cygb/Mb/GbE subtree because their orthology to Cygb "
               "(Hoffmann 2010) vs co-orthology to Cygb+Mb+GbE depends on the contested 2R timing; the scored claim "
               "is lamprey Hb vs gnathostome alpha/beta-Hb/GbY = paralogs (1R, Vertebrata). Hoffmann 2012a Table S1 "
               "and Hoffmann 2010 sd01 were not received; the received item #04 is Hoffmann 2012b."),
    )
    (HERE / "family.yaml").write_text(yaml.safe_dump(meta, sort_keys=False, allow_unicode=True, width=120))
    print(f"{FID}: {n_map}/{n_in} mapped; {per_sp}")


def _embl(oid):
    from curation_tools.map_ids import run_job
    try:
        rows, _ = run_job("EMBL-GenBank-DDBJ", "UniProtKB", [oid])
    except Exception:
        return []
    return [r.get("Entry", "") for r in rows]


def _resolve(cache, sp, recs):
    # paper alignment sequences ('aln:' keys) are already in the cache; vlib fetches UniProt candidates
    vlib.resolve_by_sequence(cache, sp, recs, strict_family=True)  # TAND-tagged family -> 99 % identity


if __name__ == "__main__":
    main()
