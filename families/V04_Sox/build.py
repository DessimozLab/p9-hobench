#!/usr/bin/env python
"""V04_Sox: vertebrate Sox genes (groups B1, B2, C, D, E, F, H) from Voldoire et al. 2017 (PLoS ONE,
PMID 28738066, PMC5524304), supplementary table S6 (pone.0180936.s006.xls: per-gene x per-species
accessions, columns = sox genes incl. teleost a/b ohnologs).

Rules (CURATION_PROTOCOL v1):
* members = every cell of S6 for a panel species; mapped to the pinned UniProt panel proteome by
  exact ID mapping (Ensembl gene/protein, RefSeq) or, for spotted gar (genomic coordinates in S6),
  by Ensembl REST overlap of the coordinates with exactly one protein-coding gene, then Ensembl
  gene -> UniProt. Genomic contig/scaffold/EST/TIGR-TC identifiers are not gene models -> unmapped.
* tree (T2 classification + T1 text events):
    root            D  (groups present in all bilaterians -> split above the panel root)
    group (B1..H)   D  2R interval (lvl_lo=Gnathostomata, lvl_hi=Vertebrata; "probable consequence
                       of the two rounds of WGD"); sox17/sox32 SSD at Teleostei inside sox17.
    gene (sox1..)   S over non-teleosts + D(Teleostei, WGD) over teleost a/b copies (synteny-curated);
                    singleton teleost genes are S (paper: same ortholog retained, synteny S3 Fig).
Run:  python families/V04_Sox/build.py   (uses cache/ for network lookups)
"""
from __future__ import annotations

import csv
import re
import sys
from pathlib import Path

import pandas as pd
import yaml

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))  # vlib.py: shared Panel V helpers (identical copy in each V0x family)
import vlib  # noqa: E402
from vlib import mk, species_s_tree  # noqa: E402

FID = "V04_Sox"
SRC = "Voldoire2017"
PAPER = vlib.ROOT / "literature/papers/PMID28738066"
S6 = vlib.ROOT / "literature/raw/B1/supp/sox_Voldoire2017__PMC5524304/pone.0180936.s006.xls"
TEXT = (PAPER / "fulltext.txt").read_text()

SPECIES = {"Elephant shark": "CALMI", "Mouse": "MOUSE", "Human": "HUMAN", "Coelacanth": "LATCH",
           "Spotted gar": "LEPOC", "Zebrafish": "DANRE", "Salmon": "SALSA", "Trout": "ONCMY",
           "Medaka": "ORYLA",
           # not in panel V
           "Cavefish": None, "Cod": None, "Tilapia": None, "Platy": None, "Tetraodon": None,
           "Fugu": None, "Stickleback": None}
TELEOSTS = {"DANRE", "ORYLA", "SALSA", "ONCMY", "ESOLU", "CORBA"}
ENS_SPECIES = {"LEPOC": "lepisosteus_oculatus"}

# group of each sox gene: from the S6 block headers (SOXB1, SOXB2, ...)
Q = {  # quotes (verbatim from fulltext.txt)
    "groups_bilaterian": "The seven other groups mentioned above are present in all bilaterians, generally "
                         "represented by a single sox gene per group in invertebrates",
    "2R": "vertebrate genomes contain for all groups, except group H represented only by one gene, two to four "
          "sox genes, as a probable consequence of the two rounds of WGD that occurred at the base of vertebrates",
    "nineteen": "These nineteen sox genes are also present as singletons in the three other Sarcopterygii "
                "analyzed (coelacanth, mouse and human).",
    "genomewide": "Sox genes were searched genome-wide (see Materials and methods section) in the eleven most "
                  "complete teleost genomes (Fig 1). The same analyses were performed in five non-teleost "
                  "species (elephant shark, mouse, human, coelacanth and spotted gar)",
    "3R": "On the 10 sox genes duplicated through the teleost-specific WGD (sox1, 4, 6, 8, 9, 10, 11, 14, 19 "
          "and 21, Fig 1)",
    "3R_ancestor": "phylogenetic reconstructions showed that both paralogs are derived from a duplication event "
                   "that occurred in the common ancestor of all studied teleost species (Fig 3)",
    "sox32": "the so-called sox32 gene in teleosts is a duplicate of sox17 that probably arose through SSD: "
             "both genes are localized on the same chromosome (S2 Fig)",
    "sox32_anc": "The identification of sox17/32 in all teleost genomes suggests that the event of SSD took "
                 "place in a common ancestor of the analyzed species.",
    "same_copy": "we observed that the same copy (i.e. the same ortholog) has been kept after the WGD (no "
                 "reciprocal gene loss) (S3 Fig)",
    "synteny_a_b": "Soxa and soxb paralogs are respectively in blue and orange.",
    "shark_sox19": "the elephant shark, a cartilaginous fish, seems to have lost sox19",
}
for k, q in Q.items():
    assert vlib.quote_in(TEXT, q), f"quote {k} not found in fulltext"

GENE_ID_RE = re.compile(r"^(ENS[A-Z]*[GPT]\d{11}|NM_\d+|Q[0-9A-Z]{5}_[A-Z]+)$")


def parse_s6():
    raw = pd.read_excel(S6, header=None, dtype=str).fillna("")
    rows = []
    blocks = [i for i in range(len(raw)) if any(str(c).startswith("SOX") for c in raw.iloc[i])]
    for b in blocks:
        hdr_grp = raw.iloc[b]
        hdr = raw.iloc[b + 1]
        # group label per column: last SOX label at or left of column
        grp_of, cur = {}, None
        for j in range(raw.shape[1]):
            if str(hdr_grp[j]).startswith("SOX"):
                cur = str(hdr_grp[j]).strip()
            grp_of[j] = cur
        r = b + 2
        while r < len(raw) and not any(str(c).startswith("SOX") for c in raw.iloc[r]):
            row = raw.iloc[r]
            sp_label = str(row[0]).strip()
            for j in range(1, raw.shape[1]):
                col = str(hdr[j]).strip()
                if not col:  # species-label column of the right-hand block
                    if str(row[j]).strip():
                        sp_label = str(row[j]).strip()
                    continue
                cell = str(row[j]).strip()
                if not cell:
                    continue
                rows.append(dict(row_xls=r + 1, col_xls=j + 1, group=grp_of[j], column=col,
                                 species_label=sp_label.strip(), cell=cell))
            r += 1
    return rows


def gene_and_copy(column: str):
    """'sox1/sox1a' -> ('sox1', 'a', True); 'sox1b' -> ('sox1','b',False); 'sox2' -> ('sox2','',True);
    'sox15/sox19a' -> ('sox19','a',True) (mammalian SOX15 = sox19 per S6 header)."""
    parts = [p.strip() for p in column.split("/")]
    teleost_name = parts[-1]
    m = re.match(r"^(sox\d+)([ab]?)$", teleost_name)
    assert m, column
    return m.group(1), m.group(2), (len(parts) > 1 or not m.group(2))


def split_cell(cell: str):
    """Salmonid cells hold several copies '(a1) x, (a2) y'. Returns [(copy_tag, id)]."""
    if "(" in cell:
        out = []
        for m in re.finditer(r"\(([^)]*)\)\s*([^,(]+)", cell):
            out.append((m.group(1).strip(), m.group(2).strip()))
        return out
    return [("", cell.strip())]


def main():
    st = vlib.species_tree()
    cache = vlib.JsonCache(HERE / "cache" / "lookups.json")
    rows = parse_s6()
    recs = []
    for r in rows:
        code = SPECIES[r["species_label"].strip()]
        gene, copy, nonteleost_col = gene_and_copy(r["column"])
        for ctag, oid in split_cell(r["cell"]):
            recs.append(dict(r, species=code, gene=gene, copy=copy, nonteleost_col=nonteleost_col,
                             salmonid_copy=ctag, original_id=oid))
    # ---------------------------------------------------------------- mapping
    for x in recs:
        code, oid = x["species"], x["original_id"]
        x.update(gene_id="", status="", reason="", mapping_method="", mapping_identity="", mapping_coverage="")
        if code is None:
            x.update(status="not_in_panel", reason=f"{x['species_label']} not in panel V")
            continue
        m = re.match(r"^Chr:([^:]+):(\d+):(\d+)$", oid)
        region = f"{m.group(1)}:{m.group(2)}-{m.group(3)}" if m else None
        if vlib.id_type(oid) or region:
            x.update(vlib.map_one(cache, code, oid, region))
        else:
            x.update(status="unmapped",
                     reason="identifier is a genomic contig/scaffold/EST/TIGR-TC, not a gene model; no exact "
                            "mapping possible (sequence mapping would fail the coverage rule)")
    for code in sorted({x["species"] for x in recs if x["species"]}):
        vlib.resolve_by_sequence(cache, code, [x for x in recs if x["species"] == code])
    # a single paper gene mapped to the same accession twice -> ambiguous (e.g. one Ensembl gene in two cells)
    seen: dict[str, list] = {}
    for x in recs:
        if x["gene_id"]:
            seen.setdefault(x["gene_id"], []).append(x)
    for acc, xs in seen.items():
        if len(xs) > 1:
            for x in xs:
                x.update(status="ambiguous", reason=f"accession {acc} assigned to {len(xs)} S6 cells",
                         gene_id="")
    mapped = [x for x in recs if x["status"] == "mapped"]

    # ---------------------------------------------------------------- tree
    groups: dict[str, dict[str, dict[str, list]]] = {}
    for x in mapped:
        g = groups.setdefault(x["group"], {}).setdefault(x["gene"], {"non": [], "a": [], "b": [], "single": []})
        if x["species"] not in TELEOSTS:
            g["non"].append((x["gene_id"], x["species"]))
        elif x["copy"]:
            g[x["copy"]].append((x["gene_id"], x["species"]))
        else:
            g["single"].append((x["gene_id"], x["species"]))
    nid = iter(range(1, 10000))
    assertions = []

    def add_assert(node_id, typ, target, tier, locator, qkey):
        assertions.append(dict(assertion_id=f"{FID}.A{len(assertions) + 1:03d}", type=typ,
                               target=node_id if typ == "event" else target, note=target,
                               tier=tier, source_id=SRC, locator=locator, quote=Q[qkey]))

    def gene_tree(gene, g):
        teleo = None
        if g["a"] or g["b"]:
            ca = species_s_tree(st, g["a"], {"tier": "T2", "src": SRC}) if g["a"] else None
            cb = species_s_tree(st, g["b"], {"tier": "T2", "src": SRC}) if g["b"] else None
            if ca is not None and cb is not None:
                i = f"{FID}.D{next(nid)}"
                teleo = mk("D", [ca, cb], id=i, lvl="Teleostei", tag="WGD", tier="T1", src=SRC)
                add_assert(i, "event", f"{gene}a|{gene}b split = teleost 3R (lvl=Teleostei)", "T1",
                           "Results 'Evolutionary fate...' para 1; S6 columns", "3R")
                add_assert(i, "event", f"{gene}a|{gene}b ancestor = common ancestor of studied teleosts",
                           "T1", "Results 'Origin of sox gene duplicates' para 1", "3R_ancestor")
            else:
                teleo = ca if ca is not None else cb
        all_single = g["single"]
        kids = []
        # non-teleosts + teleost subtree under species-tree S nodes
        if teleo is not None:
            # combine: S-tree over non-teleosts + singles, with teleost block attached at Neopterygii level
            base = g["non"] + all_single
            if base:
                outer = species_s_tree(st, base + [("__TELEO__", "DANRE")], {"tier": "T2", "src": SRC})
                for nd in outer.preorder():
                    for k, c in enumerate(nd.children):
                        if c.name == "__TELEO__":
                            nd.children[k] = teleo
                            teleo.parent = nd
                if outer.name == "__TELEO__":
                    outer = teleo
                return outer
            return teleo
        return species_s_tree(st, g["non"] + all_single, {"tier": "T2", "src": SRC})

    group_nodes = []
    for grp in sorted(groups):
        gene_nodes = []
        for gene in sorted(groups[grp], key=lambda s: int(s[3:])):
            g = groups[grp][gene]
            if gene == "sox17":
                continue
            node = gene_tree(gene, g)
            if node is not None and not node.is_leaf:
                node.name = gene
            gene_nodes.append(node)
        if grp == "SOXF":
            # sox17 / sox32: SSD in the teleost ancestor (S2 Fig); non-teleost sox17 orthologous to both
            g17 = groups[grp].get("sox17", {"non": [], "single": []})
            g32 = groups[grp].get("sox32", {"non": [], "single": []})
            t17 = [x for x in g17["single"]]
            t32 = [x for x in g32["single"]]
            gene_nodes = [n for n in gene_nodes if n is None or n.name != "sox32"]
            i = f"{FID}.D{next(nid)}"
            teleo = mk("D", [species_s_tree(st, t17, {"tier": "T2", "src": SRC}) if t17 else None,
                             species_s_tree(st, t32, {"tier": "T2", "src": SRC}) if t32 else None],
                       id=i, lvl="Teleostei", tag="SSD", tier="T1", src=SRC)
            add_assert(i, "event", "sox17|sox32 = SSD in teleost common ancestor (lvl=Teleostei)", "T1",
                       "Results 'Origin of sox gene duplicates' para 1", "sox32")
            add_assert(i, "event", "sox17|sox32 SSD in common ancestor of analyzed teleosts", "T1",
                       "Results 'Origin of sox gene duplicates' para 1", "sox32_anc")
            base = g17["non"]
            outer = species_s_tree(st, base + [("__TELEO__", "DANRE")], {"tier": "T2", "src": SRC})
            for nd in outer.preorder():
                for k, c in enumerate(nd.children):
                    if c.name == "__TELEO__":
                        nd.children[k] = teleo
            outer.name = "sox17_32"
            gene_nodes.append(outer)
        gene_nodes = [n for n in gene_nodes if n is not None]
        if len(gene_nodes) > 1:
            i = f"{FID}.D{next(nid)}"
            node = mk("D", gene_nodes, id=i, lvl_lo="Gnathostomata", lvl_hi="Vertebrata", tag="WGD",
                      tier="T1", src=SRC)
            node.name = grp
            add_assert(i, "event", f"{grp} paralogs from 2R (interval Gnathostomata..Vertebrata)", "T1",
                       "Introduction para 4", "2R")
        else:
            node = gene_nodes[0]
        group_nodes.append(node)
    i = f"{FID}.D{next(nid)}"
    root = mk("D", group_nodes, id=i, lvl_lo="Chordata", lvl_hi="ABOVE_RANGE", tier="T1", src=SRC)
    add_assert(i, "event", "Sox groups B1,B2,C,D,E,F,H split before Bilateria (above panel root)", "T1",
               "Introduction para 4", "groups_bilaterian")
    add_assert("", "completeness", "CALMI,MOUSE,HUMAN,LATCH,LEPOC,DANRE,ORYLA,SALSA,ONCMY", "T1",
               "Results para 1", "genomewide")

    # ---------------------------------------------------------------- outputs
    members = pd.DataFrame([dict(gene_id=x["gene_id"], species=x["species"], original_id_in_paper=x["original_id"],
                                 label=x["gene"] + x["copy"], mapping_method=x["mapping_method"],
                                 mapping_detail=x.get("mapping_detail", ""),
                                 mapping_identity=x["mapping_identity"], mapping_coverage=x["mapping_coverage"],
                                 source=SRC, locator=f"S6 row {x['row_xls']} col {x['col_xls']} ({x['column']})")
                            for x in mapped])
    members.to_csv(HERE / "members.tsv", sep="\t", index=False)
    rep = pd.DataFrame([dict(original_id=x["original_id"], species_label=x["species_label"],
                             species=x["species"] or "", label=x["gene"] + x["copy"] +
                             (f"({x['salmonid_copy']})" if x["salmonid_copy"] else ""),
                             status=x["status"], gene_id=x["gene_id"], reason=x["reason"],
                             locator=f"S6 row {x['row_xls']} col {x['col_xls']}") for x in recs])
    rep.to_csv(HERE / "mapping_report.tsv", sep="\t", index=False)
    (HERE / "reference_tree.nhx").write_text(vlib.nhx(root) + "\n")
    pd.DataFrame(assertions).to_csv(HERE / "assertions.tsv", sep="\t", index=False)

    inpanel = rep[rep.species != ""]
    n_in, n_map = len(inpanel), int((inpanel.status == "mapped").sum())
    gene_level = inpanel[~inpanel.reason.str.startswith("identifier is a genomic")]
    per_sp = inpanel.groupby("species").apply(lambda d: f"{int((d.status == 'mapped').sum())}/{len(d)}").to_dict()
    complete = [sp for sp in ("HUMAN", "MOUSE", "LATCH", "LEPOC", "DANRE", "ORYLA", "CALMI")
                if sp in per_sp and (lambda a, b: b > 0 and a / b >= 0.9)(*map(int, per_sp[sp].split("/")))]
    meta = dict(
        id=FID, name="Sox transcription factors (groups B1, B2, C, D, E, F, H; SRY excluded)", panel="V",
        phenomena=["WGD", "SSD", "LOSS", "DEEP"],
        descriptors=dict(n_members=len(mapped), n_species=int(members.species.nunique()),
                         mapping=f"{n_map}/{n_in} in-panel paper entries mapped "
                                 f"({n_map / n_in:.0%}); {n_map}/{int((gene_level.status != 'x').sum())} "
                                 f"of gene-level identifiers",
                         per_species_mapped=per_sp),
        sources=[dict(id=SRC, doi="10.1371/journal.pone.0180936", pmid="28738066",
                      title="Expansion by whole genome duplication and evolution of the sox gene family in teleost fish",
                      year=2017, source_method="synteny",
                      table="S6 Table (pone.0180936.s006.xls)", locator="S6 Table; Results; Fig 1-3; S2 Fig",
                      quote=Q["3R"],
                      files={str(S6.relative_to(vlib.ROOT)): vlib.sha256(S6),
                             str((PAPER / 'fulltext.txt').relative_to(vlib.ROOT)): vlib.sha256(PAPER / 'fulltext.txt')})],
        root_level="Chordata", root_single_copy=False,
        complete_species=complete,
        decoys=[], distinct_from=[],
        tier_summary=dict(T1_events=sum(a["type"] == "event" for a in assertions), T2="gene-level S polytomies "
                          "(S6 columns; synteny-curated a/b ohnologs)"),
        protocol_version=1,
        notes=("2R placement is an interval (lvl_lo=Gnathostomata, lvl_hi=Vertebrata) because the paper says "
               "'probable consequence of the two rounds of WGD'; no lamprey in source. Salmonid (salmon, trout) and "
               "elephant-shark entries of S6 are genomic scaffold/contig IDs from pre-2016 assemblies and cannot be "
               "mapped by exact ID; they are listed as unmapped (so Ss4R copies are absent from the tree). "
               "S6 header for SOXE reads 'sox8/sox8a' but non-teleost sox8 sits in the 'sox8b' column; irrelevant to "
               "the tree (non-teleosts are pre-3R), labels follow the column header. Group splits are placed at the "
               "panel root (Chordata) with lvl_hi=ABOVE_RANGE (present in all bilaterians)."),
    )
    (HERE / "family.yaml").write_text(yaml.safe_dump(meta, sort_keys=False, allow_unicode=True, width=120))
    print(f"{FID}: {len(mapped)} mapped / {n_in} in-panel entries; per species {per_sp}; complete {complete}")


if __name__ == "__main__":
    main()
