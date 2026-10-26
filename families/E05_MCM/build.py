#!/usr/bin/env python
"""E05_MCM: build members.tsv, mapping_report.tsv, reference_tree.nhx, assertions.tsv, family.yaml.

Source: Liu Y, Richards TA, Aves SJ (2009) Ancient diversification of eukaryotic MCM DNA replication proteins.
BMC Evol Biol 9:60 (PMID 19292915, PMC2667178, CC BY). Paralogue labels and accessions are read from the vector PDF
trees Additional files 1-4 (extract.py: one block of leaves per paralogue, caption order, checked against the printed
paralogue labels). Mapping (prep_mapping.py caches): NCBI accessions -> fetched sequence -> same-species DIAMOND
(protocol thresholds); JGI / CyanidioschyzonDB / GiardiaDB ids -> exact ORF-name match in the panel header_map
(unique hit). Every mapped target must hit the paper's human MCMs with e <= 1e-20 (paper's BLAST cut-off; QC of
family membership only, never used for the paralogue label).
Run from project root: python families/E05_MCM/build.py
"""
import csv, hashlib, sys
from collections import Counter, defaultdict
from pathlib import Path
import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[2]
FAM = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(FAM))
from extract import PAPER, all_leaves, header_index, id_pattern  # noqa: E402
from hobench.newick import Node, write_newick  # noqa: E402

SRC = "Liu2009"
TXT_REL = f"{PAPER}/fulltext.txt"
SUPP_REL = [f"{PAPER}/raw/supp/1471-2148-9-60-S{i}.pdf" for i in (1, 2, 3, 4)] + [f"{PAPER}/raw/supp/1471-2148-9-60-S7.txt"]
ID_MIN, COV_MIN, MARGIN, QC_E = 98.0, 0.80, 2.0, 1e-20
Q = dict(
    LECA="demonstrating that Mcm2-9 were formed by seven gene duplication events before the last common ancestor of the eukaryotes",
    ARC="Our analyses therefore suggest that the eukaryotic Mcm2-9 may be derived by gene duplication events from a single, archaeal-like, ancestral MCM",
    ARC1="we identified a single highly conserved MCM homologue in all the Archaea studied (except Methanococcus jannaschii which possesses four MCM proteins)",
    NG="In Naegleria and Xenopus there were two proteins that grouped with Mcm3",
    RECENT="providing the only examples of recent MCM duplication events among the taxa investigated",
    CLASS="Both phylogenetic analyses classified the MCM homologues as belonging to specific paralogue families",
    CTRL="Mcm2-7 protein paralogues were present in all eukaryote genomes studied",
    LOSS="Together this suggests at least five separate cases of paired Mcm8 and 9 loss during the eukaryotic radiation",
    M89LECA="Mcm8 and Mcm9 were also present in the LCEA",
    DM="Drosophila was the only eukaryote analysed which has Mcm8 but lacks Mcm9",
    COMPLETE="These multiple BLAST searches using paralogous genes as search seeds were conducted to ensure that all available paralogues were sampled from each of the target genomes.",
)
PARALOGS = ["Mcm2", "Mcm3", "Mcm4", "Mcm5", "Mcm6", "Mcm7", "Mcm8", "Mcm9"]


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main():
    leaves = [d for d in all_leaves() if d["code"]]
    allacc = Counter((d["acc"].split("(")[0], d["paralog"]) for d in all_leaves())
    acc_paralogs = defaultdict(set)
    for (a, p) in allacc:
        acc_paralogs[a].add(p)
    qm = pd.read_csv(FAM / "cache/query_meta.tsv", sep="\t", dtype=str).fillna("")
    qmeta = {(r.code, r.orig): r for r in qm.itertuples()}
    fw = pd.read_csv(FAM / "cache/fwd.tsv", sep="\t")
    rv = pd.read_csv(FAM / "cache/rev.tsv", sep="\t")
    fwg = {k: g.sort_values(["bitscore", "pident"], ascending=False) for k, g in fw.groupby("qseqid")}
    rvb = {t: g.sort_values(["bitscore", "pident"], ascending=False) for t, g in rv.groupby("qseqid")}
    qc = pd.read_csv(FAM / "cache/qc.tsv", sep="\t")
    qc_e = {r.qseqid: (float(r.evalue), r.sseqid) for r in qc.itertuples()}
    hidx = {}
    im = pd.read_csv(FAM / "cache/idmap.tsv", sep="\t", dtype=str).fillna("")
    idm = defaultdict(list)
    for x in im.itertuples():
        idm[x.original_id].append(x.gene_name)
    # isoform groups: same species + same paralogue + several leaves (human Mcm7/Mcm8 isoforms)
    group = defaultdict(list)
    for d in leaves:
        group[(d["code"], d["paralog"])].append(d["acc"])

    rep = []
    for d in leaves:
        code, par, acc = d["code"], d["paralog"], d["acc"]
        a0 = acc.split("(")[0]
        loc = f"{d['file']} ({d['figure']}), {par} block, leaf '{d['species']} [{acc}]' y={d['y']}"
        base = dict(original_id=acc, species=code, label=par, locator=loc)
        if code == "HALSA":
            base["label"] = "archaeal_MCM"
        if len(acc_paralogs[a0]) > 1:
            rep.append({**base, "gene_id": "", "status": "ambiguous",
                        "reason": f"accession printed under several paralogues in Additional files 1-4: {sorted(acc_paralogs[a0])} (source error)"})
            continue
        tok, db = id_pattern(code, acc)
        if tok:
            if code not in hidx:
                hidx[code] = header_index(code)
            hits = sorted(hidx[code][0].get(tok, ()))
            if len(hits) == 1:
                rep.append({**base, "gene_id": hits[0], "status": "mapped",
                            "reason": f"{db} -> UniProt ORF/gene name {tok.replace('#', '')} (unique, {hidx[code][1].relative_to(ROOT)})",
                            "_m": dict(mapping_method="id_mapping_orf_name", mapping_identity="", mapping_coverage="",
                                       mapping_margin="", rbh="")})
            else:
                rep.append({**base, "gene_id": "", "status": "unmapped" if not hits else "ambiguous",
                            "reason": f"{db}: ORF name {tok.replace('#', '')} has {len(hits)} panel entries {hits}"})
            continue
        r = qmeta.get((code, acc))
        if r is None or not r.qid:
            if acc.startswith("BDEG_"):
                rep.append({**base, "gene_id": "", "status": "not_in_panel",
                            "reason": "source genome = B. dendrobatidis JEL423 (Broad BDEG_ locus tags); panel BATDJ = strain "
                                      "JAM81 (UP000007241); no JEL423 sequence or cross-strain identifier is retrievable"})
            else:
                rep.append({**base, "gene_id": "", "status": "unmapped", "reason": "no NCBI record / identifier route"})
            continue
        g = fwg.get(r.qid)
        if g is None or not len(g):
            rep.append({**base, "gene_id": "", "status": "unmapped",
                        "reason": f"no same-species DIAMOND hit ({r.ncbi_record}: {r.ncbi_header})"})
            continue
        b = g.iloc[0]
        others = g[g.sseqid != b.sseqid]
        runner = float(others.iloc[0].pident) if len(others) else 0.0
        margin = float(b.pident) - runner
        qcov, tcov = b.length / b.qlen, b.length / b.slen
        rg = rvb.get(b.sseqid)
        sibs = {f"{code}|{x}" for x in group[(code, par)]}
        rbh = bool(rg is not None and rg.iloc[0].sseqid in sibs and
                   (rg[rg.bitscore == rg.iloc[0].bitscore].sseqid.isin(sibs)).all())
        info = (f"{r.route} {r.ncbi_record}; best {b.sseqid} {b.pident:.1f}% qcov {qcov:.2f} tcov {tcov:.2f} "
                f"margin {margin:.1f} rbh {rbh}")
        if b.pident >= ID_MIN and min(qcov, tcov) >= COV_MIN and margin >= MARGIN and rbh:
            rep.append({**base, "gene_id": b.sseqid, "status": "mapped", "reason": info,
                        "_m": dict(mapping_method="sequence_diamond" if r.route == "ncbi_accession" else "sequence_diamond_via_old_locus_tag",
                                   mapping_identity=f"{b.pident:.1f}", mapping_coverage=f"{min(qcov, tcov):.3f}",
                                   mapping_margin=f"{margin:.1f}", rbh="1")})
        elif b.pident >= ID_MIN and min(qcov, tcov) >= COV_MIN:
            rep.append({**base, "gene_id": "", "status": "ambiguous", "reason": f"fails margin/RBH; {info}"})
        elif (b.pident >= ID_MIN and qcov >= 0.90 and margin >= MARGIN and rbh and "isoform" in r.ncbi_header
              and len(group[(code, par)]) == 1):
            # the paper used a shorter RefSeq isoform: map by gene name (UniProt ID mapping of the RefSeq accession ->
            # one gene name; that gene name has a single panel entry), corroborated by >=98% identity over >=90% of query
            gn = {x for x in idm.get(r.ncbi_record.split(".")[0], []) if x}
            if code not in hidx:
                hidx[code] = header_index(code)
            hits = sorted(hidx[code][0].get(next(iter(gn)), ())) if len(gn) == 1 else []
            if len(hits) == 1 and hits[0] == b.sseqid:
                rep.append({**base, "gene_id": hits[0], "status": "mapped",
                            "reason": f"UniProt ID mapping {r.ncbi_record} -> gene {next(iter(gn))} (unique panel entry); {info} (source = shorter isoform)",
                            "_m": dict(mapping_method="gene_name_unique", mapping_identity=f"{b.pident:.1f}",
                                       mapping_coverage=f"{qcov:.3f}", mapping_margin=f"{margin:.1f}", rbh="1")})
            else:
                rep.append({**base, "gene_id": "", "status": "unmapped", "reason": f"isoform; gene-name route not unique {sorted(gn)} {hits}; {info}"})
        else:
            rep.append({**base, "gene_id": "", "status": "unmapped", "reason": f"below thresholds; {info} ({r.ncbi_header})"})

    # isoform leaves of one gene (same species + paralogue) whose sibling leaf is mapped to the same target
    mapped_t = {(x["species"], x["label"]): x for x in rep if x["status"] == "mapped"}
    for x in rep:
        sib = mapped_t.get((x["species"], x["label"]))
        if x["status"] == "unmapped" and sib and f"best {sib['gene_id']} " in x["reason"] and "isoform" in x["reason"]:
            x.update(status="excluded", reason=x["reason"] + f"; isoform of the gene mapped via {sib['original_id']} ({sib['gene_id']}; protocol 0.6)")

    # homology QC (paper's BLAST cut-off) on every mapped target
    for x in rep:
        if x["status"] == "mapped":
            e = qc_e.get(f"{x['species']}|{x['gene_id']}")
            if e is None or e[0] > QC_E:
                x.update(status="excluded", reason=x["reason"] + f"; target fails MCM homology QC (best hit to paper's human MCMs: {e}) -> source accession error")
                x.pop("_m", None)
    # several source leaves -> one panel entry
    byg = defaultdict(list)
    for i, x in enumerate(rep):
        if x["status"] == "mapped":
            byg[x["gene_id"]].append(i)
    for gid, ii in byg.items():
        if len(ii) < 2:
            continue
        keys = {(rep[i]["species"], rep[i]["label"]) for i in ii}
        if len(keys) == 1:  # isoforms of one gene listed separately by the paper
            for i in ii[1:]:
                rep[i].update(status="excluded", reason=rep[i]["reason"] + f"; isoform of the same gene as {rep[ii[0]]['original_id']} ({gid}; protocol 0.6)")
                rep[i].pop("_m", None)
        else:
            for i in ii:
                rep[i].update(status="ambiguous", reason=rep[i]["reason"] + f"; {gid} also target of another paralogue leaf")
                rep[i].pop("_m", None)

    mem = [dict(gene_id=x["gene_id"], species=x["species"], original_id_in_paper=x["original_id"], label=x["label"],
                **x["_m"], source=SRC, locator=x["locator"]) for x in rep if x["status"] == "mapped"]
    members = pd.DataFrame(mem).sort_values(["label", "species", "gene_id"]).reset_index(drop=True)
    members.to_csv(FAM / "members.tsv", sep="\t", index=False, quoting=csv.QUOTE_MINIMAL)
    rp = pd.DataFrame([{k: x.get(k, "") for k in ("original_id", "species", "label", "gene_id", "status", "reason", "locator")} for x in rep])
    rp = rp.sort_values(["label", "species", "original_id"]).reset_index(drop=True)
    rp.to_csv(FAM / "mapping_report.tsv", sep="\t", index=False)

    # ------------------------------------------------------------- tree
    cnt = Counter()

    def node(nid, ev, tier, **kw):
        cnt[(ev, tier)] += 1
        return Node(nid, features={"id": nid, "ev": ev, "tier": tier, "src": SRC, **kw})

    root = node("MCM_ROOT", "S", "T1")
    euk = node("MCM_LECA_DUPS", "D", "T1", lvl_lo="Eukaryota", lvl_hi="Asgard_Eukaryota")
    for _, m in members[members.label == "archaeal_MCM"].iterrows():
        root.add_child(Node(m.gene_id))
    root.add_child(euk)
    for par in PARALOGS:
        mm = members[members.label == par]
        if mm.empty:
            continue
        pn = node(f"{par}_LECA", "S", "T2")
        for sp, g in mm.groupby("species"):
            if len(g) == 1:
                pn.add_child(Node(g.gene_id.iloc[0]))
            elif sp == "NAEGR" and par == "Mcm3":
                dn = node("Mcm3_NAEGR_dup", "D", "T1", tag="LSE")  # no lvl: LCA mapping = NAEGR (species-specific)
                for gid in sorted(g.gene_id):
                    dn.add_child(Node(gid))
                pn.add_child(dn)
            else:
                raise SystemExit(f"unexpected multi-copy {sp} {par}")
        if len(pn.children) == 1 and pn.children[0].children:
            pn = pn.children[0]
        euk.add_child(pn)
    (FAM / "reference_tree.nhx").write_text(write_newick(root, features=["id", "ev", "lvl", "lvl_lo", "lvl_hi", "tier", "src", "tag"], lengths=False) + "\n")

    # ------------------------------------------------------------- assertions
    inp = rp[rp.status.isin(["mapped", "unmapped", "ambiguous"])]  # excluded / not_in_panel are not in-panel members
    frac = inp.groupby("species").status.apply(lambda s: (s == "mapped").mean())
    # a species whose source leaf was excluded by the homology QC (paper accession error) is missing a real member
    qc_fail = set(rp[rp.reason.str.contains("homology QC")].species)
    complete = sorted(c for c in frac[frac >= 0.9].index if c != "HALSA" and c not in qc_fail)
    loc_res = f"{TXT_REL} Results"
    rows = [
        dict(assertion_id="E05_A1", type="event", target="MCM_LECA_DUPS", tier="T1", source_id=SRC, locator=f"{TXT_REL} Abstract (Results)", quote=Q["LECA"]),
        dict(assertion_id="E05_A2", type="event", target="MCM_LECA_DUPS", tier="T1", source_id=SRC, locator=f"{loc_res}, archaeal MCM paragraph", quote=Q["ARC"]),
        dict(assertion_id="E05_A3", type="event", target="MCM_ROOT", tier="T1", source_id=SRC, locator=f"{loc_res}, archaeal MCM paragraph", quote=Q["ARC1"]),
        dict(assertion_id="E05_A4", type="event", target="MCM_LECA_DUPS", tier="T1", source_id=SRC, locator=f"{TXT_REL} Discussion, Gain and loss of Mcm8 and Mcm9", quote=Q["M89LECA"]),
        dict(assertion_id="E05_A5", type="event", target="Mcm3_NAEGR_dup", tier="T1", source_id=SRC, locator=f"{loc_res}, first paragraph", quote=Q["NG"]),
        dict(assertion_id="E05_A6", type="event", target="Mcm3_NAEGR_dup", tier="T1", source_id=SRC, locator=f"{loc_res}, first paragraph", quote=Q["RECENT"]),
        dict(assertion_id="E05_A7", type="tag", target="paralogue labels Mcm2-Mcm9 (Additional files 1-4)", tier="T2", source_id=SRC, locator=f"{loc_res}, phylogeny paragraph", quote=Q["CLASS"]),
        dict(assertion_id="E05_A8", type="tag", target="CONTROL (Mcm2-7 no loss)", tier="T1", source_id=SRC, locator=f"{TXT_REL} Abstract (Results)", quote=Q["CTRL"]),
        dict(assertion_id="E05_A9", type="tag", target="LOSS (Mcm8/Mcm9)", tier="T1", source_id=SRC, locator=f"{TXT_REL} Discussion, Gain and loss of Mcm8 and Mcm9", quote=Q["LOSS"]),
        dict(assertion_id="E05_A10", type="tag", target="LOSS (Drosophila Mcm9)", tier="T1", source_id=SRC, locator=f"{TXT_REL} Discussion, Mcm8/9 function paragraph", quote=Q["DM"]),
    ]
    for k, c in enumerate(complete, 1):
        rows.append(dict(assertion_id=f"E05_C{k}", type="completeness", target=c, tier="T1", source_id=SRC,
                         locator=f"{TXT_REL} Methods, Identification of MCM homologues", quote=Q["COMPLETE"]))
    with open(FAM / "assertions.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader(); w.writerows(rows)

    files = [TXT_REL] + SUPP_REL + ["families/E05_MCM/cache/query_meta.tsv", "families/E05_MCM/cache/query.fasta",
                                    "families/E05_MCM/cache/fwd.tsv", "families/E05_MCM/cache/rev.tsv",
                                    "families/E05_MCM/cache/qc.tsv", "families/E05_MCM/cache/idmap.tsv"]
    meta = dict(
        id="E05_MCM", name="MCM2-9 replicative helicase family (Mcm2-7 + Mcm8/9 LECA paralogues; archaeal MCM outgroup)",
        panel="E", phenomena=["DEEP", "CONTROL", "LOSS", "LSE"],
        descriptors=dict(n_members=int(len(members)), n_species=int(members.species.nunique()),
                         mapped_fraction=round(float((inp.status == "mapped").mean()), 3),
                         status_counts={k: int(v) for k, v in rp.status.value_counts().items()},
                         per_species_mapped_fraction={k: round(float(v), 3) for k, v in frac.items()},
                         paralogues_per_species={c: sorted(members[members.species == c].label.unique().tolist()) for c in sorted(members.species.unique())}),
        sources=[dict(id=SRC, doi="10.1186/1471-2148-9-60", pmid="19292915", pmcid="PMC2667178",
                      title="Ancient diversification of eukaryotic MCM DNA replication proteins", year=2009,
                      source_method="phylogeny", table="Additional files 1-4 (PDF trees with accessions), Additional file 7 (Newick)",
                      locator="Abstract; Results; Discussion", quote=Q["LECA"],
                      files=[dict(path=p, sha256=sha(ROOT / p)) for p in files])],
        root_level="Archaea_Eukaryota", root_single_copy=False, complete_species=complete, decoys=[], distinct_from=[],
        tier_summary={f"{e}_{t}": int(v) for (e, t), v in sorted(cnt.items())}, protocol_version=1,
        completeness_quote=Q["COMPLETE"],
        notes=("Root S (T1) joins the single Halobacterium MCM (paper: one MCM per archaeon except M. jannaschii) with the "
               "eukaryotic Mcm2-9 set; MCM_LECA_DUPS is one multifurcating D node for the seven pre-LECA duplications "
               "(branching order unresolved in the paper: 'Both phylogenetic analyses failed to resolve the branching "
               "relationship'), lvl_lo=Eukaryota, lvl_hi=Asgard_Eukaryota because the sampled Euryarchaeota and TACK "
               "(Sulfolobus, Thermofilum) carry a single MCM; no Asgard archaea were sampled. Per paralogue: <=1 gene per "
               "species -> S polytomy (T2 rule b; no two panel-E species share a WGD). Naegleria has two Mcm3 (T1 recent "
               "duplication, D at NAEGR). Source errors handled: accessions printed under two paralogues (N. crassa "
               "XP_965670 Mcm4/Mcm7, T. brucei XP_951609 Mcm3/Mcm7, M. brevicollis Monbr1:32709 Mcm8/Mcm9) -> ambiguous; "
               "accessions that resolve to non-MCM proteins (S. cerevisiae NP_010883 = SPF1; E. histolytica XP_754207 = "
               "Aspergillus protein) -> excluded by homology QC. Table 1 lists Nematostella Mcm8/9 as absent, while "
               "Additional file 4, Additional file 7 (NvMcm8/NvMcm9) and the Methods animal Mcm8/9 list include them; "
               "the figure labels are used. Panel archaea other than HALSA (SACS2, THEKO, NITMS, LOKOS) and all bacteria "
               "were not sampled by the paper; their MCMs are not members (incompleteness, not absence)."),
    )
    (FAM / "family.yaml").write_text(yaml.safe_dump(meta, sort_keys=False, allow_unicode=True, width=120))


if __name__ == "__main__":
    main()
