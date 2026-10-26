#!/usr/bin/env python
"""E02_ARF: build members.tsv, mapping_report.tsv, reference_tree.nhx, assertions.tsv, family.yaml.

Source: Vargova R, Wideman JG, Derelle R, Klimes V, Kahn RA, Dacks JB, Elias M (2021) A eukaryote-wide perspective on
the diversity and evolution of the ARF GTPase protein family. Genome Biol Evol 13:evab157 (PMID 34247240).
Supplementary table 1 (2,286 ARF family genes, 115 species blocks, orthologous-group 'Annotation', 'Conclusively
annotated' flag, accessions) and Supplementary dataset 1 (protein FASTA). Network/DIAMOND results are cached in cache/
by prep_mapping.py; this script is offline.  Run from the project root:  python families/E02_ARF/build.py
"""
import csv, hashlib, re, sys
from collections import defaultdict
from pathlib import Path
import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[2]
FAM = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(FAM))
from prep_mapping import load_table, load_dataset, row_seq, accessions, ZIP_REL  # noqa: E402
from seqmap import best_hits, ID_MIN, COV_MIN  # noqa: E402
from hobench.newick import Node, write_newick  # noqa: E402
from hobench.speciestree import SpeciesTree  # noqa: E402

FID = "E02_ARF"
SRC = "Vargova2021"
PANEL = "E"
TXT_REL = "literature/papers/PMID34247240/fulltext.txt"
TXT = f"{TXT_REL}"
# --------------------------------------------------------------------------- quotes (verbatim, fulltext.txt)
Q = dict(
    LECA=("The analyses presented above indicate that the LECA possessed at least 15 ARF family genes; Arf1 and 6, "
          "Arl1, 2, 3, 5, 6, 8, 13, 16, 17, and 18, Arfrp1, Sar1, and SarB."),
    SRB="In addition, it certainly encoded SRβ, excluded from our ScrollSaw analysis",
    SARB=("multiple lines of evidence suggest they constitute a separate ancient paralog of their own, which we call "
          "SarB"),
    ARF4=("A common origin of Arf1 and Arf4 groups was already reported (Li et al. 2004; Manolea et al. 2010), but our "
          "analysis placed this event before the divergence of ichthyosporeans to the common ancestor of Holozoa"),
    ARF4b="that is in the holozoan stem itself (Arf4)",
    VERT=("the vertebrate ARF family complement has been expanded by duplications of Arf1 and Arf4, yielding the "
          "well-known two groups of paralogs (Arf1, 2, 3 vs. Arf4 and 5)."),
    ARL4="Together with multiple duplications of Arl4, vertebrates are thus endowed with a battery of lineage-specific paralogs",
    TRIM23="its specific relationship to Arf1 is obvious from the virtually identical exon–intron structure",
    ARFB="Hence, we leave the origin of the ArfB group as unresolved.",
    HOLO=("Altogether we could recognize seven groups that apparently originated after the split of the holozoan "
          "lineage from their relatives (Holomycota)"),
    ANN=("Sequences annotated based on the phylogenetic analyses or hmmscan searches are marked by asterisk (*) in the "
         "column “Conclusively annotated” in the supplementary table 1"),
    LSE="although gene duplications specific to terminal organismal lineages are common in the ARF family",
    LOSS=("Arl1, Arl5, and Arfrp1, all associated with the Golgi apparatus, have been preserved in all main eukaryote "
          "lineages sampled, but have been lost from some more terminal branches."),
    COMPLETE="ARF family GTPases were searched in genome and/or transcriptome assemblies from 114 eukaryotic species",
    RAS=("the content of the superfamily is pigeonholed into five major families (Ras, Rho, Rab, Ran, Arf/Sar; "
         "ColicelLi 2004; Rojas et al. 2012)"),
)
LECA = ["Arf1", "Arf6", "Arl1", "Arl2", "Arl3", "Arl5", "Arl6", "Arl8", "Arl13", "Arl16", "Arl17", "Arl18", "Arfrp1",
        "Sar1", "SarB", "SRb"]


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def norm_label(a):
    a = a.strip()
    return {"SRB": "SRb"}.get(a, a)


def main():
    st = SpeciesTree.from_newick(ROOT / f"panels/{PANEL}/species_tree.nwk")
    holozoa = set(st.clade_species("Holozoa"))
    t = load_table()
    ds = load_dataset()
    p = t[t.species != ""].copy()
    idm = pd.read_csv(FAM / "cache/idmap.tsv", sep="\t", dtype=str, keep_default_na=False)
    idmap = defaultdict(set)
    for r in idm.itertuples(index=False):
        if r.uniprot_acc:
            idmap[r.original_id].add(r.uniprot_acc)
    accs = {}
    for code in p.species.unique():
        accs[code] = {l[1:].split()[0] for l in open(ROOT / f"panels/{PANEL}/proteomes/{code}.fa") if l.startswith(">")}
    fw = pd.read_csv(FAM / "cache/fw.tsv", sep="\t")
    rv = pd.read_csv(FAM / "cache/rv.tsv", sep="\t")
    hits = {}
    for code in p.species.unique():
        hits[code] = best_hits(fw[fw.species == code], rv[rv.species == code])

    # identical copies flagged by the authors ("identical to X", full length) -> exclude both
    ident = set()
    for r in p.to_dict("records"):
        m = re.match(r"identical to (\S+)\s*$", r["Notes"].strip())
        if m:
            ident |= {r["Abbreviation"].lower(), m.group(1).lower()}

    rep = []
    for r in p.to_dict("records"):
        code, ab = r["species"], r["Abbreviation"]
        label = norm_label(r["Annotation"])
        concl = r["Conclusively annotated"] == "*"
        loc = (f"Supplementary_table_1.xlsx (in {ZIP_REL}) row {r['xlsx_row']}: {ab}; Annotation={r['Annotation']}; "
               f"Conclusively annotated={'*' if concl else '-'}; basis={r['Annotation basis']}; acc={r['Accession']}")
        base = dict(original_id=ab, species=code, label=label if concl else f"({label})", locator=loc)
        if "pseudogene" in r["Notes"].lower():
            rep.append({**base, "gene_id": "", "status": "excluded", "reason": "putative pseudogene (paper note)"}); continue
        if ab.lower() in ident:
            rep.append({**base, "gene_id": "", "status": "excluded", "reason": "identical copies (paper note: " + (r["Notes"] or "referenced as identical by another row") + ")"}); continue
        if r["Number of genes"] not in ("", "1"):
            rep.append({**base, "gene_id": "", "status": "excluded", "reason": f"row represents {r['Number of genes']} identical genes"}); continue
        cands, how = set(), []
        for a in accessions(r):
            h = (idmap.get(a, set()) | idmap.get(a.split(".")[0], set())) & accs[code]
            if h:
                cands |= h; how.append(f"UniProt ID mapping of {a}")
        b = hits[code].get(ab)
        if b is None and not row_seq(r, ds):
            seqnote = "no sequence in dataset 1"
        elif b is None:
            seqnote = "no DIAMOND hit"
        else:
            seqnote = (f"seq best hit {b['t']} {b['pid']:.1f}% id, cov q{b['qcov']:.0%}/t{b['tcov']:.0%}, margin "
                       f"{b['margin']:.1f}, rbh {b['rbh']}")
        if len(cands) == 1:
            g = next(iter(cands))
            if b and b["t"] != g and b["pid"] >= ID_MIN and b["qcov"] >= COV_MIN:
                rep.append({**base, "gene_id": "", "status": "ambiguous", "reason": f"ID mapping -> {g} but {seqnote}"}); continue
            rep.append({**base, "gene_id": g, "status": "mapped", "reason": "; ".join(how) + "; " + seqnote,
                        "_m": dict(mapping_method="id_mapping",
                                   mapping_identity=f"{b['pid']:.1f}" if b and b["t"] == g else "",
                                   mapping_coverage=f"{min(b['qcov'], b['tcov']):.2f}" if b and b["t"] == g else "")})
        elif len(cands) > 1:
            if b and b["ok"] and b["t"] in cands:
                rep.append({**base, "gene_id": b["t"], "status": "mapped", "reason": "; ".join(how) + f" -> {len(cands)} entries; resolved by {seqnote}",
                            "_m": dict(mapping_method="id_mapping+sequence", mapping_identity=f"{b['pid']:.1f}",
                                       mapping_coverage=f"{min(b['qcov'], b['tcov']):.2f}")})
            else:
                rep.append({**base, "gene_id": "", "status": "ambiguous", "reason": "ID maps to several panel entries: " + ",".join(sorted(cands))})
        elif b and b["ok"]:
            rep.append({**base, "gene_id": b["t"], "status": "mapped", "reason": "same-species sequence mapping: " + seqnote,
                        "_m": dict(mapping_method="sequence", mapping_identity=f"{b['pid']:.1f}",
                                   mapping_coverage=f"{min(b['qcov'], b['tcov']):.2f}")})
        elif b and b["pid"] >= ID_MIN and b["qcov"] >= COV_MIN and b["tcov"] >= COV_MIN:
            rep.append({**base, "gene_id": "", "status": "ambiguous", "reason": "no ID route; " + seqnote + " (fails RBH/margin)"})
        else:
            rep.append({**base, "gene_id": "", "status": "unmapped", "reason": "no ID route; " + seqnote})
    # same UniProt entry from >1 rows
    byg = defaultdict(list)
    for i, r in enumerate(rep):
        if r["status"] == "mapped":
            byg[r["gene_id"]].append(i)
    for g, ii in byg.items():
        if len(ii) > 1:
            for i in ii:
                rep[i].update(status="ambiguous", reason=rep[i]["reason"] + f"; {g} also mapped by another paper row")
                rep[i].pop("_m", None)
    mem = [dict(gene_id=r["gene_id"], species=r["species"], original_id_in_paper=r["original_id"], label=r["label"],
                **r["_m"], source=SRC, locator=r["locator"]) for r in rep if r["status"] == "mapped"]
    members = pd.DataFrame(mem).sort_values(["label", "species", "gene_id"]).reset_index(drop=True)
    members.to_csv(FAM / "members.tsv", sep="\t", index=False, quoting=csv.QUOTE_MINIMAL)
    pd.DataFrame([{k: r.get(k, "") for k in ("original_id", "species", "label", "gene_id", "status", "reason", "locator")}
                  for r in rep]).to_csv(FAM / "mapping_report.tsv", sep="\t", index=False)

    # ------------------------------------------------------------------- tree
    ids = defaultdict(int)
    asserted = []  # (node_id, quote key)

    def nid(base):
        base = re.sub(r"[^A-Za-z0-9_]", "_", base)
        ids[base] += 1
        return base if ids[base] == 1 else f"{base}_{ids[base]}"

    def mk(name, ev, kids, tier="T2", **feat):
        kids = [k for k in kids if k is not None]
        if not kids:
            return None
        if len(kids) == 1:
            return kids[0]
        i = nid(name)
        nd = Node(i, features={"ev": ev, "id": i, "tier": tier, "src": SRC, **feat})
        for k in kids:
            nd.add_child(k)
        return nd

    lab_genes = defaultdict(list)  # conclusive label -> [(gene, species, orig)]
    loose = []
    for m in members.itertuples(index=False):
        if m.label.startswith("(") or m.label == "unannotated":
            loose.append(m.gene_id)
        else:
            lab_genes[m.label].append((m.gene_id, m.species, m.original_id_in_paper))

    # quoted species-specific duplications (HUMAN is the only vertebrate in panel E)
    SPECIES_D = {("Arf1", "HUMAN"): "VERT", ("Arf4", "HUMAN"): "VERT", ("Arl4", "HUMAN"): "ARL4"}

    def label_group(label, genes, name=None):
        """T2 rule b/c: <=1 gene per species -> S polytomy; species with >1 gene -> U, unless a quoted species-level
        duplication (SPECIES_D) resolves them (D at that species)."""
        name = name or label
        by = defaultdict(list)
        for g, s, _ in genes:
            by[s].append(g)
        s_kids, u_kids = [], []
        for s in sorted(by):
            gs = sorted(by[s])
            if len(gs) == 1:
                s_kids.append(Node(gs[0]))
            elif (label, s) in SPECIES_D:
                d = mk(f"{name}_{s}_dup", "D", [Node(g) for g in gs], tier="T1")  # level = its mapping (leaf branch {s})
                asserted.append((d.name, SPECIES_D[(label, s)]))
                s_kids.append(d)
            else:
                u_kids += [Node(g) for g in gs]
        s_node = mk(f"{name}_S", "S", s_kids) if len(s_kids) > 1 else (s_kids[0] if s_kids else None)
        if u_kids:
            return mk(f"{name}_grp", "U", [s_node] + u_kids)
        return s_node

    def arf1_complex():
        a1 = lab_genes.get("Arf1", [])
        non = [x for x in a1 if x[1] not in holozoa]
        hol = [x for x in a1 if x[1] in holozoa]
        a4 = lab_genes.get("Arf4", [])
        hol_d = None
        if hol and a4:
            hol_d = mk("Arf1_Arf4_Holozoa", "D", [label_group("Arf1", hol, "Arf1_holozoa"), label_group("Arf4", a4)],
                       tier="T1", lvl="Holozoa")
            asserted.append((hol_d.name, "ARF4")); asserted.append((hol_d.name, "ARF4b"))
        else:
            hol_d = mk("Arf1_Arf4_holo", "U", [label_group("Arf1", hol, "Arf1_holozoa") if hol else None,
                                               label_group("Arf4", a4) if a4 else None])
        tr = label_group("TRIM23", lab_genes.get("TRIM23", [])) if lab_genes.get("TRIM23") else None
        return mk("Arf1_complex", "U", [label_group("Arf1", non, "Arf1_nonholozoa") if non else None, hol_d, tr])

    leca_kids = []
    for lab in LECA:
        if lab == "Arf1":
            leca_kids.append(arf1_complex())
        elif lab_genes.get(lab):
            leca_kids.append(label_group(lab, lab_genes[lab]))
    leca = mk("ARF_LECA", "D", leca_kids, tier="T1", lvl_lo="Eukaryota", lvl_hi="ABOVE_RANGE")
    asserted += [(leca.name, "LECA")] + [(leca.name, k) for k, lab in (("SRB", "SRb"), ("SARB", "SarB")) if lab_genes.get(lab)]
    used = set(LECA) | {"Arf4", "TRIM23"}
    others = [label_group(lab, lab_genes[lab]) for lab in sorted(lab_genes) if lab not in used]
    root = mk("ARF_root", "U", [leca] + others + [Node(g) for g in sorted(loose)])
    (FAM / "reference_tree.nhx").write_text(
        write_newick(root, features=["ev", "lvl", "lvl_lo", "lvl_hi", "id", "tier", "src", "tag"], lengths=False) + "\n")

    # ------------------------------------------------------------------- assertions
    rows = []
    for k, (node, qk) in enumerate(asserted, 1):
        rows.append(dict(assertion_id=f"{FID}_A{k}", type="event", target=node, tier="T1", source_id=SRC,
                         locator=f"{TXT_REL} (Results)", quote=Q[qk]))
    n = len(rows)
    extra = [("tag", "DEEP (LECA paralogs)", "LECA"), ("tag", "LSE", "LSE"), ("tag", "LOSS", "LOSS"),
             ("tag", "T2 rule: only 'Conclusively annotated' labels used for S/U groups; others hang from root U", "ANN"),
             ("tag", "ArfB group placed at root U (origin unresolved)", "ARFB"),
             ("tag", "Holozoa-specific groups (Arl4, 10, 15, 19, 20) placed at root U (source paralog not resolved)", "HOLO"),
             ("tag", "TRIM23 placed inside the Arf1 complex (U)", "TRIM23"),
             ("tag", "homology with Rab (E03): Arf/Sar and Rab are Ras-superfamily families; no distinct_from declared", "RAS")]
    for typ, tgt, qk in extra:
        n += 1
        rows.append(dict(assertion_id=f"{FID}_A{n}", type=typ, target=tgt, tier="T1" if qk in ("LECA", "LSE", "LOSS") else "T2",
                         source_id=SRC, locator=f"{TXT_REL}", quote=Q[qk]))
    rp = pd.read_csv(FAM / "mapping_report.tsv", sep="\t", dtype=str, keep_default_na=False)
    inp = rp[rp.status.isin(["mapped", "unmapped", "ambiguous"])]
    frac = inp.groupby("species").status.apply(lambda s: (s == "mapped").mean())
    complete = sorted(frac[frac >= 0.9].index)
    for k, c in enumerate(complete, 1):
        rows.append(dict(assertion_id=f"{FID}_C{k}", type="completeness", target=c, tier="T1", source_id=SRC,
                         locator=f"{TXT_REL} Materials and Methods", quote=Q["COMPLETE"]))
    with open(FAM / "assertions.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader(); w.writerows(rows)

    # ------------------------------------------------------------------- family.yaml
    cnt = defaultdict(int)

    def walk(nd):
        if nd.children:
            cnt[f"{nd.features.get('ev')}_{nd.features.get('tier')}"] += 1
            for c in nd.children:
                walk(c)
    walk(root)
    panel_euk = sorted(st.clade_species("Eukaryota"))
    present = sorted(p.species.unique())
    meta = dict(
        id=FID, name="ARF family GTPases (Arf, Arl, Sar, SRbeta, TRIM23; eukaryote-wide, LECA paralogs)", panel=PANEL,
        collection="E_small_GTPases", phenomena=["DEEP", "LOSS", "LSE", "PROM"],
        descriptors=dict(n_members=int(len(members)), n_species=int(members.species.nunique()),
                         mapped_fraction=round(float((inp.status == "mapped").mean()), 3),
                         n_rows_panel_species=int(len(rp)),
                         status_counts={k: int(v) for k, v in rp.status.value_counts().items()},
                         per_species_mapped_fraction={k: round(float(v), 3) for k, v in frac.items()},
                         labels={k: int(v) for k, v in members.label.value_counts().sort_index().items()}),
        sources=[dict(id=SRC, doi="10.1093/gbe/evab157", pmid="34247240",
                      title="A Eukaryote-Wide Perspective on the Diversity and Evolution of the ARF GTPase Protein Family",
                      year=2021, journal="Genome Biology and Evolution", source_method="phylogeny",
                      table="Supplementary table 1", locator="Results; Fig. 1, 3, 5; Supplementary table 1, dataset 1",
                      quote=Q["LECA"],
                      files=[dict(path=pth, sha256=sha(ROOT / pth)) for pth in (ZIP_REL, TXT_REL)])],
        root_level="Eukaryota", root_single_copy=False,
        complete_species=complete, decoys=[], distinct_from=[],
        homologous_families=[dict(family="E03_Rab", relation="Ras superfamily (homologous small GTPases); deep merges not penalised", src=SRC)],
        tier_summary=dict(sorted(cnt.items())),
        protocol_version=1,
        notes=("Journal is Genome Biol Evol (not MBE). Tree: root U; ARF_LECA = multifurcating T1 D over the 16 LECA "
               "paralogs (15 + SRbeta; interval [Eukaryota, ABOVE_RANGE]). Arf1 child = U 'Arf1_complex' holding "
               "non-holozoan Arf1, a T1 D at Holozoa (holozoan Arf1 vs Arf4, holozoan-stem origin of Arf4) and TRIM23 "
               "(derived from Arf1). Vertebrate Arf1/Arf4/Arl4 duplications -> D at the HUMAN leaf branch (only "
               "vertebrate in panel E). Labels: only rows marked 'Conclusively annotated' (*) are grouped (S if <=1 "
               "gene per species, multi-copy species under U; T2 rules b/c; no two panel-E species share a WGD, so "
               "rule-b pseudo-orthology exception does not apply); non-conclusive and 'unannotated' rows hang from the "
               "root U. Holozoa/Metazoa-specific groups (Arl4, 10, 15, 19, 20), ArfB, Arf1/6, ArlRhodo, Arl11/14, "
               "Arl13L stay at root U (source paralog not resolved in the paper). Excluded: pseudogenes, identical "
               "copies (paper notes / 'Number of genes' > 1). Mapping: RefSeq/GenBank -> UniProt ID mapping, else "
               "same-species DIAMOND of the paper's protein (dataset 1) with protocol thresholds. Panel-E eukaryotes "
               "absent from the source: " + ", ".join(c for c in panel_euk if c not in present) +
               " (source has Phytophthora sojae, Tetrahymena/Paramecium, no microsporidia). Homology: Arf/Sar and Rab "
               "are homologous Ras-superfamily small GTPases (no distinct_from with E03_Rab)."),
    )
    (FAM / "family.yaml").write_text(yaml.safe_dump(meta, sort_keys=False, allow_unicode=True, width=120))


if __name__ == "__main__":
    main()
