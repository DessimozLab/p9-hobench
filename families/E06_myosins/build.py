#!/usr/bin/env python
"""E06_myosins: build members.tsv, mapping_report.tsv, reference_tree.nhx, assertions.tsv, family.yaml.

Source: Kollmar & Muehlhausen 2017 BMC Evol Biol 17:211 (PMID 28870165): figshare 4565155 motor-domain alignment
(class-encoded CyMoBase names <abbrev>Myo<class><variant>, Mhc = class 2, Myo_<X> = orphan) + Additional file 2
(species abbreviations). Mapping: protocol v1.1 domain rule (motor domain only) from cache/ (prep_mapping.py).
Run from project root: python families/E06_myosins/build.py
"""
import csv, hashlib, sys
from collections import defaultdict
from pathlib import Path
import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[2]
FAM = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(FAM))
from prep_mapping import ABMAP, load_seqs, parse  # noqa: E402
from hobench.newick import Node, write_newick  # noqa: E402

SRC = "Kollmar2017"
ALN_REL = "literature/papers/PMID28870165/raw/figshare_4565155/Myosin_motordomain_alignment.fas"
XLS_REL = "literature/papers/PMID28870165/raw/supp/12862_2017_1056_MOESM2_ESM.xls"
TXT_REL = "literature/papers/PMID28870165/fulltext.txt"
Q_LECA1 = ("The intron pattern analysis and the taxonomic distribution of the classes suggest two myosins in the last "
           "eukaryotic common ancestor, a class-1 prototype and another myosin, which is most likely the ancestor of all "
           "other myosin classes.")
Q_LECA2 = ("The phylogenetic analyses and gene structure comparisons strongly suggest only two myosins in the LECA, a "
           "class-1 myosin and a second myosin, which is not present in extant species anymore.")
Q_HGT = ("The sparse distribution of class-2 and class-4 myosins outside their major lineages contradicts their presence "
         "in the last eukaryotic common ancestor but instead strongly suggests early eukaryote-eukaryote horizontal gene transfer.")
Q_COMPLETE = "In addition, we performed TBLASTN searches in the respective genomes with myosins from different classes."
Q_CLASS = "Thus, we date the origin of class-22 back to the last common ancestor of the Amorphea."
ID_MIN, QCOV_MIN, QLEN_MIN, MARGIN = 98.0, 0.90, 250, 2.0
HGT_CLASSES = {"2", "4"}
PANEL_EUK = ["HUMAN", "DROME", "NEMVE", "AMPQE", "MNELE", "TRIAD", "MONBE", "CAPO3", "YEAST", "NEUCR", "BATDJ", "ENCCU",
             "DICDI", "ENTH1", "ARATH", "CHLRE", "CYAME", "THAPS", "PHYRM", "PLAF7", "TRYB2", "NAEGR", "GIAIC", "THETR",
             "ICHMU"]


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main():
    seqs = load_seqs()
    fw = pd.read_csv(FAM / "cache/fwd.tsv", sep="\t")
    rv = pd.read_csv(FAM / "cache/rev.tsv", sep="\t")
    rvb = rv.sort_values(["qseqid", "bitscore", "pident"], ascending=[True, False, False]).groupby("qseqid")
    rv_best = {t: (g.iloc[0].sseqid, (g.bitscore == g.iloc[0].bitscore).sum()) for t, g in rvb}
    fwg = {k: g.sort_values(["bitscore", "pident"], ascending=False) for k, g in fw.groupby("qseqid")}
    rep = []
    for name in sorted(seqs):
        p = parse(name)
        if not p or p[0] not in ABMAP:
            continue
        ab, cls, var = p
        code = ABMAP[ab]
        label = "orphan" if cls == "orphan" else f"Myo{cls}"
        loc = f"Myosin_motordomain_alignment.fas >{name} (abbreviation {ab}; class {cls})"
        base = dict(original_id=name, species=code, label=label, locator=loc)
        if not seqs[name]:
            rep.append({**base, "gene_id": "", "status": "unmapped", "reason": "empty sequence in alignment"}); continue
        g = fwg.get(name)
        if g is None or not len(g):
            rep.append({**base, "gene_id": "", "status": "unmapped", "reason": "no DIAMOND hit in panel proteome"}); continue
        b = g.iloc[0]
        others = g[g.sseqid != b.sseqid]
        margin = float(b.pident) - (float(others.iloc[0].pident) if len(others) else 0.0)
        qcov = min(b.length / b.qlen, 1.0)
        rb = rv_best.get(b.sseqid)
        rbh = bool(rb and rb[0] == name and rb[1] == 1)
        info = f"best {b.sseqid} {b.pident:.1f}% qcov {qcov:.2f} qlen {b.qlen} tcov {b.length / b.slen:.2f} margin {margin:.1f} rbh {rbh}"
        if b.pident >= ID_MIN and qcov >= QCOV_MIN and b.qlen >= QLEN_MIN and margin >= MARGIN and rbh:
            rep.append({**base, "gene_id": b.sseqid, "status": "mapped", "reason": "motor-domain mapping (protocol v1.1); " + info,
                        "_m": dict(mapping_method="sequence_domain_diamond", mapping_identity=f"{b.pident:.1f}",
                                   mapping_coverage="", mapping_qcov=f"{qcov:.3f}", mapping_qlen=str(int(b.qlen)),
                                   mapping_coverage_query=f"{qcov:.3f}", mapping_margin=f"{margin:.1f}", rbh="1")})
        elif b.pident >= ID_MIN and qcov >= QCOV_MIN:
            why = ("identical/near-identical motor domains (protocol 0.6)" if b.pident >= 99.5 and margin < MARGIN
                   else "fails margin/RBH/length")
            rep.append({**base, "gene_id": "", "status": "ambiguous", "reason": f"{why}; {info}"})
        else:
            rep.append({**base, "gene_id": "", "status": "unmapped", "reason": f"below thresholds; {info}"})
    byg = defaultdict(list)
    for i, r in enumerate(rep):
        if r["status"] == "mapped":
            byg[r["gene_id"]].append(i)
    for gid, ii in byg.items():
        if len(ii) > 1:
            for i in ii:
                rep[i].update(status="ambiguous", reason=rep[i]["reason"] + f"; {gid} also best hit of another source sequence")
                rep[i].pop("_m", None)
    mem = [dict(gene_id=r["gene_id"], species=r["species"], original_id_in_paper=r["original_id"], label=r["label"],
                **r["_m"], source=SRC, locator=r["locator"]) for r in rep if r["status"] == "mapped"]
    members = pd.DataFrame(mem).sort_values(["label", "species", "gene_id"]).reset_index(drop=True)
    members.to_csv(FAM / "members.tsv", sep="\t", index=False, quoting=csv.QUOTE_MINIMAL)
    rp = pd.DataFrame([{k: r.get(k, "") for k in ("original_id", "species", "label", "gene_id", "status", "reason", "locator")} for r in rep])
    rp.to_csv(FAM / "mapping_report.tsv", sep="\t", index=False)

    # ------------------------------------------------------------- tree
    cnt = defaultdict(int)

    def mk(nid, ev, tier="T2", **f):
        cnt[(ev, tier)] += 1
        return Node(nid, features={"id": nid, "ev": ev, "tier": tier, "src": SRC, **f})

    def terminal(cls):
        sub = members[members.label == cls]
        if len(sub) == 1:
            return Node(sub.gene_id.iloc[0])
        multi = sub.species.value_counts().max() > 1
        nd = mk(f"{cls}_class", "U" if (multi or cls[3:] in HGT_CLASSES) else "S")
        for gid in sorted(sub.gene_id):
            nd.add_child(Node(gid))
        return nd

    classes = sorted({c for c in members.label if c not in ("orphan",)}, key=lambda c: int(c[3:]))
    other = [terminal(c) for c in classes if c != "Myo1"]
    kids = []
    if "Myo1" in classes:
        kids.append(terminal("Myo1"))
    if len(other) > 1:
        o = mk("non_class1_myosins", "U")
        for c in other:
            o.add_child(c)
        kids.append(o)
    else:
        kids += other
    leca = mk("MYO_LECA", "D", tier="T1", lvl_lo="Eukaryota", lvl_hi="ABOVE_RANGE")
    for k in kids:
        leca.add_child(k)
    orph = sorted(members[members.label == "orphan"].gene_id)
    if orph:
        root = mk("MYO_root", "U")
        root.add_child(leca)
        for g in orph:
            root.add_child(Node(g))
    else:
        root = leca
    (FAM / "reference_tree.nhx").write_text(write_newick(root, features=["id", "ev", "lvl", "lvl_lo", "lvl_hi", "tier", "src", "tag"], lengths=False) + "\n")

    # ------------------------------------------------------------- assertions
    inp = rp[rp.status.isin(["mapped", "unmapped", "ambiguous"])]
    frac = inp.groupby("species").status.apply(lambda s: (s == "mapped").mean())
    complete = sorted(frac[frac >= 0.9].index)
    rows = [dict(assertion_id="E06_A1", type="event", target="MYO_LECA", tier="T1", source_id=SRC,
                 locator=f"{TXT_REL} Abstract, Results", quote=Q_LECA1),
            dict(assertion_id="E06_A2", type="event", target="MYO_LECA", tier="T1", source_id=SRC,
                 locator=f"{TXT_REL} Conclusions", quote=Q_LECA2),
            dict(assertion_id="E06_A3", type="tag", target="Myo2_class,Myo4_class (encoded U: HGT risk)", tier="T1",
                 source_id=SRC, locator=f"{TXT_REL} Abstract", quote=Q_HGT),
            dict(assertion_id="E06_A4", type="tag", target="lineage-specific classes (example; not encoded as D)", tier="T1",
                 source_id=SRC, locator=f"{TXT_REL} Results (class-22)", quote=Q_CLASS)]
    for k, c in enumerate(complete, 1):
        rows.append(dict(assertion_id=f"E06_C{k}", type="completeness", target=c, tier="T1", source_id=SRC,
                         locator=f"{TXT_REL} Methods", quote=Q_COMPLETE))
    with open(FAM / "assertions.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader(); w.writerows(rows)

    absent = [c for c in PANEL_EUK if c not in set(ABMAP.values())]
    meta = dict(
        id="E06_myosins", name="Myosin superfamily (class-1 vs all other classes; LECA)", panel="E",
        collection="E_cytoskeleton_motors", phenomena=["DEEP", "LSE", "LOSS", "PROM", "HGT"],
        descriptors=dict(n_members=int(len(members)), n_species=int(members.species.nunique()),
                         mapped_fraction=round(float((inp.status == "mapped").mean()), 3),
                         status_counts={k: int(v) for k, v in rp.status.value_counts().items()},
                         per_species_mapped_fraction={k: round(float(v), 3) for k, v in frac.items()},
                         classes_present=sorted(set(members.label), key=lambda c: (c == "orphan", int(c[3:]) if c[3:].isdigit() else 0))),
        sources=[dict(id=SRC, doi="10.1186/s12862-017-1056-2", pmid="28870165",
                      title="Myosin repertoire expansion coincides with eukaryotic diversification in the Mesoproterozoic era",
                      year=2017, source_method="phylogeny",
                      table="figshare 10.6084/m9.figshare.4565155 Myosin_motordomain_alignment.fas; Additional file 2",
                      locator="Abstract; Results; Conclusions", quote=Q_LECA1,
                      files=[dict(path=p, sha256=sha(ROOT / p)) for p in (ALN_REL, XLS_REL, TXT_REL)])],
        root_level="Eukaryota", root_single_copy=False, complete_species=complete, decoys=[], distinct_from=[],
        tier_summary={f"{e}_{t}": int(v) for (e, t), v in sorted(cnt.items())}, protocol_version=1.1,
        completeness_quote=Q_COMPLETE,
        notes=("T1 D node MYO_LECA (lvl_lo=Eukaryota, lvl_hi=ABOVE_RANGE): the paper infers only two LECA myosins, a "
               "class-1 prototype and a second prototype ancestral to all other classes; so class-1 vs every other class "
               "is a pre-LECA split. All other classes are post-LECA, mostly lineage-specific inventions whose parent "
               "class is not resolved, so the non-class-1 node is U (no class-vs-class level asserted; e.g. class-22 "
               "dated to Amorphea but its sister class is not given). Class-2 and class-4: HGT suggested by the paper -> "
               "U even when single-copy. Other classes: <=1 gene per species -> S polytomy (T2 rule b; no two panel-E "
               "species share a WGD), else U (rule c). Orphans (Myo_<X>) hang from the U root. Mapping: protocol v1.1 "
               "domain rule (motor domains ~690 aa; >=98% id over >=90% of query, RBH, >=2-point margin). Kinesin-myosin "
               "homology (shared motor fold): no distinct_from with E01. Panel-E eukaryotes absent from the source: "
               + ", ".join(absent) + ". Strains: YEAST=Sc_c, BATDJ=Bad_b, PLAF7=Pf_a. Sebe-Pedros 2014 (GBE 6:290) "
               "was read but not used: its DataS1 headers carry no class labels."),
    )
    (FAM / "family.yaml").write_text(yaml.safe_dump(meta, sort_keys=False, allow_unicode=True, width=120))


if __name__ == "__main__":
    main()
