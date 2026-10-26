#!/usr/bin/env python
"""E04_tubulins: build members.tsv, mapping_report.tsv, reference_tree.nhx, assertions.tsv, family.yaml.

Source: Findeisen et al. 2014 GBE 6:2274 (PMID 25169981): tubulin_alignment.fas (3,286 named full-length tubulins;
name = <species abbreviation>Tub<class><copy>, class 1-6 = alpha..zeta) + DatasetS2.xls (species abbreviations).
Mapping: same-species DIAMOND (cache/fwd.tsv, cache/rev.tsv from prep_mapping.py) with protocol thresholds.
Run from project root: python families/E04_tubulins/build.py
"""
import csv, hashlib, re, sys
from collections import defaultdict
from pathlib import Path
import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[2]
FAM = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(FAM))
from prep_mapping import ABMAP, load_seqs, parse  # noqa: E402
from hobench.newick import Node, write_newick  # noqa: E402

SRC = "Findeisen2014"
ZIP_REL = "literature/papers/PMID25169981/raw/supp/supp_evu187_suppl_data.zip"
TXT_REL = "literature/papers/PMID25169981/fulltext.txt"
GREEK = {"1": "alpha", "2": "beta", "3": "gamma", "4": "delta", "5": "epsilon", "6": "zeta"}
Q_LECA = ("Members of all subfamilies are present in all major kingdoms of the eukaryotes implying that an ancestor of "
          "each subfamily must have been present in the last common ancestor of the eukaryotes")
Q_SIX = ("The consensus of the generated phylogenetic trees, based on different tree reconstruction methods and varying "
         "data sets, shows six major eukaryotic tubulin subfamilies.")
Q_NUM = "Thus, Tub1 is used equivalent to α-tubulin or Tub3 equivalent to γ-tubulin."
Q_COMPLETE = ("Tubulin genes have been identified in iterated TBLASTN searches of the completed or almost completed "
              "genomes of 504 species")
Q_FTSZ = ("The proposed archaeal homologs grouped together with the bacterial TubZ proteins as sister clade to the FtsZ "
          "proteins indicating that tubulins are unique to eukaryotes.")
ID_MIN, COV_MIN, MARGIN = 98.0, 0.80, 2.0
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
    rv_best = {t: (g.iloc[0].sseqid, g.iloc[0].bitscore, (g.bitscore == g.iloc[0].bitscore).sum()) for t, g in rvb}
    fwg = {k: g.sort_values(["bitscore", "pident"], ascending=False) for k, g in fw.groupby("qseqid")}
    sp = pd.read_csv(ROOT / "panels/E/species.tsv", sep="\t", dtype=str).set_index("code")
    exact = {}  # species -> seq -> [acc]
    for code in set(ABMAP.values()):
        d, name_, buf = defaultdict(list), None, []
        for l in list(open(ROOT / sp.loc[code, "fasta_path"])) + [">END\n"]:
            if l.startswith(">"):
                if name_:
                    d["".join(buf).upper()].append(name_)
                name_, buf = l[1:].split()[0], []
            else:
                buf.append(l.strip())
        exact[code] = d
    src_groups = defaultdict(list)  # (species, seq) -> source names
    for name in sorted(seqs):
        p = parse(name)
        if p and p[0] in ABMAP:
            src_groups[(ABMAP[p[0]], seqs[name])].append(name)
    rep = []
    for name in sorted(seqs):
        p = parse(name)
        if not p or p[0] not in ABMAP:
            continue
        ab, cls, suf = p
        code = ABMAP[ab]
        label = f"Tub{cls}"
        loc = f"tubulin_alignment.fas >{name} (DatasetS2 abbreviation {ab}; class Tub{cls} = {GREEK[cls]})"
        base = dict(original_id=name, species=code, label=label, locator=loc)
        tg = exact[code].get(seqs[name], [])
        grp = src_groups[(code, seqs[name])]
        if len(tg) == 1:
            if grp[0] != name:
                rep.append({**base, "gene_id": "", "status": "excluded",
                            "reason": f"identical protein copy of {grp[0]} (one panel entry {tg[0]}; protocol 0.6)"}); continue
            rep.append({**base, "gene_id": tg[0], "status": "mapped",
                        "reason": "exact full-length sequence identity to a unique panel entry" + (f" (identical source copies: {','.join(grp[1:])})" if len(grp) > 1 else ""),
                        "_m": dict(mapping_method="exact_sequence", mapping_identity="100.0", mapping_coverage="1.000",
                                   mapping_margin="", rbh="")})
            continue
        if len(tg) > 1:
            rep.append({**base, "gene_id": "", "status": "excluded",
                        "reason": f"identical to {len(tg)} panel entries {','.join(sorted(tg))} (identical copies; protocol 0.6)"}); continue
        g = fwg.get(name)
        if g is None or not len(g):
            rep.append({**base, "gene_id": "", "status": "unmapped", "reason": "no DIAMOND hit in panel proteome"}); continue
        b = g.iloc[0]
        others = g[g.sseqid != b.sseqid]
        runner = float(others.iloc[0].pident) if len(others) else 0.0
        margin = float(b.pident) - runner
        qcov, tcov = b.length / b.qlen, b.length / b.slen
        rb = rv_best.get(b.sseqid)
        rbh = bool(rb and rb[0] == name and rb[2] == 1)
        info = f"best {b.sseqid} {b.pident:.1f}% qcov {qcov:.2f} tcov {tcov:.2f} margin {margin:.1f} rbh {rbh}"
        if b.pident >= ID_MIN and min(qcov, tcov) >= COV_MIN and margin >= MARGIN and rbh:
            rep.append({**base, "gene_id": b.sseqid, "status": "mapped", "reason": info,
                        "_m": dict(mapping_method="sequence_diamond", mapping_identity=f"{b.pident:.1f}",
                                   mapping_coverage=f"{min(qcov, tcov):.3f}", mapping_margin=f"{margin:.1f}", rbh="1")})
        elif b.pident >= ID_MIN and min(qcov, tcov) >= COV_MIN:
            why = "identical/near-identical copies in source or proteome (excluded per protocol 0.6)" if (b.pident >= 99.5 and margin < MARGIN) else "fails margin/RBH"
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
    root = Node("TUB_LECA", features={"id": "TUB_LECA", "ev": "D", "lvl_lo": "Eukaryota", "lvl_hi": "ABOVE_RANGE",
                                       "tier": "T1", "src": SRC})
    cnt = defaultdict(int); cnt[("D", "T1")] += 1
    for cls in sorted(members.label.unique()):
        sub = members[members.label == cls]
        if len(sub) == 1:
            root.add_child(Node(sub.gene_id.iloc[0])); continue
        multi = sub.species.value_counts().max() > 1
        ev = "U" if multi else "S"
        nd = Node(f"{cls}_{GREEK[cls[-1]]}", features={"id": f"{cls}_{GREEK[cls[-1]]}", "ev": ev, "tier": "T2", "src": SRC})
        cnt[(ev, "T2")] += 1
        if multi:  # rule c: whole class unresolved, but genes of single-copy species could still be grouped -> keep flat U
            pass
        for gid in sorted(sub.gene_id):
            nd.add_child(Node(gid))
        root.add_child(nd)
    (FAM / "reference_tree.nhx").write_text(write_newick(root, features=["id", "ev", "lvl", "lvl_lo", "lvl_hi", "tier", "src", "tag"], lengths=False) + "\n")

    # ------------------------------------------------------------- assertions
    inp = rp[rp.status.isin(["mapped", "unmapped", "ambiguous"])]
    frac = inp.groupby("species").status.apply(lambda s: (s == "mapped").mean())
    complete = sorted(frac[frac >= 0.9].index)
    rows = [dict(assertion_id="E04_A1", type="event", target="TUB_LECA", tier="T1", source_id=SRC,
                 locator=f"{TXT_REL} Discussion, para 1 (fig. 7)", quote=Q_LECA),
            dict(assertion_id="E04_A2", type="event", target="TUB_LECA", tier="T1", source_id=SRC,
                 locator=f"{TXT_REL} Discussion, para 1", quote=Q_SIX),
            dict(assertion_id="E04_A3", type="tag", target="class labels Tub1-Tub6", tier="T2", source_id=SRC,
                 locator=f"{TXT_REL} Materials and Methods, 'Important note' on naming", quote=Q_NUM)]
    for k, c in enumerate(complete, 1):
        rows.append(dict(assertion_id=f"E04_C{k}", type="completeness", target=c, tier="T1", source_id=SRC,
                         locator=f"{TXT_REL} Results, first paragraph", quote=Q_COMPLETE))
    with open(FAM / "assertions.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader(); w.writerows(rows)

    absent = [c for c in PANEL_EUK if c not in set(ABMAP.values())]
    meta = dict(
        id="E04_tubulins", name="Tubulin superfamily (alpha, beta, gamma, delta, epsilon, zeta; LECA paralogs)", panel="E",
        collection="E_cytoskeleton_motors", phenomena=["DEEP", "LOSS", "LSE"],
        descriptors=dict(n_members=int(len(members)), n_species=int(members.species.nunique()),
                         mapped_fraction=round(float((inp.status == "mapped").mean()), 3),
                         status_counts={k: int(v) for k, v in rp.status.value_counts().items()},
                         per_species_mapped_fraction={k: round(float(v), 3) for k, v in frac.items()},
                         classes_per_species={c: sorted(members[members.species == c].label.unique().tolist()) for c in sorted(members.species.unique())}),
        sources=[dict(id=SRC, doi="10.1093/gbe/evu187", pmid="25169981",
                      title="Six subgroups and extensive recent duplications characterize the evolution of the eukaryotic tubulin protein family",
                      year=2014, source_method="phylogeny", table="supplementary tubulin_alignment.fas, DatasetS2.xls",
                      locator="Discussion para 1; Methods naming note; fig. 7", quote=Q_LECA,
                      files=[dict(path=p, sha256=sha(ROOT / p)) for p in (ZIP_REL, TXT_REL)])],
        root_level="Eukaryota", root_single_copy=False, complete_species=complete, decoys=[], distinct_from=[],
        tier_summary={f"{e}_{t}": int(v) for (e, t), v in sorted(cnt.items())}, protocol_version=1,
        completeness_quote=Q_COMPLETE,
        notes=("Root D node TUB_LECA (T1; lvl_lo=Eukaryota, lvl_hi=ABOVE_RANGE) separates the six classes Tub1-6 "
               "(alpha..zeta): the paper infers an ancestor of each subfamily in the LCA of eukaryotes. Class numbers: "
               "the paper states numbers replace Greek letters (Tub1=alpha, Tub3=gamma quoted; Tub2=beta from 'beta-tubulins "
               "Tub2D and Tub2E'); Tub4/5/6 = delta/epsilon/zeta follows the Greek order (consistent with Table 1 counts "
               "and human Tub4=chr17 TUBD1, Tub5=chr6 TUBE1 in Suppl. Table S3); only the class split matters for the tree. "
               "Within a class: <=1 mapped gene per species -> S polytomy (T2 rule b; no two panel-E species share a WGD); "
               ">1 in any species -> U (rule c; alpha/beta). FtsZ/TubZ/artubulin and Prosthecobacter BtubA/B are outgroups "
               "in the paper and not members; the panel has FtsZ in bacteria/archaea (paper: '" + Q_FTSZ + "'), but no "
               "FtsZ family is curated, so no distinct_from is declared. Panel-E eukaryotes absent from the source: "
               + ", ".join(absent) + ". Strains: YEAST=S288c (Sc_c), BATDJ=JAM81 (Bad_b), PLAF7=3D7 (Pf_a); GIAIC: source "
               "has Giardia ATCC 50581 (assemblage B/GS) vs panel WB (assemblage A) -> mostly unmapped (<98%). Sequences "
               "were annotated by the authors from genome assemblies, so some differ from UniProt gene models (<98%)."),
    )
    (FAM / "family.yaml").write_text(yaml.safe_dump(meta, sort_keys=False, allow_unicode=True, width=120))


if __name__ == "__main__":
    main()
