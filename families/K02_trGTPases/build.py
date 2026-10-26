#!/usr/bin/env python
"""Build K02_trGTPases (panel K): translation elongation factors EF-Tu/EF-1alpha vs EF-G/EF-2, the LUCA paralog pair.

Sources
  Baldauf1996  Baldauf SL, Palmer JD, Doolittle WF (1996) PNAS 93:7749 (PMID 8755547): abstract, T1 event
               (EF-Tu and EF-G are products of an ancient duplication predating all extant lineages).
  Iwabe1989    Iwabe N et al. (1989) PNAS 86:9355 (PMID 2531898): abstract, T1 event (same duplication).
  Creti1994    Creti R et al. (1994) PNAS 91:3255 (PMID 8159735): abstract; EF-1alpha(Tu) and EF-2(G) are the
               cross-domain (global-phylogeny) sets, i.e. archaeal EF-1alpha = EF-Tu group, archaeal EF-2 = EF-G group.
  COG2020      NCBI COG 2020: membership of COG0050 (EF-Tu, TufA), COG5256 (EF-1alpha, TEF1), COG0480 (EF-G, FusA)
               in the panel-K genomes (see literature/papers/COG2020/{extract_panelK,seqmap_panelK}.py).
Rules: same mapping rules as K01 (exact sha256 or unique DIAMOND >=98 % id, >=80 % cov both, 2-point margin; class 3
excluded; two COG proteins on one UniProt entry -> one member). Tree:
  root D (T1; no explicit level: defaults to [LCA mapping = panel root, ABOVE_RANGE] since root_single_copy=false)
    EF-Tu group: single-copy species -> S polytomy (T2 rule b); species with >1 member -> genes attached to a U
                 node over the group (T2 rule c, no source resolves the copies)
    EF-G group : idem.
Writes only into this directory.
"""
from __future__ import annotations

import csv
import glob
import hashlib
from collections import defaultdict
from pathlib import Path

import yaml

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
PANEL = "K"
CP = "literature/papers/COG2020"
COGCSV = f"{CP}/derived/cog20_panelK.csv"
COGFA = f"{CP}/derived/cog20_panelK.fa"
PANEL_SPECIES = ["ECOLI", "PSEAE", "BRADU", "HELPY", "GEOSL", "BACSU", "MYCTU", "CHLTR", "SYNY3", "DEIRA", "THEMA",
                 "AQUAE", "BACTN", "LEPIN", "SACS2", "KORCO", "NITMS", "THEKO", "METJA", "METAC", "HALSA"]
IN_COG = set(PANEL_SPECIES) - {"BRADU", "GEOSL"}
GROUPS = {"EF-Tu": {"COG0050": "EF-Tu", "COG5256": "EF-Tu"},   # COG5256 = archaeal EF-1alpha(Tu); COG id in locator
          "EF-G": {"COG0480": "EF-G"}}                         # COG0480 = EF-G / archaeal EF-2(G)
Q = {
    "baldauf": ("The genes for the protein synthesis elongation factors Tu (EF-Tu) and G (EF-G) are the products of an "
                "ancient gene duplication, which appears to predate the divergence of all extant organismal lineages."),
    "iwabe": ("we compared a pair of duplicated genes, elongation factors Tu and G, and the alpha and beta subunits of "
              "ATPase, which are thought to have diverged by gene duplication before divergence of the primary kingdoms."),
    "creti": ("Global phylogenies were inferred by alternative tree-making methods from available EF-2(G) sequence data "
              "and contrasted with phylogenies constructed from the more conserved but shorter EF-1 alpha(Tu) sequences."),
    "cog": "Comma-delimited plain text file assigning proteins to COGs",
}


def sha256(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def read_fasta(p) -> dict[str, str]:
    d, k = {}, None
    for ln in open(p):
        if ln.startswith(">"):
            k = ln[1:].split()[0]
            d[k] = []
        elif k:
            d[k].append(ln.strip())
    return {k: "".join(v).upper().replace("*", "") for k, v in d.items()}


def panel_fasta(code: str) -> Path:
    p = ROOT / f"panels/{PANEL}/proteomes/{code}.fa"
    if p.exists():
        return p
    for r in csv.DictReader(open(ROOT / "literature/raw/C/panels_verified.tsv"), delimiter="\t"):
        if r["panel"] == PANEL and r["mnemonic"] == code:
            return Path(glob.glob(str(ROOT / f"data/proteomes/{r['upid']}_*/proteome.fasta"))[0])
    raise SystemExit(f"no proteome for {code}")


def main() -> None:
    cog2grp = {c: g for g, d in GROUPS.items() for c in d}
    cogseq = read_fasta(ROOT / COGFA)
    rows = [r for r in csv.reader(open(ROOT / COGCSV)) if r[6] in cog2grp]
    seqmap = {}
    for code in IN_COG:
        p = ROOT / f"{CP}/derived/seqmap/{code}.tsv"
        for r in csv.DictReader(open(p), delimiter="\t"):
            seqmap[(code, r["query_id"])] = r
    fcache: dict[str, dict] = {}
    report, by_target = [], defaultdict(list)
    for r in rows:
        code, pid, cog, cls = r[-1], r[2], r[6], r[8]
        grp, label = cog2grp[cog], GROUPS[cog2grp[cog]][cog]
        loc = f"cog20_panelK.csv: {r[0]},{r[1]},{pid},{cog},class={cls}"
        rep = dict(original_id=pid, species=code, label=label, gene_id="", status="", reason="", locator=loc)
        report.append(rep)
        if cls == "3":
            rep.update(status="excluded", reason="COG membership class 3 (partial match on protein and COG profile)")
            continue
        if code not in fcache:
            fcache[code] = read_fasta(panel_fasta(code))
        exact = [k for k, s in fcache[code].items() if s == cogseq[pid]]
        sm = seqmap.get((code, pid))
        if len(exact) == 1:
            g, meth, ident, cov = exact[0], "exact_sequence", 100.0, 1.0
        elif len(exact) > 1:
            rep.update(status="ambiguous", reason=f"identical sequence in {len(exact)} panel entries")
            continue
        elif sm and sm["status"] == "unique" and float(sm["mapping_identity"]) >= 98 and \
                min(float(sm["query_cov"]), float(sm["target_cov"])) >= 0.8:
            g, meth = sm["target_id"], "diamond_sequence"
            ident, cov = float(sm["mapping_identity"]), min(float(sm["query_cov"]), float(sm["target_cov"]))
        else:
            why = (f"best hit {sm['target_id']} status={sm['status']} identity={sm['mapping_identity']} "
                   f"qcov={sm['query_cov']} tcov={sm['target_cov']}") if sm else "no sequence hit"
            rep.update(status="unmapped", reason=f"below protocol thresholds: {why}")
            continue
        rep.update(status="mapped", gene_id=g, reason=meth)
        by_target[(code, g)].append((rep, dict(gene_id=g, species=code, original_id_in_paper=pid, label=label,
                                               group=grp, mapping_method=meth, mapping_identity=f"{ident:.1f}",
                                               mapping_coverage=f"{cov:.3f}", source="COG2020", locator=loc)))
    members = []
    for (code, g), lst in sorted(by_target.items()):
        lst.sort(key=lambda x: (x[1]["mapping_method"] != "exact_sequence", x[1]["original_id_in_paper"]))
        if len({x[1]["group"] for x in lst}) > 1:
            raise SystemExit(f"{g} assigned to both groups")
        members.append(lst[0][1])
        for rep, _ in lst[1:]:
            rep.update(status="excluded", reason=f"identical/collapsed copy: maps to the same UniProt entry {g} as "
                                                 f"{lst[0][1]['original_id_in_paper']} (protocol rule 6)")
    members.sort(key=lambda m: (m["group"], PANEL_SPECIES.index(m["species"]), m["gene_id"]))

    def w_tsv(name, cols, data):
        with open(HERE / name, "w", newline="") as fh:
            wr = csv.DictWriter(fh, cols, delimiter="\t", extrasaction="ignore", lineterminator="\n")
            wr.writeheader()
            wr.writerows(data)

    for code in sorted({"BRADU", "GEOSL"}):
        report.append(dict(original_id=f"{code}:EF-Tu/EF-G", species=code, label="EF-Tu (COG0050); EF-G (COG0480)",
                           gene_id="", status="unmapped", reason="genome not in COG 2020; no identifier in any source",
                           locator="cog-20.org.csv (absent)"))
    w_tsv("mapping_report.tsv", ["original_id", "species", "label", "gene_id", "status", "reason", "locator"],
          sorted(report, key=lambda r: (PANEL_SPECIES.index(r["species"]), r["label"], r["original_id"])))
    w_tsv("members.tsv", ["gene_id", "species", "original_id_in_paper", "label", "mapping_method", "mapping_identity",
                          "mapping_coverage", "source", "locator"], members)

    subtrees, multi = [], {}
    for grp, gid in (("EF-Tu", "TU"), ("EF-G", "G")):
        ms = [m for m in members if m["group"] == grp]
        cnt = defaultdict(int)
        for m in ms:
            cnt[m["species"]] += 1
        single = [m["gene_id"] for m in ms if cnt[m["species"]] == 1]
        extra = [m["gene_id"] for m in ms if cnt[m["species"]] > 1]
        multi[grp] = sorted({m["species"] for m in ms if cnt[m["species"]] > 1})
        s = f"({','.join(single)})[&&NHX:ev=S:id={gid}_orth:tier=T2:src=Creti1994]"
        if extra:
            s = f"({s},{','.join(extra)})[&&NHX:ev=U:id={gid}_grp:tier=T2:src=COG2020]"
        subtrees.append(s)
    tree = (f"({subtrees[0]},{subtrees[1]})"
            f"[&&NHX:ev=D:id=TU_G_dup:tag=DEEP:tier=T1:src=Baldauf1996];\n")
    (HERE / "reference_tree.nhx").write_text(tree)

    tu = [m["gene_id"] for m in members if m["group"] == "EF-Tu"][:3]
    gg = [m["gene_id"] for m in members if m["group"] == "EF-G"][:3]
    assertions = [
        dict(assertion_id="K02_A1", type="event", target="TU_G_dup", tier="T1", source_id="Baldauf1996",
             locator="fulltext.txt Abstract, sentence 1", quote=Q["baldauf"]),
        dict(assertion_id="K02_A2", type="event", target="TU_G_dup", tier="T1", source_id="Iwabe1989",
             locator="fulltext.txt Abstract", quote=Q["iwabe"]),
        dict(assertion_id="K02_A3", type="tag", target="TU_orth,G_orth", tier="T2", source_id="Creti1994",
             locator="abstract.txt (Europe PMC abstract), sentence 2", quote=Q["creti"]),
        dict(assertion_id="K02_A4", type="tag", target="TU_orth,G_orth", tier="T2", source_id="COG2020",
             locator="raw/Readme.2020-11-25.txt (cog-20.cog.csv); COG0050/COG5256/COG0480 in cog-20.def.tab",
             quote=Q["cog"]),
    ]
    w_tsv("assertions.tsv", ["assertion_id", "type", "target", "tier", "source_id", "locator", "quote"], assertions)

    def files(lst):
        return [dict(path=p, sha256=sha256(ROOT / p)) for p in lst]

    n_map = sum(r["status"] == "mapped" for r in report)
    n_den = sum(r["status"] != "excluded" for r in report)
    meta = dict(
        id="K02_trGTPases",
        name="Translation elongation factors EF-Tu/EF-1alpha vs EF-G/EF-2 (LUCA paralogs)",
        panel=PANEL,
        phenomena=["DEEP", "CONTROL"],
        descriptors=dict(n_members=len(members), n_species=len({m["species"] for m in members}),
                         n_duplications=1, oldest_event="above panel root (pre-LUCA)",
                         multicopy_species=multi, mapped_fraction=round(n_map / n_den, 3)),
        sources=[
            dict(id="Baldauf1996", doi="10.1073/pnas.93.15.7749", pmid="8755547",
                 title="The root of the universal tree and the origin of eukaryotes based on elongation factor phylogeny",
                 year=1996, source_method="phylogeny", files=files(["literature/papers/PMID8755547/raw/fulltext.txt"])),
            dict(id="Iwabe1989", doi="10.1073/pnas.86.23.9355", pmid="2531898",
                 title="Evolutionary relationship of archaebacteria, eubacteria, and eukaryotes inferred from "
                       "phylogenetic trees of duplicated genes", year=1989, source_method="phylogeny",
                 files=files(["literature/papers/PMID2531898/raw/fulltext.txt"])),
            dict(id="Creti1994", doi="10.1073/pnas.91.8.3255", pmid="8159735",
                 title="Evolution of translational elongation factor (EF) sequences: reliability of global phylogenies "
                       "inferred from EF-1 alpha(Tu) and EF-2(G) proteins", year=1994, source_method="phylogeny",
                 files=files(["literature/papers/PMID8159735/abstract.txt"])),
            dict(id="COG2020", doi="10.1093/nar/gkaa1018", pmid="33167031",
                 title="COG database update: focus on microbial diversity, model organisms, and widespread pathogens",
                 year=2021, source_method="expert_classification",
                 files=files([COGCSV, COGFA, f"{CP}/raw/Readme.2020-11-25.txt"]
                             + sorted(str(p.relative_to(ROOT)) for p in (ROOT / CP / "derived/seqmap").glob("*.tsv")))),
        ],
        root_single_copy=False,
        complete_species=[],
        decoys=[],
        distinct_from=[],
        tier_summary=dict(T1_events=1, T2_nodes=2, T3=0),
        protocol_version=1,
        notes=("Root D (T1): EF-Tu vs EF-G duplication predates LUCA -> the two groups are never in one HOG within the "
               "panel. Within each group, single-copy species form an S polytomy (T2; archaeal EF-1alpha/EF-2 are the "
               "cross-domain counterparts of EF-Tu/EF-G per Creti 1994). Species with several COG members (EF-Tu: "
               "ECOLI tufA/tufB, distinct UniProt entries; EF-G: EF-G paralogs in several bacteria) hang from a U "
               "node (no source resolves them). Identical tuf copies collapsed to one UniProt entry (PSEAE, DEIRA, "
               "AQUAE) count once. IF2, LepA, SelB, RF3, BipA/TypA, aIF2gamma not included (no source in hand "
               "asserts their placement). BRADU/GEOSL absent from COG 2020 -> no members."),
    )
    (HERE / "family.yaml").write_text(yaml.safe_dump(meta, sort_keys=False, allow_unicode=True, width=120))
    print(f"K02: {len(members)} members, mapped {n_map}/{n_den}, multicopy {multi}")


if __name__ == "__main__":
    main()
