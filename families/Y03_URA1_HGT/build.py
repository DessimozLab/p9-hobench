#!/usr/bin/env python
"""Build Y03_URA1_HGT (panel Y): dihydroorotate dehydrogenase (DHOD). The bacterial-type family-1A DHOD URA1 was
acquired by HGT from Lactobacillales into the Saccharomycetaceae; the native eukaryotic mitochondrial family-2 DHOD
(URA9-type; S. pombe ura3) is retained by the other yeasts and, in several Saccharomycetaceae, alongside URA1.

Sources
  Hall2005        Hall C, Brachat S, Dietrich FS (2005) Eukaryot Cell 4:1102 (PMID 15947202; PMC, not OA subset):
                  T1 HGT event, per-species DHOD type statements, GenPept accessions of resequenced DHODs.
  Gojkovic2004    Gojkovic Z et al. (2004) Mol Genet Genomics 271:387 (PMID 15014982; closed): abstract only.
  Bouwknegt2021   Bouwknegt J et al. (2021) Fungal Biol Biotechnol 8:10 (PMID 34656184, OA CC BY): Additional file 1
                  (Ura9 = Class-II DHOD orthologs, UniProt accessions; HMMER + proteinortho), Additional file 3
                  (RAxML-NG tree with bootstrap) -> family-2 members + T3 clades (BS >= 70).
  YGOBv7          YGOB v7 Pillars.tab (synteny-curated pillar of URA1, line 4856) -> family-1A members (T2).
Tree
  root H (T1, donor=OUTSIDE: Lactobacillales are not in panel Y; tag HGT)
    URA1 / family 1A (recipient=1): U[ S( S(KLE pair), S(ZT pair) ), post-WGD genes ] (T2, YGOB pillar; post-WGD
        genes left U: pillar slots are not track-assigned, YGOB v1 scores URA1 as non-syntenic '!', and the hybrid
        origin of the WGD lineage makes ZT-vs-post-WGD / KLE-vs-post-WGD merge levels tree-dependent)
    URA9 / family 2 (native): U over all members, except nested clades with bootstrap >= 70 in Bouwknegt 2021
        Additional file 3 restricted to panel members (species-overlap rule -> S, T3; post-WGD genes attached directly
        to one S node are grouped under U).
Writes only into this directory.
"""
from __future__ import annotations

import csv
import glob
import hashlib
import re
from collections import defaultdict
from pathlib import Path

import yaml

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
PANEL = "Y"
SC = HERE / "source_cache"
HALL = "literature/papers/PMID15947202/raw/fulltext.txt"
GOJ = "literature/papers/PMID15014982/abstract.txt"
BK = "literature/papers/PMID34656184"
BK_IDS = f"{BK}/raw/supp/40694_2021_117_MOESM1_ESM.txt"
BK_TREE = f"{BK}/raw/supp/40694_2021_117_MOESM3_ESM.txt"
YGOB = "literature/raw/B1/ygob/Pillars.TRUNCATED_wayback_20210728.tab"
YGOB_LINE = 4856
YGOBV1 = "literature/raw/B1/ygob/wayback/v1_Synteny_Klactis.tab"
YGOB_COLS = ["VANPO", "TETPH", "HENBL", "NAUDC", "NAUCC", "KAZNA", "KAZAF", "CANGA", "SACUV", "SAKUD", "SAMIK",
             "YEAST", "Anc", "ZYGRC", "TORDC", "KLULA", "ASHGO", "ERECY", "LAKLU", "LACTH", "LAWAL"]  # pos. 1 + pre-WGD
POST_WGD = {"YEAST", "CANGA", "NAUCC", "NAUDC", "KAZAF", "VANPO", "TETPH"}
KLE = {"KLULA", "LACTH", "ASHGO", "ERECY"}
ZT = {"ZYGRC", "TORDC"}
# Hall 2005 GenPept accessions (Results, 'URA1 gene phylogeny ...'); species as named in the paper -> panel code.
# 'K. lactis (GenPept AAR1773)' is a typo in the paper: no such accession; the K. lactis DHOD deposited in the same
# series is AAQ01773 (NCBI title 'dihydroorotate dehydrogenase [Kluyveromyces lactis]') -> used, flagged.
HALL_ACC = [("AAP39962", "C. albicans", "CANAL"), ("AAQ57200", "S. castellii", "NAUCC"),
            ("AAR17523", "C. glabrata", "CANGA"), ("AAQ01773", "K. lactis", "KLULA"),
            ("AAR17522", "K. marxianus", "KLUMA"), ("AAQ01777", "S. bayanus", "-"), ("AAQ01778", "S. mikatae", "-"),
            ("AAQ01779", "S. kluyveri (URA1)", "-"), ("AAQ57201", "S. kluyveri (URA9)", "-"),
            ("AAQ01780", "S. paradoxus", "-"), ("AAP04498", "E. faecalis", "-"), ("AAQ01774", "L. mesenteroides", "-"),
            ("AAQ01776", "L. lactis subsp. cremoris", "-"), ("AAQ01775", "L. lactis subsp. hordniae", "-")]
Q = {
    "hgt": ("The DHOD-encoding gene found in S. cerevisiae, URA1 (YKL216W), appears to have entered the "
            "Saccharomycetaceae after the divergence of the S. cerevisiae lineage from the Candida albicans lineage"),
    "donor": "URA1 gene phylogeny is consistent with horizontal transfer from lactic acid bacteria.",
    "donor2": ("This gene appears to have come from the Lactobacillales, and following its acquisition the endogenous "
               "eukaryotic DHOD gene was lost."),
    "types": ("Closely related Saccharomyces species including S. bayanus, S. paradoxus, and S. castellii have the type "
              "1a gene, whereas K. marxianus and C. glabrata have the type 2 enzyme and K. lactis, K. waltii, and S. "
              "kluyveri have both family 1- and family 2-type enzymes."),
    "both": ("The presence of two DHOD genes in K. lactis, K. waltii, and S. kluyveri supports the argument that URA1 is "
             "horizontally transferred and suggests that an ancestor of S. cerevisiae had both DHOD genes."),
    "loss": "It is likely that the family 2 DHOD was lost at some time after the divergence of C. glabrata and S. cerevisiae.",
    "ashgo": "the family 2 DHOD from A. gossypii can complement the uracil auxotrophy",
    "euk2": "Eukaryotes typically have the family 2 DHOD enzyme.",
    "goj1": ("While the majority of yeasts possess a mitochondrial DHODase, Saccharomyces cerevisiae has a cytoplasmatic "
             "enzyme, whose activity is independent of the presence of oxygen."),
    "goj2": "From the phylogenetic point of view, this enzyme is closely related to a bacterial DHODase from Lactococcus lactis.",
    "bk1": ("as query yielded 331 fungal and 73 bacterial Ura9 orthologs (Additional file 1). These Ura9 orthologs were "
            "used to calculate a maximum-likelihood phylogenetic tree"),
    "bk2": "Ura9-orthologous amino-acid sequences were then subjected to multiple sequence alignment",
    "ygob": ("TPHA0K00200 TBLA0B10140 --- NCAS0D03520 --- --- --- --- Suva_11.4 Skud_11.4 Smik_11.7 YKL216W --- "
             "ZYRO0D00396g TDEL0F05620 KLLA0C19360g --- --- SAKL0D15268g KLTH0F00308g Kwal_27.9662"),
    "ygobv1": "???!!? KLLA0C19360g undef undef undef YKL216W Scas_652.19 undef",
}
BS_MIN = 70


def sha256(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def panel_species() -> list[str]:
    p = ROOT / f"panels/{PANEL}/species.tsv"
    if p.exists():
        return [r["code"] for r in csv.DictReader(open(p), delimiter="\t")]
    return [r["mnemonic"] for r in csv.DictReader(open(ROOT / "literature/raw/C/panels_verified.tsv"), delimiter="\t")
            if r["panel"] == PANEL and r["tier"] == "core"]


def panel_fasta(code: str) -> Path:
    p = ROOT / f"panels/{PANEL}/proteomes/{code}.fa"
    if p.exists():
        return p
    for r in csv.DictReader(open(ROOT / "literature/raw/C/panels_verified.tsv"), delimiter="\t"):
        if r["panel"] == PANEL and r["mnemonic"] == code:
            return Path(glob.glob(str(ROOT / f"data/proteomes/{r['upid']}_*/proteome.fasta"))[0])
    raise SystemExit(f"no proteome for {code}")


def fasta_ids(p: Path) -> set[str]:
    return {ln[1:].split()[0] for ln in open(p) if ln.startswith(">")}


def parse_newick(s: str):
    """Minimal Newick parser -> nested (children list | leaf name, support)."""
    pos = 0

    def node():
        nonlocal pos
        if s[pos] == "(":
            pos += 1
            kids = [node()]
            while s[pos] == ",":
                pos += 1
                kids.append(node())
            assert s[pos] == ")"
            pos += 1
            m = re.match(r"[^:,();]*", s[pos:])
            lab = m.group(0)
            pos += len(lab)
            if s[pos] == ":":
                pos += len(re.match(r":[^,();]*", s[pos:]).group(0))
            return (kids, float(lab) if lab else None)
        m = re.match(r"[^:,();]+", s[pos:])
        name = m.group(0)
        pos += len(name)
        if s[pos] == ":":
            pos += len(re.match(r":[^,();]*", s[pos:]).group(0))
        return (name, None)

    return node()


def supported_clades(tree, keep: set[str]) -> dict[frozenset, float]:
    """Clades (restricted to `keep`) with bootstrap >= BS_MIN and >= 2 kept leaves, proper subsets of all kept
    -> highest support among the tree nodes that induce that kept set."""
    out: dict[frozenset, float] = {}

    def rec(n):
        lab, sup = n
        if isinstance(lab, str):
            return {lab} & keep
        s = set()
        for c in lab:
            s |= rec(c)
        if sup is not None and sup >= BS_MIN and len(s) >= 2:
            k = frozenset(s)
            out[k] = max(sup, out.get(k, 0))
        return s

    allk = rec(tree)
    return {c: v for c, v in out.items() if c != allk}


def main() -> None:
    species = panel_species()
    pset = set(species)
    ids = {c: fasta_ids(panel_fasta(c)) for c in species}
    acc2sp = {a: c for c, s in ids.items() for a in s}
    report, members = [], {}

    def add(gid, sp, orig, label, typ, meth, ident, cov, src, loc):
        if gid in members:  # same gene from a second source: keep first, note the corroboration
            members[gid]["locator"] += f" | also {src}: {orig}"
            return
        members[gid] = dict(gene_id=gid, species=sp, original_id_in_paper=orig, label=label, type=typ,
                            mapping_method=meth, mapping_identity=ident, mapping_coverage=cov, source=src,
                            locator=loc)

    # --- family 1A (URA1) from the YGOB pillar
    row = open(ROOT / YGOB).read().splitlines()[YGOB_LINE - 1].split("\t")
    loci = {r["query_locus"]: r for r in csv.DictReader(open(SC / "uniprot_ygob_ura1_loci.tsv"), delimiter="\t")}
    for col, code in enumerate(YGOB_COLS):
        g = row[col]
        if g == "---" or code == "Anc":
            continue
        loc = f"Pillars.TRUNCATED_wayback_20210728.tab:line {YGOB_LINE} col {col + 1}"
        rep = dict(original_id=g, species=code, label="URA1 (family 1A DHOD)", gene_id="", status="", reason="",
                   locator=loc, source="YGOBv7")
        report.append(rep)
        if code not in pset:
            rep.update(status="not_in_panel", reason="species not in panel Y")
            continue
        u = loci.get(g)
        if not u or u["Entry"] not in ids[code]:
            rep.update(status="unmapped", reason="no UniProt entry with this locus name in the panel proteome")
            continue
        rep.update(status="mapped", gene_id=u["Entry"], reason="gene_name_exact (UniProt locus/ORF name, 1 hit)")
        add(u["Entry"], code, g, "URA1 (family 1A DHOD)", "1A", "gene_name_exact", "100.0", "1.000", "YGOBv7", loc)

    # --- Hall 2005 GenPept accessions (same-species sequence mapping)
    for acc, sp_name, code in HALL_ACC:
        loc = "fulltext.txt Results 'URA1 gene phylogeny is consistent with horizontal transfer'"
        rep = dict(original_id=acc, species=code if code != "-" else sp_name, label="DHOD (GenPept, Hall 2005)",
                   gene_id="", status="", reason="", locator=loc, source="Hall2005")
        report.append(rep)
        if code == "-" or code not in pset:
            rep.update(status="not_in_panel", reason=f"{sp_name}: species not in panel Y")
            continue
        sm = list(csv.DictReader(open(SC / f"seqmap_{code}.tsv"), delimiter="\t"))[0]
        ok = sm["status"] in ("exact", "unique") and float(sm["mapping_identity"]) >= 98 and \
            min(float(sm["query_cov"]), float(sm["target_cov"])) >= 0.8
        if not ok:
            rep.update(status="unmapped", reason=f"below thresholds: {sm['status']} {sm['mapping_identity']}")
            continue
        g = sm["target_id"]
        # type: from the Hall 2005 per-species statement (Q["types"]): S. castellii 1a; C. glabrata, K. marxianus 2;
        # K. lactis has both -> resolved by identity with the YGOB URA1 pillar gene; C. albicans 2 (Bouwknegt 2021).
        typ = {"NAUCC": "1A", "CANGA": "2", "KLUMA": "2", "CANAL": "2"}.get(code)
        if code == "KLULA":
            typ = members[g]["type"] if g in members else None
        if typ is None or (g in members and members[g]["type"] != typ):
            raise SystemExit(f"{acc}: DHOD type not determined by a source")
        note = " (paper prints 'AAR1773'; typo, see build.py)" if acc == "AAQ01773" else ""
        rep.update(status="mapped", gene_id=g, reason=f"{sm['status']} sequence match{note}")
        lab = "URA1 (family 1A DHOD)" if typ == "1A" else "URA9-type (family 2 DHOD)"
        add(g, code, acc, lab, typ, "exact_sequence" if sm["status"] == "exact" else "diamond_sequence",
            f"{float(sm['mapping_identity']):.1f}", f"{min(float(sm['query_cov']), float(sm['target_cov'])):.3f}",
            "Hall2005", f"{loc}; GenPept {acc}{note}")

    # --- Bouwknegt 2021 Ura9 (family 2) orthologs: UniProt accessions in panel proteomes
    bk = [ln[1:].split()[0] for ln in open(ROOT / BK_IDS) if ln.startswith(">")]
    n_out = 0
    for acc in bk:
        code = acc2sp.get(acc)
        if not code:
            n_out += 1
            continue
        loc = "Additional file 1 (40694_2021_117_MOESM1_ESM.txt) FASTA header"
        report.append(dict(original_id=acc, species=code, label="Ura9 (Class-II DHOD)", gene_id=acc, status="mapped",
                           reason="exact UniProt accession in panel proteome", locator=loc, source="Bouwknegt2021"))
        if acc in members and members[acc]["type"] != "2":
            raise SystemExit(f"{acc}: type conflict")
        add(acc, code, acc, "URA9-type (family 2 DHOD)", "2", "uniprot_accession", "100.0", "1.000",
            "Bouwknegt2021", loc)
    report.append(dict(original_id=f"{n_out} accessions", species="-", label="Ura9 (Class-II DHOD)", gene_id="",
                       status="not_in_panel", reason="accession not in any panel-Y proteome (other fungi/bacteria)",
                       locator="Additional file 1", source="Bouwknegt2021"))

    mem = sorted(members.values(), key=lambda m: (m["type"], species.index(m["species"]), m["gene_id"]))
    for t in ("1A", "2"):
        sp = [m["species"] for m in mem if m["type"] == t]
        if len(sp) != len(set(sp)):
            raise SystemExit(f"type {t}: >1 gene in a species: {sp}")

    # --- tree
    def nhx(**kw):
        return "[&&NHX:" + ":".join(f"{k}={v}" for k, v in kw.items()) + "]"

    a = {m["species"]: m["gene_id"] for m in mem if m["type"] == "1A"}
    kle = [a[s] for s in species if s in a and s in KLE]
    zt = [a[s] for s in species if s in a and s in ZT]
    post = [a[s] for s in species if s in a and s in POST_WGD]
    other = [a[s] for s in species if s in a and s not in KLE | ZT | POST_WGD]
    if other:
        raise SystemExit(f"URA1 in unexpected species: {other}")

    def grp(genes, ev, nid, tier, src):
        if len(genes) == 1:
            return genes[0]
        return f"({','.join(genes)}){nhx(ev=ev, id=nid, tier=tier, src=src)}"

    # valid on both panel trees (main: ZT sister to post-WGD; alt: KLE sister to post-WGD, i.e. the hybrid-parent
    # ambiguity of the WGD lineage): KLE and ZT URA1 are asserted orthologs (S, LCA Saccharomycetaceae in both trees);
    # post-WGD URA1 genes stay unresolved against everything (U) - pillar slots carry no track and the post-WGD
    # lineage may derive from either parent (allopolyploid WGD).
    pre = f"({grp(kle, 'S', 'URA1_KLE', 'T2', 'YGOBv7')},{grp(zt, 'S', 'URA1_ZT', 'T2', 'YGOBv7')})" \
          f"{nhx(ev='S', id='URA1_preWGD', tier='T2', src='YGOBv7')}"
    ura1 = f"({','.join([pre] + post)}){nhx(ev='U', id='URA1_root', tier='T2', src='YGOBv7', recipient=1)}"
    b = {m["gene_id"] for m in mem if m["type"] == "2"}
    clades = supported_clades(parse_newick(open(ROOT / BK_TREE).read().strip().rstrip(";")), b)
    gsp = {m["gene_id"]: m["species"] for m in mem}
    cl = sorted(clades, key=len)
    for x in cl:
        for y in cl:
            if x < y or y < x or not (x & y):
                continue
            if x != y:
                raise SystemExit("supported clades not laminar")
    used = set().union(*cl) if cl else set()
    t3 = 0

    def enc(genes: frozenset, top: bool) -> str:
        """Nested encoding of supported clades. One gene per species => species-overlap rule gives S; inside an S
        node, >=2 directly attached post-WGD genes go under a U node (possible pseudo-orthologs, not synteny-curated)."""
        nonlocal t3
        inner = [c for c in cl if c < genes and not any(c < d < genes for d in cl)]
        rest = sorted(genes - set().union(*inner) if inner else genes, key=lambda g: species.index(gsp[g]))
        kids = [enc(c, False) for c in sorted(inner, key=lambda c: -len(c))]
        if top:
            return f"({','.join(kids + rest)}){nhx(ev='U', id='URA9_root', tier='T2', src='Bouwknegt2021')}"
        pw = [g for g in rest if gsp[g] in POST_WGD]
        rest = [g for g in rest if gsp[g] not in POST_WGD]
        if len(pw) >= 2:
            kids.append(f"({','.join(pw)}){nhx(ev='U', id=f'URA9_postWGD_{t3 + 1}', tier='T3', src='Bouwknegt2021')}")
        else:
            rest += pw
        t3 += 1
        return (f"({','.join(kids + rest)})"
                f"{nhx(ev='S', id=f'URA9_T3_{t3}_BS{int(clades[genes])}', tier='T3', src='Bouwknegt2021')}")

    ura9 = enc(frozenset(b), True)
    tree = f"({ura1},{ura9}){nhx(ev='H', id='DHOD_HGT', tier='T1', src='Hall2005', donor='OUTSIDE', tag='HGT')};\n"
    (HERE / "reference_tree.nhx").write_text(tree)

    def w_tsv(name, cols, data):
        with open(HERE / name, "w", newline="") as fh:
            wr = csv.DictWriter(fh, cols, delimiter="\t", extrasaction="ignore", lineterminator="\n")
            wr.writeheader()
            wr.writerows(data)

    w_tsv("members.tsv", ["gene_id", "species", "original_id_in_paper", "label", "mapping_method", "mapping_identity",
                          "mapping_coverage", "source", "locator"], mem)
    w_tsv("mapping_report.tsv", ["original_id", "species", "label", "gene_id", "status", "reason", "locator", "source"],
          report)
    ura1_ids = ",".join(m["gene_id"] for m in mem if m["type"] == "1A")
    ura9_ids = ",".join(m["gene_id"] for m in mem if m["type"] == "2")
    t3_ids = ",".join(sorted(used))
    A = [
        ("event", "DHOD_HGT", "T1", "Hall2005", "fulltext.txt Introduction (summary of results)", Q["hgt"]),
        ("event", "DHOD_HGT", "T1", "Hall2005", "fulltext.txt Results heading", Q["donor"]),
        ("event", "DHOD_HGT", "T1", "Hall2005", "fulltext.txt Introduction", Q["donor2"]),
        ("support", "DHOD_HGT", "T3", "Gojkovic2004", "abstract.txt (Europe PMC abstract)", Q["goj2"]),
        ("tag", "DHOD_HGT", "T1", "Gojkovic2004", "abstract.txt (Europe PMC abstract)", Q["goj1"]),
        ("tag", f"type 1A (URA1): {ura1_ids}", "T1", "Hall2005", "fulltext.txt Results 'URA1 gene phylogeny'",
         Q["types"]),
        ("tag", f"type 2 (URA9-type): {ura9_ids}", "T1", "Hall2005", "fulltext.txt Results 'URA1 gene phylogeny'",
         Q["types"]),
        ("tag", "KLULA has both types", "T1", "Hall2005", "fulltext.txt Results 'Analysis of gene order'", Q["both"]),
        ("tag", "CANGA family 2 retained", "T1", "Hall2005", "fulltext.txt Results 'Analysis of gene order'",
         Q["loss"]),
        ("tag", "ASHGO family 2", "T1", "Hall2005", "fulltext.txt Introduction", Q["ashgo"]),
        ("tag", "URA9_root", "T1", "Hall2005", "fulltext.txt Introduction", Q["euk2"]),
        ("tag", f"URA1_root members {ura1_ids}", "T2", "YGOBv7", f"Pillars.TRUNCATED_wayback_20210728.tab line "
         f"{YGOB_LINE}", Q["ygob"]),
        ("tag", "URA1_root (post-WGD genes left U)", "T2", "YGOBv1", "v1_Synteny_Klactis.tab line 2061 "
         "(score ???!!?: S. cerevisiae/S. castellii positions non-syntenic)", Q["ygobv1"]),
        ("tag", f"URA9_root members {ura9_ids}", "T2", "Bouwknegt2021", "fulltext.txt Results / Fig. 1 legend",
         Q["bk1"]),
        ("tag", f"URA9_T3 clades {t3_ids}", "T3", "Bouwknegt2021",
         "Additional file 3 (40694_2021_117_MOESM3_ESM.txt): '(G8JWU1:...,Q75CE1:...)94'", Q["bk2"]),
    ]
    w_tsv("assertions.tsv", ["assertion_id", "type", "target", "tier", "source_id", "locator", "quote"],
          [dict(assertion_id=f"Y03_A{i + 1}", type=t, target=tg, tier=ti, source_id=s, locator=l, quote=q)
           for i, (t, tg, ti, s, l, q) in enumerate(A)])

    def files(lst):
        return [dict(path=p, sha256=sha256(ROOT / p)) for p in lst]

    cache = sorted(str(p.relative_to(ROOT)) for p in SC.iterdir())
    n_map = sum(r["status"] == "mapped" for r in report)
    n_den = sum(r["status"] in ("mapped", "unmapped", "ambiguous") for r in report)
    meta = dict(
        id="Y03_URA1_HGT",
        name="Dihydroorotate dehydrogenase: HGT-acquired family-1A URA1 vs native family-2 URA9-type DHOD",
        panel=PANEL,
        phenomena=["HGT", "LOSS"],
        descriptors=dict(n_members=len(mem), n_species=len({m["species"] for m in mem}),
                         n_type1A=sum(m["type"] == "1A" for m in mem), n_type2=sum(m["type"] == "2" for m in mem),
                         species_with_both=sorted(set(a) & {gsp[g] for g in b}), n_duplications=0, n_hgt=1,
                         mapped_fraction=round(n_map / n_den, 3)),
        sources=[
            dict(id="Hall2005", doi="10.1128/EC.4.6.1102-1115.2005", pmid="15947202",
                 title="Contribution of horizontal gene transfer to the evolution of Saccharomyces cerevisiae",
                 year=2005, source_method="phylogeny", locator="Introduction; Results (DHOD phylogeny, gene order)",
                 quote=Q["hgt"], files=files([HALL] + [c for c in cache if "seqmap_" in c or "genpept_" in c])),
            dict(id="Gojkovic2004", doi="10.1007/s00438-004-0995-7", pmid="15014982",
                 title="Horizontal gene transfer promoted evolution of the ability to propagate under anaerobic "
                       "conditions in yeasts", year=2004, source_method="phylogeny", locator="abstract",
                 quote=Q["goj2"], files=files([GOJ])),
            dict(id="Bouwknegt2021", doi="10.1186/s40694-021-00117-4", pmid="34656184",
                 title="Class-II dihydroorotate dehydrogenases from three phylogenetically distant fungi support "
                       "anaerobic pyrimidine biosynthesis", year=2021, source_method="phylogeny",
                 locator="Additional files 1 (ortholog FASTA) and 3 (RAxML-NG tree)", quote=Q["bk1"],
                 files=files([f"{BK}/fulltext.txt", BK_IDS, BK_TREE])),
            dict(id="YGOBv7", doi="10.1101/gr.3672305", pmid="16169922",
                 title="YGOB v7 (Aug 2012) Pillars.tab (Wayback 2021-07-28, truncated)", year=2012,
                 source_method="synteny", locator=f"Pillars.tab line {YGOB_LINE}", quote=Q["ygob"],
                 files=files([YGOB, "literature/raw/B1/ygob/README", str((SC / "uniprot_ygob_ura1_loci.tsv")
                                                                        .relative_to(ROOT))])),
            dict(id="YGOBv1", doi="10.1101/gr.3672305", pmid="16169922",
                 title="YGOB v1 (Sep 2005) Synteny_Klactis.tab", year=2005, source_method="synteny",
                 locator="v1_Synteny_Klactis.tab line 2061", quote=Q["ygobv1"],
                 files=files([YGOBV1, "literature/raw/B1/ygob/wayback/v1_README.txt"])),
        ],
        root_single_copy=False,
        complete_species=[],
        decoys=[],
        distinct_from=[],
        tier_summary=dict(T1_events=1, T1_H=1, T2_nodes=4, T3_nodes=t3, events_D=0),
        protocol_version=1,
        notes=("Root H (T1, Hall 2005; Gojkovic 2004): URA1 (family 1A DHOD) entered the Saccharomycetaceae by HGT from "
               "Lactobacillales (donor OUTSIDE panel) -> URA1 xenologs and native family-2 (URA9-type) DHODs are never "
               "in one HOG at any level. URA1 subtree = YGOB pillar (synteny): KLE pair, ZT pair and KLE-vs-ZT are S (T2); "
               "post-WGD URA1 genes (YEAST, NAUCC, TETPH) hang from a U node (valid under both panel trees) (pillar slots carry no track; YGOB v1 "
               "scores URA1 as non-syntenic in S. cerevisiae/S. castellii). URA9 subtree = Bouwknegt 2021 Ura9 "
               "orthologs + Hall 2005 accessions; U except the Eremothecium pair (BS 94, T3). Species with both "
               "types: LACTH, TETPH (KLULA has both per Hall 2005 but its family-2 gene has no identifier in any "
               "source -> only URA1 listed). No member for NAUDC, KAZAF. YGOB unit YGOB_L04856 of Y01 holds the "
               "same 7 URA1 genes (gene-ownership overlap; one of the two must drop them in the panel set)."),
    )
    (HERE / "family.yaml").write_text(yaml.safe_dump(meta, sort_keys=False, allow_unicode=True, width=120))
    print(f"Y03: {len(mem)} members ({meta['descriptors']['n_type1A']} URA1, {meta['descriptors']['n_type2']} URA9), "
          f"{meta['descriptors']['n_species']} spp; mapped {n_map}/{n_den}; T3 clades {t3}")


if __name__ == "__main__":
    main()
