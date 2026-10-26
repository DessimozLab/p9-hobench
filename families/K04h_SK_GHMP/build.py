#!/usr/bin/env python
"""Build one K04 isofunctional-decoy family (panel K). The SAME script sits in every families/K04?_*/ directory;
the family built is chosen by the directory name. It writes only into its own directory.

Biology / sources
  Omelchenko2010  Omelchenko MV, Galperin MY, Wolf YI, Koonin EV (2010) Biol Direct 5:31 (PMID 20433725, OA, CC BY):
                  non-homologous isofunctional enzymes (NISE). Additional file 2 (Table S2, final list of confirmed
                  NISE, "Two distinct structural folds") gives, per EC node, the structurally unrelated forms with
                  example UniProt entries, SCOP fold and Pfam family. Pairs used here (all "Two distinct structural
                  folds" in Table S2):
                    6.1.1.6  lysyl-tRNA synthetase class I  vs class II
                    1.15.1.1 Fe,Mn superoxide dismutase      vs Cu,Zn superoxide dismutase
                    5.4.2.1  2,3-BPG-dependent PGM          vs 2,3-BPG-independent PGM
                    2.7.1.71 P-loop shikimate kinase        vs archaeal GHMP-kinase shikimate kinase
                    4.2.1.10 type I (class I aldolase fold) vs type II (flavodoxin-like) 3-dehydroquinate dehydratase
  COG2020         NCBI COG 2020 (Galperin et al. 2021): the COG whose definition names each form
                  (cog-20.def.tab) and its membership in the panel-K genomes (derived/cog20_panelK.csv), mapped onto
                  the pinned UniProt proteomes by same-species sequence identity (derived/seqmap_K04/, exact sha256 or
                  unique DIAMOND hit >=98 % id, >=80 % coverage of both, 2-point ambiguity margin).
Encoding: members = COG members (membership class 0/1/2; class 3 = partial on both -> excluded); tree = one U node over
all members (no source in hand resolves orthology vs xenology inside these groups; several are HGT-prone and several
species have >1 copy). The decoy relation is the point of the family: `decoys` = all members of the partner family,
`distinct_from` = [{partner, at_level: ALL}] (non-homologous: must never share a HOG at any level).
BRADU and GEOSL are not in COG 2020 -> no members can be derived for them (listed in mapping_report.tsv).
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
COGDEF = f"{CP}/raw/cog-20.def.tab"
SEQMAP = f"{CP}/derived/seqmap_K04"
OM = "literature/papers/PMID20433725"
OMS2 = f"{OM}/raw/supp/1745-6150-5-31-S2.PDF"
OMS1 = f"{OM}/raw/supp/1745-6150-5-31-S1.PDF"
PANEL_SPECIES = ["ECOLI", "PSEAE", "BRADU", "HELPY", "GEOSL", "BACSU", "MYCTU", "CHLTR", "SYNY3", "DEIRA", "THEMA",
                 "AQUAE", "BACTN", "LEPIN", "SACS2", "KORCO", "NITMS", "THEKO", "METJA", "METAC", "HALSA"]
NOT_IN_COG = ["BRADU", "GEOSL"]

# family id -> (partner id, COG, COG definition quote, EC, form name, Table S2 quote of this form, pair quote)
FAM = {
    "K04a_LysRS_I": dict(
        partner="K04b_LysRS_II", cog="COG1384", cogq="COG1384\tJ\tLysyl-tRNA synthetase, class I\tLysS\tAminoacyl-tRNA synthetases\t\t1IRX",
        ec="6.1.1.6", form="Class I lysyl-tRNA synthetase (LysRS-I)",
        s2q="6.1.1.6 Lysine—tRNA ligase Q9YFT9 SYK_AERPE 1IRX Adenine nucleotide alpha hydrolase-like SupFam "
            "Class I lysyl-tRNA synthetase"),
    "K04b_LysRS_II": dict(
        partner="K04a_LysRS_I", cog="COG1190", cogq="COG1190\tJ\tLysyl-tRNA synthetase (class II)\tLysU\tAminoacyl-tRNA synthetases\t\t1BBU",
        ec="6.1.1.6", form="Class II lysyl-tRNA synthetase (LysRS-II)",
        s2q="6.1.1.6 Lysine—tRNA ligase P13030 SYK1_ECOLI 1BBU Class II aaRS and biotin synthetases SupFam Class II "
            "aminoacyl-tRNA synthetase"),
    "K04c_SOD_FeMn": dict(
        partner="K04d_SOD_CuZn", cog="COG0605", cogq="COG0605\tP\tSuperoxide dismutase\tSodA\t\t\t1AP5",
        ec="1.15.1.1", form="Fe/Mn superoxide dismutase",
        s2q="1.15.1.1 Superoxide dismutase P09157 SODF_ECOLI 1ISA FeMn superoxide dismutase SupFam Fe,Mn superoxide "
            "dismutase"),
    "K04d_SOD_CuZn": dict(
        partner="K04c_SOD_FeMn", cog="COG2032", cogq="COG2032\tP\tCu/Zn superoxide dismutase\tSodC\t\t\t1ESO",
        ec="1.15.1.1", form="Cu/Zn superoxide dismutase",
        s2q="1.15.1.1 Superoxide dismutase P53635 SODC_ECOLI 1ESO Immunoglobulin-like beta-sandwich SupFam Cu,Zn "
            "superoxide dismutase"),
    "K04e_PGM_BPGdep": dict(
        partner="K04f_PGM_BPGindep", cog="COG0588", cogq="COG0588\tG\tPhosphoglycerate mutase (BPG-dependent)\tGpmA\tGlycolysis\t\t1E58",
        ec="5.4.2.1", form="2,3-bisphosphoglycerate-dependent phosphoglycerate mutase (dPGM, GpmA)",
        s2q="5.4.2.1 Phosphoglycerate mutase P31217 GPMA_ECOLI 1E59 Phosphoglycerate mutase-like SupFam "
            "Phosphoglycerate mutase"),
    "K04f_PGM_BPGindep": dict(
        partner="K04e_PGM_BPGdep", cog="COG0696", cogq="COG0696\tG\tPhosphoglycerate mutase (BPG-independent), AlkP superfamily\tGpmI\tGlycolysis\t\t1EJJ",
        ec="5.4.2.1", form="2,3-bisphosphoglycerate-independent phosphoglycerate mutase (iPGM, GpmI)",
        s2q="5.4.2.1 Phosphoglycerate mutase P52832 GPMI_PSESM 1EJJ Alkaline phosphatase-like SupFam 2,3-"
            "Bisphosphoglycerate-independent phosphoglycerate mutase"),
    "K04g_SK_Ploop": dict(
        partner="K04h_SK_GHMP", cog="COG0703", cogq="COG0703\tE\tShikimate kinase\tAroK\tAromatic amino acid biosynthesis\t\t1E6C",
        ec="2.7.1.71", form="P-loop shikimate kinase (bacterial AroK/AroL)",
        s2q="2.7.1.71 Shikimate kinase P0A6E1 AROL_ECOLI 1SHK P-loop containing NTPases SupFam Shikimate kinase"),
    "K04h_SK_GHMP": dict(
        partner="K04g_SK_Ploop", cog="COG1685", cogq="COG1685\tE\tArchaeal shikimate kinase\tAroK2",
        ec="2.7.1.71", form="GHMP-kinase-type archaeal shikimate kinase",
        s2q="2.7.1.71 Shikimate kinase Q58835 AROK_METJA similar to 1VIS Ribosomal protein S5 domain 2-like and "
            "Ferredoxin-like SupFam GHMP kinase"),
    "K04i_DHQase_I": dict(
        partner="K04j_DHQase_II", cog="COG0710", cogq="COG0710\tE\t3-dehydroquinate dehydratase\tAroD\tAromatic amino acid biosynthesis\t\t1GQN",
        ec="4.2.1.10", form="Type I 3-dehydroquinate dehydratase (AroD, class I aldolase fold)",
        s2q="4.2.1.10 3-dehydroquinate dehydratase P05194 AROD_ECOLI 1QFE TIM beta/alpha-barrel SupFam Class I "
            "aldolase"),
    "K04j_DHQase_II": dict(
        partner="K04i_DHQase_I", cog="COG0757", cogq="COG0757\tE\t3-dehydroquinate dehydratase\tAroQ\tAromatic amino acid biosynthesis\t\t1D0I",
        ec="4.2.1.10", form="Type II 3-dehydroquinate dehydratase (AroQ, flavodoxin-like fold)",
        s2q="4.2.1.10 3-dehydroquinate dehydratase P43877 AROQ_ACTPL 1UQR Flavodoxin-like SupFam Type II "
            "3-dehydroquinate dehydratase"),
}
# Table S2 lists, for each EC node, the forms and in the 'Predicted' column "Two distinct structural folds"; the pair
# quote below spans the first form's row and that verdict (elision allowed by the provenance checker).
PAIRQ = {
    "6.1.1.6": "6.1.1.6 Lysine—tRNA ligase P13030 SYK1_ECOLI 1BBU Class II aaRS and biotin synthetases SupFam Class II "
               "aminoacyl-tRNA synthetase 7735833, 11887185, 2183178 Two distinct structural folds",
    "1.15.1.1": "1.15.1.1 Superoxide dismutase P09157 SODF_ECOLI 1ISA FeMn superoxide dismutase SupFam Fe,Mn superoxide "
                "dismutase 2447093 Two distinct structural folds",
    "5.4.2.1": "5.4.2.1 Phosphoglycerate mutase P31217 GPMA_ECOLI 1E59 Phosphoglycerate mutase-like SupFam "
               "Phosphoglycerate mutase 11038361 Two distinct structural folds",
    "2.7.1.71": "2.7.1.71 Shikimate kinase P0A6E1 AROL_ECOLI 1SHK P-loop containing NTPases SupFam Shikimate kinase "
                "3001029, 9600856, 11369852 Two distinct structural",
    "4.2.1.10": "4.2.1.10 3-dehydroquinate dehydratase P43877 AROQ_ACTPL 1UQR Flavodoxin-like SupFam Type II "
                "3-dehydroquinate dehydratase 14993670, 8170389, 2848727 Two distinct structural folds",
}
QDEF = ("The final list of NISE, obtained after a detailed case-by-case analysis, includes 140 confirmed EC nodes where "
        "the presence of two or more distinct forms is supported by a comparison of their 3D structures that revealed "
        "different folds (Additional file 2, Table S2).")
QBONA = ("Here, we focus on enzyme variants that possess (or could be inferred to possess) different structural folds "
         "and therefore appear to be evolutionarily unrelated (bona fide NISE).")
QCOG = "Comma-delimited plain text file assigning proteins to COGs"


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


_FCACHE: dict[str, dict] = {}


def derive(fid: str) -> tuple[list[dict], list[dict]]:
    """Members and mapping report of family `fid` (pure function of the stored sources)."""
    f = FAM[fid]
    cog = f["cog"]
    cogseq = read_fasta(ROOT / COGFA)
    rows = [r for r in csv.reader(open(ROOT / COGCSV)) if r[6] == cog]
    seqmap = {}
    for p in sorted((ROOT / SEQMAP).glob("*.tsv")):
        for r in csv.DictReader(open(p), delimiter="\t"):
            seqmap[(p.stem, r["query_id"])] = r
    report, by_target = [], defaultdict(list)
    for r in rows:
        code, pid, cls = r[-1], r[2], r[8]
        loc = f"cog20_panelK.csv: {r[0]},{r[1]},{pid},{cog},class={cls}"
        rep = dict(original_id=pid, species=code, label=f"{f['form']} ({cog})", gene_id="", status="", reason="",
                   locator=loc)
        report.append(rep)
        if cls == "3":
            rep.update(status="excluded", reason="COG membership class 3 (partial match on protein and COG profile)")
            continue
        if code not in _FCACHE:
            _FCACHE[code] = read_fasta(panel_fasta(code))
        exact = [k for k, s in _FCACHE[code].items() if s == cogseq[pid]]
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
        by_target[(code, g)].append((rep, dict(gene_id=g, species=code, original_id_in_paper=pid,
                                               label=f"{f['form']} ({cog})", mapping_method=meth,
                                               mapping_identity=f"{ident:.1f}", mapping_coverage=f"{cov:.3f}",
                                               source="COG2020", locator=loc)))
    members = []
    for (code, g), lst in sorted(by_target.items()):
        lst.sort(key=lambda x: (x[1]["mapping_method"] != "exact_sequence", x[1]["original_id_in_paper"]))
        members.append(lst[0][1])
        for rep, _ in lst[1:]:
            rep.update(status="excluded", reason=f"collapsed copy: maps to the same UniProt entry {g} as "
                                                 f"{lst[0][1]['original_id_in_paper']} (protocol rule 6)")
    members.sort(key=lambda m: (PANEL_SPECIES.index(m["species"]), m["gene_id"]))
    for code in NOT_IN_COG:
        report.append(dict(original_id=f"{code}:{cog}", species=code, label=f"{f['form']} ({cog})", gene_id="",
                           status="not_in_panel",
                           reason="genome not in COG 2020 (only other species of the genus); membership unknown",
                           locator="cog-20.org.csv (absent)"))
    report.sort(key=lambda r: (PANEL_SPECIES.index(r["species"]), r["original_id"]))
    return members, report


def main() -> None:
    fid = HERE.name
    if fid not in FAM:
        raise SystemExit(f"unknown K04 family {fid}")
    f = FAM[fid]
    members, report = derive(fid)
    partner_members, _ = derive(f["partner"])
    decoys = sorted(m["gene_id"] for m in partner_members)
    clash = set(decoys) & {m["gene_id"] for m in members}
    if clash:
        raise SystemExit(f"{fid}: genes in both NISE forms: {sorted(clash)}")

    def w_tsv(name, cols, data):
        with open(HERE / name, "w", newline="") as fh:
            wr = csv.DictWriter(fh, cols, delimiter="\t", extrasaction="ignore", lineterminator="\n")
            wr.writeheader()
            wr.writerows(data)

    w_tsv("mapping_report.tsv", ["original_id", "species", "label", "gene_id", "status", "reason", "locator"], report)
    w_tsv("members.tsv", ["gene_id", "species", "original_id_in_paper", "label", "mapping_method", "mapping_identity",
                          "mapping_coverage", "source", "locator"], members)
    root_id = fid.split("_", 1)[0] + "_root"
    leaves = [m["gene_id"] for m in members]
    if len(leaves) >= 2:
        tree = f"({','.join(leaves)})[&&NHX:ev=U:id={root_id}:tier=T2:src=COG2020];\n"
    else:
        tree = f"{leaves[0]};\n"
    (HERE / "reference_tree.nhx").write_text(tree)

    ida = fid.split("_", 1)[0]
    assertions = [
        dict(assertion_id=f"{ida}_A1", type="decoy", target=",".join(decoys), tier="T1",
             source_id="Omelchenko2010", locator=f"Additional file 2 (Table S2), EC {f['ec']} rows (decoys = all members of {f['partner']})",
             quote=PAIRQ[f["ec"]]),
        dict(assertion_id=f"{ida}_A2", type="distinct_from", target=f"{f['partner']} at_level=ALL", tier="T1",
             source_id="Omelchenko2010", locator="fulltext.txt Results (final NISE list = Table S2)", quote=QDEF),
        dict(assertion_id=f"{ida}_A3", type="distinct_from", target=f"{f['partner']} at_level=ALL", tier="T1",
             source_id="Omelchenko2010", locator="fulltext.txt Results (bona fide NISE)", quote=QBONA),
        dict(assertion_id=f"{ida}_A4", type="tag", target=f"{fid} form", tier="T2", source_id="Omelchenko2010",
             locator=f"Additional file 2 (Table S2), EC {f['ec']}, row of this form", quote=f["s2q"]),
        dict(assertion_id=f"{ida}_A5", type="tag", target=f"{fid} members = {f['cog']}", tier="T2",
             source_id="COG2020", locator=f"raw/cog-20.def.tab, row {f['cog']}", quote=f["cogq"]),
        dict(assertion_id=f"{ida}_A6", type="tag", target=f"{fid} members = {f['cog']}", tier="T2",
             source_id="COG2020", locator="raw/Readme.2020-11-25.txt (cog-20.cog.csv)", quote=QCOG),
    ]
    w_tsv("assertions.tsv", ["assertion_id", "type", "target", "tier", "source_id", "locator", "quote"], assertions)

    def files(lst):
        return [dict(path=p, sha256=sha256(ROOT / p)) for p in lst]

    n_map = sum(r["status"] == "mapped" for r in report)
    n_den = sum(r["status"] in ("mapped", "unmapped", "ambiguous") for r in report)
    cnt = defaultdict(int)
    for m in members:
        cnt[m["species"]] += 1
    meta = dict(
        id=fid,
        name=f"{f['form']} — non-homologous isofunctional counterpart of {f['partner']} (EC {f['ec']})",
        panel=PANEL,
        collection="K04_isofunctional_decoys",
        phenomena=["ANALOG"],
        descriptors=dict(n_members=len(members), n_species=len(cnt), n_duplications=0, ec=f["ec"], cog=f["cog"],
                         multicopy_species=sorted(s for s, c in cnt.items() if c > 1),
                         mapped_fraction=round(n_map / n_den, 3) if n_den else None),
        sources=[
            dict(id="Omelchenko2010", doi="10.1186/1745-6150-5-31", pmid="20433725",
                 title="Non-homologous isofunctional enzymes: a systematic analysis of alternative solutions in "
                       "enzyme evolution", year=2010, source_method="expert_classification",
                 locator=f"Additional file 2 (Table S2), EC {f['ec']}", quote=PAIRQ[f["ec"]],
                 files=files([f"{OM}/fulltext.txt", OMS2, OMS1])),
            dict(id="COG2020", doi="10.1093/nar/gkaa1018", pmid="33167031",
                 title="COG database update: focus on microbial diversity, model organisms, and widespread pathogens",
                 year=2021, source_method="expert_classification", locator=f"cog-20.def.tab {f['cog']}",
                 quote=f["cogq"].replace("\t", " "),
                 files=files([COGCSV, COGFA, COGDEF, f"{CP}/raw/Readme.2020-11-25.txt", f"{CP}/seqmap_K04.py"]
                             + sorted(str(p.relative_to(ROOT)) for p in (ROOT / SEQMAP).glob("*.tsv")))),
        ],
        root_single_copy=False,
        complete_species=[],
        decoys=decoys,
        distinct_from=[dict(family=f["partner"], at_level="ALL", src="Omelchenko2010")],
        tier_summary=dict(T1_decoy_assertions=3, T1_events=0, T2_nodes=1, T3=0),
        protocol_version=1,
        notes=(f"Merge-test decoy family. {f['form']} (COG {f['cog']}) and its partner {f['partner']} catalyse the same "
               f"reaction (EC {f['ec']}) with different structural folds (Omelchenko et al. 2010, Table S2: 'Two "
               f"distinct structural folds') => non-homologous: they must never share a HOG at any level "
               f"(distinct_from at_level ALL; decoys = partner members). Inside the family the tree is a single U "
               f"node (no orthology asserted: no source in hand separates orthologs from xenologs/paralogs; LysRS-I "
               f"and several SODs/SKs are HGT-prone). Members = COG 2020 membership (class 0-2) in the 19 panel-K "
               f"genomes present in COG 2020; BRADU/GEOSL absent from COG 2020 (no member data). Same pairs can be "
               f"re-declared in panel E."),
    )
    (HERE / "family.yaml").write_text(yaml.safe_dump(meta, sort_keys=False, allow_unicode=True, width=120))
    print(f"{fid}: {len(members)} members / {len(cnt)} spp, mapped {n_map}/{n_den}, decoys {len(decoys)}")


if __name__ == "__main__":
    main()
