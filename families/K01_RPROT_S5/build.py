#!/usr/bin/env python
"""Build one K01 ribosomal-protein unit (collection RPROT) for panel K.

The unit is taken from this directory's name: families/K01_RPROT_<name>/ where <name> is the bacterial r-protein name
of Yutin et al. 2012 Table S1 (or the archaeal name for archaeal-only rows). The SAME build.py is copied into every
K01_RPROT_* directory; it only writes into its own directory.

Sources (all read from stored files, nothing typed by hand except quotes/rules):
  Yutin2012  Yutin N, Puigbo P, Koonin EV, Wolf YI (2012) PLoS One 7:e36972 (PMID 22615861, PMC3353972, CC0)
             Table S1 (s001.xls) = r-protein list with COG ids; Table S3 (s003.xls) = per-genome copy-number matrix
             (main set + 2nd paralogs); Table S2 (s002.xls) = genome list.
  COG2020    NCBI COG 2020 (Galperin et al. 2021 NAR 49:D274), ftp.ncbi.nih.gov/pub/COG/COG2020/data (md5 verified);
             derived/cog20_panelK.csv + .fa = rows/sequences of the panel-K genomes (extract_panelK.py),
             derived/seqmap/<CODE>.tsv = same-species sequence cross-reference onto the UniProt 2026_03 panel proteomes
             (seqmap_panelK.py -> curation_tools.map_by_sequence).
Rules (CURATION_PROTOCOL v1):
  * members = COG 2020 members (Yutin's COG id for this r-protein) in the panel-K genomes present in COG 2020;
  * mapping accepted if exact sequence identity (sha256) or unique DIAMOND hit with >=98 % identity, >=80 % coverage
    of both sequences, no other target within 2 identity points; COG membership class 3 (partial match on both
    protein and profile) -> excluded; two COG proteins collapsing onto one UniProt entry (identical copies) -> one
    member, the other reported 'excluded (identical copy)';
  * the unit is a single-copy control only if every panel species has <=1 member AND Yutin Table S3 lists no
    2nd paralog for that genome/protein; otherwise the unit is not built (status recorded in curation_log.md);
  * tree = one S polytomy (T2 rule b: same label, <=1 gene per species), root_single_copy (Yutin: ancestral r-proteins);
  * complete_species = genomes searched genome-wide by Yutin (six-frame PSSM screen) whose Table S3 count equals the
    number of uniquely mapped members (>=1).
"""
from __future__ import annotations

import csv
import glob
import hashlib
import sys
from collections import defaultdict
from pathlib import Path

import xlrd
import yaml

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
UNIT = HERE.name.split("K01_RPROT_", 1)[1]
YP = "literature/papers/PMID22615861"
S1 = f"{YP}/raw/supp/pone.0036972.s001.xls"
S3 = f"{YP}/raw/supp/pone.0036972.s003.xls"
CP = "literature/papers/COG2020"
COGCSV = f"{CP}/derived/cog20_panelK.csv"
COGFA = f"{CP}/derived/cog20_panelK.fa"
PANEL = "K"
PANEL_SPECIES = ["ECOLI", "PSEAE", "BRADU", "HELPY", "GEOSL", "BACSU", "MYCTU", "CHLTR", "SYNY3", "DEIRA", "THEMA",
                 "AQUAE", "BACTN", "LEPIN", "SACS2", "KORCO", "NITMS", "THEKO", "METJA", "METAC", "HALSA"]
ARCHAEA = {"SACS2", "KORCO", "NITMS", "THEKO", "METJA", "METAC", "HALSA"}
# panel code -> genome name in Yutin Table S2/S3 (same strain; taxid identical to the panel proteome except ECOLI
# K-12 MG1655 = panel K-12 and HALSA 'Halobacterium_sp' = NRC-1, taxid 64091). SYNY3 is not in Yutin's genome set.
YUTIN_GENOME = {
    "ECOLI": "Escherichia_coli_K_12_substr__MG1655", "PSEAE": "Pseudomonas_aeruginosa", "BRADU": "Bradyrhizobium_japonicum",
    "HELPY": "Helicobacter_pylori_26695", "GEOSL": "Geobacter_sulfurreducens", "BACSU": "Bacillus_subtilis",
    "MYCTU": "Mycobacterium_tuberculosis_H37Rv", "CHLTR": "Chlamydia_trachomatis_D_UW_3_CX",
    "DEIRA": "Deinococcus_radiodurans", "THEMA": "Thermotoga_maritima", "AQUAE": "Aquifex_aeolicus",
    "BACTN": "Bacteroides_thetaiotaomicron_VPI-5482", "LEPIN": "Leptospira_interrogans_serovar_Lai",
    "SACS2": "Sulfolobus_solfataricus", "KORCO": "Candidatus_Korarchaeum_cryptofilum_OPF8",
    "NITMS": "Nitrosopumilus_maritimus_SCM1", "THEKO": "Thermococcus_kodakaraensis_KOD1",
    "METJA": "Methanococcus_jannaschii", "METAC": "Methanosarcina_acetivorans", "HALSA": "Halobacterium_sp",
}
IN_COG = set(PANEL_SPECIES) - {"BRADU", "GEOSL"}   # B. diazoefficiens USDA110 / G. sulfurreducens PCA not in COG 2020
EXCLUDED_UNITS = {
    "L7L12": "Table S1 assigns COG2058 (archaeal L12E/P1/P2 family) to bacterial L7/L12; in COG 2020 COG2058 has no "
             "bacterial member in the panel genomes (Yutin Table S3 lists L7L12 in every bacterium), so the COG cannot "
             "serve as the membership source for this universal unit",
    "S1": "Yutin 2012 excluded S1p from the study (varied domain architecture, ubiquitous S1-like domains)",
}
Q_SCREEN = ("position-specific scoring matrices (PSSMs) for 56 bacterial and 71 archaeal r-proteins [6], [7] "
            "(Table S1) were used to screen completely sequenced prokaryotic genomes (Table S2) translated in six frames")
Q_UNIVERSAL = ("Archaeal and bacterial ribosomes contain more than 50 proteins, including 34 that are universally "
               "conserved in the three domains of cellular life (bacteria, archaea, and eukaryotes).")
Q_COG = "NCBI COG [54] and arCOG [55] databases were used as sources of initial sets of 56 bacterial and 68 archaeal r-proteins"
Q_COG2020 = "Comma-delimited plain text file assigning proteins to COGs"


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


def table_s1() -> list[dict]:
    sh = xlrd.open_workbook(str(ROOT / S1)).sheet_by_index(0)
    out, section = [], None
    for r in range(sh.nrows):
        v = [str(x).strip() for x in sh.row_values(r)]
        for key, sec in (("Universal proteins", "universal"), ("Bacterial-only", "bacterial"), ("Archaeal-only", "archaeal")):
            if v[1].startswith(key):
                section = sec
        if not section or not (v[0] or v[1]) or v[1].startswith(("Universal", "Bacterial-only", "Archaeal-only")):
            continue
        cog = v[2]
        cog_n = ("COG" + cog[3:].lstrip("0").zfill(4)) if cog.startswith("COG") else None
        out.append(dict(row=r + 1, arch=v[0], bact=v[1], cog_raw=cog, cog=cog_n, cog_annot=v[3], section=section,
                        quote=" ".join(x for x in v[:4] if x)))
    return out


def table_s3() -> tuple[dict, dict]:
    """(main[(code, name)] -> count, second[(code, name)] -> count) for the panel genomes, with (sheet,row) locators."""
    wb = xlrd.open_workbook(str(ROOT / S3))
    inv = {v: k for k, v in YUTIN_GENOME.items()}
    main, second = {}, {}
    for shn, hr, dest in [("main_set_only_bact", 1, main), ("2nd_paralogs only_bact", 0, second),
                          ("main_set_only_arch", 1, main), ("2nd_paralogs_only_archaea", 0, second)]:
        sh = wb.sheet_by_name(shn)
        hdr = [str(x).strip() for x in sh.row_values(hr)]
        for r in range(hr + 1, sh.nrows):
            v = sh.row_values(r)
            code = inv.get(str(v[0]).strip())
            if not code:
                continue
            for i in range(1, len(hdr)):
                if hdr[i] and hdr[i] != "total":
                    c = int(float(v[i])) if v[i] not in ("", None) else 0
                    dest[(code, hdr[i])] = (c, f"{Path(S3).name} sheet '{shn}' row {r + 1} col '{hdr[i]}'")
    return main, second


def main() -> None:
    units = {(u["bact"] or u["arch"]): u for u in table_s1()}
    if UNIT not in units:
        raise SystemExit(f"unit {UNIT} not in Table S1")
    u = units[UNIT]
    if UNIT in EXCLUDED_UNITS:
        (HERE / "unit_status.txt").write_text(f"excluded_by_curator\n{EXCLUDED_UNITS[UNIT]}\n")
        for f in ("members.tsv", "reference_tree.nhx", "family.yaml", "assertions.tsv", "mapping_report.tsv"):
            (HERE / f).unlink(missing_ok=True)
        print(f"{UNIT}: NOT BUILT ({EXCLUDED_UNITS[UNIT]})")
        return
    if not u["cog"]:
        raise SystemExit(f"unit {UNIT}: no COG id in Table S1 ({u['cog_raw']!r})")
    cog = u["cog"]
    s3_main, s3_second = table_s3()
    names = {"bact": u["bact"], "arch": u["arch"]}

    def yutin_count(code):
        nm = names["arch"] if code in ARCHAEA else names["bact"]
        if code not in YUTIN_GENOME or not nm or (code, nm) not in s3_main:
            return None
        m, loc = s3_main[(code, nm)]
        s, loc2 = s3_second.get((code, nm), (0, ""))
        return m + s, loc + ("; " + loc2 if s else "")

    cogseq = read_fasta(ROOT / COGFA)
    rows = [r for r in csv.reader(open(ROOT / COGCSV)) if r[6] == cog]
    seqmap = {}
    for code in IN_COG:
        p = ROOT / f"{CP}/derived/seqmap/{code}.tsv"
        if p.exists():
            for r in csv.DictReader(open(p), delimiter="\t"):
                seqmap[(code, r["query_id"])] = r
    fasta_cache: dict[str, dict] = {}

    n_cog_sp = defaultdict(int)
    for r in rows:
        n_cog_sp[r[-1]] += 1
    report, members, by_target = [], [], defaultdict(list)
    for r in rows:
        code, pid, cls = r[-1], r[2], r[8]
        loc = f"cog20_panelK.csv: {r[0]},{r[1]},{pid},{cog},class={cls}"
        rep = dict(original_id=pid, species=code, label=f"{UNIT} ({cog})", gene_id="", status="", reason="",
                   locator=loc)
        report.append(rep)
        if u["section"] == "bacterial" and code in ARCHAEA:
            rep.update(status="excluded", reason=f"archaeal member of {cog}: Yutin 2012 Table S1 lists {UNIT} as a "
                                                 f"bacterial-only r-protein (archaeal r-proteins are listed separately)")
            continue
        yc0 = yutin_count(code)
        if cls == "3" and not (n_cog_sp[code] == 1 and yc0 and yc0[0] == 1):
            rep.update(status="excluded", reason="COG membership class 3 (partial match on protein and COG profile) "
                                                 "not corroborated by a single-copy Yutin 2012 Table S3 entry")
            continue
        if code not in fasta_cache:
            fasta_cache[code] = read_fasta(panel_fasta(code))
        tgt = fasta_cache[code]
        exact = [k for k, s in tgt.items() if s == cogseq[pid]]
        sm = seqmap.get((code, pid))
        if len(exact) == 1:
            g, meth, ident, cov = exact[0], "exact_sequence", 100.0, 1.0
        elif len(exact) > 1:
            rep.update(status="ambiguous", reason=f"identical sequence in {len(exact)} panel entries: {','.join(exact)}")
            continue
        elif sm and sm["status"] == "unique" and float(sm["mapping_identity"]) >= 98 and \
                min(float(sm["query_cov"]), float(sm["target_cov"])) >= 0.8:
            g, meth = sm["target_id"], "diamond_sequence"
            ident, cov = float(sm["mapping_identity"]), min(float(sm["query_cov"]), float(sm["target_cov"]))
        else:
            why = (f"best hit {sm['target_id']} status={sm['status']} identity={sm['mapping_identity']} "
                   f"qcov={sm['query_cov']} tcov={sm['target_cov']} ({sm['flags']})") if sm else "no sequence hit"
            rep.update(status="unmapped", reason=f"below protocol thresholds: {why}")
            continue
        rep.update(status="mapped", gene_id=g, reason=meth)
        by_target[(code, g)].append((rep, dict(gene_id=g, species=code, original_id_in_paper=pid,
                                               label=f"{UNIT} ({cog})", mapping_method=meth,
                                               mapping_identity=f"{ident:.1f}", mapping_coverage=f"{cov:.3f}",
                                               source="COG2020", locator=loc)))
    for (code, g), lst in sorted(by_target.items()):
        lst.sort(key=lambda x: (x[1]["mapping_method"] != "exact_sequence", x[1]["original_id_in_paper"]))
        members.append(lst[0][1])
        for rep, _ in lst[1:]:
            rep.update(status="excluded", reason=f"identical copy: maps to the same UniProt entry {g} as "
                                                 f"{lst[0][1]['original_id_in_paper']} (protocol rule 6)")
    # Yutin presence for genomes without COG 2020 members
    have = defaultdict(int)
    for m in members:
        have[m["species"]] += 1
    for code in PANEL_SPECIES:
        yc = yutin_count(code)
        if code not in IN_COG and yc and yc[0] > 0:
            report.append(dict(original_id=f"{YUTIN_GENOME[code]}:{UNIT}", species=code, label=f"{UNIT} ({cog})",
                               gene_id="", status="unmapped",
                               reason=f"present in Yutin 2012 Table S3 (count {yc[0]}) but genome not in COG 2020: "
                                      f"no protein identifier in any source", locator=yc[1]))
        elif code in IN_COG and yc and yc[0] > 0 and have[code] == 0 and \
                not any(r["species"] == code for r in report):
            report.append(dict(original_id=f"{YUTIN_GENOME[code]}:{UNIT}", species=code, label=f"{UNIT} ({cog})",
                               gene_id="", status="unmapped",
                               reason=f"present in Yutin 2012 Table S3 (count {yc[0]}) but no {cog} member in COG 2020",
                               locator=yc[1]))

    # single-copy test
    problems = []
    for code in PANEL_SPECIES:
        n_cog = sum(1 for r in rows if r[-1] == code)
        yc = yutin_count(code)
        if n_cog > 1:
            problems.append(f"{code}: {n_cog} COG 2020 members")
        if yc and yc[0] > 1:
            problems.append(f"{code}: Yutin Table S3 count {yc[0]} (2nd paralog)")
    status = "excluded_multicopy" if problems else "built"

    # outputs ------------------------------------------------------------------------------------------------
    def w_tsv(name, cols, data):
        with open(HERE / name, "w", newline="") as fh:
            wr = csv.DictWriter(fh, cols, delimiter="\t", extrasaction="ignore", lineterminator="\n")
            wr.writeheader()
            wr.writerows(data)

    w_tsv("mapping_report.tsv", ["original_id", "species", "label", "gene_id", "status", "reason", "locator"],
          sorted(report, key=lambda r: (PANEL_SPECIES.index(r["species"]), r["original_id"])))
    (HERE / "unit_status.txt").write_text(status + ("\n" + "\n".join(problems) if problems else "") + "\n")
    if status != "built":
        for f in ("members.tsv", "reference_tree.nhx", "family.yaml", "assertions.tsv"):
            (HERE / f).unlink(missing_ok=True)
        print(f"{UNIT}: NOT BUILT ({'; '.join(problems)})")
        return

    members.sort(key=lambda m: (PANEL_SPECIES.index(m["species"]), m["gene_id"]))
    w_tsv("members.tsv", ["gene_id", "species", "original_id_in_paper", "label", "mapping_method", "mapping_identity",
                          "mapping_coverage", "source", "locator"], members)
    rid = f"{UNIT}_root"
    leaves = ",".join(m["gene_id"] for m in members)
    (HERE / "reference_tree.nhx").write_text(f"({leaves})[&&NHX:ev=S:id={rid}:tier=T2:src=Yutin2012];\n")

    complete = []
    for code in PANEL_SPECIES:
        yc = yutin_count(code)
        if yc and yc[0] >= 1 and have[code] == yc[0] and \
                all(r["status"] == "mapped" for r in report if r["species"] == code):
            complete.append(code)
    spp = sorted({m["species"] for m in members}, key=PANEL_SPECIES.index)
    arch_in = [s for s in spp if s in ARCHAEA]
    bact_in = [s for s in spp if s not in ARCHAEA]
    assertions = [
        dict(assertion_id=f"{UNIT}_A1", type="tag", target=rid, tier="T2", source_id="Yutin2012",
             locator=f"{Path(S1).name} (Table S1) row {u['row']}", quote=u["quote"][:300]),
        dict(assertion_id=f"{UNIT}_A2", type="tag", target=rid, tier="T2", source_id="Yutin2012",
             locator="fulltext.txt Methods (r-protein sets)", quote=Q_COG),
        dict(assertion_id=f"{UNIT}_A3", type="tag", target=rid, tier="T2", source_id="COG2020",
             locator="raw/Readme.2020-11-25.txt (cog-20.cog.csv)", quote=Q_COG2020),
    ]
    if u["section"] == "universal":
        assertions.append(dict(assertion_id=f"{UNIT}_A4", type="tag", target=rid, tier="T1", source_id="Yutin2012",
                               locator="fulltext.txt Abstract", quote=Q_UNIVERSAL))
    for i, code in enumerate(complete, 1):
        assertions.append(dict(assertion_id=f"{UNIT}_C{i}", type="completeness", target=code, tier="T2",
                               source_id="Yutin2012", locator=f"fulltext.txt Results (pipeline); {yutin_count(code)[1]}",
                               quote=Q_SCREEN))
    w_tsv("assertions.tsv", ["assertion_id", "type", "target", "tier", "source_id", "locator", "quote"], assertions)

    def files(lst):
        return [dict(path=p, sha256=sha256(ROOT / p)) for p in lst]

    n_rep = len(report)
    n_map = sum(1 for r in report if r["status"] == "mapped")
    meta = dict(
        id=HERE.name,
        name=f"Ribosomal protein {UNIT} ({cog}; Yutin 2012 {u['section']} r-protein)",
        panel=PANEL,
        collection="RPROT",
        phenomena=["CONTROL"] + (["DEEP"] if u["section"] == "universal" else []),
        descriptors=dict(n_members=len(members), n_species=len(spp), n_bacteria=len(bact_in), n_archaea=len(arch_in),
                         n_duplications=0, yutin_section=u["section"], cog=cog,
                         mapped_fraction=round(n_map / n_rep, 3) if n_rep else None),
        sources=[
            dict(id="Yutin2012", doi="10.1371/journal.pone.0036972", pmid="22615861",
                 title="Phylogenomics of prokaryotic ribosomal proteins", year=2012,
                 source_method="expert_classification", files=files([S1, S3, f"{YP}/fulltext.txt"])),
            dict(id="COG2020", doi="10.1093/nar/gkaa1018", pmid="33167031",
                 title="COG database update: focus on microbial diversity, model organisms, and widespread pathogens",
                 year=2021, source_method="expert_classification",
                 files=files([COGCSV, COGFA, f"{CP}/raw/Readme.2020-11-25.txt"]
                             + sorted(str(p.relative_to(ROOT)) for p in (ROOT / CP / "derived/seqmap").glob("*.tsv")))),
        ],
        root_level=None,
        root_single_copy=True,
        complete_species=complete,
        decoys=[],
        distinct_from=[],
        tier_summary=dict(T1=sum(a["tier"] == "T1" for a in assertions), T2=sum(a["tier"] == "T2" for a in assertions),
                          T3=0, events_D=0),
        protocol_version=1,
        notes=(f"Single-copy control unit of the RPROT collection. Members = COG 2020 {cog} members in the panel-K "
               f"genomes covered by COG 2020 (BRADU, GEOSL not covered), mapped by same-species sequence identity; "
               f"Yutin 2012 Table S1 assigns {cog} to r-protein {UNIT} ({u['section']}) and Table S3 gives one copy "
               f"per genome. One S polytomy (T2 rule b)."),
    )
    meta = {k: v for k, v in meta.items() if v is not None}
    (HERE / "family.yaml").write_text(yaml.safe_dump(meta, sort_keys=False, allow_unicode=True, width=120))
    print(f"{UNIT}: built {len(members)} members, {len(spp)} species, mapped {n_map}/{n_rep}, complete {len(complete)}")


if __name__ == "__main__":
    sys.exit(main())
