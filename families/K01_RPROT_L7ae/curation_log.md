# Curation log: K01_RPROT_L7ae (collection RPROT, panel K)

Status: **abandoned (unit not single-copy / not mappable; see table)**. Curator: K-panel curation agent, 2026-10-26. Time for the whole RPROT collection: ~1.5 h.

## Read
- Yutin et al. 2012 PLoS One (PMID 22615861, PMC3353972, CC0): full text (Abstract, Results "Paralogous r-proteins in
  Bacteria", "Paralogy ... archaeal", Methods), Table S1 (r-protein list, COG ids, universal/bacterial-only/archaeal-only),
  Table S2 (genomes), Table S3 (copy-number matrix: main set + 2nd paralogs). Files under literature/papers/PMID22615861/.
- NCBI COG 2020 (Galperin et al. 2021 NAR, PMID 33167031): cog-20.cog.csv, cog-20.fa.gz, cog-20.org.csv (md5 verified
  against checksums.md5.txt); panel-K subset extracted by literature/papers/COG2020/extract_panelK.py and cross-referenced
  onto the UniProt 2026_03 panel proteomes by literature/papers/COG2020/seqmap_panelK.py (curation_tools.map_by_sequence,
  ambiguity delta 2, coverage 0.8). sha256 in literature/papers/COG2020/MANIFEST.sha256.

## Decisions
- Yutin 2012 gives no per-genome protein identifiers, only COG ids (Table S1) and copy numbers (Table S3). Members are
  therefore the COG 2020 members of Yutin's COG for this r-protein in panel genomes; Yutin Table S3 is used as an
  independent copy-number check and for complete_species. COG 2020 is an expert-curated classification (not on the
  denylist); its genome set covers 19/21 panel species: BRADU (B. diazoefficiens USDA 110) and GEOSL (G. sulfurreducens)
  are missing (COG only has B. japonicum USDA 6 and G. metallireducens = other species; never mapped across species).
  Their r-proteins appear in mapping_report.tsv as 'unmapped' (present per Yutin Table S3, no identifier in any source).
  SYNY3 is in COG 2020 but not in Yutin's genome set (no completeness claim).
- Strain identity: ECOLI K-12 MG1655 (COG/Yutin) = panel K-12; HALSA = Halobacterium sp. NRC-1 (taxid 64091) in both;
  SYNY3 COG taxid 1148 = PCC 6803 (panel substrain Kazusa 1111708).
- Mapping: exact sequence (sha256) or unique DIAMOND hit >=98 % identity, >=80 % coverage of both, 2-point margin.
  Members whose COG and UniProt gene models differ in length (coverage < 0.8) are 'unmapped' even at 100 % identity.
  Reciprocal-best-hit not computed separately: all accepted non-exact hits are >=98.9 % identity, unique, single-copy units.
- COG membership class 3 (partial) accepted only when it is the only COG member in that genome AND Yutin Table S3 has
  exactly one copy there; otherwise excluded.
- Bacterial-only units (Table S1 'Bacterial-only proteins (23)'): archaeal COG members excluded (Yutin lists the archaeal
  r-proteins separately; e.g. COG0197 = L16/L10E holds archaeal L10e, listed by Yutin as archaeal-only L10e).
- Single-copy rule (PI brief): a unit is built only if every panel species has <=1 COG member AND Yutin Table S3 lists no
  2nd paralog for that genome. Paralogous units (zinc-ribbon L31/L33/L36/L28/S14/S18 variants, S7/S12 in AQUAE, S10 in
  HALSA, S17, L7ae, S21 in BRADU/GEOSL, S11 with 2 COG members in METAC) are not built. Yutin text: "all highly diverged
  paralogs are r-proteins containing various forms of the Zn-binding motif".
- Tree: one S polytomy (T2 rule b), root_single_copy true; universal units additionally tagged DEEP with the Abstract quote
  (34 universal r-proteins). Bacterial vs archaeal members of a universal unit are asserted orthologous at the panel root
  (LUCA), per Yutin's classification of the protein as universal. Risk noted: HGT of r-proteins (e.g. S4 zinc/non-zinc
  forms, Chen 2009, cited by Yutin) is not resolved by this classification.
- Excluded by curator: L7L12 (Yutin's COG2058 has no bacterial members in COG 2020; bacterial L7/L12 is elsewhere), S1
  (excluded by Yutin). S22 (PRK10057) and S31e (no COG) and archaeal-only r-proteins (arCOG ids only) were not attempted.
- Levels: no lvl tags are needed (S polytomy only); root_level left to the tree's LCA mapping.

## Collection overview (all units attempted)
| unit | status | detail |
|---|---|---|
| L1 | built | 19 members / 19 spp; mapped 19/21 |
| L10 | built | 19 members / 19 spp; mapped 19/21 |
| L11 | built | 19 members / 19 spp; mapped 19/21 |
| L13 | built | 19 members / 19 spp; mapped 19/21 |
| L14 | built | 19 members / 19 spp; mapped 19/21 |
| L15 | built | 19 members / 19 spp; mapped 19/21 |
| L16 | built | 12 members / 12 spp; mapped 12/14 |
| L17 | built | 12 members / 12 spp; mapped 12/14 |
| L18 | built | 19 members / 19 spp; mapped 19/21 |
| L19 | built | 12 members / 12 spp; mapped 12/14 |
| L2 | built | 19 members / 19 spp; mapped 19/21 |
| L20 | built | 12 members / 12 spp; mapped 12/14 |
| L21 | built | 11 members / 11 spp; mapped 11/14 (<80%: not headline) |
| L22 | built | 18 members / 18 spp; mapped 18/21 |
| L23 | built | 19 members / 19 spp; mapped 19/21 |
| L24 | built | 18 members / 18 spp; mapped 18/21 |
| L25 | built | 12 members / 12 spp; mapped 12/14 |
| L27 | built | 12 members / 12 spp; mapped 12/14 |
| L28 | excluded_multicopy | MYCTU: 3 COG 2020 members; MYCTU: Yutin Table S3 count 3 (2nd paralog) |
| L29 | built | 16 members / 16 spp; mapped 16/21 (<80%: not headline) |
| L3 | built | 19 members / 19 spp; mapped 19/21 |
| L30 | built | 16 members / 16 spp; mapped 16/18 |
| L31 | excluded_multicopy | ECOLI: 2 COG 2020 members; ECOLI: Yutin Table S3 count 2 (2nd paralog); PSEAE: 2 COG 2020 members; PSEAE: Yutin Table S3 count 2 (2nd paralog); BACSU: 2 COG 2020 members; BACSU: Yutin Table S3 count 2 (2nd paralog) |
| L32 | built | 11 members / 11 spp; mapped 11/14 (<80%: not headline) |
| L33 | excluded_multicopy | BACSU: 2 COG 2020 members; BACSU: Yutin Table S3 count 3 (2nd paralog); MYCTU: 2 COG 2020 members; MYCTU: Yutin Table S3 count 2 (2nd paralog) |
| L34 | built | 12 members / 12 spp; mapped 12/14 |
| L35 | built | 12 members / 12 spp; mapped 12/14 |
| L36 | excluded_multicopy | ECOLI: 2 COG 2020 members; ECOLI: Yutin Table S3 count 2 (2nd paralog); PSEAE: 2 COG 2020 members; PSEAE: Yutin Table S3 count 2 (2nd paralog) |
| L4 | built | 19 members / 19 spp; mapped 19/21 |
| L5 | built | 19 members / 19 spp; mapped 19/21 |
| L6 | built | 19 members / 19 spp; mapped 19/21 |
| L7L12 | excluded_by_curator | Table S1 assigns COG2058 (archaeal L12E/P1/P2 family) to bacterial L7/L12; in COG 2020 COG2058 has no bacterial member in the panel genomes (Yutin Table S3 lists L7L12 in every bacterium), so the COG cannot serve as the membership source for this universal unit |
| L7ae | excluded_multicopy | BACSU: 2 COG 2020 members; BACSU: Yutin Table S3 count 2 (2nd paralog) |
| L9 | built | 12 members / 12 spp; mapped 12/14 |
| S1 | excluded_by_curator | Yutin 2012 excluded S1p from the study (varied domain architecture, ubiquitous S1-like domains) |
| S10 | excluded_multicopy | HALSA: 2 COG 2020 members; HALSA: Yutin Table S3 count 2 (2nd paralog) |
| S11 | excluded_multicopy | METAC: 2 COG 2020 members |
| S12 | excluded_multicopy | AQUAE: 2 COG 2020 members; AQUAE: Yutin Table S3 count 2 (2nd paralog) |
| S13 | built | 19 members / 19 spp; mapped 19/21 |
| S14 | excluded_multicopy | BACSU: 2 COG 2020 members; BACSU: Yutin Table S3 count 2 (2nd paralog); MYCTU: 2 COG 2020 members; MYCTU: Yutin Table S3 count 2 (2nd paralog); DEIRA: 2 COG 2020 members |
| S15 | built | 19 members / 19 spp; mapped 19/21 |
| S16 | built | 11 members / 11 spp; mapped 11/14 (<80%: not headline) |
| S17 | excluded_multicopy | BACTN: 2 COG 2020 members; BACTN: Yutin Table S3 count 2 (2nd paralog); METAC: 2 COG 2020 members; METAC: Yutin Table S3 count 2 (2nd paralog) |
| S18 | excluded_multicopy | MYCTU: 2 COG 2020 members; MYCTU: Yutin Table S3 count 2 (2nd paralog) |
| S19 | built | 19 members / 19 spp; mapped 19/21 |
| S2 | built | 19 members / 19 spp; mapped 19/21 |
| S20 | built | 12 members / 12 spp; mapped 12/14 |
| S21 | excluded_multicopy | BRADU: Yutin Table S3 count 2 (2nd paralog); GEOSL: Yutin Table S3 count 2 (2nd paralog) |
| S3 | built | 19 members / 19 spp; mapped 19/21 |
| S4 | built | 19 members / 19 spp; mapped 19/21 |
| S5 | built | 19 members / 19 spp; mapped 19/21 |
| S6 | built | 12 members / 12 spp; mapped 12/14 |
| S7 | excluded_multicopy | AQUAE: 2 COG 2020 members; AQUAE: Yutin Table S3 count 2 (2nd paralog) |
| S8 | built | 19 members / 19 spp; mapped 19/21 |
| S9 | built | 19 members / 19 spp; mapped 19/21 |

## Problems
- Units with <80 % uniquely mapped in-panel members (L21, L29, L32, S16) stay in the directory but fail the headline
  criterion (gene-model length differences between COG/RefSeq and UniProt, or no COG member although Yutin lists one).
