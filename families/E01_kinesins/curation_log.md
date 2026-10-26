# E01_kinesins curation log
status: curated (headline-eligible: 82.2% of in-panel rows uniquely mapped)

## Read
- Wickstead, Gull & Richards 2010 BMC Evol Biol 10:110 (PMID 20423470, OA): full text (literature/papers/PMID20423470/fulltext.txt),
  Additional file 3 (S3.XLS, 1,626 rows: Analysis ID, genome-project ID, GenBank acc, Family, Subfamily, exclusion flags),
  Additional file 7 (proteome sources/versions), Additional file 8 (DatasetS1.fas, motor-domain alignment, 1,458 seqs).

## Decisions
- Panel E species in the source (16 with members): HUMAN DROME NEMVE TRIAD MONBE YEAST NEUCR ENCCU DICDI ARATH CYAME THAPS
  PLAF7 TRYB2 NAEGR GIAIC; in source but 0 mapped: BATDJ (source strain JEL423, Broad BDEG_ ids), CHLRE (JGI v3 ids),
  ENTH1 (TIGR 2005 ids) - only motor-domain sequences, see below.
- Absent from source: AMPQE, MNELE, CAPO3, PHYRM, THETR, ICHMU. NB: source code 'Capsp' = Capitella sp. I (annelid,
  Add. file 7), NOT Capsaspora; Physo = Phytophthora sojae (not P. ramorum) -> not mapped.
- Mapping: (1) GenBank/RefSeq accession -> UniProt ID mapping (cache/idmap.tsv, UniProt 2026_03), restricted to the
  species' panel proteome; (2) exact gene-model name in the panel header (JGI DRAFT names e.g. TRIADDRAFT_<id>,
  NEMVEDRAFT_v1g<id>, MONBRDRAFT_<id>, THAPSDRAFT_<id>, NAEGRDRAFT_<id>, CYME_<id>, NCU, AGI, Tb ids, GL50803_ or the
  Analysis-ID gene symbol) with a single hit; (3) every ID mapping cross-checked against the paper's motor domain (DIAMOND,
  cache/motor_vs_panel.tsv): a >=98% motor hit to a different entry -> ambiguous. Rows with only a motor-domain sequence
  (~400 aa) cannot meet the >=80% target-coverage rule and are reported unmapped even at 100% local identity (most CHLRE,
  BATDJ, ENTH1, NAEGR losses). Excluded: 'ND' divergent rows, 'near_id' rows, duplicate rows hitting the same entry.
- Tree: U root; D node KIN_MAR_LECA (T1, lvl_lo=Eukaryota, lvl_hi=ABOVE_RANGE) over the 10 MAR families
  (Kinesin-1, 2, 3, 4/10, 5, 8, 9, 13, 14, 17); inside Kinesin-9, D node K9A_9B_LECA (9A and 9B both in the MAR).
  Non-MAR families, X groups, parenthesised labels and 'ungrouped' hang from the U root. Subfamily splits undated -> U.
  Terminal groups: <=1 gene/species -> S (T2 rule b; no two panel-E species share a WGD); else U (rule c).
- complete_species: genome-wide HMMER search quoted + >=90% of the species' in-panel rows mapped. Caveat: 2007 proteomes.
- Kinesin-myosin: homologous (P-loop motor fold); no distinct_from declared.

## Problems
- 1 ARATH row ambiguous; NAEGR 18/31 unmapped (JGI 'best proteins' ids not equal to UniProt DRAFT ids).
- Time: ~1.5 h.
