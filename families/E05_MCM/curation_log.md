# E05_MCM curation log

Status: **curated** (protocol v1/v1.1). Curator: agent (panel E batch: E05, E09, E08). Time: ~1.5 h.

## Read
- Liu, Richards & Aves 2009 BMC Evol Biol 9:60 (PMID 19292915, PMC2667178, CC BY): full text (Europe PMC XML ->
  `literature/papers/PMID19292915/fulltext.txt`), Additional files 1-4 (vector PDF trees with one accession per leaf),
  5/7 (Newick, labels only, no accessions), 6 (genome sources, .doc).
- Survey entries B1 C6 / B2 E7.

## Extraction (extract.py)
- Additional files 1-4 captions give the paralogues per figure (S1 Mcm4,Mcm7; S2 Mcm5,Mcm6; S3 Mcm3,Mcm2; S4 Mcm8,
  Mcm9, archaea). Leaves are read with pypdf page coordinates, sorted top-down, split into blocks by vertical gaps
  (>12 pt), archaea split by the species list in Methods, S2 split where its first species re-occurs. The paralogue
  label printed beside each clade must fall within its block (asserted). Human accessions confirm the order
  (NP_005905 MCM4, NP_005907 MCM7, NP_006730 MCM5, NP_005906 MCM6, NP_002379 MCM3, NP_004517 MCM2, NP_877954 MCM8,
  NP_694987 MCM9).
- Additional file 7 (Newick) was used only as a cross-check of species x paralogue presence.

## Mapping
- NCBI accessions -> efetch sequence -> same-species DIAMOND (>=98% id, >=80% both cov, RBH, margin >=2).
- JGI (Naegr1, Triad1, Nemve1, Monbr1, Phyra1_1, Thaps3), CyanidioschyzonDB, GiardiaDB ids -> exact unique ORF-name
  match in panel header_map (`id_mapping_orf_name`). Chlre3 ids -> NCBI old locus tag CHLREDRAFT_<n> -> v5 record
  sequence -> DIAMOND.
- Human MCM9: paper used the 391-aa isoform 2 -> UniProt ID mapping gives gene MCM9 -> unique panel entry Q9NXL9
  (gene_name_unique; 98.5% identity over 100% of query).
- QC: every mapped target must hit the paper's human MCMs at e <= 1e-20 (paper's BLAST cut-off); membership only.
- Result: 113 mapped / 137 in-panel leaves (82.5%); 18 unmapped (THAPS 8, PHYRM 6: Thaps3/Phyra1_1 model ids absent
  from UniProt; NEMVE 2; ARATH Mcm8 92% id; ENTH1 Mcm8 wrong accession), 6 ambiguous (accession printed under two
  paralogues: NEUCR XP_965670, TRYB2 XP_951609, MONBE Monbr1:32709), 3 excluded (YEAST NP_010883 = SPF1, failed QC;
  2 human isoform duplicates), 8 not_in_panel (BATDJ: paper used strain JEL423 Broad BDEG_ loci, panel is JAM81;
  no sequences retrievable).

## Tree / decisions
- Root S (T1, Archaea_Eukaryota): Halobacterium MCM (Q9HNA5) vs eukaryotic set; paper: single MCM per archaeon
  (except M. jannaschii) and Mcm2-9 derived from a single archaeal-like MCM (hedged "may be derived").
- MCM_LECA_DUPS: one multifurcating D (T1) for the seven pre-LECA duplications, lvl_lo=Eukaryota,
  lvl_hi=Asgard_Eukaryota (sampled Euryarchaeota + TACK have one MCM; Asgard unsampled). Order of duplications is
  unresolved in the paper, so one multifurcation; Mcm8/9 sister grouping only 68% BS -> not used.
- Per paralogue S polytomy (T2 rule b; panel E has no WGD shared by two species).
- Naegleria Mcm3 x2: D (T1, "recent MCM duplication", LCA mapping NAEGR; verifier rejects species as lvl).
- Conflict: Table 1 marks Nematostella Mcm8/9 absent, but Additional file 4/7 and the Methods animal list include
  it; figure labels used (NEMVE Mcm8 unmapped anyway; Mcm9 A7S1E9 mapped).
- complete_species: >=90% mapped and no QC-excluded leaf (YEAST lost MCM3 to a source error -> not complete).
- Validation: hobench validate_family OK; verify_provenance --rebuild PASS (3 warnings: two T1 quotes without gene
  labels; S polytomy Mcm3_LECA has a NAEGR D child).
