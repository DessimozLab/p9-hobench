# V04_Sox curation log (protocol v1.1)
Status: **validated, NOT headline** (46% of in-panel S6 entries mapped; 94% of gene-level identifiers).

Read: Voldoire et al. 2017 PLoS ONE (PMID 28738066, CC BY) full text + S6 Table (pone.0180936.s006.xls,
literature/raw/B1/supp/sox_Voldoire2017__PMC5524304/, sha256 in family.yaml).

Decisions
- Members = every S6 cell for a panel species (CALMI, HUMAN, MOUSE, LATCH, LEPOC, DANRE, ORYLA, SALSA, ONCMY).
  Mapping (shared helper vlib.py): Ensembl/RefSeq ID -> UniProt ID mapping -> panel proteome; spotted-gar S6 entries are
  LepOcu1 coordinates -> Ensembl REST overlap (exactly one protein-coding gene) -> UniProt; zebrafish (RefSeq-based
  UniProt proteome, no Ensembl xrefs) -> same-species sequence mapping of the Ensembl canonical translation / UniParc
  versions of retired ENSDARP IDs (>=98% id, >=80% cov both, 2-point margin).
- Unmapped by necessity: all elephant-shark (AAVX01 WGS contigs), salmon and trout (pre-2016 scaffold IDs; one salmon
  scaffold is given for 6 different genes) and TIGR-TC/EST entries -> listed in mapping_report.tsv. Ss4R copies are
  therefore absent from the tree.
- Tree: group splits (B1,B2,C,D,E,F,H) at panel root, lvl_hi=ABOVE_RANGE ("present in all bilaterians");
  within-group paralogs = 2R as interval Gnathostomata..Vertebrata (paper: "probable consequence of the two rounds");
  teleost a/b = 3R at Teleostei (T1: 10 genes named in text); sox17/sox32 = SSD at Teleostei (T1, S2 Fig text);
  same-label genes = S polytomies (synteny-curated per text, "same copy kept").
- complete_species: HUMAN, MOUSE, LATCH, DANRE (genome-wide search stated; >=90% mapped).
Problems: S6 SOXE header says 'sox8/sox8a' but non-teleost sox8 sits in the 'sox8b' column (irrelevant to the tree).
Mouse Sox2 (97.8% id of Ensembl canonical vs P48432) and 2 zebrafish/2 medaka retired IDs remain unmapped.
Time: ~1.5 h.
