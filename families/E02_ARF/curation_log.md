# E02_ARF curation log
status: curated, protocol v1 (headline-eligible: 91.8% of in-panel rows uniquely mapped; 325 members, 22 species)

## Read
- Vargová et al. 2021, **Genome Biol Evol** 13:evab157 (not MBE), PMID 34247240, PMC8358228, OA CC BY.
  Fetched with `curation_tools.fetch_paper` -> literature/papers/PMID34247240/ (fulltext.txt, meta.json; supplement zip
  sha256 31e8c197… under raw/supp/). Used: Supplementary table 1 (2,286 genes in 115 species blocks; Annotation,
  Annotation basis, Conclusively annotated, Number of genes, Accession, Notes) and Supplementary dataset 1 (2,054 proteins).
- Text read: Results (LECA complement, SarB, Holozoa-specific paralogs, vertebrate expansions, ArfB, ArlRhodo), Methods.

## Decisions
- Panel-E species in source (22): HUMAN DROME NEMVE TRIAD MNELE AMPQE MONBE CAPO3 YEAST NEUCR BATDJ DICDI ENTH1 THETR
  GIAIC TRYB2 NAEGR ARATH CHLRE CYAME PLAF7 THAPS. Absent: ENCCU, PHYRM (source has P. sojae), ICHMU.
- Mapping: RefSeq/GenBank accession -> UniProt ID mapping (cache/idmap.tsv, release 2026_03) restricted to the species'
  panel proteome, cross-checked by DIAMOND of the paper protein; else same-species sequence mapping (>=98% id, >=80% cov
  both, RBH, margin >=2; seqmap.py). Conflicts -> ambiguous; same UniProt entry from 2 rows -> ambiguous.
  Excluded: HsaArf2 (pseudogene), author-flagged identical copies (both partners), rows with 'Number of genes' > 1.
- Tree (T1 events, all quoted in assertions.tsv):
  - ARF_LECA: multifurcating D, lvl_lo=Eukaryota / lvl_hi=ABOVE_RANGE, over the 16 LECA paralogs ("at least 15 … Arf1
    and 6, Arl1, 2, 3, 5, 6, 8, 13, 16, 17, and 18, Arfrp1, Sar1, and SarB" + SRβ). SarB has no conclusive panel member.
  - Arf1 child = U 'Arf1_complex': non-holozoan Arf1 | D at Holozoa (holozoan Arf1 vs Arf4; "holozoan stem itself
    (Arf4)") | TRIM23 (derived from Arf1, exon-intron quote).
  - Vertebrate duplications of Arf1, Arf4, Arl4 -> D nodes over the HUMAN copies (no lvl: level = HUMAN leaf branch,
    the only vertebrate in panel E; verifier does not accept species codes as lvl). Verifier warns "S polytomy with >1
    member in HUMAN" for their parent S nodes - expected, the HUMAN copies sit under the quoted D.
- T2: only rows with 'Conclusively annotated' = * are grouped; <=1 gene/species -> S polytomy, multi-copy species -> U
  (rules b/c; no two panel-E species share a WGD). Non-conclusive, 'unannotated', and groups whose source paralog the
  paper leaves unresolved (ArfB, Arf1/6, ArlRhodo, Arl4, Arl10, Arl11/14, Arl13L, Arl15, Arl19, Arl20) hang from the U root.
- Homology note: Arf/Sar and Rab are Ras-superfamily small GTPases (quote in assertions) -> no distinct_from with E03_Rab;
  deep merges between E02 and E03 must not be penalised (SYNTHESIS B6).
- complete_species (17): genome-wide search quoted ("searched in genome and/or transcriptome assemblies from 114
  eukaryotic species") + >=90% of the species' in-panel rows mapped. Caveat: 2018 assemblies.

## Problems
- MNELE 13/19 and AMPQE 13/20 mapped (new MNELE assembly UP001729289; AMPQE fused RefSeq models XP_0198567xx).
- NEMVE 17/23: JGI numeric IDs only, several paper models longer/shorter than UniProt entries (<80% coverage).
- Verifier warning A6: SRβ quote vs label 'SRb' (spelling only).

Time: ~1.2 h.
