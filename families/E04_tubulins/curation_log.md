# E04_tubulins curation log
status: curated (headline-eligible: 86.3% of in-panel source sequences uniquely mapped)

## Read
- Findeisen et al. 2014 GBE 6:2274 (PMID 25169981, OA): fulltext.txt, Table 1, supplementary zip
  (tubulin_alignment.fas: 3,286 full-length named tubulins; DatasetS2.xls: 504 species + abbreviations;
  SupplementaryMaterial.pdf: Suppl. Table S3 human tubulin loci used only to sanity-check Tub4/Tub5 = TUBD1/TUBE1).

## Decisions
- Names = <abbrev>Tub<class><copy>; class 1-6 = alpha..zeta (paper: numbers replace Greek letters; Tub1=alpha, Tub3=gamma,
  Tub2=beta quoted; 4/5/6 by Greek order, consistent with Table 1 counts). Only the class split is used in the tree.
- 24/25 panel-E eukaryotes in the source (MNELE absent). Strain choice: YEAST=Sc_c (S288c), BATDJ=Bad_b (JAM81),
  PLAF7=Pf_a (3D7). GIAIC source strain ATCC 50581 (assemblage B) differs from panel WB -> 1/2 mapped.
- Mapping (same species only): exact full-length sequence identity to a unique panel entry (exact_sequence), else DIAMOND
  >=98% id, >=80% coverage both, RBH, >=2-point margin. Identical protein copies (in source or proteome) excluded (protocol 0.6),
  one representative kept when a single panel entry exists. Near-identical alpha/beta paralogs that are not exact stay ambiguous.
- Tree: root D TUB_LECA (T1, lvl_lo=Eukaryota, lvl_hi=ABOVE_RANGE) over Tub1..Tub6 (quote: ancestor of each subfamily in
  LCA of eukaryotes). Classes with >1 gene in any species -> U (alpha, beta, gamma); delta, epsilon, zeta -> S polytomy (T2 rule b).
- FtsZ/TubZ/artubulins/BtubA-B were outgroups in the paper; not members; no distinct_from (no curated FtsZ family; quote kept in notes).
- complete_species: TBLASTN of whole genome assemblies (quoted) + >=90% mapped.

## Problems
- alpha/beta/gamma give no S assertions (multi-copy species -> rule c). Time ~45 min.
