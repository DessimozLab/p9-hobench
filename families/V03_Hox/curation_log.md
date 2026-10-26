# V03_Hox curation log (protocol v1.1)
Status: **validated, headline-eligible** (98% of nomenclature entries uniquely mapped; 5 species).

Read: Holland et al. 2007 (PMID 17963489; Table 2, gene-family definition), Kumar et al. 2024 review (PMID 39079985;
2R -> HoxA-D, teleost 7-8 clusters by TS-3R), Pascual-Anaya et al. 2012 (PMID 23198682; posterior PG orthology to
invertebrates unclear), Session et al. 2016 (PMID 27762356; L/S subgenomes), Kikuchi et al. 2025 (PMID 41268905;
49 zebrafish hox genes), Braasch et al. 2016 gar genome (PMID 26950095; text + Supplementary Table 23).

Decisions
- Members from nomenclature authorities (downloaded 2026-10-26, sha256 in family.yaml): HGNC (HUMAN), MGI (MOUSE),
  ZFIN (DANRE; hox a/b = 3R clusters), Xenbase (XENTR; XENLA .L/.S). Mapping via the authorities' UniProt
  cross-references (unique accession in the panel proteome), else NCBI GeneID -> UniProt.
- NOT used: Braasch 2016 gar/medaka-teleost ohnolog tables (Suppl. Table 23 / Source Data 4) - their orthologs come
  from "phylogenetic [Vilella 2009 EnsemblCompara GeneTrees] and conserved synteny analyses" -> Compara is denylisted.
  Gar Hox (Suppl. Fig. 12) is a figure only. Medaka/salmon/shark/coelacanth/lamprey/amphioxus: no OA machine-readable
  expert Hox tables found (Amores 1998, Kuraku & Meyer 2009, Kurosawa 2006, Mungpakdee 2008 not OA).
- Tree: root = Holland gene families Hox1, Hox2, Hox3, Hox4, Hox5, Hox6-8, Hox9-13 (single genes in the bilaterian
  ancestor -> lvl_lo=Chordata, lvl_hi=ABOVE_RANGE); PG splits inside Hox6-8 / Hox9-13 = pre-2R tandem duplications
  (T2, lvl_lo=Vertebrata, lvl_hi=ABOVE_RANGE); clusters A-D per PG = 2R (interval Gnathostomata..Vertebrata);
  zebrafish a/b = 3R (Teleostei); X. laevis .L/.S = allotetraploidy (D on the XENLA branch, tag ALLO); same symbol
  across species = S. PG rule cross-checked against Holland Table 2 (assert in build.py).
- complete_species: HUMAN (Holland: 39 HOX genes), DANRE (Kikuchi: 49 hox genes); all mapped.
Problems: 4 XENLA genes (hoxb2.L, hoxb8.L, hoxc8.L, hoxd1.S) unmapped/ambiguous. All 2R/PG events lie above the members'
LCA (Euteleostomi) -> cross-cluster pairs score as 'never in one HOG within range'.
Time: ~1.5 h.
