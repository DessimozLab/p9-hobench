# E07b_adaptins_medium_small curation log
status: curated, protocol v1 (mapped fraction E07a 0.72, E07b 0.72, E07c 0.81; E07a/E07b below 80% -> not headline)

## Read
- Hirst et al. 2014 eLife 3:e02866 (PMID 24867644, OA) via fetch_paper: fulltext.txt, Figure 3-source data 1 (1,413 rows:
  section = subunit, annotation, accession, length), Figure 1-source data 1/2 (reverse-HHpred hits; beta/gamma-COP among
  large, delta/zeta-COP among medium/small homologues). supp1.zip alignments are masked (not usable for mapping).

## Decisions
- Split into three families by subunit homology class (large; medium+small; COPI/TSET scaffold): classes are not asserted
  homologous to each other ("The large adaptor subunits share sequence and structure homology, as do the medium and small
  subunits"). Same build.py in all three dirs (family chosen by directory name); prep_mapping.py run from E07a only.
- Tree: T1 D (lvl_lo=Eukaryota, lvl_hi=ABOVE_RANGE) over subunit groups: AP-1..5 split pre-LECA (quote), AP/COPI common
  ancestry pre-LECA, TSET in LECA and distinct, large-type and small/medium-type duplications pre-LECA (Fig. 4 legend),
  TTRAY1/2 and alpha/beta'-COP duplications pre-LECA. AP-1/AP-2 beta (1B/2B labels) pooled into a U group.
  Groups: <=1 gene/species -> S, else U. Unparseable annotation -> loose under U root.
- Mapping: RefSeq/GenBank/UniProt/SGD -> UniProt ID mapping; gene-model names in panel header maps (AGI, DDB_G, EHI_,
  PF3D7_, Tb927, AMSG_, CAOG_ zero-padded, Cre->CHLRE_..v5, NCU, CYME_, GL50803_, JGI protein id -> <ORG>DRAFT_<id>) with
  a +-20% length check against the table's Length column; else NCBI sequence of the RefSeq id with an exact same-species
  match (or unique >=99.5%, cov >=0.8 both). Giardia GL50581 (assemblage B) = not_in_panel.
- complete_species empty: table lists sequences used for phylogenies, AP/COPI sampling not exhaustive.

## Problems
- Unmapped: Broad BDEG_ (BATDJ strain JEL423), old PlasmoDB ids (PF13_0062...), Tb11 / some Tb927 ids not in headers,
  JGI EST/model names (estExt_..., e_gw...), JGI v3/v4 numeric CHLRE ids; DICDI TSET genes absent from panel proteome.
