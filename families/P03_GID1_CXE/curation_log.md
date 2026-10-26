# P03_GID1_CXE — curation log

Status: **curated, headline-eligible** (87.9 % of in-panel members uniquely mapped; 94 members, 9 species). ~35 min.

## Source
Yoshida et al. 2018 PNAS 115:E7844 (PMID 30068603, PMC6099883, OA; journal verified = PNAS). Read: Introduction,
Results (GID1 establishment; diversification in angiosperms), Discussion, Methods (sequence collection, phylogeny).
Files: Dataset S2 (sd02.xlsx; 259 named genes: CXEs of At/Os/Sm/Pp from Marshall et al. 2003 + 169 GID1s from 66
species), Dataset S1 (sd01; CXE/GID1 alignment, 4 species) and Dataset S3 (sd03; GID1 alignment). Ungapped
alignment rows are the query sequences. Fig. S1 (CXE clades I-VII) is an image only -> CXE clade structure not used.
source_method: phylogeny.

## Mapping
Panel species: ARATH, ORYSJ, SELMO, PHYPA (CXE + GID1), SOYBN, PHAVU, VITVI, AMBTC, MAIZE (GID1 only).
Sequence (exact / >=98 %, >=80 % cov, RBH, margin 2) then exact ids: AGI; SELMODRAFT_<JGI id>; Phvul.* ->
PHAVU_*g (UniProt GN); LOC_Os -> RAP (MSU7 locus overlap, see P02); XP_ via UniProt ID mapping.
Unmapped: 11 rice CXEs (MSU7 models with no RAP gene at the locus / divergent RAP models), PpCXE2 (90 % id),
maize GID1-2. Bug fixed during curation: query ids with spaces were truncated by DIAMOND (plib._qid).
complete_species: none (Dataset S2 is not a complete CXE inventory; GID1 sets are by reciprocal BLAST only).

## Tree
- D_GID1_origin (T1, lo Tracheophyta, hi ABOVE_RANGE): GID1 vs all tracheophyte CXEs listed (incl. the GID1-like
  clade-IV CXEs). Basis: "GID1 appeared after the divergence of vascular plants from the moss lineage" (no GID1 in
  Physcomitrella/Marchantia) and "GID1 was derived from one specific CXE group". The progenitor CXE lineage is not
  identified, so only "separate at Tracheophyta and below" is asserted; moss CXEs/GID1-likes stay under the root U.
- GID1: U[ Selaginella GID1-1/-2 (L-type; multicopy, monophyly not stated), S_Magnoliopsida[ Amborella (BA),
  S_Mesangiospermae[ monocot M-type (U: rice 1 + maize 1 mapped; maize copies 'might be caused by recent genome
  duplication'), D_GID1_AB ] ] ]. S nodes T2, quotes: Amborella and gymnosperms one copy; "default copy number in
  the early stage of seed plants might have been one".
- D_GID1_AB (T1, lvl Rosids): A/B split "just before or after the establishment of core eudicots"; panel P has no
  basal eudicot, so the whole admissible interval is the single panel branch into Rosids (VITVI has both types).
  The gamma hexaploidy lies on the same branch; the paper does not attribute the split to it -> no WGD tag.
- Within A/B types (Arabidopsis GID1a/c, soybean, Phaseolus copies): U.

## Phenomenon
Nested family: GID1 is a receptor lineage inside the CXE enzyme family; a tool must keep GID1 and CXEs in one HOG
only at Embryophyta or older, and separate them at Tracheophyta.
