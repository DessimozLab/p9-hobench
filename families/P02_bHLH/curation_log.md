# P02_bHLH — curation log

Status: **curated, headline-eligible** (82.5 % of in-panel paper members uniquely mapped). ~45 min.

## Source
Pires & Dolan 2010 MBE 27:862 (PMID 19942615, PMC2839125, OA). Read: abstract, Methods, all Results sections.
Files: File009 (566 rows species/protein/locus), File011 (FASTA, full-length), msp288_2.pdf = Suppl. Fig. S2 whose
text layer contains the table "Subfamily / Proteins" (544 plant proteins -> 26 subfamilies + Orphans; pypdf;
line-wrap artefacts 'bH LH' repaired; every protein assigned exactly once, checked in build.py).
source_method: phylogeny (ML/NJ/MP; subfamily support values in family.yaml descriptors).

## Mapping
Panel species: ARATH (158), ORYSJ (173), SELMO (103), PHYPA (98), CHLRE (3). Sequence first (exact or >=98 % id,
>=80 % cov both, RBH, >=2-point margin), then exact identifiers: AGI; LOC_Os -> RAP-DB gene by MSU7 locus
coordinates (rice.uga.edu all.locus_brief_info.7.0, sha256 in family.yaml; MSU7 pseudomolecules = IRGSP-1.0) +
unique same-strand >=50 % overlap with an Ensembl Plants RAP gene (cached) -> UniProt GN; PHYPADRAFT_/CHLREDRAFT_
-> UniProt GN (no hits: current UniProt PHYPA/CHLRE models use v3/v5 identifiers).
Rule added: when several paper proteins have the *identical* sequence and that sequence is exactly one panel
protein (8 Selaginella pairs = haplotype copies in the v1 assembly), the lowest-numbered is mapped and the others
are `excluded` (protocol 0.6, identical copies) — otherwise all would be "ambiguous" although the gene is certain.
Result: ARATH 100 %, ORYSJ 90 %, SELMO 94 %, PHYPA 32 % (v1.1 models vs v3 UniProt: coverage/identity failures),
CHLRE 1/3. complete_species ARATH, ORYSJ, SELMO (genome-wide sets; >=90 % mapped).

## Tree
- D_subfamilies_Embryophyta (lo Embryophyta, hi ABOVE_RANGE): the 20 subfamilies with moss members (computed from
  the S2 table; matches the paper's count of 20).
- D_subfamilies_Tracheophyta (lo Tracheophyta): Vb, X, XIII, XIV (lycophyte + angiosperm, no moss).
- D_subfamilies_angiosperm (lo Mesangiospermae): Ib(2), XV. Note: the text says "3 include only angiosperm
  proteins"; the S2 table gives 2 (III(d+e) has moss but no lycophyte members, so it counts among the 20).
- Upper bounds ABOVE_RANGE (inter-subfamily relationships unresolved; "basal nodes often have low support").
- Pairs across the three age groups: U (parent subfamily of younger subfamilies unknown); orphans and CHLRE: U;
  within subfamilies: U.
Earlier draft nested the three D nodes; rejected because hobench requires a D level >= the LCA of its genes, and
nesting would imply an unsupported order of subfamily origins.

## Notes
W:b warnings: T1 quotes give subfamily counts, not names (acceptable).
