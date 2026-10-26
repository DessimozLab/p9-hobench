# P01_homeobox — curation log

Status: **curated, NOT headline-eligible** (62.4 % of in-panel paper members uniquely mapped; protocol 0.5 needs 80 %).
Kept for secondary / class-boundary analyses (per K01 precedent: built units below 80 % stay in the directory).
Curator: Panel-P curation agent, 2026-10-26, ~1 h.

## Source
Mukherjee, Brocchieri & Bürglin 2009 MBE 26:2775 (PMID 19734295, PMC2775110, CC BY-NC, OA). Read: abstract, Methods,
Results (all class sections), Discussion (evolutionary history, rooting). Data: Supplementary table 1
(`literature/raw/B2/supp/PMC2775110/msp201_5.html`, 902 rows: name, accession, alternative names, species code,
class, protein sequence). Tree PDFs (msp201_1-4) not used (figures; T3fig not needed).
source_method: phylogeny (ML class trees + codomain architecture).

## Mapping
Species codes At/Os/Os_j/Zm/Pp/Sm/Cr/Ol (genome-wide searches) and Vv/Gm/Br/Pv (EST/GSS rows only).
Cascade: (1) same-species sequence (exact; or >=98 % id, >=80 % cov of both, RBH within the paper's species set,
>=2-point margin; DIAMOND --very-sensitive, hits cached in cache/diamond_*.tsv); (2) exact identifiers with a single
panel hit: SELMODRAFT_<Selmo1 id>, OSTLU_<id>, AGI, LOC_Os/RAP loci (gene2acc), GenBank/RefSeq protein accessions via
UniProt ID mapping (cache/uniprot_idmap.json). Several paper rows on one UniProt entry -> all ambiguous.
Result: 319 mapped / 511 in-panel non-excluded rows. ARATH 99 %, OSTLU 100 %, ORYSJ 85 %, SELMO 71 %, MAIZE 45 %,
PHYPA 17 %, CHLRE 0 %. Reason for the low rates: the paper's moss, maize, Chlamydomonas and part of the Selaginella
proteins are the authors' own homology-based re-predictions on 2007-era assemblies (Methods); current UniProt models
differ in exon structure (coverage < 80 % or identity 84-97 %), and those rows carry only scaffold coordinates (no
liftable IDs). Not fixable within the protocol.
Excluded (47): pseudogenes; rows whose label denotes a codomain protein without homeodomain (START ONLY, PLINC ONLY,
MEKHLA ONLY, SAWADEE ONLY, PH-START, RhoGAP-START, TM-START, START-HD-SAD, 'no HD'), SAWADEE-like, and CHLRE
'DDT_WSD motif' (Table 2 counts one Cr DDT gene).
complete_species: ARATH, OSTLU only (genome-wide search + >=90 % mapped).

## Tree (T1 events quoted in assertions.tsv)
- D_TALE_nonTALE (lo Viridiplantae, hi ABOVE_RANGE): TALE and non-TALE present in the plant/animal ancestor.
- D_KNOX_BEL (Viridiplantae..ABOVE): KNOX and BEL in red + green algae, present in the plant LCA.
- D_KNOXI_KNOXII (Embryophyta..ABOVE): both present in the moss/vascular ancestor. Only angiosperm genes carry the
  KNOX I/II sub-label; other KNOX genes sit in the KNOX U node.
- D_DDT_WOX_PINTOX (Viridiplantae..ABOVE): the three classes were present in the green-plant LCA.
- D_HDZIP_I_IV (Embryophyta..Streptophyta): HD-ZIP classes absent from Chlorophyta/Rhodophyta, present before the
  moss/vascular split ("possibly before transition to land"); upper bound = the streptophyte stem (paper's
  lineage-specific expansion interpretation, stated as more plausible than early duplication + loss).
- S nodes at Viridiplantae for DDT and WOX (chlorophyte gene vs land-plant genes), T2, quote "each of these classes was
  originally represented by one gene". WOX loss in Chlamydomonas tagged (quote).
- U: inter-class relations among PLINC/PHD/NDX/LD/SAWADEE/HD-ZIP and DDT/WOX/PINTOX ("the trees do not provide
  reliable information on the evolutionary relations among most of the classes"); within classes (paralogs, WGD
  lineages, no synteny curation: protocol 3 rule 2b/2c); unclassified TALE genes.
No WGD event asserted (the paper does not date within-class duplications).

## Problems / notes
- Verifier warnings W:b on 3 T1 quotes (they name classes by group, e.g. "TALE and non-TALE"), acceptable.
- `plib.py` is a copy of families/_shared_P/plib.py so that `verify_provenance --rebuild` works in isolation.
