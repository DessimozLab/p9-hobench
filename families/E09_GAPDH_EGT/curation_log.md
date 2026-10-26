# E09_GAPDH_EGT curation log

Status: **abandoned** (fewer than 4 panel species with source-listed members for the EGT event). Time: ~40 min.

## Read
- Figge et al. 1999 MBE 16:429 (requested/49.pdf, sha256 in requested/RECEIVED.tsv): gap1 ~ eukaryotic GapC, gap2 ~
  GapAB; quote: "eukaryotes acquired their nuclear genes for GapC and GapAB via endosymbiotic gene transfer from the
  antecedents of mitochondria and chloroplasts". Members only as figure leaves / GenBank ids of non-panel bacteria.
- Stanley et al.? 2013 PLoS One (PMID 23936198, CC BY): Data S1 = 351 trimmed GAPDH sequences (~325 aa) + FastTree.
  Panel species in the green group (Gap2 + GapA/B): SYNY3 (Gap2), ARATH (GapA x2, GapB), CHLRE (GapA). No
  Cyanidioschyzon (red algae only Galdieria GapA). Quote (Fig. 2): "Gap2 of cyanobacteria became the GapA of the
  Archaeplastida, with a further gene duplication event in land plants giving rise to GapB."
- PMID 42482530 (2026 New Phytol, CC BY): GapB present in the LCA of Viridiplantae (contests Petersen 2006
  Streptophyta origin); Table S2 Chlorophyta only.
- Martin & Cerff 2017 Protoplasma (PMID 28265765): review, no member ids.

## Why abandoned
- The quotable EGT (cyanobacterial Gap2 -> Archaeplastida GapA/B, H node donor=SYNY3) has source-listed members in
  only 3 panel species (SYNY3, ARATH, CHLRE). Red algal GapAB is stated (Figge 1999) but no source lists a
  C. merolae sequence; THAPS uses GapC1 (plastid-retargeted cytosolic copy), not GapA.
- Extending to cytosolic GapC (HUMAN, YEAST, DICDI, ...) would add only unasserted (U) members: the GapC
  mitochondrial-donor hypothesis has explicit exceptions (Euglena/trypanosome GapC), no panel donor member is
  source-listed (no BRADU/ECOLI gap1 in any source), and eukaryotic GapC copy numbers (GapC1/2, GapCp, TDH1-3,
  GAPDHS) are not resolved in an OA source.
- A minimal 3-species version (H[donor=SYNY3](Gap2, GapA/B clade); D GapA/GapB lvl_lo=ARATH, lvl_hi=Viridiplantae;
  CHLRE GapA under U) is feasible from PMID 23936198 Data S1 (domain-rule mapping, v1.1) if the PI waives the
  4-species rule.
