# Y02_MH2015_ohnologs — curation log

Status: **abandoned as a standalone family** (kept as a documented cross-check of Y01). Time ~30 min.

Source: Marcet-Houben & Gabaldón 2015, PLoS Biol 13:e1002220 (PMID26252497, PMC4529251, CC BY), fetched with
curation_tools (OA XML + 24 supplements). S1 Data `pbio.1002220.s001.xls` (sha256 7830b3b7…), sheet `Fig1C`: 365 rows =
245 unique *S. cerevisiae* ohnolog pairs with an "Inferred topology" label ('A – A' 160 pairs, 'A - AB' 81, both 4).

Why no family:
1. Only *S. cerevisiae* identifiers are given (1 panel species; protocol needs ≥4). The per-gene reduced-phylome trees
   (seed + one ZT + one KLE species + outgroups) were linked from genome.crg.es/~mmarcet/yeast_hybrids/phylome_table.htm;
   the page is archived (2015-08-15) but the tree files are not, and the live site redirects.
2. The ohnolog pairs are not independent of Y01: the paper restricts to "conserved pairs of WGD-ohnologs as defined in
   the Yeast Gene Order Browser (YGOB) database" (Results).
3. The topology labels 'A – A' / 'A - AB' are not defined verbatim in the text or legends. From the Fig 1C caption
   (two populations: duplication mapped at the WGD branch vs at the pre-KLE branch) and the Fig 4A definitions (A = with
   ZT, B = with KLE, C = sister to both) the likely reading is A–A = both copies with ZT, A–AB = one copy outside ZT+KLE;
   but this is an inference, and the label does not say *which* copy is which, so no per-copy ZT placement can be
   asserted (protocol §0.2 ⇒ ZT vs post-WGD stays U in Y01).

Cross-check with Y01 (`build.py` → `crosscheck.tsv`):
- 245/245 MH2015 pairs are S. cerevisiae ohnolog pairs in the complete YGOB v1 Pillars.tab.
- 135/245 fall in the truncated v7 pillars used by Y01: **135/135 are encoded as the same YEAST ohnolog pair (D[WGD])
  in Y01; 0 discordant.** The other 110 are outside the truncated part of the v7 file.
- Caveat for the A2 encoding: 85/245 pairs (35 %) carry 'A - AB', i.e. MH2015's trees place the split of these
  ohnologs *before* the KLE divergence ("more ancestral duplications are prominent with >50% of the duplications being
  mapped to the two nodes preceding the WGD (pre-ZT and pre-KLE)"). Y01's interval [ZT+WGD, Saccharomycetaceae]
  scores a tool that splits such pairs above Saccharomycetaceae as wrong. Consider widening lvl_hi to Saccharomycetaceae
  (merge up to its parent) if the PI accepts MH2015's pre-KLE signal as real rather than an artefact.
