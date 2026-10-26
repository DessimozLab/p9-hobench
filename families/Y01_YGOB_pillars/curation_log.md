# Y01_YGOB_pillars — curation log (collection YGOB, panel Y)

Status: **v1 built and validated (0 errors / 0 warnings on all units)** against a *provisional* panel-Y species tree
(panels/Y not yet published; level roles are derived by LCA in build.py, so a rebuild with
`--species-tree panels/Y/species_tree.nwk` only renames levels). Rebuild: `python families/Y01_YGOB_pillars/build.py`.

## Sources read (sha256 in collection.yaml)
- YGOB v7 (Aug 2012) `Pillars.tab`, Wayback 2021-07-28 = 2025-02-14 capture, both truncated by the archive at 1,048,576 B
  (sha256 e8e12eda…). 5,173 complete pillar lines used; the partial line 5,174 (`Kpol_1032.86\tTPHA0N01030\t`) dropped.
  Order of the 5,173 lines is the file's, i.e. roughly the first third of the ancestral genome, not a random sample.
- YGOB v1 (Sep 2005; the release described in Byrne & Wolfe 2005) from Wayback, complete: README,
  `Synteny_{Klactis,Agossypii,Kwaltii,Skluyveri}.tab` (Track A/B per outgroup locus with 6-char synteny score),
  Pillars.tab, classes.tab, genome tabs (Proteins.fsa truncated, unused). Saved in literature/raw/B1/ygob/wayback/.
- Byrne & Wolfe 2005 (PMID16169922, PMC HTML text), Scannell et al. 2007 (PMID17494770, PMC HTML text; supplements
  behind PMC proof-of-work wall), Marcet-Houben & Gabaldón 2015 (PMID26252497, OA XML).

## Key facts / decisions
1. **v7 Pillars.tab carries no track information**: README "The order of the columns and the genes is ARBITRARY, and does
   NOT indicate presence on an A or B track"; paper: "the two slots for genes from a post-WGD species are not
   preassigned to particular tracks". Hence cross-species post-WGD pairs from v7 alone are U (pseudo-orthology risk,
   Scannell 2007: "Almost half of their surviving single-copy genes are not orthologs but paralogs").
2. **Tracks from YGOB v1 synteny tables** (README rule: same track ⇒ orthologs if score robust 1/0; A vs B ⇒ paralogs)
   resolve YEAST × CANGA only (v1 has S. cer, C. gla, S. castellii; S. castellii `Scas_` contigs of Cliften 2003 cannot
   be mapped to NAUCC). Transferred to a v7 pillar only if (a) the v1 row is anchored on the same KLLA/AGOS gene,
   (b) its YEAST+CANGA gene set equals the v7 pillar's, (c) Scer A/B and Cgl A/B score positions are all robust (1/0),
   (d) all robust rows (up to 4 outgroup tables) agree. Same-track ⇒ S(YEAST,CANGA); different tracks ⇒ D[WGD].
   Single copies on different tracks = **pseudo-orthologs ⇒ D tagged WGD,LOSS** (hidden-paralogy test).
3. **WGD as allopolyploidisation (SYNTHESIS A2 / ME-2)**: every WGD D has lvl_lo=PostWGD, lvl_hi=ZT+WGD ancestor
   ⇒ merge ∈ [ZT+WGD, Saccharomycetaceae]. Tree per unit:
   `S@Saccharomycetaceae( S(KLE genes), U( S(ZT genes), D_WGD/S track structure for YEAST-CANGA, per-species D_WGD for
   other two-copy species, single copies ) )`. KLE vs post-WGD = orthologs at Saccharomycetaceae; ZT vs post-WGD
   unasserted. Caveat for the PI: MH2015 report that most ohnolog gene trees place the split on the *pre-KLE* branch
   ("more ancestral duplications are prominent with >50% ... mapped to the two nodes preceding the WGD (pre-ZT and
   pre-KLE)"); the A2 interval excludes a merge above Saccharomycetaceae. Y02 will quantify this.
4. Pre-WGD species in a pillar (≤1 per species, synteny-curated) ⇒ S polytomy per clade (T2 rule b).
5. Unit = one v7 pillar with ≥4 panel species mapped, ≥80 % of in-panel genes mapped, ≥1 asserted node.
   2,700 lines skipped (<4 panel species: mostly singleton/lineage-specific pillars), 41 for mapping <80 %, 12 with no
   asserted pair. complete_species = [] (truncated file ⇒ no genome-wide claim).
6. Tags: WGD on all ohnolog D nodes; LOSS on pseudo-ortholog D nodes and as unit phenomenon when a pillar mixes 1- and
   2-copy post-WGD species; ALLO as unit phenomenon (MH2015 quote). **CONV not tagged**: no stored source states
   conversion for specific RP pairs (Evangelisti & Conant 2010 not fetched).

## Mapping (UniProt 2026_03 ordered-locus/ORF names; source_cache/README)
Exact name match within species only (normalisation: case, `_`; Kpol `.`→`p`; AGOS_ prefix optional).
Accessions carrying several YGOB loci (identical copies merged in one entry) are excluded. Within units (excl. SACUV)
>98 % mapped; unmapped mostly KLULA (88) and VANPO (76) loci absent from current UniProt.
**SACUV: 0 mapped** — YGOB `Suva_` names refer to the Scannell 2011 assembly; UniProt UP001162085 (GCA_947243795) has no
cross-reference and the archived YGOB AA.fsa is truncated, so same-species sequence mapping is impossible. Listed as
unmapped in mapping_report.tsv. Not in panel: S. kudriavzevii, S. mikatae, L. kluyveri, L. waltii.

## Layout
units/<YGOB_Anc_x_y>/ {family.yaml, members.tsv, reference_tree.nhx} (each loadable by hobench load_family);
collection-level collection.yaml, members.tsv, mapping_report.tsv, assertions.tsv (rule rows with quotes + one row per
asserted node with file:line locators), build.py.

## Time
~2 h (incl. discovering the complete YGOB v1 track tables on Wayback).
