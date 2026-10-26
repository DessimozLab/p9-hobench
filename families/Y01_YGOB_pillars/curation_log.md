# Y01_YGOB_pillars — curation log (collection YGOB, panel Y)

Status: **v1 built and validated against the final panel Y (panels/Y/species_tree.nwk): 1,977 units, 0 errors / 0 warnings.**

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

## Verification (final panel Y)
- hobench validate_family: 1,977/1,977 units OK (0 errors, 0 warnings) on panels/Y/species_tree.nwk; also valid on a
  tree with non-monophyletic KLE (build derives level roles by LCA; WGD lvl_hi = child of the youngest
  post-WGD/KLE ancestor present in the unit).
- scripts/verify_provenance.py: **PASS (0 errors, 0 warnings; checks a–f)** on an assembled single-family view
  (family.yaml from collection.yaml + all unit trees under one root U node). The script requires family.yaml at
  families/<ID>/, which this collection layout deliberately omits so that hobench load_families does not load 26k genes
  as one family; PI to decide how collections are routed. README quotes are matched in `README.txt` / `v1_README.txt`
  (verbatim copies of the Wayback README files, created because the verifier only searches known text suffixes).
- Panel Y has no HENBL, KAZNA, SACUV (YGOB columns → not_in_panel). Final: 26,173 members, 13 species, 1,591 WGD D
  nodes (1,341 single-species ohnolog pairs, 201 YEAST/CANGA track-resolved, 49 pseudo-ortholog D tagged LOSS),
  1,189 same-track YEAST–CANGA S nodes, unique mapping 98.8 %.

## PI decision (after Y02): widened WGD interval
All 1,591 WGD D nodes now carry lvl_lo=PostWGD, lvl_hi=Saccharomycetaceae ⇒ ohnolog merge ∈ [ZT_PostWGD,
parent(Saccharomycetaceae)], because MH2015 find that a substantial fraction of ohnolog pairs coalesce before the KLE
divergence, consistent with an allopolyploid origin: "We found that the pre-KLE duplication peak was more apparent in
the subset of duplications leading to conserved pairs of ohnologs, which indicates that this ancestral duplication peak
is indeed related to the observed WGD paralogous blocks" and "the main result that the majority of ohnologs have
apparently diverged before the expected WGD" (PMID26252497; Y02 cross-check: 85/245 pairs labelled 'A - AB').
Consequence: in the 525+6 units that contain a WGD D node, KLE vs post-WGD pairs are no longer asserted (root U;
KLE vs ZT kept as S at Saccharomycetaceae); units without a WGD D keep KLE vs post-WGD orthologs at Saccharomycetaceae.
Re-validated: 1,977 units, 0 errors; verify_provenance (assembled view) PASS, 0 errors/0 warnings.

## PI decision: URA1 pillar excluded
Unit YGOB_L04856 (YKL216W/URA1, DHOD acquired by HGT) is dropped via EXCLUDED_UNITS in build.py; it is encoded with its HGT history in families/Y03_URA1_HGT. Now 1,976 units; validate 0 errors; verify_provenance PASS.

## Scannell et al. 2006 (Nature, PMID16421571) delivered — no extension possible
Files (requested/02_*, sha256 in requested/RECEIVED.tsv) copied to literature/papers/PMID16421571/raw/ (git-ignored; text extracted with pypdf). Supp1.xls = GO-category table only; supp2–5 = notes (phylogeny, timing, model, GO); main text gives class counts, no per-locus lists (0 Scas_/CAGL0 identifiers in any file). Its Methods describe exactly the YGOB track scoring already used here (v1 Synteny_*.tab: "scored up to 18 times: on tracks A and B in each of the three post-WGD species"). The only S. castellii track data (YGOB v1, Scas_ contig names of Cliften 2003) cannot be mapped to NAUCC (no sequences: v1 Proteins.fsa archive truncated, S. cerevisiae only). Hidden-paralogy assertions therefore remain YEAST–CANGA only; no change to units.
