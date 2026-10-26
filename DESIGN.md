# DESIGN — HOBench: a literature-based benchmark for hierarchical orthology inference

Version: v1 (after review round 1, T+1h05). Review rounds: see `reviews/`. Changes are logged in §12 with review ids.

## 1. Goal and scope (fixed by the brief)
Build (i) a curated compendium of reference gene families whose evolutionary histories are resolved in publications,
mapped onto (ii) panels of annotated species proteomes, plus (iii) metrics specific to *hierarchical* orthology (HOGs ≡
reconciled gene trees), (iv) a Python wrapper that runs HOG-producing tools (OrthoFinder 3.1.5, FastOMA 0.5.0) and
computes the metrics, and (v) an analysis of data difficulty (method performance + agreement + phenomenon tags).

## 2. Core concepts and formal definitions
- **Panel**: species set + rooted, binary reference species tree T with named internal nodes = *levels*
  (`panels/<P>/species_tree.nwk`, `levels.tsv` with NCBI/literature names and citations; ids `^[A-Za-z0-9_-]{3,}$`).
- **HOG at level L**: set of genes descending from one ancestral gene in the last common ancestor of the species under L.
- **Reference family F**: members M_F (UniProt accessions in pinned panel proteomes) + rooted, possibly multifurcating
  reference tree R_F with node events S, D, U (unresolved; no assertion), H (HGT acquisition) and provenance per node.
- **Merge level** m(a,b): youngest level at which a and b are in the same HOG. **Reference merge levels are intervals**
  [lo, hi] (hi may be ∞) [A1]: S node → [S0, S0] (S0 = species LCA); D node with `lvl` → [parent(lvl)]; with
  `lvl_lo`/`lvl_hi` → [parent(lvl_lo), parent(lvl_hi)]; default for D without level = [parent(LCA mapping), tightest
  ancestor bound]; cross-copy pairs of a duplication at the family root → [parent(X), ∞] unless `root_single_copy: true`.
  A prediction (point value, lowest level where both share a HOG, ∞ if never) is correct iff it lies in the interval;
  signed error = distance to the nearest bound (+ = too old = over-split; − = too young = over-merged). Exact-level and
  interval pairs are reported separately.
- **Theory**: the matrix of merge levels of a reconciled tree is a symbolic ultrametric; orthology relations of a
  reconciled tree are cographs (Hellmuth et al. 2013; Lafond & El-Mabrouk 2014). This is why consistency checks on the
  reference (validator) are necessary and why per-level HOG partitions and merge levels are equivalent views [A8].

### 2.1 Assertion tiers (every scored relation traces to one)
- **T1 — text-asserted**: explicit statement in the paper, quote stored and machine-verified.
- **T2 — classification-derived**: curated labels compiled by the rules of `design/CURATION_PROTOCOL.md` §3. The rule
  "≤1 gene per species ⇒ S" is **not** applied across WGD-descended lineages unless the set is synteny-curated
  (YGOB, Xenbase L/S, Hox clusters) — reciprocal ohnolog loss creates pseudo-orthologs [A3].
- **T3 — tree-derived** (published topology, support ≥70% BS / ≥0.95 PP); **T3fig** — read from figures, double extraction.
Headline = T1+T2; T3/T3fig reported separately.

### 2.2 Special histories
- **Allopolyploid WGD (yeast)** [A2]: hybrid origin (Marcet-Houben & Gabaldón 2015) ⇒ ohnolog pairs apart at PostWGD and
  below, merge ∈ [ZT+WGD ancestor, Saccharomycetaceae]; post-WGD vs KLE copies: orthologs at Saccharomycetaceae;
  post-WGD vs ZT: unasserted unless placed by a source. YGOB pseudo-orthologs are scored as paralogs.
- **Delayed rediploidisation (salmonid LORe)** [A5]: one D per salmonid lineage; requires *Coregonus* + salmon + trout + pike.
- **Homeologs (*X. laevis*)** [A6]: D on the XENLA leaf branch, tag ALLO.
- **HGT/EGT** [A4]: HOG semantics define xenolog levels: apart at every recipient-lineage level; merge at LCA(donor,
  recipient) if the donor lineage is in the panel, else ∞. Encoded as an `H` node (interval semantics); U as fallback.
- **Gene conversion** [A7]: tag CONV (yeast RP ohnologs, tufA/tufB, OPN1LW/MW); excluded from the headline; identical
  copies excluded as unmappable.
- **Contested event timing** (2R vs cyclostomes, plant ζ/ε/ρ/α/β) [A6]: events placed by panel branch with intervals.

## 3. Species panels (UniProt Reference Proteomes, release **2026_03 pinned**; canonical one-protein-per-gene FTP files; CC BY 4.0) [D1, D2]
| Panel | Scope | Spp | Key events covered | Data mode |
|---|---|---|---|---|
| Y | Saccharomycetaceae + outgroups | 16 | yeast WGD (hybrid origin), reciprocal loss, HGT (URA1, GAL) | full |
| K | Prokaryotes (14 bacteria, 7 archaea) | 21 | LUCA paralogs, HGT, single-copy controls, analogous-enzyme decoys | full |
| E | Eukaryota-wide + 5 archaea + 4 bacteria | 34 | LECA paralogs, long branches (Microsporidia, Giardia), EGT | subset |
| V | Vertebrata (+amphioxus, Ciona; + trout, pike, *Coregonus*) | 17 | 1R/2R, 3R, 4R (AORe/LORe), *X. laevis* allotetraploidy | subset |
| P | Viridiplantae (+ *Phaseolus* if compute allows) | 13–14 | ζ/ε, γ, ρ, α/β, Brassica triplication, soybean, maize | subset |
Panel D (Drosophila) dropped [D8]. Per panel: `panels/<P>/species.tsv` (code, UP id, taxid, release, sha256, n proteins),
`levels.tsv` (level id, NCBI taxid/name or literature name, contested flag, citation), binary named `species_tree.nwk`
(identical input to both tools and to scoring) and `species_tree_alt.nwk` (alternative at a contested node), generated by
`scripts/build_panels.py`. Validation: tree leaves = species codes = FASTA basenames = member species.
Rule: panel species are chosen to maximise overlap with species sampled by the source papers.

## 4. Reference families (≈25–30 target; ≥3 per panel; diversity over count)
Selection: machine-readable literature data; curated sources only (**denylist** [C3]: OMA, OrthoFinder, Ensembl
Compara/HCOP, Compara-seeded names (VGNC), OHNOLOGS, PANTHER, TreeFam, eggNOG, OrthoDB, PLAZA, LOC names, KRAB-ZNF ages);
diversity across clade × depth × rate × complexity × phenomenon. Each source records `source_method`.
- **Y**: YGOB pillars [collection, T2, synteny]; Marcet-Houben & Gabaldón 2015 ohnologs [collection, T1/T3]; URA1 HGT
  (Hall 2005) [T1]; *S. pombe* GAL-cluster HGT (Slot & Rokas 2010) [T1] [D9].
- **K**: ribosomal proteins restricted to proteins single-copy in every panel species [C2]; translation GTPases
  (EF-Tu/EF-G and other trGTPases; LUCA paralogs) [T1]; nitrogenase Nif/Anf/Vnf/Bch (if ≥4 panel species) [T2];
  **non-homologous isofunctional enzymes** as merge decoys (Omelchenko 2010: LysRS I vs II, FBP aldolase I vs II, SOD
  Cu/Zn vs Fe/Mn) [C1] (β-lactamase A/C/D are homologous — removed).
- **E**: kinesins, myosins, tubulins, ARF/Arl/Sar, Rab, adaptins/COPI/TSET, MCM2–9 (control), B-family DNA pol (FUS),
  GAPDH EGT (GapA/B vs GapC) [D9; dynein dropped].
- **V**: globins, opsins, Hox, Sox; teleost 3R ohnologs (Braasch 2016) [collection]; salmonid 4R AORe/LORe [collection];
  *X. laevis* L/S [collection, ≤500 units]; KRAB-ZNF (integrity metrics only); chemokines (if time).
- **P**: plant homeobox, bHLH, GID1/CXE, MIPs, NLR restricted to RNL/ZAR1/TNL loss [D9], MADS (if quotable).
- **Cross-panel**: SwissTree reconciled trees mapped onto panels with ≥4 panel species (T1-equivalent) [D10].
Homology declarations [B6]: Rab/Arf, kinesin/myosin, EF-Tu/EF-G, tubulin/FtsZ, MCM are homologous — merges at or above
their documented split are not penalised; contamination counts only declared decoys or `distinct_from` relations.
Fusion families declare an anchor domain; pairs with the second parent family are unasserted.
**Genome-scale collections** (YGOB, Braasch 3R, MH2015, *X. laevis* L/S, ribosomal proteins) are scored as a separate
track with two-stage aggregation (unit → collection) [B10].

## 5. Curation pipeline (binding protocol: `design/CURATION_PROTOCOL.md`) [C4–C9]
Agents curate blind to tool outputs; members come from committed per-family `build.py` parsers (no typed IDs); every
label/event/tag/decoy has a stored, machine-verified quote; mapping thresholds ≥98% identity (≥99% XENLA, salmonids, TAND),
≥80% coverage of both, reciprocal best hit, ≥2-point margin, same species only [D3]; ≥80% unique mapping for headline;
`scripts/verify_provenance.py` must pass; blind verifier re-derives min(all, 30) seeded assertions + all D-levels and decoys
(kill rule ≥2/30) [C6]; pre-freeze within-family DIAMOND QC (members with no in-family hit / best hit in another family
flagged) [C8]; freeze at T+12 h by `families/MANIFEST.sha256` + git tag `refs-frozen-v1`; metrics frozen by tag
`metrics-frozen-v1`; post-freeze fixes only via `families/ERRATA.md` with before/after numbers [C9].

## 6. Metrics (implemented in `hobench/`)
**Four primary endpoints** (pre-registered; Holm correction across OF-vs-FastOMA comparisons) [B8]:
1. **Per-level HOG F1, gene-balanced (B-cubed)** — depth profile + mean over levels [B2].
2. **Node-balanced merge-level accuracy** — mean over asserted reference nodes of the fraction of their pairs whose
   predicted merge level falls in the reference interval [A1, B2].
3. **Two-sided event recovery** per literature duplication: copies apart at the event level AND ≥θ of each copy's
   reference-together pairs together (θ ∈ {0.5, 0.8}); plus spurious-duplication count inside reference S clades [B3].
4. **Family-membership F1** over declared-complete (family, species) cells [B5].
Secondary/exploratory: pair-pooled versions; signed level error, also as fraction of family depth, and a ±1-level
tolerance variant [B12]; within-family orthology P/R (from merge levels and from OF `Orthologues/`; not QfO-comparable)
[B13]; root-level OrthoBench/OF3 measures (P/R/F, % fission, % fusion, % exact, dispersion entropy) [B7];
contamination = decoy genes /(own + decoy) with unlabelled homologs reported separately [B5, B6]; xenolog separation;
per-gene/per-species correctness; baseline-normalised skill (s−b)/(1−b), b = max(singletons, one-HOG, mean relabel) [B4].
Statistics [B9, B11]: unit = family (collections = one unit after cluster bootstrap); macro averages; paired family
bootstrap (same resampled IDs across methods, stratified by panel; seeds in `hobench/seeds.py`) + Wilcoxon signed-rank;
no CI when n < 5; n reported everywhere; phenomenon effects estimated within families (e.g. WGD-tagged D vs S nodes)
with cluster bootstrap. No composite score. Difficulty = 1 − mean skill across methods; agreement = AMI per level.

## 7. Experiments
- **E1 Main**: OrthoFinder 3.1.5 (`-s` binary panel tree) and FastOMA 0.5.0 (same tree; precomputed OMAmer hogmaps) per
  panel; variants: OF `-M dendroblast`, FastOMA without user tree (if supported) [E7]; baselines: OF flat orthogroups,
  all-singletons, one-HOG-per-family, random relabel; positive control (reference → perfect). FastOMA OrthoXML read with
  duplications at the enclosing group's declared level (primary) and pyHam-late placement (sensitivity) [E1]; levels
  mapped by clade species set, hard error on unknown levels [E2, E3].
- **Data modes** [D5, D6]: full proteomes for Y and K; family-neighbourhood subsets for E, V, P: DIAMOND `--more-sensitive`
  e ≤ 1e-5, fixed `--dbsize`, from raw paper member sequences (frozen before runs); + capped seeded sample of
  domain-only hits (promiscuous decoys) + seeded background fraction; ≤70k proteins per panel; `mode` column in every
  result row, never pooled. **Validation (pre-registered)**: per tool on Y and K, subset vs full: paired per-family Δ of
  primary endpoints within ±0.05 (TOST, 90% paired CI) AND per-pair merge-level identity; otherwise E/V/P results are
  reported as conditional on subset mode.
- **E2 Robustness** (Y, full mode, seed 1; contamination seeds 2–3 if on schedule) [E5]: (a) fragmentation of 10% of
  family members; (c) contamination (20 proteins of a distant species into one proteome); (d) alternative species tree at
  a contested node; (e) taxon removal (pre-WGD outgroups); (b) fusion = COULD. One clean replicate per tool as noise floor.
  Δ on matched pairs, split affected/unaffected genes, paired CI. Fragment scoring rule fixed before running: fragments
  inherit the parent's relations; the parent is credited if either fragment is placed correctly (also report both).
- **OMAmer leakage experiment dropped** (infeasible in budget) [E6]; circularity handled by the denylist; reported as a
  limitation (FastOMA's OMAmer database is built from OMA HOGs, which may include our panel species).

## 8. Software deliverables
`hobench/` (formats, parsers, metrics, controls, freeze, CLI), `curation_tools/`, `wrapper/` (tool runs with manifests,
subset builder, perturbations), `families/`, `panels/`, `scripts/` (build steps, verification, audit), `tests/`.
Release "HOBench-core" [F2]: exact input FASTAs (full + subset), trees and manifests with checksums (archive, not git);
references exported as per-level OrthoXML (U marked) and flattened pairs [F3]; family cards from family.yaml; coverage
table HOBench vs OrthoBench vs SwissTree (clade × depth × phenomenon) [F4].

## 9. Reproducibility
Git-tracked code, references, result tables; DATA.md (source, version, checksum, licence, redistribution, script); seeds
fixed (`hobench/seeds.py`, dataset seed 20261026); run manifests (commands, versions, input/DB/tree checksums) asserted at
scoring; every REPORT.md number generated by a script into `results/` with commit hash and checked by a fresh agent [C11];
metric tests reviewed/extended by a different agent incl. mutation checks [C10].

## 10. Risks and mitigations
| Risk | Mitigation |
|---|---|
| Reference incompleteness (no negative set) | score only asserted relations; membership FPs only in declared-complete cells |
| Identifier drift / mis-mapping | strict sequence mapping, per-family mapping report, ≥80% unique mapping for headline |
| Agent hallucination of labels/IDs/quotes | parsers not typing; mechanical quote/ID verification; blind 30-item re-derivation |
| Circularity (tool-derived truth) | denylist + `scripts/audit_sources.py`; OMAmer leakage stated as limitation |
| Bias toward well-studied families | per-tag/per-panel reporting with n; within-family phenomenon contrasts |
| Subset-mode bias | pre-registered equivalence test on Y and K |
| Compute / account limits | full mode only Y,K; ≤5 concurrent agents; /tmp workdirs; omamer ≤2 parallel |

## 11. Timeline (revised after losing T+1h05→T+4h55 to an account rate limit)
T+5–7: hobench v1 metrics, panels/trees/hogmaps, wrapper; Y/K full runs start. T+5–12: curation in waves of ≤4 agents
(Y, K first: vertical slice of ≥2 frozen families in Y and K scored end-to-end by T+9). T+12: freeze. T+9–15: subsets +
E/V/P runs. T+12–17: scoring, verification, robustness (Y). T+17–21: analysis, follow-ups, review of results.
T+21–23: REPORT.md.

## 12. Change log
- v0 — initial draft.
- v1 — round-1 review (reviews/round1/SYNTHESIS.md): interval reference merge levels (A1); hybrid yeast WGD (A2); T2
  WGD restriction (A3); HGT semantics (A4); salmonid LORe + *Coregonus* (A5); branch-placed contested events (A6); CONV
  stratum (A7); theory grounding (A8); metric fixes and four primary endpoints (B1–B13); β-lactamase decoys replaced by
  analogous enzymes (C1); single-copy r-proteins only (C2); denylist (C3); mechanical provenance + blind verification +
  freeze (C4–C9); pinned release, strict mapping, panel tables, DATA.md (D1–D4); subset design + pre-registered validation
  (D5, D6); OMArk covariate (D7); panel D dropped (D8); family list changes (D9); SwissTree (D10); FastOMA level reading
  and level mapping (E1–E3); critical-path runs (E4); narrowed robustness (E5); E3 leakage dropped (E6); method variants
  (E7); release/export/coverage table (F2–F4); literature corrections (F5).
