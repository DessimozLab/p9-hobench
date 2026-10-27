# HOBench — a literature-based benchmark for hierarchical orthology inference

**Status:** final report (2026-10-27, T+19h).
Every number below is produced by a committed script into `results/`; file pointers are given in brackets.
Process record: `RESEARCH_TREE.md` (decisions D0–D11), `LESSONS.md`, `reviews/` (2 review rounds), `results/ERRATA.md`.

## Summary

We built HOBench, a benchmark that scores *hierarchical* orthology (HOGs ≡ reconciled gene trees) against evolutionary
histories asserted in the literature, rather than against flat orthogroups or pairwise orthology.

- **Data.** 5 species panels (102 panel slots, 90 unique species; UniProt Reference Proteomes 2026_03, CC BY 4.0): yeasts Y (16), prokaryotes K
  (21), eukaryotes-wide E (34), vertebrates V (17), plants P (14), each with a binary, literature-resolved species tree
  whose internal nodes are named taxonomic levels and an alternative topology at a contested node.
- **References.** 32 curated catalogue rows (13 headline-scored items + 3 controls; the rest secondary or decoys): ~31.8k
  genes, 85 text-asserted (T1) and 2,138 classification-derived (T2) duplication events, 1 HGT, plus 5 pairs of
  non-homologous isofunctional enzymes as merge decoys. Two genome-scale collections: YGOB yeast WGD pillars
  (1,976 units; hidden-paralogy test) and *Xenopus laevis* L/S homeologs (498 units after erratum). Every scored relation has a stored
  source and a machine-verified quote; members were extracted by committed parsers (no typed IDs); 22 of 23 verification
  groups passed blind re-derivation (kill rule: ≥2/30 discrepancies; V01 failed with 2 → secondary); 0 denylist (tool-derived truth) violations; references frozen by hash.
- **Metrics.** A unified "merge level" formalism with interval-valued literature references, from which we compute four
  pre-registered primary endpoints: per-level HOG B-cubed F1 (depth profile), node-balanced merge-level accuracy,
  two-sided duplication-event recovery, and family-membership F1 — always shown next to trivial baselines and as skill.
- **Tools.** A Python wrapper (`wrapper/`) runs OrthoFinder 3.1.5 and FastOMA 0.5.0 with identical inputs and species
  trees, records manifests, and `scripts/evaluate_panel.py` scores any OrthoXML / OrthoFinder / TSV output.
- **Findings.** (i) In one pooled confirmatory test over 13 headline items (5 panels; the difference is driven mainly by
  panels E and V), OrthoFinder scores higher than FastOMA on per-level B-cubed F1 (+0.054, Holm p = 0.035); merge accuracy (+0.134, Holm p = 0.058), event recovery (+0.285, n = 11) and membership F1
  (+0.315, n = 5) point the same way but are not significant. (ii) Both tools fall into the yeast WGD hidden-paralogy trap:
  of 49 synteny-certified pseudo-ortholog pairs (reciprocal ohnolog loss), OrthoFinder classifies 4% and FastOMA 12%
  correctly, versus 87–91% on 1,189 true ortholog pairs. (iii) Whole-genome duplications are hard for both tools in
  opposite ways: OrthoFinder fails to keep ohnolog copies apart, FastOMA fails to keep each copy's orthologs together.
  (iv) FastOMA over-splits vertebrate homeologs/3R paralogs (Hox 2R/3R events: 8–9% recovered vs 77–82% for OrthoFinder)
  and over-merges duplication pairs on the eukaryote panel (24.4% of D-pairs, mostly at Eukaryota). (v) Annotation fragmentation is the dominant robustness
  threat (FastOMA merge accuracy −0.114 vs OrthoFinder −0.050 when both fragments must be placed correctly); contamination,
  a species-tree change at a contested node, and outgroup removal change primary endpoints by ≤0.005 (even with 200
  contaminants, ~half of which both tools absorb into reference HOGs as spurious paralogs). (vi) OrthoFinder's
  DendroBLAST mode has the best pooled merge accuracy (0.886) and event recovery (0.617) of three methods (not on every
  panel: on V it is below OrthoFinder default). (vii) We found no significant association between family difficulty
  and tool agreement (ρ = −0.41, p = 0.19, n = 12).
- **Main caveats.** Small number of independent families per panel (descriptive per panel; one pooled confirmatory test);
  three large panels scored in a validated subset mode; FastOMA's OMAmer database derives from OMA HOGs (possible
  circularity in its favour on OMA species); about 11 h were lost to three account usage-limit outages, forcing scope cuts.

## 1. Approach

1. **Literature & data survey** (4 parallel agents; all citations verified against Europe PMC) → `LITERATURE.md`,
   `literature/{A,B1,B2,C}_*.md`. Key gap: no real-data benchmark scores HOGs per taxonomic level or literature events
   (QfO flattens HOGs to pairs; OrthoBench is a single level; SwissTree is pairwise-scored).
2. **Design through critique.** `DESIGN.md` v0 → round 1 (7 independent reviewers: external PI, prior work, statistics,
   data, algorithms/compute, AI-agent failure modes, molecular evolution; synthesis in `reviews/round1/SYNTHESIS.md`) →
   v1 → round 2 (PI+statistics, data/agent audit; `reviews/round2/`). Round 1 found real defects (point-valued references
   mis-scored correct WGD predictions; trivial predictors scored well; β-lactamase classes wrongly called unrelated;
   FastOMA level names silently renamed) — all fixed before scoring.
3. **Curation** under a binding protocol (`design/CURATION_PROTOCOL.md`): blind to tool outputs, parsers not typing,
   quote-verified provenance, strict same-species ID/sequence mapping (≥98% identity, ≥80% coverage; ≥99% for
   homeologs/salmonids; v1.1 domain rule for motor-domain-only sources), ≥80% unique mapping for headline, blind 30-item
   re-derivation with a kill rule, denylist of tool-derived sources (OMA, OrthoFinder, Ensembl Compara/HCOP, VGNC,
   OHNOLOGS, PANTHER, TreeFam, eggNOG, OrthoDB, PLAZA), freeze by `families/MANIFEST.sha256` + tag `refs-frozen-v1`.
4. **Human-in-the-loop data access.** `requested/REQUESTED.md` ranked 50 paywalled papers/data files; 49 were delivered
   (`requested/RECEIVED.tsv`); used where they contained machine-readable IDs (Session 2016 Xenopus homeologs; others
   lacked per-gene IDs: Scannell 2006, Manning 2002 for non-human kinomes, Parenicová 2003 beyond Arabidopsis).
5. **Tool runs → scoring → controls → analysis**, with every metric defect discovered by a control logged with
   before/after values (`results/ERRATA.md` E1–E5).

## 2. Data

### 2.1 Panels (`panels/<P>/`, `DATA.md`)
| Panel | Species | Proteins | Mode | Events covered |
|---|---|---|---|---|
| Y yeasts | 16 | 85k | full (+subset for validation) | allopolyploid WGD, reciprocal loss, HGT (URA1) |
| K prokaryotes | 21 | 68k | full (+subset) | LUCA paralogs, single-copy controls, analogous-enzyme decoys |
| E eukaryotes + 9 prokaryote outgroups | 34 | 340k → 25.9k subset | subset | LECA paralogs, long branches |
| V vertebrates (+amphioxus, Ciona, trout, pike) | 17 | 448k → 70k subset | subset | 2R, 3R, *X. laevis* allotetraploidy, tandem Hox |
| P plants (+Phaseolus) | 14 | 418k → 25.4k subset | subset | plant WGDs, HGT (MIP NIPs), nested families |
Species trees: pruned/resolved from Shen 2018 (yeasts), Coleman 2021 (prokaryotes), Burki 2020 (eukaryotes), fish/plant
phylogenomics; contested nodes and citations in `panels/<P>/levels.tsv`. Subsets: DIAMOND neighbourhood of frozen
literature members + capped domain-only decoys + seeded background (`panels/<P>/subset/manifest.json`).

### 2.2 Reference families (`results/catalogue/families.md`, `coverage.md`, `abandoned.md`)
Headline: Y01 YGOB pillars (collection, 1,976 units; WGD/LOSS/ALLO), Y03 URA1 HGT, E01 kinesins (Wickstead 2010),
E02 ARF/Arl/Sar (Vargová 2021), E03 Rab (Klöpper 2012), E04 tubulins (Findeisen 2014), E07c adaptin scaffolds
(Hirst 2014), V02 visual opsins (Lagman 2013), V03 Hox clusters (curated nomenclature), V07 *X. laevis* L/S homeologs
(Session 2016; collection, 498 units), P02 bHLH (Pires & Dolan 2010), P03 GID1 within CXE (Yoshida 2018), P04 MIPs
(Danielson & Johanson 2008). Decoys (contamination/merge test only, no events): K04a–j analogous enzymes (Omelchenko
2010). Controls: K01 ribosomal proteins
(40 single-copy units, Yutin 2012), K02 EF-Tu/EF-G (LUCA paralogs), E05 MCM2–9 (Liu 2009). Secondary (mapping < 80%
or failed verification): E06 myosins, E07a/b, P01 plant homeobox, V01 globins, V04 Sox. Abandoned with reasons
(e.g. Braasch 2016 teleost 3R tables derive from Ensembl Compara → denylisted; nitrogenase/GAPDH-EGT/GAL-HGT had
< 4 panel species; salmonid 4R dropped for budget; *Coregonus* proteome downloaded empty).

Phenomenon coverage (catalogue rows): WGD 5, ALLO 3, LOSS/hidden paralogy 17, tandem/LSE 9, HGT 3, deep (LECA/LUCA) 17,
promiscuous domains 4, analogous decoys 10, controls 3. Event depth (headline T1+T2): 2,136 recent, 38 intermediate,
5 deep, 13 LECA/LUCA. **Gaps:** no family with literature-resolved ILS, fusion/fission, or a dedicated fast-evolving
lineage reached headline; HGT is represented by one headline family.

## 3. Metrics (`hobench/metrics.py`; `design/core_spec.md`; DESIGN §2, §6)
- **Merge level** m(a,b): youngest species-tree level at which genes a,b share a HOG. Literature references give an
  **interval** [lo, hi] (S node → species LCA; duplication on branch X → parent(X); contested timing, allopolyploidy, HGT
  and family-root duplications → intervals, possibly unbounded). A prediction is correct iff its merge level is inside.
  This one quantity encodes orthology, paralogy, duplication timing and per-level HOG membership (theory: merge
  matrices of reconciled trees are symbolic ultrametrics).
- **Primary endpoints** (pre-registered): (1) per-level HOG B-cubed F1 (gene-balanced), averaged over levels with a depth
  profile; (2) node-balanced merge-level accuracy; (3) two-sided event recovery (copies apart at the event, each copy's
  orthologs together; testable internal-branch events only); (4) family-membership F1 over declared-complete species,
  counting only reference-asserted separations and declared decoys as errors.
- **Secondary:** OrthoBench/OF3 root-level measures (P/R/F, fission/fusion/exact %), within-family orthology P/R,
  signed level error (over-split vs over-merge), ±1-level tolerance, contamination by decoys, in-paralog recovery
  (leaf events), S/D-balanced merge accuracy (post-hoc), baseline-normalised skill.
- **Statistics:** unit = family (a collection = one item after cluster aggregation); per-panel results descriptive;
  one confirmatory paired permutation test per primary endpoint (OF vs FastOMA, pooled headline items, Holm over 4).
- **Validation:** 95 unit tests with hand-derived answers; mutation check 10/10 mutants killed
  (`design/MUTATION_REPORT.md`); positive control perfect on every family; negative controls (all-singletons, one big
  HOG, random relabel) at their analytic values; control-detected defects fixed and logged (E1 untestable events, E2 flat
  projection, E4 leaf events, E5 unresolved-root membership).

## 4. Software
`hobench/` (formats, parsers for OrthoXML/OrthoFinder/TSV, metrics, controls, stats, freeze), `curation_tools/` (paper
fetching, ID extraction/mapping, sequence cross-referencing, proteome fetching), `wrapper/` (`python -m wrapper.run
--tool orthofinder|fastoma --panel P --mode full|subset`; subset builder; perturbation generator), `scripts/`
(evaluation, confirmatory, equivalence, catalogue, analysis, audit, verification), `README.md` quickstart.

## 5. Results

### 5.1 Main table (headline items; mean over items; skill vs best trivial baseline) [`results/analysis/main_table.md`]
| Panel (mode) | Method | Level B-cubed F1 | Merge acc. | Event rec. | Membership F1 |
|---|---|---|---|---|---|
| Y (full), n=2 | OrthoFinder | 0.928 | 0.879 | 0.060 | NA |
|  | FastOMA | 0.920 | 0.942 | 0.134 | NA |
|  | best trivial baseline | 0.790 | 0.872 | 0.000 | NA |
| E (subset), n=5 (membership 3) | OrthoFinder | 0.981 | 0.729 | 0.400 | 0.870 |
|  | FastOMA | 0.906 | 0.531 | 0.025 | 0.563 |
|  | best trivial baseline | 0.849 | 0.838 (one big HOG) | 0.000 | 0.180 |
| P (subset), n=3 (events 2, membership 1) | OrthoFinder | 0.994 | 0.907 | 0.667 | 0.493 |
|  | FastOMA | 0.997 | 0.867 | 0.667 | 0.049 |
| V (subset), n=3 (events 2, membership 1) | OrthoFinder | 0.987 | 0.960 | 0.907 | 0.721 |
|  | FastOMA | 0.881 | 0.707 | 0.204 | 0.511 |
K has no headline item with asserted pairs (K02 is a control; K04 decoys test contamination only).
Note: on E, *one big HOG* beats both tools on node-balanced merge accuracy, because most asserted nodes are speciation
polytomies and both tools over-split speciation pairs (FastOMA 66.7%, OrthoFinder 20.4% of S-pairs). Skill scores and the
post-hoc S/D-balanced variant are reported for this reason.

### 5.2 Confirmatory test (OrthoFinder − FastOMA, pooled headline items, paired permutation, Holm) [`results/confirmatory.tsv`]
| Endpoint | n | Mean difference | Holm p |
|---|---|---|---|
| Per-level B-cubed F1 | 13 | +0.054 | 0.035 |
| Node merge accuracy | 13 | +0.134 | 0.058 |
| Event recovery | 11 | +0.285 | 0.12 |
| Membership F1 | 5 | +0.315 | 0.12 |

### 5.3 Phenomena [`results/analysis/phenomena.md`]
- **Hidden paralogy (yeast WGD, reciprocal loss):** pseudo-ortholog pairs correctly called paralogs — OrthoFinder 0.041
  [0.000, 0.105], FastOMA 0.122 [0.038, 0.224] (n = 49); same-track orthologs 0.910 / 0.865 (n = 1,189).
- **WGD events** (243 testable, V+Y): recovery OF 0.22, FastOMA 0.24; OF fails on separation (copies merged too young;
  separation 0.36), FastOMA on cohesion (0.61). Within units, WGD D-nodes score far below S-nodes (YGOB: OF −0.68,
  FastOMA −0.21; n = 532 units).
- **Hox 2R / 3R** (V03): events recovered OF 0.77 / 0.82, FastOMA 0.08 / 0.09, flat orthogroups 0 / 0 (n = 13 / 11).
- **Allopolyploid homeologs** (*X. laevis*, 333 units in subset): L/S pairs correctly placed OF 0.988, FastOMA 0.429
  (FastOMA over-splits 55.6%); homeolog vs speciation nodes within units: FastOMA −0.26 [−0.29, −0.23], OF −0.005.
- **Ancient paralogs** (E): FastOMA merges 24.4% of pairs under duplication nodes too early (24.0% at Eukaryota; OF 0.1%).
- **LUCA paralogs (K02 EF-Tu/EF-G, control):** neither tool recovers the event; failure is cohesion (OF splits bacterial
  vs archaeal EF-Tu; FastOMA builds no LUCA-level EF-G group).

### 5.4 Error anatomy and difficulty [`results/analysis/error_anatomy.md`, `difficulty.md`]
- P02 bHLH: FastOMA membership precision 0.959 but recall 0.025 — the 435 members are scattered over 186 groups
  (median 1 member).
- V: FastOMA errors are over-splitting (over-merging ≤ 1.7% at every level); zebrafish 3R nodes 91% over-split.
- Difficulty (1 − mean skill over defined endpoints, 3 methods) vs 3-way agreement (mean pairwise AMI): Spearman
  ρ = −0.41 (p = 0.19, n = 12 headline items; an earlier ρ = −0.58 used n = 20 including non-headline families) — a
  non-significant trend. Hardest headline families: E02 ARF, E01 kinesins, P02 bHLH (short/promiscuous domains, deep
  paralogy).
- Third method (OrthoFinder `-M dendroblast`, pooled over E, P, V, Y): merge accuracy 0.886 vs 0.846 (OF default) vs
  0.712 (FastOMA), n = 13; event recovery 0.617 vs 0.479 vs 0.194, n = 11 [`results/analysis/main_table.md`].

### 5.5 Subset-mode validation (pre-registered equivalence, TOST ±0.05) [`results/equivalence.md`]
On Y (subset = 54% of proteome) and K (6%), subset − full deltas lie within ±0.05 (90% paired CI) for all primary
endpoints and tools (e.g. Y OrthoFinder level F1 −0.0003 [−0.0008, +0.0002], n = 1,976 units), except K event recovery
(n = 1 event; flat baseline flips) which is inconclusive. Pair merge-level identity subset vs full: 0.961–0.996.
Equivalence was shown only on Y and K; E/P/V results are reported **conditional on subset mode** (not tested there).

### 5.6 Robustness (Y, full mode; 1,977 units; seed 1) [`results/analysis/robustness.md`, `results/robustness/`]
Noise floor: an OrthoFinder replicate is bit-identical to the original run (0 differences over 4,156 unit-metrics), so
nonzero deltas are attributable to the perturbation (FastOMA was not replicated).
| Perturbation | OrthoFinder Δ merge acc. | FastOMA Δ merge acc. | Comment |
|---|---|---|---|
| Fragmentation of 10% of members (2,615 genes → 2 pieces), credit "both fragments" | −0.050 [−0.057, −0.044] | −0.114 [−0.122, −0.106] | affected genes lose ~25–28% correctness |
| Fragmentation, credit "either fragment" (pre-registered primary rule) | −0.003 [−0.008, +0.002] | −0.024 | FastOMA's *unaffected* genes also degrade (−0.015) |
| Contamination (20 human proteins into *K. lactis*) | +0.000 | −0.000 | only 8–13 units touched |
| Contamination, follow-up (200 human homologs of Y members into *K. lactis*) | +0.0004 [−0.0003, +0.0012] | −0.0010 [−0.0036, +0.0015] | 99–105 units touched; see below |
| Alternative species tree (1 contested node) | −0.001 | −0.004 | level F1 −0.003 both |
| Taxon removal (ZT clade: 2 pre-WGD species) | +0.001 | −0.003 | — |
Contaminant fate (follow-up, `results/Y_contam200_full/contaminant_placement_*.tsv`): 107/200 (FastOMA) and 103/200
(OrthoFinder) contaminants were placed in a HOG together with yeast reference genes — FastOMA mostly at
Saccharomycetaceae (119 placements at that level; i.e. as a deep paralog), OrthoFinder mostly at the shallow
Kluyveromyces–Eremothecium level (79; i.e. as a recent in-paralog). Neither tool isolates foreign sequences, but
the relations among the genuine reference genes are unaffected; contamination therefore shows up as spurious
lineage-specific duplications rather than as broken orthology — detectable only with dedicated QC (e.g. OMArk).

Interpretation: within this design, annotation fragmentation is the main robustness threat and FastOMA is more sensitive
to it; the alternative-tree and taxon-removal perturbations were mild (one node, one clade) and a 10× stronger contamination
follow-up still left reference relations intact; harsher species-tree errors are a follow-up.

## 6. Limitations
1. **Few independent families per panel** (2–5 headline items): per-panel numbers are descriptive; only the pooled
   confirmatory test is inferential, and it has n = 5–13.
2. **Reference incompleteness:** only asserted relations are scored; many families are partly unresolved (U), and
   declared-complete species are few, so membership F1 is testable for only 5 items.
3. **Selection bias** toward well-studied, machine-readable families; phenomena without headline families: ILS,
   gene fusion/fission, fast-evolving lineages; HGT has one headline family.
4. **Circularity risks:** FastOMA places proteins with OMAmer against OMA HOGs (not testable within budget); K01/K04 use
   NCBI COG membership to locate the papers' proteins (shares COG's clustering assumptions); YGOB is synteny-based.
5. **Interval hedges:** the yeast WGD upper bound (allopolyploid model) makes most YGOB events untestable (202/1,592
   testable); V01 globins failed blind verification because of interval widening.
6. **Subset mode** for E/V/P is validated only on Y/K; the V subset hit the 70k cap (seeded with 300 sampled V07 units;
   333 of 498 units were fully present and scored).
7. **Two tools** (+ OrthoFinder DendroBLAST and flat variants); single runs (OrthoFinder shown deterministic; FastOMA
   not replicated); robustness only on Y with mild perturbations; no OMA standalone, no simulation-based comparison.
8. **Process:** three account-limit outages (~11 h in total) forced cuts (salmonid 4R, chemokines, KRAB-ZNF, SwissTree
   extension — partially started, unused — MADS, PolB, kinome).

## 7. Follow-ups and remaining ideas
- Add SwissTree reconciled trees and salmonid AORe/LORe (Robertson 2017; requires a non-empty *Coregonus* proteome).
- Request the full YGOB Pillars.tab (only #01 not delivered) to extend hidden-paralogy tests to all post-WGD species.
- More methods (OMA standalone, SonicParanoid2/Broccoli as flat baselines, GeneRax-reconciled trees) to make agreement
  analyses informative; precomputed OMA/Compara HOGs as methods (never truth).
- Full-proteome runs for V/E/P to remove the subset caveat; OMAmer database without panel species for FastOMA.
- Dated species trees (TimeTree) to express level errors in Myr; ILS-aware tolerance.
- Literature-mining at scale: the curation toolkit + protocol generalise to more families (kinomes with KinBase files,
  MADS with per-species tables, NLR ancestral lineages from Shao 2016).

## 8. Reproducibility
Code, references, run manifests, compressed HOG outputs and all result tables are in git; `DATA.md` lists source,
release, checksum and licence for every dataset; seeds are fixed (`hobench/seeds.py`, dataset seed 20261026); references
and metrics are frozen by hash/tag with errata logs; `README.md` gives the end-to-end commands.
