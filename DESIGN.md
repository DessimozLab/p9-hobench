# DESIGN — HOBench: a literature-based benchmark for hierarchical orthology inference

Version: v0 (PI draft, T+0h50). Review rounds: see `reviews/`. Changes after review are logged in §12.

## 1. Goal and scope (fixed by the brief)
Build (i) a curated compendium of reference gene families whose evolutionary histories are resolved in publications,
mapped onto (ii) panels of annotated species proteomes, plus (iii) metrics specific to *hierarchical* orthology (HOGs ≡
reconciled gene trees), (iv) a Python wrapper that runs HOG-producing tools (OrthoFinder 3.1.5, FastOMA 0.5.0) and
computes the metrics, and (v) an analysis of data difficulty (method performance + agreement + phenomenon tags).

## 2. Core concepts and formal definitions
- **Panel**: species set + rooted reference species tree T with named internal nodes = *levels*.
- **HOG at level L**: set of genes descending from one ancestral gene in the last common ancestor of the species under L.
- **Reference family F**: members M_F (UniProt accessions in panel proteomes) + a rooted, possibly multifurcating
  reference tree R_F with node events S (speciation), D (duplication, optional explicit level, tags WGD/TAND/…),
  U (unresolved; no assertion), H (HGT acquisition edge; see §2.2). Each node carries provenance (source id).
- **Merge level** m(a,b): the youngest level at which a and b are in the same HOG (reference: S → species LCA S0;
  D at level D → parent_T(D); prediction: lowest level where they share a HOG). Levels older than the family's
  taxonomic range (LCA of member species) collapse to ∞ ("never together within range"). This single quantity encodes
  orthology (m = S0), paralogy, and duplication timing.

### 2.1 Assertion tiers (every scored relation traces to one)
- **T1 — text-asserted**: an explicit statement/figure in the paper ("X and Y are ohnologs from the TGD", "Kinesin-2 present in LECA").
- **T2 — classification-derived**: the paper's curated labels/nomenclature (subfamily tables, YGOB pillars, gene names in
  curated nomenclatures) compiled into S/D/U structure by a documented rule (e.g. members of one labelled orthology group
  with ≤1 gene per species → S polytomy; labelled subfamilies known to pre-date level L → D at L).
- **T3 — tree-derived**: published tree topology, only clades with support ≥70% BS / ≥0.95 PP, events labelled by
  species-overlap rule; flagged as weaker (our reconciliation of their tree).
Headline numbers use T1+T2; T3 reported separately.

### 2.2 HGT
HGT families are encoded as two reference sub-trees (donor-lineage, recipient-lineage) joined by an `H` node; pairs across
an H node are asserted as "not orthologous at the recipient's level" (xenologs) and scored by a dedicated metric
(xenolog separation) rather than by merge level, because HOG semantics do not define a correct level for xenologs.

## 3. Species panels (survey C; UniProt reference proteomes 2026_03 / QfO 2026_02; CC BY 4.0)
| Panel | Scope | Spp | Proteins | Key events covered |
|---|---|---|---|---|
| Y | Saccharomycetaceae + outgroups | 16 | 85k | yeast WGD (hybrid origin), reciprocal loss, HGT (URA1) |
| V | Vertebrata (+amphioxus, Ciona) | 14 (+2: trout, pike for salmonid 4R) | ~370k | 1R/2R, 3R, 4R, *X. laevis* allotetraploidy |
| P | Viridiplantae | 13 | 389k | ζ/ε, γ, ρ, α/β, Brassica triplication, soybean, maize |
| E | Eukaryota-wide + 5 archaea + 4 bacteria | 34 | 340k | LECA paralogs, long branches (Microsporidia, Giardia), archaeal origins |
| K | Prokaryotes | 21 | 68k | LUCA paralogs, HGT, single-copy controls |
| D (stretch) | 12 *Drosophila* + outgroup | 13 | ~190k | recent LSE/tandem P450s, ILS-prone short internodes |
Species trees: NCBI topology, polytomies resolved from cited phylogenomic studies; contested nodes recorded and used in the
species-tree perturbation experiment (§7). **Rule:** panel species are chosen to maximise overlap with the species sampled
by the source papers (genome-wide surveys), because reference completeness is only known for sampled species.

## 4. Reference families (initial target ≈ 30, ≥4 per panel, chosen for diversity)
Selection criteria: machine-readable literature data; curated (not tool-derived; OMA/OrthoFinder/Compara/PANTHER-derived
sources excluded for orthology truth); diversity across clade × depth × rate × complexity × phenomenon.
- **Y**: YGOB pillars (WGD ohnologs + single-copy pillars; genome-scale collection) [T2]; Marcet-Houben 2015 ohnolog pairs
  [T1/T3]; URA1 HGT (Hall 2005) [T1]; GAL1/GAL3 [T1].
- **V**: globins [T1/T2]; opsins [T2]; Hox clusters [T2 nomenclature]; Sox (Voldoire 2017) [T2]; teleost 3R ohnologs
  (Braasch 2016) [T2, genome-scale]; salmonid 4R AORe/LORe (Robertson 2017) [T1/T3]; KRAB-ZNF [T2, merging stress];
  chemokines (Nomiyama 2013) [T2]; *X. laevis* L/S homeologs (nomenclature .L/.S) [T2].
- **P**: plant homeobox classes (Mukherjee 2009) [T2]; NLR (RefPlantNLR) [T2/T3]; bHLH (Pires & Dolan 2010) [T2];
  GID1 within CXE (Yoshida 2018) [T2]; MIPs (Danielson & Johanson 2008) [T3]; MADS-box [T2, if curatable].
- **E**: kinesins (Wickstead 2010) [T2]; myosins [T2]; ARF/Arl/Sar (Vargová 2021) [T2]; Rab (Klöpper 2012) [T2]; tubulins
  (Findeisen 2014) [T2]; MCM2–9 (Liu 2009) [T2/T3, control]; dynein HCs [T3]; adaptins/COPI (Hirst 2014) [T2];
  B-family DNA pol (Tahirov 2009) [T3]; kinome (KinBase) [T2].
- **K**: ribosomal proteins (Yutin 2012; single-copy negative control) [T2]; EF-Tu/EF-G (LUCA paralogs) [T1];
  nitrogenase Nif/Anf/Vnf/Bch (Garcia 2020) [T2]; β-lactamase classes A/B/C/D (non-homologous classes ⇒ merge decoys) [T2].
Each family gets phenomenon tags {WGD, TAND, LSE, LOSS/hidden-paralogy, HGT, FAST, PROM (promiscuous domains), FUS,
DEEP, CONTROL} and descriptors (depth of oldest/youngest asserted event, size, mean pairwise identity, #duplications).

## 5. Curation pipeline (agents; blind to tool outputs)
1. `curation_tools.fetch_paper` → OA full text + supplements (Europe PMC), checksums in `literature/papers/<id>/meta.json`.
2. Extract identifiers + labels (supplement tables/trees/text) → `families/<ID>/raw_members.tsv` with provenance (file, row, quote).
3. Map to panel proteomes: UniProt ID mapping; else same-species sequence match (exact, or ≥95% identity & ≥80% coverage,
   unambiguous). Unmapped members are recorded, never silently dropped.
4. Compile reference tree (S/D/U/H) with explicit level placement of literature events; `hobench` validator must pass.
5. **Independent verification**: a different agent re-derives a random 10% of assertions + all T1 events from the source and
   reports discrepancies; families with >5% discrepancy are fixed or demoted.
6. Curators never see OrthoFinder/FastOMA results (curation frozen — git tag `refs-frozen-v1` — before any scoring).

## 6. Metrics (implemented in `hobench/metrics.py`; see design/core_spec.md)
Primary (per family, per level):
1. **Per-level HOG pair P/R/F1** over asserted pairs (depth profile).
2. **Merge-level accuracy** + signed level error (over-split vs over-merge rates).
3. **Event recovery** for every literature duplication (recovered / misplaced by k levels / missed); WGD events aggregated.
4. **Orthology P/R** (QfO-comparable, derived from merge levels) and three-class relation accuracy (ortholog / in-paralog / out-paralog).
5. **Family integrity**: split excess (predicted top-level HOGs minus reference-implied), contamination by members of other
   reference families or declared non-homologous decoys (merge), coverage (members placed in any HOG).
6. **Per-gene / per-species correctness** (fraction of a gene's asserted relations correct).
7. **Xenolog separation** for HGT families.
Robustness metrics (§7): Δ of metrics 1–5 under controlled perturbation.
Aggregation: macro over families; per tag, per panel, per depth bin; 95% CI by bootstrap over families (fixed seeds in
`hobench/seeds.py`); paired method comparisons by paired bootstrap on family scores. No single composite score.
Difficulty: per-family mean error across methods; inter-method agreement (ARI/VI of HOG partitions per level on family genes).

## 7. Experiments
- **E1 Main**: OrthoFinder 3.1.5 and FastOMA 0.5.0 on each panel with the same reference species tree; baselines:
  OrthoFinder flat orthogroups (no hierarchy), all-singletons, one-HOG-per-family-at-root, random (positive/negative controls).
- **Data mode**: full proteomes for Y and K; *family-neighbourhood subsets* (all panel proteins with DIAMOND e ≤ 1e-5 to any
  reference member, plus a fixed random background of 1,000 proteins/species, seed 20261026) for V, P, E (compute: full V/P/E
  ≈ 4–8 h per tool-run). Validation control: run both modes on Y and K and compare metrics (must agree within CI; else
  report the bias).
- **E2 Robustness** (on subsets): (a) annotation fragmentation — split 10% of family members into two fragments;
  (b) gene fusion — fuse 5% of members with a random non-family protein; (c) contamination — insert 20 proteins from a
  distant species into one proteome; (d) species-tree error — swap a contested node / collapse to NCBI polytomies;
  (e) taxon sampling — remove pre-WGD outgroups (Y) / gar (V). Metric deltas vs clean run on unaffected + affected genes.
- **E3 OMAmer leakage check**: FastOMA with a reduced OMAmer database (clade-level or LUCA minus panel species, if feasible).

## 8. Software deliverables
`hobench/` (formats, parsers, metrics, controls, CLI), `curation_tools/`, `wrapper/` (`hobench run --tool orthofinder|fastoma
--panel V --mode subset` → normalized predictions + metrics), `families/` (reference data), `panels/` (species lists, trees,
manifests), `scripts/` (build steps), `tests/` (toy truths, controls). Work dirs on local NVMe /tmp; RAM cap 46 GB.

## 9. Reproducibility
Git-tracked code, references, results tables; DATA.md (source, version, checksum, licence per dataset); fixed seeds; one
`make`-style driver script per stage; no restricted data committed (PDFs, Phytozome). Every metric validated on toy data
with hand-derived answers before use.

## 10. Risks and mitigations
| Risk | Mitigation |
|---|---|
| Reference incompleteness (no negative set) | Score only asserted relations; precision penalties only within declared-complete species |
| Identifier drift / mis-mapping | sequence-based mapping with identity/coverage thresholds; mapping report per family |
| Agent hallucination of labels/IDs | provenance per assertion; independent 10% re-derivation; validator |
| Circularity (tool-derived "truth") | exclude OMA/OF/Compara-derived sources; OMAmer leakage experiment |
| Bias toward well-studied families | report per-tag/per-panel with n; do not over-generalize |
| Compute | subsets + validation; /tmp workdirs; omamer ≤2 parallel |

## 11. Timeline (T = hours since START)
T1–2.5 design reviews; T2–8 curation (parallel) + proteome download + subset building; T4–14 tool runs; T8–16 scoring,
verification, robustness; T16–21 analysis, follow-ups; T21–23 REPORT.md.

## 12. Change log
- v0 — initial draft.
