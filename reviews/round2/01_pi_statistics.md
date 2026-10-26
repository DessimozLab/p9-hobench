# Round 2 — external PI / statistics review (reviewer 01)

## (a) Round-1 must-items
Applied correctly: A1, A2, B1–B3, B5, B8, E1 (declared and late readings are identical), E2, E3 and C2.

Not done:
- **C9**: there is no MANIFEST, freeze tag or ERRATA. The slice was scored with `dirty=True`.
- **C6**: there are no `answers.tsv` for Y01, E06, E07a, E07b, P01, P04, V01 or V04.
- **C3**: `audit_sources.py` is missing.
- **D6**: the subset-vs-full equivalence test has not been run.
- **B5**: `complete_species` is empty for all of Y, K02, E03, E07*, P03, P04, V01, V02 and V07.

## (b) Endpoints on real data
1. **Event recovery is degenerate in Y.** all_singletons scores 0.861 on YGOB units, beating OF (0.24) and FastOMA (0.61); cohesion is vacuous without together-pairs. In K, the endpoint is a single event (K02) on which every tool scores 0.
2. **Level F1 and merge accuracy are near-saturated by lumping.** On YGOB, one-HOG scores 0.964/0.886, against OF 0.977/0.900 and FastOMA 0.967/0.900. On RPROT, one-HOG scores 1.0/1.0. S nodes reward lumping and D intervals [parent(X), ∞] reward splitting, so the raw means largely measure the reference's S:D mix.
3. **Membership F1 is precision-blind** wherever a family has no complete cells or decoys. One-HOG scores 0.999 on YGOB but 0.047 on RPROT.
4. **of_flat equals all_singletons on YGOB**, and scores below it on merge accuracy (0.057). This is a projection or ID bug.
5. **The K headline (n=2) is half a control (RPROT).** The FastOMA-vs-OF gap there (0.697 vs 0.976) is a control-track finding, not a "prokaryote" finding.

## (c) Coverage achieved vs claims
About 20 headline items, counting each collection as one, against a target of 25–30. Y and K have only 2–3 items each.

Tag coverage:
- **0 curated families**: ILS, fusion, fast evolution, contamination, species-tree error, annotation error, alignment error, EGT (E09 has no yaml).
- **1 family each**: CONV (excluded by design), NESTED, RETRO, SSD.
- **HGT**: Y03 (one event), P04, and E06 (non-headline).
- **Well covered**: WGD, LOSS, DEEP, PROM.

Cheapest closures:
- **FAST (≤1 h, no tool runs)**: pre-register a rate covariate built from reference sequences only (MAFFT+FastTree tree length per unit) over 1,976 YGOB units, the RPROT units and the E families. Add a per-gene long-branch stratum in E (Microsporidia, Giardia).
- **Fusion (~2 h wall)**: run E2(b) in Y with the existing perturb.py. Fuse members of two YGOB units in one species, and score spurious merging of the two units. Fix the scoring rule before the run.
- **Contamination, species-tree error, annotation error**: E2(c), (d) and (a) cover these. Label them "perturbations of real data", not curated phenomena. Also report decoy contamination for K04 and for the domain-only decoys.
- **ILS**: there is no cheap legitimate family. Report spurious duplications on single-copy units only as "gene/species-tree discordance".
- **EGT**: add E09 only if it passes C6 before the freeze.

## (d) Statistics at small n
- **Per panel (n=2–9)**: descriptive only. The summary currently prints min–max inside a column labelled "CI".
- **Confirmatory**: for each primary endpoint, run one exact paired sign-flip permutation test of OF − FastOMA over all headline items, stratified by panel, with Holm correction over 4 and the median paired difference as the effect. Wilcoxon cannot reach p<0.0125 below n=8, so per-panel tests are impossible. Pool full- and subset-mode items only if D6 passes.
- **Within-collection unit bootstraps** (YGOB, RPROT, robustness) are conditional on one source.
- **Robustness noise floor**: needed for both tools. Only OF rep2 exists; run FastOMA rep2 or demonstrate that it is deterministic.

## Concrete changes
**MUST**
- M1. Testable-event rule: each copy needs ≥2 genes from ≥2 species. Report n testable. Regression test: singletons and one-HOG score exactly 0 on event recovery.
- M2. Show every primary endpoint next to all trivial baselines and as skill (s−b)/(1−b) with the per-family b. Add S/D-balanced merge accuracy, mean(acc_S, acc_D), on which both trivial baselines score ≤~0.5. Log it in ERRATA as a post-slice secondary, not a primary.
- M3. Set membership F1 to NA for families with no complete cell and no decoy, and report its n. Rename `one_hog_per_family` to "one big HOG".
- M4. Gate: any method identical to all_singletons on >50% of units is not scored until explained.
- M5. Put CONTROL families (RPROT, K02, E05) in a controls track, out of the panel headline means.
- M6. Before final scoring, complete C9 (with an ERRATA entry stating that the Y/K slice was seen), C6 for all headline families, C3 and D6. Families without C6 become non-headline.
- M7. Adopt the plan in (d) and fix the CI labelling.

**SHOULD**
- The rate covariate, the long-branch stratum and the Y fusion perturbation.
- Check why FastOMA cohesion is 0 for K02 EF-Tu and Y03 URA1.

**COULD**
- E09 EGT, and Y02 as a cross-check collection.

## (e) The REPORT must not claim
- That ILS, fusion, alignment error or fast evolution are covered by curated families.
- "95 families". Give items, collections and units separately.
- Per-panel rankings, or any CI or p-value at n<5.
- Phenomenon-general conclusions drawn from a single collection.
- That subset mode is equivalent to full mode before D6 passes.
- That high level F1 or merge accuracy means quality, without the baseline columns.
- That the Y and K references are blind to tool outputs, unless they were frozen before any further edit.
