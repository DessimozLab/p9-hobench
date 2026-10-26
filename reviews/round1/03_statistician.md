# Round 1, reviewer 03: Statistician

Checked on toy families (/tmp/claude-212988/sim*.py). "KIN" is 10 single-copy subfamilies × 5 species, joined by a duplication at the family root. "RIB" is one single-copy family.

1. **ST-1 (MUST): F1 is dropped, not scored 0, when a method abstains.** `metrics._prf` (l.177–181) returns `None` when tp+fp = 0 or tp+fn = 0, and `macro_average` then silently excludes that family.
   - all-singletons: every F1 is `None`, so the macro mean covers 0 families.
   - one-big-HOG on toy F1 (pseudo-orthology): Teleostei has tp=0, fp=1, so its F1 is `None` and `level_mean_f1 = 1.0`.

   *Fix:* F1 = 2tp/(2tp+fp+fn) whenever tp+fp+fn > 0.

2. **ST-2 (MUST): The code does not collapse old levels to ∞, as §2 requires.** `predicted_merge_matrix` does not collapse predicted levels above the family level. So a correct HOG at the family level that is merged at the panel root is scored "−inf, over-merged". OrthoFinder N0/orthogroups always do this.

   *Fix:* set Mpred = INF wherever pred ⊄ family level. Document that ∞ means "not together at or below the family level".

3. **ST-3 (MUST): Never merging wins on families with a root duplication.** Pairs across the root duplication grow quadratically. On KIN, all-singletons scores:
   - merge accuracy 0.918
   - gene_mean_frac_correct 0.918
   - paralog_f1 0.947
   - **event recovery 1.0**

  

   *Fixes:*
   - (a) Make the headline merge accuracy **node-balanced**: average, over the asserted reference nodes, the accuracy of the pairs whose LCA is that node. Keep the pair-pooled value as secondary.
   - (b) Make event recovery **two-sided**: the copies must be separated at D, and ≥ θ of each copy's reference-together pairs must be together at D. For ∞ events, only the second condition applies.
   - (c) Report a per-family **baseline-normalized skill**: (s − b) / (1 − b), where b = max(singletons, big-HOG, mean relabel over 10 seeds).

4. **ST-4 (MUST): Precision cannot see over-merging with non-members.** One-big-HOG scores 1.0 on every pair metric for RIB. `contamination` (l.433) divides by `hog_size`, which includes unlabelled genes, so a giant HOG in an 85k proteome scores ≈ 0. The value also differs between subset and full mode. The "declared-complete" precision promised in §10 is not implemented.

   *Fixes:*
   - Add `complete_species` to `family.yaml`.
   - Add **membership P/R/F1** at the family level, counting non-members from complete species in the cover HOGs as false positives.
   - Compute contamination as contam/(own+contam), and report unlabelled genes separately.

5. **ST-5 (MUST): Genome-scale collections must be split into units.** YGOB as one family is about 55k genes, so the n² int64 matrices (about 25 GB each) will exhaust memory. It would also weigh the same as one globin family.

   *Fixes:*
   - Split YGOB pillars, the 3R pairs and the *X. laevis* L/S pairs into units that carry a `collection` id.
   - Aggregate in two stages: the mean over units with a cluster-bootstrap CI, then each collection counts as one unit in the macro mean.
   - `events_by_tag` (l.522–532) pools events across families with no CI, so YGOB would dominate WGD. Weight it by family or collection and bootstrap by cluster.

6. **ST-6 (MUST): Paired inference and CI rules.**
   - (a) `bootstrap_ci` resamples a filtered value vector, so it cannot pair methods. Resample **family IDs**, using the same indices for every method, on the common family set.
   - (b) Implement the paired comparison that §6 promises: a paired bootstrap of the mean difference plus a Wilcoxon signed-rank test, stratified by panel.
   - (c) For n = 1, the current CI has zero width. Give no CI below n = 5 (report the range), and show n everywhere.
   - (d) Name the primary endpoints in §6: node-balanced merge accuracy, pooled level F1, two-sided event recovery and membership F1. Apply Holm correction to the OrthoFinder vs FastOMA comparisons on these. Label everything else exploratory.

7. **ST-7 (SHOULD): `split_excess` is signed and unnormalized.** Opposite errors in different families cancel, and big families dominate. *Fix:* report over-split and over-merge separately, as fractions of members.

8. **ST-8 (SHOULD): Phenomenon effects are underpowered per family tag.** With about 30 families there are only 3–8 per tag, and tags are confounded with panel. *Fix:* estimate the effects within families, e.g. accuracy on pairs split by WGD-tagged D nodes vs pairs split by S nodes. Use a cluster bootstrap by family.

9. **ST-9 (SHOULD): Robustness (E2) needs replicates and a noise floor.**
   - Run the random perturbations with SEEDS[0:3].
   - Run the clean baseline twice (different thread counts) to measure tool nondeterminism.
   - Compute Δ on matched asserted pairs, split into affected and unaffected genes, with a paired CI. Call a method "robust" only if |Δ| is within the noise.
   - Fix the scoring of fragments before running: each fragment inherits its parent's position, and both are scored.

10. **ST-10 (SHOULD): Difficulty is confounded.** With 2 methods, the mean error mostly reflects size, depth and the share of ∞ pairs.
    - *Fix:* define difficulty from the skill score in ST-3c.
    - Report Spearman ρ against the descriptors, with permutation p-values.
    - Measure agreement with AMI rather than ARI, counting unplaced genes as singletons.

11. **ST-11 (SHOULD): Overlapping CIs do not establish equivalence.** The mode check in §7 ("agree within CI") relies on this. *Fix:* compute paired per-family full − subset differences, and check them against a pre-set margin (|Δ| ≤ 0.03) with TOST or a 90% paired CI.

12. **ST-12 (COULD): Smaller items.**
    - Ties for the modal merge level (l.331) depend on insertion order; break them deterministically.
     - Note that within-family "orthology P/R" is not QfO-comparable, since its negatives are only paralogs.
