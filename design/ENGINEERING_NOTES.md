# Engineering notes

## hobench core package (v0.2: review round 1 implemented)

Code: `hobench/`. Tests: `tests/` (toy data with hand-derived answers, plus excerpts of real FastOMA 0.5.0 and
OrthoFinder 3.1.5 output in `tests/fixtures/`). Run `python -m pytest tests/`.

Module map (public API):

| Module | Contents |
|---|---|
| `speciestree` | `SpeciesTree` |
| `reference` | `load_family(ies)`, `write_family`, `validate_family`, `validate_panel`, `reconcile`, `tree_from_classification` |
| `predictions` | `Prediction`, `read_orthoxml`, `read_orthofinder`, `read_tsv` |
| `metrics` | `evaluate_family`, `family_integrity`, `reference_matrices`, `representative_matrix`, `interval_correct`, `interval_error`, `panel_context`, `family_scalars`, `SCALARS`, `PRIMARY` |
| `stats` | `aggregate`, `mean_ci`, `bootstrap_indices`, `compare_methods`, `holm`, `ami`, `agreement` |
| `controls` | `reference_prediction`, `all_singletons`, `one_big_hog`, `random_relabel`, `swap_leaves`, `perturbed_reference_prediction` |
| `panel` | `load_panel`, `validate_members` |
| `freeze` | `python -m hobench.freeze` |
| `evaluate` | `python -m hobench.evaluate` |

Every module docstring states its semantics.

### Reference encoding (A1, A2, A4)
- **Merge intervals.** Reference merge levels are intervals `[lo, hi]`; `hi` may be ∞. For a pair whose LCA node in
  the reference tree is:
  - **S:** `lo = hi = lca(species)`.
  - **D:** `[parent(Dlo), parent(Dhi)]`, where the duplication level is `lvl` (exact) or `lvl_lo`/`lvl_hi`.
    `lvl_hi=ABOVE_RANGE` gives ∞.
  - **H:** `lca(donor, recipient clade)`, or ∞ when `donor=OUTSIDE`.
  - **U:** not asserted.
- **Default duplication interval.** It runs from the LCA mapping to the *tightest ancestor bound*:
  - Below a speciation, the bound is the child (toward the D node) of `lca(mapping, closest sibling lineage)`, strictly
    below the speciation. For a binary S consistent with the species tree this is the S node's mapping.
  - Below a D or H node, the bound is that node's `lvl_hi`.
  - At the top of the tree it is `ABOVE_RANGE`. With `root_single_copy: true` it is strictly below `root_level`, or
    below the parent of the root node's level.
- **Root duplications.** A duplication with no constraining non-U ancestor always gets `hi = ∞`, even with an explicit
  `lvl`, unless `root_single_copy` is set (SYNTHESIS A1 / AC-1). Never merging the copies is then also correct.
- **Default family level.** If `root_level` is absent, the family level is the oldest of: the root's level, every D/H
  `lvl_lo`, and every finite D/H merge `hi`. This was found on YGOB units: a U root over an interval D otherwise scored
  integrity at a level where nothing is asserted.
- **H nodes.** Exactly one child carries `recipient=1`, and `donor=<level|OUTSIDE>` is required. The donor must not be
  nested with the recipient clade. Pairs between transferred genes and genes outside the recipient clade are excluded
  from within-family orthology. If the donor relation has no source, curators use U (AC-10 fallback).
- **Allopolyploid yeast WGD (A2).** Expressed as `D[lvl_lo=PostWGD, lvl_hi=ZT_WGD]`, with ZT-lineage genes under a U.
  Tested in `test_reference.py` and `test_metrics.py`.
- **Collections.** Layout is `families/<ID>/collection.yaml` plus `units/<unit>/`. `load_families` loads each unit as a
  family and inherits `collection` and `panel` from `collection.yaml`; a conflicting value is an error.
- **Validation errors.** All of these are errors, not warnings:
  - unknown or missing `ev`;
  - `lvl` together with `lvl_lo`/`lvl_hi`;
  - `lvl_hi` younger than `lvl_lo`;
  - bound violations;
  - S children overlapping in the species tree (taking their youngest level into account);
  - level keys on S or U nodes;
  - bad `distinct_from`, `complete_species` or `root_single_copy` values;
  - decoys that are also members.

### Scoring (B1–B3, B5–B7, B12, B13)
- **Correctness.** A predicted merge level `p` is correct iff `lo ⊆ p ⊆ hi`. The signed error is `depth(bound) −
  depth(p)` to the nearest bound: `+` means over-split (predicted too old), `−` means over-merged (too young). Exact
  and interval pairs are reported separately (`accuracy_exact_pairs`, `accuracy_interval_pairs`).
- **Primary endpoints** (`metrics.PRIMARY`, B8):
  - `level_bcubed_f1`: mean over levels ≤ the family level of gene-balanced B-cubed F1. Self is counted in both
    clusters. A pair is skipped at the levels that lie inside its interval.
  - `merge_accuracy_node`: mean over asserted reference nodes of the accuracy of the pairs separated at that node.
  - `event_recovery_0.5`: two-sided event recovery (θ = 0.8 is reported as a sensitivity).
  - `membership_f1`.
- **Event recovery.** An event (D or H) is recovered at θ iff:
  - ≥ θ of its cross-copy pairs are inside the merge interval, i.e. apart below it (`separation`) and merged by its
    upper bound (`merge_side`); and
  - every copy keeps ≥ θ of its reference-together pairs together at the duplication level (`cohesion`).

  Single-gene copies make cohesion vacuous. Separation alone is therefore *not* sufficient, otherwise all-singletons
  would "recover" F1's duplication. The modal predicted level breaks ties by smallest |error|, then by level index.
- **B1.** F1 = 2tp/(2tp+fp+fn) whenever anything is counted. It is never dropped.
- **Spurious duplications.** For each maximal reference clade containing only speciation nodes, the count at each
  level L is (#predicted HOGs among the clade's genes below L) − 1, summed over levels. The rate divides by the number
  of opportunities.
- **Integrity at the family level** (B5–B7, OrthoBench and OF3 definitions; formulas in the `family_integrity`
  docstring):
  - **RefOGs** are the reference-derived HOGs at the family level. Unplaced members join if there is a single RefOG;
    otherwise they are excluded and counted.
  - **PredOGs** are each member's predicted HOG at the highest level ≤ the family level. A singleton counts as missing,
    except for a single-gene RefOG, where being alone is exact.
  - **False positives** are only *judgeable* genes:
    - members of another RefOG of the same family;
    - declared decoys;
    - members of a `distinct_from` family, when the HOG's level ⊆ `at_level` (`ALL` means at any level);
    - non-members from `complete_species`.

    Other genes, including undeclared other reference families, are *unlabelled* and reported only (B6).
  - **Measures reported:**
    - membership P/R/F1 (OF3 overlap-weighted per RefOG, then size-weighted);
    - OrthoBench pairwise P/R/F (`benchmark.py`, q_even);
    - % fission, % fusion, % exact, % missing genes;
    - normalised entropy;
    - contamination = contam/(own+contam);
    - unsigned `frac_members_split_off` and `frac_over_merged`.

    The signed `split_excess` and the "≥90% cover" rule are removed.
- **B12 and B13.** `mean_abs_depth_frac_error` is |error| / (family depth). `accuracy_tol1_*` treats ±1-level errors as
  correct (a sensitivity analysis). The orthology metric is renamed `within_family_orthology`; its negatives are
  within-family paralogs, so it is not comparable to QfO.
- **Raises (`ScoringError`).** Scoring raises when a predicted gene sits at a level its species is not below, when the
  prediction and the reference disagree on a gene's species, or when a complete species was not counted by the parser.

### Aggregation and inference (B4, B9, B10)
- **Items.** A family is one item, and a collection is one item. A collection's value is the mean over its units, with
  a unit-bootstrap CI.
- **Tracks.** `aggregate()` reports three tracks: `families` (ordinary families only), `collections` (one entry per
  collection) and `combined`. `events_by_tag` averages per item first, so a collection with many events counts once.
- **Bootstrap.** Percentile bootstrap over item IDs, stratified by panel. `bootstrap_indices` is deterministic in (IDs,
  strata, seed), so methods scored on the same items get paired resamples. No CI is given when n < 5; min and max are
  reported instead, and n is always reported.
- **Comparisons.** `compare_methods` gives the paired mean difference with a CI plus a Wilcoxon signed-rank p-value.
  `holm` adjusts across comparisons.
- **Skill.** Skill = (s − b)/(1 − b), where b = max(all-singletons, one-big-HOG, mean of the within-species-relabelled
  positive control over `SEEDS[:10]`). It is computed per family for the primary metrics (12 extra scorings per panel;
  disable with `--no-skill`).
- **Agreement.** `agreement()` computes per-level AMI on reference members, using arithmetic normalisation; it matches
  scikit-learn.
- **Scale.** Y01 YGOB (1,977 units, 26k genes): positive control 10 s, scoring 7.5 s, scoring with skill about 2.5 min,
  peak memory 375 MB. A 5,000-unit synthetic collection scores in about 15 s.

### Parsers (E1–E3)
- **OrthoXML levels (E2).** Level names in a file (`taxonId` resolved through `<taxonomy>`, or `TaxRange`) are mapped
  to our levels **by clade species set**. The clades come from the tool's tree: `--tool-species-tree`, e.g. FastOMA's
  `species_tree_checked.nwk`; else the file's `<taxonomy>`; else our tree, by name. Any name that cannot be mapped is a
  hard error. This makes FastOMA's renaming of `KE` to `internal_0` harmless. A group with no level annotation gets
  the MRCA of its members (pyHam).
- **OrthoXML duplication placement (E1).**
  - `declared` (primary): a `paralogGroup` is a duplication just below its group's declared level, and each copy gets
    its own chain of induced HOGs.
  - `late` (sensitivity, pyHam): the duplication sits at parent(MRCA(copy levels)).

  `meta['dup_placement_ambiguous']` counts the groups where the two readings differ: 37 on the 4-yeast calibration
  run. Skipped levels between non-paralog children are induced as in pyHam's `_resolve_missing_levels`. A group at the
  same level as its child, or declared at a species, is spliced into its parent.
- **OrthoXML strictness.** Strict by default. These raise:
  - a declared level that is not above all members;
  - a duplication older than its group;
  - a geneRef to an unknown gene;
  - a gene in two groups;
  - an unknown species.

  `strict=False` repairs these and records warnings; it is for sensitivity runs only.
- **OrthoFinder (E3).**
  - `N*.tsv` nodes are mapped to our levels by exact clade.
  - OrthoFinder 3.x writes no `N0.tsv`, so the root comes from `Orthogroups/Orthogroups.tsv`, and each HOG is checked
    to lie inside its OG.
  - Every one of our levels must match an OrthoFinder clade, so OrthoFinder must be run with `-s <panel tree>`.
    `--allow-induced-levels` (perturbed-tree runs only) induces missing levels from the smallest containing clade.
  - The `map:` ID normaliser is strict; `map_or_auto:` passes unmapped IDs through explicitly.

### Reproducibility and curation guards (C4, C9)
- **Member checks (`panel.validate_members`).** These are errors:
  - a gene missing from the sha256-verified pinned FASTA of its species;
  - a sequence-based mapping (`mapping_method` matching diamond, mmseqs, blast or sequence) with identity below 98%
    (`mapping_identity` is a percentage), or below 99% for XENLA, salmonids, `strict_mapping` species, families tagged
    TAND, or families with `strict_mapping: true`;
  - coverage below 0.80 (`mapping_coverage`, a fraction; it is the minimum over both sequences, or given as
    `mapping_qcov`/`mapping_tcov`);
  - unmapped or ambiguous methods.

  Any other method name is treated as an exact ID mapping and is exempt from the thresholds.
- **`evaluate`.**
  - Requires `--panel-dir` unless `--no-fasta-check` is given explicitly; the choice is recorded in the output.
  - Records the hobench git commit and dirty flag, plus sha256 hashes of the panel tree, the prediction, a hash of the
    families tree, the manifest's hash and its differences, `species.tsv`, and the tool tree.
- **`python -m hobench.freeze`.** Writes `families/MANIFEST.sha256`; `--check` verifies the files against it.

### Open items and deviations (for PI)
- **Representative grouping is a choice.** Integrity RefOGs and the positive control choose one merge level inside
  each interval: the youngest one consistent with descendant duplications. U pairs are grouped at the oldest level
  implied inside the U. This choice cannot change pair scores, but it does define the RefOGs when an interval
  straddles the family level.
- **Event recovery is stricter than ST-3b.** Cohesion alone is vacuous for single-gene copies, so recovery also
  requires the cross-copy pairs to fall inside the merge interval.
- **Not done:**
  - B11, within-family phenomenon effects. The per-node accuracy (`nodes` in the JSON) makes it computable.
  - B13's second source for OrthoFinder orthology (`Orthologues/`). This is COULD.
  - "Collapsed duplications" in the dose-response control.
  - Myr-based tolerance (B12). This is COULD.

### Change log
- v0.1: initial package. It used pyHam-late duplication placement, name-based level mapping, lenient parsers, exact
  merge levels with a ∞ collapse at the family root, `split_excess`, and an unpaired bootstrap.
- v0.2: review round 1, as described above. The signed-error labels in core_spec.md were corrected: + means
  over-split.
