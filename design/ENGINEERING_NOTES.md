# Engineering notes

## hobench core package (metrics library, v0.1)

Code: `hobench/`, tests: `tests/` (toy data with hand-derived answers + excerpts of real
FastOMA 0.5.0 / OrthoFinder 3.1.5 output in `tests/fixtures/`). Run: `python -m pytest tests/`.

### Representation
- Levels are species-tree node indices; merge level `INF = -1` (never in one HOG),
  `NOT_ASSERTED = -2` (pair below a `U` node). Species leaves are levels too (a lineage-specific
  duplication can sit at a species), but predicted HOGs are stored for internal levels only.
- Reference merge levels are written per reference-tree node over products of child leaf
  sets (O(n²) array writes); predicted merge levels are filled level by level, deepest first,
  over predicted HOG blocks. A 2000-member family (2M pairs) evaluates in ~1–2 s.
- Parsers keep memberships only for reference members and decoys (`keep_genes`) but count
  all genes in HOG sizes (for contamination). Memory stays small on whole-proteome output.

### OrthoXML level induction (follows pyHam 2.0, `parsers.py`)
1. Level of an `orthologGroup`: `taxonId` → `<taxonomy>` name (OrthoXML 0.5), else the
   `TaxRange` property, else the MRCA of its children (one child level → its parent).
2. `paralogGroup` = one duplication (directly nested paralogGroups are merged into one, as in
   pyHam). Duplication level `M = parent(MRCA(levels of the copies))`. If `M` is below the
   enclosing group, an induced HOG at `M` holds the copies. Each copy gets its own chain of
   induced HOGs down to its level, so copies are never together below `M`. This places the
   duplication **as late as the file allows**. FastOMA writes `TaxRange` = MRCA of members, not
   the level at which the subHOG was inferred, so that information is lost in the file. We read
   the file as written, like pyHam.
3. Skipped levels for non-paralog children: children lying strictly below the same child
   level `C` share one induced HOG at `C` (no duplication was recorded between them). Children
   sitting exactly at `C` stay separate HOGs, as in pyHam's `_resolve_missing_levels`.
4. A child group at the same level as its parent, or a group declared at a species, is
   spliced into its parent.
5. Deviations from pyHam (lenient by default, `strict=True` restores pyHam's errors):
   - A declared level that is not above all members falls back to the MRCA, with a warning.
   - A duplication predating its enclosing group is placed just below that group.
   - A gene in two groups keeps its first assignment and is counted in
     `meta['membership_conflicts']`.
   - Reading the pyHam source (not tested by running pyHam), its `_resolve_missing_levels`
     appears able to group two copies of one duplication that lie strictly below the same
     intermediate level, which would undo the duplication. We always keep copies apart.

### OrthoFinder
- `N*.tsv` nodes are mapped to our levels by **exact clade species set** using
  `Species_Tree/SpeciesTree_rooted_node_labels.txt`.
- OrthoFinder 3.1.5 writes **no `N0.tsv`**, so root-level groups are taken from
  `Orthogroups/Orthogroups.tsv`. Every N* HOG is checked to lie inside the OG named in its
  `OG` column (`meta['n_hogs_not_nested_in_og']`; 0 on the 4-yeast test run).
- Our levels with no matching OrthoFinder clade (different topology or rooting) are induced
  from the smallest OrthoFinder clade containing them, restricted to their species. This
  assumes no duplication in between. Induced levels are listed in `meta['induced_levels']`.
- Gene IDs: the default normaliser `auto` maps `sp|ACC|NAME` / `tr|ACC|NAME` to `ACC` and
  takes the first token of FASTA headers; `map:<tsv>` handles anything else.

### Decisions and deviations from core_spec.md (for PI review)
- **Sign labels.** We keep the spec's sign: error = depth(ref) − depth(pred), so + means
  predicted too old. The spec's parenthetical labels are swapped, though. A predicted merge
  that is too old means genes are kept apart for too long, which is **over-splitting**. Too
  young means **over-merging** (pseudo-orthologs called orthologs). The code and outputs use
  `frac_over_split` (+, including `+inf` = never merged) and `frac_over_merged` (−, including
  `-inf`). Pairs where both merge levels are ∞ count as exact agreement ("0").
- **Family level.** Equal to `root_level` in family.yaml if given, else the level of the
  reference-tree root (an explicit D `lvl`, or the LCA mapping). A duplication whose parent
  level lies above the family level gives merge = ∞. This is the spec's "no parent level
  inside the family" case.
- **Default D level.** A D node without `lvl` sits at its LCA mapping (the youngest
  consistent level). Events report `level_explicit`.
- **Extra validation.** Children of an S node must map to disjoint clades, taking explicit
  duplication levels into account. A D `lvl` must lie strictly below an S ancestor's mapping
  and within a D ancestor's `lvl`. Without these checks, "same HOG at L" stops being
  transitive. The positive control reports `inconsistent_pairs` as a further check.
- **Metric 2** (per-level F1) uses only levels at or below the family level. Pairs above it
  add no new reference information; ∞-pair mistakes still show up in metric 1. "Same HOG at
  L" is read directly from the predicted HOGs, so non-nested predictions are scored as given.
- **Metric 3** uses cross-species pairs only. Within a reference family every non-ortholog
  pair is a paralog, so the paralog class is the complement. A prediction that never merges a
  pair counts as "not ortholog".
- **Metric 4** reports `frac_correct`, `recovered` (≥ 0.5), the modal predicted merge level
  and `misplacement` (signed level distance from the modal predicted merge to the reference
  merge), plus fractions over-split and over-merged.
- **Metric 5.** Each member's "root-level HOG" is its HOG at the highest level ≤ the family
  level, and members with no HOG count as singletons. Contamination counts genes of *other*
  reference families plus declared `decoys`; unlabelled genes are reported separately. Some
  families legitimately have more than one HOG at the family level (a duplication at the
  family root). We therefore also compute the reference's own expected values
  (`ref_n_hogs_cover`, `ref_coverage`) and report **`split_excess`** = predicted − expected.
- **U nodes in the positive control.** Genes under a `U` node are grouped at the oldest
  level implied inside it. The pairs remain unasserted, so this choice cannot affect scores.
- **Random relabelling** permutes gene labels within species. A global permutation would put
  genes at levels their species are not below.
- **Dose-response.** k random swaps of two leaves of the same species anywhere in the panel
  (so genes can move between families, which also exercises contamination). "Collapsed
  duplications" is not implemented yet.
- **Bootstrap.** Percentile CI of the macro mean, resampling families. The seed comes from
  `hobench/seeds.py::SEEDS` (fixed list; append, never edit).
