# Round 1 — Algorithms & compute (AC)

Evidence scripts: /tmp/claude-212988/*.py.

1. **AC-1 MUST — Root-duplication "∞ collapse" penalises correct predictions.** `reference.reconcile`: without `root_level`, `family_level` = root D's `lvl`, so `merge = parent(lvl)` becomes INF. Counterexample: `((DANRE_a,ORYLA_a)S,(DANRE_b,ORYLA_b)S)D[lvl=Teleostei]`. The correct prediction (copies apart at Teleostei, one HOG at Euteleostomi) scores merge accuracy 0.33, and the TGD event scores frac_correct 0 / "-inf". Never merging the copies scores 1.0. This hits the headline families (3R/4R ohnologs, *X. laevis* L/S, 2R) whenever no outgroup member is curated. **Fix:** make reference merges *intervals* [lo, hi]. For cross-copy pairs of a root D at level X: [parent(X), ∞]. It narrows only if the curator asserts a single ancestral gene at `root_level` (new field). Correct iff the predicted merge lies in it; signed error = distance to the nearest bound. ∞ alone only if X is the panel root; update F4 test and spec.

2. **AC-2 MUST — Ambiguous duplication levels (losses, missing `lvl`).** A D without `lvl` defaults to its LCA mapping, the *youngest* possible level. Example: WGD copies retained only in YEAST get merge parent(YEAST), so a tool placing the duplication at the WGD is scored over-split. **Fix:** use the same interval mechanism. `lvl_min` = `lvl`, else the mapping. `lvl_max` = an explicit `lvl_max`, else the tightest ancestor bound already computed in `validate_family` (child-toward of the S ancestor, the D ancestor's level, or `root_level`/panel root). Report interval pairs separately; never score them as exact.

3. **AC-3 MUST — The OrthoXML duplication placement contradicts FastOMA's own output.** ENGINEERING_NOTES says FastOMA writes TaxRange = MRCA of members. The calibration file disproves this: 49 orthologGroups are declared at Saccharomycetaceae with members in only one subclade. In 43 of about 580 paralogGroups (7%), pyHam's "late" placement M = parent(MRCA(copies)) falls strictly *below* the enclosing group. E.g. two YEAST copies in a Saccharomycetaceae group: FastOMA did not merge them at the WGD level, but the late rule does, penalising FastOMA on WGD ohnologs with losses (the Y panel's core phenomenon). **Fix:** for FastOMA, place the duplication at the enclosing group's declared level (primary). pyHam-late = sensitivity; count `dup_placement_ambiguous` per family; add a toy test.

4. **AC-4 MUST — FastOMA silently renames internal nodes.** `FastOMA/check_input.py:add_internal_node_prune` renames names that are short (<3 chars), duplicated or contain characters outside `[A-Za-z0-9_-]`. In calibration `KE` became `internal_0`, and 4,200 TaxRanges fell back to MRCA inference silently. **Fix:** (a) have `SpeciesTree` / the panel builder enforce `^[A-Za-z0-9_-]{3,}$` for level names. (b) Have `read_orthoxml` take `species_tree_checked.nwk`, build `level_map` *by clade species set* (as the OrthoFinder parser does), and error if `taxrange_unknown > 0`. (c) Required for E2(d), where names can denote different clades.

5. **AC-5 MUST — OrthoFinder configuration.** Calibration ran without `-s` (OF inferred and rooted its own tree). Production must pass `-s <panel tree>` with leaf names equal to the FASTA basenames. OF 3.1.5 rejects trees whose root has more than 2 children (`run/species_info.py:324`), and internal polytomies are untested. Give the tools a deterministically binarised tree. In main runs assert `meta['induced_levels'] == {}`. For E2(d), score only levels shared by both trees, plus the level-free orthology P/R. N0: the code is correct; in 3.1.5 `Orthogroups.tsv` is built from `WorkingDirectory/N0.ids.tsv` (`orthologues.py:690`).

6. **AC-6 MUST — Critical path and compute budget (12 cores / 46 GB).** V/P/E subsets need curated members (≈T8); runs are planned T4–14. Plan:
   - Start Y and K full runs as soon as proteomes are available.
   - Precompute OMAmer hogmaps once for all panels: ≈1.2 M proteins at ~350 q/s with 2 processes ≈ 1 h. Pass via FastOMA `--hogmap_in`; never co-run omamer with other >10 GB jobs.
   - Freeze subsets from raw paper member lists at a fixed time; members missing from the subset count as absent.
   - Estimates (scaled from calibration, NVMe workdirs): per run, OF ≈ 0.5 h and FastOMA ≈ 0.75 h for Y/K full or V/P/E subsets ≤ 70k proteins; Y/K subset validation ≈ 1 h.
   - E2 on **Y full only**, 1 seed: (d) and (e) reuse OF `-ft`/`-b` (minutes) and the hogmaps. (a)–(c) rerun with re-search of only the modified proteomes, ≈ 4 h.
   - Total ≈ 12–13 h serial, ≈ 8 h wall with OF ∥ FastOMA (6 CPUs each).
   - **Cut:** panel D, E2 on subsets, and E3 as specified. A leakage-free OMAmer DB (OMA HDF5 + `omamer mkdb`) is infeasible. Replace E3 with one FastOMA run on a smaller prebuilt clade DB, if downloadable, labelled "DB sensitivity".

7. **AC-7 SHOULD — Subset size control.** Promiscuous domains (Rab/ARF, kinome, KRAB-ZNF, NLR→LRR-RLKs) and genome-scale sets (3R, L/S) make neighbourhoods near-full proteomes. Compute sizes now. Require ≥ 50% query coverage for neighbourhood hits and add a seeded sample of domain-only hits as decoys. Sample genome-scale families to ≤ 300 groups per panel (seed from `seeds.py`). Target ≤ 70k proteins per panel.

8. **AC-8 SHOULD — Subset validation statistic.** "Agree within CI" with ~5 families/panel has little power. Report per-pair identity of predicted merge levels (subset vs full) plus metric Δs; the bias is tool-asymmetric (OMAmer is per-protein, OF's MCL is context-dependent).

9. **AC-9 SHOULD — Noise floor.** Run one replicate of Y per tool. Compare E2 deltas to this variation.

10. **AC-10 SHOULD — HGT is not implemented.** `EVENTS = ("S","D","U")`: no H node. Encode H as U for now. Score xenologs only where well-defined (xenolog vs a native non-orthologous copy, e.g. URA1 vs URA9, via decoys); otherwise report merge levels descriptively.

11. **AC-11 COULD — Run manifest.** Record command lines, OMAmer DB checksum and the species-tree checksum used by tools and reference; assert equality at scoring.
