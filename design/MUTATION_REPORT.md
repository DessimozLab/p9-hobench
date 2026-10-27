# Mutation report (scripts/mutation_check.py)

Each mutant is applied to a temp copy of `hobench/`; `pytest -q tests` is run on the copy. Killed = suite fails.

| Mutant | Killed before | Killed by (before) | Killed after | Note |
|---|---|---|---|---|
| a sign flip of signed level error | yes | test_interval_correct_and_error (+5) | yes | |
| b swap S/D in `reference_matrices` | yes | 24 tests | yes | |
| c ABOVE_RANGE/INF treated as finite root | yes, but only incidentally | `test_large_family_is_fast` alone (a performance test); AC-1 tests did NOT catch it | yes (+ `test_nested_duplication_under_root_duplication_has_infinite_upper_bound`) | weak kill -> new test. Note `_parent_or_inf(ABOVE)` is reachable only for a D nested directly under a root D with unbounded Dhi (root D itself is forced to INF at reference.py:508). |
| d interval check `lo<=x<=hi` -> `x==lo` | yes | 4 tests | yes | |
| e event testability rule removed | yes | test_event_testability (sole) | yes | |
| f F1 dropped on empty P/R denominator | yes | test_pseudo_orthology_one_big_hog (sole) | yes | |
| g OrthoXML dup placement default declared->late | yes | test_paralog_group_dup_placement_declared_vs_late (+1) | yes | |
| h flat projection disabled | yes | test_flat_tsv_projection (sole) | yes | |
| i membership FP counts other-family genes outside complete cells | yes | 3 tests | yes | |
| j bootstrap pairing broken (b misaligned) | yes | test_paired_comparison_and_holm (sole) | yes | |

No mutant survived the original suite; one (c) was killed only by an unrelated test and now has a direct test in
`tests/test_mutation_gaps.py` (derivation in its docstring; passes on real code, fails on the mutant).
Untested caveat: mutants were chosen from the brief; equivalent-mutant analysis not done beyond (c).

## Control gate (scripts/controls_gate.py; families Y03_URA1_HGT, E02_ARF, V03_Hox, each under its own panel tree)

| Family | positive | all-singletons recovered | one-big-HOG recovered |
|---|---|---|---|
| Y03_URA1_HGT | perfect (1/1 events) | 0/1 | 0/1 |
| E02_ARF | perfect (5/5) | 0/5 | 3/5 |
| V03_Hox | perfect (61/61) | 0/61 | 34/61 |

Finding (deviation from "one-big-HOG = 0"): every one-big-HOG recovery is a duplication placed AT a species leaf
(in-paralogs). One-big puts same-species pairs at parent(species), which equals the reference merge, and single-gene
copies make cohesion vacuous, so these events count as recovered (and testable). Gate therefore asserts: singletons = 0
everywhere; one-big recovers nothing except species-leaf duplications. Decide whether leaf-level duplications with
single-gene copies should be non-testable or whether one-big-HOG legitimately scores them (it does get those pairs right).
