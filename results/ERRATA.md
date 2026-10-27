# Results errata

## 2026-10-27 — E1 (testable events) and E2 (flat baseline projection)

Code: hobench v0.2.1. Rationale and rules are in design/ENGINEERING_NOTES.md ("Errata after the Y/K vertical slice").
- **E1.** Event recovery now counts *testable* events only; the rule depends on the reference alone. An event is testable
  iff its cross-copy merge interval is bounded above within the family level, or one of its copies has a
  reference-together pair. Y has 202 testable and 1,390 untestable events; K has 1 testable and 0 untestable.
- **E2.** of_flat orthogroups are now projected to every level as "no duplications". Before, they were read as root-level
  HOGs only, which behaved like singletons below the root.

Primary endpoints, before → after: mean over items, a collection counts as one item, scope = all families.
Re-run with `scripts/evaluate_panel.py --panel {Y,K} --mode full`. Everything else is unchanged.

### Panel Y (items)

| method | variant | level_bcubed_f1 | merge_accuracy_node | event_recovery_0.5 | membership_f1 |
|---|---|---|---|---|---|
| positive_control | reference | 1.000 | 1.000 | 1.000 | 1.000 |
| all_singletons |  | 0.695 | 0.128 | 0.431 → **0.000** | 0.001 |
| fastoma | default_declared | 0.920 | 0.942 | 0.306 → **0.134** | 0.918 |
| fastoma | default_late | 0.920 | 0.942 | 0.305 → **0.139** | 0.918 |
| of_flat | default | 0.686 → **0.790** | 0.057 → **0.872** | 0.431 → **0.007** | 0.337 → **0.833** |
| of_flat | dendroblast | 0.686 → **0.790** | 0.057 → **0.872** | 0.431 → **0.010** | 0.337 → **0.833** |
| one_hog_per_family |  | 0.790 | 0.872 | 0.000 | 0.836 |
| orthofinder | default | 0.928 | 0.879 | 0.122 → **0.060** | 0.831 |
| orthofinder | dendroblast | 0.928 | 0.880 | 0.127 → **0.072** | 0.831 |
| random_relabel | mean10 | 0.694 | 0.129 | 0.431 → **0.000** | 0.216 |

### Panel K (items)

| method | variant | level_bcubed_f1 | merge_accuracy_node | event_recovery_0.5 | membership_f1 |
|---|---|---|---|---|---|

YGOB unit level (mean over 1,976 units), event_recovery_0.5 before → after:

| method | before | after |
|---|---|---|
| all_singletons | 0.861 | 0.000 |
| random_relabel | 0.861 | 0.000 |
| FastOMA declared | 0.613 | 0.269 |
| FastOMA late | 0.610 | 0.279 |
| OrthoFinder default | 0.244 | 0.119 |
| OrthoFinder dendroblast | 0.254 | 0.144 |
| of_flat default | 0.861 | 0.015 |

## E3 (evaluation layer, round-2 decisions R2-2/3/5/7; scripts/evaluate_panel.py, scripts/confirmatory.py)
Re-run of Y_full and K_full (hobench with E1 testable-event rule and E2 flat projection). Changes: baselines and skill beside every
endpoint; `one_hog_per_family` renamed `one_big_hog`; membership F1 NA without complete cells/decoys; tracks control/collection/family
(headline = family+collection only; control = K01 RPROT, K02, E05 MCM, reported separately); sanity gate vs all-singletons;
summaries descriptive (no CIs, relabel min-max); runs with manifest status != done refused; S/D-balanced merge accuracy added
as a POST-RESULTS secondary (mean of S-node and D-node accuracy; not pre-registered).
Before (e0f02f4, all items pooled) -> after (headline track), Y: OF-default event recovery 0.122 -> 0.060; FastOMA-declared 0.306 -> 0.134;
all-singletons 0.431 -> 0.000 (E1); of_flat merge acc 0.057 -> 0.872 and bcubed 0.686 -> 0.790 (E2: now equals one big HOG at root-level families).
Membership F1: Y now NA (was 0.83/0.92); K headline n=1 (K04 decoys): OF-default 0.949, FastOMA 0.897 (K pooled before: 0.94/0.78).
K: K01 RPROT, K02 moved to control track; no K headline item has asserted pairs, so K level/merge/event headline = NA.
Confirmatory (results/confirmatory.tsv, Y+K pooled, OF-default minus FastOMA-declared, Holm over 4): level n=2 p=0.50; merge n=2 p=1;
event n=2 p=1; membership n=10 p=0.45; all Holm p=1.0. n is far below what a sign-flip test can resolve (min p with n=2 is 0.5).

## E4 (PI decision, 2026-10-27 ~00:05 UTC, before any E/P/V scoring was seen)
Controls gate (scripts/controls_gate.py, design/MUTATION_REPORT.md) showed one-big-HOG "recovers" duplications placed on a
species leaf with single-gene copies (E02 3/5, V03 34/61): those events only test in-paralog grouping, which lumping
satisfies correctly. Rule: headline event recovery is computed over **internal-branch** testable events only; leaf-level
events are reported separately as "in-paralog recovery" (secondary). Y/K values unaffected (no leaf events) — to be
confirmed in the re-run.
