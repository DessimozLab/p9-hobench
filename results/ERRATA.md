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
E4 implementation (scripts/evaluate_panel.py, evaluation layer): per-event `event_branch` (leaf = event level is a species, else internal) in events.tsv;
`event_recovery_*` now over internal-branch testable events; `in_paralog_recovery_0.5` (leaf events) is a secondary column.
Y/K full re-run: all four primary endpoint means identical to before in every scope/method (0 differences); no leaf events in Y/K.
Subset mode (`--mode subset`): only families whose members are all in panels/<P>/subset/members.tsv are scored; excluded ids in results/<P>_subset/excluded_families.txt.

## E5 (PI decision, 2026-10-27 ~01:00 UTC; found by the analysis agent in P02 error anatomy)
Membership F1 built RefOGs from the family-level reference partition; when the reference root is **U** (unresolved) the
whole family became one RefOG, i.e. the metric demanded an assertion the reference does not make (DESIGN §2: only
asserted relations are scored). Rule: under a U root, each child subtree is its own RefOG (as under a D root); under a
U node below the root likewise if the U node is at the family level. Before/after values for all panels recorded below
by the re-run. This changes scores for every method equally in definition; it was triggered by a FastOMA case (P02).

### E5 implementation and before → after
Implemented in commits d2c56a6 and the follow-up refinement. A U node at the top of the tree (or at the family level) does
not join its children into one RefOG. A U node also asserts no separation, so another RefOG's members count as false
positives only if the reference asserts the two apart at the family level. A singleton RefOG is scored only if it is
asserted apart from every other RefOG, and a family with no eligible RefOG gets membership = None (untestable,
excluded). Applying the literal rule alone was degenerate on K, where U roots sit directly over single genes:
all-singletons membership became 1.000 and FastOMA/OrthoFinder fell to 0.35/0.26 on K full. Only membership_f1
(n in parentheses = families scored) depends on RefOGs. Re-run with `scripts/evaluate_panel.py`; summary.tsv means.

#### Y_full (scope headline)

| method | variant | membership_f1 before | after |
|---|---|---|---|
| positive_control | reference | – (n=0) | – (n=0) |
| all_singletons |  | – (n=0) | – (n=0) |
| fastoma | default_declared | – (n=0) | – (n=0) |
| fastoma | default_late | – (n=0) | – (n=0) |
| of_flat | default | – (n=0) | – (n=0) |
| of_flat | dendroblast | – (n=0) | – (n=0) |
| one_big_hog |  | – (n=0) | – (n=0) |
| orthofinder | default | – (n=0) | – (n=0) |
| orthofinder | dendroblast | – (n=0) | – (n=0) |
| random_relabel | mean10 | – (n=0) | – (n=0) |
| random_relabel | seeds10 | – (n=) | – (n=) |

#### K_full (scope headline)

| method | variant | membership_f1 before | after |
|---|---|---|---|
| positive_control | reference | 1.000 (n=10) | – (n=0) |
| all_singletons |  | 0.000 (n=10) | – (n=0) |
| fastoma | default_declared | 0.897 (n=10) | – (n=0) |
| fastoma | default_late | 0.897 (n=10) | – (n=0) |
| of_flat | default | 0.949 (n=10) | – (n=0) |
| of_flat | dendroblast | 0.956 (n=10) | – (n=0) |
| one_big_hog |  | 0.753 (n=10) | – (n=0) |
| orthofinder | default | 0.949 (n=10) | – (n=0) |
| orthofinder | dendroblast | 0.956 (n=10) | – (n=0) |
| random_relabel | mean10 | 0.221 (n=10) | – (n=0) |

#### E_subset (scope headline)

| method | variant | membership_f1 before | after |
|---|---|---|---|
| positive_control | reference | 1.000 (n=3) | 1.000 (n=3) |
| all_singletons |  | 0.075 (n=3) | 0.000 (n=3) |
| fastoma | default_declared | 0.574 (n=3) | 0.563 (n=3) |
| fastoma | default_late | 0.574 (n=3) | 0.563 (n=3) |
| of_flat | default | 0.289 (n=3) | 0.152 (n=3) |
| one_big_hog |  | 0.049 (n=3) | 0.014 (n=3) |
| orthofinder | default | 0.741 (n=3) | 0.870 (n=3) |
| random_relabel | mean10 | 0.120 (n=3) | 0.180 (n=3) |
| random_relabel | seeds10 | 0.120 (n=) | 0.180 (n=) |

#### P_subset (scope headline)

| method | variant | membership_f1 before | after |
|---|---|---|---|
| positive_control | reference | 1.000 (n=1) | 1.000 (n=1) |
| all_singletons |  | 0.000 (n=1) | 0.000 (n=1) |
| fastoma | default_declared | 0.018 (n=1) | 0.049 (n=1) |
| fastoma | default_late | 0.018 (n=1) | 0.049 (n=1) |
| of_flat | default | 0.364 (n=1) | 0.493 (n=1) |
| one_big_hog |  | 0.697 (n=1) | 0.550 (n=1) |
| orthofinder | default | 0.364 (n=1) | 0.493 (n=1) |
| random_relabel | mean10 | 0.395 (n=1) | 0.277 (n=1) |
| random_relabel | seeds10 | 0.395 (n=) | 0.277 (n=) |

#### K_subset (scope headline)

| method | variant | membership_f1 before | after |
|---|---|---|---|
| positive_control | reference | 1.000 (n=10) | – (n=0) |
| all_singletons |  | 0.000 (n=10) | – (n=0) |
| fastoma | default_declared | 0.910 (n=10) | – (n=0) |
| fastoma | default_late | 0.910 (n=10) | – (n=0) |
| of_flat | default | – (n=) | – (n=0) |
| one_big_hog |  | 0.753 (n=10) | – (n=0) |
| orthofinder | default | 0.953 (n=10) | – (n=0) |
| random_relabel | mean10 | 0.221 (n=10) | – (n=0) |

