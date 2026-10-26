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
