# Panel P / subset: evaluation summary

hobench git 0b5dd0a216 (dirty=True); 4 families (+0 excluded: members absent from subset) (0 control-track, 1 non-headline); families tree hash 8c1d8bb8a7fc.
Descriptive only (mean over items; a collection is one item; no CIs). Headline = tracks family+collection, headline-eligible; control track (K01 RPROT, K02, E05) reported separately. Trivial baselines are rows of every table; skill = (s-b)/(1-b), b = per-family max(singletons, one big HOG, mean relabel). membership F1 is NA for families without complete cells or decoys. Confirmatory test: scripts/confirmatory.py.

## scope: headline

| method | variant | level_bcubed_f1 (n) | skill | merge_accuracy_node (n) | skill | event_recovery_0.5 (n) | skill | membership_f1 (n) | skill | S/D-bal merge acc (post-hoc) |
|---|---|---|---|---|---|---|---|---|---|---|
| positive_control | reference | 1.000 (3) | 1.000 | 1.000 (3) | 1.000 | 1.000 (2) | 1.000 | 1.000 (1) | 1.000 | 1.000 |
| all_singletons |  | 0.941 (3) | -0.261 | 0.550 (3) | -0.250 | 0.000 (2) | 0.000 | 0.000 (1) | -2.296 | 0.250 |
| one_big_hog |  | 0.148 (3) | -12.667 | 0.233 (3) | -0.167 | 0.000 (2) | 0.000 | 0.697 (1) | 0.000 | 0.500 |
| random_relabel | mean10 | 0.908 (3) | -0.455 | 0.616 (3) | -0.035 | 0.000 (2) | 0.000 | 0.395 (1) | -0.994 | 0.363 |
| fastoma | default_declared | 0.997 (3) | 0.918 | 0.867 (3) | 0.667 | 0.667 (2) | 0.667 | 0.018 (1) | -2.235 | 0.688 |
| fastoma | default_late | 0.997 (3) | 0.918 | 0.867 (3) | 0.667 | 0.667 (2) | 0.667 | 0.018 (1) | -2.235 | 0.688 |
| of_flat | default | 0.481 (3) | -7.990 | 0.579 (3) | -0.039 | 0.167 (2) | 0.167 | 0.364 (1) | -1.096 | 0.550 |
| orthofinder | default | 0.994 (3) | 0.920 | 0.907 (3) | 0.767 | 0.667 (2) | 0.667 | 0.364 (1) | -1.096 | 0.842 |

random relabel, min-max of item-mean over 10 seeds: level_bcubed_f1 0.901446-0.912248; merge_accuracy_node 0.607556-0.629788; event_recovery_0.5 0-0; membership_f1 0.373164-0.419816

## scope: control

| method | variant | level_bcubed_f1 (n) | skill | merge_accuracy_node (n) | skill | event_recovery_0.5 (n) | skill | membership_f1 (n) | skill | S/D-bal merge acc (post-hoc) |
|---|---|---|---|---|---|---|---|---|---|---|
| positive_control | reference | NA (0) | NA | NA (0) | NA | NA (0) | NA | NA (0) | NA | NA |
| all_singletons |  | NA (0) | NA | NA (0) | NA | NA (0) | NA | NA (0) | NA | NA |
| one_big_hog |  | NA (0) | NA | NA (0) | NA | NA (0) | NA | NA (0) | NA | NA |
| random_relabel | mean10 | NA (0) | NA | NA (0) | NA | NA (0) | NA | NA (0) | NA | NA |
| fastoma | default_declared | NA (0) | NA | NA (0) | NA | NA (0) | NA | NA (0) | NA | NA |
| fastoma | default_late | NA (0) | NA | NA (0) | NA | NA (0) | NA | NA (0) | NA | NA |
| of_flat | default | NA (0) | NA | NA (0) | NA | NA (0) | NA | NA (0) | NA | NA |
| orthofinder | default | NA (0) | NA | NA (0) | NA | NA (0) | NA | NA (0) | NA | NA |

random relabel, min-max of item-mean over 10 seeds: 

## Sanity gate (tool identical to all-singletons on >50% of units)

- fastoma default_declared: 0.0% of units identical to all-singletons
- fastoma default_late: 0.0% of units identical to all-singletons
- of_flat default: 0.0% of units identical to all-singletons
- orthofinder default: 0.0% of units identical to all-singletons

## Positive control / degenerate baselines (per-family min/mean/max)

- positive_control level_bcubed_f1: mean 1.0000, min 1.0000, max 1.0000, n=4
- positive_control merge_accuracy_node: mean 1.0000, min 1.0000, max 1.0000, n=4
- positive_control event_recovery_0.5: mean 1.0000, min 1.0000, max 1.0000, n=3
- positive_control membership_f1: mean 1.0000, min 1.0000, max 1.0000, n=2
- all_singletons level_bcubed_f1: mean 0.9491, min 0.8720, max 1.0000, n=4
- all_singletons merge_accuracy_node: mean 0.5554, min 0.2500, max 1.0000, n=4
- all_singletons event_recovery_0.5: mean 0.0000, min 0.0000, max 0.0000, n=3
- all_singletons membership_f1: mean 0.0016, min 0.0000, max 0.0031, n=2
- one_big_hog level_bcubed_f1: mean 0.1398, min 0.0452, max 0.3196, n=4
- one_big_hog merge_accuracy_node: mean 0.2464, min 0.0000, max 0.5000, n=4
- one_big_hog event_recovery_0.5: mean 0.0000, min 0.0000, max 0.0000, n=3
- one_big_hog membership_f1: mean 0.4754, min 0.2542, max 0.6966, n=2
- random_relabel level_bcubed_f1: mean 0.8898, min 0.8365, max 0.9159, n=4
- random_relabel merge_accuracy_node: mean 0.5920, min 0.3984, max 0.9838, n=4
- random_relabel event_recovery_0.5: mean 0.0000, min 0.0000, max 0.0000, n=3
- random_relabel membership_f1: mean 0.2961, min 0.1972, max 0.3949, n=2

POSITIVE CONTROL: PERFECT on all primary endpoints.
