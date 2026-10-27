# Panel P / subset: evaluation summary

hobench git 3b738caa34 (dirty=True); 4 families (+0 excluded: members absent from subset) (0 control-track, 1 non-headline); families tree hash 6a93f4150bb2.
Descriptive only (mean over items; a collection is one item; no CIs). Headline = tracks family+collection, headline-eligible; control track (K01 RPROT, K02, E05) reported separately. Trivial baselines are rows of every table; skill = (s-b)/(1-b), b = per-family max(singletons, one big HOG, mean relabel). membership F1 is NA for families without complete cells or decoys. Confirmatory test: scripts/confirmatory.py.

## scope: headline

| method | variant | level_bcubed_f1 (n) | skill | merge_accuracy_node (n) | skill | event_recovery_0.5 (n) | skill | membership_f1 (n) | skill | S/D-bal merge acc (post-hoc) |
|---|---|---|---|---|---|---|---|---|---|---|
| positive_control | reference | 1.000 (3) | 1.000 | 1.000 (3) | 1.000 | 1.000 (2) | 1.000 | 1.000 (1) | 1.000 | 1.000 |
| all_singletons |  | 0.941 (3) | -0.260 | 0.550 (3) | -0.250 | 0.000 (2) | 0.000 | 0.000 (1) | -1.224 | 0.250 |
| one_big_hog |  | 0.148 (3) | -12.665 | 0.233 (3) | -0.167 | 0.000 (2) | 0.000 | 0.550 (1) | 0.000 | 0.500 |
| random_relabel | mean10 | 0.916 (3) | -0.339 | 0.616 (3) | -0.035 | 0.000 (2) | 0.000 | 0.277 (1) | -0.607 | 0.360 |
| fastoma | default_declared | 0.997 (3) | 0.918 | 0.867 (3) | 0.667 | 0.667 (2) | 0.667 | 0.049 (1) | -1.116 | 0.688 |
| fastoma | default_late | 0.997 (3) | 0.918 | 0.867 (3) | 0.667 | 0.667 (2) | 0.667 | 0.049 (1) | -1.116 | 0.688 |
| of_flat | default | 0.481 (3) | -7.988 | 0.579 (3) | -0.039 | 0.167 (2) | 0.167 | 0.493 (1) | -0.127 | 0.550 |
| of_flat | dendroblast | 0.489 (3) | -7.733 | 0.578 (3) | -0.041 | 0.167 (2) | 0.167 | 0.493 (1) | -0.127 | 0.549 |
| orthofinder | default | 0.994 (3) | 0.920 | 0.907 (3) | 0.767 | 0.667 (2) | 0.667 | 0.493 (1) | -0.127 | 0.842 |
| orthofinder | dendroblast | 1.000 (3) | 1.000 | 1.000 (3) | 1.000 | 1.000 (2) | 1.000 | 0.493 (1) | -0.127 | 1.000 |

random relabel, min-max of item-mean over 10 seeds: level_bcubed_f1 0.909316-0.919967; merge_accuracy_node 0.608804-0.631333; event_recovery_0.5 0-0; membership_f1 0.260647-0.295063

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
| of_flat | dendroblast | NA (0) | NA | NA (0) | NA | NA (0) | NA | NA (0) | NA | NA |
| orthofinder | default | NA (0) | NA | NA (0) | NA | NA (0) | NA | NA (0) | NA | NA |
| orthofinder | dendroblast | NA (0) | NA | NA (0) | NA | NA (0) | NA | NA (0) | NA | NA |

random relabel, min-max of item-mean over 10 seeds: 

## Sanity gate (tool identical to all-singletons on >50% of units)

- fastoma default_declared: 0.0% of units identical to all-singletons
- fastoma default_late: 0.0% of units identical to all-singletons
- of_flat default: 0.0% of units identical to all-singletons
- of_flat dendroblast: 0.0% of units identical to all-singletons
- orthofinder default: 0.0% of units identical to all-singletons
- orthofinder dendroblast: 0.0% of units identical to all-singletons

## Positive control / degenerate baselines (per-family min/mean/max)

- positive_control level_bcubed_f1: mean 1.0000, min 1.0000, max 1.0000, n=4
- positive_control merge_accuracy_node: mean 1.0000, min 1.0000, max 1.0000, n=4
- positive_control event_recovery_0.5: mean 1.0000, min 1.0000, max 1.0000, n=3
- positive_control membership_f1: mean 1.0000, min 1.0000, max 1.0000, n=2
- all_singletons level_bcubed_f1: mean 0.9491, min 0.8720, max 1.0000, n=4
- all_singletons merge_accuracy_node: mean 0.5554, min 0.2500, max 1.0000, n=4
- all_singletons event_recovery_0.5: mean 0.0000, min 0.0000, max 0.0000, n=3
- all_singletons membership_f1: mean 0.0000, min 0.0000, max 0.0000, n=2
- one_big_hog level_bcubed_f1: mean 0.1398, min 0.0452, max 0.3196, n=4
- one_big_hog merge_accuracy_node: mean 0.2464, min 0.0000, max 0.5000, n=4
- one_big_hog event_recovery_0.5: mean 0.0000, min 0.0000, max 0.0000, n=3
- one_big_hog membership_f1: mean 0.4704, min 0.3906, max 0.5503, n=2
- random_relabel level_bcubed_f1: mean 0.9025, min 0.8611, max 0.9199, n=4
- random_relabel merge_accuracy_node: mean 0.5957, min 0.3980, max 0.9851, n=4
- random_relabel event_recovery_0.5: mean 0.0000, min 0.0000, max 0.0000, n=3
- random_relabel membership_f1: mean 0.2499, min 0.2226, max 0.2772, n=2

POSITIVE CONTROL: PERFECT on all primary endpoints.
