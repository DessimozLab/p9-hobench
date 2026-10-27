# Panel K / full: evaluation summary

hobench git c304c5e4b8 (dirty=True); 51 families (+0 excluded: members absent from subset) (41 control-track, 4 non-headline); families tree hash 8c1d8bb8a7fc.
Descriptive only (mean over items; a collection is one item; no CIs). Headline = tracks family+collection, headline-eligible; control track (K01 RPROT, K02, E05) reported separately. Trivial baselines are rows of every table; skill = (s-b)/(1-b), b = per-family max(singletons, one big HOG, mean relabel). membership F1 is NA for families without complete cells or decoys. Confirmatory test: scripts/confirmatory.py.

## scope: headline

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

## scope: control

| method | variant | level_bcubed_f1 (n) | skill | merge_accuracy_node (n) | skill | event_recovery_0.5 (n) | skill | membership_f1 (n) | skill | S/D-bal merge acc (post-hoc) |
|---|---|---|---|---|---|---|---|---|---|---|
| positive_control | reference | 1.000 (2) | 1.000 | 1.000 (2) | 1.000 | 1.000 (1) | 1.000 | 1.000 (1) | 1.000 | 1.000 |
| all_singletons |  | 0.611 (2) | 0.000 | 0.167 (2) | -1.000 | 0.000 (1) | 0.000 | 0.000 (1) | -0.100 | 0.500 |
| one_big_hog |  | 0.764 (2) | -0.675 | 0.833 (2) | 0.000 | 0.000 (1) | 0.000 | 0.047 (1) | -0.048 | 0.500 |
| random_relabel | mean10 | 0.617 (2) | -0.001 | 0.181 (2) | -0.967 | 0.000 (1) | 0.000 | 0.090 (1) | 0.000 | 0.505 |
| fastoma | default_declared | 0.969 (2) | 0.933 | 0.676 (2) | -0.034 | 0.000 (1) | 0.000 | 0.818 (1) | 0.802 | 0.742 |
| fastoma | default_late | 0.969 (2) | 0.933 | 0.676 (2) | -0.034 | 0.000 (1) | 0.000 | 0.818 (1) | 0.802 | 0.742 |
| of_flat | default | 0.991 (2) | 0.949 | 0.893 (2) | 0.431 | 0.000 (1) | 0.000 | 0.985 (1) | 0.984 | 0.858 |
| of_flat | dendroblast | 0.996 (2) | 0.982 | 0.904 (2) | 0.497 | 0.000 (1) | 0.000 | 0.985 (1) | 0.984 | 0.874 |
| orthofinder | default | 0.980 (2) | 0.873 | 0.883 (2) | 0.371 | 0.000 (1) | 0.000 | 0.985 (1) | 0.984 | 0.843 |
| orthofinder | dendroblast | 0.985 (2) | 0.907 | 0.894 (2) | 0.436 | 0.000 (1) | 0.000 | 0.985 (1) | 0.984 | 0.859 |

random relabel, min-max of item-mean over 10 seeds: level_bcubed_f1 0.612697-0.624494; merge_accuracy_node 0.173102-0.187408; event_recovery_0.5 0-0; membership_f1 0.0881279-0.0947852

## Sanity gate (tool identical to all-singletons on >50% of units)

- fastoma default_declared: 19.6% of units identical to all-singletons
- fastoma default_late: 19.6% of units identical to all-singletons
- of_flat default: 19.6% of units identical to all-singletons
- of_flat dendroblast: 19.6% of units identical to all-singletons
- orthofinder default: 19.6% of units identical to all-singletons
- orthofinder dendroblast: 19.6% of units identical to all-singletons

## Positive control / degenerate baselines (per-family min/mean/max)

- positive_control level_bcubed_f1: mean 1.0000, min 1.0000, max 1.0000, n=41
- positive_control merge_accuracy_node: mean 1.0000, min 1.0000, max 1.0000, n=41
- positive_control event_recovery_0.5: mean 1.0000, min 1.0000, max 1.0000, n=1
- positive_control membership_f1: mean 1.0000, min 1.0000, max 1.0000, n=40
- all_singletons level_bcubed_f1: mean 0.5103, min 0.4838, max 0.7176, n=41
- all_singletons merge_accuracy_node: mean 0.0081, min 0.0000, max 0.3333, n=41
- all_singletons event_recovery_0.5: mean 0.0000, min 0.0000, max 0.0000, n=1
- all_singletons membership_f1: mean 0.0000, min 0.0000, max 0.0000, n=40
- one_big_hog level_bcubed_f1: mean 0.9885, min 0.5271, max 1.0000, n=41
- one_big_hog merge_accuracy_node: mean 0.9919, min 0.6667, max 1.0000, n=41
- one_big_hog event_recovery_0.5: mean 0.0000, min 0.0000, max 0.0000, n=1
- one_big_hog membership_f1: mean 0.0472, min 0.0442, max 0.0493, n=40
- random_relabel level_bcubed_f1: mean 0.5223, min 0.4971, max 0.7174, n=41
- random_relabel merge_accuracy_node: mean 0.0262, min 0.0106, max 0.3443, n=41
- random_relabel event_recovery_0.5: mean 0.0000, min 0.0000, max 0.0000, n=1
- random_relabel membership_f1: mean 0.0904, min 0.0714, max 0.1277, n=40

POSITIVE CONTROL: PERFECT on all primary endpoints.
