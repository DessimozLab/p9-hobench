# Panel Y / full: evaluation summary

hobench git 99ae136671 (dirty=True); 1977 families (+0 excluded: members absent from subset) (0 control-track, 0 non-headline); families tree hash 6a93f4150bb2.
Descriptive only (mean over items; a collection is one item; no CIs). Headline = tracks family+collection, headline-eligible; control track (K01 RPROT, K02, E05) reported separately. Trivial baselines are rows of every table; skill = (s-b)/(1-b), b = per-family max(singletons, one big HOG, mean relabel). membership F1 is NA for families without complete cells or decoys. Confirmatory test: scripts/confirmatory.py.

## scope: headline

| method | variant | level_bcubed_f1 (n) | skill | merge_accuracy_node (n) | skill | event_recovery_0.5 (n) | skill | membership_f1 (n) | skill | S/D-bal merge acc (post-hoc) |
|---|---|---|---|---|---|---|---|---|---|---|
| positive_control | reference | 1.000 (2) | 1.000 | 1.000 (2) | 1.000 | 1.000 (2) | 1.000 | NA (0) | NA | 1.000 |
| all_singletons |  | 0.695 (2) | -0.808 | 0.128 (2) | -2.947 | 0.000 (2) | 0.000 | NA (0) | NA | 0.500 |
| one_big_hog |  | 0.790 (2) | -0.398 | 0.872 (2) | -0.084 | 0.000 (2) | 0.000 | NA (0) | NA | 0.500 |
| random_relabel | mean10 | 0.694 (2) | -0.808 | 0.129 (2) | -2.947 | 0.000 (2) | 0.000 | NA (0) | NA | 0.500 |
| fastoma | default_declared | 0.920 (2) | 0.470 | 0.942 (2) | 0.633 | 0.134 (2) | 0.134 | NA (0) | NA | 0.756 |
| fastoma | default_late | 0.920 (2) | 0.468 | 0.942 (2) | 0.628 | 0.139 (2) | 0.139 | NA (0) | NA | 0.754 |
| fastoma | rep2_declared | 0.920 (2) | 0.469 | 0.942 (2) | 0.632 | 0.134 (2) | 0.134 | NA (0) | NA | 0.755 |
| fastoma | rep2_late | 0.920 (2) | 0.467 | 0.941 (2) | 0.627 | 0.139 (2) | 0.139 | NA (0) | NA | 0.753 |
| of_flat | default | 0.790 (2) | -0.360 | 0.872 (2) | -0.049 | 0.007 (2) | 0.007 | NA (0) | NA | 0.539 |
| of_flat | dendroblast | 0.790 (2) | -0.358 | 0.872 (2) | -0.047 | 0.010 (2) | 0.010 | NA (0) | NA | 0.541 |
| orthofinder | default | 0.928 (2) | 0.425 | 0.879 (2) | 0.031 | 0.060 (2) | 0.060 | NA (0) | NA | 0.596 |
| orthofinder | dendroblast | 0.928 (2) | 0.425 | 0.880 (2) | 0.047 | 0.072 (2) | 0.072 | NA (0) | NA | 0.603 |
| orthofinder | rep2 | 0.928 (2) | 0.425 | 0.879 (2) | 0.031 | 0.060 (2) | 0.060 | NA (0) | NA | 0.596 |

random relabel, min-max of item-mean over 10 seeds: level_bcubed_f1 0.692366-0.694716; merge_accuracy_node 0.127503-0.128844; event_recovery_0.5 0-0

## scope: control

| method | variant | level_bcubed_f1 (n) | skill | merge_accuracy_node (n) | skill | event_recovery_0.5 (n) | skill | membership_f1 (n) | skill | S/D-bal merge acc (post-hoc) |
|---|---|---|---|---|---|---|---|---|---|---|
| positive_control | reference | NA (0) | NA | NA (0) | NA | NA (0) | NA | NA (0) | NA | NA |
| all_singletons |  | NA (0) | NA | NA (0) | NA | NA (0) | NA | NA (0) | NA | NA |
| one_big_hog |  | NA (0) | NA | NA (0) | NA | NA (0) | NA | NA (0) | NA | NA |
| random_relabel | mean10 | NA (0) | NA | NA (0) | NA | NA (0) | NA | NA (0) | NA | NA |
| fastoma | default_declared | NA (0) | NA | NA (0) | NA | NA (0) | NA | NA (0) | NA | NA |
| fastoma | default_late | NA (0) | NA | NA (0) | NA | NA (0) | NA | NA (0) | NA | NA |
| fastoma | rep2_declared | NA (0) | NA | NA (0) | NA | NA (0) | NA | NA (0) | NA | NA |
| fastoma | rep2_late | NA (0) | NA | NA (0) | NA | NA (0) | NA | NA (0) | NA | NA |
| of_flat | default | NA (0) | NA | NA (0) | NA | NA (0) | NA | NA (0) | NA | NA |
| of_flat | dendroblast | NA (0) | NA | NA (0) | NA | NA (0) | NA | NA (0) | NA | NA |
| orthofinder | default | NA (0) | NA | NA (0) | NA | NA (0) | NA | NA (0) | NA | NA |
| orthofinder | dendroblast | NA (0) | NA | NA (0) | NA | NA (0) | NA | NA (0) | NA | NA |
| orthofinder | rep2 | NA (0) | NA | NA (0) | NA | NA (0) | NA | NA (0) | NA | NA |

random relabel, min-max of item-mean over 10 seeds: 

## Sanity gate (tool identical to all-singletons on >50% of units)

- fastoma default_declared: 1.9% of units identical to all-singletons
- fastoma default_late: 1.9% of units identical to all-singletons
- fastoma rep2_declared: 1.8% of units identical to all-singletons
- fastoma rep2_late: 1.8% of units identical to all-singletons
- of_flat default: 0.4% of units identical to all-singletons
- of_flat dendroblast: 0.4% of units identical to all-singletons
- orthofinder default: 0.5% of units identical to all-singletons
- orthofinder dendroblast: 0.6% of units identical to all-singletons
- orthofinder rep2: 0.5% of units identical to all-singletons

## Positive control / degenerate baselines (per-family min/mean/max)

- positive_control level_bcubed_f1: mean 1.0000, min 1.0000, max 1.0000, n=1977
- positive_control merge_accuracy_node: mean 1.0000, min 1.0000, max 1.0000, n=1977
- positive_control event_recovery_0.5: mean 1.0000, min 1.0000, max 1.0000, n=202
- all_singletons level_bcubed_f1: mean 0.6052, min 0.4606, max 1.0000, n=1977
- all_singletons merge_accuracy_node: mean 0.1141, min 0.0000, max 1.0000, n=1977
- all_singletons event_recovery_0.5: mean 0.0000, min 0.0000, max 0.0000, n=202
- one_big_hog level_bcubed_f1: mean 0.9638, min 0.6152, max 1.0000, n=1977
- one_big_hog merge_accuracy_node: mean 0.8859, min 0.0000, max 1.0000, n=1977
- one_big_hog event_recovery_0.5: mean 0.0000, min 0.0000, max 0.0000, n=202
- random_relabel level_bcubed_f1: mean 0.6054, min 0.4606, max 1.0000, n=1977
- random_relabel merge_accuracy_node: mean 0.1145, min 0.0000, max 1.0000, n=1977
- random_relabel event_recovery_0.5: mean 0.0000, min 0.0000, max 0.0000, n=202

POSITIVE CONTROL: PERFECT on all primary endpoints.
