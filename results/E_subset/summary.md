# Panel E / subset: evaluation summary

hobench git 3b738caa34 (dirty=True); 9 families (+0 excluded: members absent from subset) (1 control-track, 3 non-headline); families tree hash 6a93f4150bb2.
Descriptive only (mean over items; a collection is one item; no CIs). Headline = tracks family+collection, headline-eligible; control track (K01 RPROT, K02, E05) reported separately. Trivial baselines are rows of every table; skill = (s-b)/(1-b), b = per-family max(singletons, one big HOG, mean relabel). membership F1 is NA for families without complete cells or decoys. Confirmatory test: scripts/confirmatory.py.

## scope: headline

| method | variant | level_bcubed_f1 (n) | skill | merge_accuracy_node (n) | skill | event_recovery_0.5 (n) | skill | membership_f1 (n) | skill | S/D-bal merge acc (post-hoc) |
|---|---|---|---|---|---|---|---|---|---|---|
| positive_control | reference | 1.000 (5) | 1.000 | 1.000 (5) | 1.000 | 1.000 (5) | 1.000 | 1.000 (3) | 1.000 | 1.000 |
| all_singletons |  | 0.849 (5) | 0.000 | 0.146 (5) | -5.600 | 0.000 (5) | 0.000 | 0.000 (3) | -0.223 | 0.363 |
| one_big_hog |  | 0.217 (5) | -10.447 | 0.838 (5) | 0.000 | 0.000 (5) | 0.000 | 0.014 (3) | -0.206 | 0.603 |
| random_relabel | mean10 | 0.835 (5) | -0.210 | 0.147 (5) | -5.590 | 0.000 (5) | 0.000 | 0.180 (3) | 0.000 | 0.363 |
| fastoma | default_declared | 0.906 (5) | 0.371 | 0.531 (5) | -3.153 | 0.025 (5) | 0.025 | 0.563 (3) | 0.466 | 0.628 |
| fastoma | default_late | 0.904 (5) | 0.342 | 0.531 (5) | -3.153 | 0.025 (5) | 0.025 | 0.563 (3) | 0.466 | 0.628 |
| of_flat | default | 0.358 (5) | -9.673 | 0.758 (5) | -0.716 | 0.000 (5) | 0.000 | 0.152 (3) | -0.040 | 0.626 |
| of_flat | dendroblast | 0.353 (5) | -9.709 | 0.757 (5) | -0.719 | 0.000 (5) | 0.000 | 0.152 (3) | -0.040 | 0.625 |
| orthofinder | default | 0.981 (5) | 0.871 | 0.729 (5) | -1.235 | 0.400 (5) | 0.400 | 0.870 (3) | 0.840 | 0.805 |
| orthofinder | dendroblast | 0.983 (5) | 0.905 | 0.795 (5) | -0.715 | 0.625 (5) | 0.625 | 0.907 (3) | 0.886 | 0.858 |

random relabel, min-max of item-mean over 10 seeds: level_bcubed_f1 0.833067-0.836406; merge_accuracy_node 0.145601-0.148992; event_recovery_0.5 0-0; membership_f1 0.172638-0.189956

## scope: control

| method | variant | level_bcubed_f1 (n) | skill | merge_accuracy_node (n) | skill | event_recovery_0.5 (n) | skill | membership_f1 (n) | skill | S/D-bal merge acc (post-hoc) |
|---|---|---|---|---|---|---|---|---|---|---|
| positive_control | reference | 1.000 (1) | 1.000 | 1.000 (1) | 1.000 | 1.000 (1) | 1.000 | 1.000 (1) | 1.000 | 1.000 |
| all_singletons |  | 0.486 (1) | -0.043 | 0.000 (1) | -10.000 | 0.000 (1) | 0.000 | 0.000 (1) | -0.303 | 0.000 |
| one_big_hog |  | 0.358 (1) | -0.303 | 0.909 (1) | 0.000 | 0.000 (1) | 0.000 | 0.233 (1) | 0.000 | 0.750 |
| random_relabel | mean10 | 0.507 (1) | 0.000 | 0.023 (1) | -9.746 | 0.000 (1) | 0.000 | 0.136 (1) | -0.125 | 0.034 |
| fastoma | default_declared | 0.901 (1) | 0.800 | 0.752 (1) | -1.733 | 0.000 (1) | 0.000 | 0.621 (1) | 0.506 | 0.654 |
| fastoma | default_late | 0.901 (1) | 0.800 | 0.752 (1) | -1.733 | 0.000 (1) | 0.000 | 0.621 (1) | 0.506 | 0.654 |
| of_flat | default | 0.358 (1) | -0.303 | 0.909 (1) | 0.000 | 0.000 (1) | 0.000 | 0.983 (1) | 0.977 | 0.750 |
| of_flat | dendroblast | 0.358 (1) | -0.303 | 0.909 (1) | 0.000 | 0.000 (1) | 0.000 | 0.983 (1) | 0.977 | 0.750 |
| orthofinder | default | 1.000 (1) | 1.000 | 1.000 (1) | 1.000 | 1.000 (1) | 1.000 | 0.983 (1) | 0.977 | 1.000 |
| orthofinder | dendroblast | 0.997 (1) | 0.995 | 0.898 (1) | -0.125 | 1.000 (1) | 1.000 | 0.976 (1) | 0.969 | 0.938 |

random relabel, min-max of item-mean over 10 seeds: level_bcubed_f1 0.504128-0.515157; merge_accuracy_node 0.0174643-0.0302245; event_recovery_0.5 0-0; membership_f1 0.125464-0.14518

## Sanity gate (tool identical to all-singletons on >50% of units)

- fastoma default_declared: 0.0% of units identical to all-singletons
- fastoma default_late: 0.0% of units identical to all-singletons
- of_flat default: 0.0% of units identical to all-singletons
- of_flat dendroblast: 0.0% of units identical to all-singletons
- orthofinder default: 0.0% of units identical to all-singletons
- orthofinder dendroblast: 0.0% of units identical to all-singletons

## Positive control / degenerate baselines (per-family min/mean/max)

- positive_control level_bcubed_f1: mean 1.0000, min 1.0000, max 1.0000, n=9
- positive_control merge_accuracy_node: mean 1.0000, min 1.0000, max 1.0000, n=9
- positive_control event_recovery_0.5: mean 1.0000, min 1.0000, max 1.0000, n=9
- positive_control membership_f1: mean 1.0000, min 1.0000, max 1.0000, n=5
- all_singletons level_bcubed_f1: mean 0.8189, min 0.4861, max 0.9826, n=9
- all_singletons merge_accuracy_node: mean 0.1154, min 0.0000, max 0.2500, n=9
- all_singletons event_recovery_0.5: mean 0.0000, min 0.0000, max 0.0000, n=9
- all_singletons membership_f1: mean 0.0000, min 0.0000, max 0.0000, n=5
- one_big_hog level_bcubed_f1: mean 0.2189, min 0.0745, max 0.5352, n=9
- one_big_hog merge_accuracy_node: mean 0.8655, min 0.7500, max 0.9333, n=9
- one_big_hog event_recovery_0.5: mean 0.0000, min 0.0000, max 0.0000, n=9
- one_big_hog membership_f1: mean 0.0562, min 0.0074, max 0.2325, n=5
- random_relabel level_bcubed_f1: mean 0.8113, min 0.5074, max 0.9707, n=9
- random_relabel merge_accuracy_node: mean 0.1218, min 0.0231, max 0.2514, n=9
- random_relabel event_recovery_0.5: mean 0.0000, min 0.0000, max 0.0000, n=9
- random_relabel membership_f1: mean 0.1958, min 0.1339, max 0.3021, n=5

POSITIVE CONTROL: PERFECT on all primary endpoints.
