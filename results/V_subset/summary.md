# Panel V / subset: evaluation summary

hobench git 3b738caa34 (dirty=True); 337 families (+165 excluded: members absent from subset) (0 control-track, 2 non-headline); families tree hash 6a93f4150bb2.
Descriptive only (mean over items; a collection is one item; no CIs). Headline = tracks family+collection, headline-eligible; control track (K01 RPROT, K02, E05) reported separately. Trivial baselines are rows of every table; skill = (s-b)/(1-b), b = per-family max(singletons, one big HOG, mean relabel). membership F1 is NA for families without complete cells or decoys. Confirmatory test: scripts/confirmatory.py.

## scope: headline

| method | variant | level_bcubed_f1 (n) | skill | merge_accuracy_node (n) | skill | event_recovery_0.5 (n) | skill | membership_f1 (n) | skill | S/D-bal merge acc (post-hoc) |
|---|---|---|---|---|---|---|---|---|---|---|
| positive_control | reference | 1.000 (3) | 1.000 | 1.000 (3) | 1.000 | 1.000 (2) | 1.000 | 1.000 (1) | 1.000 | 1.000 |
| all_singletons |  | 0.517 (3) | -0.007 | 0.016 (3) | -7.685 | 0.000 (2) | -0.002 | 0.000 (1) | -0.453 | 0.050 |
| one_big_hog |  | 0.479 (3) | -0.643 | 0.924 (3) | 0.000 | 0.000 (2) | -0.002 | 0.312 (1) | 0.000 | 0.801 |
| random_relabel | mean10 | 0.512 (3) | -0.018 | 0.024 (3) | -7.577 | 0.002 (2) | 0.000 | 0.115 (1) | -0.286 | 0.056 |
| fastoma | default_declared | 0.881 (3) | 0.816 | 0.707 (3) | -0.905 | 0.204 (2) | 0.202 | 0.511 (1) | 0.290 | 0.629 |
| fastoma | default_late | 0.882 (3) | 0.818 | 0.709 (3) | -0.886 | 0.222 (2) | 0.221 | 0.511 (1) | 0.290 | 0.632 |
| of_flat | default | 0.664 (3) | -0.115 | 0.903 (3) | -0.271 | 0.056 (2) | 0.054 | 0.721 (1) | 0.595 | 0.792 |
| of_flat | dendroblast | 0.664 (3) | -0.115 | 0.902 (3) | -0.271 | 0.056 (2) | 0.054 | 0.721 (1) | 0.595 | 0.791 |
| orthofinder | default | 0.987 (3) | 0.962 | 0.960 (3) | 0.514 | 0.907 (2) | 0.907 | 0.721 (1) | 0.595 | 0.963 |
| orthofinder | dendroblast | 0.969 (3) | 0.898 | 0.928 (3) | 0.069 | 0.759 (2) | 0.759 | 0.721 (1) | 0.595 | 0.931 |

random relabel, min-max of item-mean over 10 seeds: level_bcubed_f1 0.508436-0.518737; merge_accuracy_node 0.0175803-0.0328821; event_recovery_0.5 0-0.0185185; membership_f1 0.104784-0.127032

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

- fastoma default_declared: 4.7% of units identical to all-singletons
- fastoma default_late: 4.7% of units identical to all-singletons
- of_flat default: 0.0% of units identical to all-singletons
- of_flat dendroblast: 0.3% of units identical to all-singletons
- orthofinder default: 0.3% of units identical to all-singletons
- orthofinder dendroblast: 0.3% of units identical to all-singletons

## Positive control / degenerate baselines (per-family min/mean/max)

- positive_control level_bcubed_f1: mean 1.0000, min 1.0000, max 1.0000, n=337
- positive_control merge_accuracy_node: mean 1.0000, min 1.0000, max 1.0000, n=337
- positive_control event_recovery_0.5: mean 1.0000, min 1.0000, max 1.0000, n=4
- positive_control membership_f1: mean 1.0000, min 1.0000, max 1.0000, n=2
- all_singletons level_bcubed_f1: mean 0.5008, min 0.4611, max 0.7166, n=337
- all_singletons merge_accuracy_node: mean 0.0006, min 0.0000, max 0.1364, n=337
- all_singletons event_recovery_0.5: mean 0.0000, min 0.0000, max 0.0000, n=4
- all_singletons membership_f1: mean 0.0000, min 0.0000, max 0.0000, n=2
- one_big_hog level_bcubed_f1: mean 0.9909, min 0.0722, max 1.0000, n=337
- one_big_hog merge_accuracy_node: mean 0.9979, min 0.6818, max 1.0000, n=337
- one_big_hog event_recovery_0.5: mean 0.0000, min 0.0000, max 0.0000, n=4
- one_big_hog membership_f1: mean 0.2196, min 0.1275, max 0.3117, n=2
- random_relabel level_bcubed_f1: mean 0.5020, min 0.4416, max 0.7198, n=337
- random_relabel merge_accuracy_node: mean 0.0026, min 0.0000, max 0.1589, n=337
- random_relabel event_recovery_0.5: mean 0.0009, min 0.0000, max 0.0037, n=4
- random_relabel membership_f1: mean 0.1257, min 0.1149, max 0.1366, n=2

POSITIVE CONTROL: PERFECT on all primary endpoints.
