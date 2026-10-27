# Panel E / subset: evaluation summary

hobench git 0b5dd0a216 (dirty=False); 9 families (+0 excluded: members absent from subset) (1 control-track, 4 non-headline); families tree hash 8c1d8bb8a7fc.
Descriptive only (mean over items; a collection is one item; no CIs). Headline = tracks family+collection, headline-eligible; control track (K01 RPROT, K02, E05) reported separately. Trivial baselines are rows of every table; skill = (s-b)/(1-b), b = per-family max(singletons, one big HOG, mean relabel). membership F1 is NA for families without complete cells or decoys. Confirmatory test: scripts/confirmatory.py.

## scope: headline

| method | variant | level_bcubed_f1 (n) | skill | merge_accuracy_node (n) | skill | event_recovery_0.5 (n) | skill | membership_f1 (n) | skill | S/D-bal merge acc (post-hoc) |
|---|---|---|---|---|---|---|---|---|---|---|
| positive_control | reference | 1.000 (2) | 1.000 | 1.000 (2) | 1.000 | 1.000 (2) | 1.000 | 1.000 (2) | 1.000 | 1.000 |
| all_singletons |  | 0.863 (2) | 0.000 | 0.132 (2) | -6.250 | 0.000 (2) | 0.000 | 0.101 (2) | -0.044 | 0.329 |
| one_big_hog |  | 0.137 (2) | -12.755 | 0.848 (2) | 0.000 | 0.000 (2) | 0.000 | 0.042 (2) | -0.119 | 0.629 |
| random_relabel | mean10 | 0.832 (2) | -0.510 | 0.136 (2) | -6.220 | 0.000 (2) | 0.000 | 0.129 (2) | -0.014 | 0.327 |
| fastoma | default_declared | 0.883 (2) | 0.220 | 0.484 (2) | -3.841 | 0.031 (2) | 0.031 | 0.534 (2) | 0.447 | 0.579 |
| fastoma | default_late | 0.881 (2) | 0.183 | 0.484 (2) | -3.841 | 0.031 (2) | 0.031 | 0.535 (2) | 0.447 | 0.579 |
| orthofinder | default | 0.978 (2) | 0.847 | 0.747 (2) | -1.369 | 0.500 (2) | 0.500 | 0.748 (2) | 0.705 | 0.810 |

random relabel, min-max of item-mean over 10 seeds: level_bcubed_f1 0.826558-0.835374; merge_accuracy_node 0.132884-0.142501; event_recovery_0.5 0-0; membership_f1 0.123498-0.135269

## scope: control

| method | variant | level_bcubed_f1 (n) | skill | merge_accuracy_node (n) | skill | event_recovery_0.5 (n) | skill | membership_f1 (n) | skill | S/D-bal merge acc (post-hoc) |
|---|---|---|---|---|---|---|---|---|---|---|
| positive_control | reference | 1.000 (1) | 1.000 | 1.000 (1) | 1.000 | 1.000 (1) | 1.000 | 1.000 (1) | 1.000 | 1.000 |
| all_singletons |  | 0.486 (1) | -0.058 | 0.000 (1) | -10.000 | 0.000 (1) | 0.000 | 0.000 (1) | -0.303 | 0.000 |
| one_big_hog |  | 0.358 (1) | -0.321 | 0.909 (1) | 0.000 | 0.000 (1) | 0.000 | 0.233 (1) | 0.000 | 0.750 |
| random_relabel | mean10 | 0.514 (1) | 0.000 | 0.030 (1) | -9.671 | 0.000 (1) | 0.000 | 0.147 (1) | -0.111 | 0.043 |
| fastoma | default_declared | 0.901 (1) | 0.797 | 0.752 (1) | -1.733 | 0.000 (1) | 0.000 | 0.621 (1) | 0.506 | 0.654 |
| fastoma | default_late | 0.901 (1) | 0.797 | 0.752 (1) | -1.733 | 0.000 (1) | 0.000 | 0.621 (1) | 0.506 | 0.654 |
| orthofinder | default | 1.000 (1) | 1.000 | 1.000 (1) | 1.000 | 1.000 (1) | 1.000 | 0.983 (1) | 0.977 | 1.000 |

random relabel, min-max of item-mean over 10 seeds: level_bcubed_f1 0.510248-0.520778; merge_accuracy_node 0.0255337-0.0348721; event_recovery_0.5 0-0; membership_f1 0.139873-0.155577

## Sanity gate (tool identical to all-singletons on >50% of units)

- fastoma default_declared: 0.0% of units identical to all-singletons
- fastoma default_late: 0.0% of units identical to all-singletons
- orthofinder default: 0.0% of units identical to all-singletons

## Positive control / degenerate baselines (per-family min/mean/max)

- positive_control level_bcubed_f1: mean 1.0000, min 1.0000, max 1.0000, n=9
- positive_control merge_accuracy_node: mean 1.0000, min 1.0000, max 1.0000, n=9
- positive_control event_recovery_0.5: mean 1.0000, min 1.0000, max 1.0000, n=9
- positive_control membership_f1: mean 1.0000, min 1.0000, max 1.0000, n=5
- all_singletons level_bcubed_f1: mean 0.8189, min 0.4861, max 0.9826, n=9
- all_singletons merge_accuracy_node: mean 0.1154, min 0.0000, max 0.2500, n=9
- all_singletons event_recovery_0.5: mean 0.0000, min 0.0000, max 0.0000, n=9
- all_singletons membership_f1: mean 0.0579, min 0.0000, max 0.1785, n=5
- one_big_hog level_bcubed_f1: mean 0.2189, min 0.0745, max 0.5352, n=9
- one_big_hog merge_accuracy_node: mean 0.8655, min 0.7500, max 0.9333, n=9
- one_big_hog event_recovery_0.5: mean 0.0000, min 0.0000, max 0.0000, n=9
- one_big_hog membership_f1: mean 0.1111, min 0.0214, max 0.2325, n=5
- random_relabel level_bcubed_f1: mean 0.8032, min 0.5142, max 0.9590, n=9
- random_relabel merge_accuracy_node: mean 0.1242, min 0.0299, max 0.2538, n=9
- random_relabel event_recovery_0.5: mean 0.0000, min 0.0000, max 0.0000, n=9
- random_relabel membership_f1: mean 0.1290, min 0.0803, max 0.1549, n=5

POSITIVE CONTROL: PERFECT on all primary endpoints.
