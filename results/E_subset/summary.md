# Panel E / subset: evaluation summary

hobench git 5e46b018e2 (dirty=True); 9 families (+0 excluded: members absent from subset) (1 control-track, 3 non-headline); families tree hash 8c1d8bb8a7fc.
Descriptive only (mean over items; a collection is one item; no CIs). Headline = tracks family+collection, headline-eligible; control track (K01 RPROT, K02, E05) reported separately. Trivial baselines are rows of every table; skill = (s-b)/(1-b), b = per-family max(singletons, one big HOG, mean relabel). membership F1 is NA for families without complete cells or decoys. Confirmatory test: scripts/confirmatory.py.

## scope: headline

| method | variant | level_bcubed_f1 (n) | skill | merge_accuracy_node (n) | skill | event_recovery_0.5 (n) | skill | membership_f1 (n) | skill | S/D-bal merge acc (post-hoc) |
|---|---|---|---|---|---|---|---|---|---|---|
| positive_control | reference | 1.000 (5) | 1.000 | 1.000 (5) | 1.000 | 1.000 (5) | 1.000 | 1.000 (3) | 1.000 | 1.000 |
| all_singletons |  | 0.849 (5) | 0.000 | 0.146 (5) | -5.600 | 0.000 (5) | 0.000 | 0.075 (3) | -0.059 | 0.363 |
| one_big_hog |  | 0.217 (5) | -10.447 | 0.838 (5) | 0.000 | 0.000 (5) | 0.000 | 0.049 (3) | -0.095 | 0.603 |
| random_relabel | mean10 | 0.823 (5) | -0.409 | 0.149 (5) | -5.574 | 0.000 (5) | 0.000 | 0.120 (3) | -0.010 | 0.360 |
| fastoma | default_declared | 0.906 (5) | 0.371 | 0.531 (5) | -3.153 | 0.025 (5) | 0.025 | 0.574 (3) | 0.499 | 0.628 |
| fastoma | default_late | 0.904 (5) | 0.342 | 0.531 (5) | -3.153 | 0.025 (5) | 0.025 | 0.574 (3) | 0.499 | 0.628 |
| of_flat | default | 0.358 (5) | -9.673 | 0.758 (5) | -0.716 | 0.000 (5) | 0.000 | 0.289 (3) | 0.180 | 0.626 |
| orthofinder | default | 0.981 (5) | 0.871 | 0.729 (5) | -1.235 | 0.400 (5) | 0.400 | 0.741 (3) | 0.701 | 0.805 |

random relabel, min-max of item-mean over 10 seeds: level_bcubed_f1 0.819188-0.826205; merge_accuracy_node 0.14592-0.154424; event_recovery_0.5 0-0; membership_f1 0.116867-0.124755

## scope: control

| method | variant | level_bcubed_f1 (n) | skill | merge_accuracy_node (n) | skill | event_recovery_0.5 (n) | skill | membership_f1 (n) | skill | S/D-bal merge acc (post-hoc) |
|---|---|---|---|---|---|---|---|---|---|---|
| positive_control | reference | 1.000 (1) | 1.000 | 1.000 (1) | 1.000 | 1.000 (1) | 1.000 | 1.000 (1) | 1.000 | 1.000 |
| all_singletons |  | 0.486 (1) | -0.058 | 0.000 (1) | -10.000 | 0.000 (1) | 0.000 | 0.000 (1) | -0.303 | 0.000 |
| one_big_hog |  | 0.358 (1) | -0.321 | 0.909 (1) | 0.000 | 0.000 (1) | 0.000 | 0.233 (1) | 0.000 | 0.750 |
| random_relabel | mean10 | 0.514 (1) | 0.000 | 0.030 (1) | -9.671 | 0.000 (1) | 0.000 | 0.147 (1) | -0.111 | 0.043 |
| fastoma | default_declared | 0.901 (1) | 0.797 | 0.752 (1) | -1.733 | 0.000 (1) | 0.000 | 0.621 (1) | 0.506 | 0.654 |
| fastoma | default_late | 0.901 (1) | 0.797 | 0.752 (1) | -1.733 | 0.000 (1) | 0.000 | 0.621 (1) | 0.506 | 0.654 |
| of_flat | default | 0.358 (1) | -0.321 | 0.909 (1) | 0.000 | 0.000 (1) | 0.000 | 0.983 (1) | 0.977 | 0.750 |
| orthofinder | default | 1.000 (1) | 1.000 | 1.000 (1) | 1.000 | 1.000 (1) | 1.000 | 0.983 (1) | 0.977 | 1.000 |

random relabel, min-max of item-mean over 10 seeds: level_bcubed_f1 0.510248-0.520778; merge_accuracy_node 0.0255337-0.0348721; event_recovery_0.5 0-0; membership_f1 0.139873-0.155577

## Sanity gate (tool identical to all-singletons on >50% of units)

- fastoma default_declared: 0.0% of units identical to all-singletons
- fastoma default_late: 0.0% of units identical to all-singletons
- of_flat default: 0.0% of units identical to all-singletons
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
