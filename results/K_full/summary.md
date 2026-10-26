# Panel K / full: evaluation summary

hobench git 7149dd29c6 (dirty=True); 51 families (41 control-track, 4 non-headline); families tree hash aeea5fbf3546.
Descriptive only (mean over items; a collection is one item; no CIs). Headline = tracks family+collection, headline-eligible; control track (K01 RPROT, K02, E05) reported separately. Trivial baselines are rows of every table; skill = (s-b)/(1-b), b = per-family max(singletons, one big HOG, mean relabel). membership F1 is NA for families without complete cells or decoys. Confirmatory test: scripts/confirmatory.py.

## scope: headline

| method | variant | level_bcubed_f1 (n) | skill | merge_accuracy_node (n) | skill | event_recovery_0.5 (n) | skill | membership_f1 (n) | skill | S/D-bal merge acc (post-hoc) |
|---|---|---|---|---|---|---|---|---|---|---|
| positive_control | reference | NA (0) | NA | NA (0) | NA | NA (0) | NA | 1.000 (1) | 1.000 | NA |
| all_singletons |  | NA (0) | NA | NA (0) | NA | NA (0) | NA | 0.000 (1) | -4.595 | NA |
| one_big_hog |  | NA (0) | NA | NA (0) | NA | NA (0) | NA | 0.753 (1) | 0.000 | NA |
| random_relabel | mean10 | NA (0) | NA | NA (0) | NA | NA (0) | NA | 0.221 (1) | -3.333 | NA |
| fastoma | default_declared | NA (0) | NA | NA (0) | NA | NA (0) | NA | 0.897 (1) | 0.595 | NA |
| fastoma | default_late | NA (0) | NA | NA (0) | NA | NA (0) | NA | 0.897 (1) | 0.595 | NA |
| of_flat | default | NA (0) | NA | NA (0) | NA | NA (0) | NA | 0.949 (1) | 0.730 | NA |
| of_flat | dendroblast | NA (0) | NA | NA (0) | NA | NA (0) | NA | 0.956 (1) | 0.772 | NA |
| orthofinder | default | NA (0) | NA | NA (0) | NA | NA (0) | NA | 0.949 (1) | 0.730 | NA |
| orthofinder | dendroblast | NA (0) | NA | NA (0) | NA | NA (0) | NA | 0.956 (1) | 0.772 | NA |

random relabel, min-max of item-mean over 10 seeds: membership_f1 0.206552-0.24145

## scope: control

| method | variant | level_bcubed_f1 (n) | skill | merge_accuracy_node (n) | skill | event_recovery_0.5 (n) | skill | membership_f1 (n) | skill | S/D-bal merge acc (post-hoc) |
|---|---|---|---|---|---|---|---|---|---|---|
| positive_control | reference | 1.000 (2) | 1.000 | 1.000 (2) | 1.000 | 1.000 (1) | 1.000 | 1.000 (1) | 1.000 | 1.000 |
| all_singletons |  | 0.611 (2) | 0.000 | 0.167 (2) | -1.000 | 0.000 (1) | 0.000 | 0.000 (1) | -0.101 | 0.500 |
| one_big_hog |  | 0.764 (2) | -0.675 | 0.833 (2) | 0.000 | 0.000 (1) | 0.000 | 0.047 (1) | -0.049 | 0.500 |
| random_relabel | mean10 | 0.618 (2) | -0.001 | 0.182 (2) | -0.966 | 0.000 (1) | 0.000 | 0.091 (1) | 0.000 | 0.504 |
| fastoma | default_declared | 0.969 (2) | 0.933 | 0.676 (2) | -0.034 | 0.000 (1) | 0.000 | 0.818 (1) | 0.802 | 0.742 |
| fastoma | default_late | 0.969 (2) | 0.933 | 0.676 (2) | -0.034 | 0.000 (1) | 0.000 | 0.818 (1) | 0.802 | 0.742 |
| of_flat | default | 0.991 (2) | 0.949 | 0.893 (2) | 0.431 | 0.000 (1) | 0.000 | 0.985 (1) | 0.984 | 0.858 |
| of_flat | dendroblast | 0.996 (2) | 0.982 | 0.904 (2) | 0.497 | 0.000 (1) | 0.000 | 0.985 (1) | 0.984 | 0.874 |
| orthofinder | default | 0.980 (2) | 0.873 | 0.883 (2) | 0.371 | 0.000 (1) | 0.000 | 0.985 (1) | 0.984 | 0.843 |
| orthofinder | dendroblast | 0.985 (2) | 0.907 | 0.894 (2) | 0.436 | 0.000 (1) | 0.000 | 0.985 (1) | 0.984 | 0.859 |

random relabel, min-max of item-mean over 10 seeds: level_bcubed_f1 0.612807-0.624325; merge_accuracy_node 0.172614-0.188326; event_recovery_0.5 0-0; membership_f1 0.0884645-0.0951757

## Sanity gate (tool identical to all-singletons on >50% of units)

- fastoma default_declared: 0.0% of units identical to all-singletons
- fastoma default_late: 0.0% of units identical to all-singletons
- of_flat default: 0.0% of units identical to all-singletons
- of_flat dendroblast: 0.0% of units identical to all-singletons
- orthofinder default: 0.0% of units identical to all-singletons
- orthofinder dendroblast: 0.0% of units identical to all-singletons

## Positive control / degenerate baselines (per-family min/mean/max)

- positive_control level_bcubed_f1: mean 1.0000, min 1.0000, max 1.0000, n=41
- positive_control merge_accuracy_node: mean 1.0000, min 1.0000, max 1.0000, n=41
- positive_control event_recovery_0.5: mean 1.0000, min 1.0000, max 1.0000, n=1
- positive_control membership_f1: mean 1.0000, min 1.0000, max 1.0000, n=50
- all_singletons level_bcubed_f1: mean 0.5103, min 0.4838, max 0.7176, n=41
- all_singletons merge_accuracy_node: mean 0.0081, min 0.0000, max 0.3333, n=41
- all_singletons event_recovery_0.5: mean 0.0000, min 0.0000, max 0.0000, n=1
- all_singletons membership_f1: mean 0.0000, min 0.0000, max 0.0000, n=50
- one_big_hog level_bcubed_f1: mean 0.9885, min 0.5271, max 1.0000, n=41
- one_big_hog merge_accuracy_node: mean 0.9919, min 0.6667, max 1.0000, n=41
- one_big_hog event_recovery_0.5: mean 0.0000, min 0.0000, max 0.0000, n=1
- one_big_hog membership_f1: mean 0.1884, min 0.0442, max 0.9333, n=50
- random_relabel level_bcubed_f1: mean 0.5225, min 0.4971, max 0.7173, n=41
- random_relabel merge_accuracy_node: mean 0.0270, min 0.0106, max 0.3447, n=41
- random_relabel event_recovery_0.5: mean 0.0000, min 0.0000, max 0.0000, n=1
- random_relabel membership_f1: mean 0.1172, min 0.0736, max 0.3212, n=50

POSITIVE CONTROL: PERFECT on all primary endpoints.
