# Panel Y / subset: evaluation summary

hobench git d2c56a6326 (dirty=True); 1976 families (+1 excluded: members absent from subset) (0 control-track, 0 non-headline); families tree hash 8c1d8bb8a7fc.
Descriptive only (mean over items; a collection is one item; no CIs). Headline = tracks family+collection, headline-eligible; control track (K01 RPROT, K02, E05) reported separately. Trivial baselines are rows of every table; skill = (s-b)/(1-b), b = per-family max(singletons, one big HOG, mean relabel). membership F1 is NA for families without complete cells or decoys. Confirmatory test: scripts/confirmatory.py.

## scope: headline

| method | variant | level_bcubed_f1 (n) | skill | merge_accuracy_node (n) | skill | event_recovery_0.5 (n) | skill | membership_f1 (n) | skill | S/D-bal merge acc (post-hoc) |
|---|---|---|---|---|---|---|---|---|---|---|
| positive_control | reference | 1.000 (2) | 1.000 | 1.000 (2) | 1.000 | 1.000 (2) | 1.000 | NA (0) | NA | 1.000 |
| all_singletons |  | 0.695 (2) | -0.808 | 0.128 (2) | -2.947 | 0.000 (2) | 0.000 | NA (0) | NA | 0.500 |
| one_big_hog |  | 0.790 (2) | -0.398 | 0.872 (2) | -0.084 | 0.000 (2) | 0.000 | NA (0) | NA | 0.500 |
| random_relabel | mean10 | 0.695 (2) | -0.806 | 0.129 (2) | -2.947 | 0.000 (2) | 0.000 | NA (0) | NA | 0.500 |
| fastoma | default_declared | 0.920 (2) | 0.469 | 0.942 (2) | 0.633 | 0.134 (2) | 0.134 | NA (0) | NA | 0.756 |
| fastoma | default_late | 0.920 (2) | 0.468 | 0.942 (2) | 0.630 | 0.139 (2) | 0.139 | NA (0) | NA | 0.755 |
| of_flat | default | 0.789 (2) | -0.367 | 0.871 (2) | -0.056 | 0.007 (2) | 0.007 | NA (0) | NA | 0.533 |
| orthofinder | default | 0.927 (2) | 0.424 | 0.878 (2) | 0.029 | 0.057 (2) | 0.057 | NA (0) | NA | 0.595 |

random relabel, min-max of item-mean over 10 seeds: level_bcubed_f1 0.694637-0.694782; merge_accuracy_node 0.12861-0.128928; event_recovery_0.5 0-0

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

- fastoma default_declared: 1.6% of units identical to all-singletons
- fastoma default_late: 1.6% of units identical to all-singletons
- of_flat default: 0.4% of units identical to all-singletons
- orthofinder default: 0.5% of units identical to all-singletons

## Positive control / degenerate baselines (per-family min/mean/max)

- positive_control level_bcubed_f1: mean 1.0000, min 1.0000, max 1.0000, n=1976
- positive_control merge_accuracy_node: mean 1.0000, min 1.0000, max 1.0000, n=1976
- positive_control event_recovery_0.5: mean 1.0000, min 1.0000, max 1.0000, n=202
- all_singletons level_bcubed_f1: mean 0.6052, min 0.4606, max 1.0000, n=1976
- all_singletons merge_accuracy_node: mean 0.1142, min 0.0000, max 1.0000, n=1976
- all_singletons event_recovery_0.5: mean 0.0000, min 0.0000, max 0.0000, n=202
- one_big_hog level_bcubed_f1: mean 0.9638, min 0.6152, max 1.0000, n=1976
- one_big_hog merge_accuracy_node: mean 0.8858, min 0.0000, max 1.0000, n=1976
- one_big_hog event_recovery_0.5: mean 0.0000, min 0.0000, max 0.0000, n=202
- random_relabel level_bcubed_f1: mean 0.6055, min 0.4606, max 1.0000, n=1976
- random_relabel merge_accuracy_node: mean 0.1146, min 0.0000, max 1.0000, n=1976
- random_relabel event_recovery_0.5: mean 0.0000, min 0.0000, max 0.0000, n=202

POSITIVE CONTROL: PERFECT on all primary endpoints.
