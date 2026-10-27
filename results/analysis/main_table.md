# HOBench main table (headline track)

Source: results/<P>_<mode>/per_family.tsv (headline = track family+collection, headline-eligible; controls K01/K02/E05 excluded).
Value = mean over items (a collection = one item = mean of its units); n = items with a defined value; skill = (s-b)/(1-b), b = per-family max of trivial baselines (mean over items). NA = no item defines the endpoint. Descriptive only; no CIs (n<5 almost everywhere). Panels: Y full, K full, E subset, P subset

## Y (full)

| method | panels | level_bcubed_f1 (n) | skill | merge_accuracy_node (n) | skill | event_recovery_0.5 (n) | skill | membership_f1 (n) | skill |
|---|---|---|---|---|---|---|---|---|---|
| orthofinder/default | Y | 0.928 (2) | 0.425 | 0.879 (2) | 0.031 | 0.060 (2) | 0.060 | NA (0) | NA |
| orthofinder/dendroblast | Y | 0.928 (2) | 0.425 | 0.880 (2) | 0.047 | 0.072 (2) | 0.072 | NA (0) | NA |
| fastoma/default_declared | Y | 0.920 (2) | 0.470 | 0.942 (2) | 0.633 | 0.134 (2) | 0.134 | NA (0) | NA |
| fastoma/default_late | Y | 0.920 (2) | 0.468 | 0.942 (2) | 0.628 | 0.139 (2) | 0.139 | NA (0) | NA |
| of_flat/default | Y | 0.790 (2) | -0.360 | 0.872 (2) | -0.049 | 0.007 (2) | 0.007 | NA (0) | NA |
| of_flat/dendroblast | Y | 0.790 (2) | -0.358 | 0.872 (2) | -0.047 | 0.010 (2) | 0.010 | NA (0) | NA |
| all_singletons | Y | 0.695 (2) | -0.808 | 0.128 (2) | -2.947 | 0.000 (2) | 0.000 | NA (0) | NA |
| one_big_hog | Y | 0.790 (2) | -0.398 | 0.872 (2) | -0.084 | 0.000 (2) | 0.000 | NA (0) | NA |
| random_relabel/mean10 | Y | 0.694 (2) | -0.808 | 0.129 (2) | -2.947 | 0.000 (2) | 0.000 | NA (0) | NA |
| positive_control/reference | Y | 1.000 (2) | 1.000 | 1.000 (2) | 1.000 | 1.000 (2) | 1.000 | NA (0) | NA |

## K (full)

| method | panels | level_bcubed_f1 (n) | skill | merge_accuracy_node (n) | skill | event_recovery_0.5 (n) | skill | membership_f1 (n) | skill |
|---|---|---|---|---|---|---|---|---|---|
| orthofinder/default | K | NA (0) | NA | NA (0) | NA | NA (0) | NA | 0.949 (10) | 0.730 |
| orthofinder/dendroblast | K | NA (0) | NA | NA (0) | NA | NA (0) | NA | 0.956 (10) | 0.772 |
| fastoma/default_declared | K | NA (0) | NA | NA (0) | NA | NA (0) | NA | 0.897 (10) | 0.595 |
| fastoma/default_late | K | NA (0) | NA | NA (0) | NA | NA (0) | NA | 0.897 (10) | 0.595 |
| of_flat/default | K | NA (0) | NA | NA (0) | NA | NA (0) | NA | 0.949 (10) | 0.730 |
| of_flat/dendroblast | K | NA (0) | NA | NA (0) | NA | NA (0) | NA | 0.956 (10) | 0.772 |
| all_singletons | K | NA (0) | NA | NA (0) | NA | NA (0) | NA | 0.000 (10) | -4.595 |
| one_big_hog | K | NA (0) | NA | NA (0) | NA | NA (0) | NA | 0.753 (10) | 0.000 |
| random_relabel/mean10 | K | NA (0) | NA | NA (0) | NA | NA (0) | NA | 0.221 (10) | -3.333 |
| positive_control/reference | K | NA (0) | NA | NA (0) | NA | NA (0) | NA | 1.000 (10) | 1.000 |

## E (subset)

| method | panels | level_bcubed_f1 (n) | skill | merge_accuracy_node (n) | skill | event_recovery_0.5 (n) | skill | membership_f1 (n) | skill |
|---|---|---|---|---|---|---|---|---|---|
| orthofinder/default | E | 0.981 (5) | 0.871 | 0.729 (5) | -1.235 | 0.400 (5) | 0.400 | 0.741 (3) | 0.701 |
| fastoma/default_declared | E | 0.906 (5) | 0.371 | 0.531 (5) | -3.153 | 0.025 (5) | 0.025 | 0.574 (3) | 0.499 |
| fastoma/default_late | E | 0.904 (5) | 0.342 | 0.531 (5) | -3.153 | 0.025 (5) | 0.025 | 0.574 (3) | 0.499 |
| of_flat/default | E | 0.358 (5) | -9.673 | 0.758 (5) | -0.716 | 0.000 (5) | 0.000 | 0.289 (3) | 0.180 |
| all_singletons | E | 0.849 (5) | 0.000 | 0.146 (5) | -5.600 | 0.000 (5) | 0.000 | 0.075 (3) | -0.059 |
| one_big_hog | E | 0.217 (5) | -10.447 | 0.838 (5) | 0.000 | 0.000 (5) | 0.000 | 0.049 (3) | -0.095 |
| random_relabel/mean10 | E | 0.823 (5) | -0.409 | 0.149 (5) | -5.574 | 0.000 (5) | 0.000 | 0.120 (3) | -0.010 |
| positive_control/reference | E | 1.000 (5) | 1.000 | 1.000 (5) | 1.000 | 1.000 (5) | 1.000 | 1.000 (3) | 1.000 |

## P (subset)

| method | panels | level_bcubed_f1 (n) | skill | merge_accuracy_node (n) | skill | event_recovery_0.5 (n) | skill | membership_f1 (n) | skill |
|---|---|---|---|---|---|---|---|---|---|
| orthofinder/default | P | 0.994 (3) | 0.920 | 0.907 (3) | 0.767 | 0.667 (2) | 0.667 | 0.364 (1) | -1.096 |
| fastoma/default_declared | P | 0.997 (3) | 0.918 | 0.867 (3) | 0.667 | 0.667 (2) | 0.667 | 0.018 (1) | -2.235 |
| fastoma/default_late | P | 0.997 (3) | 0.918 | 0.867 (3) | 0.667 | 0.667 (2) | 0.667 | 0.018 (1) | -2.235 |
| of_flat/default | P | 0.481 (3) | -7.990 | 0.579 (3) | -0.039 | 0.167 (2) | 0.167 | 0.364 (1) | -1.096 |
| all_singletons | P | 0.941 (3) | -0.261 | 0.550 (3) | -0.250 | 0.000 (2) | 0.000 | 0.000 (1) | -2.296 |
| one_big_hog | P | 0.148 (3) | -12.667 | 0.233 (3) | -0.167 | 0.000 (2) | 0.000 | 0.697 (1) | 0.000 |
| random_relabel/mean10 | P | 0.908 (3) | -0.455 | 0.616 (3) | -0.035 | 0.000 (2) | 0.000 | 0.395 (1) | -0.994 |
| positive_control/reference | P | 1.000 (3) | 1.000 | 1.000 (3) | 1.000 | 1.000 (2) | 1.000 | 1.000 (1) | 1.000 |

## Pooled headline means (all panels above, each in its analysis mode; items pooled, unweighted)

| method | panels | level_bcubed_f1 (n) | skill | merge_accuracy_node (n) | skill | event_recovery_0.5 (n) | skill | membership_f1 (n) | skill |
|---|---|---|---|---|---|---|---|---|---|
| orthofinder/default | EKPY | 0.974 (10) | 0.783 | 0.812 (10) | -0.509 | 0.384 (9) | 0.384 | 0.862 (14) | 0.593 |
| orthofinder/dendroblast | KY | 0.928 (2) | 0.425 | 0.880 (2) | 0.047 | 0.072 (2) | 0.072 | 0.956 (10) | 0.772 |
| fastoma/default_declared | EKPY | 0.936 (10) | 0.514 | 0.714 (10) | -1.463 | 0.192 (9) | 0.192 | 0.765 (14) | 0.373 |
| fastoma/default_late | EKPY | 0.935 (10) | 0.498 | 0.714 (10) | -1.464 | 0.193 (9) | 0.193 | 0.765 (14) | 0.373 |
| of_flat/default | EKPY | 0.481 (10) | -7.229 | 0.727 (10) | -0.417 | 0.039 (9) | 0.039 | 0.766 (14) | 0.482 |
| of_flat/dendroblast | KY | 0.790 (2) | -0.358 | 0.872 (2) | -0.047 | 0.010 (2) | 0.010 | 0.956 (10) | 0.772 |
| all_singletons | EKPY | 0.846 (10) | -0.237 | 0.264 (10) | -3.822 | 0.000 (9) | 0.000 | 0.016 (14) | -3.459 |
| one_big_hog | EKPY | 0.311 (10) | -8.707 | 0.663 (10) | -0.056 | 0.000 (9) | 0.000 | 0.598 (14) | -0.020 |
| random_relabel/mean10 | EKPY | 0.823 (10) | -0.508 | 0.285 (10) | -3.760 | 0.000 (9) | 0.000 | 0.212 (14) | -2.453 |
| positive_control/reference | EKPY | 1.000 (10) | 1.000 | 1.000 (10) | 1.000 | 1.000 (9) | 1.000 | 1.000 (14) | 1.000 |

Note: modes are never pooled within a panel; pooling across panels mixes full (Y,K) and subset (E,P,V) modes and is descriptive only (confirmatory test: scripts/confirmatory.py).
