# HOBench main table (headline track)

Source: results/<P>_<mode>/per_family.tsv (headline = track family+collection, headline-eligible; controls K01/K02/E05 excluded).
Value = mean over items (a collection = one item = mean of its units); n = items with a defined value; skill = (s-b)/(1-b), b = per-family max of trivial baselines (mean over items). NA = no item defines the endpoint. Descriptive only; no CIs (n<5 almost everywhere). Panels: Y full, K full, E subset, P subset, V subset

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

## E (subset)

| method | panels | level_bcubed_f1 (n) | skill | merge_accuracy_node (n) | skill | event_recovery_0.5 (n) | skill | membership_f1 (n) | skill |
|---|---|---|---|---|---|---|---|---|---|
| orthofinder/default | E | 0.981 (5) | 0.871 | 0.729 (5) | -1.235 | 0.400 (5) | 0.400 | 0.870 (3) | 0.840 |
| fastoma/default_declared | E | 0.906 (5) | 0.371 | 0.531 (5) | -3.153 | 0.025 (5) | 0.025 | 0.563 (3) | 0.466 |
| fastoma/default_late | E | 0.904 (5) | 0.342 | 0.531 (5) | -3.153 | 0.025 (5) | 0.025 | 0.563 (3) | 0.466 |
| of_flat/default | E | 0.358 (5) | -9.673 | 0.758 (5) | -0.716 | 0.000 (5) | 0.000 | 0.152 (3) | -0.040 |
| all_singletons | E | 0.849 (5) | 0.000 | 0.146 (5) | -5.600 | 0.000 (5) | 0.000 | 0.000 (3) | -0.223 |
| one_big_hog | E | 0.217 (5) | -10.447 | 0.838 (5) | 0.000 | 0.000 (5) | 0.000 | 0.014 (3) | -0.206 |
| random_relabel/mean10 | E | 0.835 (5) | -0.210 | 0.147 (5) | -5.590 | 0.000 (5) | 0.000 | 0.180 (3) | 0.000 |
| positive_control/reference | E | 1.000 (5) | 1.000 | 1.000 (5) | 1.000 | 1.000 (5) | 1.000 | 1.000 (3) | 1.000 |

## P (subset)

| method | panels | level_bcubed_f1 (n) | skill | merge_accuracy_node (n) | skill | event_recovery_0.5 (n) | skill | membership_f1 (n) | skill |
|---|---|---|---|---|---|---|---|---|---|
| orthofinder/default | P | 0.994 (3) | 0.920 | 0.907 (3) | 0.767 | 0.667 (2) | 0.667 | 0.493 (1) | -0.127 |
| fastoma/default_declared | P | 0.997 (3) | 0.918 | 0.867 (3) | 0.667 | 0.667 (2) | 0.667 | 0.049 (1) | -1.116 |
| fastoma/default_late | P | 0.997 (3) | 0.918 | 0.867 (3) | 0.667 | 0.667 (2) | 0.667 | 0.049 (1) | -1.116 |
| of_flat/default | P | 0.481 (3) | -7.988 | 0.579 (3) | -0.039 | 0.167 (2) | 0.167 | 0.493 (1) | -0.127 |
| all_singletons | P | 0.941 (3) | -0.260 | 0.550 (3) | -0.250 | 0.000 (2) | 0.000 | 0.000 (1) | -1.224 |
| one_big_hog | P | 0.148 (3) | -12.665 | 0.233 (3) | -0.167 | 0.000 (2) | 0.000 | 0.550 (1) | 0.000 |
| random_relabel/mean10 | P | 0.916 (3) | -0.339 | 0.616 (3) | -0.035 | 0.000 (2) | 0.000 | 0.277 (1) | -0.607 |
| positive_control/reference | P | 1.000 (3) | 1.000 | 1.000 (3) | 1.000 | 1.000 (2) | 1.000 | 1.000 (1) | 1.000 |

## V (subset)

| method | panels | level_bcubed_f1 (n) | skill | merge_accuracy_node (n) | skill | event_recovery_0.5 (n) | skill | membership_f1 (n) | skill |
|---|---|---|---|---|---|---|---|---|---|
| orthofinder/default | V | 0.987 (3) | 0.962 | 0.960 (3) | 0.514 | 0.907 (2) | 0.907 | 0.721 (1) | 0.595 |
| fastoma/default_declared | V | 0.881 (3) | 0.816 | 0.707 (3) | -0.905 | 0.204 (2) | 0.202 | 0.511 (1) | 0.290 |
| fastoma/default_late | V | 0.882 (3) | 0.818 | 0.709 (3) | -0.886 | 0.222 (2) | 0.221 | 0.511 (1) | 0.290 |
| of_flat/default | V | 0.664 (3) | -0.115 | 0.903 (3) | -0.271 | 0.056 (2) | 0.054 | 0.721 (1) | 0.595 |
| all_singletons | V | 0.517 (3) | -0.007 | 0.016 (3) | -7.685 | 0.000 (2) | -0.002 | 0.000 (1) | -0.453 |
| one_big_hog | V | 0.479 (3) | -0.643 | 0.924 (3) | 0.000 | 0.000 (2) | -0.002 | 0.312 (1) | 0.000 |
| random_relabel/mean10 | V | 0.512 (3) | -0.018 | 0.024 (3) | -7.577 | 0.002 (2) | 0.000 | 0.115 (1) | -0.286 |
| positive_control/reference | V | 1.000 (3) | 1.000 | 1.000 (3) | 1.000 | 1.000 (2) | 1.000 | 1.000 (1) | 1.000 |

## Pooled headline means (all panels above, each in its analysis mode; items pooled, unweighted)

| method | panels | level_bcubed_f1 (n) | skill | merge_accuracy_node (n) | skill | event_recovery_0.5 (n) | skill | membership_f1 (n) | skill |
|---|---|---|---|---|---|---|---|---|---|
| orthofinder/default | EPVY | 0.977 (13) | 0.816 | 0.846 (13) | -0.323 | 0.479 (11) | 0.479 | 0.765 (5) | 0.597 |
| orthofinder/dendroblast | Y | 0.928 (2) | 0.425 | 0.880 (2) | 0.047 | 0.072 (2) | 0.072 | NA (0) | NA |
| fastoma/default_declared | EPVY | 0.923 (13) | 0.569 | 0.712 (13) | -1.361 | 0.194 (11) | 0.194 | 0.450 (5) | 0.115 |
| fastoma/default_late | EPVY | 0.923 (13) | 0.556 | 0.713 (13) | -1.359 | 0.198 (11) | 0.198 | 0.450 (5) | 0.115 |
| of_flat/default | EPVY | 0.523 (13) | -5.936 | 0.768 (13) | -0.391 | 0.042 (11) | 0.041 | 0.334 (5) | 0.069 |
| of_flat/dendroblast | Y | 0.790 (2) | -0.358 | 0.872 (2) | -0.047 | 0.010 (2) | 0.010 | NA (0) | NA |
| all_singletons | EPVY | 0.770 (13) | -0.195 | 0.206 (13) | -4.524 | 0.000 (11) | -0.000 | 0.000 (5) | -0.469 |
| one_big_hog | EPVY | 0.350 (13) | -7.241 | 0.724 (13) | -0.046 | 0.000 (11) | -0.000 | 0.181 (5) | -0.124 |
| random_relabel/mean10 | EPVY | 0.758 (13) | -0.307 | 0.224 (13) | -4.461 | 0.000 (11) | 0.000 | 0.186 (5) | -0.179 |
| positive_control/reference | EPVY | 1.000 (13) | 1.000 | 1.000 (13) | 1.000 | 1.000 (11) | 1.000 | 1.000 (5) | 1.000 |

Note: modes are never pooled within a panel; pooling across panels mixes full (Y,K) and subset (E,P,V) modes and is descriptive only (confirmatory test: scripts/confirmatory.py).
