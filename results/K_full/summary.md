# Panel K / full: evaluation summary

hobench git 86b9ad7ec5 (dirty=True); 51 families, 4 non-headline; families tree hash a62dac6b1989.
Mean over items (collection = one item), scope=all families; CI = stratified bootstrap, none for n<5 (then min-max). Details: summary.tsv (also scope=headline and per-collection tracks).

| method | variant | level_bcubed_f1 | merge_accuracy_node | event_recovery_0.5 | membership_f1 | n items |
|---|---|---|---|---|---|---|
| positive_control | reference | 1.000 (1.00-1.00) | 1.000 (1.00-1.00) | 1.000 (1.00-1.00) | 1.000 (1.00-1.00) | 2 |
| all_singletons |  | 0.611 (0.51-0.72) | 0.167 (0.00-0.33) | 0.000 (0.00-0.00) | 0.000 (0.00-0.00) | 2 |
| one_hog_per_family |  | 0.764 (0.53-1.00) | 0.833 (0.67-1.00) | 0.000 (0.00-0.00) | 0.492 (0.05-0.75) | 2 |
| random_relabel | mean10 | 0.618 (0.52-0.72) | 0.182 (0.02-0.34) | 0.000 (0.00-0.00) | 0.144 (0.09-0.22) | 2 |
| fastoma | default_declared | 0.969 (0.96-0.98) | 0.676 (0.66-0.70) | 0.000 (0.00-0.00) | 0.780 (0.62-0.90) | 2 |
| fastoma | default_late | 0.969 (0.96-0.98) | 0.676 (0.66-0.70) | 0.000 (0.00-0.00) | 0.780 (0.62-0.90) | 2 |
| of_flat | default | 0.634 (0.53-0.73) | 0.408 (0.31-0.51) | 0.000 (0.00-0.00) | 0.723 (0.64-0.88) | 2 |
| of_flat | dendroblast | 0.634 (0.53-0.73) | 0.408 (0.31-0.51) | 0.000 (0.00-0.00) | 0.729 (0.64-0.89) | 2 |
| orthofinder | default | 0.980 (0.96-1.00) | 0.883 (0.79-0.98) | 0.000 (0.00-0.00) | 0.939 (0.88-0.99) | 2 |
| orthofinder | dendroblast | 0.985 (0.97-1.00) | 0.894 (0.81-0.98) | 0.000 (0.00-0.00) | 0.945 (0.89-0.99) | 2 |

Unit-level (families track + each collection's units; mean over units):

| method | variant | collection | level_bcubed_f1 | merge_accuracy_node | event_recovery_0.5 | membership_f1 | n units |
|---|---|---|---|---|---|---|---|
| positive_control | reference | K04_isofunctional_decoys | NA | NA | NA | 1.000 | 0 |
| positive_control | reference | RPROT | 1.000 | 1.000 | NA | 1.000 | 40 |
| all_singletons |  | K04_isofunctional_decoys | NA | NA | NA | 0.000 | 0 |
| all_singletons |  | RPROT | 0.505 | 0.000 | NA | 0.000 | 40 |
| one_hog_per_family |  | K04_isofunctional_decoys | NA | NA | NA | 0.753 | 0 |
| one_hog_per_family |  | RPROT | 1.000 | 1.000 | NA | 0.047 | 40 |
| random_relabel | mean10 | K04_isofunctional_decoys | NA | NA | NA | 0.221 | 0 |
| random_relabel | mean10 | RPROT | 0.518 | 0.019 | NA | 0.091 | 40 |
| fastoma | default_declared | K04_isofunctional_decoys | NA | NA | NA | 0.897 | 0 |
| fastoma | default_declared | RPROT | 0.956 | 0.697 | NA | 0.818 | 40 |
| fastoma | default_late | K04_isofunctional_decoys | NA | NA | NA | 0.897 | 0 |
| fastoma | default_late | RPROT | 0.956 | 0.697 | NA | 0.818 | 40 |
| of_flat | default | K04_isofunctional_decoys | NA | NA | NA | 0.649 | 0 |
| of_flat | default | RPROT | 0.534 | 0.306 | NA | 0.638 | 40 |
| of_flat | dendroblast | K04_isofunctional_decoys | NA | NA | NA | 0.656 | 0 |
| of_flat | dendroblast | RPROT | 0.534 | 0.306 | NA | 0.638 | 40 |
| orthofinder | default | K04_isofunctional_decoys | NA | NA | NA | 0.949 | 0 |
| orthofinder | default | RPROT | 0.996 | 0.976 | NA | 0.985 | 40 |
| orthofinder | dendroblast | K04_isofunctional_decoys | NA | NA | NA | 0.956 | 0 |
| orthofinder | dendroblast | RPROT | 0.996 | 0.976 | NA | 0.985 | 40 |

## Sanity checks (per-family min/max over families with a defined value)

- positive_control level_bcubed_f1: mean 1.0000, min 1.0000, max 1.0000, n=41
- positive_control merge_accuracy_node: mean 1.0000, min 1.0000, max 1.0000, n=41
- positive_control event_recovery_0.5: mean 1.0000, min 1.0000, max 1.0000, n=1
- positive_control membership_f1: mean 1.0000, min 1.0000, max 1.0000, n=51
- all_singletons level_bcubed_f1: mean 0.5103, min 0.4838, max 0.7176, n=41
- all_singletons merge_accuracy_node: mean 0.0081, min 0.0000, max 0.3333, n=41
- all_singletons event_recovery_0.5: mean 0.0000, min 0.0000, max 0.0000, n=1
- all_singletons membership_f1: mean 0.0000, min 0.0000, max 0.0000, n=51
- one_hog_per_family level_bcubed_f1: mean 0.9885, min 0.5271, max 1.0000, n=41
- one_hog_per_family merge_accuracy_node: mean 0.9919, min 0.6667, max 1.0000, n=41
- one_hog_per_family event_recovery_0.5: mean 0.0000, min 0.0000, max 0.0000, n=1
- one_hog_per_family membership_f1: mean 0.1980, min 0.0442, max 0.9333, n=51
- random_relabel level_bcubed_f1: mean 0.5225, min 0.4971, max 0.7173, n=41
- random_relabel merge_accuracy_node: mean 0.0270, min 0.0106, max 0.3447, n=41
- random_relabel event_recovery_0.5: mean 0.0000, min 0.0000, max 0.0000, n=1
- random_relabel membership_f1: mean 0.1172, min 0.0736, max 0.3212, n=51

POSITIVE CONTROL: PERFECT on all primary endpoints.
