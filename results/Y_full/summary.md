# Panel Y / full: evaluation summary

hobench git 4eb0772ad5 (dirty=True); 1977 families, 0 non-headline; families tree hash afa5a53647f2.
Mean over items (collection = one item), scope=all families; CI = stratified bootstrap, none for n<5 (then min-max). Details: summary.tsv (also scope=headline and per-collection tracks).

| method | variant | level_bcubed_f1 | merge_accuracy_node | event_recovery_0.5 | membership_f1 | n items |
|---|---|---|---|---|---|---|
| positive_control | reference | 1.000 (1.00-1.00) | 1.000 (1.00-1.00) | 1.000 (1.00-1.00) | 1.000 (1.00-1.00) | 2 |
| all_singletons |  | 0.695 (0.61-0.78) | 0.128 (0.11-0.14) | 0.000 (0.00-0.00) | 0.001 (0.00-0.00) | 2 |
| one_hog_per_family |  | 0.790 (0.62-0.96) | 0.872 (0.86-0.89) | 0.000 (0.00-0.00) | 0.836 (0.67-1.00) | 2 |
| random_relabel | mean10 | 0.694 (0.61-0.78) | 0.129 (0.11-0.14) | 0.000 (0.00-0.00) | 0.216 (0.15-0.29) | 2 |
| fastoma | default_declared | 0.920 (0.87-0.97) | 0.942 (0.90-0.98) | 0.134 (0.00-0.27) | 0.918 (0.89-0.95) | 2 |
| fastoma | default_late | 0.920 (0.87-0.97) | 0.942 (0.90-0.98) | 0.139 (0.00-0.28) | 0.918 (0.89-0.95) | 2 |
| of_flat | default | 0.790 (0.62-0.96) | 0.872 (0.86-0.89) | 0.007 (0.00-0.01) | 0.833 (0.67-0.99) | 2 |
| of_flat | dendroblast | 0.790 (0.62-0.96) | 0.872 (0.86-0.89) | 0.010 (0.00-0.02) | 0.833 (0.67-0.99) | 2 |
| orthofinder | default | 0.928 (0.88-0.98) | 0.879 (0.86-0.90) | 0.060 (0.00-0.12) | 0.831 (0.67-0.99) | 2 |
| orthofinder | dendroblast | 0.928 (0.88-0.98) | 0.880 (0.86-0.90) | 0.072 (0.00-0.14) | 0.831 (0.67-0.99) | 2 |

Unit-level (families track + each collection's units; mean over units):

| method | variant | collection | level_bcubed_f1 | merge_accuracy_node | event_recovery_0.5 | membership_f1 | n units |
|---|---|---|---|---|---|---|---|
| positive_control | reference | YGOB | 1.000 | 1.000 | 1.000 | 1.000 | 1976 |
| all_singletons |  | YGOB | 0.605 | 0.114 | 0.000 | 0.001 | 1976 |
| one_hog_per_family |  | YGOB | 0.964 | 0.886 | 0.000 | 0.999 | 1976 |
| random_relabel | mean10 | YGOB | 0.605 | 0.115 | 0.000 | 0.146 | 1976 |
| fastoma | default_declared | YGOB | 0.967 | 0.900 | 0.269 | 0.947 | 1976 |
| fastoma | default_late | YGOB | 0.967 | 0.900 | 0.279 | 0.946 | 1976 |
| of_flat | default | YGOB | 0.964 | 0.887 | 0.015 | 0.992 | 1976 |
| of_flat | dendroblast | YGOB | 0.964 | 0.888 | 0.020 | 0.992 | 1976 |
| orthofinder | default | YGOB | 0.977 | 0.900 | 0.119 | 0.988 | 1976 |
| orthofinder | dendroblast | YGOB | 0.977 | 0.904 | 0.144 | 0.988 | 1976 |

## Sanity checks (per-family min/max over families with a defined value)

- positive_control level_bcubed_f1: mean 1.0000, min 1.0000, max 1.0000, n=1977
- positive_control merge_accuracy_node: mean 1.0000, min 1.0000, max 1.0000, n=1977
- positive_control event_recovery_0.5: mean 1.0000, min 1.0000, max 1.0000, n=202
- positive_control membership_f1: mean 1.0000, min 1.0000, max 1.0000, n=1977
- all_singletons level_bcubed_f1: mean 0.6052, min 0.4606, max 1.0000, n=1977
- all_singletons merge_accuracy_node: mean 0.1141, min 0.0000, max 1.0000, n=1977
- all_singletons event_recovery_0.5: mean 0.0000, min 0.0000, max 0.0000, n=202
- all_singletons membership_f1: mean 0.0010, min 0.0000, max 1.0000, n=1977
- one_hog_per_family level_bcubed_f1: mean 0.9638, min 0.6152, max 1.0000, n=1977
- one_hog_per_family merge_accuracy_node: mean 0.8859, min 0.0000, max 1.0000, n=1977
- one_hog_per_family event_recovery_0.5: mean 0.0000, min 0.0000, max 0.0000, n=202
- one_hog_per_family membership_f1: mean 0.9991, min 0.2222, max 1.0000, n=1977
- random_relabel level_bcubed_f1: mean 0.6054, min 0.4606, max 1.0000, n=1977
- random_relabel merge_accuracy_node: mean 0.1146, min 0.0000, max 1.0000, n=1977
- random_relabel event_recovery_0.5: mean 0.0000, min 0.0000, max 0.0000, n=202
- random_relabel membership_f1: mean 0.1464, min 0.0952, max 1.0000, n=1977

POSITIVE CONTROL: PERFECT on all primary endpoints.
