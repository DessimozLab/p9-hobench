# Robustness: Y alttree (full mode), perturbed vs clean Y full

Perturbation: alttree, seed 20261026; panel Y_alttree. Descriptive: mean over scored units (YGOB units + Y03; control track excluded), delta = perturbed - clean; CI = 95% bootstrap over units. See scripts/evaluate_robustness.py docstring for the rules.

## fastoma_default_declared (induced levels: 1; shared: 14)

| credit | metric | n units | mean delta | 95% CI |
|---|---|---|---|---|
| - | level_bcubed_f1 | 1977 | -0.0029 | [-0.0046, -0.0013] |
| - | merge_accuracy_node | 1977 | -0.0035 | [-0.0070, +0.0001] |
| - | event_recovery_0.5 | 202 | +0.0050 | [-0.0248, +0.0347] |
| - | gene_frac_correct_affected | 1977 | -0.0050 | [-0.0095, -0.0005] |
| - | gene_frac_correct_unaffected | 1 | +0.0000 | NA |

## of_flat_default (induced levels: 1; shared: 14)

| credit | metric | n units | mean delta | 95% CI |
|---|---|---|---|---|
| - | level_bcubed_f1 | 1977 | -0.0037 | [-0.0040, -0.0033] |
| - | merge_accuracy_node | 1977 | -0.0002 | [-0.0005, +0.0000] |
| - | event_recovery_0.5 | 202 | +0.0000 | [+0.0000, +0.0000] |
| - | gene_frac_correct_affected | 1977 | -0.0003 | [-0.0008, +0.0000] |
| - | gene_frac_correct_unaffected | 1 | +0.0000 | NA |

## orthofinder_default (induced levels: 1; shared: 14)

| credit | metric | n units | mean delta | 95% CI |
|---|---|---|---|---|
| - | level_bcubed_f1 | 1977 | -0.0028 | [-0.0035, -0.0022] |
| - | merge_accuracy_node | 1977 | -0.0012 | [-0.0029, +0.0004] |
| - | event_recovery_0.5 | 202 | +0.0000 | [-0.0149, +0.0149] |
| - | gene_frac_correct_affected | 1977 | -0.0021 | [-0.0044, -0.0002] |
| - | gene_frac_correct_unaffected | 1 | +0.0000 | NA |

Noise floor (OF default vs rep2 on clean Y): results/robustness/noise_floor_orthofinder.tsv

