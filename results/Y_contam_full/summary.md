# Robustness: Y contam (full mode), perturbed vs clean Y full

Perturbation: contaminate, seed 20261026; panel Y_contam. Descriptive: mean over scored units (YGOB units + Y03; control track excluded), delta = perturbed - clean; CI = 95% bootstrap over units. See scripts/evaluate_robustness.py docstring for the rules.

Contaminant placement (fastoma_default_declared): 12/20 contaminants share a HOG with reference genes at some level; 11/20 at the root level (details: contaminant_placement_fastoma_default_declared.tsv).

## fastoma_default_declared

| credit | metric | n units | mean delta | 95% CI |
|---|---|---|---|---|
| - | level_bcubed_f1 | 1977 | -0.0001 | [-0.0013, +0.0012] |
| - | merge_accuracy_node | 1977 | -0.0002 | [-0.0027, +0.0023] |
| - | event_recovery_0.5 | 202 | +0.0000 | [+0.0000, +0.0000] |
| - | gene_frac_correct_unaffected | 1967 | -0.0018 | [-0.0050, +0.0017] |
| - | gene_frac_correct_affected | 13 | -0.0196 | [-0.0587, +0.0000] |

Contaminant placement (of_flat_default): 8/20 contaminants share a HOG with reference genes at some level; 8/20 at the root level (details: contaminant_placement_of_flat_default.tsv).

## of_flat_default

| credit | metric | n units | mean delta | 95% CI |
|---|---|---|---|---|
| - | level_bcubed_f1 | 1977 | +0.0000 | [+0.0000, +0.0000] |
| - | merge_accuracy_node | 1977 | +0.0000 | [+0.0000, +0.0000] |
| - | event_recovery_0.5 | 202 | +0.0000 | [+0.0000, +0.0000] |
| - | gene_frac_correct_unaffected | 1969 | +0.0000 | [+0.0000, +0.0000] |
| - | gene_frac_correct_affected | 8 | +0.0000 | [+0.0000, +0.0000] |

Contaminant placement (orthofinder_default): 8/20 contaminants share a HOG with reference genes at some level; 8/20 at the root level (details: contaminant_placement_orthofinder_default.tsv).

## orthofinder_default

| credit | metric | n units | mean delta | 95% CI |
|---|---|---|---|---|
| - | level_bcubed_f1 | 1977 | -0.0000 | [-0.0000, +0.0000] |
| - | merge_accuracy_node | 1977 | +0.0001 | [+0.0000, +0.0003] |
| - | event_recovery_0.5 | 202 | +0.0050 | [+0.0000, +0.0149] |
| - | gene_frac_correct_unaffected | 1969 | +0.0000 | [+0.0000, +0.0001] |
| - | gene_frac_correct_affected | 8 | +0.0000 | [+0.0000, +0.0000] |

Noise floor (OF default vs rep2 on clean Y): results/robustness/noise_floor_orthofinder.tsv

