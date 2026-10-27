# Robustness: Y contam (full mode), perturbed vs clean Y full

Perturbation: contaminate, seed 20261026; panel Y_contam200. Descriptive: mean over scored units (YGOB units + Y03; control track excluded), delta = perturbed - clean; CI = 95% bootstrap over units. See scripts/evaluate_robustness.py docstring for the rules.

Contaminant placement (fastoma_default_declared): 133/200 contaminants share a HOG with reference genes at some level; 107/200 at the root level (details: contaminant_placement_fastoma_default_declared.tsv).

## fastoma_default_declared

| credit | metric | n units | mean delta | 95% CI |
|---|---|---|---|---|
| - | level_bcubed_f1 | 1977 | -0.0005 | [-0.0017, +0.0006] |
| - | merge_accuracy_node | 1977 | -0.0010 | [-0.0036, +0.0015] |
| - | event_recovery_0.5 | 202 | +0.0000 | [+0.0000, +0.0000] |
| - | gene_frac_correct_unaffected | 1879 | -0.0014 | [-0.0045, +0.0017] |
| - | gene_frac_correct_affected | 105 | +0.0011 | [-0.0080, +0.0094] |

Contaminant placement (orthofinder_default): 103/200 contaminants share a HOG with reference genes at some level; 103/200 at the root level (details: contaminant_placement_orthofinder_default.tsv).

## orthofinder_default

| credit | metric | n units | mean delta | 95% CI |
|---|---|---|---|---|
| - | level_bcubed_f1 | 1977 | -0.0001 | [-0.0003, +0.0001] |
| - | merge_accuracy_node | 1977 | +0.0004 | [-0.0003, +0.0012] |
| - | event_recovery_0.5 | 202 | +0.0000 | [+0.0000, +0.0000] |
| - | gene_frac_correct_unaffected | 1883 | +0.0001 | [+0.0000, +0.0004] |
| - | gene_frac_correct_affected | 99 | +0.0044 | [-0.0123, +0.0227] |

Noise floor (OF default vs rep2 on clean Y): results/robustness/noise_floor_orthofinder.tsv

