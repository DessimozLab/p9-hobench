# Robustness: Y frag (full mode), perturbed vs clean Y full

Perturbation: fragment, seed 20261026; panel Y_frag. Descriptive: mean over scored units (YGOB units + Y03; control track excluded), delta = perturbed - clean; CI = 95% bootstrap over units. See scripts/evaluate_robustness.py docstring for the rules.

## fastoma_default_declared

| credit | metric | n units | mean delta | 95% CI |
|---|---|---|---|---|
| either | level_bcubed_f1 | 1977 | -0.0112 | [-0.0135, -0.0090] |
| either | merge_accuracy_node | 1977 | -0.0238 | [-0.0298, -0.0178] |
| either | event_recovery_0.5 | 202 | +0.0149 | [-0.0248, +0.0545] |
| either | gene_frac_correct_affected | 1409 | -0.0158 | [-0.0268, -0.0050] |
| either | gene_frac_correct_unaffected | 1977 | -0.0154 | [-0.0220, -0.0094] |
| both | level_bcubed_f1 | 1977 | -0.0452 | [-0.0484, -0.0420] |
| both | merge_accuracy_node | 1977 | -0.1139 | [-0.1219, -0.1065] |
| both | event_recovery_0.5 | 202 | -0.0347 | [-0.0743, +0.0050] |
| both | gene_frac_correct_affected | 1409 | -0.2755 | [-0.2953, -0.2574] |
| both | gene_frac_correct_unaffected | 1977 | -0.0477 | [-0.0546, -0.0402] |

## of_flat_default

| credit | metric | n units | mean delta | 95% CI |
|---|---|---|---|---|
| either | level_bcubed_f1 | 1977 | +0.0015 | [+0.0009, +0.0022] |
| either | merge_accuracy_node | 1977 | +0.0047 | [+0.0029, +0.0067] |
| either | event_recovery_0.5 | 202 | +0.0000 | [-0.0149, +0.0149] |
| either | gene_frac_correct_affected | 1409 | +0.0415 | [+0.0320, +0.0514] |
| either | gene_frac_correct_unaffected | 1977 | +0.0034 | [+0.0012, +0.0054] |
| both | level_bcubed_f1 | 1977 | -0.0128 | [-0.0144, -0.0112] |
| both | merge_accuracy_node | 1977 | -0.0336 | [-0.0379, -0.0295] |
| both | event_recovery_0.5 | 202 | -0.0099 | [-0.0248, +0.0000] |
| both | gene_frac_correct_affected | 1409 | -0.2186 | [-0.2395, -0.1977] |
| both | gene_frac_correct_unaffected | 1977 | -0.0212 | [-0.0244, -0.0180] |

## orthofinder_default

| credit | metric | n units | mean delta | 95% CI |
|---|---|---|---|---|
| either | level_bcubed_f1 | 1977 | -0.0028 | [-0.0043, -0.0014] |
| either | merge_accuracy_node | 1977 | -0.0033 | [-0.0084, +0.0019] |
| either | event_recovery_0.5 | 202 | +0.0297 | [-0.0149, +0.0792] |
| either | gene_frac_correct_affected | 1409 | +0.0207 | [+0.0088, +0.0327] |
| either | gene_frac_correct_unaffected | 1977 | -0.0018 | [-0.0077, +0.0041] |
| both | level_bcubed_f1 | 1977 | -0.0209 | [-0.0232, -0.0186] |
| both | merge_accuracy_node | 1977 | -0.0505 | [-0.0575, -0.0438] |
| both | event_recovery_0.5 | 202 | -0.0297 | [-0.0693, +0.0099] |
| both | gene_frac_correct_affected | 1409 | -0.2526 | [-0.2761, -0.2305] |
| both | gene_frac_correct_unaffected | 1977 | -0.0306 | [-0.0372, -0.0239] |

Noise floor (OF default vs rep2 on clean Y): results/robustness/noise_floor_orthofinder.tsv

