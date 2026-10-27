# Robustness: Y alttree (full mode), perturbed vs clean Y full

Perturbation: alttree, seed 20261026; panel Y_alttree. Descriptive: mean over scored units (YGOB units + Y03; control track excluded), delta = perturbed - clean; CI = 95% bootstrap over units. See scripts/evaluate_robustness.py docstring for the rules.

## orthofinder_default (induced levels: 1; shared: 14)

| credit | metric | n units | mean delta | 95% CI |
|---|---|---|---|---|
| - | level_bcubed_f1 | 1977 | -0.0028 | [-0.0034, -0.0022] |
| - | merge_accuracy_node | 1977 | -0.0012 | [-0.0031, +0.0004] |
| - | event_recovery_0.5 | 202 | +0.0000 | [-0.0149, +0.0149] |
| - | gene_frac_correct_affected | 1977 | -0.0021 | [-0.0043, -0.0002] |
| - | gene_frac_correct_unaffected | 1 | +0.0000 | NA |

Noise floor: OF rep2 replicate not available yet.

