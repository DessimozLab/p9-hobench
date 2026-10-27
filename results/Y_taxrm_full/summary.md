# Robustness: Y taxrm (full mode), perturbed vs clean Y full

Perturbation: remove, seed 20261026; panel Y_taxrm. Descriptive: mean over scored units (YGOB units + Y03; control track excluded), delta = perturbed - clean; CI = 95% bootstrap over units. See scripts/evaluate_robustness.py docstring for the rules.

## fastoma_default_declared

| credit | metric | n units | mean delta | 95% CI |
|---|---|---|---|---|
| - | level_bcubed_f1 | 1972 | -0.0006 | [-0.0020, +0.0010] |
| - | merge_accuracy_node | 1972 | -0.0028 | [-0.0068, +0.0015] |
| - | event_recovery_0.5 | 202 | +0.0050 | [-0.0248, +0.0396] |
| - | gene_frac_correct_affected | 1963 | -0.0057 | [-0.0108, -0.0008] |
| - | gene_frac_correct_unaffected | 1959 | -0.0046 | [-0.0099, +0.0010] |

## of_flat_default

| credit | metric | n units | mean delta | 95% CI |
|---|---|---|---|---|
| - | level_bcubed_f1 | 1972 | -0.0002 | [-0.0007, +0.0002] |
| - | merge_accuracy_node | 1972 | -0.0004 | [-0.0020, +0.0011] |
| - | event_recovery_0.5 | 202 | +0.0000 | [+0.0000, +0.0000] |
| - | gene_frac_correct_affected | 1963 | -0.0010 | [-0.0032, +0.0014] |
| - | gene_frac_correct_unaffected | 1959 | -0.0005 | [-0.0017, +0.0006] |

## orthofinder_default

| credit | metric | n units | mean delta | 95% CI |
|---|---|---|---|---|
| - | level_bcubed_f1 | 1972 | -0.0008 | [-0.0019, +0.0003] |
| - | merge_accuracy_node | 1972 | +0.0011 | [-0.0039, +0.0060] |
| - | event_recovery_0.5 | 202 | +0.0050 | [-0.0396, +0.0495] |
| - | gene_frac_correct_affected | 1963 | +0.0016 | [-0.0049, +0.0077] |
| - | gene_frac_correct_unaffected | 1959 | -0.0020 | [-0.0055, +0.0011] |

