# Subset vs full equivalence (pre-registered; TOST margin +-0.05, 90% paired bootstrap CI)

Delta = subset - full, paired per reference family/unit (all tracks; collections are bootstrapped over their units). Equivalent = 90% CI inside +-0.05. Pair identity = fraction of member pairs with identical predicted merge level.

| panel | method | variant | metric | n | mean delta | 90% CI | equivalent | pair identity (n pairs) |
|---|---|---|---|---|---|---|---|---|
| Y | fastoma | default_declared | event_recovery_0.5 | 202 | +0.0000 | [+0.0000, +0.0000] | YES |  |
| Y | fastoma | default_declared | level_bcubed_f1 | 1976 | +0.0002 | [-0.0006, +0.0010] | YES | 0.9869 (165393) |
| Y | fastoma | default_declared | merge_accuracy_node | 1976 | +0.0000 | [-0.0022, +0.0022] | YES |  |
| Y | fastoma | default_late | event_recovery_0.5 | 202 | +0.0000 | [+0.0000, +0.0000] | YES |  |
| Y | fastoma | default_late | level_bcubed_f1 | 1976 | +0.0002 | [-0.0005, +0.0010] | YES |  |
| Y | fastoma | default_late | merge_accuracy_node | 1976 | +0.0002 | [-0.0020, +0.0023] | YES |  |
| Y | of_flat | default | event_recovery_0.5 | 202 | +0.0000 | [-0.0099, +0.0099] | YES |  |
| Y | of_flat | default | level_bcubed_f1 | 1976 | -0.0004 | [-0.0007, -0.0002] | YES | 0.9956 (165393) |
| Y | of_flat | default | merge_accuracy_node | 1976 | -0.0015 | [-0.0023, -0.0008] | YES |  |
| Y | orthofinder | default | event_recovery_0.5 | 202 | -0.0050 | [-0.0347, +0.0248] | YES |  |
| Y | orthofinder | default | level_bcubed_f1 | 1976 | -0.0003 | [-0.0008, +0.0002] | YES | 0.9785 (165393) |
| Y | orthofinder | default | merge_accuracy_node | 1976 | -0.0006 | [-0.0030, +0.0018] | YES |  |
| K | fastoma | default_declared | event_recovery_0.5 | 1 | +0.0000 | [+0.0000, +0.0000] | YES |  |
| K | fastoma | default_declared | level_bcubed_f1 | 41 | -0.0013 | [-0.0029, +0.0000] | YES | 0.9959 (6898) |
| K | fastoma | default_declared | membership_f1 | 50 | -0.0017 | [-0.0055, +0.0003] | YES |  |
| K | fastoma | default_declared | merge_accuracy_node | 41 | -0.0003 | [-0.0009, +0.0000] | YES |  |
| K | fastoma | default_late | event_recovery_0.5 | 1 | +0.0000 | [+0.0000, +0.0000] | YES |  |
| K | fastoma | default_late | level_bcubed_f1 | 41 | -0.0013 | [-0.0029, +0.0000] | YES |  |
| K | fastoma | default_late | membership_f1 | 50 | -0.0017 | [-0.0055, +0.0003] | YES |  |
| K | fastoma | default_late | merge_accuracy_node | 41 | -0.0003 | [-0.0009, +0.0000] | YES |  |
| K | of_flat | default | event_recovery_0.5 | 1 | +1.0000 | [+1.0000, +1.0000] | no |  |
| K | of_flat | default | level_bcubed_f1 | 41 | +0.0021 | [+0.0004, +0.0047] | YES | 0.9700 (6898) |
| K | of_flat | default | membership_f1 | 50 | +0.0073 | [-0.0028, +0.0205] | YES |  |
| K | of_flat | default | merge_accuracy_node | 41 | +0.0210 | [+0.0044, +0.0450] | YES |  |
| K | orthofinder | default | event_recovery_0.5 | 1 | +1.0000 | [+1.0000, +1.0000] | no |  |
| K | orthofinder | default | level_bcubed_f1 | 41 | +0.0015 | [-0.0002, +0.0041] | YES | 0.9613 (6898) |
| K | orthofinder | default | membership_f1 | 50 | +0.0073 | [-0.0028, +0.0201] | YES |  |
| K | orthofinder | default | merge_accuracy_node | 41 | +0.0201 | [+0.0037, +0.0441] | YES |  |

Notes: K event recovery has n=1 event (K02 TU_G_dup), so a delta of +1.0 is a single flip, not an estimate; treat as inconclusive. Y membership F1 is NA (no complete cells/decoys). Pair identity is computed once per tool (shown on the level_bcubed_f1 row). FastOMA variants: declared/late; dendroblast and OF variants other than default have no subset runs.
