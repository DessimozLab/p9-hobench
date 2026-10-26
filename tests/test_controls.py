"""Positive / negative / dose-response controls."""
import numpy as np
import pytest

from hobench import controls
from hobench.metrics import PRIMARY, evaluate_family, family_scalars, panel_context
from hobench.stats import aggregate
from hobench.reference import validate_family
from hobench.seeds import SEEDS


def _all(fams, st, pred):
    ctx = panel_context(fams, st)
    return [evaluate_family(f, st, pred, ctx) for f in fams]


def test_positive_control_is_perfect(st, fams):
    pred = controls.reference_prediction(fams, st)
    assert pred.meta["inconsistent_pairs"] == 0
    assert pred.check(st) == []
    for r in _all(fams, st, pred):
        sc = family_scalars(r)
        for k in ("merge_accuracy_node", "merge_accuracy_pairs", "level_bcubed_f1", "level_pair_f1_pooled",
                  "membership_f1", "membership_precision", "membership_recall", "orthobench_f",
                  "gene_mean_frac_correct", "gene_frac_perfect"):
            assert sc[k] == 1.0, (r["id"], k, sc[k])
        assert sc["wf_ortholog_f1"] in (1.0, None)
        assert sc["contamination"] == 0.0 and sc["fission_pct"] == 0.0 and sc["fusion_pct"] == 0.0
        assert sc["exact_pct"] == 100.0 and sc["missing_genes_pct"] == 0.0 and sc["entropy"] == 0.0
        assert sc["event_recovery_0.5"] in (1.0, None) and sc["event_recovery_0.8"] in (1.0, None)
        assert sc["spurious_dup_rate"] in (0.0, None)
        assert set(r["pair_merge"]["error_histogram"]) == {"0"}


def test_negative_controls_known_values(st, fams):
    """Per-family merge accuracies (derivations in test_metrics.py).
    all-singletons, pair-pooled: F1 0, F2 0, F3 0, F4 4/6; node-balanced: F4 = (root D 1 + 2 S
    nodes 0)/3. one-big-HOG, pooled: F1 9/10, F2 11/15, F3 1, F4 2/6; node-balanced:
    F1 3/4, F2 4/5, F3 1, F4 2/3 (root D wrong, both S right). Membership F1: singletons 0
    everywhere (B1); one-big 1 except F4 (two RefOGs fused: P = 1/2, R = 1 -> 2/3)."""
    rs = _all(fams, st, controls.all_singletons(fams))
    assert [r["pair_merge"]["accuracy_pairs"] for r in rs] == pytest.approx([0, 0, 0, 4 / 6])
    assert [r["pair_merge"]["accuracy_node_balanced"] for r in rs] == pytest.approx([0, 0, 0, 1 / 3])
    assert [r["integrity"]["membership_f1"] for r in rs] == [0.0] * 4
    assert all(r["within_family_orthology"]["ortholog"]["recall"] == 0.0 for r in rs)
    assert all(not e["recovered_0.5"] for r in rs for e in r["events"])
    rb = _all(fams, st, controls.one_big_hog(fams, st))
    assert [r["pair_merge"]["accuracy_pairs"] for r in rb] == pytest.approx([0.9, 11 / 15, 1.0, 2 / 6])
    assert [r["pair_merge"]["accuracy_node_balanced"] for r in rb] == pytest.approx([0.75, 0.8, 1.0, 2 / 3])
    assert [r["integrity"]["membership_f1"] for r in rb] == pytest.approx([1, 1, 1, 2 / 3])
    assert rb[3]["integrity"]["fusion_pct"] == 100.0
    agg = aggregate(rb, n_boot=200)
    assert agg["combined"]["merge_accuracy_pairs"]["mean"] == pytest.approx((0.9 + 11 / 15 + 1 + 2 / 6) / 4)
    assert agg["combined"]["merge_accuracy_pairs"]["ci_low"] is None   # n = 4 < 5: no CI (B9)


def test_random_relabel(st, fam_by_id):
    """Panel = F2 alone; positive control relabelled within species. HUMAN and MOUSE each have
    2 genes (x, y) -> 4 equally likely permutations:
      identity                  -> 15/15 correct
      swap human only           -> HxMx, HyMy, HxMy, HyMx wrong -> 11/15
      swap mouse only           -> 11/15
      swap both (copies swapped consistently) -> 15/15
    so each relabelling scores 1 or 11/15 and the expectation is (2 + 22/15)/4 = 13/15."""
    f = fam_by_id["F2"]
    pos = controls.reference_prediction([f], st)
    ctx = panel_context([f], st)
    accs = [evaluate_family(f, st, controls.random_relabel(pos, seed=s), ctx)["pair_merge"]["accuracy_pairs"]
            for s in range(400)]
    assert set(np.round(accs, 9)) <= {1.0, round(11 / 15, 9)}
    assert np.mean(accs) == pytest.approx(13 / 15, abs=0.03)
    # deterministic for a fixed seed; HOG sizes preserved
    a, b = controls.random_relabel(pos, SEEDS[0]), controls.random_relabel(pos, SEEDS[0])
    assert a.levels == b.levels and a.sizes == pos.sizes


def test_swap_leaves_invariants(st, fams):
    pert = controls.swap_leaves(fams, k=5, seed=SEEDS[0])
    before = sorted(g for f in fams for g in f.member_ids)
    after = sorted(g for f in pert for g in f.member_ids)
    assert before == after
    for f in pert:
        assert sorted(f.tree.leaf_names()) == sorted(f.member_ids)
        assert validate_family(f, st).ok          # same-species swaps keep trees consistent
    # originals untouched
    assert fams[0].tree.leaf_names() == ["HUMAN_1", "MOUSE_1", "CHICK_1", "DANRE_a", "ORYLA_b"]
    assert controls.swap_leaves(fams, 0, SEEDS[0])[1].tree.leaf_names() == fams[1].tree.leaf_names()


def test_dose_response_monotone(st, fams):
    """Mean (over the fixed seed list) macro node-balanced merge accuracy decreases with the number of
    random same-species leaf swaps; k = 0 is the positive control."""
    ctx = panel_context(fams, st)
    means = []
    for k in (0, 1, 2, 4, 8, 16):
        vals = []
        for s in SEEDS:
            p = controls.perturbed_reference_prediction(fams, st, k, s)
            vals.append(np.mean([evaluate_family(f, st, p, ctx)["pair_merge"]["accuracy_node_balanced"] for f in fams]))
        means.append(np.mean(vals))
    assert means[0] == 1.0
    assert all(a > b for a, b in zip(means, means[1:])), means
