"""Positive / negative / dose-response controls."""
import numpy as np
import pytest

from hobench import controls
from hobench.metrics import aggregate, evaluate_family, family_scalars, panel_context
from hobench.reference import validate_family
from hobench.seeds import SEEDS


def _all(fams, st, pred):
    ctx = panel_context(fams)
    return [evaluate_family(f, st, pred, ctx) for f in fams]


def test_positive_control_is_perfect(st, fams):
    pred = controls.reference_prediction(fams, st)
    assert pred.meta["inconsistent_pairs"] == 0
    assert pred.check(st) == []
    for r in _all(fams, st, pred):
        sc = family_scalars(r)
        assert sc["merge_accuracy"] == 1.0, r["id"]
        assert sc["level_mean_f1"] == 1.0 and sc["level_pooled_f1"] == 1.0
        assert sc["ortholog_f1"] == 1.0
        assert sc["gene_mean_frac_correct"] == 1.0 and sc["gene_frac_perfect"] == 1.0
        assert sc["contamination"] == 0.0 and sc["split_excess"] == 0
        assert sc["event_recovery_rate"] in (1.0, None)
        assert all(e["misplacement"] == 0 for e in r["events"])
        assert set(r["pair_merge"]["error_histogram"]) == {"0"}


def test_negative_controls_known_values(st, fams):
    """Per-family merge accuracies (derivations in test_metrics.py):
    all-singletons: F1 0, F2 0, F3 0, F4 4/6 (only the ∞ pairs are right);
    one-big-HOG:    F1 9/10, F2 11/15, F3 1 (3/3 asserted), F4 2/6."""
    rs = _all(fams, st, controls.all_singletons(fams))
    assert [r["pair_merge"]["accuracy"] for r in rs] == pytest.approx([0, 0, 0, 4 / 6])
    assert all(r["orthology"]["ortholog"]["recall"] == 0.0 for r in rs)
    assert all(r["integrity"]["coverage"] == 0.0 for r in rs)
    rb = _all(fams, st, controls.one_big_hog(fams, st))
    assert [r["pair_merge"]["accuracy"] for r in rb] == pytest.approx([0.9, 11 / 15, 1.0, 2 / 6])
    assert all(r["orthology"]["ortholog"]["recall"] == 1.0 for r in rb)
    agg = aggregate(rb, n_boot=200)
    assert agg["macro"]["merge_accuracy"]["mean"] == pytest.approx((0.9 + 11 / 15 + 1 + 2 / 6) / 4)


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
    ctx = panel_context([f])
    accs = [evaluate_family(f, st, controls.random_relabel(pos, seed=s), ctx)["pair_merge"]["accuracy"]
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
    """Mean (over the fixed seed list) macro merge accuracy decreases with the number of
    random same-species leaf swaps; k = 0 is the positive control."""
    ctx = panel_context(fams)
    means = []
    for k in (0, 1, 2, 4, 8, 16):
        vals = []
        for s in SEEDS:
            p = controls.perturbed_reference_prediction(fams, st, k, s)
            vals.append(np.mean([evaluate_family(f, st, p, ctx)["pair_merge"]["accuracy"] for f in fams]))
        means.append(np.mean(vals))
    assert means[0] == 1.0
    assert all(a > b for a, b in zip(means, means[1:])), means
