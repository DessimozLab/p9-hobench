"""Aggregation, paired inference and agreement (SYNTHESIS B4, B9, B10)."""
import math
import time

import numpy as np
import pytest

from hobench import controls
from hobench.evaluate import evaluate_panel
from hobench.stats import (aggregate, agreement, ami, bootstrap_indices, compare_methods, holm, items,
                           mean_ci)

from conftest import make_family


def fake(fid, value, panel="P", collection=None, tags=(), events=()):
    """Minimal family result carrying one scalar (node-balanced merge accuracy)."""
    return {"id": fid, "panel": panel, "collection": collection, "phenomena": list(tags),
            "pair_merge": {"accuracy_node_balanced": value}, "events": list(events)}


def ev(rec, tags=("WGD",)):
    return {"tags": list(tags), "recovered_0.5": rec, "recovered_0.8": rec, "frac_correct": float(rec)}


def test_ami_hand_values():
    """[0,0,1,1] vs [0,1,0,1]: MI = 0. EMI: each of the 4 cells has n_ij ~ Hypergeom with
    P(0)=1/6, P(1)=4/6, P(2)=1/6; term = n_ij/4 ln(n_ij) -> E = (1/6)(1/2)ln2 per cell, ln2/3 total.
    AMI = (0 - ln2/3)/(ln2 - ln2/3) = -1/2."""
    assert ami([0, 0, 1, 1], [0, 0, 1, 1]) == 1.0
    assert ami([0, 0, 1, 1], ["b", "b", "a", "a"]) == 1.0
    assert ami([0, 0, 1, 1], [0, 1, 0, 1]) == pytest.approx(-0.5)
    try:
        from sklearn.metrics import adjusted_mutual_info_score
    except ImportError:
        return
    rng = np.random.default_rng(1)
    for _ in range(5):
        a, b = rng.integers(0, 4, 30), rng.integers(0, 5, 30)
        assert ami(a, b) == pytest.approx(adjusted_mutual_info_score(a, b), abs=1e-9)


def test_bootstrap_indices_stratified_and_deterministic():
    ids = [f"f{i}" for i in range(10)]
    strata = {x: ("A" if int(x[1:]) < 4 else "B") for x in ids}
    s1, idx1 = bootstrap_indices(ids, 50, 7, strata)
    s2, idx2 = bootstrap_indices(list(reversed(ids)), 50, 7, strata)
    assert s1 == s2 and (idx1 == idx2).all()          # same IDs -> same resamples (paired)
    a_pos = {k for k, x in enumerate(s1) if strata[x] == "A"}
    for row in idx1:
        assert sum(1 for k in row if k in a_pos) == 4   # stratum sizes kept


def test_mean_ci_rules():
    small = mean_ci({"a": 0.2, "b": 0.6, "c": None})
    assert small == {"mean": 0.4, "ci_low": None, "ci_high": None, "n": 2, "min": 0.2, "max": 0.6}
    big = mean_ci({str(i): float(i) for i in range(10)}, n_boot=500)
    assert big["n"] == 10 and 0 <= big["ci_low"] <= 4.5 <= big["ci_high"] <= 9


def test_two_stage_collections_and_event_weighting():
    """3 ordinary families (1, 0, 0.5) + collection C with 4 units (1, 1, 1, 0 -> 0.75).
    families track = 0.5; collection C = 0.75 (n = 4 units, no CI); combined (C = one item) =
    (1 + 0 + 0.5 + 0.75)/4 = 0.5625. WGD events: family F1 has 1 recovered event, C has 9 not
    recovered -> per-item means 1 and 0 -> 0.5 (pooling would give 0.1)."""
    res = [fake("F1", 1.0, events=[ev(True)]), fake("F2", 0.0), fake("F3", 0.5)]
    res += [fake(f"U{i}", v, collection="C", events=[ev(False)] * (3 if i < 3 else 0))
            for i, v in enumerate([1, 1, 1, 0])]
    agg = aggregate(res, n_boot=100)
    m = "merge_accuracy_node"
    assert agg["families"][m]["mean"] == pytest.approx(0.5)
    assert agg["collections"]["C"]["n_units"] == 4 and agg["collections"]["C"][m]["mean"] == pytest.approx(0.75)
    assert agg["combined"][m]["mean"] == pytest.approx(0.5625) and agg["combined"][m]["n"] == 4
    assert agg["events_by_tag"]["WGD"]["recovery_rate_0.5"] == pytest.approx(0.5)
    assert agg["events_by_tag"]["WGD"]["n_events"] == 10 and agg["events_by_tag"]["WGD"]["n_items"] == 2
    assert set(items(res)) == {"F1", "F2", "F3", "collection:C"}


def test_paired_comparison_and_holm():
    a = [fake(f"F{i}", 0.5 + 0.1 * (i % 3), panel="P" if i < 4 else "Q") for i in range(8)]
    b = [fake(f"F{i}", 0.4 + 0.1 * (i % 3), panel="P" if i < 4 else "Q") for i in range(8)]
    c = compare_methods(a, b, "merge_accuracy_node", n_boot=200)
    # every paired difference is exactly 0.1 -> CI collapses to 0.1
    assert c["n"] == 8 and c["mean_diff"] == pytest.approx(0.1)
    assert c["ci_low"] == pytest.approx(0.1) and c["ci_high"] == pytest.approx(0.1)
    assert c["wilcoxon_p"] < 0.05
    same = compare_methods(a, a, "merge_accuracy_node", n_boot=50)
    assert same["mean_diff"] == 0 and same["wilcoxon_p"] == 1.0
    few = compare_methods(a[:3], b[:3], "merge_accuracy_node")
    assert few["n"] == 3 and few["ci_low"] is None and few["min_diff"] == pytest.approx(0.1)
    # Holm: sorted p (0.01, 0.03, 0.04) * (3, 2, 1) = (0.03, 0.06, 0.04) -> monotone (0.03, 0.06, 0.06)
    h = holm({"a": 0.01, "b": 0.04, "c": 0.03, "d": None})
    assert h["a"] == pytest.approx(0.03) and h["c"] == pytest.approx(0.06) and h["b"] == pytest.approx(0.06)
    assert h["d"] is None


def test_agreement_positive_vs_relabelled(st, fams):
    pos = controls.reference_prediction(fams, st)
    assert agreement(fams[1], st, pos, pos)["mean"] == 1.0
    big = controls.one_big_hog(fams, st)
    ag = agreement(fams[1], st, pos, big)
    assert ag["per_level"]["Mammalia"] == pytest.approx(ami([0, 0, 1, 1], [0, 0, 0, 0]))


def test_collection_of_5000_units_is_cheap(st):
    """B10: a genome-scale collection of 5,000 small units is scored unit by unit (compact
    results), never as one big family."""
    fams = [make_family(f"U{i}", f"((HUMAN_{i},MOUSE_{i})[&&NHX:ev=S],CHICK_{i})[&&NHX:ev=S];", collection="YGOB")
            for i in range(5000)]
    t0 = time.time()
    pred = controls.reference_prediction(fams, st, check=False)
    res = evaluate_panel(st, fams, pred, n_boot=100, with_skill=False)
    dt = time.time() - t0
    assert res["aggregate"]["collections"]["YGOB"]["n_units"] == 5000
    assert res["aggregate"]["combined"]["merge_accuracy_node"]["n"] == 1
    assert res["aggregate"]["collections"]["YGOB"]["merge_accuracy_node"]["mean"] == 1.0
    assert "genes" not in res["families"][0] and "nodes" not in res["families"][0]
    assert dt < 180, dt
