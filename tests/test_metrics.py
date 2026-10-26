"""Metrics on the toy families (see conftest.py) with hand-derived expected values.

Depths: Euteleostomi 0, Amniota 1, Teleostei 1, Mammalia 2, CHICK/DANRE/ORYLA 2, HUMAN/MOUSE 3.
Signed error = depth(ref merge) - depth(pred merge); + = predicted too old (over-split),
- = predicted too young (over-merged).
"""
import numpy as np
import pytest

from hobench import controls
from hobench.metrics import (INF, NOT_ASSERTED, aggregate, bootstrap_ci, evaluate_family,
                             family_scalars, panel_context, predicted_merge_matrix, hog_codes,
                             reference_merge_matrix)
from hobench.predictions import Prediction
from hobench.reference import reconcile

from conftest import make_family


def _eval(f, st, pred, fams):
    return evaluate_family(f, st, pred, panel_context(fams))


def test_reference_merge_matrix_F2(st, fam_by_id):
    rf = reconcile(fam_by_id["F2"], st)
    M = reference_merge_matrix(rf)
    g = rf.gene_index
    name = lambda a, b: st.names[M[g[a], g[b]]]
    assert name("HUMAN_x", "MOUSE_x") == "Mammalia"      # S node: S0 = lca(HUMAN, MOUSE)
    assert name("HUMAN_x", "HUMAN_y") == "Amniota"       # D at Mammalia -> parent = Amniota
    assert name("HUMAN_x", "MOUSE_y") == "Amniota"
    assert name("MOUSE_y", "CHICK_g") == "Amniota"       # S, S0 = Amniota
    assert name("CHICK_g", "DANRE_g") == "Euteleostomi"
    rf3 = reconcile(fam_by_id["F3"], st)
    M3 = reference_merge_matrix(rf3)
    g3 = rf3.gene_index
    assert M3[g3["HUMAN_u1"], g3["MOUSE_u"]] == NOT_ASSERTED
    assert st.names[M3[g3["HUMAN_u1"], g3["CHICK_u"]]] == "Amniota"
    rf4 = reconcile(fam_by_id["F4"], st)
    M4 = reference_merge_matrix(rf4)
    assert M4[rf4.gene_index["HUMAN_r1"], rf4.gene_index["CHICK_r2"]] == INF


def test_pseudo_orthology_one_big_hog(st, fams, fam_by_id):
    """F1 under a species-tree-shaped prediction (everything merges at its species LCA).

    10 pairs. All amniote/fish pairs merge at Euteleostomi, amniote pairs at their LCA: correct.
    DANRE_a-ORYLA_b: ref merge = parent(Teleostei) = Euteleostomi (depth 0), predicted at
    S0 = Teleostei (depth 1) -> error 0 - 1 = -1 (pseudo-orthologs called orthologs).
    accuracy = 9/10. Orthology (all 10 pairs cross-species): ref orthologs 9, predicted 10,
    TP 9, FP 1 -> P = 0.9, R = 1, F1 = 1.8/1.9. Paralogs: ref 1, predicted 0 -> recall 0.
    """
    r = _eval(fam_by_id["F1"], st, controls.one_big_hog(fams, st), fams)
    pm = r["pair_merge"]
    assert r["n_asserted_pairs"] == 10
    assert pm["accuracy"] == pytest.approx(0.9)
    assert pm["error_histogram"] == {"-1": 1, "0": 9}
    assert pm["frac_over_merged"] == pytest.approx(0.1) and pm["frac_over_split"] == 0
    o = r["orthology"]
    assert (o["n_ref_orthologs"], o["n_ref_paralogs"]) == (9, 1)
    assert o["ortholog"]["precision"] == pytest.approx(0.9) and o["ortholog"]["recall"] == 1.0
    assert o["ortholog"]["f1"] == pytest.approx(1.8 / 1.9)
    assert o["paralog"]["recall"] == 0.0 and o["paralog"]["precision"] is None
    # Event: 1 cross-copy pair, predicted at Teleostei instead of Euteleostomi -> not recovered,
    # misplaced by -1 level (too young).
    (e,) = r["events"]
    assert e["level"] == "Teleostei" and e["ref_merge"] == "Euteleostomi" and e["tags"] == ["WGD"]
    assert e["frac_correct"] == 0.0 and not e["recovered"]
    assert e["modal_pred_merge"] == "Teleostei" and e["misplacement"] == -1
    # Per level: at Teleostei the only pair (DANRE_a, ORYLA_b) is predicted together but should
    # be apart -> FP=1, precision 0, recall undefined.
    lv = {d["level"]: d for d in r["per_level"]}
    assert lv["Teleostei"]["fp"] == 1 and lv["Teleostei"]["precision"] == 0.0 and lv["Teleostei"]["recall"] is None
    assert lv["Euteleostomi"]["f1"] == 1.0
    # Contamination: the predicted HOG at Euteleostomi holds all 19 panel genes, 5 are F1's.
    integ = r["integrity"]
    assert integ["n_hogs_cover"] == 1 and integ["ref_n_hogs_cover"] == 1 and integ["split_excess"] == 0
    assert integ["contamination"] == pytest.approx(14 / 19) and integ["coverage"] == 1.0


def test_lineage_specific_duplication_one_big_hog(st, fams, fam_by_id):
    """F2 (6 genes, 15 pairs) under the species-tree-shaped prediction.

    Wrong pairs: the 4 cross-copy mammal pairs (HxHy, HxMy, HyMx, MxMy): ref Amniota (depth 1),
    predicted Mammalia (depth 2; same-species pairs merge at parent(HUMAN) = Mammalia) -> -1.
    accuracy = 11/15.
    Per level: Mammalia has 6 pairs, ref-together = {HxMx, HyMy}, predicted all 6 together:
    P = 2/6, R = 1, F1 = 0.5. Amniota: 10 pairs all together in both -> 1. Euteleostomi: 15 -> 1.
    Teleostei has no pair. mean F1 = (0.5+1+1)/3; pooled TP = 2+10+15 = 27, FP = 4 -> F1 = 54/58.
    Orthology over 13 cross-species pairs: ref orthologs 11, predicted 13 -> P = 11/13, R = 1.
    Gene level: each mammal gene has 5 pairs, 3 correct (0.6); CHICK_g and DANRE_g 5/5.
    """
    r = _eval(fam_by_id["F2"], st, controls.one_big_hog(fams, st), fams)
    assert r["pair_merge"]["accuracy"] == pytest.approx(11 / 15)
    assert r["pair_merge"]["error_histogram"] == {"-1": 4, "0": 11}
    lv = {d["level"]: d for d in r["per_level"]}
    assert set(lv) == {"Mammalia", "Amniota", "Euteleostomi"}
    assert (lv["Mammalia"]["tp"], lv["Mammalia"]["fp"], lv["Mammalia"]["fn"]) == (2, 4, 0)
    assert lv["Mammalia"]["f1"] == pytest.approx(0.5)
    assert r["per_level_summary"]["mean_f1"] == pytest.approx(2.5 / 3)
    assert r["per_level_summary"]["pooled"]["f1"] == pytest.approx(54 / 58)
    o = r["orthology"]
    assert o["n_pairs"] == 13 and o["n_ref_orthologs"] == 11
    assert o["ortholog"]["precision"] == pytest.approx(11 / 13) and o["ortholog"]["recall"] == 1.0
    (e,) = r["events"]
    assert e["n_pairs"] == 4 and e["frac_correct"] == 0.0 and e["modal_pred_merge"] == "Mammalia"
    assert e["misplacement"] == -1
    genes = {d["gene"]: d for d in r["genes"]}
    assert genes["HUMAN_x"]["frac_correct"] == pytest.approx(0.6) and genes["CHICK_g"]["frac_correct"] == 1.0
    assert r["species"]["HUMAN"] == {"n_genes": 2, "mean_frac_correct": pytest.approx(0.6),
                                     "n_perfect": 0, "n_majority": 2}
    assert r["genes_summary"]["mean_frac_correct"] == pytest.approx((4 * 0.6 + 2) / 6)
    assert r["genes_summary"]["frac_perfect"] == pytest.approx(2 / 6)


def test_polytomy_pairs_excluded(st, fams, fam_by_id):
    """F3: of the 6 pairs only the 3 with CHICK_u are asserted (merge Amniota). Two predictions
    that differ only inside the U polytomy must score identically."""
    f = fam_by_id["F3"]
    a, b = Prediction(), Prediction()
    for p in (a, b):
        for g in ("HUMAN_u1", "HUMAN_u2", "MOUSE_u", "CHICK_u"):
            p.add("Amniota", "h", g)
            p.add("Euteleostomi", "e", g)
    a.add("Mammalia", "m1", "HUMAN_u1"); a.add("Mammalia", "m1", "MOUSE_u"); a.add("Mammalia", "m2", "HUMAN_u2")
    for g in ("HUMAN_u1", "HUMAN_u2", "MOUSE_u"):
        b.add("Mammalia", f"s_{g}", g)
    ra, rb = _eval(f, st, a, fams), _eval(f, st, b, fams)
    assert ra["n_pairs_total"] == 6 and ra["n_asserted_pairs"] == 3
    for r in (ra, rb):
        assert r["pair_merge"]["accuracy"] == 1.0
        assert r["orthology"]["ortholog"]["f1"] == 1.0
    assert {d["level"] for d in ra["per_level"]} == {"Amniota"}  # Mammalia: no asserted pair
    # with no prediction, the 3 asserted pairs are all '+inf' errors
    rs = _eval(f, st, controls.all_singletons(fams), fams)
    assert rs["pair_merge"]["error_histogram"] == {"+inf": 3}


def test_root_duplication_infinite_merge(st, fams, fam_by_id):
    """F4: Hr1-Cr1 and Hr2-Cr2 merge at Amniota; the 4 cross-copy pairs never merge (INF).
    all-singletons gets the 4 INF pairs right: 4/6; one-big-HOG gets only the 2 orthologs: 2/6,
    the other 4 are '-inf' (merged although they never should). The reference itself has 2
    HOGs at the family level (Amniota), so the expected split is 2."""
    f = fam_by_id["F4"]
    rs = _eval(f, st, controls.all_singletons(fams), fams)
    assert rs["pair_merge"]["accuracy"] == pytest.approx(4 / 6)
    assert rs["pair_merge"]["error_histogram"] == {"0": 4, "+inf": 2}
    assert rs["events"][0]["recovered"] and rs["events"][0]["misplacement"] == 0
    assert rs["integrity"]["n_hogs_cover"] == 4 and rs["integrity"]["ref_n_hogs_cover"] == 2
    assert rs["integrity"]["split_excess"] == 2 and rs["integrity"]["coverage"] == 0.0
    rb = _eval(f, st, controls.one_big_hog(fams, st), fams)
    assert rb["pair_merge"]["accuracy"] == pytest.approx(2 / 6)
    assert rb["pair_merge"]["error_histogram"] == {"-inf": 4, "0": 2}
    assert rb["events"][0]["misplacement"] == "-inf"
    # one big HOG at Amniota = 16 amniote genes of the panel, 4 of them F4 -> 12/16
    assert rb["integrity"]["split_excess"] == -1 and rb["integrity"]["contamination"] == pytest.approx(0.75)


def test_over_split_sign(st, fams, fam_by_id):
    """F2 where the prediction keeps CHICK_g out of the amniote HOG: CHICK_g pairs (4 mammal
    genes) are predicted to merge at Euteleostomi (depth 0) instead of Amniota (depth 1):
    error 1 - 0 = +1 (too old, over-split)."""
    ref = controls.reference_prediction(fams, st)
    p = Prediction(levels={lv: dict(m) for lv, m in ref.levels.items()}, sizes=ref.sizes,
                   gene_species=ref.gene_species)
    p.levels["Amniota"]["CHICK_g"] = "elsewhere"
    r = _eval(fam_by_id["F2"], st, p, fams)
    assert r["pair_merge"]["error_histogram"] == {"0": 11, "+1": 4}
    assert r["pair_merge"]["frac_over_split"] == pytest.approx(4 / 15)


def test_predicted_merge_uses_lowest_level(st, fam_by_id):
    """Non-nested input: two genes together at Euteleostomi and at Mammalia but not at
    Amniota -> merge_pred is the lowest level (Mammalia)."""
    rf = reconcile(fam_by_id["F2"], st)
    p = Prediction()
    p.add("Mammalia", "m", "HUMAN_x"); p.add("Mammalia", "m", "MOUSE_x")
    p.add("Amniota", "a1", "HUMAN_x"); p.add("Amniota", "a2", "MOUSE_x")
    p.add("Euteleostomi", "e", "HUMAN_x"); p.add("Euteleostomi", "e", "MOUSE_x")
    M = predicted_merge_matrix(rf, hog_codes(rf, p))
    g = rf.gene_index
    assert st.names[M[g["HUMAN_x"], g["MOUSE_x"]]] == "Mammalia"
    assert M[g["HUMAN_x"], g["CHICK_g"]] == INF
    assert p.check(st) and "not nested" in p.check(st)[0]


def test_aggregate_and_bootstrap(st, fams):
    pred = controls.one_big_hog(fams, st)
    ctx = panel_context(fams)
    res = [evaluate_family(f, st, pred, ctx) for f in fams]
    agg = aggregate(res, n_boot=500)
    # macro mean of per-family accuracies 0.9, 11/15, 1, 2/6
    assert agg["macro"]["merge_accuracy"]["mean"] == pytest.approx((0.9 + 11 / 15 + 1 + 2 / 6) / 4)
    lo, hi = agg["macro"]["merge_accuracy"]["ci_low"], agg["macro"]["merge_accuracy"]["ci_high"]
    assert 2 / 6 <= lo <= hi <= 1.0
    assert agg["per_tag"]["WGD"]["merge_accuracy"]["mean"] == pytest.approx(0.9)
    assert agg["events_by_tag"]["WGD"] == {"n_events": 1, "recovery_rate": 0.0, "mean_frac_correct": 0.0}
    assert bootstrap_ci([0.5, 0.5, 0.5]) == (0.5, 0.5)
    assert bootstrap_ci([0.1, 0.9], seed=1) == bootstrap_ci([0.1, 0.9], seed=1)  # deterministic
    assert family_scalars(res[0])["merge_accuracy"] == pytest.approx(0.9)


def test_large_family_is_fast(st):
    """A 2000-member family (~2M pairs) evaluates quickly thanks to block/vectorised code."""
    import time
    from hobench.newick import Node
    sp = ["HUMAN", "MOUSE", "CHICK", "DANRE", "ORYLA"]
    root = Node(features={"ev": "D"})        # 400 copies duplicated at the root level
    for c in range(400):
        clade = Node(features={"ev": "S"})
        amn = clade.add_child(Node(features={"ev": "S"}))
        mam = amn.add_child(Node(features={"ev": "S"}))
        mam.add_child(Node(f"HUMAN_{c}")); mam.add_child(Node(f"MOUSE_{c}"))
        amn.add_child(Node(f"CHICK_{c}"))
        tel = clade.add_child(Node(features={"ev": "S"}))
        tel.add_child(Node(f"DANRE_{c}")); tel.add_child(Node(f"ORYLA_{c}"))
        root.add_child(clade)
    from hobench.newick import write_newick
    f = make_family("BIG", write_newick(root, features=["ev"]))
    assert len(f.member_ids) == 2000
    t0 = time.time()
    pos = controls.reference_prediction([f], st)
    r = evaluate_family(f, st, pos, panel_context([f]), per_gene=False)
    t1 = time.time()
    assert r["n_asserted_pairs"] == 2000 * 1999 // 2
    assert r["pair_merge"]["accuracy"] == 1.0
    assert r["events"][0]["n_pairs"] == 2000 * 1999 // 2 - 400 * 10
    r2 = evaluate_family(f, st, controls.one_big_hog([f], st), panel_context([f]), per_gene=False)
    t2 = time.time()
    # one-big-HOG: only the 400*10 within-copy pairs are right
    assert r2["pair_merge"]["accuracy"] == pytest.approx(4000 / (2000 * 1999 / 2))
    assert t1 - t0 < 60 and t2 - t1 < 60, (t1 - t0, t2 - t1)
