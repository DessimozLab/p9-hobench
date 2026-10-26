"""Metrics on the toy families (see conftest.py) with hand-derived expected values.

Depths: Euteleostomi 0, Amniota 1, Teleostei 1, Mammalia 2, CHICK/DANRE/ORYLA 2, HUMAN/MOUSE 3.
Signed error = depth(nearest bound) - depth(pred); + = predicted too old (over-split),
- = predicted too young (over-merged). Reference merges are intervals [lo, hi].
"""
import numpy as np
import pytest

from hobench import controls
from hobench.metrics import (INF, NOT_ASSERTED, ScoringError, evaluate_family, family_scalars, hog_codes,
                             interval_correct, interval_error, panel_context, predicted_merge_matrix,
                             reference_matrices)
from hobench.predictions import Prediction
from hobench.reference import reconcile

from conftest import make_family


def _eval(f, st, pred, fams):
    return evaluate_family(f, st, pred, panel_context(fams, st))


def _lv(r):
    return {d["level"]: d for d in r["per_level"]}


def test_interval_correct_and_error(st):
    i = st.idx
    lo = np.array([i("Amniota")] * 5 + [INF])
    hi = np.array([i("Euteleostomi")] * 4 + [INF, INF])
    p = np.array([i("Amniota"), i("Euteleostomi"), i("Mammalia"), INF, INF, i("Amniota")])
    # in [Amniota, Euteleostomi]: Amniota yes, Euteleostomi yes, Mammalia no (-1),
    # INF no (+inf, hi finite); [Amniota, INF] contains INF; exact INF vs Amniota -> -inf
    assert interval_correct(st, p, lo, hi).tolist() == [True, True, False, False, True, False]
    assert interval_error(st, p, lo, hi).tolist() == [0, 0, -1, np.inf, 0, -np.inf]
    # above a finite hi: [Mammalia, Mammalia] predicted at Euteleostomi -> 2 - 0 = +2
    assert interval_error(st, np.array([i("Euteleostomi")]), np.array([i("Mammalia")]),
                          np.array([i("Mammalia")])).tolist() == [2]


def test_reference_matrices_F2_F4(st, fam_by_id):
    rf = reconcile(fam_by_id["F2"], st)
    R = reference_matrices(rf)
    g = rf.gene_index
    nm = lambda a, b: (st.names[R.lo[g[a], g[b]]], st.names[R.hi[g[a], g[b]]])
    assert nm("HUMAN_x", "MOUSE_x") == ("Mammalia", "Mammalia")      # S: S0
    assert nm("HUMAN_x", "HUMAN_y") == ("Amniota", "Amniota")        # D at Mammalia (not root): exact
    assert nm("CHICK_g", "DANRE_g") == ("Euteleostomi", "Euteleostomi")
    rf4 = reconcile(fam_by_id["F4"], st)
    R4 = reference_matrices(rf4)
    a, b = rf4.gene_index["HUMAN_r1"], rf4.gene_index["CHICK_r2"]
    assert st.names[R4.lo[a, b]] == "Euteleostomi" and R4.hi[a, b] == INF     # root D: [parent(X), inf]
    rf3 = reconcile(fam_by_id["F3"], st)
    R3 = reference_matrices(rf3)
    assert R3.lo[rf3.gene_index["HUMAN_u1"], rf3.gene_index["MOUSE_u"]] == NOT_ASSERTED


def test_pseudo_orthology_one_big_hog(st, fams, fam_by_id):
    """F1 under the species-tree-shaped prediction (everything merges at its species LCA).

    10 pairs; only DANRE_a-ORYLA_b is wrong: ref exact parent(Teleostei) = Euteleostomi (depth 0),
    predicted Teleostei (depth 1) -> -1. Pair-pooled 9/10. Node-balanced: nodes root S (6 pairs,
    1.0), Amniota S (2, 1.0), Mammalia S (1, 1.0), D (1, 0.0) -> 3/4.
    Per level (pairs with S0 below L): Mammalia 1 TP; Amniota 3 TP; Euteleostomi 10 TP;
    Teleostei: the fish pair is FP -> pair F1 = 2*0/(0+1+0) = 0 (B1: not dropped). B-cubed at
    Teleostei: each fish gene P = (0+1)/(1+1) = 1/2, R = (0+1)/(0+1) = 1 -> F = 2/3; other levels 1.
    mean B-cubed = (1+1+1+2/3)/4 = 11/12; mean pair F1 = 3/4; pooled TP 14, FP 1 -> 28/29.
    Within-family orthology (all 10 cross-species): TP 9, FP 1 -> P 0.9, R 1.
    """
    r = _eval(fam_by_id["F1"], st, controls.one_big_hog(fams, st), fams)
    pm = r["pair_merge"]
    assert pm["accuracy_pairs"] == pytest.approx(0.9) and pm["accuracy_node_balanced"] == pytest.approx(0.75)
    assert pm["error_histogram"] == {"-1": 1, "0": 9}
    assert pm["frac_over_merged"] == pytest.approx(0.1) and pm["frac_over_split"] == 0
    assert pm["accuracy_exact_pairs"] == pytest.approx(0.9) and pm["accuracy_interval_pairs"] is None
    lv = _lv(r)
    assert lv["Teleostei"]["f1"] == 0.0 and lv["Teleostei"]["fp"] == 1
    assert lv["Teleostei"]["bcubed_f1"] == pytest.approx(2 / 3)
    assert r["per_level_summary"]["mean_bcubed_f1"] == pytest.approx(11 / 12)
    assert r["per_level_summary"]["mean_pair_f1"] == pytest.approx(0.75)
    assert r["per_level_summary"]["pooled"]["f1"] == pytest.approx(28 / 29)
    o = r["within_family_orthology"]
    assert o["ortholog"]["precision"] == pytest.approx(0.9) and o["ortholog"]["recall"] == 1.0
    (e,) = r["events"]
    assert e["frac_correct"] == 0.0 and e["separation"] == 0.0 and not e["recovered_0.5"]
    assert e["modal_pred_merge"] == "Teleostei" and e["misplacement"] == -1
    # Integrity: one RefOG (all 5 together at Euteleostomi); the big HOG also holds 14 genes of
    # the other toy families, but they are undeclared (homologous?) -> unlabelled, not FP (B6).
    it = r["integrity"]
    assert it["n_refogs"] == 1 and it["membership_f1"] == 1.0 and it["contamination"] == 0.0
    assert it["n_other_reference_family_genes"] == 14 and it["exact_pct"] == 100.0


def test_declared_decoys_and_complete_species(st, fam_by_id, fams):
    """B5/B6 on F1 + one-big HOG at Euteleostomi (19 panel genes, 5 of F1).
    distinct_from F2 at ALL levels: the 6 F2 genes are contaminants: contamination 6/(5+6);
    membership precision = 5/11 (overlap-weighted), recall 1 -> F1 = 2(5/11)/(16/11) = 10/16.
    distinct_from F2 at Amniota: the HOG is at Euteleostomi ⊄ Amniota -> no contamination.
    complete_species [HUMAN]: HUMAN genes in the HOG = 7 (1 of F1) -> 6 FP for precision,
    contamination unchanged (only declared decoys count)."""
    nwk = "(((HUMAN_1,MOUSE_1)[&&NHX:ev=S],CHICK_1)[&&NHX:ev=S],(DANRE_a,ORYLA_b)[&&NHX:ev=D:lvl=Teleostei])[&&NHX:ev=S];"
    others = [f for f in fams if f.id != "F1"]
    for extra, contam, f1, fus in [
        ({"distinct_from": [{"family": "F2", "at_level": "ALL", "src": "S1"}]}, 6 / 11, 10 / 16, 100.0),
        ({"distinct_from": [{"family": "F2", "at_level": "Amniota", "src": "S1"}]}, 0.0, 1.0, 0.0),
        ({"complete_species": ["HUMAN"]}, 0.0, 10 / 16, 0.0),
        ({"decoys": ["HUMAN_x", "HUMAN_y"]}, 2 / 7, 2 * (5 / 7) / (1 + 5 / 7), 100.0),
    ]:
        f1fam = make_family("F1", nwk, **extra)
        panel = [f1fam] + others
        r = evaluate_family(f1fam, st, controls.one_big_hog(panel, st), panel_context(panel, st))
        it = r["integrity"]
        assert it["contamination"] == pytest.approx(contam), extra
        assert it["membership_f1"] == pytest.approx(f1), extra
        assert it["fusion_pct"] == fus, extra


def test_lineage_specific_duplication_one_big_hog(st, fams, fam_by_id):
    """F2 (6 genes, 15 pairs). Wrong: the 4 cross-copy mammal pairs (ref Amniota, pred Mammalia)
    -> -1; pooled 11/15. Node-balanced: root S 1, Amniota S 1, D 0, two Mammalia S 1 -> 4/5.
    Mammalia level: 6 pairs, ref-together {HxMx, HyMy}, all predicted together: pair F1 = 4/8;
    B-cubed each mammal gene P = (1+1)/(3+1) = 1/2, R = 1 -> F = 2/3. Amniota, Euteleostomi: 1.
    mean B-cubed = (1+1+2/3)/3 = 8/9. Event: separation 0, cohesion 1 (each copy's mammal pair
    together at Mammalia), not recovered. Gene level: mammal genes 3/5, CHICK_g and DANRE_g 5/5."""
    r = _eval(fam_by_id["F2"], st, controls.one_big_hog(fams, st), fams)
    assert r["pair_merge"]["accuracy_pairs"] == pytest.approx(11 / 15)
    assert r["pair_merge"]["accuracy_node_balanced"] == pytest.approx(0.8)
    assert r["pair_merge"]["error_histogram"] == {"-1": 4, "0": 11}
    lv = _lv(r)
    assert set(lv) == {"Mammalia", "Amniota", "Euteleostomi"}
    assert (lv["Mammalia"]["tp"], lv["Mammalia"]["fp"], lv["Mammalia"]["fn"]) == (2, 4, 0)
    assert lv["Mammalia"]["f1"] == pytest.approx(0.5) and lv["Mammalia"]["bcubed_f1"] == pytest.approx(2 / 3)
    assert r["per_level_summary"]["mean_bcubed_f1"] == pytest.approx(8 / 9)
    (e,) = r["events"]
    assert e["separation"] == 0.0 and e["cohesion_min"] == 1.0 and not e["recovered_0.5"]
    assert e["modal_pred_merge"] == "Mammalia" and e["misplacement"] == -1
    genes = {d["gene"]: d for d in r["genes"]}
    assert genes["HUMAN_x"]["frac_correct"] == pytest.approx(0.6) and genes["CHICK_g"]["frac_correct"] == 1.0
    assert r["genes_summary"]["mean_frac_correct"] == pytest.approx((4 * 0.6 + 2) / 6)
    # tolerance ±1 (B12): the four -1 errors become correct; depth fraction = 1/family depth 3
    assert r["pair_merge"]["accuracy_tol1_pairs"] == 1.0
    assert r["family_depth"] == 3 and r["pair_merge"]["mean_abs_depth_frac_error"] == pytest.approx(4 / 15 / 3)


def test_polytomy_pairs_excluded(st, fams, fam_by_id):
    """F3: of the 6 pairs only the 3 with CHICK_u are asserted (merge Amniota). Two predictions
    that differ only inside the U polytomy score identically."""
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
        assert r["pair_merge"]["accuracy_node_balanced"] == 1.0
        assert r["within_family_orthology"]["ortholog"]["f1"] == 1.0
    assert ra["per_level_summary"] == rb["per_level_summary"]
    rs = _eval(f, st, controls.all_singletons(fams), fams)
    assert rs["pair_merge"]["error_histogram"] == {"+inf": 3}


def test_root_duplication_interval_AC1(st):
    """Regression for AC-1: TGD toy ((DANRE_a,ORYLA_a)S,(DANRE_b,ORYLA_b)S)D[lvl=Teleostei] with no
    outgroup member. Cross-copy pairs: [parent(Teleostei)=Euteleostomi, inf]. The correct
    prediction (copies apart at Teleostei, together at Euteleostomi) and the never-merge prediction
    both score 1.0 (old code: 0.33 vs 1.0). Merging the copies at Teleostei scores the 2
    orthologous pairs only: node-balanced (0 + 1 + 1)/3."""
    f = make_family("TGD", "((DANRE_a,ORYLA_a)[&&NHX:ev=S],(DANRE_b,ORYLA_b)[&&NHX:ev=S])[&&NHX:ev=D:lvl=Teleostei];")
    good, never, merged = Prediction(), Prediction(), Prediction()
    for p in (good, never):
        p.add("Teleostei", "a", "DANRE_a"); p.add("Teleostei", "a", "ORYLA_a")
        p.add("Teleostei", "b", "DANRE_b"); p.add("Teleostei", "b", "ORYLA_b")
    for g in ("DANRE_a", "ORYLA_a", "DANRE_b", "ORYLA_b"):
        good.add("Euteleostomi", "e", g)
        merged.add("Teleostei", "t", g)
    ctx = panel_context([f], st)
    for p in (good, never):
        r = evaluate_family(f, st, p, ctx)
        assert r["pair_merge"]["accuracy_pairs"] == 1.0 and r["n_interval_pairs"] == 4
        assert r["events"][0]["recovered_0.5"] and r["events"][0]["recovered_0.8"]
    r = evaluate_family(f, st, merged, ctx)
    assert r["pair_merge"]["accuracy_node_balanced"] == pytest.approx(2 / 3)
    assert r["pair_merge"]["error_histogram"] == {"-1": 4, "0": 2}


def test_all_singletons_cannot_win_root_duplication_ST3(st):
    """Regression for ST-3: KIN-like toy, 3 single-copy subfamilies x 4 species joined by a root
    duplication. All-singletons gets every cross-subfamily pair right ([.., inf]) but fails each
    subfamily's 3 S nodes: node-balanced accuracy = 1/4 (root D right, 3 S nodes wrong... each
    subfamily has 3 S nodes) -> 1/(1+9) = 0.1; event not recovered (cohesion 0)."""
    sub = lambda k: (f"(((HUMAN_{k},MOUSE_{k})[&&NHX:ev=S],CHICK_{k})[&&NHX:ev=S],DANRE_{k})[&&NHX:ev=S]")
    f = make_family("KIN", "(" + ",".join(sub(k) for k in "abc") + ")[&&NHX:ev=D];")
    r = evaluate_family(f, st, controls.all_singletons([f]), panel_context([f], st))
    assert r["pair_merge"]["accuracy_pairs"] == pytest.approx(48 / 66)      # 3*16 cross pairs of 66
    assert r["pair_merge"]["accuracy_node_balanced"] == pytest.approx(0.1)
    assert not r["events"][0]["recovered_0.5"] and r["events"][0]["cohesion_min"] == 0.0
    assert r["integrity"]["membership_f1"] == 0.0     # B1: scored 0, not dropped


def test_hgt_transfer_scoring(st):
    """A4: DANRE_x acquired from the CHICK lineage; vs native teleost genes merge exactly at
    lca(CHICK, DANRE) = Euteleostomi; xenolog pairs excluded from within-family orthology."""
    nwk = ("((DANRE_x[&&NHX:recipient=1],(DANRE_n,ORYLA_n)[&&NHX:ev=S])[&&NHX:ev=H:donor=CHICK:tag=HGT],"
           "CHICK_c)[&&NHX:ev=S];")
    f = make_family("HGT", nwk)
    good, naive = Prediction(), Prediction()
    for p in (good, naive):
        for g in ("DANRE_x", "DANRE_n", "ORYLA_n", "CHICK_c"):
            p.add("Euteleostomi", "e", g)
    good.add("Teleostei", "n", "DANRE_n"); good.add("Teleostei", "n", "ORYLA_n"); good.add("Teleostei", "x", "DANRE_x")
    for g in ("DANRE_x", "DANRE_n", "ORYLA_n"):
        naive.add("Teleostei", "t", g)          # xenolog merged with natives at Teleostei
    ctx = panel_context([f], st)
    rg, rn = evaluate_family(f, st, good, ctx), evaluate_family(f, st, naive, ctx)
    assert rg["pair_merge"]["accuracy_pairs"] == 1.0 and rg["events"][0]["recovered_0.5"]
    # naive: DANRE_x-DANRE_n (same species -> merge Teleostei) and DANRE_x-ORYLA_n (Teleostei) wrong (-1)
    assert rn["pair_merge"]["error_histogram"] == {"-1": 2, "0": 4}
    # 3 xenolog pairs (DANRE_x vs the 3 others) are excluded from the orthology classification
    assert rg["within_family_orthology"]["n_xenolog_pairs_excluded"] == 3
    assert rg["within_family_orthology"]["n_pairs"] == 3   # Dn-On, Dn-Cc, On-Cc


def test_yeast_allopolyploid_interval_A2():
    """A2: ohnologs apart at PostWGD, merge anywhere in [ZTWGD, Saccharomycetaceae] is correct;
    merging at PostWGD (pseudo-orthology) is -1. Per-level F1 skips ZTWGD for ohnolog pairs."""
    from hobench.speciestree import SpeciesTree
    yst = SpeciesTree.from_newick("(((YEAST,CANGA)PostWGD,ZYGRO)ZTWGD,(KLULA,ASHGO)KLE)Saccharomycetaceae;")
    f = make_family("YWGD", "(((YEAST_1,CANGA_1)[&&NHX:ev=S],(YEAST_2,CANGA_2)[&&NHX:ev=S])"
                            "[&&NHX:ev=D:lvl_lo=PostWGD:lvl_hi=ZTWGD],(KLULA_1,ASHGO_1)[&&NHX:ev=S])[&&NHX:ev=S];")
    ctx = panel_context([f], yst)
    copies = {"YEAST_1": "a", "CANGA_1": "a", "YEAST_2": "b", "CANGA_2": "b"}
    for merge_at in ("ZTWGD", "Saccharomycetaceae"):
        p = Prediction()
        for g, c in copies.items():
            p.add("PostWGD", c, g)
        if merge_at == "ZTWGD":
            for g in copies:
                p.add("ZTWGD", "z", g)
        else:
            for g, c in copies.items():
                p.add("ZTWGD", c, g)
        p.add("KLE", "k", "KLULA_1"); p.add("KLE", "k", "ASHGO_1")
        for g in list(copies) + ["KLULA_1", "ASHGO_1"]:
            p.add("Saccharomycetaceae", "r", g)
        r = evaluate_family(f, yst, p, ctx)
        assert r["pair_merge"]["accuracy_pairs"] == 1.0, merge_at
        assert _lv(r)["ZTWGD"]["n_skipped_interval"] == 4    # the 4 ohnolog pairs, inside the interval
    bad = Prediction()
    for g in copies:
        bad.add("PostWGD", "x", g)
    r = evaluate_family(f, yst, bad, ctx)
    # cross-ohnolog pairs: YEAST_1-YEAST_2 merge parent(YEAST)=PostWGD (depth 2) vs lo ZTWGD (1): -1
    assert r["pair_merge"]["error_histogram"]["-1"] == 4


def test_predicted_merge_uses_lowest_level(st, fam_by_id):
    rf = reconcile(fam_by_id["F2"], st)
    p = Prediction()
    p.add("Mammalia", "m", "HUMAN_x"); p.add("Mammalia", "m", "MOUSE_x")
    p.add("Amniota", "a1", "HUMAN_x"); p.add("Amniota", "a2", "MOUSE_x")
    p.add("Euteleostomi", "e", "HUMAN_x"); p.add("Euteleostomi", "e", "MOUSE_x")
    M = predicted_merge_matrix(rf, hog_codes(rf, p))
    g = rf.gene_index
    assert st.names[M[g["HUMAN_x"], g["MOUSE_x"]]] == "Mammalia"
    assert M[g["HUMAN_x"], g["CHICK_g"]] == INF
    assert p.nesting_problems(st) and not p.fatal_problems(st)


def test_no_silent_fallbacks(st, fams, fam_by_id):
    """Scoring raises on a gene placed at a level its species is not below, and on a species
    mismatch between prediction and reference."""
    p = Prediction()
    p.add("Teleostei", "t", "HUMAN_x"); p.add("Teleostei", "t", "DANRE_g")
    with pytest.raises(ScoringError, match="not below"):
        _eval(fam_by_id["F2"], st, p, fams)
    q = Prediction(gene_species={"HUMAN_x": "MOUSE"})
    with pytest.raises(ScoringError, match="different species"):
        _eval(fam_by_id["F2"], st, q, fams)


def test_spurious_duplications(st, fams, fam_by_id):
    """F1's pure-S clade (HUMAN_1, MOUSE_1, CHICK_1): levels Mammalia (2 genes) and Amniota (3)
    -> opportunities (2-1)+(3-1) = 3 (Euteleostomi is above the clade). All-singletons splits
    everything: 1 + 2 = 3 spurious -> rate 1; the one-big HOG: 0."""
    s = _eval(fam_by_id["F1"], st, controls.all_singletons(fams), fams)["spurious_duplications"]
    assert (s["n_clades"], s["count"], s["opportunities"], s["rate"]) == (1, 3, 3, 1.0)
    b = _eval(fam_by_id["F1"], st, controls.one_big_hog(fams, st), fams)["spurious_duplications"]
    assert b["count"] == 0


def test_orthobench_measures_against_reference_implementation(st):
    """B7: membership/OrthoBench measures for one RefOG of 4 genes split 2+1 with 1 missing.
    Reference formulas (Open_Orthobench benchmark.py, q_even; OF3 supplement S1):
      overlaps t = (2, 1), missing 1, n = 4.
      TP = C(2,2)+C(1,2) = 1; FN = [(4-2)*2 + (4-1)*1 + 1*(4-1)]/2 = 5; FP = 0
      -> normalised by n-1 = 3: P = 1, R = 1/(1+5) = 1/6.
      OF3 recall RC = (2*2/4 + 1*1/4)/(2+1) = 5/12; precision 1; fission: 2 PredOGs -> 100%.
      Entropy: p = (1/2, 1/4, 1/4 missing) -> H = 1.5 bits / log2(4) = 0.75.
    """
    f = make_family("OB", "(((HUMAN_1,MOUSE_1)[&&NHX:ev=S],CHICK_1)[&&NHX:ev=S],DANRE_1)[&&NHX:ev=S];")
    p = Prediction()
    p.add("Euteleostomi", "A", "HUMAN_1"); p.add("Euteleostomi", "A", "MOUSE_1")
    p.add("Amniota", "B", "CHICK_1"); p.add("Amniota", "B", "OTHER_CHICK")  # size-2 PredOG, unlabelled partner
    it = evaluate_family(f, st, p, panel_context([f], st))["integrity"]
    assert it["orthobench_precision"] == 1.0 and it["orthobench_recall"] == pytest.approx(1 / 6)
    assert it["membership_recall"] == pytest.approx(5 / 12) and it["membership_precision"] == 1.0
    assert it["fission_pct"] == 100.0 and it["missing_genes_pct"] == 25.0
    assert it["entropy"] == pytest.approx(0.75) and it["exact_pct"] == 0.0
    assert it["n_unlabelled"] == 1
    # cross-check with a direct transcription of benchmark.py's pairwise counts
    ref, preds = {"HUMAN_1", "MOUSE_1", "CHICK_1", "DANRE_1"}, [{"HUMAN_1", "MOUSE_1"}, {"CHICK_1", "OTHER_CHICK"}]
    tp = fp = fn = 0.0
    not_present = set(ref)
    for P in preds:
        ov = len(ref & P)
        if ov:
            not_present -= P
            tp += ov * (ov - 1) / 2
            fn += (len(ref) - ov) * ov
    fn = (fn + len(not_present) * (len(ref) - 1)) / 2
    assert tp / (tp + fn) == pytest.approx(it["orthobench_recall"])


def test_large_family_is_fast(st):
    """A 2000-member family (~2M pairs) evaluates quickly."""
    import time
    from hobench.newick import Node, write_newick
    root = Node(features={"ev": "D"})
    for c in range(400):
        clade = Node(features={"ev": "S"})
        amn = clade.add_child(Node(features={"ev": "S"}))
        mam = amn.add_child(Node(features={"ev": "S"}))
        mam.add_child(Node(f"HUMAN_{c}")); mam.add_child(Node(f"MOUSE_{c}"))
        amn.add_child(Node(f"CHICK_{c}"))
        tel = clade.add_child(Node(features={"ev": "S"}))
        tel.add_child(Node(f"DANRE_{c}")); tel.add_child(Node(f"ORYLA_{c}"))
        root.add_child(clade)
    f = make_family("BIG", write_newick(root, features=["ev"]))
    t0 = time.time()
    pos = controls.reference_prediction([f], st)
    r = evaluate_family(f, st, pos, panel_context([f], st), per_gene=False)
    t1 = time.time()
    assert r["n_asserted_pairs"] == 2000 * 1999 // 2
    assert r["pair_merge"]["accuracy_pairs"] == 1.0
    r2 = evaluate_family(f, st, controls.one_big_hog([f], st), panel_context([f], st), per_gene=False)
    t2 = time.time()
    # the root D maps to the panel root: cross-copy merge = [parent(root), inf] = exactly inf, so
    # one-big-HOG gets only the 400*10 within-copy pairs right
    assert r2["pair_merge"]["accuracy_pairs"] == pytest.approx(4000 / (2000 * 1999 / 2))
    assert t1 - t0 < 60 and t2 - t1 < 60, (t1 - t0, t2 - t1)


def test_u_root_over_interval_duplication_positive_control(st):
    """Regression (found on YGOB units): root U over a D whose lvl_lo lies above the U's LCA
    mapping. U((HUMAN_a,HUMAN_b)D[lvl_lo=Mammalia,lvl_hi=Amniota], MOUSE_c): the D has no
    constraining ancestor -> merge [Amniota, inf]; the family level is raised to the D's lvl_lo
    (Mammalia), where the reference fixes no grouping across the U: 3 singleton RefOGs. The
    positive control (everyone alone at Mammalia) is perfect, including integrity."""
    f = make_family("UR", "((HUMAN_a,HUMAN_b)[&&NHX:ev=D:lvl_lo=Mammalia:lvl_hi=Amniota],MOUSE_c)[&&NHX:ev=U];")
    rf = reconcile(f, st)
    assert st.names[rf.family_level] == "Mammalia"
    pos = controls.reference_prediction([f], st)
    r = evaluate_family(f, st, pos, panel_context([f], st))
    sc = family_scalars(r)
    assert r["integrity"]["n_refogs"] == 3
    assert sc["merge_accuracy_node"] == 1.0 and sc["membership_f1"] == 1.0 and sc["exact_pct"] == 100.0
