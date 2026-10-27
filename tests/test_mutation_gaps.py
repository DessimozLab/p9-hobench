"""Tests added after scripts/mutation_check.py (see design/MUTATION_REPORT.md): each kills a mutant
that the original suite caught only incidentally (not by a test about the mutated behaviour)."""
import numpy as np
import pytest

from hobench.metrics import evaluate_family, panel_context, reference_matrices
from hobench.predictions import Prediction
from hobench.reference import INF, reconcile
from conftest import make_family


def test_nested_duplication_under_root_duplication_has_infinite_upper_bound(st):
    """Mutant c: ABOVE_RANGE (an unbounded Dhi inherited from an enclosing root duplication) treated as
    the finite root level instead of inf.
    Family: ((DANRE_a1,DANRE_a2)D, (DANRE_b,ORYLA_b)S)D[lvl_lo=Teleostei,lvl_hi=ABOVE_RANGE]; the inner D has no
    level: lo = its mapping DANRE, hi = bound = enclosing D's Dhi = ABOVE_RANGE (unbounded -> inf).  Inner pair interval = [parent(DANRE)=Teleostei,
    inf].  Prediction: only DANRE_a1 + ORYLA_b in one HOG at Teleostei; everything else unplaced, so
    the pair (DANRE_a1, DANRE_a2) is never merged = inf, inside [Teleostei, inf] -> correct.
    If ABOVE were the finite root, inf would lie above hi -> incorrect."""
    f = make_family("NEST", "((DANRE_a1,DANRE_a2)[&&NHX:ev=D],"
                            "(DANRE_b,ORYLA_b)[&&NHX:ev=S])[&&NHX:ev=D:lvl_lo=Teleostei:lvl_hi=ABOVE_RANGE];")
    rf = reconcile(f, st)
    R = reference_matrices(rf)
    i, j = rf.genes.index("DANRE_a1"), rf.genes.index("DANRE_a2")
    assert R.lo[i, j] == st.idx("Teleostei") and R.hi[i, j] == INF
    p = Prediction()
    p.add("Teleostei", "x", "ORYLA_b"); p.add("Teleostei", "x", "DANRE_a1")
    from hobench.metrics import predicted_merge_matrix, hog_codes, interval_correct
    Mp = predicted_merge_matrix(rf, hog_codes(rf, p))
    assert Mp[i, j] == INF
    assert interval_correct(st, np.array([Mp[i, j]]), R.lo[[i], [j]], R.hi[[i], [j]])[0]
