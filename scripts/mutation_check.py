#!/usr/bin/env python
"""Mutation check for hobench: apply each mutation to a temporary copy of hobench/ (+tests/),
run `pytest -q tests` there, and record whether the suite FAILS (killed, good) or passes
(survived = test gap).  Usage: python scripts/mutation_check.py [--tests tests/test_x.py ...]
[--json out.json].  The real hobench/ is never modified."""
import argparse, json, shutil, subprocess, sys, tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent

# id -> (description, file, [(old, new, count)])  ; each patch asserts `old` occurs `count` times
MUTANTS = {
 "a_sign_flip": ("flip sign of signed level error", "metrics.py", [
    ("err[below_lo] = (d[_safe(lo)] - d[_safe(p)])[below_lo]", "err[below_lo] = (d[_safe(p)] - d[_safe(lo)])[below_lo]", 1),
    ("err[above_fin] = (d[_safe(hi)] - d[_safe(p)])[above_fin]", "err[above_fin] = (d[_safe(p)] - d[_safe(hi)])[above_fin]", 1)]),
 "b_swap_S_D": ("swap S and D handling in reference_matrices", "metrics.py", [
    ('        if ev.ev == "S":\n            val_lo = val_hi = lambda A, B: st.lca_idx', '        if ev.ev == "D":\n            val_lo = val_hi = lambda A, B: st.lca_idx', 1),
    ('        elif ev.ev in ("D", "H"):\n            val_lo = (lambda v', '        elif ev.ev in ("S", "H"):\n            val_lo = (lambda v', 1)]),
 "c_no_inf": ("drop INF/ABOVE_RANGE: treat as finite root", "reference.py", [
    ("    if lvl == ABOVE:\n        return INF\n    p = int(st.parent_idx[lvl])\n    return p if p >= 0 else INF",
     "    if lvl == ABOVE:\n        return 0\n    p = int(st.parent_idx[lvl])\n    return p if p >= 0 else 0", 1)]),
 "d_interval_eq": ("interval check lo<=x<=hi -> x==lo", "metrics.py", [
    ("return np.where(lf, ~below_lo & ~above_hi, ~pf)", "return np.where(lf, pf & (p == lo), ~pf)", 1)]),
 "e_no_testable": ("event testability rule removed", "metrics.py", [
    ('d["testable"] = bool(bounded or any(x is not None for x in cohesions))', 'd["testable"] = True', 1)]),
 "f_f1_drop": ("F1 dropped when precision or recall denominator empty", "metrics.py", [
    ("    return 2 * tp / d if d else None", "    return 2 * tp / d if (tp + fp and tp + fn) else None", 1)]),
 "g_dup_late": ("OrthoXML dup_placement default declared -> late", "predictions.py", [
    ('dup_placement: str = "declared",', 'dup_placement: str = "late",', 1)]),
 "h_no_flat": ("flat TSV projection disabled", "predictions.py", [
    ("            if flat:\n                if not species", "            if False:\n                if not species", 1)]),
 "i_fp_outside_complete": ("membership FP counts other-family genes outside complete cells", "metrics.py", [
    ('return (d["n_members"] - d["own"][r_]) + d["contam"] + d["complete_nonmembers"]',
     'return (d["n_members"] - d["own"][r_]) + d["contam"] + d["complete_nonmembers"] + d["other_family"]', 1)]),
 "j_unpaired_bootstrap": ("paired bootstrap broken (method b items misaligned)", "stats.py", [
    ('diffs = {k: ia[k]["scalars"][metric] - ib[k]["scalars"][metric] for k in common}',
     'diffs = {k: ia[k]["scalars"][metric] - ib[common[-1 - n]]["scalars"][metric] for n, k in enumerate(common)}', 1)]),
}


def apply(src_dir: Path, fname: str, patches):
    p = src_dir / fname
    s = p.read_text()
    for old, new, cnt in patches:
        assert s.count(old) == cnt, f"{fname}: patch target found {s.count(old)}x (expected {cnt}): {old[:60]!r}"
        s = s.replace(old, new)
    p.write_text(s)


def run(mid, tests):
    desc, fname, patches = MUTANTS[mid]
    with tempfile.TemporaryDirectory(prefix=f"mut_{mid}_") as td:
        td = Path(td)
        shutil.copytree(ROOT / "hobench", td / "hobench", ignore=shutil.ignore_patterns("__pycache__"))
        shutil.copytree(ROOT / "tests", td / "tests", ignore=shutil.ignore_patterns("__pycache__"))
        for e in ROOT.iterdir():   # everything else (wrapper/, config/, families/, ...) read-only via symlink
            if e.name not in ("hobench", "tests", ".git") and not (td / e.name).exists():
                (td / e.name).symlink_to(e)
        apply(td / "hobench", fname, patches)
        r = subprocess.run([sys.executable, "-m", "pytest", "-q", "-p", "no:cacheprovider", *tests],
                           cwd=td, capture_output=True, text=True)
        failed = [l for l in r.stdout.splitlines() if l.startswith("FAILED")]
        if len(failed) > 1: failed = [f"{len(failed)} tests fail; first: " + failed[0][7:90]]
        assert r.returncode in (0, 1), f"harness error (rc={r.returncode}): {r.stdout[-400:]}"
        return r.returncode == 1, (failed[0][:110] if failed else r.stdout.strip().splitlines()[-1])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tests", nargs="*", default=["tests"])
    ap.add_argument("--only", nargs="*")
    ap.add_argument("--json")
    a = ap.parse_args()
    res = {}
    for mid in a.only or MUTANTS:
        killed, note = run(mid, a.tests)
        res[mid] = {"desc": MUTANTS[mid][0], "killed": killed, "note": note}
        print(f"{mid:24s} {'KILLED' if killed else 'SURVIVED'}  {note}", flush=True)
    if a.json:
        Path(a.json).write_text(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
