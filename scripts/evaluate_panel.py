#!/usr/bin/env python
"""Score every run under runs/<P>/<mode>/ plus baselines/controls against all reference
families of panel P.  Usage: python scripts/evaluate_panel.py --panel Y --mode full

Writes results/<P>_<mode>/{per_family,per_level,events,summary}.tsv, summary.md, provenance.json.
Methods: each run dir (tool/variant from its manifest.json; FastOMA twice: dup placement
'declared' = primary, 'late' = sensitivity) and controls (positive_control, all_singletons,
one_big_hog, random_relabel = mean over seeds.SEEDS[:10] of the relabelled positive control).
Families without family.yaml (abandoned units) are not loaded.  Non-headline families
(unique mapping < 80 %, `headline: false`, or a curation_log status line saying NOT headline)
are scored and flagged in column `headline`.  Metrics/references are never altered here.
"""
from __future__ import annotations

import argparse, copy, csv, json, re, sys, tarfile, time
from collections import defaultdict
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from hobench import controls                                              # noqa: E402
from hobench.evaluate import _git_state, _path_hash, N_RELABEL_SEEDS      # noqa: E402
from hobench.freeze import MANIFEST, check_manifest, manifest_entries, tree_hash  # noqa: E402
from hobench.ids import make_normalizer                                   # noqa: E402
from hobench.metrics import (PRIMARY, SCALARS, evaluate_family, family_scalars,   # noqa: E402
                             panel_context)
from hobench.panel import load_panel, sha256_file, validate_members       # noqa: E402
from hobench.predictions import read_prediction                           # noqa: E402
from hobench.reference import load_families, reconcile, validate_panel    # noqa: E402
from hobench.seeds import SEEDS                                           # noqa: E402
from hobench.speciestree import SpeciesTree                               # noqa: E402
from hobench.stats import aggregate, compare_methods, holm                # noqa: E402

PF_METRICS = PRIMARY + ["in_paralog_recovery_0.5", "merge_accuracy_sd_posthoc", "event_recovery_0.8", "merge_accuracy_pairs", "merge_accuracy_tol1", "spurious_dup_rate",
                        "membership_precision", "membership_recall", "orthobench_f", "frac_over_split",
                        "frac_over_merged", "level_pair_f1_pooled", "mean_abs_depth_frac_error", "wf_ortholog_f1"]
N_OF = {"level_bcubed_f1": "n_asserted_pairs", "merge_accuracy_node": "n_asserted_pairs", "membership_f1": "n_members"}
HEAVY = ("nodes", "per_level", "genes", "species", "per_gene")


# ---------------------------------------------------------------- headline flag
CATALOGUE = {}
with open(ROOT / "results" / "catalogue" / "families.tsv") as _fh:    # single source of truth for headline/track
    for _r in csv.DictReader(_fh, delimiter="\t"):
        CATALOGUE[_r["id"]] = _r
COLL_IDS = [k for k, r in CATALOGUE.items() if int(r["n_units"]) > 1]


def is_collection_unit(f) -> bool:
    """True collections (many units, one item): YGOB, RPROT.  Panel-level `collection` labels (E_cytoskeleton_motors,
    K04_isofunctional_decoys, ...) do not make their families one item."""
    return bool(f.collection) and any(f.collection in cid for cid in COLL_IDS)


def headline_info(f) -> tuple[bool, str]:
    """(headline, reason) from results/catalogue/families.tsv.  'n: no events' (K04 decoys) stays in scoring as
    headline (membership/contamination).  Collection units: per-unit unique mapping < 80 % (mapping_report.tsv)."""
    row = CATALOGUE.get(f.id)
    if row is not None:
        h = row["headline"]
        return (True, "") if h.startswith("y") or h.startswith("n: no events") else (False, h)
    if f.meta.get("headline") is False:
        return False, "headline:false"
    p = Path(f.path) if f.path else None
    if p is not None:
        for cand in (p / "mapping_report.tsv", p.parent.parent / "mapping_report.tsv"):
            if cand.exists():
                rows = list(csv.DictReader(open(cand), delimiter="\t"))
                if rows and "unit" in rows[0]:
                    rows = [r for r in rows if r["unit"] == f.id]
                rows = [r for r in rows if r.get("status") in ("mapped", "unmapped", "ambiguous", "multi", "failed")]
                if rows:
                    fr = sum(r["status"] == "mapped" for r in rows) / len(rows)
                    if fr < 0.8:
                        return False, f"mapped {fr:.0%}"
                break
    return True, ""


# ---------------------------------------------------------------- scoring
SCALARS["in_paralog_recovery_0.5"] = ("events_summary", "in_paralog_recovery_0.5")   # leaf-branch events (secondary, ERRATA E4)
SCALARS["merge_accuracy_sd_posthoc"] = ("pair_merge", "sd_balanced_posthoc")   # POST-RESULTS secondary (ERRATA)
CONTROL_COLLECTIONS = {"RPROT"}
CONTROL_IDS = {"K02_trGTPases", "E05_MCM"}


def track_of(f) -> str:
    """control = K01 RPROT units, K02_trGTPases, E05_MCM; collection = YGOB; family = rest."""
    if f.id in CATALOGUE:
        return CATALOGUE[f.id]["track"]
    if f.collection in CONTROL_COLLECTIONS or f.id in CONTROL_IDS or f.id.startswith("K01_RPROT"):
        return "control"
    return "collection" if f.collection == "YGOB" else "family"


def score(st, rfs, pred, ctx):   # items: r['collection'] kept only for true collection units
    """-> list of slim result dicts (scalars intact) + per_level/event rows."""
    out = []
    for rf in rfs:
        r = evaluate_family(rf, st, pred, ctx, per_gene=False, compact=False)
        r["_levels"] = r["per_level"]
        if not is_collection_unit(rf.family):
            r["collection"] = None
        # ERRATA E4: headline event recovery = internal-branch testable events only; leaf events -> in-paralog recovery
        es = r["events_summary"]
        for e in r["events"]:
            e["event_branch"] = "leaf" if st.is_species(e["level"]) else "internal"
        for br, suffix in (("internal", ""), ("leaf", None)):
            sel = [e for e in r["events"] if e["event_branch"] == br and e.get("testable", True)]
            for th in (0.5, 0.8):
                v = float(np.mean([e[f"recovered_{th}"] for e in sel])) if sel else None
                if suffix is None:
                    es[f"in_paralog_recovery_{th}"] = v
                else:
                    es[f"recovery_rate_{th}"] = v
            es["n_internal_testable" if suffix == "" else "n_leaf_testable"] = len(sel)
        acc = defaultdict(list)
        for nd in r["nodes"]:
            if nd.get("accuracy") is not None:
                acc[nd["ev"]].append(nd["accuracy"])
        r["pair_merge"]["sd_balanced_posthoc"] = (0.5 * (np.mean(acc["S"]) + np.mean(acc["D"]))
                                                  if acc["S"] and acc["D"] else None)
        if not (rf.family.complete_species or rf.family.decoys):   # R2-3: membership F1 undefined without complete cells / decoys
            for k in ("membership_f1", "membership_precision", "membership_recall"):
                r["integrity"][k] = None
        for k in HEAVY:
            r.pop(k, None)
        out.append(r)
    return out


def set_path(d, path, v):
    for k in path[:-1]:
        d = d.setdefault(k, {})
    d[path[-1]] = v


def mean_results(runs):
    base = copy.deepcopy(runs[0])
    for i, b in enumerate(base):
        scs = [family_scalars(r[i]) for r in runs]
        for name, path in SCALARS.items():
            xs = [s[name] for s in scs if s.get(name) is not None]
            set_path(b, path, float(np.mean(xs)) if xs else None)
    return base


def set_skill(results, base):
    for r in results:
        sc = family_scalars(r)
        r["skill"] = {}
        for m in PRIMARY:
            b = base[r["id"]][m]
            r["skill"][m] = None if sc[m] is None or b is None or b >= 1 else (sc[m] - b) / (1 - b)


def baselines(ctrl_single, ctrl_big, relabel_runs):
    out = {}
    for i, r in enumerate(ctrl_single):
        d = {}
        for m in PRIMARY:
            rv = [family_scalars(x[i])[m] for x in relabel_runs]
            rv = [x for x in rv if x is not None]
            c = [x for x in (family_scalars(r)[m], family_scalars(ctrl_big[i])[m], float(np.mean(rv)) if rv else None)
                 if x is not None]
            d[m] = max(c) if c else None
        out[r["id"]] = d
    return out


# ---------------------------------------------------------------- predictions
def load_run(run: Path, st, kw, cache: Path, placements):
    man = json.loads((run / "manifest.json").read_text())
    if man.get("status") != "done":
        raise SystemExit(f"REFUSING to score {run}: manifest status {man.get('status')!r} != 'done'")
    tool, variant = man["tool"], man.get("variant", "default")
    res = []
    if tool == "fastoma":
        for dp in placements:
            pred = read_prediction("orthoxml", run / "FastOMA_HOGs.orthoxml.gz", st, id_attr="protId",
                                   tool_species_tree=run / "species_tree_checked.nwk", dup_placement=dp, **kw)
            res.append((tool, f"{variant}_{dp}", pred, run / "FastOMA_HOGs.orthoxml.gz", man))
    elif tool == "of_flat":
        res.append((tool, variant, read_prediction("tsv", run / "hogs.tsv.gz", st, flat=True, **kw), run / "hogs.tsv.gz", man))
    elif tool == "orthofinder":
        d = cache / run.parent.parent.name / run.parent.name / run.name
        if not (d / "Phylogenetic_Hierarchical_Orthogroups").exists():
            d.mkdir(parents=True, exist_ok=True)
            with tarfile.open(run / "of_core.tar.gz") as t:
                t.extractall(d)
        pred = read_prediction("orthofinder", d / "Phylogenetic_Hierarchical_Orthogroups", st,
                               of_species_tree=run / "SpeciesTree_rooted_node_labels.txt", **kw)
        res.append((tool, variant, pred, run / "of_core.tar.gz", man))
    else:
        raise SystemExit(f"unknown tool {tool} in {run}")
    return res


def fnum(x):
    return "" if x is None else f"{x:.6g}"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--panel", required=True)
    ap.add_argument("--mode", default="full")
    ap.add_argument("--families", default=str(ROOT / "families"))
    ap.add_argument("--runs", default=None)
    ap.add_argument("--out", default=None)
    ap.add_argument("--cache", default="/tmp/hobench_eval_cache")
    ap.add_argument("--n-boot", type=int, default=2000)
    ap.add_argument("--all-metrics", action="store_true", help="per_family.tsv with every SCALAR (large)")
    a = ap.parse_args()
    P, mode = a.panel, a.mode
    rundir = Path(a.runs or ROOT / "runs" / P / mode)
    out = Path(a.out or ROOT / "results" / f"{P}_{mode}")
    out.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    log = lambda *x: print(f"[{time.time()-t0:6.0f}s]", *x, file=sys.stderr, flush=True)

    ptree = ROOT / "panels" / P / "species_tree.nwk"
    st = SpeciesTree.from_newick(str(ptree))
    fams = load_families(a.families, panel=P)
    excluded = []
    if mode == "subset":   # score only families whose members are all present in the subset proteomes
        sub = {r["id"] for r in csv.DictReader(open(ROOT / "panels" / P / "subset" / "members.tsv"), delimiter="\t")}
        excluded = [f.id for f in fams if not set(f.member_ids) <= sub]
        fams = [f for f in fams if f.id not in set(excluded)]
        log(f"subset mode: {len(fams)} families scored, {len(excluded)} excluded (members absent from subset)")
        (out / "excluded_families.txt").write_text("\n".join(excluded) + "\n")
    rep = validate_panel(fams, st)
    panel = load_panel(ROOT / "panels" / P)
    for e in panel.verify_checksums():
        rep.errors.append(e)
    for f in fams:
        rep.extend(validate_members(f, panel))
    for w in rep.warnings[:20]:
        log("WARNING", w)
    if not rep.ok:
        print(rep, file=sys.stderr)
        return 1
    hl = {f.id: headline_info(f) for f in fams}
    tk = {f.id: track_of(f) for f in fams}
    keep, cs = set(), set()
    for f in fams:
        keep |= set(f.member_ids) | set(f.decoys)
        cs |= set(f.complete_species)
    kw = {"normalize": make_normalizer("auto"), "keep_genes": keep, "count_species": cs}
    ctx = panel_context(fams, st)
    rfs = [reconcile(f, st, check=False) for f in fams]
    log(f"{len(fams)} families ({sum(not h[0] for h in hl.values())} non-headline)")

    methods = []   # (method, variant, kind, results, source)
    # controls
    pos = controls.reference_prediction(fams, st, check=False)
    pos_res = score(st, rfs, pos, ctx); log("positive control")
    sing = score(st, rfs, controls.all_singletons(fams), ctx); log("singletons")
    big = score(st, rfs, controls.one_big_hog(fams, st), ctx); log("one-HOG")
    rel_runs = []
    for s in SEEDS[:N_RELABEL_SEEDS]:
        rel_runs.append(score(st, rfs, controls.random_relabel(pos, s), ctx)); log("relabel", s)
    rel = mean_results(rel_runs)
    base = baselines(sing, big, rel_runs)
    methods += [("positive_control", "reference", "control", pos_res, None),
                ("all_singletons", "", "control", sing, None),
                ("one_big_hog", "", "control", big, None),
                ("random_relabel", f"mean{N_RELABEL_SEEDS}", "control", rel, None)]
    # tool runs
    prov_runs = {}
    for run in sorted(p for p in rundir.iterdir() if (p / "manifest.json").exists()):
        for tool, variant, pred, src, man in load_run(run, st, kw, Path(a.cache), ["declared", "late"]):
            fatal = pred.fatal_problems(st)
            if fatal:
                raise SystemExit(f"{run.name}: prediction inconsistent with panel tree: {fatal[:5]}")
            miss = sorted(ctx.count_species - pred.count_species)
            if miss:
                raise SystemExit(f"{run.name}: no per-HOG counts for complete species {miss}")
            res = score(st, rfs, pred, ctx)
            methods.append((tool, variant, "tool", res, run))
            meta = {k: v for k, v in pred.meta.items() if k != "warnings"}
            prov_runs[f"{tool}/{variant}"] = {
                "run_dir": str(run), "prediction_sha256": _path_hash(Path(src)),
                "manifest_sha256": sha256_file(run / "manifest.json"),
                "manifest_git_commit": man.get("git_commit"), "prediction_meta": meta,
                "n_warnings": len(pred.meta.get("warnings", [])), "nesting_problems": len(pred.nesting_problems(st)),
                "n_genes_dropped_not_in_reference": None}
            log("scored", tool, variant)
    for m in methods:
        set_skill(m[3], base)

    # ---------------------------------------------------------------- tables
    pf_keys = list(SCALARS) if a.all_metrics else PF_METRICS
    pfm = pf_keys + [f"skill_{m}" for m in PRIMARY]
    with open(out / "per_family.tsv", "w", newline="") as fh, \
            open(out / "per_level.tsv", "w", newline="") as fl, open(out / "events.tsv", "w", newline="") as fe:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        wl = csv.writer(fl, delimiter="\t", lineterminator="\n")
        we = csv.writer(fe, delimiter="\t", lineterminator="\n")
        w.writerow(["family", "collection", "track", "method", "variant", "metric", "value", "n", "headline", "headline_note"])
        wl.writerow(["family", "collection", "method", "variant", "level", "depth", "n_pairs", "n_genes", "tp", "fp", "fn",
                     "pair_f1", "bcubed_precision", "bcubed_recall", "bcubed_f1"])
        we.writerow(["family", "collection", "method", "variant", "event", "ev", "level", "tags", "n_copies", "n_pairs",
                     "frac_correct", "separation", "cohesion_min", "merge_side", "frac_over_split", "frac_over_merged",
                     "modal_pred_merge", "event_branch", "testable", "recovered_0.5", "recovered_0.8"])
        for meth, var, kind, res, _ in methods:
            for r in res:
                sc = family_scalars(r)
                coll, (hd, note) = r.get("collection") or "", hl[r["id"]]
                for m in pfm:
                    if sc.get(m) is None:
                        continue
                    nn = r.get(N_OF.get(m, ""), None)
                    if m.startswith("event") or m == "spurious_dup_rate":
                        nn = r["events_summary"]["n_internal_testable"] if m.startswith("event") else r["spurious_duplications"]["opportunities"]
                    elif nn is None:
                        nn = r["n_asserted_pairs"] if "merge" in m or "level" in m else r["n_members"]
                    w.writerow([r["id"], coll, tk[r["id"]], meth, var, m, fnum(sc[m]), nn, int(hd), note])
                for L in r["_levels"]:
                    wl.writerow([r["id"], coll, meth, var, L["level"], L["depth"], L["n_pairs"], L["n_genes"], L["tp"], L["fp"],
                                 L["fn"], fnum(L["f1"]), fnum(L["bcubed_precision"]), fnum(L["bcubed_recall"]),
                                 fnum(L["bcubed_f1"])])
                for e in r["events"]:
                    we.writerow([r["id"], coll, meth, var, e["event"], e["ev"], e["level"], ",".join(e["tags"]),
                                 e["n_copies"], e["n_pairs"], fnum(e["frac_correct"]), fnum(e["separation"]),
                                 fnum(e["cohesion_min"]), fnum(e["merge_side"]), fnum(e["frac_over_split"]),
                                 fnum(e["frac_over_merged"]), e["modal_pred_merge"], e["event_branch"], int(e.get("testable", True)), int(e["recovered_0.5"]),
                                 int(e["recovered_0.8"])])
    log("per_* written")

    # ---- summary (descriptive only: means over items, min/max, n; no CIs -- confirmatory test is scripts/confirmatory.py)
    from hobench.stats import items as stat_items
    smetrics = PRIMARY + [f"skill_{m}" for m in PRIMARY] + ["in_paralog_recovery_0.5", "merge_accuracy_sd_posthoc", "event_recovery_0.8"]
    scopes = {"headline": lambda r: tk[r["id"]] != "control" and hl[r["id"]][0],
              "noncontrol_all": lambda r: tk[r["id"]] != "control",
              "control": lambda r: tk[r["id"]] == "control"}

    def describe(res, scope):
        its = stat_items([r for r in res if scopes[scope](r)])
        out_ = {}
        for m in smetrics:
            xs = [v["scalars"][m] for v in its.values() if v["scalars"].get(m) is not None]
            out_[m] = (float(np.mean(xs)) if xs else None, min(xs) if xs else None, max(xs) if xs else None, len(xs))
        return out_

    rows, desc = [], {}
    for meth, var, kind, res, _ in methods:
        for scope in scopes:
            d = desc[(meth, var, scope)] = describe(res, scope)
            for m in smetrics:
                mean_, lo, hi, n_ = d[m]
                rows.append([scope, meth, var, kind, m, fnum(mean_), fnum(lo), fnum(hi), n_, "", ""])
    # relabel: min-max over the 10 seeds of the item-mean
    for scope in scopes:
        for m in smetrics:
            vals = []
            for rr in rel_runs:
                v = describe(rr, scope)[m][0]
                if v is not None:
                    vals.append(v)
            if vals:
                rows.append([scope, "random_relabel", f"seeds{N_RELABEL_SEEDS}", "control", m, fnum(np.mean(vals)), "", "", "",
                             fnum(min(vals)), fnum(max(vals))])
    # sanity gate: tool identical to all-singletons on >50% of units (all primary metrics equal)
    gate = []
    sing_sc = {r["id"]: family_scalars(r) for r in sing}
    for meth, var, kind, res, _ in methods:
        if kind != "tool":
            continue
        same = 0
        for r in res:
            sc, sg = family_scalars(r), sing_sc[r["id"]]
            same += all((sc[m] is None and sg[m] is None) or (sc[m] is not None and sg[m] is not None and abs(sc[m] - sg[m]) < 1e-9)
                        for m in PRIMARY)
        fr = same / len(res)
        gate.append((meth, var, fr, fr > 0.5))
        rows.append(["gate", meth, var, "tool", "frac_units_identical_to_all_singletons", fnum(fr), "", "", len(res), "", ""])
    with open(out / "summary.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["scope", "method", "variant", "kind", "metric", "mean", "min", "max", "n", "seed_min", "seed_max"])
        w.writerows(rows)

    # provenance
    fdir = Path(a.families)
    frozen = (fdir / MANIFEST).exists()
    prov = {"hobench_git": _git_state(ROOT), "script": "scripts/evaluate_panel.py", "panel": P, "mode": mode,
            "panel_tree_sha256": sha256_file(ptree), "panel_species_tsv_sha256": sha256_file(ROOT / "panels" / P / "species.tsv"),
            "families_tree_hash": tree_hash(manifest_entries(fdir)), "families_manifest_present": frozen,
            "families_manifest_differences": check_manifest(fdir) if frozen else None,
            "n_families": len(fams), "n_excluded_subset": len(excluded), "n_non_headline": sum(not h[0] for h in hl.values()), "n_control_track": sum(t == "control" for t in tk.values()),
            "seeds_relabel": SEEDS[:N_RELABEL_SEEDS], "n_boot": a.n_boot, "runs": prov_runs}
    (out / "provenance.json").write_text(json.dumps(prov, indent=1, default=str))

    # summary.md
    def cell(d, m):
        mean_, lo, hi, n_ = d[m]
        return "NA" if mean_ is None else f"{mean_:.3f}"
    md = [f"# Panel {P} / {mode}: evaluation summary", "",
          f"hobench git {prov['hobench_git']['commit'][:10]} (dirty={prov['hobench_git']['dirty']}); {len(fams)} families "
          f"(+{len(excluded)} excluded: members absent from subset) ({sum(t == 'control' for t in tk.values())} control-track, {prov['n_non_headline']} non-headline); "
          f"families tree hash {prov['families_tree_hash'][:12]}.",
          "Descriptive only (mean over items; a collection is one item; no CIs). Headline = tracks family+collection, "
          "headline-eligible; control track (K01 RPROT, K02, E05) reported separately. Trivial baselines are rows of every table; "
          "skill = (s-b)/(1-b), b = per-family max(singletons, one big HOG, mean relabel). membership F1 is NA for families "
          "without complete cells or decoys. Confirmatory test: scripts/confirmatory.py.", ""]
    for scope in ("headline", "control"):
        md += [f"## scope: {scope}", "",
               "| method | variant | " + " | ".join(f"{m} (n) | skill" for m in PRIMARY) + " | S/D-bal merge acc (post-hoc) |",
               "|---|---|" + "---|" * (2 * len(PRIMARY) + 1)]
        for meth, var, kind, res, _ in methods:
            d = desc[(meth, var, scope)]
            md.append(f"| {meth} | {var} | " + " | ".join(f"{cell(d, m)} ({d[m][3]}) | {cell(d, 'skill_' + m)}" for m in PRIMARY)
                      + f" | {cell(d, 'merge_accuracy_sd_posthoc')} |")
        rl = [r for r in rows if r[0] == scope and r[1] == "random_relabel" and r[2].startswith("seeds") and r[4] in PRIMARY]
        md += ["", "random relabel, min-max of item-mean over 10 seeds: " +
               "; ".join(f"{r[4]} {r[9]}-{r[10]}" for r in rl), ""]
    md += ["## Sanity gate (tool identical to all-singletons on >50% of units)", ""]
    md += [f"- {m} {v}: {fr:.1%} of units identical to all-singletons" + ("  **FLAGGED**" if f else "") for m, v, fr, f in gate]
    md += ["", "## Positive control / degenerate baselines (per-family min/mean/max)", ""]
    bad = []
    for meth, var, kind, res, _ in methods:
        if kind != "control":
            continue
        for m in PRIMARY:
            xs = [family_scalars(r)[m] for r in res if family_scalars(r)[m] is not None]
            if xs:
                md.append(f"- {meth} {m}: mean {np.mean(xs):.4f}, min {min(xs):.4f}, max {max(xs):.4f}, n={len(xs)}")
                if meth == "positive_control" and min(xs) < 1 - 1e-9:
                    bad.append((m, [r["id"] for r in res if (family_scalars(r)[m] or 1) < 1 - 1e-9][:10]))
    md += ["", "POSITIVE CONTROL: " + ("PERFECT on all primary endpoints." if not bad else f"NOT perfect: {bad}")]
    (out / "summary.md").write_text("\n".join(md) + "\n")
    log("done ->", out)
    print("\n".join(md[:30]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
