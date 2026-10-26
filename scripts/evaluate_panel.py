#!/usr/bin/env python
"""Score every run under runs/<P>/<mode>/ plus baselines/controls against all reference
families of panel P.  Usage: python scripts/evaluate_panel.py --panel Y --mode full

Writes results/<P>_<mode>/{per_family,per_level,events,summary}.tsv, summary.md, provenance.json.
Methods: each run dir (tool/variant from its manifest.json; FastOMA twice: dup placement
'declared' = primary, 'late' = sensitivity) and controls (positive_control, all_singletons,
one_hog_per_family, random_relabel = mean over seeds.SEEDS[:10] of the relabelled positive control).
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

PF_METRICS = PRIMARY + ["event_recovery_0.8", "merge_accuracy_pairs", "merge_accuracy_tol1", "spurious_dup_rate",
                        "membership_precision", "membership_recall", "orthobench_f", "frac_over_split",
                        "frac_over_merged", "level_pair_f1_pooled", "mean_abs_depth_frac_error", "wf_ortholog_f1"]
N_OF = {"level_bcubed_f1": "n_asserted_pairs", "merge_accuracy_node": "n_asserted_pairs", "membership_f1": "n_members"}
HEAVY = ("nodes", "per_level", "genes", "species", "per_gene")


# ---------------------------------------------------------------- headline flag
def headline_info(f) -> tuple[bool, str]:
    """(headline, reason).  Non-headline: yaml headline:false; mapping_report unique-mapped < 80 %;
    curation_log status line (first 6 lines) containing 'not headline'."""
    if f.meta.get("headline") is False:
        return False, "headline:false"
    p = Path(f.path) if f.path else None
    if p is None:
        return True, ""
    for cand in (p / "mapping_report.tsv", p.parent.parent / "mapping_report.tsv"):
        if cand.exists():
            try:
                rows = list(csv.DictReader(open(cand), delimiter="\t"))
                if "unit" in rows[0]:
                    rows = [r for r in rows if r["unit"] == f.id]
                rows = [r for r in rows if r.get("status") in ("mapped", "unmapped", "ambiguous", "multi", "failed")]
                if rows:
                    fr = sum(r["status"] == "mapped" for r in rows) / len(rows)
                    if fr < 0.8:
                        return False, f"mapped {fr:.0%}"
            except Exception:
                pass
            break
    log = p / "curation_log.md"
    if log.exists():
        head = " ".join(log.read_text().splitlines()[:6]).lower()
        if re.search(r"\bnot headline", head) and "superseded" not in head.split("not headline")[0][-0:]:
            return False, "log status: not headline"
    return True, ""


# ---------------------------------------------------------------- scoring
def score(st, rfs, pred, ctx):
    """-> list of slim result dicts (scalars intact) + per_level/event rows."""
    out = []
    for rf in rfs:
        r = evaluate_family(rf, st, pred, ctx, per_gene=False, compact=False)
        r["_levels"] = r["per_level"]
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
    tool, variant = man["tool"], man.get("variant", "default")
    res = []
    if tool == "fastoma":
        for dp in placements:
            pred = read_prediction("orthoxml", run / "FastOMA_HOGs.orthoxml.gz", st, id_attr="protId",
                                   tool_species_tree=run / "species_tree_checked.nwk", dup_placement=dp, **kw)
            res.append((tool, f"{variant}_{dp}", pred, run / "FastOMA_HOGs.orthoxml.gz", man))
    elif tool == "of_flat":
        res.append((tool, variant, read_prediction("tsv", run / "hogs.tsv.gz", st, **kw), run / "hogs.tsv.gz", man))
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
                ("one_hog_per_family", "", "control", big, None),
                ("random_relabel", f"mean{N_RELABEL_SEEDS}", "control", rel, None)]
    del rel_runs
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
        w.writerow(["family", "collection", "method", "variant", "metric", "value", "n", "headline", "headline_note"])
        wl.writerow(["family", "collection", "method", "variant", "level", "depth", "n_pairs", "n_genes", "tp", "fp", "fn",
                     "pair_f1", "bcubed_precision", "bcubed_recall", "bcubed_f1"])
        we.writerow(["family", "collection", "method", "variant", "event", "ev", "level", "tags", "n_copies", "n_pairs",
                     "frac_correct", "separation", "cohesion_min", "merge_side", "frac_over_split", "frac_over_merged",
                     "modal_pred_merge", "recovered_0.5", "recovered_0.8"])
        for meth, var, kind, res, _ in methods:
            for r in res:
                sc = family_scalars(r)
                coll, (hd, note) = r.get("collection") or "", hl[r["id"]]
                for m in pfm:
                    if sc.get(m) is None:
                        continue
                    nn = r.get(N_OF.get(m, ""), None)
                    if m.startswith("event") or m == "spurious_dup_rate":
                        nn = r["events_summary"]["n_events"] if m.startswith("event") else r["spurious_duplications"]["opportunities"]
                    elif nn is None:
                        nn = r["n_asserted_pairs"] if "merge" in m or "level" in m else r["n_members"]
                    w.writerow([r["id"], coll, meth, var, m, fnum(sc[m]), nn, int(hd), note])
                for L in r["_levels"]:
                    wl.writerow([r["id"], coll, meth, var, L["level"], L["depth"], L["n_pairs"], L["n_genes"], L["tp"], L["fp"],
                                 L["fn"], fnum(L["f1"]), fnum(L["bcubed_precision"]), fnum(L["bcubed_recall"]),
                                 fnum(L["bcubed_f1"])])
                for e in r["events"]:
                    we.writerow([r["id"], coll, meth, var, e["event"], e["ev"], e["level"], ",".join(e["tags"]),
                                 e["n_copies"], e["n_pairs"], fnum(e["frac_correct"]), fnum(e["separation"]),
                                 fnum(e["cohesion_min"]), fnum(e["merge_side"]), fnum(e["frac_over_split"]),
                                 fnum(e["frac_over_merged"]), e["modal_pred_merge"], int(e["recovered_0.5"]),
                                 int(e["recovered_0.8"])])
    log("per_* written")

    # summary: aggregate per method on all families and on headline-only
    smetrics = PRIMARY + [f"skill_{m}" for m in PRIMARY] + ["event_recovery_0.8", "spurious_dup_rate"]
    rows = []
    aggs = {}
    for meth, var, kind, res, _ in methods:
        for scope, sel in (("all", res), ("headline", [r for r in res if hl[r["id"]][0]])):
            ag = aggregate(sel, n_boot=a.n_boot)
            aggs[(meth, var, scope)] = ag
            tracks = [("combined", ag["combined"]), ("families", ag["families"])] + \
                     [(f"collection:{c}", {k: v for k, v in d.items() if k != "n_units"}) for c, d in ag["collections"].items()]
            for tname, tr in tracks:
                for m in smetrics:
                    d = tr[m]
                    rows.append([scope, tname, meth, var, kind, m, fnum(d["mean"]), fnum(d["ci_low"]), fnum(d["ci_high"]),
                                 fnum(d.get("min")), fnum(d.get("max")), d["n"]])
    # paired comparisons OF vs FastOMA (declared), Holm over the 4 primary
    def get(tool, var):
        for m in methods:
            if m[0] == tool and m[1] == var:
                return m[3]
    cmp_rows = []
    of_d, fo_d = get("orthofinder", "default"), get("fastoma", "default_declared")
    if of_d and fo_d:
        for scope in ("all", "headline"):
            sel = (lambda rs: rs) if scope == "all" else (lambda rs: [r for r in rs if hl[r["id"]][0]])
            cs_ = {m: compare_methods(sel(of_d), sel(fo_d), m, n_boot=a.n_boot) for m in PRIMARY}
            hp = holm({m: c["wilcoxon_p"] for m, c in cs_.items()})
            for m, c in cs_.items():
                cmp_rows.append([scope, "paired:orthofinder_default-fastoma_default_declared", "", "", "comparison", m,
                                 fnum(c["mean_diff"]), fnum(c["ci_low"]), fnum(c["ci_high"]), fnum(c["min_diff"]),
                                 fnum(c["max_diff"]), c["n"]])
                cmp_rows.append([scope, "paired:wilcoxon_p/holm", "", "", "comparison", m, fnum(c["wilcoxon_p"]),
                                 fnum(hp[m]), "", "", "", c["n"]])
    with open(out / "summary.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["scope", "track", "method", "variant", "kind", "metric", "mean", "ci_low", "ci_high", "min", "max", "n"])
        w.writerows(rows + cmp_rows)

    # provenance
    fdir = Path(a.families)
    frozen = (fdir / MANIFEST).exists()
    prov = {"hobench_git": _git_state(ROOT), "script": "scripts/evaluate_panel.py", "panel": P, "mode": mode,
            "panel_tree_sha256": sha256_file(ptree), "panel_species_tsv_sha256": sha256_file(ROOT / "panels" / P / "species.tsv"),
            "families_tree_hash": tree_hash(manifest_entries(fdir)), "families_manifest_present": frozen,
            "families_manifest_differences": check_manifest(fdir) if frozen else None,
            "n_families": len(fams), "n_non_headline": sum(not h[0] for h in hl.values()),
            "seeds_relabel": SEEDS[:N_RELABEL_SEEDS], "n_boot": a.n_boot, "runs": prov_runs}
    (out / "provenance.json").write_text(json.dumps(prov, indent=1, default=str))

    # summary.md
    L = [f"# Panel {P} / {mode}: evaluation summary", "",
         f"hobench git {prov['hobench_git']['commit'][:10]} (dirty={prov['hobench_git']['dirty']}); "
         f"{len(fams)} families, {prov['n_non_headline']} non-headline; families tree hash {prov['families_tree_hash'][:12]}.",
         "Mean over items (collection = one item), scope=all families; CI = stratified bootstrap, none for n<5 "
         "(then min-max). Details: summary.tsv (also scope=headline and per-collection tracks).", "",
         "| method | variant | " + " | ".join(PRIMARY) + " | n items |", "|---|---|" + "---|" * (len(PRIMARY) + 1)]
    for meth, var, kind, res, _ in methods:
        c = aggs[(meth, var, "all")]["combined"]
        cells = []
        for m in PRIMARY:
            d = c[m]
            if d["mean"] is None:
                cells.append("NA")
            elif d["ci_low"] is not None:
                cells.append(f"{d['mean']:.3f} [{d['ci_low']:.3f},{d['ci_high']:.3f}]")
            else:
                cells.append(f"{d['mean']:.3f} ({d['min']:.2f}-{d['max']:.2f})")
        L.append(f"| {meth} | {var} | " + " | ".join(cells) + f" | {c[PRIMARY[0]]['n']} |")
    L += ["", "Unit-level (families track + each collection's units; mean over units):", ""]
    L += ["| method | variant | collection | " + " | ".join(PRIMARY) + " | n units |", "|---|---|---|" + "---|" * (len(PRIMARY) + 1)]
    for meth, var, kind, res, _ in methods:
        ag = aggs[(meth, var, "all")]
        for c, d in ag["collections"].items():
            L.append(f"| {meth} | {var} | {c} | " + " | ".join("NA" if d[m]["mean"] is None else f"{d[m]['mean']:.3f}" for m in PRIMARY)
                     + f" | {d[PRIMARY[0]]['n']} |")
    # sanity
    L += ["", "## Sanity checks (per-family min/max over families with a defined value)", ""]
    bad = []
    for meth, var, kind, res, _ in methods:
        if kind != "control":
            continue
        for m in PRIMARY:
            xs = [family_scalars(r)[m] for r in res if family_scalars(r)[m] is not None]
            if xs:
                L.append(f"- {meth} {m}: mean {np.mean(xs):.4f}, min {min(xs):.4f}, max {max(xs):.4f}, n={len(xs)}")
                if meth == "positive_control" and min(xs) < 1 - 1e-9:
                    bad.append((m, [r["id"] for r in res if (family_scalars(r)[m] or 1) < 1 - 1e-9][:10]))
    L.append("")
    L.append("POSITIVE CONTROL: " + ("PERFECT on all primary endpoints." if not bad else f"NOT perfect: {bad}"))
    (out / "summary.md").write_text("\n".join(L) + "\n")
    log("done ->", out)
    print("\n".join(L[:12 + 7]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
