#!/usr/bin/env python
"""Perturbation-aware scoring of E2 robustness runs (Y, full mode) against the clean Y references (DESIGN sec. 7, E2).

  python scripts/evaluate_robustness.py --pert frag|contam|taxrm|alttree [--run orthofinder_default ...]

Prediction runs: runs/Y_<pert>/full/<run>/ (manifest status must be 'done').  Scored on the unperturbed Y
reference families with the clean Y tree (taxrm: reference pruned to remaining species; alttree: prediction lifted to
the Y levels, level B-cubed F1 averaged over levels SHARED with the alternative tree only).  Per family and primary
endpoint: clean, perturbed, delta; plus gene-level split (affected vs unaffected genes, mean per-gene fraction of
pairs correct) and, for fragments, both credit rules.
 frag    : `<acc>_f1/_f2` stand in for the parent (all relations inherited).  credit=either: per-family best of the
           two single-fragment readings (parent credited if either fragment is placed correctly); credit=both: per-family
           worst of the two (both fragments correct).  Affected genes = fragmented parents.
 contam  : contaminant IDs are excluded from reference scoring; their placement (HOGs shared with reference genes) is
           written to results/Y_contam_full/contaminant_placement.tsv.  Affected genes = reference genes sharing a
           root-level HOG with a contaminant.
 taxrm   : removed species' genes are dropped from the references (tree leaves, members); the prediction (pruned tree) is
           lifted to Y levels (a removed level = copy of the pruned level with the same remaining species, or empty).
           Affected genes = remaining species under the parent of the removed lineage's LCA.
 alttree : prediction parsed on the alternative tree, lifted to Y levels (identical clades copied, others induced from the
           smallest containing level restricted to the clade).  Affected = species in clades that differ between the trees.
Noise floor: |OF default - OF rep2| on clean Y if a replicate exists (runs/Y/full/orthofinder_*rep2*), written to
results/robustness/noise_floor_orthofinder.tsv.
Outputs: results/robustness/<pert>_<tool>.tsv, results/Y_<pert>_full/summary.md.
"""
from __future__ import annotations

import argparse, copy, csv, json, sys
from collections import defaultdict
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / "scripts"))
import evaluate_panel as ep                                                      # noqa: E402
from hobench.ids import make_normalizer                                         # noqa: E402
from hobench.metrics import PRIMARY, evaluate_family, family_scalars, panel_context   # noqa: E402
from hobench.predictions import Prediction                                      # noqa: E402
from hobench.reference import load_families, reconcile                          # noqa: E402
from hobench.seeds import SEEDS                                                 # noqa: E402
from hobench.speciestree import SpeciesTree                                     # noqa: E402

METRICS = PRIMARY
BOOT = 2000


# ---------------------------------------------------------------- prediction transforms
def rebuild(pred: Prediction, levels: dict[str, dict[str, str]]) -> Prediction:
    new = Prediction(meta=dict(pred.meta), count_species=set(pred.count_species), gene_species=dict(pred.gene_species))
    for lv, m in levels.items():
        for g, h in m.items():
            new.add(lv, h, g, species=pred.gene_species.get(g))
    return new


def relabel(pred, mapping, drop=()):
    """gene -> parent gene (mapping), genes in drop removed; first assignment wins."""
    drop = set(drop)
    lv = {}
    for L, m in pred.levels.items():
        d = {}
        for g, h in m.items():
            if g in drop:
                continue
            d.setdefault(mapping.get(g, g), h)
        lv[L] = d
    out = rebuild(pred, lv)
    out.gene_species = {mapping.get(g, g): s for g, s in pred.gene_species.items() if g not in drop}
    return out


def lift(pred: Prediction, src: SpeciesTree, dst: SpeciesTree, keep_species=None):
    """Re-express a prediction made on tree `src` on the levels of `dst` (see module docstring). Returns
    (prediction, shared_levels, induced_levels)."""
    keep = set(keep_species) if keep_species is not None else set(dst.species)
    src_clades = {}
    for L in src.names:
        if not src.is_species(L):
            src_clades[L] = frozenset(src.clade_species(L))
    by_set = {}
    for L, c in src_clades.items():
        by_set.setdefault(c, []).append(L)
    out, shared, induced = {}, [], {}
    for L in dst.names:
        if dst.is_species(L):
            continue
        C = frozenset(dst.clade_species(L)) & keep
        if len(C) < 2:
            continue
        cand = [m for m in by_set.get(C, []) if m in pred.levels]
        if cand:
            out[L] = pred.levels[cand[0]]
            if frozenset(dst.clade_species(L)) == C:
                shared.append(L)
            continue
        sup = sorted((m for m, c in src_clades.items() if C <= c and m in pred.levels), key=lambda m: len(src_clades[m]))
        if not sup:
            continue
        m = sup[0]
        out[L] = {g: h for g, h in pred.levels[m].items() if pred.gene_species.get(g) in C}
        induced[L] = m
    return rebuild(pred, out), shared, induced


def prune_family(f, drop_sp):
    """Reference family restricted to species not in drop_sp (leaves removed, unary nodes spliced)."""
    g = copy.deepcopy(f)
    gsp = f.gene_species
    root = g.tree
    for nd in list(root.postorder()):
        if nd.is_leaf and gsp.get(nd.name) in drop_sp and nd.parent is not None:
            nd.parent.children.remove(nd)
    changed = True
    while changed:
        changed = False
        for nd in list(root.postorder()):
            if nd.parent is None:
                continue
            if not nd.children and gsp.get(nd.name) is None:      # emptied internal node (also 'leaf' by shape)
                nd.parent.children.remove(nd); changed = True
            elif len(nd.children) == 1:
                p, c = nd.parent, nd.children[0]
                p.children[p.children.index(nd)] = c; c.parent = p; changed = True
        if len(root.children) == 1 and root.children[0].children:
            root = root.children[0]; root.parent = None; changed = True
    g.tree = root
    keep_g = {n.name for n in root.leaves()}
    g.members = g.members[g.members["gene_id"].astype(str).isin(keep_g)].reset_index(drop=True)
    g.meta = dict(g.meta)
    g.meta["complete_species"] = [s for s in f.complete_species if s not in drop_sp]
    return g


# ---------------------------------------------------------------- scoring helpers
def score_pred(st, rfs, pred, ctx):
    out = {}
    for rf in rfs:
        r = evaluate_family(rf, st, pred, ctx, per_gene=True, compact=False)
        genes = {x["gene"]: x["frac_correct"] for x in r.get("genes", [])}
        sc = family_scalars(r)
        out[r["id"]] = ({m: sc.get(m) for m in METRICS}, genes)
        if not (rf.family.complete_species or rf.family.decoys):
            out[r["id"]][0]["membership_f1"] = None
    return out


def combine(a, b, fn):
    out = {}
    for fid in a:
        sa, ga = a[fid]; sb, gb = b[fid]
        s = {m: (None if sa[m] is None or sb[m] is None else fn(sa[m], sb[m])) for m in METRICS}
        g = {k: (None if ga[k] is None or gb.get(k) is None else fn(ga[k], gb[k])) for k in ga}
        out[fid] = (s, g)
    return out


def gene_split(genes, affected):
    a = [v for k, v in genes.items() if k in affected and v is not None]
    u = [v for k, v in genes.items() if k not in affected and v is not None]
    return (float(np.mean(a)) if a else None, float(np.mean(u)) if u else None, len(a), len(u))


def fmt(x):
    return "" if x is None else f"{x:.6g}"


def boot_ci(d, rng):
    d = np.asarray(d, float)
    if len(d) < 5:
        return None, None
    bm = d[rng.integers(0, len(d), size=(BOOT, len(d)))].mean(axis=1)
    return float(np.percentile(bm, 2.5)), float(np.percentile(bm, 97.5))


# ---------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pert", required=True, choices=["frag", "contam", "taxrm", "alttree"])
    ap.add_argument("--run", nargs="*", default=None, help="run dir names (default: all done runs in runs/Y_<pert>/full)")
    ap.add_argument("--cache", default="/tmp/hobench_eval_cache")
    ap.add_argument("--name", default=None, help="panel/run suffix if not equal to --pert (e.g. contam200)")
    a = ap.parse_args()
    pert = a.pert
    pname = a.name or pert
    rdir = ROOT / "runs" / f"Y_{pname}" / "full"
    pdir = ROOT / "panels" / f"Y_{pname}"
    pj = json.loads((pdir / "perturbation.json").read_text())
    out_rob = ROOT / "results" / "robustness"; out_rob.mkdir(parents=True, exist_ok=True)
    out_dir = ROOT / "results" / f"Y_{pname}_full"; out_dir.mkdir(parents=True, exist_ok=True)

    st = SpeciesTree.from_newick(str(ROOT / "panels" / "Y" / "species_tree.nwk"))
    st_p = SpeciesTree.from_newick(str(pdir / "species_tree.nwk"))
    fams = load_families(str(ROOT / "families"), panel="Y")
    hl = {f.id: ep.headline_info(f) for f in fams}
    tk = {f.id: ep.track_of(f) for f in fams}
    # perturbation bookkeeping
    frag = {fr["id"]: c["orig_id"] for c in pj.get("changes", []) if pert == "frag" for fr in c["fragments"]}
    contam = {c["new_id"]: c for c in pj.get("changes", []) if pert == "contam"}
    removed = set(pj.get("removed_species", []))
    keep_sp = set(st.species) - removed
    keep = set(); 
    for f in fams:
        keep |= set(f.member_ids) | set(f.decoys)
    kw = {"normalize": make_normalizer("auto"), "keep_genes": keep | set(frag) | set(contam),
          "count_species": set(c for f in fams for c in f.complete_species) - removed}
    kw_clean = {"normalize": make_normalizer("auto"), "keep_genes": keep, "count_species": set(c for f in fams for c in f.complete_species)}

    # references / context
    if pert == "taxrm":
        fams_s = [prune_family(f, removed) for f in fams]
        fams_s = [f for f in fams_s if len(f.members) >= 2 and f.tree.children]   # drop families emptied by removal
    else:
        fams_s = fams
    ctx = panel_context(fams_s, st)
    rfs = [reconcile(f, st, check=False) for f in fams_s]

    # affected genes
    all_sp = {f.id: f.gene_species for f in fams}
    if pert == "frag":
        affected_by = {f.id: {g for g in f.member_ids if g in set(frag.values())} for f in fams}
    elif pert == "taxrm":
        lca = st.lca(*removed); par = st.parent(lca) or lca
        zone = set(st.clade_species(par)) - removed
        affected_by = {f.id: {g for g, s in all_sp[f.id].items() if s in zone} for f in fams}
    elif pert == "alttree":
        sp = set()
        for _, c in pj["clades"].get("only_parent", []) + pj["clades"].get("only_new", []):
            sp |= set(c)
        affected_by = {f.id: {g for g, s in all_sp[f.id].items() if s in sp} for f in fams}
    else:
        affected_by = None   # contam: per run (depends on placement)

    runs = sorted(p for p in rdir.iterdir() if (p / "manifest.json").exists()) if rdir.exists() else []
    if a.run:
        runs = [r for r in runs if r.name in a.run]
    clean_cache = {}
    rng = np.random.default_rng(SEEDS[0])
    md = [f"# Robustness: Y {pert} (full mode), perturbed vs clean Y full", "",
          f"Perturbation: {pj['perturbation']}, seed {pj.get('seed')}; panel {pdir.name}. Descriptive: mean over scored units "
          "(YGOB units + Y03; control track excluded), delta = perturbed - clean; CI = 95% bootstrap over units. "
          "See scripts/evaluate_robustness.py docstring for the rules.", ""]
    nf = None
    for run in runs:
        man = json.loads((run / "manifest.json").read_text())
        if man.get("status") != "done":
            print(f"skip {run.name}: manifest status {man.get('status')!r}", file=sys.stderr); continue
        for tool, var, pred, src, _ in ep.load_run(run, st_p, kw, Path(a.cache), ["declared"]):
            name = f"{tool}_{var}"
            if pert == "frag":
                preds = {}
                for which in ("f1", "f2"):
                    mp = {fid: par for fid, par in frag.items() if fid.endswith("_" + which)}
                    drop = [fid for fid in frag if not fid.endswith("_" + which)]
                    preds[which] = relabel(pred, mp, drop)
                variants = {"f1": preds["f1"], "f2": preds["f2"]}
            elif pert == "contam":
                variants = {"": relabel(pred, {}, list(contam))}
            elif pert == "taxrm":
                variants = {"": lift(pred, st_p, st, keep_sp)[0]}
            else:
                lp, shared, induced = lift(pred, st_p, st)
                variants = {"": lp}
            sc = {k: score_pred(st, rfs, p, ctx) for k, p in variants.items()}
            if pert == "frag":
                perturbed = {"either": combine(sc["f1"], sc["f2"], max), "both": combine(sc["f1"], sc["f2"], min)}
            else:
                perturbed = {"": sc[""]}
            # alttree: level F1 over shared levels only
            if pert == "alttree":
                shared = set(shared)
                for rf in rfs:
                    r = evaluate_family(rf, st, variants[""], ctx, per_gene=False, compact=False)
                    xs = [l["bcubed_f1"] for l in r["per_level"] if l["level"] in shared and l["bcubed_f1"] is not None]
                    perturbed[""][r["id"]][0]["level_bcubed_f1"] = float(np.mean(xs)) if xs else None
            # clean
            if (tool, var) not in clean_cache:
                crun = ROOT / "runs" / "Y" / "full" / run.name
                cp = [x for x in ep.load_run(crun, st, kw_clean, Path(a.cache), ["declared"]) if x[0] == tool and x[1] == var][0][2]
                cpred = cp
                if pert == "taxrm":    # clean prediction restricted to remaining species, on the pruned references
                    cpred = relabel(cp, {}, [g for g, s in cp.gene_species.items() if s in removed])
                clean_cache[(tool, var)] = score_pred(st, rfs, cpred, ctx)
            clean = clean_cache[(tool, var)]
            # contaminant placement + affected genes
            aff = affected_by
            if pert == "contam":
                aff, place = {f.id: set() for f in fams}, []
                g2f = {g: f.id for f in fams for g in f.member_ids}
                roots = [L for L in pred.levels if st_p.depth_of(L) == 0] or [min(pred.levels, key=st_p.depth_of)]
                for cid, c in contam.items():
                    joined = {}
                    jl = ""
                    for L in sorted(pred.levels, key=lambda x: -st_p.depth_of(x)):     # deepest first
                        h = pred.levels[L].get(cid)
                        if h is None:
                            continue
                        refs = [g for g in pred.hog_members(L, h) if g in g2f]
                        if refs and not joined:
                            jl = L; joined = {g2f[g] for g in refs}
                            jn = len(refs); jsp = sorted({pred.gene_species.get(g, "?") for g in refs})
                    rl = roots[0]; hr = pred.levels[rl].get(cid)
                    rfam = defaultdict(int)
                    if hr is not None:
                        for g in pred.hog_members(rl, hr):
                            if g in g2f:
                                rfam[g2f[g]] += 1
                                aff[g2f[g]].add(g)
                    place.append([cid, c["source_id"], tool, var, jl, ";".join(sorted(joined)), len(joined),
                                  ";".join(f"{k}:{v}" for k, v in sorted(rfam.items())[:8]), sum(rfam.values())])
                with open(out_dir / f"contaminant_placement_{name}.tsv", "w") as fh:
                    w = csv.writer(fh, delimiter="\t", lineterminator="\n")
                    w.writerow(["contaminant", "source_id", "tool", "variant", "deepest_level_sharing_ref_HOG", "families_at_that_level",
                                "n_families", "root_level_families:n_genes", "n_ref_genes_in_root_HOG"])
                    w.writerows(place)
                md += [f"Contaminant placement ({name}): {sum(1 for p in place if p[6] > 0)}/{len(place)} contaminants share a HOG "
                       f"with reference genes at some level; {sum(1 for p in place if p[8] > 0)}/{len(place)} at the root level "
                       f"(details: contaminant_placement_{name}.tsv).", ""]
            # rows
            rows, acc = [], defaultdict(list)
            for credit, pc in perturbed.items():
                for fid in clean:
                    cs, cg = clean[fid]; ps, pg = pc[fid]
                    coll = "YGOB" if tk[fid] == "collection" else ""
                    base = [fid, tk[fid], int(hl[fid][0]), pert, tool, var, credit]
                    for m in METRICS:
                        d = None if cs[m] is None or ps[m] is None else ps[m] - cs[m]
                        rows.append(base + [m, fmt(cs[m]), fmt(ps[m]), fmt(d), "", ""])
                        if d is not None and tk[fid] != "control":
                            acc[(credit, m)].append(d)
                    ca, cu, na, nu = gene_split(cg, aff[fid]); pa, pu, _, _ = gene_split(pg, aff[fid])
                    for lab, c_, p_ in (("gene_frac_correct_affected", ca, pa), ("gene_frac_correct_unaffected", cu, pu)):
                        d = None if c_ is None or p_ is None else p_ - c_
                        rows.append(base + [lab, fmt(c_), fmt(p_), fmt(d), na if "aff" in lab and "un" not in lab else nu, ""])
                        if d is not None and tk[fid] != "control":
                            acc[(credit, lab)].append(d)
            with open(out_rob / f"{pname}_{name}.tsv", "w", newline="") as fh:
                w = csv.writer(fh, delimiter="\t", lineterminator="\n")
                w.writerow(["family", "track", "headline", "pert", "tool", "variant", "credit", "metric", "clean", "perturbed",
                            "delta", "n_genes", "note"])
                w.writerows(rows)
            md += [f"## {name}" + (f" (induced levels: {len(induced)}; shared: {len(shared)})" if pert == "alttree" else ""), "",
                   "| credit | metric | n units | mean delta | 95% CI |", "|---|---|---|---|---|"]
            for (credit, m), ds in acc.items():
                lo, hi = boot_ci(ds, rng)
                md.append(f"| {credit or '-'} | {m} | {len(ds)} | {np.mean(ds):+.4f} | " + ("NA" if lo is None else f"[{lo:+.4f}, {hi:+.4f}]") + " |")
            md.append("")
    # noise floor
    reps = sorted((ROOT / "runs" / "Y" / "full").glob("orthofinder_*rep2*"))
    if reps and (reps[0] / "manifest.json").exists() and json.loads((reps[0] / "manifest.json").read_text()).get("status") == "done":
        tool, var, p2, _, _ = ep.load_run(reps[0], st, kw_clean, Path(a.cache), ["declared"])[0]
        s2 = score_pred(st, [reconcile(f, st, check=False) for f in fams], p2, panel_context(fams, st)) if pert != "taxrm" else None
        if s2 is not None:
            base = clean_cache.get(("orthofinder", "default"))
            if base:
                with open(out_rob / "noise_floor_orthofinder.tsv", "w") as fh:
                    fh.write("family\tmetric\tdefault\trep2\tdelta\n")
                    for fid in base:
                        for m in METRICS:
                            x, y = base[fid][0][m], s2[fid][0][m]
                            fh.write(f"{fid}\t{m}\t{fmt(x)}\t{fmt(y)}\t{fmt(None if x is None or y is None else y - x)}\n")
                md += ["Noise floor (OF default vs rep2 on clean Y): results/robustness/noise_floor_orthofinder.tsv", ""]
    else:
        md += ["Noise floor: OF rep2 replicate not available yet.", ""]
    (out_dir / "summary.md").write_text("\n".join(md) + "\n")
    print("\n".join(md[:40]))


if __name__ == "__main__":
    main()
