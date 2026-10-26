"""Evaluate a predicted set of HOGs against the reference families of a panel.

Usage::

    python -m hobench.evaluate --panel-tree panel.nwk --families families/ \\
        --panel-dir panels/Y --prediction FastOMA_HOGs.orthoxml --format orthoxml \\
        --tool-species-tree species_tree_checked.nwk --out results.json

Writes ``<out>.json`` (per-family results, skill scores, two-stage aggregates with
bootstrap CIs, provenance: git commit + data hashes) and ``<out>.tsv`` (one row per family).
Inconsistent inputs raise; nothing is silently repaired. Validation against the pinned panel
FASTAs (``--panel-dir``) is mandatory unless ``--no-fasta-check`` is given explicitly (and is
then recorded in the output).
"""
from __future__ import annotations

import argparse
import csv
import json
import subprocess
import sys
from pathlib import Path
from typing import Optional, Sequence

import numpy as np

from . import controls
from .freeze import MANIFEST, check_manifest, manifest_entries, tree_hash
from .ids import make_normalizer
from .metrics import PRIMARY, SCALARS, ScoringError, evaluate_family, family_scalars, panel_context
from .panel import load_panel, sha256_file, validate_members
from .predictions import Prediction, read_prediction
from .reference import ReferenceFamily, load_families, reconcile, validate_panel
from .seeds import SEEDS
from .speciestree import SpeciesTree
from .stats import aggregate

N_RELABEL_SEEDS = 10


def baseline_scores(families: Sequence[ReferenceFamily], st: SpeciesTree, ctx=None,
                    metrics: Sequence[str] = PRIMARY) -> dict[str, dict[str, Optional[float]]]:
    """B4: per family and metric, b = max(all-singletons, one-big-HOG, mean over 10 seeds of
    the within-species relabelled positive control)."""
    ctx = ctx or panel_context(families, st)
    rfs = [reconcile(f, st, check=False) for f in families]

    def scores(pred: Prediction) -> dict[str, dict]:
        return {rf.id: family_scalars(evaluate_family(rf, st, pred, ctx, per_gene=False, compact=True)) for rf in rfs}

    single = scores(controls.all_singletons(families))
    big = scores(controls.one_big_hog(families, st))
    pos = controls.reference_prediction(families, st, check=False)
    rel = [scores(controls.random_relabel(pos, s)) for s in SEEDS[:N_RELABEL_SEEDS]]
    out: dict[str, dict] = {}
    for rf in rfs:
        d = {}
        for m in metrics:
            rv = [r[rf.id][m] for r in rel if r[rf.id][m] is not None]
            cands = [x for x in (single[rf.id][m], big[rf.id][m], float(np.mean(rv)) if rv else None) if x is not None]
            d[m] = max(cands) if cands else None
        out[rf.id] = d
    return out


def skill(score: Optional[float], base: Optional[float]) -> Optional[float]:
    """(s − b)/(1 − b); None if undefined (missing values or b ≥ 1)."""
    if score is None or base is None or base >= 1.0:
        return None
    return (score - base) / (1.0 - base)


def evaluate_panel(st: SpeciesTree, families: Sequence[ReferenceFamily], pred: Prediction,
                   n_boot: int = 2000, seed: int = SEEDS[0], per_gene: bool = True,
                   with_skill: bool = True) -> dict:
    fatal = pred.fatal_problems(st)
    if fatal:
        raise ScoringError("prediction is inconsistent with the panel species tree:\n  " + "\n  ".join(fatal[:20]))
    ctx = panel_context(families, st)
    missing_counts = sorted(ctx.count_species - pred.count_species)
    if missing_counts:
        raise ScoringError(f"prediction did not record per-HOG counts for complete species {missing_counts}")
    results = []
    for f in families:
        results.append(evaluate_family(reconcile(f, st, check=False), st, pred, ctx,
                                       per_gene=per_gene, compact=bool(f.collection)))
    if with_skill:
        base = baseline_scores(families, st, ctx)
        for r in results:
            sc = family_scalars(r)
            r["baseline"] = base[r["id"]]
            r["skill"] = {m: skill(sc[m], base[r["id"]][m]) for m in PRIMARY}
    return {"families": results, "aggregate": aggregate(results, n_boot=n_boot, seed=seed),
            "prediction_meta": {k: v for k, v in pred.meta.items() if k != "warnings"},
            "prediction_warnings": pred.meta.get("warnings", [])[:200],
            "prediction_nesting_problems": pred.nesting_problems(st)[:50]}


def _read_map(path: Optional[str]) -> Optional[dict[str, str]]:
    if not path:
        return None
    out = {}
    with open(path) as fh:
        for row in csv.reader(fh, delimiter="\t"):
            if row and not row[0].startswith("#") and len(row) >= 2:
                out[row[0]] = row[1]
    return out


def _git_state(cwd: Path) -> dict:
    def run(*cmd):
        try:
            return subprocess.run(cmd, cwd=cwd, capture_output=True, text=True, check=True).stdout.strip()
        except (subprocess.CalledProcessError, FileNotFoundError):
            return None
    commit = run("git", "rev-parse", "HEAD")
    status = run("git", "status", "--porcelain", "--untracked-files=no")
    return {"commit": commit, "dirty": bool(status) if status is not None else None}


def _path_hash(p: Path) -> str:
    if p.is_file():
        return sha256_file(p)
    entries = [(sha256_file(f), f.relative_to(p).as_posix()) for f in sorted(p.rglob("*"))
               if f.is_file() and (f.suffix in (".tsv", ".txt", ".orthoxml", ".xml") or f.name.startswith("N"))]
    return tree_hash(entries)


def provenance(args, families_dir: Path) -> dict:
    fam_entries = manifest_entries(families_dir)
    frozen = (families_dir / MANIFEST).exists()
    diffs = check_manifest(families_dir) if frozen else None
    out = {
        "hobench_git": _git_state(Path(__file__).resolve().parent),
        "panel_tree_sha256": sha256_file(args.panel_tree),
        "prediction_sha256": _path_hash(Path(args.prediction)),
        "families_tree_hash": tree_hash(fam_entries),
        "families_manifest_present": frozen,
        "families_manifest_sha256": sha256_file(families_dir / MANIFEST) if frozen else None,
        "families_manifest_differences": diffs,
    }
    if args.panel_dir:
        out["panel_species_tsv_sha256"] = sha256_file(Path(args.panel_dir) / "species.tsv")
    if getattr(args, "tool_species_tree", None):
        out["tool_species_tree_sha256"] = sha256_file(args.tool_species_tree)
    return out


def write_tsv(res: dict, path: Path) -> None:
    keys = list(SCALARS) + [f"skill_{m}" for m in PRIMARY]
    cols = ["family", "name", "collection", "phenomena", "family_level", "n_members", "n_asserted_pairs"] + keys
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(cols)
        for r in res["families"]:
            sc = family_scalars(r)
            w.writerow([r["id"], r["name"], r.get("collection") or "", ",".join(r["phenomena"]), r["family_level"],
                        r["n_members"], r["n_asserted_pairs"]]
                       + ["" if sc.get(k) is None else f"{sc[k]:.6g}" for k in keys])


def main(argv: Optional[Sequence[str]] = None) -> int:
    ap = argparse.ArgumentParser(prog="python -m hobench.evaluate", description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--panel-tree", required=True, help="rooted Newick species tree with named levels")
    ap.add_argument("--families", required=True, help="directory containing <FAMILY_ID>/family.yaml")
    ap.add_argument("--panel", default=None, help="only families whose family.yaml 'panel' equals this")
    ap.add_argument("--panel-dir", default=None, help="panels/<P> with species.tsv (pinned FASTAs)")
    ap.add_argument("--no-fasta-check", action="store_true",
                    help="explicitly skip checking members against the pinned panel FASTAs")
    ap.add_argument("--prediction", required=True, help="OrthoXML file | OrthoFinder results/HOG dir | TSV")
    ap.add_argument("--format", required=True, choices=["orthoxml", "orthofinder", "tsv"])
    ap.add_argument("--id-normalizer", default="auto",
                    help="identity|first_token|uniprot|auto|map:<tsv>, '+'-combinable (default auto)")
    ap.add_argument("--species-map", help="TSV: species name in tool output -> species code in the tree")
    ap.add_argument("--level-map", help="TSV (orthoxml): level name in file -> level name in the tool tree")
    ap.add_argument("--tool-species-tree", help="(orthoxml) species tree given to the tool, e.g. FastOMA "
                                                "species_tree_checked.nwk; levels are mapped by clade")
    ap.add_argument("--dup-placement", default="declared", choices=["declared", "late"],
                    help="(orthoxml) duplication placement: declared level (primary) or pyHam-late (sensitivity)")
    ap.add_argument("--of-species-tree", help="OrthoFinder SpeciesTree_rooted_node_labels.txt (auto-detected)")
    ap.add_argument("--allow-induced-levels", action="store_true",
                    help="(orthofinder) induce unmatched levels; only for perturbed species-tree runs")
    ap.add_argument("--flat", action="store_true",
                    help="(tsv) single flat partition: project each group to every level (no duplications)")
    ap.add_argument("--orthoxml-id-attr", default="protId", help="gene attribute holding the ID (protId)")
    ap.add_argument("--out", required=True, help="output path (.json or .tsv; both are written)")
    ap.add_argument("--n-boot", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=SEEDS[0])
    ap.add_argument("--no-per-gene", action="store_true", help="omit per-gene rows from the JSON")
    ap.add_argument("--no-skill", action="store_true", help="skip baseline-normalised skill (12 extra scorings)")
    a = ap.parse_args(argv)

    st = SpeciesTree.from_newick(a.panel_tree)
    families = load_families(a.families, panel=a.panel)
    if not families:
        print(f"no families found in {a.families}" + (f" for panel {a.panel}" if a.panel else ""), file=sys.stderr)
        return 2
    rep = validate_panel(families, st)
    if a.panel_dir:
        panel = load_panel(a.panel_dir)
        for e in panel.verify_checksums():
            rep.errors.append(e)
        for f in families:
            rep.extend(validate_members(f, panel))
    elif not a.no_fasta_check:
        print("ERROR: --panel-dir is required (or pass --no-fasta-check explicitly)", file=sys.stderr)
        return 2
    for wmsg in rep.warnings:
        print("WARNING:", wmsg, file=sys.stderr)
    if not rep.ok:
        print(rep, file=sys.stderr)
        return 1
    keep, count_species = set(), set()
    for f in families:
        keep.update(f.member_ids)
        keep.update(f.decoys)
        count_species.update(f.complete_species)
    kw: dict = {"normalize": make_normalizer(a.id_normalizer), "keep_genes": keep, "count_species": count_species}
    if a.format == "orthoxml":
        kw.update(species_map=_read_map(a.species_map), level_map=_read_map(a.level_map),
                  id_attr=a.orthoxml_id_attr, tool_species_tree=a.tool_species_tree, dup_placement=a.dup_placement)
    elif a.format == "orthofinder":
        kw.update(species_map=_read_map(a.species_map), of_species_tree=a.of_species_tree,
                  allow_induced_levels=a.allow_induced_levels)
    elif a.flat:
        kw["flat"] = True
    pred = read_prediction(a.format, a.prediction, st, **kw)
    res = evaluate_panel(st, families, pred, n_boot=a.n_boot, seed=a.seed, per_gene=not a.no_per_gene,
                         with_skill=not a.no_skill)
    res["meta"] = {"panel_tree": a.panel_tree, "families_dir": a.families, "panel": a.panel,
                   "panel_dir": a.panel_dir, "fasta_check": not a.no_fasta_check or bool(a.panel_dir),
                   "prediction": a.prediction, "format": a.format, "id_normalizer": a.id_normalizer,
                   "dup_placement": a.dup_placement if a.format == "orthoxml" else None,
                   "allow_induced_levels": a.allow_induced_levels, "n_boot": a.n_boot, "seed": a.seed,
                   "provenance": provenance(a, Path(a.families))}
    out = Path(a.out)
    stem = out.with_suffix("") if out.suffix in (".json", ".tsv") else out
    stem.parent.mkdir(parents=True, exist_ok=True)
    jpath, tpath = stem.with_suffix(".json"), stem.with_suffix(".tsv")
    jpath.write_text(json.dumps(res, indent=1, default=str))
    write_tsv(res, tpath)
    comb = res["aggregate"]["combined"]
    print(f"{len(families)} families evaluated -> {jpath}, {tpath}", file=sys.stderr)
    for k in PRIMARY:
        m = comb[k]
        if m["mean"] is not None:
            ci = f"[{m['ci_low']:.4f}, {m['ci_high']:.4f}]" if m["ci_low"] is not None else "(no CI, n<5)"
            print(f"  {k:24s} {m['mean']:.4f}  {ci}  n={m['n']}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
