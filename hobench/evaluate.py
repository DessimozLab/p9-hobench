"""Evaluate a predicted set of HOGs against the reference families of a panel.

Usage::

    python -m hobench.evaluate --panel-tree panel.nwk --families families/ \\
        --prediction FastOMA_HOGs.orthoxml --format orthoxml --out results.json

Writes ``<out>.json`` (full per-family results + aggregates with bootstrap CIs) and
``<out>.tsv`` (one row per family with the headline scalars). ``--out`` may end in ``.json``
or ``.tsv``; the other file is written next to it.
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path
from typing import Optional, Sequence

from .ids import make_normalizer
from .metrics import SCALARS, aggregate, evaluate_family, family_scalars, panel_context
from .predictions import Prediction, read_prediction
from .reference import ReferenceFamily, load_families, reconcile, validate_panel
from .seeds import SEEDS
from .speciestree import SpeciesTree


def evaluate_panel(st: SpeciesTree, families: Sequence[ReferenceFamily], pred: Prediction,
                   n_boot: int = 2000, seed: int = SEEDS[0], per_gene: bool = True) -> dict:
    ctx = panel_context(families)
    results = [evaluate_family(reconcile(f, st, check=False), st, pred, ctx, per_gene=per_gene)
               for f in families]
    return {"families": results, "aggregate": aggregate(results, n_boot=n_boot, seed=seed),
            "prediction_meta": {k: v for k, v in pred.meta.items() if k != "warnings"},
            "prediction_warnings": pred.meta.get("warnings", [])[:200],
            "prediction_problems": pred.check(st)[:50]}


def _read_map(path: Optional[str]) -> Optional[dict[str, str]]:
    if not path:
        return None
    out = {}
    with open(path) as fh:
        for row in csv.reader(fh, delimiter="\t"):
            if row and not row[0].startswith("#") and len(row) >= 2:
                out[row[0]] = row[1]
    return out


def write_tsv(res: dict, path: Path) -> None:
    cols = ["family", "name", "phenomena", "family_level", "n_members", "n_asserted_pairs"] + list(SCALARS)
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(cols)
        for r in res["families"]:
            sc = family_scalars(r)
            w.writerow([r["id"], r["name"], ",".join(r["phenomena"]), r["family_level"], r["n_members"],
                        r["n_asserted_pairs"]] + ["" if sc[k] is None else f"{sc[k]:.6g}" for k in SCALARS])


def main(argv: Optional[Sequence[str]] = None) -> int:
    ap = argparse.ArgumentParser(prog="python -m hobench.evaluate", description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--panel-tree", required=True, help="rooted Newick species tree with named levels")
    ap.add_argument("--families", required=True, help="directory containing <FAMILY_ID>/family.yaml")
    ap.add_argument("--panel", default=None, help="only families whose family.yaml 'panel' equals this")
    ap.add_argument("--prediction", required=True, help="OrthoXML file | OrthoFinder results/HOG dir | TSV")
    ap.add_argument("--format", required=True, choices=["orthoxml", "orthofinder", "tsv"])
    ap.add_argument("--id-normalizer", default="auto",
                    help="identity|first_token|uniprot|auto|map:<tsv>, '+'-combinable (default auto)")
    ap.add_argument("--species-map", help="TSV: species name in tool output -> species code in the tree")
    ap.add_argument("--level-map", help="TSV (orthoxml): level name in file -> level name in the tree")
    ap.add_argument("--of-species-tree", help="OrthoFinder SpeciesTree_rooted_node_labels.txt (auto-detected)")
    ap.add_argument("--orthoxml-id-attr", default="protId", help="gene attribute holding the ID (protId)")
    ap.add_argument("--out", required=True, help="output path (.json or .tsv; both are written)")
    ap.add_argument("--n-boot", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=SEEDS[0])
    ap.add_argument("--no-per-gene", action="store_true", help="omit per-gene rows from the JSON")
    ap.add_argument("--skip-validation", action="store_true")
    a = ap.parse_args(argv)

    st = SpeciesTree.from_newick(a.panel_tree)
    families = load_families(a.families, panel=a.panel)
    if not families:
        print(f"no families found in {a.families}" + (f" for panel {a.panel}" if a.panel else ""), file=sys.stderr)
        return 2
    if not a.skip_validation:
        rep = validate_panel(families, st)
        for wmsg in rep.warnings:
            print("WARNING:", wmsg, file=sys.stderr)
        if not rep.ok:
            print(rep, file=sys.stderr)
            return 1
    keep = set()
    for f in families:
        keep.update(f.member_ids)
        keep.update(f.decoys)
    kw: dict = {"normalize": make_normalizer(a.id_normalizer), "keep_genes": keep}
    if a.format == "orthoxml":
        kw.update(species_map=_read_map(a.species_map), level_map=_read_map(a.level_map),
                  id_attr=a.orthoxml_id_attr)
    elif a.format == "orthofinder":
        kw.update(species_map=_read_map(a.species_map), of_species_tree=a.of_species_tree)
    pred = read_prediction(a.format, a.prediction, st, **kw)
    res = evaluate_panel(st, families, pred, n_boot=a.n_boot, seed=a.seed, per_gene=not a.no_per_gene)
    res["meta"] = {"panel_tree": a.panel_tree, "families_dir": a.families, "panel": a.panel,
                   "prediction": a.prediction, "format": a.format, "id_normalizer": a.id_normalizer,
                   "n_boot": a.n_boot, "seed": a.seed}
    out = Path(a.out)
    stem = out.with_suffix("") if out.suffix in (".json", ".tsv") else out
    stem.parent.mkdir(parents=True, exist_ok=True)
    jpath, tpath = stem.with_suffix(".json"), stem.with_suffix(".tsv")
    jpath.write_text(json.dumps(res, indent=1, default=str))
    write_tsv(res, tpath)
    macro = res["aggregate"]["macro"]
    print(f"{len(families)} families evaluated -> {jpath}, {tpath}", file=sys.stderr)
    for k in ("merge_accuracy", "level_mean_f1", "ortholog_f1", "event_recovery_rate", "contamination"):
        m = macro[k]
        if m["mean"] is not None:
            print(f"  {k:24s} {m['mean']:.4f}  [{m['ci_low']:.4f}, {m['ci_high']:.4f}]  n={m['n_families']}",
                  file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
