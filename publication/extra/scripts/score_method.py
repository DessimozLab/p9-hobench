#!/usr/bin/env python
"""Score your own hierarchical-orthology output on a HOBench panel, next to the paper's methods and baselines.

Run your method on the panel proteomes (panels/<P>/proteomes/<CODE>.fa for --mode full, panels/<P>/subset/<CODE>.fa
for --mode subset) with the panel species tree (panels/<P>/species_tree.nwk), then:

  python scripts/score_method.py --panel K --mode full --name mytool --format orthoxml --input my_hogs.orthoxml.gz
  python scripts/score_method.py --panel E --mode subset --name mytool --format orthofinder \
         --input Results_X/Phylogenetic_Hierarchical_Orthogroups --species-tree Results_X/Species_Tree/SpeciesTree_rooted_node_labels.txt
  python scripts/score_method.py --panel Y --mode full --name mytool --format tsv --input hogs.tsv   # level, hog_id, gene

Formats: OrthoXML 0.3/0.5 (TaxRange levels are matched to panel levels by species content; gene ids taken from the
`protId` attribute unless --id-attr is given), OrthoFinder 3 HOG directory, or a TSV with columns level, hog_id, gene
(level = a level id from panels/<P>/levels.tsv; add --flat for a single root-level partition). Gene ids may be UniProt
accessions or UniProt-style headers (sp|ACC|NAME).

Output: results_user/<name>/<P>_<mode>/summary.md (primary endpoints with trivial baselines and skill), per_family.tsv,
per_level.tsv, events.tsv. With --compare (default) the paper's OrthoFinder and FastOMA runs are scored alongside.
"""
from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--panel", required=True, choices=["Y", "K", "E", "V", "P"])
    ap.add_argument("--mode", required=True, choices=["full", "subset"])
    ap.add_argument("--name", required=True, help="method name used in result tables")
    ap.add_argument("--format", required=True, choices=["orthoxml", "orthofinder", "tsv"])
    ap.add_argument("--input", required=True, help="file (orthoxml, tsv) or directory (orthofinder)")
    ap.add_argument("--species-tree", default=None, help="tool's own species tree, if its level names differ")
    ap.add_argument("--id-attr", default="protId", help="OrthoXML gene attribute holding protein ids")
    ap.add_argument("--flat", action="store_true", help="TSV holds only root-level groups (no hierarchy)")
    ap.add_argument("--no-compare", action="store_true", help="do not score the paper's runs alongside")
    ap.add_argument("--n-boot", type=int, default=200)
    a = ap.parse_args()

    runs = ROOT / "results_user" / a.name / "runs" / a.panel / a.mode
    if runs.exists():
        shutil.rmtree(runs)
    run = runs / f"{a.name}_default"
    run.mkdir(parents=True)
    src = Path(a.input).resolve()
    dst = run / src.name
    dst.symlink_to(src)
    man = {"tool": "user", "name": a.name, "variant": "default", "status": "done", "format": a.format,
           "file": src.name, "id_attr": a.id_attr, "flat": a.flat, "panel": a.panel, "mode": a.mode,
           "registered": datetime.now(timezone.utc).isoformat()}
    if a.species_tree:
        t = Path(a.species_tree).resolve()
        (run / "tool_species_tree.nwk").symlink_to(t)
        man["species_tree"] = "tool_species_tree.nwk"
    (run / "manifest.json").write_text(json.dumps(man, indent=1))
    if not a.no_compare:
        paper = ROOT / "runs" / a.panel / a.mode
        for d in sorted(paper.glob("*")) if paper.exists() else []:
            if (d / "manifest.json").exists():
                (runs / d.name).symlink_to(d.resolve())
    out = ROOT / "results_user" / a.name / f"{a.panel}_{a.mode}"
    cmd = [sys.executable, str(ROOT / "scripts" / "evaluate_panel.py"), "--panel", a.panel, "--mode", a.mode,
           "--runs", str(runs), "--out", str(out), "--n-boot", str(a.n_boot)]
    print(" ".join(cmd), flush=True)
    rc = subprocess.call(cmd, cwd=ROOT)
    if rc == 0:
        print(f"\nresults: {out}/summary.md")
    return rc


if __name__ == "__main__":
    sys.exit(main())
