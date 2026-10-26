"""End-to-end CLI: python -m hobench.evaluate."""
import csv
import json
import subprocess
import sys
from pathlib import Path

import pytest

from hobench import controls
from hobench.evaluate import main
from hobench.reference import write_family

from conftest import TOY_TREE
from test_orthoxml import YEAST_TREE, snf3_family

FIX = Path(__file__).parent / "fixtures"
ROOT = Path(__file__).parent.parent


def test_cli_tsv_prediction(st, fams, tmp_path):
    fdir = tmp_path / "families"
    for f in fams:
        write_family(fdir / f.id, f.meta, f.members, f.tree)
    (tmp_path / "panel.nwk").write_text(TOY_TREE)
    pred_path = tmp_path / "onebig.tsv"
    controls.one_big_hog(fams, st).write_tsv(pred_path)
    out = tmp_path / "out" / "results.json"
    rc = main(["--panel-tree", str(tmp_path / "panel.nwk"), "--families", str(fdir), "--prediction",
               str(pred_path), "--format", "tsv", "--out", str(out), "--n-boot", "200"])
    assert rc == 0
    res = json.loads(out.read_text())
    assert [r["id"] for r in res["families"]] == ["F1", "F2", "F3", "F4"]
    assert res["aggregate"]["macro"]["merge_accuracy"]["mean"] == pytest.approx((0.9 + 11 / 15 + 1 + 2 / 6) / 4)
    rows = list(csv.DictReader(open(out.with_suffix(".tsv")), delimiter="\t"))
    assert float(rows[0]["merge_accuracy"]) == pytest.approx(0.9)
    assert rows[3]["family"] == "F4" and float(rows[3]["split_excess"]) == -1


def test_cli_orthofinder_real_fixture_subprocess(tmp_path):
    fam = snf3_family()
    write_family(tmp_path / "families" / "SNF3", fam.meta, fam.members, fam.tree)
    (tmp_path / "yeast.nwk").write_text(YEAST_TREE)
    out = tmp_path / "of_results.tsv"
    cmd = [sys.executable, "-m", "hobench.evaluate", "--panel-tree", str(tmp_path / "yeast.nwk"),
           "--families", str(tmp_path / "families"), "--prediction",
           str(FIX / "orthofinder_of3" / "Results_fixture"), "--format", "orthofinder", "--out", str(out),
           "--n-boot", "100", "--no-per-gene"]
    p = subprocess.run(cmd, cwd=ROOT, capture_output=True, text=True)
    assert p.returncode == 0, p.stderr
    res = json.loads(out.with_suffix(".json").read_text())
    assert res["families"][0]["pair_merge"]["accuracy"] == 1.0
    assert "genes" not in res["families"][0]


def test_cli_rejects_invalid_reference(tmp_path, fams):
    bad = fams[0]
    bad.tree.features["ev"] = "Q"               # invalid event code -> validation error, exit 1
    write_family(tmp_path / "families" / "F1", bad.meta, bad.members, bad.tree)
    (tmp_path / "panel.nwk").write_text(TOY_TREE)
    (tmp_path / "p.tsv").write_text("level\thog_id\tgene\n")
    rc = main(["--panel-tree", str(tmp_path / "panel.nwk"), "--families", str(tmp_path / "families"),
               "--prediction", str(tmp_path / "p.tsv"), "--format", "tsv", "--out", str(tmp_path / "r.json")])
    assert rc == 1
