"""End-to-end CLI (python -m hobench.evaluate), panel FASTA checks (C4) and freezing (C9)."""
import csv
import hashlib
import json
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

from hobench import controls
from hobench.evaluate import main
from hobench.freeze import check_manifest, main as freeze_main, write_manifest
from hobench.panel import load_panel, validate_members
from hobench.reference import write_family

from conftest import TOY_TREE, make_family
from test_orthoxml import YEAST_TREE, snf3_family

FIX = Path(__file__).parent / "fixtures"
ROOT = Path(__file__).parent.parent


def make_panel(tmp_path: Path, fams, extra_ids=()) -> Path:
    """panels/TOY/species.tsv + one FASTA per species holding every member (UniProt-style
    headers) plus ``extra_ids``."""
    pdir = tmp_path / "panels" / "TOY"
    (pdir / "proteomes").mkdir(parents=True)
    by_sp: dict[str, list[str]] = {}
    for f in fams:
        for g, s in f.gene_species.items():
            by_sp.setdefault(s, []).append(g)
    for g in extra_ids:
        by_sp.setdefault(g.split("_")[0], []).append(g)
    rows = []
    for s in ["HUMAN", "MOUSE", "CHICK", "DANRE", "ORYLA"]:
        fa = pdir / "proteomes" / f"{s}.fa"
        fa.write_text("".join(f">sp|{g}|{g}_X some description\nMSEQ\n" for g in by_sp.get(s, [])))
        rows.append({"code": s, "upid": "UP0", "taxid": "1", "mnemonic": s,
                     "fasta_path": f"panels/TOY/proteomes/{s}.fa",
                     "sha256": hashlib.sha256(fa.read_bytes()).hexdigest()})
    pd.DataFrame(rows).to_csv(pdir / "species.tsv", sep="\t", index=False)
    return pdir


def setup_panel(tmp_path, fams, st):
    fdir = tmp_path / "families"
    for f in fams:
        write_family(fdir / f.id, f.meta, f.members, f.tree)
    (tmp_path / "panel.nwk").write_text(TOY_TREE)
    pred_path = tmp_path / "onebig.tsv"
    controls.one_big_hog(fams, st).write_tsv(pred_path)
    return fdir, pred_path, make_panel(tmp_path, fams)


def test_cli_tsv_prediction(st, fams, tmp_path):
    fdir, pred_path, pdir = setup_panel(tmp_path, fams, st)
    out = tmp_path / "out" / "results.json"
    args = ["--panel-tree", str(tmp_path / "panel.nwk"), "--families", str(fdir), "--prediction",
            str(pred_path), "--format", "tsv", "--out", str(out), "--n-boot", "200"]
    assert main(args) == 2                         # --panel-dir required unless explicitly skipped
    assert main(args + ["--panel-dir", str(pdir)]) == 0
    res = json.loads(out.read_text())
    assert [r["id"] for r in res["families"]] == ["F1", "F2", "F3", "F4"]
    comb = res["aggregate"]["combined"]
    assert comb["merge_accuracy_pairs"]["mean"] == pytest.approx((0.9 + 11 / 15 + 1 + 2 / 6) / 4)
    assert comb["merge_accuracy_pairs"]["n"] == 4 and comb["merge_accuracy_pairs"]["ci_low"] is None
    prov = res["meta"]["provenance"]
    assert prov["panel_tree_sha256"] == hashlib.sha256((tmp_path / "panel.nwk").read_bytes()).hexdigest()
    assert prov["prediction_sha256"] and prov["families_tree_hash"] and "hobench_git" in prov
    assert prov["families_manifest_present"] is False
    # B4: the one-big-HOG is itself a baseline -> skill <= 0 for every primary metric
    for r in res["families"]:
        assert all(v is None or v <= 1e-12 for v in r["skill"].values()), r["skill"]
    rows = list(csv.DictReader(open(out.with_suffix(".tsv")), delimiter="\t"))
    assert float(rows[0]["merge_accuracy_pairs"]) == pytest.approx(0.9)
    assert "skill_merge_accuracy_node" in rows[0]
    # explicit skip is allowed and recorded
    assert main(args + ["--no-fasta-check", "--no-skill"]) == 0
    assert json.loads(out.read_text())["meta"]["fasta_check"] is False


def test_cli_orthofinder_real_fixture_subprocess(tmp_path):
    fam = snf3_family()
    write_family(tmp_path / "families" / "SNF3", fam.meta, fam.members, fam.tree)
    (tmp_path / "yeast.nwk").write_text(YEAST_TREE)
    out = tmp_path / "of_results.tsv"
    cmd = [sys.executable, "-m", "hobench.evaluate", "--panel-tree", str(tmp_path / "yeast.nwk"),
           "--families", str(tmp_path / "families"), "--prediction",
           str(FIX / "orthofinder_of3" / "Results_fixture"), "--format", "orthofinder", "--out", str(out),
           "--n-boot", "100", "--no-per-gene", "--no-fasta-check"]
    p = subprocess.run(cmd, cwd=ROOT, capture_output=True, text=True)
    assert p.returncode == 0, p.stderr
    res = json.loads(out.with_suffix(".json").read_text())
    assert res["families"][0]["pair_merge"]["accuracy_pairs"] == 1.0
    assert "genes" not in res["families"][0]


def test_cli_rejects_invalid_reference(tmp_path, fams, st):
    bad = fams[0]
    bad.tree.features["ev"] = "Q"               # invalid event code -> validation error, exit 1
    write_family(tmp_path / "families" / "F1", bad.meta, bad.members, bad.tree)
    (tmp_path / "panel.nwk").write_text(TOY_TREE)
    (tmp_path / "p.tsv").write_text("level\thog_id\tgene\n")
    rc = main(["--panel-tree", str(tmp_path / "panel.nwk"), "--families", str(tmp_path / "families"),
               "--prediction", str(tmp_path / "p.tsv"), "--format", "tsv", "--out", str(tmp_path / "r.json"),
               "--no-fasta-check"])
    assert rc == 1


def test_validate_members_against_panel(tmp_path, fams):
    """C4: absence from the pinned FASTA and mapping-threshold violations are errors."""
    pdir = make_panel(tmp_path, fams)
    panel = load_panel(pdir)
    assert panel.verify_checksums() == []
    assert validate_members(fams[0], panel).ok
    f = make_family("M", "(HUMAN_1,CHICK_1)[&&NHX:ev=S];")
    f.members.loc[0, "gene_id"] = "HUMAN_404"
    assert any("not found in the pinned HUMAN FASTA" in e for e in validate_members(f, panel).errors)

    def with_mapping(method, ident, cov, sp="HUMAN", **meta):
        g = make_family("M", f"({sp}_1,CHICK_1)[&&NHX:ev=S];", **meta)
        g.members["mapping_method"] = method
        g.members["mapping_identity"] = ident
        g.members["mapping_coverage"] = cov
        return validate_members(g, panel)

    assert with_mapping("UniProtKB_AC-ID", "", "").ok                 # exact ID mapping: exempt
    assert with_mapping("diamond", "98.5", "0.9").ok
    assert any("identity 97.0% < 98.0%" in e for e in with_mapping("diamond", "97.0", "0.9").errors)
    assert any("< 99.0%" in e for e in with_mapping("diamond", "98.5", "0.9", phenomena=["TAND"]).errors)
    assert any("coverage 0.70 < 0.8" in e for e in with_mapping("mmseqs", "99.5", "0.7").errors)
    assert any("fraction" in e for e in with_mapping("diamond", "99.5", "85").errors)
    assert with_mapping("diamond", "99.5", "1.0/0.9").ok                # "query/target" form: min is checked
    assert any("coverage 0.70 < 0.8" in e for e in with_mapping("diamond", "99.5", "1.0/0.7").errors)
    assert any("not an accepted mapping" in e for e in with_mapping("unmapped", "", "").errors)
    # stricter species (XENLA) need >= 99% (panel species.tsv strict_mapping column also works)
    sp_tab = pd.read_csv(pdir / "species.tsv", sep="\t", dtype=str)
    sp_tab["strict_mapping"] = ["1", "0", "0", "0", "0"]
    sp_tab.to_csv(pdir / "species.tsv", sep="\t", index=False)
    assert any("< 99.0%" in e for e in validate_members(_fam_with("diamond", "98.5", "0.9"), load_panel(pdir)).errors)
    # checksum mismatch detected
    (pdir / "proteomes" / "HUMAN.fa").write_text(">sp|X|X\nM\n")
    assert any("sha256 mismatch for HUMAN" in e for e in load_panel(pdir).verify_checksums())


def _fam_with(method, ident, cov):
    g = make_family("M", "(HUMAN_1,CHICK_1)[&&NHX:ev=S];")
    g.members["mapping_method"] = method
    g.members["mapping_identity"] = ident
    g.members["mapping_coverage"] = cov
    return g


def test_freeze_manifest(tmp_path, fams, st):
    fdir, pred_path, pdir = setup_panel(tmp_path, fams, st)
    assert freeze_main(["--families", str(fdir), "--check"]) == 1      # no manifest yet
    assert freeze_main(["--families", str(fdir)]) == 0
    assert check_manifest(fdir) == []
    assert freeze_main(["--families", str(fdir), "--check"]) == 0
    out = tmp_path / "r.json"
    assert main(["--panel-tree", str(tmp_path / "panel.nwk"), "--families", str(fdir), "--prediction",
                 str(pred_path), "--format", "tsv", "--out", str(out), "--panel-dir", str(pdir),
                 "--n-boot", "100", "--no-skill"]) == 0
    prov = json.loads(out.read_text())["meta"]["provenance"]
    assert prov["families_manifest_present"] and prov["families_manifest_differences"] == []
    (fdir / "F1" / "family.yaml").write_text((fdir / "F1" / "family.yaml").read_text() + "\nnotes: edited\n")
    assert check_manifest(fdir) == ["changed: F1/family.yaml"]
