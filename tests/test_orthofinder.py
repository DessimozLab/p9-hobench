"""OrthoFinder Phylogenetic_Hierarchical_Orthogroups parsing."""
import shutil
from pathlib import Path

import pytest

from hobench.metrics import evaluate_family, panel_context
from hobench.predictions import PredictionError, read_orthofinder
from hobench.speciestree import SpeciesTree

from test_orthoxml import YEAST_TREE, snf3_family

FIX = Path(__file__).parent / "fixtures"
HEADER = ["HOG", "OG", "Gene Tree Parent Clade", "CHICK", "DANRE", "HUMAN", "MOUSE", "ORYLA"]


def write_of(tmp_path: Path, of_tree: str, files: dict[str, list[list[str]]]) -> Path:
    res = tmp_path / "Results_X"
    (res / "Species_Tree").mkdir(parents=True)
    (res / "Species_Tree" / "SpeciesTree_rooted_node_labels.txt").write_text(of_tree)
    d = res / "Phylogenetic_Hierarchical_Orthogroups"
    d.mkdir()
    for name, rows in files.items():
        (d / f"{name}.tsv").write_text("\n".join("\t".join(r) for r in [HEADER] + rows) + "\n")
    return res


def dec(*genes):  # OrthoFinder keeps the FASTA header's first token
    return ", ".join(f"sp|{g}|{g}_X" for g in genes)


# F2 as OrthoFinder 2.5 would report it (N0 file present). OrthoFinder node labels:
# N0 = Euteleostomi, N1 = Amniota, N2 = Teleostei, N3 = Mammalia.
OF_TREE = "(((HUMAN:1,MOUSE:1)N3:1,CHICK:1)N1:1,(DANRE:1,ORYLA:1)N2:1)N0;"
F2_FILES = {
    "N0": [["N0.HOG0000000", "OG0000000", "n0", dec("CHICK_g"), dec("DANRE_g"), dec("HUMAN_x", "HUMAN_y"),
            dec("MOUSE_x", "MOUSE_y"), ""]],
    "N1": [["N1.HOG0000000", "OG0000000", "n1", dec("CHICK_g"), "", dec("HUMAN_x", "HUMAN_y"),
            dec("MOUSE_x", "MOUSE_y"), ""]],
    "N2": [["N2.HOG0000000", "OG0000000", "n2", "", dec("DANRE_g"), "", "", ""]],
    "N3": [["N3.HOG0000000", "OG0000000", "n3", "", "", dec("HUMAN_x"), dec("MOUSE_x"), ""],
           ["N3.HOG0000001", "OG0000000", "n4", "", "", dec("HUMAN_y"), dec("MOUSE_y"), ""]],
}


def test_of25_style_with_n0(st, tmp_path, fams, fam_by_id):
    res = write_of(tmp_path, OF_TREE, F2_FILES)
    pred = read_orthofinder(res, st)
    assert pred.meta["of_to_ours"] == {"N0": "Euteleostomi", "N1": "Amniota", "N2": "Teleostei", "N3": "Mammalia"}
    assert pred.hog("Mammalia", "HUMAN_x") == "N3.HOG0000000"      # decorated ID normalised
    assert pred.hog("Mammalia", "HUMAN_y") == "N3.HOG0000001"
    assert pred.gene_species["MOUSE_y"] == "MOUSE"
    r = evaluate_family(fam_by_id["F2"], st, pred, panel_context(fams))
    assert r["pair_merge"]["accuracy"] == 1.0 and r["events"][0]["recovered"]
    # path may also point at the HOG directory or one N*.tsv file
    p2 = read_orthofinder(res / "Phylogenetic_Hierarchical_Orthogroups" / "N1.tsv", st)
    assert p2.levels == pred.levels


def test_topology_mismatch_induces_level(st, tmp_path):
    """OrthoFinder's tree groups MOUSE with CHICK, so our Mammalia {HUMAN, MOUSE} has no
    matching OrthoFinder clade; it is induced from the smallest OrthoFinder clade containing it
    (N1 = {HUMAN, MOUSE, CHICK} = our Amniota), restricted to HUMAN and MOUSE genes."""
    of_tree = "((HUMAN:1,(MOUSE:1,CHICK:1)N3:1)N1:1,(DANRE:1,ORYLA:1)N2:1)N0;"
    files = {"N0": F2_FILES["N0"], "N1": F2_FILES["N1"], "N2": F2_FILES["N2"],
             "N3": [["N3.HOG0000000", "OG0000000", "n3", dec("CHICK_g"), "", "", dec("MOUSE_x", "MOUSE_y"), ""]]}
    pred = read_orthofinder(write_of(tmp_path, of_tree, files), st)
    assert pred.meta["unmatched_of_nodes"] == ["N3"]
    assert pred.meta["induced_levels"] == {"Mammalia": "N1"}
    assert pred.hog("Mammalia", "HUMAN_x") == pred.hog("Mammalia", "MOUSE_y") == "Mammalia<N1.HOG0000000"
    assert pred.hog("Mammalia", "CHICK_g") is None
    assert pred.hog_size("Mammalia", "Mammalia<N1.HOG0000000") == 4


def test_of3_real_fixture_root_from_orthogroups():
    """Excerpt of a real OrthoFinder 3.1.5 run (4 yeasts): there is no N0.tsv; root-level
    groups come from Orthogroups/Orthogroups.tsv. SNF3/RGT2 (WGD ohnologs) must come out exactly
    as the reference: separate at the post-WGD clade (N1), together at the root."""
    st = SpeciesTree.from_newick(YEAST_TREE)
    pred = read_orthofinder(FIX / "orthofinder_of3" / "Results_fixture", st)
    assert pred.meta["of_to_ours"] == {"N1": "Saccharomycetaceae_WGD", "N2": "internal_0", "N0": "Saccharomycetaceae"}
    assert pred.meta["root_from_orthogroups"].endswith("Orthogroups.tsv")
    assert pred.meta["n_hogs_not_nested_in_og"] == 0 and pred.meta["induced_levels"] == {}
    assert pred.hog("Saccharomycetaceae", "P10870") == "OG0000077"
    assert pred.hog("Saccharomycetaceae_WGD", "P10870") == "N1.HOG0000101"
    assert pred.hog("Saccharomycetaceae_WGD", "Q12300") == "N1.HOG0000102"
    assert pred.check(st) == []
    fam = snf3_family()
    r = evaluate_family(fam, st, pred, panel_context([fam]))
    assert r["pair_merge"]["accuracy"] == 1.0 and r["events"][0]["recovered"]


def test_of3_nesting_violation_detected(tmp_path):
    st = SpeciesTree.from_newick(YEAST_TREE)
    res = tmp_path / "R"
    shutil.copytree(FIX / "orthofinder_of3" / "Results_fixture", res)
    n1 = res / "Phylogenetic_Hierarchical_Orthogroups" / "N1.tsv"
    n1.write_text(n1.read_text().replace("N1.HOG0000835\tOG0000665", "N1.HOG0000835\tOG0000077"))
    pred = read_orthofinder(res, st)
    assert pred.meta["n_hogs_not_nested_in_og"] == 1
    (res / "Orthogroups" / "Orthogroups.tsv").unlink()
    with pytest.raises(PredictionError, match="root level"):
        read_orthofinder(res, st)


def test_species_mismatch_is_an_error(tmp_path):
    st = SpeciesTree.from_newick("((HUMAN,MOUSE)Mammalia,CHICK)Amniota;")
    with pytest.raises(PredictionError, match="not in our species tree"):
        read_orthofinder(write_of(tmp_path, OF_TREE, F2_FILES), st)
