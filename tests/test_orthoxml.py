"""OrthoXML parsing and pyHam-style level induction.

Toy species tree: (((HUMAN,MOUSE)Mammalia,CHICK)Amniota,(DANRE,ORYLA)Teleostei)Euteleostomi.
"""
from pathlib import Path

import pytest

from hobench.metrics import evaluate_family, panel_context
from hobench.predictions import PredictionError, read_orthoxml
from hobench.speciestree import SpeciesTree

from conftest import make_family

FIX = Path(__file__).parent / "fixtures"
SPECIES = ["HUMAN", "MOUSE", "CHICK", "DANRE", "ORYLA"]


def orthoxml(tmp_path, genes: list[str], groups: str, taxonomy: str = "", version="0.3") -> Path:
    """Write an OrthoXML file. Gene <gene id> = position, protId = decorated UniProt-style
    'sp|<gene>|<gene>_X' (normalised back to <gene> by the default 'auto' normaliser).
    In ``groups`` refer to genes as {GENE_NAME}."""
    ids = {g: str(i + 1) for i, g in enumerate(genes)}
    sp_xml = []
    for s in SPECIES:
        gs = [g for g in genes if g.split("_")[0] == s]
        if gs:
            sp_xml.append(f'<species name="{s}" NCBITaxId="0"><database name="d" version="1"><genes>'
                          + "".join(f'<gene id="{ids[g]}" protId="sp|{g}|{g}_X"/>' for g in gs)
                          + "</genes></database></species>")
    body = groups
    for g, i in ids.items():
        body = body.replace("{" + g + "}", f'<geneRef id="{i}"/>')
    xml = (f'<?xml version="1.0" encoding="utf-8"?>\n<orthoXML xmlns="http://orthoXML.org/2011/" '
           f'version="{version}" origin="test" originVersion="1">' + "".join(sp_xml) + taxonomy
           + f"<groups>{body}</groups></orthoXML>")
    p = tmp_path / "test.orthoxml"
    p.write_text(xml)
    return p


def tr(level: str) -> str:
    return f'<property name="TaxRange" value="{level}"/>'


def same(p, level, a, b) -> bool:
    ha, hb = p.hog(level, a), p.hog(level, b)
    return ha is not None and ha == hb


def test_flat_group_induces_all_levels(st, tmp_path, fams, fam_by_id):
    """A single group at Euteleostomi listing five genes directly (levels Amniota, Mammalia,
    Teleostei skipped). No paralogGroup -> no duplication: genes below the same child level
    are grouped, recursively: Amniota {H,M,C}, Mammalia {H,M}, Teleostei {D,O}."""
    genes = ["HUMAN_1", "MOUSE_1", "CHICK_1", "DANRE_a", "ORYLA_b"]
    p = orthoxml(tmp_path, genes, f"<orthologGroup id='g1'>{tr('Euteleostomi')}" +
                 "".join("{" + g + "}" for g in genes) + "</orthologGroup>")
    pred = read_orthoxml(p, st)
    assert set(pred.levels) == {"Euteleostomi", "Amniota", "Mammalia", "Teleostei"}
    assert len(set(pred.levels["Euteleostomi"].values())) == 1
    assert sorted(pred.levels["Amniota"]) == ["CHICK_1", "HUMAN_1", "MOUSE_1"]
    assert same(pred, "Amniota", "HUMAN_1", "CHICK_1") and same(pred, "Mammalia", "HUMAN_1", "MOUSE_1")
    assert same(pred, "Teleostei", "DANRE_a", "ORYLA_b")
    assert pred.check(st) == []
    assert pred.gene_species["DANRE_a"] == "DANRE"
    # on F1 this equals the species-tree-shaped prediction: 9/10 pairs (pseudo-orthologs wrong)
    r = evaluate_family(fam_by_id["F1"], st, pred, panel_context(fams, st))
    assert r["pair_merge"]["accuracy_pairs"] == pytest.approx(0.9)


def test_paralog_group_dup_at_enclosing_level(st, tmp_path):
    """Group at Amniota: CHICK_g + paralogGroup(HUMAN_x, MOUSE_y). The copies' levels are
    HUMAN and MOUSE, MRCA Mammalia, so the duplication is at parent(Mammalia) = Amniota =
    the enclosing group: the two copies are separate HOGs at Mammalia (each its own chain);
    all three genes are together at Amniota."""
    p = orthoxml(tmp_path, ["HUMAN_x", "MOUSE_y", "CHICK_g"],
                 f"<orthologGroup>{tr('Amniota')}{{CHICK_g}}<paralogGroup>{{HUMAN_x}}{{MOUSE_y}}"
                 f"</paralogGroup></orthologGroup>")
    pred = read_orthoxml(p, st)
    assert same(pred, "Amniota", "HUMAN_x", "MOUSE_y") and same(pred, "Amniota", "HUMAN_x", "CHICK_g")
    assert pred.hog("Mammalia", "HUMAN_x") is not None and pred.hog("Mammalia", "MOUSE_y") is not None
    assert not same(pred, "Mammalia", "HUMAN_x", "MOUSE_y")
    assert "Euteleostomi" not in pred.levels


def test_paralog_group_dup_placement_declared_vs_late(st, tmp_path):
    """Group declared at Euteleostomi: DANRE_g, MOUSE_z, CHICK_c and paralogGroup(HUMAN_x, HUMAN_y).

    declared (E1, primary): the duplication is just below the declared level Euteleostomi, so
    each copy has its own chain: HUMAN_x and HUMAN_y are apart at Amniota and Mammalia (merge
    Euteleostomi); MOUSE_z and CHICK_c (non-paralogs strictly below Amniota) share an Amniota HOG.
    late (pyHam): copies are both HUMAN -> duplication at parent(HUMAN) = Mammalia: an induced
    Mammalia HOG holds both; MOUSE_z gets its own Mammalia HOG (pyHam); Amniota groups all four.
    The readings differ -> dup_placement_ambiguous = 1."""
    p = orthoxml(tmp_path, ["HUMAN_x", "HUMAN_y", "MOUSE_z", "CHICK_c", "DANRE_g"],
                 f"<orthologGroup>{tr('Euteleostomi')}{{DANRE_g}}<paralogGroup>{{HUMAN_x}}{{HUMAN_y}}"
                 f"</paralogGroup>{{MOUSE_z}}{{CHICK_c}}</orthologGroup>")
    d = read_orthoxml(p, st)
    assert d.meta["dup_placement"] == "declared" and d.meta["dup_placement_ambiguous"] == 1
    assert not same(d, "Amniota", "HUMAN_x", "HUMAN_y") and not same(d, "Mammalia", "HUMAN_x", "HUMAN_y")
    assert same(d, "Amniota", "MOUSE_z", "CHICK_c") and not same(d, "Amniota", "HUMAN_x", "MOUSE_z")
    assert same(d, "Euteleostomi", "HUMAN_x", "HUMAN_y")
    late = read_orthoxml(p, st, dup_placement="late")
    assert same(late, "Mammalia", "HUMAN_x", "HUMAN_y") and not same(late, "Mammalia", "HUMAN_x", "MOUSE_z")
    assert same(late, "Amniota", "HUMAN_x", "MOUSE_z") and same(late, "Amniota", "HUMAN_y", "CHICK_c")
    assert late.sizes["Amniota"][late.hog("Amniota", "HUMAN_x")] == 4
    assert late.meta["dup_placement_ambiguous"] == 1


def test_clade_level_mapping_E2(tmp_path):
    """FastOMA renames short level names (KE -> internal_0). Levels are mapped by clade species
    set via the file's <taxonomy> (or the tool tree), so our 'KE' is found; an unknown TaxRange
    is a hard error."""
    our = SpeciesTree.from_newick("((YEAST,CANGA)PostWGD,(KLULA,ASHGO)KE)Sacch;")
    pred = read_orthoxml(FIX / "fastoma_excerpt.orthoxml", our)
    assert "clade sets" in pred.meta["level_mapping"]
    assert same(pred, "KE", "Q75A29", "Q6CL74") and same(pred, "PostWGD", "P10870", "Q6FNU3")
    bad = tmp_path / "bad.orthoxml"
    bad.write_text((FIX / "fastoma_excerpt.orthoxml").read_text().replace('value="internal_0"', 'value="Foo"')
                   .replace('taxonId="5"', ""))
    with pytest.raises(PredictionError, match="cannot be mapped"):
        read_orthoxml(bad, our)
    # explicit tool tree (species_tree_checked.nwk style) on a 0.3-style file without taxonomy
    tool = tmp_path / "tool.nwk"
    tool.write_text("((YEAST,CANGA)Saccharomycetaceae_WGD,(KLULA,ASHGO)internal_0)Saccharomycetaceae;")
    nt = tmp_path / "notax.orthoxml"
    txt = (FIX / "fastoma_excerpt.orthoxml").read_text()
    import re
    txt = re.sub(r"<taxonomy>.*?</taxonomy>", "", txt, flags=re.S)
    txt = re.sub(r' taxonId="\d+"', "", txt)
    nt.write_text(txt)
    p2 = read_orthoxml(nt, our, tool_species_tree=tool)
    assert same(p2, "KE", "Q75A29", "Q6CL74")
    with pytest.raises(PredictionError, match="cannot be mapped"):
        read_orthoxml(nt, our)   # names internal_0 / Saccharomycetaceae_WGD are not in our tree


def test_nested_groups_and_v05_taxonomy(st, tmp_path, fams, fam_by_id):
    """F2 encoded the way FastOMA writes it (OrthoXML 0.5: taxonId referring to <taxonomy>,
    no TaxRange property), with Teleostei skipped for DANRE_g. Must reproduce F2 exactly."""
    taxonomy = ('<taxonomy><taxon id="1" name="Euteleostomi"><taxon id="2" name="Amniota">'
                '<taxon id="3" name="Mammalia"><taxon id="4" name="HUMAN"/><taxon id="5" name="MOUSE"/></taxon>'
                '<taxon id="6" name="CHICK"/></taxon><taxon id="7" name="Teleostei"><taxon id="8" name="DANRE"/>'
                '<taxon id="9" name="ORYLA"/></taxon></taxon></taxonomy>')
    groups = ('<orthologGroup id="HOG:1_1" taxonId="1">{DANRE_g}<orthologGroup id="HOG:1_2" taxonId="2">'
              '{CHICK_g}<paralogGroup><orthologGroup id="HOG:1.1a_3" taxonId="3">{HUMAN_x}{MOUSE_x}'
              '</orthologGroup><orthologGroup id="HOG:1.1b_3" taxonId="3">{HUMAN_y}{MOUSE_y}</orthologGroup>'
              '</paralogGroup></orthologGroup></orthologGroup>')
    genes = ["HUMAN_x", "MOUSE_x", "HUMAN_y", "MOUSE_y", "CHICK_g", "DANRE_g"]
    pred = read_orthoxml(orthoxml(tmp_path, genes, groups, taxonomy, version="0.5"), st)
    r = evaluate_family(fam_by_id["F2"], st, pred, panel_context(fams, st))
    assert r["pair_merge"]["accuracy_pairs"] == 1.0 and r["events"][0]["recovered_0.8"]


def test_level_from_mrca_when_unannotated(st, tmp_path):
    """No TaxRange/taxonId: level = MRCA of children (pyHam). {HUMAN_1, MOUSE_1} -> Mammalia;
    a group with a single species child -> its parent level."""
    p = orthoxml(tmp_path, ["HUMAN_1", "MOUSE_1", "CHICK_1", "DANRE_1", "ORYLA_1"],
                 "<orthologGroup>{HUMAN_1}{MOUSE_1}</orthologGroup>"
                 "<orthologGroup>{CHICK_1}</orthologGroup>"
                 "<orthologGroup>{DANRE_1}{ORYLA_1}</orthologGroup>")
    pred = read_orthoxml(p, st)
    assert same(pred, "Mammalia", "HUMAN_1", "MOUSE_1") and pred.hog("Amniota", "HUMAN_1") is None
    assert pred.hog("Amniota", "CHICK_1") is not None   # single-species group -> parent(CHICK)
    assert same(pred, "Teleostei", "DANRE_1", "ORYLA_1")


def test_conflicts_and_unknown_species(st, tmp_path):
    p = orthoxml(tmp_path, ["HUMAN_1", "DANRE_1"],
                 f"<orthologGroup id='bad'>{tr('Mammalia')}{{HUMAN_1}}{{DANRE_1}}</orthologGroup>")
    with pytest.raises(PredictionError, match="declared level"):
        read_orthoxml(p, st)                       # strict by default: no silent repair
    pred = read_orthoxml(p, st, strict=False)      # sensitivity use: repaired + recorded
    assert same(pred, "Euteleostomi", "HUMAN_1", "DANRE_1") and pred.meta["warnings"]
    st2 = SpeciesTree.from_newick("((HUMAN,MOUSE)Mammalia,CHICK)Amniota;")
    with pytest.raises(PredictionError, match="species not in the species tree"):
        read_orthoxml(p, st2)
    with pytest.raises(PredictionError, match="does not refer to a gene"):
        read_orthoxml(p, st2, ignore_unknown_species=True)


def test_keep_genes_and_sizes(st, tmp_path):
    genes = ["HUMAN_1", "MOUSE_1", "CHICK_1"]
    p = orthoxml(tmp_path, genes, f"<orthologGroup>{tr('Amniota')}" + "".join("{" + g + "}" for g in genes)
                 + "</orthologGroup>")
    pred = read_orthoxml(p, st, keep_genes={"HUMAN_1"})
    assert list(pred.levels["Amniota"]) == ["HUMAN_1"]
    h = pred.hog("Amniota", "HUMAN_1")
    assert pred.hog_size("Amniota", h) == 3 and pred.hog_members("Amniota", h) == ["HUMAN_1"]


# ------------------------------------------------------------------ real FastOMA output
YEAST_TREE = "((YEAST,CANGA)Saccharomycetaceae_WGD,(KLULA,ASHGO)internal_0)Saccharomycetaceae;"


def snf3_family():
    """SNF3/RGT2 glucose sensors: ohnologs from the yeast WGD (on the branch leading to the
    post-WGD clade). Reference tree (UniProt accessions):
    root S( D[lvl=Saccharomycetaceae_WGD]( S(SNF3 P10870, CANGA Q6FNU3), S(RGT2 Q12300, CANGA Q6FQS7) ),
            S(ASHGO Q75A29, KLULA Q6CL74) ).
    Cross-ohnolog pairs merge at parent(WGD clade) = Saccharomycetaceae."""
    import pandas as pd
    from hobench.newick import parse_newick
    from hobench.reference import ReferenceFamily
    t = parse_newick("(((P10870,Q6FNU3)[&&NHX:ev=S],(Q12300,Q6FQS7)[&&NHX:ev=S])"
                     "[&&NHX:ev=D:lvl=Saccharomycetaceae_WGD:tag=WGD],(Q75A29,Q6CL74)[&&NHX:ev=S])[&&NHX:ev=S];")
    sp = {"P10870": "YEAST", "Q12300": "YEAST", "Q6FNU3": "CANGA", "Q6FQS7": "CANGA",
          "Q75A29": "ASHGO", "Q6CL74": "KLULA"}
    m = pd.DataFrame({"gene_id": list(sp), "species": list(sp.values())})
    return ReferenceFamily("SNF3", t, m, meta={"id": "SNF3", "name": "SNF3/RGT2", "panel": "yeast",
                                               "phenomena": ["WGD"], "sources": [{"id": "S1", "doi": "x"}]})


def test_fastoma_excerpt():
    st = SpeciesTree.from_newick(YEAST_TREE)
    pred = read_orthoxml(FIX / "fastoma_excerpt.orthoxml", st)
    assert pred.check(st) == [] and pred.meta["warnings"] == []
    # SNF3 group: ohnolog pairs separate at the WGD clade, together at the root
    assert same(pred, "Saccharomycetaceae_WGD", "P10870", "Q6FNU3")
    assert not same(pred, "Saccharomycetaceae_WGD", "P10870", "Q12300")
    assert same(pred, "Saccharomycetaceae", "P10870", "Q12300")
    assert same(pred, "internal_0", "Q75A29", "Q6CL74")
    # in-paralog group: paralogGroup of two KLULA genes directly in a group declared at the root
    # (with a YEAST gene). declared (E1): duplication just below the root -> the KLULA copies are
    # apart at internal_0; late (pyHam): duplication at parent(KLULA) = internal_0 -> together.
    k1, k2, y = "Q6CR14", "Q6CVW9", "P22580"   # KLULA, KLULA, YEAST (AMDY_YEAST)
    assert not same(pred, "internal_0", k1, k2) and same(pred, "Saccharomycetaceae", k1, k2)
    assert pred.meta["dup_placement_ambiguous"] == 1
    assert same(read_orthoxml(FIX / "fastoma_excerpt.orthoxml", st, dup_placement="late"), "internal_0", k1, k2)
    assert same(pred, "Saccharomycetaceae", k1, y) and pred.hog("internal_0", y) is None
    assert pred.hog("Saccharomycetaceae_WGD", y) is not None   # skipped level induced (singleton)
    fam = snf3_family()
    r = evaluate_family(fam, st, pred, panel_context([fam], st))
    assert r["pair_merge"]["accuracy_pairs"] == 1.0 and r["events"][0]["recovered_0.8"]


def test_gene_in_two_root_groups(st, tmp_path):
    p = orthoxml(tmp_path, ["HUMAN_1", "MOUSE_1", "CHICK_1"],
                 f"<orthologGroup>{tr('Mammalia')}{{HUMAN_1}}{{MOUSE_1}}</orthologGroup>"
                 f"<orthologGroup>{tr('Amniota')}{{HUMAN_1}}{{CHICK_1}}</orthologGroup>")
    with pytest.raises(PredictionError, match="two HOGs"):
        read_orthoxml(p, st)
    pred = read_orthoxml(p, st, strict=False)
    assert pred.meta["membership_conflicts"] == 1
    assert same(pred, "Mammalia", "HUMAN_1", "MOUSE_1")   # first assignment kept
