import pandas as pd
import pytest

from hobench.newick import write_newick
from hobench.reference import (INF, NOT_ASSERTED, ReferenceError, load_families, reconcile,
                               tree_from_classification, validate_family, validate_panel, write_family)

from conftest import make_family


def test_toy_families_valid(st, fams):
    for f in fams:
        rep = validate_family(f, st)
        assert rep.ok, str(rep)
    assert validate_panel(fams, st).ok


def test_reconcile_levels_and_merges(st, fam_by_id):
    # F1: the D node at Teleostei merges at parent(Teleostei) = Euteleostomi (= family level)
    rf = reconcile(fam_by_id["F1"], st)
    assert st.names[rf.family_level] == "Euteleostomi"
    d = [e for e in rf.events if e.ev == "D"][0]
    assert st.names[d.level] == "Teleostei" and d.level_explicit and d.tags == ["WGD"] and d.src == "S1"
    assert st.names[d.merge] == "Euteleostomi"
    # F4: D without lvl -> level = LCA mapping = Amniota = family level -> merge = INF
    rf4 = reconcile(fam_by_id["F4"], st)
    d4 = [e for e in rf4.events if e.ev == "D"][0]
    assert st.names[d4.level] == "Amniota" and not d4.level_explicit and d4.merge == INF
    # the same family with root_level=Euteleostomi: merge becomes finite (Euteleostomi)
    f4b = make_family("F4b", write_newick(fam_by_id["F4"].tree, features=["ev"]), root_level="Euteleostomi")
    d4b = [e for e in reconcile(f4b, st).events if e.ev == "D"][0]
    assert st.names[d4b.merge] == "Euteleostomi"
    # F3: U node -> not asserted
    rf3 = reconcile(fam_by_id["F3"], st)
    assert [e.merge for e in rf3.events if e.ev == "U"] == [NOT_ASSERTED]


@pytest.mark.parametrize("nwk,msg", [
    # two human genes whose LCA is a speciation node
    ("((HUMAN_1,HUMAN_2)[&&NHX:ev=S],CHICK_1)[&&NHX:ev=S];", "both from species 'HUMAN'"),
    # duplication level below the LCA mapping of its genes (Amniota)
    ("(HUMAN_1,CHICK_1)[&&NHX:ev=D:lvl=Mammalia];", "must be equal or ancestral"),
    # S children overlap once the explicit duplication level is taken into account
    ("((HUMAN_a,HUMAN_b)[&&NHX:ev=D:lvl=Amniota],CHICK_c)[&&NHX:ev=S];", "overlap in the species tree"),
    # nested duplication older than the duplication above it
    ("(((HUMAN_a,HUMAN_b)[&&NHX:ev=D:lvl=Amniota],MOUSE_c)[&&NHX:ev=D:lvl=Mammalia],CHICK_d)[&&NHX:ev=S];",
     "constraint from an ancestor"),
    ("((HUMAN_1,MOUSE_1),CHICK_1)[&&NHX:ev=S];", "has no 'ev'"),
    ("((HUMAN_1,MOUSE_1)[&&NHX:ev=X],CHICK_1)[&&NHX:ev=S];", "must be S, D or U"),
    ("(HUMAN_1,CHICK_1)[&&NHX:ev=D:lvl=Primates];", "not a node of the species tree"),
    ("(HUMAN_1,CHICK_1)[&&NHX:ev=S:src=S9];", "src='S9'"),
    ("((HUMAN_1)[&&NHX:ev=S],CHICK_1)[&&NHX:ev=S];", "single child"),
])
def test_validation_errors(st, nwk, msg):
    rep = validate_family(make_family("BAD", nwk), st)
    assert not rep.ok
    assert any(msg in e for e in rep.errors), rep.errors
    with pytest.raises(ReferenceError):
        reconcile(make_family("BAD", nwk), st)


def test_member_table_errors(st):
    f = make_family("BAD", "(HUMAN_1,CHICK_1)[&&NHX:ev=S];")
    f.members.loc[0, "species"] = "PANTR"
    assert any("not in the panel species tree" in e for e in validate_family(f, st).errors)
    f = make_family("BAD", "(HUMAN_1,CHICK_1)[&&NHX:ev=S];")
    f.members = f.members[f.members.gene_id != "CHICK_1"]
    assert any("not listed in members" in e for e in validate_family(f, st).errors)
    f = make_family("BAD", "(HUMAN_1,HUMAN_1)[&&NHX:ev=D];")
    assert any("duplicated" in e for e in validate_family(f, st).errors)
    f = make_family("BAD", "(HUMAN_1,CHICK_1)[&&NHX:ev=S];", sources=[{"id": "S1"}])
    assert any("needs a 'doi' or 'pmid'" in e for e in validate_family(f, st).errors)


def test_unplaced_member_warning_and_panel_uniqueness(st, fams):
    f = make_family("W", "(HUMAN_w,CHICK_w)[&&NHX:ev=S];")
    f.members = pd.concat([f.members, pd.DataFrame([{"gene_id": "MOUSE_w", "species": "MOUSE"}])],
                          ignore_index=True)
    rep = validate_family(f, st)
    assert rep.ok and any("not placed" in w for w in rep.warnings)
    dup = make_family("F9", "(HUMAN_1,CHICK_9)[&&NHX:ev=S];")  # HUMAN_1 already in F1
    rep = validate_panel(fams + [dup], st)
    assert any("member of both 'F1' and 'F9'" in e for e in rep.errors)


def test_write_load_roundtrip(st, fams, tmp_path):
    for f in fams:
        write_family(tmp_path / f.id, f.meta, f.members, f.tree)
    loaded = load_families(tmp_path)
    assert [f.id for f in loaded] == ["F1", "F2", "F3", "F4"]
    assert validate_panel(loaded, st).ok
    for a, b in zip(fams, loaded):
        assert write_newick(a.tree, features=["ev", "lvl", "tag", "src"]) == \
               write_newick(b.tree, features=["ev", "lvl", "tag", "src"])
        assert list(a.members.gene_id) == list(b.members.gene_id)
    assert load_families(tmp_path, panel="other") == []


def test_tree_from_classification(st):
    table = pd.DataFrame([
        ("H_A1", "HUMAN", "PG1|A"), ("M_A1", "MOUSE", "PG1|A"), ("C_A1", "CHICK", "PG1|A"),
        ("H_A1b", "HUMAN", "PG1|A"),
        ("H_B1", "HUMAN", "PG1|B"), ("C_B1", "CHICK", "PG1|B"),
        ("H_A2", "HUMAN", "PG2|A"), ("M_A2", "MOUSE", "PG2|A"),
    ], columns=["gene", "species", "path"])
    groups = pd.DataFrame([{"group": "", "ev": "D", "lvl": "Amniota", "tag": ""},
                           {"group": "PG1", "ev": "D", "lvl": "", "tag": "WGD"}])
    t = tree_from_classification(table, st, groups)
    # Expected: root D(Amniota) splitting PG1 and PG2; PG1 = D(WGD) splitting A and B; each
    # terminal group follows the species tree; the two human genes of PG1|A form a U node
    # (in-paralogs, no assertion between them); PG2 has a single sub-group (A) -> collapsed.
    assert write_newick(t, features=["ev", "lvl", "tag"]) == (
        "(((((H_A1,H_A1b)[&&NHX:ev=U],M_A1)[&&NHX:ev=S],C_A1)A[&&NHX:ev=S],(H_B1,C_B1)B[&&NHX:ev=S])"
        "PG1[&&NHX:ev=D:tag=WGD],(H_A2,M_A2)A[&&NHX:ev=S])[&&NHX:ev=D:lvl=Amniota];")
    # same tree from inline per-component annotations
    table2 = table.assign(events="D|D", levels="Amniota|", tags="|WGD")
    table2.loc[table2.path.str.startswith("PG2"), ["events", "levels", "tags"]] = ["D", "Amniota", ""]
    t2 = tree_from_classification(table2, st)
    assert write_newick(t2, features=["ev", "lvl", "tag"]) == write_newick(t, features=["ev", "lvl", "tag"])
    members = table.rename(columns={"gene": "gene_id"})
    fam = make_family("CLS", write_newick(t, features=["ev", "lvl", "tag"]), root_level="Euteleostomi")
    fam.members = members
    assert validate_family(fam, st).ok
    # genes classified only to a group that has sub-groups require ev=U on that group
    bad = pd.concat([table, pd.DataFrame([("C_x", "CHICK", "PG1")], columns=table.columns)])
    with pytest.raises(ReferenceError, match="classified only down to it"):
        tree_from_classification(bad, st, groups)
    # inconsistent inline annotations
    table3 = table.assign(events="D")
    table3.loc[0, "events"] = "S"
    with pytest.raises(ReferenceError, match="inconsistent ev"):
        tree_from_classification(table3, st)
