import pandas as pd
import pytest

from hobench.newick import write_newick
from hobench.reference import (ABOVE, INF, NOT_ASSERTED, ReferenceError, load_families, reconcile,
                               tree_from_classification, validate_family, validate_panel, write_family)

from conftest import make_family


def test_toy_families_valid(st, fams):
    for f in fams:
        rep = validate_family(f, st)
        assert rep.ok, str(rep)
    assert validate_panel(fams, st).ok


def _ev(rf, ev):
    return [e for e in rf.events if e.ev == ev]


def _merge(st, e):
    n = lambda x: "inf" if x == INF else st.names[x]
    return (n(e.merge_lo), n(e.merge_hi))


def test_reconcile_levels_and_merges(st, fam_by_id):
    # F1: exact lvl=Teleostei -> merge exactly parent(Teleostei) = Euteleostomi
    rf = reconcile(fam_by_id["F1"], st)
    assert st.names[rf.family_level] == "Euteleostomi"
    (d,) = _ev(rf, "D")
    assert st.names[d.level] == "Teleostei" and d.level_source == "exact" and d.tags == ["WGD"] and d.src == "S1"
    assert _merge(st, d) == ("Euteleostomi", "Euteleostomi") and not d.is_interval
    # F4: root D without lvl -> lvl in [Amniota (mapping), ABOVE_RANGE] -> merge [Euteleostomi, inf]
    rf4 = reconcile(fam_by_id["F4"], st)
    (d4,) = _ev(rf4, "D")
    assert d4.level_source == "default" and d4.lvl_hi == ABOVE and d4.is_interval
    assert _merge(st, d4) == ("Euteleostomi", "inf")
    assert st.names[rf4.family_level] == "Amniota"
    # root_single_copy: single gene at parent(Amniota) -> exact merge Euteleostomi
    nwk4 = write_newick(fam_by_id["F4"].tree, features=["ev"])
    d4b = _ev(reconcile(make_family("F4b", nwk4, root_single_copy=True), st), "D")[0]
    assert _merge(st, d4b) == ("Euteleostomi", "Euteleostomi")
    d4c = _ev(reconcile(make_family("F4c", nwk4, root_level="Euteleostomi", root_single_copy=True), st), "D")[0]
    assert _merge(st, d4c) == ("Euteleostomi", "Euteleostomi")
    # root_level alone only sets the scope: still [Euteleostomi, inf], family level Euteleostomi
    rf4d = reconcile(make_family("F4d", nwk4, root_level="Euteleostomi"), st)
    assert _merge(st, _ev(rf4d, "D")[0]) == ("Euteleostomi", "inf") and st.names[rf4d.family_level] == "Euteleostomi"
    # F3: U node -> not asserted
    assert [e.merge_lo for e in _ev(reconcile(fam_by_id["F3"], st), "U")] == [NOT_ASSERTED]


def test_default_bound_hidden_paralogy(st):
    """AC-2: copies retained only in HUMAN, sibling CHICK. The duplication lies between HUMAN
    (mapping) and the lineage below the speciation with CHICK: lvl in [HUMAN, Mammalia], so the
    merge interval is [parent(HUMAN), parent(Mammalia)] = [Mammalia, Amniota]."""
    rf = reconcile(make_family("HP", "((HUMAN_a,HUMAN_b)[&&NHX:ev=D],CHICK_c)[&&NHX:ev=S];"), st)
    (d,) = _ev(rf, "D")
    assert (st.names[d.lvl_lo], st.names[d.lvl_hi]) == ("HUMAN", "Mammalia")
    assert _merge(st, d) == ("Mammalia", "Amniota")
    # binary S consistent with the tree: closest sibling gives the same bound as the S mapping
    rf2 = reconcile(make_family("HP2", "(((HUMAN_a,HUMAN_b)[&&NHX:ev=D],MOUSE_c)[&&NHX:ev=S],CHICK_d)[&&NHX:ev=S];"), st)
    assert _merge(st, _ev(rf2, "D")[0]) == ("Mammalia", "Mammalia")
    # multifurcating S {HUMAN},{MOUSE},{CHICK}: closest sibling of HUMAN is MOUSE -> bound below Mammalia
    rf3 = reconcile(make_family("HP3", "((HUMAN_a,HUMAN_b)[&&NHX:ev=D],MOUSE_c,CHICK_d)[&&NHX:ev=S];"), st)
    assert _merge(st, _ev(rf3, "D")[0]) == ("Mammalia", "Mammalia")
    # nested D without lvl under a D with interval: child bound = parent's lvl_hi
    rf4 = reconcile(make_family("HP4", "(((HUMAN_a,HUMAN_b)[&&NHX:ev=D],MOUSE_c)[&&NHX:ev=D:lvl_lo=Mammalia:lvl_hi=Amniota],"
                                       "DANRE_d)[&&NHX:ev=S];"), st)
    outer, inner = _ev(rf4, "D")
    assert _merge(st, outer) == ("Amniota", "Euteleostomi")
    assert (st.names[inner.lvl_lo], st.names[inner.lvl_hi]) == ("HUMAN", "Amniota")


def test_explicit_interval_and_above_range(st):
    rf = reconcile(make_family("IV", "((HUMAN_a,MOUSE_a)[&&NHX:ev=S],(HUMAN_b,MOUSE_b)[&&NHX:ev=S])"
                                     "[&&NHX:ev=D:lvl_lo=Mammalia:lvl_hi=ABOVE_RANGE];"), st)
    assert _merge(st, _ev(rf, "D")[0]) == ("Amniota", "inf")


def test_yeast_allopolyploid_expressible():
    """A2: ohnologs apart at PostWGD and below, merge in [ZTWGD, Saccharomycetaceae];
    PostWGD vs KLE orthologs at the root; ZYGRO vs ohnologs unasserted (U)."""
    from hobench.speciestree import SpeciesTree
    yst = SpeciesTree.from_newick("(((YEAST,CANGA)PostWGD,ZYGRO)ZTWGD,(KLULA,ASHGO)KLE)Saccharomycetaceae;")
    nwk = ("((((YEAST_1,CANGA_1)[&&NHX:ev=S],(YEAST_2,CANGA_2)[&&NHX:ev=S])"
           "[&&NHX:ev=D:lvl_lo=PostWGD:lvl_hi=ZTWGD:tag=WGD,ALLO],ZYGRO_1)[&&NHX:ev=U],"
           "(KLULA_1,ASHGO_1)[&&NHX:ev=S])[&&NHX:ev=S];")
    f = make_family("YWGD", nwk)
    assert validate_family(f, yst).ok, str(validate_family(f, yst))
    rf = reconcile(f, yst)
    (d,) = _ev(rf, "D")
    assert _merge(yst, d) == ("ZTWGD", "Saccharomycetaceae") and d.tags == ["WGD", "ALLO"]
    # lvl_hi above the speciation bound (the root S) is an error
    bad = make_family("YB", nwk.replace("lvl_hi=ZTWGD", "lvl_hi=Saccharomycetaceae"))
    assert any("strictly below" in e for e in validate_family(bad, yst).errors)


def test_transfer_node(st):
    """A4: H node; recipient clade = DANRE_x (xenolog from the Amniota lineage... donor CHICK).
    Merge of xenolog vs native teleost genes = lca(CHICK, DANRE) = Euteleostomi."""
    nwk = "((DANRE_x[&&NHX:recipient=1],(DANRE_n,ORYLA_n)[&&NHX:ev=S])[&&NHX:ev=H:donor=CHICK:tag=HGT],CHICK_c)[&&NHX:ev=S];"
    f = make_family("HGT", nwk)
    rep = validate_family(f, st)
    assert rep.ok, str(rep)
    (h,) = _ev(reconcile(f, st), "H")
    assert _merge(st, h) == ("Euteleostomi", "Euteleostomi") and h.recipient_child == 0
    f2 = make_family("HGT2", nwk.replace("donor=CHICK", "donor=OUTSIDE").replace("CHICK_c", "HUMAN_c"))
    assert _merge(st, _ev(reconcile(f2, st), "H")[0]) == ("inf", "inf")
    for bad, msg in [(nwk.replace("[&&NHX:recipient=1]", ""), "exactly one child"),
                     (nwk.replace(":donor=CHICK", ""), "donor=<level|OUTSIDE> is required"),
                     (nwk.replace("donor=CHICK", "donor=Teleostei"), "nested")]:
        assert any(msg in e for e in validate_family(make_family("B", bad), st).errors), msg


@pytest.mark.parametrize("nwk,msg", [
    # two human genes whose LCA is a speciation node
    ("((HUMAN_1,HUMAN_2)[&&NHX:ev=S],CHICK_1)[&&NHX:ev=S];", "both from species 'HUMAN'"),
    # duplication level below the LCA mapping of its genes (Amniota)
    ("(HUMAN_1,CHICK_1)[&&NHX:ev=D:lvl=Mammalia];", "must be equal or ancestral"),
    # S children overlap once the explicit duplication level is taken into account
    ("((HUMAN_a,HUMAN_b)[&&NHX:ev=D:lvl=Amniota],CHICK_c)[&&NHX:ev=S];", "overlap in the species tree"),
    # nested duplication older than the duplication above it
    ("(((HUMAN_a,HUMAN_b)[&&NHX:ev=D:lvl=Amniota],MOUSE_c)[&&NHX:ev=D:lvl=Mammalia],CHICK_d)[&&NHX:ev=S];",
     "constraint from D ancestor"),
    ("((HUMAN_a,HUMAN_b)[&&NHX:ev=D:lvl=HUMAN:lvl_hi=Mammalia],CHICK_c)[&&NHX:ev=S];", "not both"),
    ("((HUMAN_a,HUMAN_b)[&&NHX:ev=D:lvl_lo=Mammalia:lvl_hi=HUMAN],CHICK_c)[&&NHX:ev=S];", "must be equal or ancestral to lvl_lo"),
    ("((HUMAN_1,MOUSE_1)[&&NHX:ev=S:lvl_hi=Amniota],CHICK_1)[&&NHX:ev=S];", "only allowed on D/H"),
    ("((HUMAN_1,MOUSE_1),CHICK_1)[&&NHX:ev=S];", "has no 'ev'"),
    ("((HUMAN_1,MOUSE_1)[&&NHX:ev=X],CHICK_1)[&&NHX:ev=S];", "must be S, D, U or H"),
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
