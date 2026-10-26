import pytest

from hobench.newick import parse_newick, write_newick
from hobench.speciestree import SpeciesTree, SpeciesTreeError


def test_newick_nhx_roundtrip():
    t = parse_newick("((a:1,'b c':2)X[&&NHX:ev=S],d)R;")
    assert t.name == "R" and [c.name for c in t.children] == ["X", "d"]
    assert t.children[0].features == {"ev": "S"}
    assert t.leaf_names() == ["a", "b c", "d"]
    assert t.children[0].children[0].length == 1.0
    s = write_newick(t, features=["ev"])
    assert s == "((a:1,'b c':2)X[&&NHX:ev=S],d)R;"
    assert write_newick(parse_newick(s), features=["ev"]) == s


def test_basic_relations(st):
    assert st.lca("HUMAN", "MOUSE") == "Mammalia"
    assert st.lca("HUMAN", "CHICK", "MOUSE") == "Amniota"
    assert st.lca("HUMAN", "DANRE") == "Euteleostomi"
    assert st.lca("HUMAN") == "HUMAN"
    assert st.is_ancestor("Amniota", "HUMAN")
    assert st.is_ancestor("Amniota", "Amniota")
    assert not st.is_ancestor("Amniota", "Amniota", strict=True)
    assert not st.is_ancestor("Mammalia", "CHICK")
    assert st.parent("Mammalia") == "Amniota" and st.parent("Euteleostomi") is None
    assert st.depth_of("HUMAN") == 3 and st.depth_of("Euteleostomi") == 0
    assert st.path("HUMAN", "DANRE") == ["HUMAN", "Mammalia", "Amniota", "Euteleostomi", "Teleostei", "DANRE"]
    assert st.distance("HUMAN", "DANRE") == 5
    assert st.distance("Mammalia", "Euteleostomi") == 2
    assert st.clade_species("Amniota") == {"HUMAN", "MOUSE", "CHICK"}
    assert st.find_clade({"DANRE", "ORYLA"}) == "Teleostei"
    assert st.find_clade({"HUMAN", "CHICK"}) is None
    assert set(st.species) == {"HUMAN", "MOUSE", "CHICK", "DANRE", "ORYLA"}
    assert set(st.levels) == {"Euteleostomi", "Amniota", "Mammalia", "Teleostei"}
    assert st.ancestors("CHICK") == ["CHICK", "Amniota", "Euteleostomi"]


def test_lca_matrix_matches_bruteforce(st):
    for a in st.names:
        for b in st.names:
            anc_a = set(st.ancestors(a))
            common = [x for x in st.ancestors(b) if x in anc_a]
            assert st.names[st.lca_idx[st.idx(a), st.idx(b)]] == common[0]


def test_prune(st):
    assert write_newick(st.prune_to(["HUMAN", "CHICK", "DANRE"])) == "((HUMAN,CHICK)Amniota,DANRE)Euteleostomi;"


@pytest.mark.parametrize("nwk,msg", [
    ("((A,B)X,(C,D)X)R;", "not unique"),
    ("((A,B),(C,D)Y)R;", "unnamed"),
    ("(((A,B)Z)X,C)R;", "unary"),
])
def test_invalid_trees(nwk, msg):
    with pytest.raises(SpeciesTreeError, match=msg):
        SpeciesTree.from_newick(nwk)
