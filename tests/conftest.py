"""Shared toy data with hand-derivable answers.

Toy species tree (depth from the root in brackets)::

    Euteleostomi [0]
    ├── Amniota [1]
    │   ├── Mammalia [2]
    │   │   ├── HUMAN [3]
    │   │   └── MOUSE [3]
    │   └── CHICK [2]
    └── Teleostei [1]
        ├── DANRE [2]
        └── ORYLA [2]

Toy reference families (gene ID prefix = species code):

F1 "pseudo-orthology": teleost WGD (TGD) at Teleostei with differential loss — zebrafish kept
   copy a, medaka kept copy b. The two fish genes look like orthologs (S0 = Teleostei) but
   their LCA is the duplication, so merge = parent(Teleostei) = Euteleostomi.
F2 "lineage-specific duplication" in the Mammalia lineage (lvl=Mammalia): copies x and y in
   human and mouse; merge across copies = parent(Mammalia) = Amniota.
F3 "unresolved polytomy": a U node over two human genes and one mouse gene; only the three
   pairs with the chicken gene are asserted.
F4 "duplication at the family root" (no lvl → LCA mapping Amniota = family level): pairs
   across the two copies merge at ∞ (never in one HOG).
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))  # make `hobench` importable

import pandas as pd
import pytest

from hobench.newick import parse_newick
from hobench.reference import ReferenceFamily
from hobench.speciestree import SpeciesTree

TOY_TREE = "(((HUMAN,MOUSE)Mammalia,CHICK)Amniota,(DANRE,ORYLA)Teleostei)Euteleostomi;"

F1_NHX = ("(((HUMAN_1,MOUSE_1)[&&NHX:ev=S],CHICK_1)[&&NHX:ev=S],"
          "(DANRE_a,ORYLA_b)[&&NHX:ev=D:lvl=Teleostei:tag=WGD:src=S1])[&&NHX:ev=S];")
F2_NHX = ("((((HUMAN_x,MOUSE_x)[&&NHX:ev=S],(HUMAN_y,MOUSE_y)[&&NHX:ev=S])[&&NHX:ev=D:lvl=Mammalia],"
          "CHICK_g)[&&NHX:ev=S],DANRE_g)[&&NHX:ev=S];")
F3_NHX = "((HUMAN_u1,HUMAN_u2,MOUSE_u)[&&NHX:ev=U],CHICK_u)[&&NHX:ev=S];"
F4_NHX = "((HUMAN_r1,CHICK_r1)[&&NHX:ev=S],(HUMAN_r2,CHICK_r2)[&&NHX:ev=S])[&&NHX:ev=D];"


def make_family(fid: str, nwk: str, phenomena=(), **meta) -> ReferenceFamily:
    tree = parse_newick(nwk)
    genes = tree.leaf_names()
    members = pd.DataFrame({"gene_id": genes, "species": [g.split("_")[0] for g in genes],
                            "original_id_in_paper": genes, "label": "", "mapping_method": "toy",
                            "mapping_identity": "1.0", "source": "S1"})
    m = {"id": fid, "name": f"toy {fid}", "panel": "toy", "phenomena": list(phenomena),
         "sources": [{"id": "S1", "doi": "10.0000/toy", "figure": "Fig. 1", "quote": "toy"}]}
    m.update(meta)
    return ReferenceFamily(id=fid, tree=tree, members=members, meta=m)


@pytest.fixture
def st() -> SpeciesTree:
    return SpeciesTree.from_newick(TOY_TREE)


@pytest.fixture
def fams() -> list[ReferenceFamily]:
    return [make_family("F1", F1_NHX, ["WGD", "pseudo-orthology"]),
            make_family("F2", F2_NHX, ["lineage-specific-dup"]),
            make_family("F3", F3_NHX, ["polytomy"]),
            make_family("F4", F4_NHX, ["root-dup"])]


@pytest.fixture
def fam_by_id(fams) -> dict[str, ReferenceFamily]:
    return {f.id: f for f in fams}
