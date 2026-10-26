"""Reference (literature-derived) gene families.

On-disk format — one directory per family, ``families/<FAMILY_ID>/``:

``family.yaml``
    ``id`` (= directory name), ``name``, ``panel``, ``phenomena`` (list of tags, e.g.
    ``[WGD, HGT]``), ``sources`` (list of ``{id, doi and/or pmid, figure/table/locator, quote}``),
    ``members`` (path of the members table, default ``members.tsv``), optional ``tree``
    (default ``reference_tree.nhx``), optional ``root_level`` (the level at which the family
    is asserted to be a single ancestral gene; default = level of the reference-tree root),
    optional ``decoys`` (gene IDs asserted to be NON-homologous to this family), ``notes``.

``members.tsv``
    columns ``gene_id, species, original_id_in_paper, label, mapping_method,
    mapping_identity, source`` (only ``gene_id`` and ``species`` are required).

``reference_tree.nhx``
    rooted, possibly multifurcating tree; leaves are gene IDs; every internal node carries
    ``ev=S|D|U`` and optionally ``lvl=<level>`` (duplications: the species-tree node in whose
    ancestral lineage the duplication happened), ``tag=WGD|HGT|ILS|...`` (several tags may be
    joined with ``,``) and ``src=<source id>``.

The *reconciled view* (:class:`ReconciledFamily`) attaches species-tree indices to every
internal node: its LCA mapping, its level (explicit ``lvl`` for D, else the mapping) and the
merge level that pairs separated at that node get (see ``design/core_spec.md``).
"""
from __future__ import annotations

import csv
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, Optional

import numpy as np
import pandas as pd
import yaml

from .newick import Node, read_newick, write_newick
from .speciestree import SpeciesTree

EVENTS = ("S", "D", "U")
NHX_KEYS = ["ev", "lvl", "tag", "src"]
MEMBER_COLUMNS = ["gene_id", "species", "original_id_in_paper", "label",
                  "mapping_method", "mapping_identity", "source"]
INF = -1           # merge level "never in one HOG"
NOT_ASSERTED = -2  # pair under a U node


class ReferenceError(ValueError):
    pass


@dataclass
class ValidationReport:
    errors: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)

    @property
    def ok(self) -> bool:
        return not self.errors

    def extend(self, other: "ValidationReport") -> None:
        self.errors.extend(other.errors)
        self.warnings.extend(other.warnings)

    def raise_if_errors(self) -> None:
        if self.errors:
            raise ReferenceError("reference validation failed:\n  - " + "\n  - ".join(self.errors))

    def __str__(self) -> str:
        lines = [f"ERROR: {e}" for e in self.errors] + [f"WARNING: {w}" for w in self.warnings]
        return "\n".join(lines) if lines else "OK"


@dataclass
class ReferenceFamily:
    id: str
    tree: Node
    members: pd.DataFrame
    meta: dict = field(default_factory=dict)
    path: Optional[Path] = None

    # -- convenience accessors -----------------------------------------------------
    @property
    def name(self) -> str:
        return self.meta.get("name", self.id)

    @property
    def panel(self) -> Optional[str]:
        return self.meta.get("panel")

    @property
    def phenomena(self) -> list[str]:
        return [str(t) for t in (self.meta.get("phenomena") or [])]

    @property
    def root_level(self) -> Optional[str]:
        return self.meta.get("root_level")

    @property
    def decoys(self) -> list[str]:
        return [str(x) for x in (self.meta.get("decoys") or [])]

    @property
    def gene_species(self) -> dict[str, str]:
        return dict(zip(self.members["gene_id"].astype(str), self.members["species"].astype(str)))

    @property
    def member_ids(self) -> list[str]:
        return list(self.members["gene_id"].astype(str))


# ----------------------------------------------------------------------------- I/O
def load_family(directory: str | Path) -> ReferenceFamily:
    d = Path(directory)
    ypath = d / "family.yaml"
    if not ypath.exists():
        raise ReferenceError(f"{d}: missing family.yaml")
    meta = yaml.safe_load(ypath.read_text()) or {}
    fid = str(meta.get("id", d.name))
    mpath = d / str(meta.get("members", "members.tsv"))
    tpath = d / str(meta.get("tree", "reference_tree.nhx"))
    if not mpath.exists():
        raise ReferenceError(f"{fid}: members table {mpath} not found")
    if not tpath.exists():
        raise ReferenceError(f"{fid}: reference tree {tpath} not found")
    members = pd.read_csv(mpath, sep="\t", dtype=str, keep_default_na=False, comment=None)
    tree = read_newick(tpath.read_text())
    return ReferenceFamily(id=fid, tree=tree, members=members, meta=meta, path=d)


def load_families(root: str | Path, panel: Optional[str] = None) -> list[ReferenceFamily]:
    """Load every ``<root>/*/family.yaml`` (sorted by directory name); optionally keep only
    families whose ``panel`` field equals ``panel``."""
    out = []
    for d in sorted(Path(root).iterdir()):
        if d.is_dir() and (d / "family.yaml").exists():
            fam = load_family(d)
            if panel is None or fam.panel == panel:
                out.append(fam)
    return out


def write_family(directory: str | Path, meta: dict, members: pd.DataFrame, tree: Node) -> Path:
    """Write a family directory (family.yaml, members.tsv, reference_tree.nhx)."""
    d = Path(directory)
    d.mkdir(parents=True, exist_ok=True)
    meta = dict(meta)
    meta.setdefault("members", "members.tsv")
    meta.setdefault("tree", "reference_tree.nhx")
    (d / "family.yaml").write_text(yaml.safe_dump(meta, sort_keys=False, allow_unicode=True))
    cols = [c for c in MEMBER_COLUMNS if c in members.columns] + \
           [c for c in members.columns if c not in MEMBER_COLUMNS]
    members[cols].to_csv(d / meta["members"], sep="\t", index=False, quoting=csv.QUOTE_MINIMAL)
    (d / meta["tree"]).write_text(write_newick(tree, features=NHX_KEYS, lengths=False) + "\n")
    return d


# ------------------------------------------------------------------- reconciled view
@dataclass
class EventNode:
    """An internal node of the reference tree with its species-tree annotation."""
    node: Node
    id: str                 # node name, or an automatic id like "n12"
    ev: str                 # S, D or U
    mapping: int            # species-tree index of the LCA of its leaves' species
    level: int              # D: explicit lvl (else mapping); S/U: mapping
    level_explicit: bool
    merge: int              # D: merge level index or INF; S: -3 (per-pair S0); U: NOT_ASSERTED
    tags: list[str]
    src: Optional[str]
    child_leaves: list[np.ndarray]  # gene indices (into ReconciledFamily.genes) per child


PER_PAIR_S0 = -3


@dataclass
class ReconciledFamily:
    family: ReferenceFamily
    st: SpeciesTree
    genes: list[str]             # tree leaves, preorder
    gene_index: dict[str, int]
    species_idx: np.ndarray      # species-tree index per gene
    family_level: int            # species-tree index
    events: list[EventNode]

    @property
    def id(self) -> str:
        return self.family.id


def _split_tags(v: Optional[str]) -> list[str]:
    if not v:
        return []
    return [t.strip().upper() for t in str(v).replace(";", ",").replace("|", ",").split(",") if t.strip()]


def reconcile(fam: ReferenceFamily, st: SpeciesTree, check: bool = True) -> ReconciledFamily:
    """Attach species-tree levels to the reference tree. Raises if the family is invalid
    (when ``check``)."""
    if check:
        validate_family(fam, st).raise_if_errors()
    gsp = fam.gene_species
    leaves = fam.tree.leaves()
    genes = [lf.name for lf in leaves]
    gindex = {g: i for i, g in enumerate(genes)}
    sp_idx = np.array([st.idx(gsp[g]) for g in genes], dtype=np.int64)

    # leaf-index sets and mappings bottom-up
    leafset: dict[int, np.ndarray] = {}
    mapping: dict[int, int] = {}
    for nd in fam.tree.postorder():
        if nd.is_leaf:
            i = gindex[nd.name]
            leafset[id(nd)] = np.array([i], dtype=np.int64)
            mapping[id(nd)] = int(sp_idx[i])
        else:
            leafset[id(nd)] = np.concatenate([leafset[id(c)] for c in nd.children])
            mapping[id(nd)] = st.lca_index(mapping[id(c)] for c in nd.children)

    def node_level(nd: Node) -> tuple[int, bool]:
        if nd.is_leaf:
            return mapping[id(nd)], False
        ev = nd.features.get("ev", "U").upper()
        lvl = nd.features.get("lvl")
        if ev == "D" and lvl:
            return st.idx(lvl), True
        return mapping[id(nd)], False

    root_lvl, _ = node_level(fam.tree)
    fam_level = st.idx(fam.root_level) if fam.root_level else root_lvl

    events: list[EventNode] = []
    name_count: dict[str, int] = {}
    for nd in fam.tree.preorder():
        if not nd.is_leaf and nd.name:
            name_count[nd.name] = name_count.get(nd.name, 0) + 1
    for k, nd in enumerate(fam.tree.preorder()):
        if nd.is_leaf:
            continue
        ev = nd.features.get("ev", "U").upper()
        lvl, explicit = node_level(nd)
        if ev == "S":
            merge = PER_PAIR_S0
        elif ev == "D":
            p = int(st.parent_idx[lvl])
            merge = p if (p >= 0 and st.anc[fam_level, p]) else INF
        else:
            merge = NOT_ASSERTED
        events.append(EventNode(
            node=nd, id=(nd.name if name_count.get(nd.name) == 1 else f"{nd.name}#n{k}") if nd.name else f"n{k}",
            ev=ev, mapping=mapping[id(nd)], level=lvl,
            level_explicit=explicit, merge=merge, tags=_split_tags(nd.features.get("tag")),
            src=nd.features.get("src"), child_leaves=[leafset[id(c)] for c in nd.children]))
    return ReconciledFamily(fam, st, genes, gindex, sp_idx, fam_level, events)


# ----------------------------------------------------------------------- validation
def validate_family(fam: ReferenceFamily, st: SpeciesTree) -> ValidationReport:
    """Check metadata, members table and reference tree against the species tree."""
    r = ValidationReport()
    P = f"[{fam.id}]"
    meta = fam.meta

    # --- metadata
    if fam.path is not None and fam.path.name != fam.id:
        r.warnings.append(f"{P} family id {fam.id!r} differs from directory name {fam.path.name!r}")
    for key in ("id", "name", "panel", "sources"):
        if not meta.get(key):
            r.errors.append(f"{P} family.yaml: missing required field '{key}'")
    source_ids = set()
    for i, s in enumerate(meta.get("sources") or []):
        if not isinstance(s, dict):
            r.errors.append(f"{P} family.yaml: sources[{i}] must be a mapping")
            continue
        sid = s.get("id")
        if not sid:
            r.errors.append(f"{P} family.yaml: sources[{i}] has no 'id'")
        else:
            source_ids.add(str(sid))
        if not (s.get("doi") or s.get("pmid")):
            r.errors.append(f"{P} family.yaml: source {sid!r} needs a 'doi' or 'pmid'")
        if not any(s.get(k) for k in ("figure", "table", "locator")):
            r.warnings.append(f"{P} family.yaml: source {sid!r} has no figure/table/locator")
        if not s.get("quote"):
            r.warnings.append(f"{P} family.yaml: source {sid!r} has no supporting quote")
    if fam.root_level and fam.root_level not in st:
        r.errors.append(f"{P} family.yaml: root_level {fam.root_level!r} is not a node of the species tree")

    # --- members table
    m = fam.members
    for col in ("gene_id", "species"):
        if col not in m.columns:
            r.errors.append(f"{P} members table lacks required column '{col}'")
    if r.errors:
        return r
    ids = m["gene_id"].astype(str)
    dup = sorted(set(ids[ids.duplicated()]))
    if dup:
        r.errors.append(f"{P} members table: duplicated gene_id(s) {dup[:10]}")
    if (ids == "").any():
        r.errors.append(f"{P} members table: empty gene_id in {int((ids == '').sum())} row(s)")
    panel_species = set(st.species)
    bad_sp = sorted(set(m["species"].astype(str)) - panel_species)
    if bad_sp:
        r.errors.append(f"{P} members table: species not in the panel species tree: {bad_sp}")
    if "source" in m.columns:
        unknown_src = sorted({s for s in m["source"].astype(str) if s and s not in source_ids})
        if unknown_src and source_ids:
            r.warnings.append(f"{P} members table: source id(s) {unknown_src[:5]} not declared in family.yaml")

    # --- tree
    gsp = fam.gene_species
    leaves = fam.tree.leaf_names()
    seen: set[str] = set()
    dups = sorted({x for x in leaves if x in seen or seen.add(x)})
    if dups:
        r.errors.append(f"{P} reference tree: duplicated leaf gene ID(s) {dups[:10]}")
    not_member = sorted(set(leaves) - set(gsp))
    if not_member:
        r.errors.append(f"{P} reference tree: leaves not listed in members table: {not_member[:10]}")
    unplaced = sorted(set(gsp) - set(leaves))
    if unplaced:
        r.warnings.append(f"{P} {len(unplaced)} member(s) not placed in the reference tree "
                          f"(count for family integrity only): {unplaced[:5]}")
    if r.errors:
        return r
    if any(gsp[g] not in panel_species for g in leaves):
        return r  # already reported

    mapping: dict[int, int] = {}
    species_sets: dict[int, set[str]] = {}
    for nd in fam.tree.postorder():
        if nd.is_leaf:
            mapping[id(nd)] = st.idx(gsp[nd.name])
            species_sets[id(nd)] = {gsp[nd.name]}
            continue
        mapping[id(nd)] = st.lca_index(mapping[id(c)] for c in nd.children)
        species_sets[id(nd)] = set().union(*(species_sets[id(c)] for c in nd.children))

    def label(nd: Node) -> str:
        if nd.name:
            return repr(nd.name)
        lv = nd.leaf_names()
        return f"<{nd.features.get('ev', '?')} node over {lv[:3]}{'...' if len(lv) > 3 else ''}>"

    # explicit/implicit level of a node, for constraint checks
    def level_of(nd: Node) -> int:
        if not nd.is_leaf and nd.features.get("ev", "").upper() == "D" and nd.features.get("lvl") in st:
            return st.idx(nd.features["lvl"])
        return mapping[id(nd)]

    for nd in fam.tree.preorder():
        if nd.is_leaf:
            continue
        ev = nd.features.get("ev")
        if ev is None:
            r.errors.append(f"{P} internal node {label(nd)} has no 'ev' (must be S, D or U)")
            continue
        ev = ev.upper()
        if ev not in EVENTS:
            r.errors.append(f"{P} internal node {label(nd)} has ev={ev!r} (must be S, D or U)")
            continue
        if len(nd.children) < 2:
            r.errors.append(f"{P} internal node {label(nd)} has a single child; remove it")
        lvl = nd.features.get("lvl")
        if lvl:
            if lvl not in st:
                r.errors.append(f"{P} node {label(nd)}: lvl={lvl!r} is not a node of the species tree")
            elif ev == "D":
                mp = mapping[id(nd)]
                if not st.anc[st.idx(lvl), mp]:
                    r.errors.append(
                        f"{P} duplication {label(nd)}: lvl={lvl!r} must be equal or ancestral to the LCA "
                        f"of its genes' species ({st.names[mp]!r})")
            elif ev == "S" and st.idx(lvl) != mapping[id(nd)]:
                r.warnings.append(f"{P} speciation {label(nd)}: lvl={lvl!r} differs from its LCA mapping "
                                  f"{st.names[mapping[id(nd)]]!r}; lvl is ignored for S nodes")
        src = nd.features.get("src")
        if src and source_ids and src not in source_ids:
            r.errors.append(f"{P} node {label(nd)}: src={src!r} is not a source id in family.yaml")

        if ev == "S":
            kids = nd.children
            for a in range(len(kids)):
                for b in range(a + 1, len(kids)):
                    common = species_sets[id(kids[a])] & species_sets[id(kids[b])]
                    if common:
                        s = sorted(common)[0]
                        ga = next(x for x in kids[a].leaf_names() if gsp[x] == s)
                        gb = next(x for x in kids[b].leaf_names() if gsp[x] == s)
                        r.errors.append(
                            f"{P} speciation node {label(nd)}: genes {ga!r} and {gb!r} are both from species "
                            f"{s!r} but their LCA is this S node (same-species genes must split at a D or U node)")
                        continue
                    la, lb = level_of(kids[a]), level_of(kids[b])
                    if st.anc[la, lb] or st.anc[lb, la]:
                        r.errors.append(
                            f"{P} speciation node {label(nd)}: children {label(kids[a])} (level "
                            f"{st.names[la]!r}) and {label(kids[b])} (level {st.names[lb]!r}) overlap in the "
                            f"species tree; children of a speciation must descend into disjoint clades "
                            f"(check the lvl of duplications below it)")

    # top-down level bounds: a D's lvl must lie below an S ancestor's mapping (strictly)
    # and within a D ancestor's lvl (non-strictly); U passes its bound through.
    root_bound: Optional[tuple[int, bool]] = None
    if fam.root_level and fam.root_level in st:
        root_bound = (st.idx(fam.root_level), False)
        rl = level_of(fam.tree)
        if not st.anc[root_bound[0], rl]:
            r.errors.append(f"{P} root_level {fam.root_level!r} is not equal or ancestral to the level of the "
                            f"reference-tree root ({st.names[rl]!r})")
    stack: list[tuple[Node, Optional[tuple[int, bool]]]] = [(fam.tree, root_bound)]
    while stack:
        nd, bound = stack.pop()
        if nd.is_leaf:
            continue
        ev = nd.features.get("ev", "U").upper()
        if ev == "D" and nd.features.get("lvl") in st and bound is not None:
            lv = st.idx(nd.features["lvl"])
            b, strict = bound
            if not st.anc[b, lv] or (strict and b == lv):
                rel = "strictly below" if strict else "equal to or below"
                r.errors.append(f"{P} duplication {label(nd)}: lvl={nd.features['lvl']!r} must be {rel} "
                                f"{st.names[b]!r} (constraint from an ancestor node)")
        if ev == "S":
            child_bound = (mapping[id(nd)], True)
        elif ev == "D":
            child_bound = (level_of(nd), False)
        else:
            child_bound = bound
        stack.extend((c, child_bound) for c in nd.children)
    return r


def validate_panel(families: Iterable[ReferenceFamily], st: SpeciesTree) -> ValidationReport:
    """Validate every family, then cross-family constraints (gene IDs unique across families)."""
    r = ValidationReport()
    owner: dict[str, str] = {}
    for fam in families:
        r.extend(validate_family(fam, st))
        for g in fam.member_ids:
            if g in owner and owner[g] != fam.id:
                r.errors.append(f"gene {g!r} is a member of both {owner[g]!r} and {fam.id!r}")
            owner.setdefault(g, fam.id)
    return r


# ------------------------------------------------------- tree from a classification
def tree_from_classification(
    table: pd.DataFrame,
    st: SpeciesTree,
    groups: Optional[pd.DataFrame] = None,
    sep: str = "|",
    inparalog_event: str = "U",
) -> Node:
    """Build a reference tree from a hierarchical classification.

    ``table`` has one row per gene with columns ``gene``, ``species``, ``path`` (e.g.
    ``HoxA|HoxA1``) and optionally ``events`` / ``levels`` / ``tags`` (``sep``-separated,
    aligned with the path components: entry *i* annotates the node that splits group
    ``path[:i]`` into its sub-groups; entry 0 is the family root).

    ``groups`` (alternative/complement) has columns ``group`` (a path prefix; ``""`` = root),
    ``ev``, and optional ``lvl``, ``tag``, ``src`` — the annotation of the node splitting that
    group into its sub-groups. Unannotated splits default to ``U`` (no assertion).

    Within a terminal group, genes are assumed to be orthologs following the species tree:
    S nodes at the species-tree topology restricted to the group's species; several genes of
    one species are joined under an ``inparalog_event`` node (``U`` = no assertion between
    them, the default; ``D`` = lineage-specific duplication at that species).

    Genes classified only down to a group that has sub-groups are allowed only if that group
    is annotated ``U``.
    """
    req = {"gene", "species", "path"}
    if not req <= set(table.columns):
        raise ReferenceError(f"classification table needs columns {sorted(req)}")
    ann: dict[str, dict] = {}

    def put(group: str, key: str, value: str) -> None:
        if value is None or (isinstance(value, float) and np.isnan(value)) or str(value) == "":
            return
        d = ann.setdefault(group, {})
        if key in d and d[key] != str(value):
            raise ReferenceError(f"inconsistent {key} for group {group or '<root>'!r}: {d[key]!r} vs {value!r}")
        d[key] = str(value)

    # tree of groups
    root: dict = {"sub": {}, "genes": []}
    for row in table.itertuples(index=False):
        comps = [c for c in str(row.path).split(sep) if c != ""]
        node = root
        for c in comps:
            node = node["sub"].setdefault(c, {"sub": {}, "genes": []})
        node["genes"].append((str(row.gene), str(row.species)))
        for col, key in (("events", "ev"), ("levels", "lvl"), ("tags", "tag")):
            if col in table.columns:
                vals = str(getattr(row, col)).split(sep) if getattr(row, col) not in (None, "") else []
                for i, v in enumerate(vals):
                    if i <= len(comps):
                        put(sep.join(comps[:i]), key, v.strip())
    if groups is not None:
        for row in groups.to_dict("records"):
            g = str(row.get("group", "") or "")
            if g in ("<root>", "*"):
                g = ""
            for key in ("ev", "lvl", "tag", "src"):
                if key in row:
                    put(g, key, row[key])

    def species_subtree(genes: list[tuple[str, str]]) -> Node:
        by_sp: dict[str, list[str]] = {}
        for g, s in genes:
            by_sp.setdefault(s, []).append(g)
        pruned = st.prune_to(by_sp)

        def conv(n: Node) -> Node:
            if n.is_leaf:
                gs = by_sp[n.name]
                if len(gs) == 1:
                    return Node(gs[0])
                out = Node(features={"ev": inparalog_event})
                if inparalog_event == "D":
                    out.features["lvl"] = n.name
                for g in gs:
                    out.add_child(Node(g))
                return out
            out = Node(features={"ev": "S"})
            for c in n.children:
                out.add_child(conv(c))
            return out

        return conv(pruned)

    def build(name: str, path: str, g: dict) -> Node:
        if not g["sub"]:
            nd = species_subtree(g["genes"])
            if name and not nd.is_leaf and not nd.name:
                nd.name = name
            return nd
        a = ann.get(path, {})
        ev = a.get("ev", "U").upper()
        kids = [build(c, f"{path}{sep}{c}" if path else c, sub) for c, sub in g["sub"].items()]
        if g["genes"]:
            if ev != "U":
                raise ReferenceError(
                    f"group {path or '<root>'!r} has sub-groups and also genes classified only down to it "
                    f"({[x for x, _ in g['genes']][:3]}); annotate the group as ev=U or classify them fully")
            kids.extend(Node(x) for x, _ in g["genes"])
        if len(kids) == 1:
            return kids[0]
        nd = Node(name=name, features={"ev": ev})
        for key in ("lvl", "tag", "src"):
            if key in a:
                nd.features[key] = a[key]
        for k in kids:
            nd.add_child(k)
        return nd

    return build("", "", root)
