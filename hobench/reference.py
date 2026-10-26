"""Reference (literature-derived) gene families.

On-disk format — one directory per family, ``families/<FAMILY_ID>/`` (see
``design/CURATION_PROTOCOL.md``):

``family.yaml``
    ``id`` (= directory name), ``name``, ``panel``, optional ``collection`` (genome-scale
    collections such as YGOB pillars: each family is one *unit* of the collection),
    ``phenomena`` (tags), ``sources`` (list of ``{id, doi and/or pmid, figure/table/locator,
    quote, source_method}``), ``members`` (default ``members.tsv``), ``tree`` (default
    ``reference_tree.nhx``), ``root_level`` (family scope level; default = level of the tree
    root), ``root_single_copy`` (bool, default false: the family was a single gene at
    ``root_level`` — or, without ``root_level``, at the parent of the root node's level),
    ``complete_species`` (species whose membership is declared complete), ``decoys`` (gene IDs
    asserted NON-homologous to this family), ``distinct_from`` (``[{family, at_level, src}]``:
    merging this family with ``family`` in a HOG at a level ⊆ ``at_level`` is wrong;
    ``at_level: ALL`` = non-homologous), ``notes``.

``members.tsv``
    ``gene_id, species, original_id_in_paper, label, mapping_method, mapping_identity
    (percent), mapping_coverage (fraction, min over both sequences), source, locator``.

``reference_tree.nhx``
    leaves = gene IDs; internal nodes ``ev=S|D|U|H`` plus

    * D: ``lvl=<level>`` (exact) or ``lvl_lo=<level>`` / ``lvl_hi=<level|ABOVE_RANGE>``: the
      youngest / oldest species-tree node in whose ancestral lineage the duplication may have
      happened. Missing values default to (LCA mapping, tightest ancestor bound).
    * H (horizontal transfer): exactly one child carries ``recipient=1`` (the transferred
      clade); ``donor=<level|OUTSIDE>``.
    * optional ``tag=WGD|TAND|HGT|ALLO|CONV|...`` (comma-separated), ``src=<source id>``,
      ``id=<node id>``, ``tier=T1|T2|T3|T3fig``.

Reference merge levels are **intervals** ``[lo, hi]`` of species-tree nodes (``hi`` may be
``INF`` = never in one HOG within the panel). For a pair whose reference LCA node is

* S: ``lo = hi = lca(sp(a), sp(b))`` (orthologs);
* D with duplication level interval ``[Dlo, Dhi]``: ``[parent(Dlo), parent(Dhi)]`` (``parent``
  of the panel root, or ``Dhi = ABOVE_RANGE``, gives ``INF``). A duplication with no
  constraining (non-U) ancestor — a *root* duplication — always gets ``hi = INF`` unless
  ``root_single_copy`` is asserted, even with an explicit ``lvl`` (SYNTHESIS A1);
* H: ``lo = hi = lca(donor, mapping(recipient child))``, or ``INF`` if the donor is outside the
  panel (the xenolog is apart from native homologs at every recipient-lineage level);
* U: not asserted.

Default duplication bound (``Dhi`` when not given): the tightest constraint from the nearest
non-U ancestor — below a speciation the duplication happened after the lineage split from its
closest sibling lineage (``child_toward(lca(mapping, closest sibling), mapping)``); below a D/H
node it is at most that node's ``Dhi``; at the top of the tree it is ``ABOVE_RANGE`` unless
``root_single_copy`` asserts a single ancestral gene (then strictly below that level).
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

EVENTS = ("S", "D", "U", "H")
NHX_KEYS = ["id", "ev", "lvl", "lvl_lo", "lvl_hi", "donor", "recipient", "tag", "tier", "src"]
MEMBER_COLUMNS = ["gene_id", "species", "original_id_in_paper", "label",
                  "mapping_method", "mapping_identity", "mapping_coverage", "source", "locator"]
INF = -1           # merge level "never in one HOG"
NOT_ASSERTED = -2  # pair under a U node
PER_PAIR_S0 = -3   # S node: merge = species LCA of the pair
ABOVE = -1         # duplication-level bound "ABOVE_RANGE"
ABOVE_RANGE_TOKENS = ("ABOVE_RANGE", "OUTSIDE", "INF")
NONHOMOLOGOUS_TOKENS = ("ALL", "NONHOMOLOGOUS")


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
    def collection(self) -> Optional[str]:
        c = self.meta.get("collection")
        return str(c) if c else None

    @property
    def phenomena(self) -> list[str]:
        return [str(t) for t in (self.meta.get("phenomena") or [])]

    @property
    def root_level(self) -> Optional[str]:
        return self.meta.get("root_level")

    @property
    def root_single_copy(self) -> bool:
        v = self.meta.get("root_single_copy", False)
        if not isinstance(v, bool):
            raise ReferenceError(f"[{self.id}] root_single_copy must be true/false, got {v!r}")
        return v

    @property
    def complete_species(self) -> list[str]:
        return [str(x) for x in (self.meta.get("complete_species") or [])]

    @property
    def decoys(self) -> list[str]:
        return [str(x) for x in (self.meta.get("decoys") or [])]

    @property
    def distinct_from(self) -> list[dict]:
        return list(self.meta.get("distinct_from") or [])

    @property
    def gene_species(self) -> dict[str, str]:
        """gene_id -> species (cached per members-table object; returns a copy)."""
        key = id(self.members)
        cache = self.__dict__.get("_gsp_cache")
        if cache is None or cache[0] != key or cache[1] is not self.members:
            d = dict(zip(self.members["gene_id"].astype(str), self.members["species"].astype(str)))
            self.__dict__["_gsp_cache"] = cache = (key, self.members, d)
        return dict(cache[2])

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
    members = pd.read_csv(mpath, sep="\t", dtype=str, keep_default_na=False)
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
    """An internal node of the reference tree with its species-tree annotation.

    ``lvl_lo``/``lvl_hi``: duplication-level interval (D, H; ``lvl_hi == ABOVE`` = above the
    panel range); for S/U both equal the LCA mapping. ``merge_lo``/``merge_hi``: merge-level
    interval for pairs separated at this node (``PER_PAIR_S0`` for S, ``NOT_ASSERTED`` for U,
    ``INF`` = never in one HOG)."""
    node: Node
    id: str
    ev: str
    mapping: int
    lvl_lo: int
    lvl_hi: int
    level_source: str        # 'exact' | 'interval' | 'default' | 'mapping' (S/U) | 'transfer' (H)
    merge_lo: int
    merge_hi: int
    tags: list[str]
    src: Optional[str]
    tier: Optional[str]
    child_leaves: list[np.ndarray]   # gene indices per child
    recipient_child: Optional[int] = None   # H: index of the transferred child
    donor: Optional[str] = None             # H: donor level name or 'OUTSIDE'

    # backwards-compatible aliases
    @property
    def level(self) -> int:
        return self.lvl_lo

    @property
    def level_explicit(self) -> bool:
        return self.level_source in ("exact", "interval")

    @property
    def merge(self) -> int:
        return self.merge_lo

    @property
    def is_interval(self) -> bool:
        return self.ev in ("D", "H") and self.merge_lo != self.merge_hi


@dataclass
class ReconciledFamily:
    family: ReferenceFamily
    st: SpeciesTree
    genes: list[str]             # tree leaves, preorder
    gene_index: dict[str, int]
    species_idx: np.ndarray      # species-tree index per gene
    family_level: int            # species-tree index
    events: list[EventNode]      # preorder

    @property
    def id(self) -> str:
        return self.family.id


def _split_tags(v: Optional[str]) -> list[str]:
    if not v:
        return []
    return [t.strip().upper() for t in str(v).replace(";", ",").replace("|", ",").split(",") if t.strip()]


def _child_toward(st: SpeciesTree, anc: int, desc: int) -> int:
    """Child of ``anc`` on the path to ``desc`` (``desc`` strictly below ``anc``)."""
    cur = desc
    while st.parent_idx[cur] != anc:
        cur = int(st.parent_idx[cur])
        if cur < 0:
            raise ReferenceError("internal error: _child_toward called with non-ancestor")
    return cur


def _parent_or_inf(st: SpeciesTree, lvl: int) -> int:
    if lvl == ABOVE:
        return INF
    p = int(st.parent_idx[lvl])
    return p if p >= 0 else INF


def _annotate(fam: ReferenceFamily, st: SpeciesTree) -> tuple[ValidationReport, Optional[ReconciledFamily]]:
    """Single pass that validates the tree and, if valid, builds the reconciled view."""
    r = ValidationReport()
    P = f"[{fam.id}]"
    gsp = fam.gene_species
    leaves = fam.tree.leaves()
    genes = [lf.name for lf in leaves]
    gindex = {g: i for i, g in enumerate(genes)}
    sp_idx = np.array([st.idx(gsp[g]) for g in genes], dtype=np.int64)
    source_ids = {str(s.get("id")) for s in (fam.meta.get("sources") or []) if isinstance(s, dict) and s.get("id")}

    def label(nd: Node) -> str:
        nid = nd.features.get("id") or nd.name
        if nid:
            return repr(nid)
        lv = nd.leaf_names()
        return f"<{nd.features.get('ev', '?')} node over {lv[:3]}{'...' if len(lv) > 3 else ''}>"

    def level_token(nd: Node, key: str) -> Optional[int]:
        """Parse a level-valued NHX field; None if absent; ABOVE for ABOVE_RANGE; error -> None."""
        v = nd.features.get(key)
        if v is None or v == "":
            return None
        if v.upper() in ABOVE_RANGE_TOKENS:
            return ABOVE
        if v not in st:
            r.errors.append(f"{P} node {label(nd)}: {key}={v!r} is not a node of the species tree")
            return None
        return st.idx(v)

    # bottom-up: leaf sets, species sets, mappings
    leafset: dict[int, np.ndarray] = {}
    mapping: dict[int, int] = {}
    species_sets: dict[int, set] = {}
    for nd in fam.tree.postorder():
        if nd.is_leaf:
            i = gindex[nd.name]
            leafset[id(nd)] = np.array([i], dtype=np.int64)
            mapping[id(nd)] = int(sp_idx[i])
            species_sets[id(nd)] = {gsp[nd.name]}
        else:
            leafset[id(nd)] = np.concatenate([leafset[id(c)] for c in nd.children])
            mapping[id(nd)] = st.lca_index(mapping[id(c)] for c in nd.children)
            species_sets[id(nd)] = set().union(*(species_sets[id(c)] for c in nd.children))

    root_level = st.idx(fam.root_level) if fam.root_level and fam.root_level in st else None
    if fam.root_level and fam.root_level not in st:
        r.errors.append(f"{P} family.yaml: root_level {fam.root_level!r} is not a node of the species tree")

    # per-node event info; top-down so bounds are known
    info: dict[int, dict] = {}
    events: list[EventNode] = []
    order = list(fam.tree.preorder())
    parent_of: dict[int, Optional[Node]] = {id(n): n.parent for n in order}

    def own_lo(nd: Node) -> int:
        """Youngest level of a node (for sibling-disjointness): explicit lvl/lvl_lo or mapping."""
        if nd.is_leaf:
            return mapping[id(nd)]
        if id(nd) in info:
            return info[id(nd)]["lvl_lo"]
        return mapping[id(nd)]

    def ancestor_bound(nd: Node) -> tuple[int, bool, str]:
        """(bound level or ABOVE, strict?, description) from the nearest non-U ancestor."""
        child, a = nd, nd.parent
        while a is not None and a.features.get("ev", "").upper() == "U":
            child, a = a, a.parent
        if a is None:
            if fam.root_single_copy:
                if root_level is not None:
                    c_lvl = root_level
                else:
                    top = order[0]
                    t_lo = info[id(top)]["lvl_lo"] if id(top) in info and top is not nd else (
                        own_explicit_lo(top) if own_explicit_lo(top) is not None else mapping[id(top)])
                    c_lvl = int(st.parent_idx[t_lo]) if t_lo >= 0 else -1
                    if c_lvl < 0:
                        r.errors.append(f"{P} root_single_copy: the root node's level has no parent level")
                        return ABOVE, False, "top of tree"
                return c_lvl, True, f"root_single_copy at {st.names[c_lvl]!r}"
            return ABOVE, False, "top of tree"
        aev = a.features.get("ev", "").upper()
        if aev == "S":
            # closest sibling lineage of `child` among a's children
            m = mapping[id(nd)]
            best = None
            for sib in a.children:
                if sib is child:
                    continue
                l = int(st.lca_idx[m, own_lo(sib)])
                if best is None or st.depth[l] > st.depth[best]:
                    best = l
            return best, True, f"speciation {label(a)}"
        a_info = info.get(id(a))
        hi = a_info["lvl_hi"] if a_info else ABOVE
        return hi, False, f"{aev} ancestor {label(a)}"

    def own_explicit_lo(nd: Node) -> Optional[int]:
        f = nd.features
        if f.get("lvl") and f["lvl"] in st:
            return st.idx(f["lvl"])
        if f.get("lvl_lo") and f["lvl_lo"] in st:
            return st.idx(f["lvl_lo"])
        return None

    name_count: dict[str, int] = {}
    for nd in order:
        if not nd.is_leaf:
            nm = nd.features.get("id") or nd.name
            if nm:
                name_count[nm] = name_count.get(nm, 0) + 1

    for k, nd in enumerate(order):
        if nd.is_leaf:
            continue
        ev = nd.features.get("ev")
        if ev is None:
            r.errors.append(f"{P} internal node {label(nd)} has no 'ev' (must be S, D, U or H)")
            continue
        ev = ev.upper()
        if ev not in EVENTS:
            r.errors.append(f"{P} internal node {label(nd)} has ev={ev!r} (must be S, D, U or H)")
            continue
        if len(nd.children) < 2:
            r.errors.append(f"{P} internal node {label(nd)} has a single child; remove it")
        src = nd.features.get("src")
        if src and source_ids and src not in source_ids:
            r.errors.append(f"{P} node {label(nd)}: src={src!r} is not a source id in family.yaml")
        mp = mapping[id(nd)]
        lo = hi = mp
        source = "mapping"
        merge_lo = merge_hi = NOT_ASSERTED
        rec_child = None
        donor_name = None
        if ev in ("S", "U"):
            for key in ("lvl", "lvl_lo", "lvl_hi", "donor"):
                if nd.features.get(key):
                    if ev == "S" and key == "lvl" and nd.features[key] in st and st.idx(nd.features[key]) == mp:
                        continue
                    r.errors.append(f"{P} {ev} node {label(nd)}: '{key}' is only allowed on D/H nodes "
                                    f"(LCA mapping is {st.names[mp]!r})")
            merge_lo = merge_hi = PER_PAIR_S0 if ev == "S" else NOT_ASSERTED
        elif ev == "D":
            bound, strict, why = ancestor_bound(nd)
            lvl = level_token(nd, "lvl")
            lvl_lo = level_token(nd, "lvl_lo")
            lvl_hi = level_token(nd, "lvl_hi")
            if lvl is not None and (lvl_lo is not None or lvl_hi is not None):
                r.errors.append(f"{P} duplication {label(nd)}: give either lvl or lvl_lo/lvl_hi, not both")
            if lvl == ABOVE or lvl_lo == ABOVE:
                r.errors.append(f"{P} duplication {label(nd)}: ABOVE_RANGE is only allowed for lvl_hi")
                lvl = lvl_lo = None
            if lvl is not None:
                lo = hi = lvl
                source = "exact"
            else:
                lo = lvl_lo if lvl_lo is not None else mp
                if lvl_hi is not None:
                    hi = lvl_hi
                    source = "interval"
                else:
                    # default upper bound: tightest ancestor constraint
                    if bound == ABOVE:
                        hi = ABOVE
                    elif strict:
                        hi = _child_toward(st, bound, lo) if (bound != lo and st.anc[bound, lo]) else lo
                    else:
                        hi = bound
                    source = "interval" if lvl_lo is not None else "default"
            if not st.anc[lo, mp]:
                r.errors.append(f"{P} duplication {label(nd)}: level {st.names[lo]!r} must be equal or ancestral "
                                f"to the LCA of its genes' species ({st.names[mp]!r})")
            if hi != ABOVE and not st.anc[hi, lo]:
                r.errors.append(f"{P} duplication {label(nd)}: lvl_hi {st.names[hi]!r} must be equal or ancestral "
                                f"to lvl_lo {st.names[lo]!r}")
            # bound checks
            if bound != ABOVE:
                for nm, v in (("level" if lo == hi else "lvl_lo", lo), ("lvl_hi", hi)):
                    if v == ABOVE:
                        r.errors.append(f"{P} duplication {label(nd)}: lvl_hi=ABOVE_RANGE but the {why} bounds it "
                                        f"at {st.names[bound]!r}")
                        continue
                    ok = st.anc[bound, v] and not (strict and bound == v)
                    if not ok:
                        rel = "strictly below" if strict else "equal to or below"
                        r.errors.append(f"{P} duplication {label(nd)}: {nm} {st.names[v]!r} must be {rel} "
                                        f"{st.names[bound]!r} (constraint from {why})")
            merge_lo, merge_hi = _parent_or_inf(st, lo), _parent_or_inf(st, hi)
            if why == "top of tree":
                # a duplication with no constraining ancestor (root D): the panel need not contain
                # the single ancestral gene, so never merging the copies is also correct
                merge_hi = INF
        elif ev == "H":
            flagged = [i for i, c in enumerate(nd.children) if str(c.features.get("recipient", "")) in ("1", "true", "True")]
            if len(flagged) != 1:
                r.errors.append(f"{P} transfer {label(nd)}: exactly one child must carry recipient=1 "
                                f"(found {len(flagged)})")
            donor = nd.features.get("donor")
            if not donor:
                r.errors.append(f"{P} transfer {label(nd)}: donor=<level|OUTSIDE> is required "
                                f"(use ev=U if the donor relation is not source-supported)")
            for key in ("lvl", "lvl_lo", "lvl_hi"):
                if nd.features.get(key):
                    r.errors.append(f"{P} transfer {label(nd)}: '{key}' is not allowed on H nodes (use donor=)")
            if len(flagged) == 1 and donor:
                rec_child = flagged[0]
                rmap = mapping[id(nd.children[rec_child])]
                donor_name = donor
                if donor.upper() in ABOVE_RANGE_TOKENS:
                    merge_lo = merge_hi = INF
                    lo, hi = mp, ABOVE
                elif donor not in st:
                    r.errors.append(f"{P} transfer {label(nd)}: donor={donor!r} is not a node of the species tree")
                else:
                    d = st.idx(donor)
                    if st.anc[d, rmap] or st.anc[rmap, d]:
                        r.errors.append(f"{P} transfer {label(nd)}: donor {donor!r} and recipient clade "
                                        f"{st.names[rmap]!r} are nested; not a transfer")
                    m = int(st.lca_idx[d, rmap])
                    merge_lo = merge_hi = m
                    lo = hi = _child_toward(st, m, mp) if (m != mp and st.anc[m, mp]) else mp
            source = "transfer"
        info[id(nd)] = {"lvl_lo": lo, "lvl_hi": hi}
        nid = nd.features.get("id") or nd.name
        eid = (nid if name_count.get(nid) == 1 else f"{nid}#n{k}") if nid else f"n{k}"
        events.append(EventNode(
            node=nd, id=eid, ev=ev, mapping=mp, lvl_lo=lo, lvl_hi=hi, level_source=source,
            merge_lo=merge_lo, merge_hi=merge_hi, tags=_split_tags(nd.features.get("tag")), src=src,
            tier=nd.features.get("tier"), child_leaves=[leafset[id(c)] for c in nd.children],
            recipient_child=rec_child, donor=donor_name))

    # speciation consistency: children in disjoint clades (using their youngest levels)
    for nd in order:
        if nd.is_leaf or nd.features.get("ev", "").upper() != "S":
            continue
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
                        f"{s!r} but their LCA is this S node (same-species genes must split at a D, H or U node)")
                    continue
                la, lb = own_lo(kids[a]), own_lo(kids[b])
                if st.anc[la, lb] or st.anc[lb, la]:
                    r.errors.append(
                        f"{P} speciation node {label(nd)}: children {label(kids[a])} (level {st.names[la]!r}) and "
                        f"{label(kids[b])} (level {st.names[lb]!r}) overlap in the species tree; children of a "
                        f"speciation must descend into disjoint clades (check the lvl of duplications below it)")

    top_lo = info[id(fam.tree)]["lvl_lo"] if id(fam.tree) in info else (mapping[id(fam.tree)])
    if root_level is not None and not st.anc[root_level, top_lo]:
        r.errors.append(f"{P} root_level {fam.root_level!r} is not equal or ancestral to the level of the "
                        f"reference-tree root ({st.names[top_lo]!r})")
    if r.errors:
        return r, None
    fam_level = root_level if root_level is not None else top_lo
    return r, ReconciledFamily(fam, st, genes, gindex, sp_idx, fam_level, events)


def reconcile(fam: ReferenceFamily, st: SpeciesTree, check: bool = True) -> ReconciledFamily:
    """Attach species-tree levels and merge intervals to the reference tree. Always raises
    on an invalid tree; with ``check`` the metadata/members checks of
    :func:`validate_family` must pass too."""
    if check:
        validate_family(fam, st).raise_if_errors()
    rep, rf = _annotate(fam, st)
    if rf is None:
        rep.raise_if_errors()
    return rf


# ----------------------------------------------------------------------- validation
def _validate_meta_members(fam: ReferenceFamily, st: SpeciesTree, r: ValidationReport) -> set[str]:
    P = f"[{fam.id}]"
    meta = fam.meta
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
        if not any(s.get(k) for k in ("figure", "table", "locator", "files")):
            r.warnings.append(f"{P} family.yaml: source {sid!r} has no figure/table/locator")
        if not s.get("quote"):
            r.warnings.append(f"{P} family.yaml: source {sid!r} has no supporting quote")
    try:
        fam.root_single_copy
    except ReferenceError as e:
        r.errors.append(str(e))
    for s in fam.complete_species:
        if s not in st or not st.is_species(s):
            r.errors.append(f"{P} complete_species: {s!r} is not a species of the panel tree")
    for i, d in enumerate(fam.distinct_from):
        if not isinstance(d, dict) or not d.get("family") or not d.get("at_level"):
            r.errors.append(f"{P} distinct_from[{i}] needs 'family' and 'at_level'")
            continue
        lv = str(d["at_level"])
        if lv.upper() not in NONHOMOLOGOUS_TOKENS and lv not in st:
            r.errors.append(f"{P} distinct_from[{i}]: at_level {lv!r} is not a level (or ALL)")
        if not d.get("src"):
            r.errors.append(f"{P} distinct_from[{i}] ({d['family']}): a source id 'src' is required")
        elif source_ids and str(d["src"]) not in source_ids:
            r.errors.append(f"{P} distinct_from[{i}]: src={d['src']!r} is not a source id in family.yaml")

    m = fam.members
    for col in ("gene_id", "species"):
        if col not in m.columns:
            r.errors.append(f"{P} members table lacks required column '{col}'")
    if r.errors:
        return source_ids
    ids = m["gene_id"].astype(str)
    dup = sorted(set(ids[ids.duplicated()]))
    if dup:
        r.errors.append(f"{P} members table: duplicated gene_id(s) {dup[:10]}")
    if (ids == "").any():
        r.errors.append(f"{P} members table: empty gene_id in {int((ids == '').sum())} row(s)")
    bad_sp = sorted(set(m["species"].astype(str)) - set(st.species))
    if bad_sp:
        r.errors.append(f"{P} members table: species not in the panel species tree: {bad_sp}")
    if "source" in m.columns:
        unknown_src = sorted({s for s in m["source"].astype(str) if s and s not in source_ids})
        if unknown_src and source_ids:
            r.warnings.append(f"{P} members table: source id(s) {unknown_src[:5]} not declared in family.yaml")
    overlap = sorted(set(fam.decoys) & set(ids))
    if overlap:
        r.errors.append(f"{P} decoys are also members: {overlap[:5]}")
    return source_ids


def validate_family(fam: ReferenceFamily, st: SpeciesTree) -> ValidationReport:
    """Check metadata, members table and reference tree against the species tree."""
    r = ValidationReport()
    P = f"[{fam.id}]"
    _validate_meta_members(fam, st, r)
    if r.errors:
        return r
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
    rep, _ = _annotate(fam, st)
    r.extend(rep)
    return r


def validate_panel(families: Iterable[ReferenceFamily], st: SpeciesTree) -> ValidationReport:
    """Validate every family, then cross-family constraints: gene IDs unique across families,
    ``distinct_from`` targets exist."""
    r = ValidationReport()
    families = list(families)
    owner: dict[str, str] = {}
    ids = {f.id for f in families}
    for fam in families:
        r.extend(validate_family(fam, st))
        for g in fam.member_ids:
            if g in owner and owner[g] != fam.id:
                r.errors.append(f"gene {g!r} is a member of both {owner[g]!r} and {fam.id!r}")
            owner.setdefault(g, fam.id)
        for d in fam.distinct_from:
            if isinstance(d, dict) and d.get("family") and d["family"] not in ids:
                r.warnings.append(f"[{fam.id}] distinct_from family {d['family']!r} is not in this panel")
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

