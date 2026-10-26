"""Rooted reference species tree with named levels.

Every node of the tree has a unique name. Leaves are species codes (e.g. ``HUMAN``);
internal nodes are taxonomic *levels* (e.g. ``Mammalia``). Species are levels too in the
sense of the core formalism (a lineage-specific duplication can sit at a species), but
predicted HOGs are only ever stored for internal levels.

All nodes get an integer index (preorder). Precomputed numpy arrays make level arithmetic
vectorisable:

* ``parent_idx[i]``  parent index (-1 for root)
* ``depth[i]``       number of edges from the root
* ``anc[a, b]``      True iff ``a`` is an ancestor of ``b`` or ``a == b`` ("a ⊇ b")
* ``lca_idx[a, b]``  index of the lowest common ancestor
"""
from __future__ import annotations

from collections import Counter
from pathlib import Path
from typing import Iterable, Optional

import numpy as np

from .newick import Node, read_newick


class SpeciesTreeError(ValueError):
    pass


class SpeciesTree:
    def __init__(self, root: Node, allow_unary: bool = False):
        self.root = root
        nodes = list(root.preorder())
        names = [n.name for n in nodes]
        unnamed = sum(1 for x in names if not x)
        if unnamed:
            raise SpeciesTreeError(
                f"species tree has {unnamed} unnamed node(s); every internal node must be named "
                f"with its taxonomic level and every leaf with its species code")
        dups = sorted(x for x, c in Counter(names).items() if c > 1)
        if dups:
            raise SpeciesTreeError(f"species tree node names are not unique: {dups[:20]}")
        if not allow_unary:
            unary = [n.name for n in nodes if len(n.children) == 1]
            if unary:
                raise SpeciesTreeError(
                    f"species tree has unary internal node(s) {unary[:10]}: a level with a single child "
                    f"has the same species set as its child; remove it (or pass allow_unary=True)")
        self._nodes = nodes
        self.names: list[str] = names
        self.index: dict[str, int] = {x: i for i, x in enumerate(names)}
        n = len(nodes)
        self.n = n
        node_idx = {id(nd): i for i, nd in enumerate(nodes)}
        self.parent_idx = np.full(n, -1, dtype=np.int64)
        self.depth = np.zeros(n, dtype=np.int64)
        for i, nd in enumerate(nodes):
            if nd.parent is not None:
                p = node_idx[id(nd.parent)]
                self.parent_idx[i] = p
                self.depth[i] = self.depth[p] + 1
        self.is_leaf = np.array([nd.is_leaf for nd in nodes], dtype=bool)
        self.children_idx: list[list[int]] = [[node_idx[id(c)] for c in nd.children] for nd in nodes]

        # subtree index ranges: in preorder each subtree is a contiguous block
        self.subtree_end = np.zeros(n, dtype=np.int64)  # exclusive end
        for i in reversed(range(n)):
            ch = self.children_idx[i]
            self.subtree_end[i] = self.subtree_end[ch[-1]] if ch else i + 1
        idx = np.arange(n)
        self.anc = (idx[:, None] <= idx[None, :]) & (idx[None, :] < self.subtree_end[:, None])

        # LCA matrix: for each node v, pairs drawn from different child subtrees have lca v
        lca = np.empty((n, n), dtype=np.int64)
        for v in range(n):
            lca[v, v] = v
        for v in reversed(range(n)):
            blocks = [np.arange(c, self.subtree_end[c]) for c in self.children_idx[v]]
            # v versus its descendants
            desc = np.arange(v + 1, self.subtree_end[v])
            lca[v, desc] = v
            lca[desc, v] = v
            for a in range(len(blocks)):
                for b in range(a + 1, len(blocks)):
                    lca[np.ix_(blocks[a], blocks[b])] = v
                    lca[np.ix_(blocks[b], blocks[a])] = v
        self.lca_idx = lca

        self.species: list[str] = [names[i] for i in range(n) if self.is_leaf[i]]
        self.levels: list[str] = [names[i] for i in range(n) if not self.is_leaf[i]]
        self._clade_cache: dict[int, frozenset[str]] = {}

    # -- construction -----------------------------------------------------------------
    @classmethod
    def from_newick(cls, path_or_text: str | Path, **kw) -> "SpeciesTree":
        return cls(read_newick(path_or_text), **kw)

    # -- basic lookups ---------------------------------------------------------------
    def idx(self, name: str) -> int:
        try:
            return self.index[name]
        except KeyError:
            raise SpeciesTreeError(f"unknown species/level {name!r}") from None

    def __contains__(self, name: str) -> bool:
        return name in self.index

    @property
    def root_name(self) -> str:
        return self.names[0]

    def parent(self, name: str) -> Optional[str]:
        p = self.parent_idx[self.idx(name)]
        return None if p < 0 else self.names[p]

    def children(self, name: str) -> list[str]:
        return [self.names[c] for c in self.children_idx[self.idx(name)]]

    def depth_of(self, name: str) -> int:
        return int(self.depth[self.idx(name)])

    def is_species(self, name: str) -> bool:
        return bool(self.is_leaf[self.idx(name)])

    # -- relations ---------------------------------------------------------------------
    def is_ancestor(self, a: str, b: str, strict: bool = False) -> bool:
        """True iff level ``a`` is an ancestor of ``b`` (or equal, unless ``strict``)."""
        ia, ib = self.idx(a), self.idx(b)
        return bool(self.anc[ia, ib]) and not (strict and ia == ib)

    def lca(self, *names: str) -> str:
        if not names:
            raise SpeciesTreeError("lca() of nothing")
        return self.names[self.lca_index(self.idx(x) for x in names)]

    def lca_index(self, idxs: Iterable[int]) -> int:
        it = iter(idxs)
        cur = next(it)
        for i in it:
            cur = int(self.lca_idx[cur, i])
        return int(cur)

    def ancestors(self, name: str, include_self: bool = True) -> list[str]:
        """Path from ``name`` up to the root."""
        out = []
        i = self.idx(name)
        if not include_self:
            i = self.parent_idx[i]
        while i >= 0:
            out.append(self.names[i])
            i = self.parent_idx[i]
        return out

    def path(self, a: str, b: str) -> list[str]:
        """Node names on the path a → lca → b (inclusive)."""
        m = self.lca(a, b)
        up = self.ancestors(a)
        up = up[: up.index(m) + 1]
        down = self.ancestors(b)
        down = down[: down.index(m)]
        return up + list(reversed(down))

    def distance(self, a: str, b: str) -> int:
        """Number of species-tree edges between two levels."""
        ia, ib = self.idx(a), self.idx(b)
        m = self.lca_idx[ia, ib]
        return int(self.depth[ia] + self.depth[ib] - 2 * self.depth[m])

    def clade_species(self, name: str) -> frozenset[str]:
        i = self.idx(name)
        if i not in self._clade_cache:
            self._clade_cache[i] = frozenset(
                self.names[j] for j in range(i, self.subtree_end[i]) if self.is_leaf[j])
        return self._clade_cache[i]

    def find_clade(self, species: Iterable[str]) -> Optional[str]:
        """Name of the node whose clade is exactly this species set, or None."""
        sp = frozenset(species)
        if not sp:
            return None
        m = self.lca(*sp)
        return m if self.clade_species(m) == sp else None

    def prune_to(self, species: Iterable[str]) -> Node:
        """Return a new Node tree restricted to ``species``, with unary nodes suppressed
        (the surviving node keeps the name of the *lowest* level, i.e. the clade's LCA)."""
        keep = set(species)
        unknown = keep - set(self.species)
        if unknown:
            raise SpeciesTreeError(f"species not in tree: {sorted(unknown)}")

        def rec(i: int) -> Optional[Node]:
            if self.is_leaf[i]:
                return Node(self.names[i]) if self.names[i] in keep else None
            kids = [k for k in (rec(c) for c in self.children_idx[i]) if k is not None]
            if not kids:
                return None
            if len(kids) == 1:
                return kids[0]
            nd = Node(self.names[i])
            for k in kids:
                nd.add_child(k)
            return nd

        out = rec(0)
        if out is None:
            raise SpeciesTreeError("prune_to(): empty species set")
        return out
