"""Normalised predicted HOGs and parsers for tool outputs.

A :class:`Prediction` stores, for each *internal* level ``L`` of the panel species tree, a
mapping ``gene -> HOG id at L`` (only for genes of species below ``L``), plus the total size of
every HOG. Leaf (species) levels are implicit: every gene is its own HOG there.

To keep memory bounded on whole-proteome outputs, parsers accept ``keep_genes``: only those
genes' memberships are stored, but HOG sizes always count *all* genes (needed for the
contamination metric). The evaluator passes the union of all reference members and decoys.

Parsers
-------
* :func:`read_orthoxml` — OrthoXML 0.3/0.5 (FastOMA, OMA). Level-induction follows pyHam
  (see the function docstring and ``design/ENGINEERING_NOTES.md``).
* :func:`read_orthofinder` — OrthoFinder ``Phylogenetic_Hierarchical_Orthogroups/N*.tsv``.
* :func:`read_tsv` — simple ``level, hog_id, gene[, species]`` table.
"""
from __future__ import annotations

import csv
import itertools
import re
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Iterable, Optional

from .ids import Normalizer, make_normalizer
from .newick import read_newick
from .speciestree import SpeciesTree


class PredictionError(ValueError):
    pass


@dataclass
class Prediction:
    levels: dict[str, dict[str, str]] = field(default_factory=dict)   # level -> gene -> hog
    sizes: dict[str, dict[str, int]] = field(default_factory=dict)    # level -> hog -> n genes (all)
    gene_species: dict[str, str] = field(default_factory=dict)        # kept genes only
    meta: dict = field(default_factory=dict)
    _members: Optional[dict[str, dict[str, list[str]]]] = field(default=None, repr=False)

    def add(self, level: str, hog: str, gene: Optional[str], strict: bool = True) -> None:
        """Register a gene in a HOG; ``gene=None`` only increments the HOG size (gene not kept).

        A gene may belong to only one HOG per level. On a second, different assignment we raise
        (``strict``) or keep the first one and count it in ``meta['membership_conflicts']``."""
        self.sizes.setdefault(level, {})
        self.sizes[level][hog] = self.sizes[level].get(hog, 0) + 1
        if gene is not None:
            lv = self.levels.setdefault(level, {})
            if gene in lv and lv[gene] != hog:
                if strict:
                    raise PredictionError(f"gene {gene!r} assigned to two HOGs at level {level!r}: "
                                          f"{lv[gene]!r} and {hog!r}")
                self.meta["membership_conflicts"] = self.meta.get("membership_conflicts", 0) + 1
                return
            lv[gene] = hog
        self._members = None

    def hog(self, level: str, gene: str) -> Optional[str]:
        return self.levels.get(level, {}).get(gene)

    def hog_size(self, level: str, hog: str) -> int:
        return self.sizes.get(level, {}).get(hog, 0)

    def hog_members(self, level: str, hog: str) -> list[str]:
        """Kept genes of a HOG (all genes if the parser kept everything)."""
        if self._members is None:
            mem: dict[str, dict[str, list[str]]] = {}
            for lv, m in self.levels.items():
                d: dict[str, list[str]] = defaultdict(list)
                for g, h in m.items():
                    d[h].append(g)
                mem[lv] = dict(d)
            self._members = mem
        return self._members.get(level, {}).get(hog, [])

    def genes(self) -> set[str]:
        out: set[str] = set()
        for m in self.levels.values():
            out.update(m)
        return out

    def check(self, st: SpeciesTree) -> list[str]:
        """Return problems: genes at a level their species is not below, non-nested HOGs."""
        problems = []
        for lv, m in self.levels.items():
            if lv not in st:
                problems.append(f"level {lv!r} not in species tree")
                continue
            clade = st.clade_species(lv)
            bad = [g for g in m if g in self.gene_species and self.gene_species[g] not in clade]
            if bad:
                problems.append(f"level {lv!r}: {len(bad)} gene(s) whose species is not below it, e.g. {bad[:3]}")
        # nesting: a HOG at L maps into a single HOG at parent(L)
        for lv, m in self.levels.items():
            if lv not in st:
                continue
            par = st.parent(lv)
            if par is None or par not in self.levels:
                continue
            pm = self.levels[par]
            target: dict[str, set] = defaultdict(set)
            for g, h in m.items():
                if g in pm:
                    target[h].add(pm[g])
            split = [h for h, t in target.items() if len(t) > 1]
            if split:
                problems.append(f"level {lv!r}: {len(split)} HOG(s) not nested in a single HOG at {par!r}")
        return problems

    def write_tsv(self, path: str | Path) -> None:
        """Write the kept memberships in the simple TSV format (level, hog_id, gene, species)."""
        with open(path, "w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t", lineterminator="\n")
            w.writerow(["level", "hog_id", "gene", "species"])
            for lv in sorted(self.levels):
                for g, h in sorted(self.levels[lv].items(), key=lambda x: (x[1], x[0])):
                    w.writerow([lv, h, g, self.gene_species.get(g, "")])


def _keep(keep_genes: Optional[set[str]], g: str) -> Optional[str]:
    return g if keep_genes is None or g in keep_genes else None


# ============================================================================ OrthoXML
def _local(tag) -> str:
    return tag.rsplit("}", 1)[-1] if isinstance(tag, str) else ""


@dataclass(eq=False)
class _H:
    """A HOG node while building the induced HOG tree (orthoXML parser)."""
    level: int
    children: list = field(default_factory=list)   # of _H or _G
    dup: object = None                              # duplication key if this is a paralog copy


@dataclass(eq=False)
class _G:
    gene: Optional[str]   # normalised ID, or None if not kept
    level: int            # species index
    dup: object = None


def read_orthoxml(
    path: str | Path,
    st: SpeciesTree,
    normalize: str | Normalizer | None = "auto",
    id_attr: str = "protId",
    species_map: Optional[dict[str, str]] = None,
    level_map: Optional[dict[str, str]] = None,
    keep_genes: Optional[set[str]] = None,
    strict: bool = False,
    ignore_unknown_species: bool = False,
) -> Prediction:
    """Parse OrthoXML (0.3 / 0.5) into a :class:`Prediction`, inducing HOGs at skipped levels.

    **Level of a group.** An ``orthologGroup``'s level is, in order of preference, the name of
    its ``taxonId`` in the file's ``<taxonomy>`` (OrthoXML 0.5), its ``TaxRange`` property
    (``level_map`` can rename either), or the MRCA of its children (as pyHam does). A
    ``geneRef`` sits at its species.

    **Level induction (pyHam semantics).** A HOG exists at every level between a group's level
    and the level of the group/gene below it:

    1. *Duplications.* All copies of one ``paralogGroup`` (directly nested paralogGroups are
       one duplication, as in pyHam) duplicated on the branch just above the MRCA of the
       copies' levels: the duplication's level is ``M = parent(MRCA(copy levels))`` (pyHam's
       ``DuplicationNode.set_MRCA``). If ``M`` is below the enclosing group's level ``X``, an
       induced HOG at ``M`` holds the copies (they are together at ``M`` and above). Each copy
       then gets its own chain of induced single-lineage HOGs from below ``M`` down to its own
       level: copies are never together below ``M``. If ``M`` is above ``X`` the file is
       inconsistent (pyHam raises); with ``strict=False`` we place the duplication at ``X`` and
       record a warning.
    2. *Skipped levels (non-paralog children).* Children of a HOG at ``X`` that are not
       paralog copies and lie strictly below the same child level ``C`` of ``X`` are grouped
       into one induced HOG at ``C`` (recursively), like pyHam's ``_resolve_missing_levels``:
       no duplication was recorded between them, so they descend from one gene at ``C``.
       Children sitting exactly at ``C`` remain separate HOGs at ``C``.
    3. A child group at the same level as its parent is spliced into the parent (pyHam
       collapses such "same-level" HOGs).

    Each explicit or induced HOG node at level ``L`` yields one HOG at ``L`` whose genes are all
    genes below that node.
    """
    from lxml import etree

    norm = make_normalizer(normalize)
    species_map = species_map or {}
    level_map = level_map or {}
    pred = Prediction(meta={"format": "orthoxml", "source": str(path)})
    warnings: list[str] = []
    gene_info: dict[str, tuple[Optional[str], int]] = {}  # xml gene id -> (kept norm id, species idx)
    taxon_names: dict[str, str] = {}
    unknown_species: set[str] = set()
    current_species: Optional[str] = None
    hog_counter = itertools.count()
    stats = Counter()

    def level_idx(name: Optional[str]) -> Optional[int]:
        if name is None:
            return None
        name = level_map.get(name, name)
        return st.index.get(name)

    # ---- tree building helpers ----
    def parse_group(el) -> Optional[_H]:
        items: list = []
        props: dict[str, str] = {}

        def collect(e, dupkey):
            for ch in e:
                t = _local(ch.tag)
                if t == "geneRef":
                    gi = gene_info.get(ch.get("id"))
                    if gi is None:
                        stats["generef_unknown"] += 1
                        continue
                    items.append(_G(gi[0], gi[1], dupkey))
                elif t == "orthologGroup":
                    h = parse_group(ch)
                    if h is not None:
                        h.dup = dupkey
                        items.append(h)
                elif t == "paralogGroup":
                    collect(ch, dupkey if dupkey is not None else object())

        for ch in el:
            t = _local(ch.tag)
            if t == "property":
                props[ch.get("name")] = ch.get("value")
            elif t == "geneRef":
                gi = gene_info.get(ch.get("id"))
                if gi is None:
                    stats["generef_unknown"] += 1
                    continue
                items.append(_G(gi[0], gi[1], None))
            elif t == "orthologGroup":
                h = parse_group(ch)
                if h is not None:
                    items.append(h)
            elif t == "paralogGroup":
                collect(ch, object())
        if not items:
            return None

        child_levels = [it.level for it in items]
        lvl = None
        if el.get("taxonId") is not None:
            lvl = level_idx(taxon_names.get(el.get("taxonId")))
        if lvl is None:
            for key in ("TaxRange", "taxid"):
                if props.get(key) is not None:
                    lvl = level_idx(props[key])
                    if lvl is not None:
                        break
                    stats["taxrange_unknown"] += 1
        mrca = st.lca_index(child_levels)
        if lvl is None:
            # pyHam: a single child level -> its parent; otherwise the MRCA of the children
            if len(set(child_levels)) == 1:
                lvl = int(st.parent_idx[mrca]) if st.parent_idx[mrca] >= 0 else mrca
            else:
                lvl = mrca
        elif not st.anc[lvl, mrca]:
            msg = (f"group {el.get('id')!r}: declared level {st.names[lvl]!r} is not above all its "
                   f"members (their MRCA is {st.names[mrca]!r})")
            if strict:
                raise PredictionError(msg)
            warnings.append(msg)
            lvl = mrca
        if st.is_leaf[lvl]:
            # a group declared at a species (in-paralogs of one species): transparent, its genes
            # are spliced into the enclosing group (our leaf levels are implicit singletons)
            stats["species_level_groups"] += 1
            return _H(level=lvl, children=items)
        # same-level (or species-level) child groups: splice
        spliced = []
        for it in items:
            if isinstance(it, _H) and (it.level == lvl or st.is_leaf[it.level]):
                stats["same_level_collapsed"] += 1
                for c in it.children:
                    if c.dup is None:
                        c.dup = it.dup
                    spliced.append(c)
            else:
                spliced.append(it)
        items = spliced
        node = _H(level=lvl, children=items)
        place_duplications(node)
        fill_levels(node)
        return node

    def place_duplications(node: _H) -> None:
        by_dup: dict[object, list] = defaultdict(list)
        for it in node.children:
            if it.dup is not None:
                by_dup[it.dup].append(it)
        for key, copies in by_dup.items():
            m = st.lca_index(c.level for c in copies)
            M = int(st.parent_idx[m])
            if M < 0 or not st.anc[node.level, M]:
                msg = (f"duplication with copies at {sorted({st.names[c.level] for c in copies})} predates "
                       f"its enclosing group at {st.names[node.level]!r}; placed just below the group")
                if strict:
                    raise PredictionError(msg)
                warnings.append(msg)
                M = node.level
            if M != node.level:
                ind = _H(level=M, children=copies)
                node.children = [c for c in node.children if c.dup is not key] + [ind]
                fill_levels(ind)

    def child_toward(parent_level: int, lvl: int) -> int:
        """The child of ``parent_level`` on the path down to ``lvl``."""
        cur = lvl
        while st.parent_idx[cur] != parent_level:
            cur = int(st.parent_idx[cur])
        return cur

    def chain(item, parent_level: int):
        """Wrap ``item`` in single-lineage induced HOGs up to (excluding) ``parent_level``."""
        cur = item
        p = int(st.parent_idx[cur.level])
        while p != parent_level:
            cur = _H(level=p, children=[cur])
            p = int(st.parent_idx[p])
        return cur

    def fill_levels(node: _H) -> None:
        X = node.level
        new_children = []
        groups: dict[int, list] = defaultdict(list)
        for it in node.children:
            if int(st.parent_idx[it.level]) == X:
                new_children.append(it)
            elif it.dup is not None:
                top = chain(it, X)
                top.dup = it.dup
                new_children.append(top)
            else:
                groups[child_toward(X, it.level)].append(it)
        for C, members in groups.items():
            if len(members) == 1:
                new_children.append(chain(members[0], X))
            else:
                shared = _H(level=C, children=members)
                fill_levels(shared)
                new_children.append(shared)
        node.children = new_children

    def emit(root: _H) -> None:
        """Assign ids to every node and record memberships; returns via pred."""
        stack = [(root, False)]
        genes_below: dict[int, list] = {}
        order = []
        while stack:  # postorder
            nd, done = stack.pop()
            if isinstance(nd, _G):
                continue
            if done:
                order.append(nd)
            else:
                stack.append((nd, True))
                stack.extend((c, False) for c in nd.children if isinstance(c, _H))
        for nd in order:
            gl: list = []
            for c in nd.children:
                if isinstance(c, _G):
                    gl.append(c.gene)
                else:
                    gl.extend(genes_below.pop(id(c)))
            genes_below[id(nd)] = gl
            hid = f"H{next(hog_counter)}"
            lvname = st.names[nd.level]
            for g in gl:
                pred.add(lvname, hid, g, strict=False)
            stats["hogs"] += 1

    # ---- streaming parse ----
    context = etree.iterparse(str(path), events=("start", "end"), huge_tree=True, remove_comments=True)
    taxon_stack: list = []
    for event, el in context:
        t = _local(el.tag)
        if event == "start":
            if t == "species":
                name = el.get("name")
                current_species = species_map.get(name, name)
                if current_species not in st.index or not st.is_leaf[st.index[current_species]]:
                    unknown_species.add(name)
            elif t == "taxon":
                taxon_names[el.get("id")] = el.get("name")
            continue
        # end events
        if t == "gene":
            if current_species is None:
                continue
            raw = el.get(id_attr) or el.get("protId") or el.get("geneId") or el.get("id")
            if current_species in st.index and st.is_leaf[st.index[current_species]]:
                gid = norm(raw)
                kept = _keep(keep_genes, gid)
                gene_info[el.get("id")] = (kept, st.index[current_species])
                if kept is not None:
                    if kept in pred.gene_species and pred.gene_species[kept] != current_species:
                        warnings.append(f"gene ID {kept!r} occurs in two species")
                    pred.gene_species[kept] = current_species
            el.clear()
        elif t == "species":
            current_species = None
            el.clear()
        elif t in ("orthologGroup", "paralogGroup"):
            parent = el.getparent()
            if parent is not None and _local(parent.tag) == "groups":
                if unknown_species and not ignore_unknown_species:
                    raise PredictionError(
                        f"{path}: species not in the species tree: {sorted(unknown_species)[:10]} "
                        f"(use species_map or ignore_unknown_species=True)")
                roots = [el] if t == "orthologGroup" else [c for c in el if _local(c.tag) == "orthologGroup"]
                for r in roots:
                    h = parse_group(r)
                    if h is not None and not st.is_leaf[h.level]:
                        emit(h)
                el.clear()
                while el.getprevious() is not None:
                    del parent[0]
    pred.meta.update({"warnings": warnings, "stats": dict(stats),
                      "unknown_species": sorted(unknown_species)})
    return pred


# ========================================================================= OrthoFinder
def _find_of_files(path: Path) -> tuple[Path, Optional[Path]]:
    """Locate the HOG directory and the rooted species tree with node labels."""
    p = Path(path)
    if p.is_file():
        hog_dir = p.parent
    elif (p / "Phylogenetic_Hierarchical_Orthogroups").is_dir():
        hog_dir = p / "Phylogenetic_Hierarchical_Orthogroups"
    else:
        hog_dir = p
    tree = None
    for cand in (hog_dir.parent / "Species_Tree" / "SpeciesTree_rooted_node_labels.txt",
                 hog_dir / "SpeciesTree_rooted_node_labels.txt"):
        if cand.exists():
            tree = cand
            break
    return hog_dir, tree


def read_orthofinder(
    path: str | Path,
    st: SpeciesTree,
    of_species_tree: Optional[str | Path] = None,
    normalize: str | Normalizer | None = "auto",
    species_map: Optional[dict[str, str]] = None,
    keep_genes: Optional[set[str]] = None,
    induce_missing_levels: bool = True,
    orthogroups: Optional[str | Path] = None,
) -> Prediction:
    """Parse OrthoFinder (v2.5/3.x) ``Phylogenetic_Hierarchical_Orthogroups/N*.tsv``.

    ``path`` may be the results directory, the HOG directory, or one ``N*.tsv`` file (then the
    other files in that directory are read too). Each ``Nk.tsv`` lists the HOGs at OrthoFinder
    species-tree node ``Nk`` (columns ``HOG``, ``OG``, ``Gene Tree Parent Clade``, then one
    column per species with ``", "``-separated gene IDs). ``Nk`` is mapped to *our* level whose
    clade has exactly the same species set, using OrthoFinder's
    ``Species_Tree/SpeciesTree_rooted_node_labels.txt`` (or ``of_species_tree``).

    OrthoFinder 3.x writes no file for the root node (``N0``); root-level groups are then the
    orthogroups of ``Orthogroups/Orthogroups.tsv`` (or ``orthogroups=``; columns
    ``Orthogroup`` then one per species). We check that every ``N*`` HOG is contained in the
    orthogroup named in its ``OG`` column; violations are counted in
    ``meta['n_hogs_not_nested_in_og']``.

    Our levels without an exactly matching OrthoFinder clade (e.g. OrthoFinder's tree was
    rooted/resolved differently) are, if ``induce_missing_levels``, induced from the HOGs at
    the smallest OrthoFinder clade containing them, restricted to genes below the level (i.e.
    we assume no duplication between the two levels, as for skipped levels in OrthoXML). This
    is recorded in ``meta['induced_levels']``.
    """
    norm = make_normalizer(normalize)
    species_map = species_map or {}
    hog_dir, tree_path = _find_of_files(Path(path))
    if of_species_tree is not None:
        tree_path = Path(of_species_tree)
    if tree_path is None or not Path(tree_path).exists():
        raise PredictionError(f"OrthoFinder species tree with node labels not found near {hog_dir}; "
                              f"pass of_species_tree=")
    of_tree = read_newick(Path(tree_path).read_text())
    og_path = Path(orthogroups) if orthogroups else hog_dir.parent / "Orthogroups" / "Orthogroups.tsv"

    def sp(name: str) -> str:
        return species_map.get(name, name)

    of_clades: dict[str, frozenset[str]] = {}
    for nd in of_tree.postorder():
        if not nd.is_leaf:
            of_clades[nd.name] = frozenset(sp(x) for x in nd.leaf_names())
    unknown = sorted({sp(x) for x in of_tree.leaf_names()} - set(st.species))
    if unknown:
        raise PredictionError(f"OrthoFinder species not in our species tree: {unknown[:10]} (use species_map)")
    of_to_ours: dict[str, str] = {}
    unmatched_of = []
    for n, clade in of_clades.items():
        ours = st.find_clade(clade)
        if ours is not None:
            of_to_ours[n] = ours
        else:
            unmatched_of.append(n)

    files = sorted(hog_dir.glob("N*.tsv"), key=lambda p: (len(p.stem), p.stem))
    files = [f for f in files if re.fullmatch(r"N\d+", f.stem)]
    if not files and not og_path.exists():
        raise PredictionError(f"no N*.tsv files in {hog_dir} and no {og_path}")
    pred = Prediction(meta={"format": "orthofinder", "source": str(hog_dir), "of_to_ours": of_to_ours,
                            "unmatched_of_nodes": unmatched_of})

    def read_table(f: Path, first_species_col: int):
        """Rows as (group id, {species: [genes]}) and, for HOG files, the OG column."""
        rows, og_of = [], {}
        with open(f, newline="") as fh:
            rd = csv.reader(fh, delimiter="\t")
            header = next(rd)
            sp_cols = header[first_species_col:]
            for row in rd:
                if not row:
                    continue
                d: dict[str, list[str]] = {}
                for col, cell in zip(sp_cols, row[first_species_col:]):
                    cell = cell.strip()
                    if cell:
                        d[sp(col)] = [norm(x) for x in cell.split(", ") if x.strip()]
                rows.append((row[0], d))
                if first_species_col == 3 and len(row) > 1:
                    og_of[row[0]] = row[1]
        return rows, og_of

    of_hogs: dict[str, list[tuple[str, dict[str, list[str]]]]] = {}
    og_column: dict[str, dict[str, str]] = {}
    for f in files:
        if f.stem in of_clades:
            of_hogs[f.stem], og_column[f.stem] = read_table(f, 3)
    root_name = of_tree.name
    if root_name not in of_hogs:
        if not og_path.exists():
            raise PredictionError(f"no {root_name}.tsv in {hog_dir} and no orthogroups file {og_path} "
                                  f"for the root level")
        rows, _ = read_table(og_path, 1)
        of_hogs[root_name] = rows
        pred.meta["root_from_orthogroups"] = str(og_path)
        og_genes = {og: {g for gs in d.values() for g in gs} for og, d in rows}
        bad = 0
        for node, rws in of_hogs.items():
            if node == root_name:
                continue
            for hog, d in rws:
                og = og_column[node].get(hog)
                if og not in og_genes or not {g for gs in d.values() for g in gs} <= og_genes[og]:
                    bad += 1
        pred.meta["n_hogs_not_nested_in_og"] = bad

    def add_hogs(level: str, rows, restrict: Optional[frozenset[str]] = None, prefix: str = ""):
        for hog, d in rows:
            hid = prefix + hog
            for s, genes in d.items():
                if restrict is not None and s not in restrict:
                    continue
                for g in genes:
                    kept = _keep(keep_genes, g)
                    pred.add(level, hid, kept, strict=False)
                    if kept is not None:
                        pred.gene_species[kept] = s

    for n, rows in of_hogs.items():
        if n in of_to_ours:
            add_hogs(of_to_ours[n], rows)
    induced = {}
    if induce_missing_levels:
        mapped = {of_to_ours[n] for n in of_hogs if n in of_to_ours}
        for lv in st.levels:
            if lv in mapped:
                continue
            clade = st.clade_species(lv)
            cands = [n for n in of_hogs if of_clades[n] >= clade]
            if not cands:
                continue
            best = min(cands, key=lambda n: len(of_clades[n]))
            induced[lv] = best
            add_hogs(lv, of_hogs[best], restrict=clade, prefix=f"{lv}<")
    pred.meta["induced_levels"] = induced
    return pred


# ================================================================================ TSV
def read_tsv(
    path: str | Path,
    st: Optional[SpeciesTree] = None,
    normalize: str | Normalizer | None = "identity",
    keep_genes: Optional[set[str]] = None,
    gene_species: Optional[dict[str, str]] = None,
) -> Prediction:
    """Read the simple baseline format: tab-separated ``level, hog_id, gene[, species]``
    with an optional header line starting with ``level``. HOG ids only need to be unique
    within a level."""
    norm = make_normalizer(normalize)
    pred = Prediction(meta={"format": "tsv", "source": str(path)})
    with open(path, newline="") as fh:
        for i, row in enumerate(csv.reader(fh, delimiter="\t")):
            if not row or row[0].startswith("#"):
                continue
            if i == 0 and row[0].strip().lower() == "level":
                continue
            if len(row) < 3:
                raise PredictionError(f"{path}:{i + 1}: expected level, hog_id, gene[, species]")
            level, hog, gene = row[0], row[1], norm(row[2])
            if st is not None and level not in st:
                raise PredictionError(f"{path}:{i + 1}: unknown level {level!r}")
            kept = _keep(keep_genes, gene)
            pred.add(level, hog, kept, strict=False)
            if kept is not None and len(row) > 3 and row[3]:
                pred.gene_species[kept] = row[3]
    if gene_species:
        for g in pred.genes():
            pred.gene_species.setdefault(g, gene_species.get(g, ""))
    return pred


def read_prediction(fmt: str, path: str | Path, st: SpeciesTree, **kw) -> Prediction:
    fmt = fmt.lower()
    if fmt == "orthoxml":
        return read_orthoxml(path, st, **kw)
    if fmt == "orthofinder":
        return read_orthofinder(path, st, **kw)
    if fmt == "tsv":
        return read_tsv(path, st, **kw)
    raise PredictionError(f"unknown prediction format {fmt!r} (orthoxml|orthofinder|tsv)")
