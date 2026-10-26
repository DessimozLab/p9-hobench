"""Normalised predicted HOGs and parsers for tool outputs.

A :class:`Prediction` stores, for each *internal* level ``L`` of the panel species tree, a
mapping ``gene -> HOG id at L`` (only for genes of species below ``L``), plus the total size of
every HOG. Leaf (species) levels are implicit: every gene is its own HOG there.

To keep memory bounded on whole-proteome outputs, parsers accept ``keep_genes``: only those
genes' memberships are stored, but HOG sizes always count *all* genes; per-HOG gene counts
are also recorded for the species in ``count_species`` (declared-complete species, needed for
membership precision). The evaluator passes the union of all reference members and decoys.

Parsers (all strict by default: inconsistencies raise :class:`PredictionError`)
-------
* :func:`read_orthoxml` — OrthoXML 0.3/0.5 (FastOMA, OMA). Levels are mapped by clade
  species set; duplication placement ``declared`` (default) or ``late`` (pyHam).
* :func:`read_orthofinder` — OrthoFinder ``Phylogenetic_Hierarchical_Orthogroups/N*.tsv``
  (root level from ``Orthogroups/Orthogroups.tsv`` when there is no ``N0.tsv``).
* :func:`read_tsv` — simple ``level, hog_id, gene[, species]`` table.
"""
from __future__ import annotations

import csv
import gzip
import itertools
import re
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, Optional

from .ids import Normalizer, make_normalizer
from .newick import Node, read_newick
from .speciestree import SpeciesTree

DUP_PLACEMENTS = ("declared", "late")


class PredictionError(ValueError):
    pass


@dataclass
class Prediction:
    levels: dict[str, dict[str, str]] = field(default_factory=dict)   # level -> gene -> hog
    sizes: dict[str, dict[str, int]] = field(default_factory=dict)    # level -> hog -> n genes (all)
    gene_species: dict[str, str] = field(default_factory=dict)        # kept genes only
    meta: dict = field(default_factory=dict)
    count_species: set[str] = field(default_factory=set)              # species counted per HOG
    hog_species: dict[str, dict[str, dict[str, int]]] = field(default_factory=dict)  # level -> hog -> sp -> n
    _members: Optional[dict[str, dict[str, list[str]]]] = field(default=None, repr=False)

    def add(self, level: str, hog: str, gene: Optional[str], strict: bool = True,
            species: Optional[str] = None) -> None:
        """Register a gene in a HOG; ``gene=None`` only increments the HOG size (gene not kept).
        ``species`` (the gene's species) is needed for species in ``count_species``.

        A gene may belong to only one HOG per level: a second, different assignment raises
        (``strict``) or keeps the first one and is counted in ``meta['membership_conflicts']``."""
        if gene is not None:
            lv = self.levels.setdefault(level, {})
            if gene in lv and lv[gene] != hog:
                if strict:
                    raise PredictionError(f"gene {gene!r} assigned to two HOGs at level {level!r}: "
                                          f"{lv[gene]!r} and {hog!r}")
                self.meta["membership_conflicts"] = self.meta.get("membership_conflicts", 0) + 1
                return
            lv[gene] = hog
        self.sizes.setdefault(level, {})
        self.sizes[level][hog] = self.sizes[level].get(hog, 0) + 1
        if species is not None and species in self.count_species:
            d = self.hog_species.setdefault(level, {}).setdefault(hog, {})
            d[species] = d.get(species, 0) + 1
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
        """Return problems: unknown levels, genes at a level their species is not below
        (both fatal for scoring), and non-nested HOGs (informational; scored as given)."""
        return self.fatal_problems(st) + self.nesting_problems(st)

    def fatal_problems(self, st: SpeciesTree) -> list[str]:
        problems = []
        for lv, m in self.levels.items():
            if lv not in st or st.is_species(lv):
                problems.append(f"level {lv!r} is not an internal level of the species tree")
                continue
            clade = st.clade_species(lv)
            bad = [g for g in m if g in self.gene_species and self.gene_species[g] not in clade]
            if bad:
                problems.append(f"level {lv!r}: {len(bad)} gene(s) whose species is not below it, e.g. {bad[:3]}")
        return problems

    def nesting_problems(self, st: SpeciesTree) -> list[str]:
        problems = []
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
    species: str
    dup: object = None


def _clade_level_map(tool_root: Node, st: SpeciesTree, species_map: dict[str, str]) -> dict[str, int]:
    """Map every node name of a tool's species tree/taxonomy to our level with the same clade
    species set (names of nodes with no matching clade are absent)."""
    out: dict[str, int] = {}
    sets: dict[int, frozenset] = {}
    for nd in tool_root.postorder():
        if nd.is_leaf:
            sets[id(nd)] = frozenset([species_map.get(nd.name, nd.name)])
        else:
            sets[id(nd)] = frozenset().union(*(sets[id(c)] for c in nd.children))
        if nd.name:
            ours = st.find_clade(sets[id(nd)]) if sets[id(nd)] <= set(st.species) else None
            if ours is not None:
                out[nd.name] = st.idx(ours)
    return out


def read_orthoxml(
    path: str | Path,
    st: SpeciesTree,
    normalize: str | Normalizer | None = "auto",
    id_attr: str = "protId",
    species_map: Optional[dict[str, str]] = None,
    level_map: Optional[dict[str, str]] = None,
    keep_genes: Optional[set[str]] = None,
    count_species: Optional[set[str]] = None,
    tool_species_tree: Optional[str | Path] = None,
    dup_placement: str = "declared",
    strict: bool = True,
    ignore_unknown_species: bool = False,
) -> Prediction:
    """Parse OrthoXML (0.3 / 0.5) into a :class:`Prediction`, inducing HOGs at skipped levels.

    **Level of a group (E2).** Level *names* in the file (``taxonId`` → ``<taxonomy>`` name, or
    the ``TaxRange`` property) are mapped to our levels **by clade species set**, using the
    tool's own species tree: ``tool_species_tree`` (e.g. FastOMA's
    ``species_tree_checked.nwk``) if given, else the file's ``<taxonomy>``, else our tree (by
    name). ``level_map`` renames before lookup. Any name that cannot be mapped is a hard error.
    A group with no level annotation gets the MRCA of its children (pyHam).

    **Duplication placement (E1).** ``declared`` (default): a ``paralogGroup`` inside a group
    declared at level ``X`` is a duplication just below ``X``; each copy gets its own chain of
    induced HOGs below ``X``. ``late`` (pyHam, sensitivity analysis): the duplication sits at
    ``M = parent(MRCA(copy levels))``; if ``M`` is below ``X`` an induced HOG at ``M`` holds the
    copies. Groups where the two readings differ are counted in
    ``meta['dup_placement_ambiguous']``. Directly nested paralogGroups are one duplication.

    **Skipped levels.** Non-paralog children lying strictly below the same child level ``C``
    of a group share one induced HOG at ``C`` (no duplication was recorded between them);
    children exactly at ``C`` remain separate HOGs (pyHam ``_resolve_missing_levels``). A child
    group at the same level as its parent, or a group declared at a species, is spliced in.

    With ``strict`` (default) inconsistencies raise: a declared level not above its members,
    a duplication older than its group, a geneRef to an unknown gene, a gene in two groups.
    ``strict=False`` repairs them and records warnings (sensitivity use only).
    """
    from lxml import etree

    if dup_placement not in DUP_PLACEMENTS:
        raise PredictionError(f"dup_placement must be one of {DUP_PLACEMENTS}")
    norm = make_normalizer(normalize)
    species_map = species_map or {}
    level_map = level_map or {}
    pred = Prediction(meta={"format": "orthoxml", "source": str(path), "dup_placement": dup_placement},
                      count_species=set(count_species or ()))
    warnings: list[str] = []
    gene_info: dict[str, tuple[Optional[str], int, str]] = {}
    taxon_names: dict[str, str] = {}
    unknown_species: set[str] = set()
    current_species: Optional[str] = None
    hog_counter = itertools.count()
    stats = Counter()
    name_to_level: Optional[dict[str, int]] = None
    if tool_species_tree is not None:
        name_to_level = _clade_level_map(read_newick(Path(tool_species_tree).read_text()), st, species_map)
        pred.meta["level_mapping"] = f"clade sets from {tool_species_tree}"
    taxonomy_root: Optional[Node] = None
    taxon_stack: list[Node] = []

    def problem(msg: str) -> None:
        if strict:
            raise PredictionError(f"{path}: {msg}")
        warnings.append(msg)

    def level_idx(name: Optional[str]) -> Optional[int]:
        if name is None:
            return None
        name = level_map.get(name, name)
        if name_to_level is not None:
            if name in name_to_level:
                return name_to_level[name]
        elif name in st.index:
            return st.index[name]
        raise PredictionError(f"{path}: level {name!r} cannot be mapped to a level of the panel species tree "
                              f"(by clade species set); check the species tree given to the tool")

    def parse_group(el) -> Optional[_H]:
        items: list = []
        props: dict[str, str] = {}

        def gene_item(ch, dupkey):
            gi = gene_info.get(ch.get("id"))
            if gi is None:
                problem(f"geneRef id={ch.get('id')!r} does not refer to a gene of a panel species")
                stats["generef_unknown"] += 1
                return
            items.append(_G(gi[0], gi[1], gi[2], dupkey))

        def collect(e, dupkey):
            for ch in e:
                t = _local(ch.tag)
                if t == "geneRef":
                    gene_item(ch, dupkey)
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
                gene_item(ch, None)
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
            tname = taxon_names.get(el.get("taxonId"))
            if tname is None:
                raise PredictionError(f"{path}: group {el.get('id')!r}: taxonId {el.get('taxonId')!r} "
                                      f"is not in the file's <taxonomy>")
            lvl = level_idx(tname)
        if lvl is None:
            for key in ("TaxRange", "taxid"):
                if props.get(key) is not None:
                    lvl = level_idx(props[key])
                    break
        mrca = st.lca_index(child_levels)
        if lvl is None:
            stats["level_from_mrca"] += 1
            if len(set(child_levels)) == 1:
                lvl = int(st.parent_idx[mrca]) if st.parent_idx[mrca] >= 0 else mrca
            else:
                lvl = mrca
        elif not st.anc[lvl, mrca]:
            problem(f"group {el.get('id')!r}: declared level {st.names[lvl]!r} is not above all its "
                    f"members (their MRCA is {st.names[mrca]!r})")
            lvl = mrca
        if st.is_leaf[lvl]:
            stats["species_level_groups"] += 1
            return _H(level=lvl, children=items)
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
        node = _H(level=lvl, children=spliced)
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
                problem(f"duplication with copies at {sorted({st.names[c.level] for c in copies})} predates "
                        f"its enclosing group at {st.names[node.level]!r}")
                M = node.level
            if M != node.level:
                stats["dup_placement_ambiguous"] += 1
                if dup_placement == "late":
                    ind = _H(level=M, children=copies)
                    node.children = [c for c in node.children if c.dup is not key] + [ind]
                    fill_levels(ind)

    def child_toward(parent_level: int, lvl: int) -> int:
        cur = lvl
        while st.parent_idx[cur] != parent_level:
            cur = int(st.parent_idx[cur])
        return cur

    def chain(item, parent_level: int):
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
        stack = [(root, False)]
        genes_below: dict[int, list] = {}
        order = []
        while stack:  # postorder
            nd, done = stack.pop()
            if done:
                order.append(nd)
            else:
                stack.append((nd, True))
                stack.extend((c, False) for c in nd.children if isinstance(c, _H))
        for nd in order:
            gl: list = []
            for c in nd.children:
                if isinstance(c, _G):
                    gl.append((c.gene, c.species))
                else:
                    gl.extend(genes_below.pop(id(c)))
            genes_below[id(nd)] = gl
            hid = f"H{next(hog_counter)}"
            lvname = st.names[nd.level]
            for g, s in gl:
                pred.add(lvname, hid, g, strict=strict, species=s)
            stats["hogs"] += 1

    # .gz input is decompressed explicitly (lxml's path-based iterparse does not do it; FastOMA output is gzipped)
    src = gzip.open(path, "rb") if str(path).endswith(".gz") else str(path)
    context = etree.iterparse(src, events=("start", "end"), huge_tree=True, remove_comments=True)
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
                nd = Node(el.get("name"))
                if taxon_stack:
                    taxon_stack[-1].add_child(nd)
                taxon_stack.append(nd)
            continue
        if t == "taxon":
            nd = taxon_stack.pop()
            if not taxon_stack:
                taxonomy_root = nd
                if name_to_level is None:
                    name_to_level = _clade_level_map(taxonomy_root, st, species_map)
                    pred.meta["level_mapping"] = "clade sets from the file's <taxonomy>"
        elif t == "gene":
            if current_species is None:
                continue
            raw = el.get(id_attr) or el.get("protId") or el.get("geneId") or el.get("id")
            if current_species in st.index and st.is_leaf[st.index[current_species]]:
                gid = norm(raw)
                kept = _keep(keep_genes, gid)
                gene_info[el.get("id")] = (kept, st.index[current_species], current_species)
                if kept is not None:
                    if kept in pred.gene_species and pred.gene_species[kept] != current_species:
                        problem(f"gene ID {kept!r} occurs in two species")
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
                if name_to_level is None:
                    pred.meta["level_mapping"] = "names of the panel species tree"
                    name_to_level = None  # name lookup in our tree
                roots = [el] if t == "orthologGroup" else [c for c in el if _local(c.tag) == "orthologGroup"]
                for r in roots:
                    h = parse_group(r)
                    if h is not None and not st.is_leaf[h.level]:
                        emit(h)
                el.clear()
                while el.getprevious() is not None:
                    del parent[0]
    pred.meta.update({"warnings": warnings, "stats": dict(stats),
                      "dup_placement_ambiguous": stats.get("dup_placement_ambiguous", 0),
                      "unknown_species": sorted(unknown_species)})
    return pred


# ========================================================================= OrthoFinder
def _find_of_files(path: Path) -> tuple[Path, Optional[Path]]:
    """Locate the HOG directory and the rooted species tree with node labels (the
    orthogroups file is ``<hog_dir>/../Orthogroups/Orthogroups.tsv``)."""
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
    count_species: Optional[set[str]] = None,
    allow_induced_levels: bool = False,
    orthogroups: Optional[str | Path] = None,
) -> Prediction:
    """Parse OrthoFinder (v2.5/3.x) ``Phylogenetic_Hierarchical_Orthogroups/N*.tsv``.

    ``path`` may be the results directory, the HOG directory, or one ``N*.tsv`` file. Each
    ``Nk.tsv`` lists the HOGs at OrthoFinder node ``Nk`` (columns ``HOG``, ``OG``, ``Gene Tree
    Parent Clade``, then one column per species with ``", "``-separated gene IDs). ``Nk`` is
    mapped to *our* level with exactly the same clade species set, using
    ``Species_Tree/SpeciesTree_rooted_node_labels.txt`` (or ``of_species_tree``).

    OrthoFinder 3.x writes no file for the root node; root-level groups are then the
    orthogroups of ``Orthogroups/Orthogroups.tsv`` (or ``orthogroups=``). Every ``N*`` HOG must
    lie inside the OG named in its ``OG`` column (``meta['n_hogs_not_nested_in_og']``).

    **E3.** Every internal level of our tree must match an OrthoFinder clade (run OrthoFinder
    with ``-s <panel tree>``); otherwise this raises. With ``allow_induced_levels`` (perturbed
    species-tree experiments only) unmatched levels are induced from the smallest
    OrthoFinder clade containing them, restricted to their species, and listed in
    ``meta['induced_levels']``.
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
                            "unmatched_of_nodes": unmatched_of}, count_species=set(count_species or ()))

    def read_table(f: Path, first_species_col: int):
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

    of_hogs: dict[str, list] = {}
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
                    pred.add(level, hid, kept, strict=True, species=s)
                    if kept is not None:
                        pred.gene_species[kept] = s

    for n, rows in of_hogs.items():
        if n in of_to_ours:
            add_hogs(of_to_ours[n], rows)
    mapped = {of_to_ours[n] for n in of_hogs if n in of_to_ours}
    missing = [lv for lv in st.levels if lv not in mapped]
    induced = {}
    if missing and not allow_induced_levels:
        raise PredictionError(f"OrthoFinder output has no clade matching our level(s) {missing[:10]}; run "
                              f"OrthoFinder with -s <panel tree> (or allow_induced_levels for perturbed-tree runs)")
    for lv in missing:
        clade = st.clade_species(lv)
        cands = [n for n in of_hogs if of_clades[n] >= clade]
        if not cands:
            raise PredictionError(f"no OrthoFinder clade contains level {lv!r}")
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
    count_species: Optional[set[str]] = None,
    gene_species: Optional[dict[str, str]] = None,
    flat: bool = False,
) -> Prediction:
    """Read the simple baseline format: tab-separated ``level, hog_id, gene[, species]``
    with an optional header line starting with ``level``. HOG ids only need to be unique
    within a level. The species column is required when ``count_species`` is non-empty.

    ``flat``: the file holds a single flat partition (e.g. OrthoFinder orthogroups at the root,
    "of_flat"). Each group is projected to every internal level at or below its row level as
    the group restricted to that clade's species (one HOG per group per level = "no
    duplications" reading), id ``<hog>@<level>``. Requires ``st`` and species."""
    if flat and st is None:
        raise PredictionError("read_tsv(flat=True) needs the species tree")
    norm = make_normalizer(normalize)
    pred = Prediction(meta={"format": "tsv", "source": str(path), "flat": flat},
                      count_species=set(count_species or ()))
    with (gzip.open(path, "rt", newline="") if str(path).endswith(".gz") else open(path, newline="")) as fh:
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
            species = row[3] if len(row) > 3 and row[3] else (gene_species or {}).get(gene)
            if pred.count_species and species is None:
                raise PredictionError(f"{path}:{i + 1}: species column required (complete-species counting)")
            kept = _keep(keep_genes, gene)
            if flat:
                if not species or species not in st.index:
                    raise PredictionError(f"{path}:{i + 1}: flat projection needs a known species for {gene!r}")
                if not st.anc[st.idx(level), st.idx(species)]:
                    raise PredictionError(f"{path}:{i + 1}: species {species!r} is not below level {level!r}")
                for lv in st.ancestors(species, include_self=False):
                    pred.add(lv, hog if lv == level else f"{hog}@{lv}", kept, strict=True, species=species)
                    if lv == level:
                        break
            else:
                pred.add(level, hog, kept, strict=True, species=species)
            if kept is not None and species:
                pred.gene_species[kept] = species
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
