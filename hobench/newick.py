"""Minimal, dependency-free Newick / NHX reader and writer.

We use our own parser (rather than ete3) so that behaviour is fully specified:

* labels may be quoted (``'a b'``, with ``''`` as an escaped quote) or unquoted;
  unquoted labels are taken literally (underscores are *not* turned into spaces);
* branch lengths (``:0.1``) are parsed into ``Node.length``;
* bracketed comments are ignored except NHX comments ``[&&NHX:key=value:key=value]``,
  whose key/value pairs are stored in ``Node.features`` (values are strings);
* internal node labels are kept as ``Node.name`` (OrthoFinder ``N0``, taxon names, ...).
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterator, Optional

_SPECIAL = set("()[]':;,")


@dataclass(eq=False)
class Node:
    name: str = ""
    length: Optional[float] = None
    children: list["Node"] = field(default_factory=list)
    parent: Optional["Node"] = None
    features: dict[str, str] = field(default_factory=dict)

    # -- structure ---------------------------------------------------------------
    def add_child(self, child: "Node") -> "Node":
        child.parent = self
        self.children.append(child)
        return child

    @property
    def is_leaf(self) -> bool:
        return not self.children

    @property
    def is_root(self) -> bool:
        return self.parent is None

    def preorder(self) -> Iterator["Node"]:
        stack = [self]
        while stack:
            n = stack.pop()
            yield n
            stack.extend(reversed(n.children))

    def postorder(self) -> Iterator["Node"]:
        stack: list[tuple[Node, bool]] = [(self, False)]
        while stack:
            n, done = stack.pop()
            if done or not n.children:
                yield n
            else:
                stack.append((n, True))
                stack.extend((c, False) for c in reversed(n.children))

    def leaves(self) -> list["Node"]:
        return [n for n in self.preorder() if n.is_leaf]

    def leaf_names(self) -> list[str]:
        return [n.name for n in self.leaves()]

    def __repr__(self) -> str:  # pragma: no cover - debugging aid
        return f"Node({self.name!r}, n_children={len(self.children)})"


class NewickError(ValueError):
    pass


def _tokenize(text: str) -> Iterator[tuple[str, str]]:
    """Yield (kind, value) tokens: kind in {'(', ')', ',', ':', ';', 'label', 'comment'}."""
    i, n = 0, len(text)
    while i < n:
        c = text[i]
        if c.isspace():
            i += 1
        elif c in "(),:;":
            yield c, c
            i += 1
        elif c == "[":
            j = text.find("]", i)
            if j < 0:
                raise NewickError("unterminated '[' comment")
            yield "comment", text[i + 1 : j]
            i = j + 1
        elif c == "'":
            buf = []
            i += 1
            while True:
                if i >= n:
                    raise NewickError("unterminated quoted label")
                if text[i] == "'":
                    if i + 1 < n and text[i + 1] == "'":
                        buf.append("'")
                        i += 2
                        continue
                    i += 1
                    break
                buf.append(text[i])
                i += 1
            yield "label", "".join(buf)
        else:
            j = i
            while j < n and text[j] not in _SPECIAL and not text[j].isspace():
                j += 1
            yield "label", text[i:j]
            i = j


def _parse_nhx(comment: str) -> Optional[dict[str, str]]:
    if not comment.startswith("&&NHX"):
        return None
    out: dict[str, str] = {}
    for part in comment[len("&&NHX") :].split(":"):
        if not part:
            continue
        if "=" not in part:
            raise NewickError(f"malformed NHX field {part!r} in [{comment}]")
        k, v = part.split("=", 1)
        out[k] = v
    return out


def parse_newick(text: str) -> Node:
    """Parse a single Newick/NHX tree from a string."""
    tokens = list(_tokenize(text.strip()))
    if not tokens:
        raise NewickError("empty Newick string")
    pos = 0

    def peek() -> tuple[str, str]:
        return tokens[pos] if pos < len(tokens) else (";", ";")

    def parse_annotations(node: Node) -> None:
        nonlocal pos
        # label, then any mix of :length and [comments]
        if peek()[0] == "label":
            node.name = tokens[pos][1]
            pos += 1
        while True:
            kind = peek()[0]
            if kind == "comment":
                feats = _parse_nhx(tokens[pos][1])
                if feats:
                    node.features.update(feats)
                pos += 1
            elif kind == ":":
                pos += 1
                if peek()[0] != "label":
                    raise NewickError("expected branch length after ':'")
                try:
                    node.length = float(tokens[pos][1])
                except ValueError as e:
                    raise NewickError(f"bad branch length {tokens[pos][1]!r}") from e
                pos += 1
            else:
                break

    def parse_subtree() -> Node:
        nonlocal pos
        node = Node()
        if peek()[0] == "(":
            pos += 1
            while True:
                node.add_child(parse_subtree())
                kind = peek()[0]
                if kind == ",":
                    pos += 1
                    continue
                if kind == ")":
                    pos += 1
                    break
                raise NewickError(f"unexpected token {peek()[1]!r} inside parentheses")
        parse_annotations(node)
        return node

    root = parse_subtree()
    if peek()[0] != ";" and pos < len(tokens):
        raise NewickError(f"unexpected trailing token {peek()[1]!r}")
    return root


def read_newick(path_or_text: str | Path) -> Node:
    """Read a tree from a file path, or parse it directly if given Newick text."""
    p = Path(str(path_or_text))
    s = str(path_or_text).strip()
    if not s.startswith("(") and p.exists():
        s = p.read_text()
    return parse_newick(s)


def _quote(label: str) -> str:
    if label == "" or not any(c in _SPECIAL or c.isspace() for c in label):
        return label
    return "'" + label.replace("'", "''") + "'"


def write_newick(root: Node, features: Optional[list[str]] = None, lengths: bool = True) -> str:
    """Serialise a tree. If ``features`` is given, write those keys as NHX (in that order,
    followed by any other keys alphabetically when ``features == []`` means 'all')."""

    def fmt(n: Node) -> str:
        s = ""
        if n.children:
            s = "(" + ",".join(fmt(c) for c in n.children) + ")"
        s += _quote(n.name)
        if lengths and n.length is not None:
            s += f":{n.length:g}"
        if features is not None and n.features:
            keys = features if features else sorted(n.features)
            parts = [f"{k}={n.features[k]}" for k in keys if k in n.features and n.features[k] != ""]
            if parts:
                s += "[&&NHX:" + ":".join(parts) + "]"
        return s

    return fmt(root) + ";"


def copy_tree(root: Node) -> Node:
    """Deep copy of a tree (names, lengths, features)."""
    new_root = Node(root.name, root.length, features=dict(root.features))
    stack = [(root, new_root)]
    while stack:
        old, new = stack.pop()
        for c in old.children:
            nc = new.add_child(Node(c.name, c.length, features=dict(c.features)))
            stack.append((c, nc))
    return new_root
