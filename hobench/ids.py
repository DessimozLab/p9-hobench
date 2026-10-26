"""Gene-ID normalisation hooks.

Tool outputs often carry decorated identifiers: full FASTA headers, UniProt-style
``sp|P12345|NAME_HUMAN`` / ``tr|A0A0...|...`` strings, or versioned accessions. A
*normaliser* is any ``Callable[[str], str]`` mapping a raw identifier to the canonical gene
ID used in the reference families (``members.tsv: gene_id``).

``make_normalizer(spec)`` builds one from a short spec string:

* ``identity``       — unchanged
* ``first_token``    — first whitespace-separated token (FASTA header → ID)
* ``uniprot``        — first token, then ``db|ACC|NAME`` → ``ACC``
* ``auto``           — same as ``uniprot`` (safe on plain accessions too)
* ``map:<path>``     — two-column TSV raw→canonical (raw or its first token); an ID absent from
  the map raises ``KeyError`` (no silent fallback)
* ``map_or_auto:<path>`` — as ``map:`` but IDs absent from the map explicitly go through ``auto``
* ``strip_version``  can be appended with ``+`` (e.g. ``auto+strip_version``): drop a trailing ``.N``
"""
from __future__ import annotations

import csv
import re
from pathlib import Path
from typing import Callable

Normalizer = Callable[[str], str]

_UNIPROT_RE = re.compile(r"^(?:sp|tr)\|([^|]+)\|")
_VERSION_RE = re.compile(r"\.\d+$")


def identity(x: str) -> str:
    return x


def first_token(x: str) -> str:
    x = x.strip()
    return x.split()[0] if x else x


def uniprot(x: str) -> str:
    x = first_token(x)
    m = _UNIPROT_RE.match(x)
    return m.group(1) if m else x


def strip_version(x: str) -> str:
    return _VERSION_RE.sub("", x)


_BUILTIN: dict[str, Normalizer] = {
    "identity": identity,
    "first_token": first_token,
    "uniprot": uniprot,
    "auto": uniprot,
    "strip_version": strip_version,
}


def load_id_map(path: str | Path) -> dict[str, str]:
    out: dict[str, str] = {}
    with open(path, newline="") as fh:
        for row in csv.reader(fh, delimiter="\t"):
            if not row or row[0].startswith("#"):
                continue
            if len(row) < 2:
                raise ValueError(f"{path}: ID map rows need two columns (raw, canonical): {row}")
            out[row[0]] = row[1]
    return out


def make_normalizer(spec: str | Normalizer | None = "auto") -> Normalizer:
    if spec is None:
        return uniprot
    if callable(spec):
        return spec
    funcs: list[Normalizer] = []
    for part in spec.split("+"):
        if part.startswith("map:") or part.startswith("map_or_auto:"):
            fallback = part.startswith("map_or_auto:")
            mapping = load_id_map(part.split(":", 1)[1])

            def _m(x: str, _mapping=mapping, _fb=fallback) -> str:
                v = _mapping.get(x) or _mapping.get(first_token(x))
                if v is None:
                    if not _fb:
                        raise KeyError(f"ID {x!r} not in the ID map (use map_or_auto: to pass unmapped IDs through)")
                    v = uniprot(x)
                return v

            funcs.append(_m)
        elif part in _BUILTIN:
            funcs.append(_BUILTIN[part])
        else:
            raise ValueError(f"unknown ID normaliser {part!r}; choose from {sorted(_BUILTIN)} or map:<tsv>")

    def norm(x: str) -> str:
        for f in funcs:
            x = f(x)
        return x

    return norm
