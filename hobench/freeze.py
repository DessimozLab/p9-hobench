"""Freeze reference families: ``python -m hobench.freeze [--families families] [--check]``.

Writes ``<families>/MANIFEST.sha256`` (``sha256  relative/path`` for every file under the
families directory except the manifest itself and caches), or with ``--check`` verifies the
current files against it (exit status 1 and a list of differences on mismatch).
"""
from __future__ import annotations

import argparse
import hashlib
import sys
from pathlib import Path

from .panel import sha256_file

MANIFEST = "MANIFEST.sha256"
_SKIP_PARTS = {"__pycache__", ".ipynb_checkpoints"}


def manifest_entries(families_dir: str | Path) -> list[tuple[str, str]]:
    root = Path(families_dir)
    out = []
    for p in sorted(root.rglob("*")):
        if not p.is_file() or p.name == MANIFEST or _SKIP_PARTS & set(p.parts):
            continue
        out.append((sha256_file(p), p.relative_to(root).as_posix()))
    return out


def write_manifest(families_dir: str | Path) -> Path:
    root = Path(families_dir)
    path = root / MANIFEST
    path.write_text("".join(f"{h}  {rel}\n" for h, rel in manifest_entries(root)))
    return path


def check_manifest(families_dir: str | Path) -> list[str]:
    """Differences between the manifest and the files on disk (empty list = frozen state)."""
    root = Path(families_dir)
    path = root / MANIFEST
    if not path.exists():
        return [f"{path} does not exist"]
    want = {}
    for line in path.read_text().splitlines():
        if line.strip():
            h, rel = line.split("  ", 1)
            want[rel] = h
    have = {rel: h for h, rel in manifest_entries(root)}
    diffs = [f"changed: {k}" for k in sorted(set(want) & set(have)) if want[k] != have[k]]
    diffs += [f"missing: {k}" for k in sorted(set(want) - set(have))]
    diffs += [f"new (not in manifest): {k}" for k in sorted(set(have) - set(want))]
    return diffs


def tree_hash(entries: list[tuple[str, str]]) -> str:
    """One hash over (relative path, file hash) pairs."""
    h = hashlib.sha256()
    for fh, rel in sorted(entries, key=lambda x: x[1]):
        h.update(f"{rel}\0{fh}\n".encode())
    return h.hexdigest()


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="python -m hobench.freeze", description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--families", default="families")
    ap.add_argument("--check", action="store_true")
    a = ap.parse_args(argv)
    if a.check:
        diffs = check_manifest(a.families)
        for d in diffs:
            print(d, file=sys.stderr)
        print("OK" if not diffs else f"{len(diffs)} difference(s)", file=sys.stderr)
        return 1 if diffs else 0
    p = write_manifest(a.families)
    print(f"wrote {p} ({len(p.read_text().splitlines())} files)", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
