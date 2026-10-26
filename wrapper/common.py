"""Shared helpers for the HOBench tool wrapper: paths, panels, hashing, FASTA I/O, timing.

The wrapper is blind to reference families: nothing here reads ``families/``.
"""
from __future__ import annotations

import gzip
import hashlib
import json
import os
import re
import subprocess
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, Iterator, Optional

from hobench.newick import Node, read_newick, write_newick

ROOT = Path(__file__).resolve().parent.parent
WORK_ROOT = Path(os.environ.get("HOBENCH_WORK", "/tmp/hobench_work"))
CONDA_BIN = ROOT / "conda_envs" / "ortho" / "bin"
LEVEL_RE = re.compile(r"^[A-Za-z0-9_-]{3,}$")  # SYNTHESIS D2; also FastOMA's renaming rule
GLOBAL_SEED = 20261026  # hobench.seeds.SEEDS[0]


class WrapperError(RuntimeError):
    pass


# ----------------------------------------------------------------------------- env / hashing
def tool_env() -> dict[str, str]:
    env = dict(os.environ)
    env["PATH"] = f"{CONDA_BIN}:{env.get('PATH', '')}"
    env.setdefault("NXF_HOME", str(ROOT / "tools" / "nxf_home"))
    env.setdefault("NXF_OFFLINE", "true")
    return env


def sha256_file(path: Path, chunk: int = 1 << 22) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        while True:
            b = fh.read(chunk)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def sha256_cached(path: Path) -> str:
    """sha256 of a big file (OMAmer DB), cached next to it keyed by size+mtime."""
    path = Path(path).resolve()
    st = path.stat()
    key = f"{st.st_size}:{int(st.st_mtime)}"
    cache = path.with_name(path.name + ".sha256.json")
    try:
        d = json.loads(cache.read_text())
        if d.get("key") == key:
            return d["sha256"]
    except (OSError, ValueError, KeyError):
        pass
    digest = sha256_file(path)
    try:
        cache.write_text(json.dumps({"key": key, "sha256": digest}))
    except OSError:
        pass
    return digest


def stable_unit(seed: int, *parts: str) -> float:
    """Deterministic, order-independent pseudo-random number in [0,1) for an item.
    Selecting items with stable_unit < f gives nested samples as f grows."""
    h = hashlib.sha256(("|".join([str(seed), *parts])).encode()).digest()
    return int.from_bytes(h[:8], "big") / 2.0**64


# ----------------------------------------------------------------------------- FASTA
def read_fasta(path: Path) -> Iterator[tuple[str, str]]:
    """Yield (header_without_>, sequence)."""
    opener = gzip.open if str(path).endswith(".gz") else open
    header, seq = None, []
    with opener(path, "rt") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq)
                header, seq = line[1:], []
            elif line:
                seq.append(line.strip())
    if header is not None:
        yield header, "".join(seq)


def fasta_id(header: str) -> str:
    """The identifier tools see: first whitespace token of the header."""
    return header.split()[0] if header.strip() else ""


def accession(raw_id: str) -> str:
    """``sp|ACC|NAME`` -> ``ACC``; plain IDs unchanged."""
    m = re.match(r"^(?:sp|tr)\|([^|]+)\|", raw_id)
    return m.group(1) if m else raw_id


def write_fasta(path: Path, records: Iterable[tuple[str, str]], width: int = 60) -> int:
    n = 0
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as fh:
        for header, seq in records:
            fh.write(f">{header}\n")
            for i in range(0, len(seq), width):
                fh.write(seq[i:i + width] + "\n")
            n += 1
    return n


# ----------------------------------------------------------------------------- trees
def tree_clades(root: Node) -> dict[str, frozenset[str]]:
    """internal-node name -> leaf set."""
    out = {}
    for n in root.postorder():
        if not n.is_leaf:
            out[n.name] = frozenset(n.leaf_names())
    return out


def validate_species_tree(root: Node, codes: Iterable[str], require_binary: bool = True) -> list[str]:
    """Return a list of problems (empty == OK). Checks: leaves == species codes, internal
    names present/unique/regex-safe (so FastOMA does not rename them), binary."""
    problems = []
    codes = set(codes)
    leaves = root.leaf_names()
    if len(leaves) != len(set(leaves)):
        problems.append("duplicate leaf names")
    if set(leaves) != codes:
        problems.append(f"leaves != species codes: only_tree={sorted(set(leaves) - codes)} "
                        f"only_panel={sorted(codes - set(leaves))}")
    names = []
    for n in root.preorder():
        if n.is_leaf:
            continue
        if not n.name or not LEVEL_RE.match(n.name):
            problems.append(f"internal node name {n.name!r} not matching {LEVEL_RE.pattern}")
        names.append(n.name)
        if require_binary and len(n.children) != 2:
            problems.append(f"node {n.name!r} has {len(n.children)} children (binary tree required)")
    dup = {x for x in names if names.count(x) > 1}
    if dup:
        problems.append(f"duplicate internal names {sorted(dup)}")
    return problems


def prune_tree(root: Node, keep: set[str]) -> Node:
    """Return a copy restricted to leaves in ``keep``; unary nodes are removed (the child keeps
    its own name), so the surviving levels are exactly those whose clade is still split."""
    def rec(n: Node) -> Optional[Node]:
        if n.is_leaf:
            return Node(n.name, n.length) if n.name in keep else None
        kids = [k for k in (rec(c) for c in n.children) if k is not None]
        if not kids:
            return None
        if len(kids) == 1:
            k = kids[0]
            if k.length is not None and n.length is not None:
                k.length += n.length
            return k
        new = Node(n.name, n.length)
        for k in kids:
            new.add_child(k)
        return new
    out = rec(root)
    if out is None:
        raise WrapperError("pruning removed every leaf")
    out.length = None
    return out


# ----------------------------------------------------------------------------- panels
@dataclass
class Panel:
    name: str
    dir: Path
    species: list[dict[str, str]]  # rows of species.tsv (must contain 'code')
    fastas: dict[str, Path]  # code -> FASTA (mode-dependent)
    tree_path: Path
    mode: str = "full"
    notes: list[str] = field(default_factory=list)

    @property
    def codes(self) -> list[str]:
        return [r["code"] for r in self.species]


def read_species_tsv(path: Path) -> list[dict[str, str]]:
    rows = []
    with open(path) as fh:
        header = None
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if header is None:
                header = [h.strip().lower() for h in parts]
                if "code" not in header:
                    raise WrapperError(f"{path}: header lacks a 'code' column: {header}")
                continue
            rows.append(dict(zip(header, parts)))
    return rows


def panel_dir(panel: str, panels_root: Optional[Path] = None) -> Path:
    return (panels_root or ROOT / "panels") / panel


def load_panel(panel: str, mode: str = "full", panels_root: Optional[Path] = None,
               species_tree: Optional[Path] = None) -> Panel:
    d = panel_dir(panel, panels_root)
    sp_tsv = d / "species.tsv"
    if not sp_tsv.exists():
        raise WrapperError(f"{sp_tsv} missing")
    species = read_species_tsv(sp_tsv)
    sub = {"full": "proteomes", "subset": "subset"}[mode]
    fastas = {}
    for r in species:
        c = r["code"]
        cands = [d / sub / f"{c}.fa", d / sub / f"{c}.fasta"]
        p = next((x for x in cands if x.exists()), None)
        if p is None:
            raise WrapperError(f"no FASTA for {c} in {d / sub}")
        fastas[c] = p
    tree = Path(species_tree) if species_tree else d / "species_tree.nwk"
    if not tree.exists():
        raise WrapperError(f"species tree {tree} missing")
    return Panel(panel, d, species, fastas, tree, mode)


# ----------------------------------------------------------------------------- timing
_TIME_KEYS = {
    "Elapsed (wall clock) time (h:mm:ss or m:ss)": "wall_clock",
    "User time (seconds)": "user_s",
    "System time (seconds)": "sys_s",
    "Maximum resident set size (kbytes)": "max_rss_kb",
    "Percent of CPU this job got": "cpu_percent",
    "Exit status": "exit_status",
}


def parse_time_v(text: str) -> dict:
    out = {}
    for line in text.splitlines():
        line = line.strip()
        for k, v in _TIME_KEYS.items():
            if line.startswith(k + ":"):
                out[v] = line[len(k) + 1:].strip()
    for k in ("user_s", "sys_s"):
        if k in out:
            out[k] = float(out[k])
    if "max_rss_kb" in out:
        out["max_rss_kb"] = int(out["max_rss_kb"])
    if "exit_status" in out:
        out["exit_status"] = int(out["exit_status"])
    return out


def run_timed(cmd: list[str], log: Path, cwd: Optional[Path] = None, env=None) -> dict:
    """Run ``cmd`` under /usr/bin/time -v, stream stdout+stderr to ``log``.
    Note: GNU time's max RSS is the largest single (waited-for) descendant, not the sum."""
    tfile = log.with_suffix(".time")
    full = ["/usr/bin/time", "-v", "-o", str(tfile)] + [str(c) for c in cmd]
    t0 = time.time()
    with open(log, "w") as fh:
        fh.write("# " + " ".join(map(str, cmd)) + "\n")
        fh.flush()
        rc = subprocess.call(full, stdout=fh, stderr=subprocess.STDOUT, cwd=cwd, env=env or tool_env())
    res = {"cmd": [str(c) for c in cmd], "returncode": rc, "wall_s": round(time.time() - t0, 1),
           "log": str(log)}
    if tfile.exists():
        res.update(parse_time_v(tfile.read_text()))
    return res


def capture(cmd: list[str]) -> str:
    try:
        p = subprocess.run(cmd, capture_output=True, text=True, env=tool_env(), timeout=120)
        return re.sub(r"\x1b\[[0-9;]*m|\x1b\]8;[^\x1b]*\x1b\\", "", (p.stdout + p.stderr)).strip()
    except (OSError, subprocess.TimeoutExpired) as e:  # pragma: no cover
        return f"ERROR: {e}"


def git_commit() -> str:
    try:
        return subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT, capture_output=True,
                              text=True).stdout.strip()
    except OSError:  # pragma: no cover
        return ""


def write_json(path: Path, obj) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(obj, indent=1, sort_keys=False, default=str) + "\n")
    tmp.replace(path)


def log(msg: str) -> None:
    print(f"[{time.strftime('%Y-%m-%dT%H:%M:%S')}] {msg}", file=sys.stderr, flush=True)


def write_tree(root: Node, path: Path) -> None:
    path.write_text(write_newick(root, lengths=False) + "\n")


__all__ = [n for n in dir() if not n.startswith("_")] + ["read_newick"]
