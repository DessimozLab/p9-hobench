"""Robustness perturbation generators (DESIGN §7 E2; SYNTHESIS E5; Y panel, full mode).

    python -m wrapper.perturb fragment    --panel Y --genes genes.txt [--fraction 0.1]
    python -m wrapper.perturb contaminate --panel Y --target CODE --source-fasta F --n 20 [--source-label HUMAN]
    python -m wrapper.perturb alttree     --panel Y [--tree panels/Y/species_tree_alt.nwk]
    python -m wrapper.perturb remove      --panel Y --taxa CODE1,CODE2

Each writes a self-contained panel ``panels/<P>_<name>/`` (species.tsv, proteomes/<CODE>.fa,
species_tree.nwk) and ``perturbation.json`` listing every change, which scoring needs:
fragments (orig id -> two new ids + coordinates), inserted contaminants (new id -> source id),
tree clades lost/gained, removed species/proteins/levels. Unchanged proteomes are relative
symlinks to the parent panel; changed ones are new files. Seeded with 20261026 by default.
"""
from __future__ import annotations

import argparse
import math
import os
import shutil
import sys
import time
from pathlib import Path

from . import common as C
from hobench.newick import read_newick

MIN_FRAG = 30  # aa; FastOMA drops sequences < 20 aa


def _new_panel(panel: C.Panel, name: str, panels_root) -> Path:
    out = C.panel_dir(f"{panel.name}_{name}", panels_root)
    if out.exists():
        shutil.rmtree(out)
    (out / "proteomes").mkdir(parents=True)
    return out


def _link(src: Path, dst: Path) -> None:
    os.symlink(os.path.relpath(src.resolve(), dst.parent), dst)


def _write_species(out: Path, rows: list[dict]) -> None:
    cols = list(rows[0].keys())
    with open(out / "species.tsv", "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(r.get(c, "") for c in cols) + "\n")


def _base_manifest(panel: C.Panel, kind: str, seed: int, args) -> dict:
    return {"perturbation": kind, "parent_panel": panel.name, "seed": seed,
            "created": time.strftime("%Y-%m-%dT%H:%M:%S%z"), "git_commit": C.git_commit(),
            "argv": sys.argv, "parent_tree_sha256": C.sha256_file(panel.tree_path),
            "parent_proteomes": {c: C.sha256_file(p) for c, p in panel.fastas.items()}}


def _finish(out: Path, man: dict) -> None:
    man["proteomes"] = {p.stem: C.sha256_file(p) for p in sorted((out / "proteomes").glob("*.fa"))}
    man["species_tree_sha256"] = C.sha256_file(out / "species_tree.nwk")
    C.write_json(out / "perturbation.json", man)
    C.log(f"perturbed panel written: {out}")


def _read_ids(path: Path) -> list[str]:
    return [ln.split()[0] for ln in path.read_text().splitlines() if ln.strip() and not ln.startswith("#")]


# ----------------------------------------------------------------------------- (a) fragmentation
def fragment(panel: C.Panel, genes: list[str], fraction: float, seed: int, out: Path, man: dict) -> None:
    wanted = set(genes)
    eligible, not_found, too_short = [], set(wanted), []
    recs = {}
    for code, p in panel.fastas.items():
        recs[code] = list(C.read_fasta(p))
        for h, s in recs[code]:
            rid = C.fasta_id(h)
            key = rid if rid in wanted else C.accession(rid) if C.accession(rid) in wanted else None
            if key is None:
                continue
            not_found.discard(key)
            (eligible if len(s) >= 2 * MIN_FRAG else too_short).append((code, rid))
    k = math.ceil(fraction * len(eligible)) if eligible else 0
    chosen = set(sorted(eligible, key=lambda x: C.stable_unit(seed, "fragment", x[1]))[:k])
    chosen_ids = {rid for _, rid in chosen}
    changes = []
    for code, p in panel.fastas.items():
        dst = out / "proteomes" / f"{code}.fa"
        if not any(c == code for c, _ in chosen):
            _link(p, dst)
            continue
        new = []
        for h, s in recs[code]:
            rid = C.fasta_id(h)
            if rid not in chosen_ids:
                new.append((h, s))
                continue
            L = len(s)
            u = C.stable_unit(seed, "breakpoint", rid)
            b = int(round(L * (0.3 + 0.4 * u)))
            b = min(max(b, MIN_FRAG), L - MIN_FRAG)
            acc = C.accession(rid)
            f1, f2 = f"{acc}_f1", f"{acc}_f2"
            new.append((f"{f1} fragment=1/2 of={rid} range=1-{b}", s[:b]))
            new.append((f"{f2} fragment=2/2 of={rid} range={b + 1}-{L}", s[b:]))
            changes.append({"species": code, "orig_id": rid, "orig_len": L, "breakpoint": b,
                            "fragments": [{"id": f1, "start": 1, "end": b}, {"id": f2, "start": b + 1, "end": L}]})
        C.write_fasta(dst, new)
    man.update(fraction=fraction, n_genes_listed=len(wanted), n_eligible=len(eligible),
               genes_not_found=sorted(not_found), too_short=[r for _, r in too_short],
               n_fragmented=len(changes), changes=changes)


# ----------------------------------------------------------------------------- (c) contamination
def contaminate(panel: C.Panel, target: str, src_fasta: Path, n: int, label: str, seed: int,
                ids: list[str] | None, out: Path, man: dict) -> None:
    if target not in panel.fastas:
        raise C.WrapperError(f"target {target} not in panel {panel.name}")
    src = {C.fasta_id(h): (h, s) for h, s in C.read_fasta(src_fasta)}
    if ids:
        pick = [i for i in ids if i in src or any(C.accession(k) == i for k in src)]
        pick = [i if i in src else next(k for k in src if C.accession(k) == i) for i in pick]
    else:
        pick = sorted((k for k, (h, s) in src.items() if len(s) >= 50),
                      key=lambda r: C.stable_unit(seed, "contaminant", r))[:n]
    changes = []
    for code, p in panel.fastas.items():
        dst = out / "proteomes" / f"{code}.fa"
        if code != target:
            _link(p, dst)
            continue
        new = list(C.read_fasta(p))
        for rid in pick:
            nid = f"CONTAM_{label}_{C.accession(rid)}"
            new.append((f"{nid} contaminant_from={label} source_id={rid}", src[rid][1]))
            changes.append({"species": target, "new_id": nid, "source_id": rid, "source_header": src[rid][0],
                            "len": len(src[rid][1])})
        C.write_fasta(dst, new)
    man.update(target=target, source_fasta=str(src_fasta), source_sha256=C.sha256_file(src_fasta),
               source_label=label, n_requested=n, n_inserted=len(changes), changes=changes)


# ----------------------------------------------------------------------------- (d) alt tree / (e) removal
def tree_diff(a, b, restrict: set[str] | None = None) -> dict:
    """Compare clades by species set (parent clades restricted to ``restrict`` if given; a
    restricted clade claimed by several parent nodes is attributed to the youngest one)."""
    ca, cb = C.tree_clades(a), C.tree_clades(b)
    if restrict is not None:
        ca = {k: v & restrict for k, v in ca.items()}
    sa = {}
    for k, v in sorted(ca.items(), key=lambda kv: -len(C.tree_clades(a)[kv[0]])):
        if len(v) >= 2:
            sa[v] = k  # younger (smaller original clade) overwrites older
    sb = {v: k for k, v in cb.items()}
    return {"shared": sorted([sa[c], sb[c]] for c in set(sa) & set(sb)),
            "only_parent": sorted([sa[c], sorted(c)] for c in set(sa) - set(sb)),
            "only_new": sorted([sb[c], sorted(c)] for c in set(sb) - set(sa))}


def alttree(panel: C.Panel, tree: Path, out: Path, man: dict) -> None:
    alt = read_newick(tree)
    probs = C.validate_species_tree(alt, panel.codes)
    if probs:
        raise C.WrapperError(f"alternative tree invalid: {probs}")
    for code, p in panel.fastas.items():
        _link(p, out / "proteomes" / f"{code}.fa")
    shutil.copyfile(tree, out / "species_tree.nwk")
    man.update(alt_tree=str(tree), alt_tree_sha256=C.sha256_file(tree),
               clades=tree_diff(read_newick(panel.tree_path), alt))


def remove(panel: C.Panel, taxa: list[str], out: Path, man: dict) -> list[dict]:
    bad = [t for t in taxa if t not in panel.fastas]
    if bad:
        raise C.WrapperError(f"unknown taxa {bad}")
    keep = [c for c in panel.codes if c not in taxa]
    base = read_newick(panel.tree_path)
    pruned = C.prune_tree(base, set(keep))
    C.write_tree(pruned, out / "species_tree.nwk")
    for code in keep:
        _link(panel.fastas[code], out / "proteomes" / f"{code}.fa")
    removed_prot = {t: [C.fasta_id(h) for h, _ in C.read_fasta(panel.fastas[t])] for t in taxa}
    man.update(removed_species=taxa, removed_protein_counts={t: len(v) for t, v in removed_prot.items()},
               clades=tree_diff(base, pruned, set(keep)),
               removed_levels=sorted(set(C.tree_clades(base)) - set(C.tree_clades(pruned))))
    return [r for r in panel.species if r["code"] in keep]


# ----------------------------------------------------------------------------- main
def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="python -m wrapper.perturb", description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="kind", required=True)
    for k in ("fragment", "contaminate", "alttree", "remove"):
        s = sub.add_parser(k)
        s.add_argument("--panel", required=True)
        s.add_argument("--name", default=None, help="suffix of the new panel (default: frag/contam/alttree/taxrm)")
        s.add_argument("--seed", type=int, default=C.GLOBAL_SEED)
        s.add_argument("--panels-root", type=Path, default=None)
        if k == "fragment":
            s.add_argument("--genes", type=Path, required=True, help="gene IDs/accessions, one per line")
            s.add_argument("--fraction", type=float, default=0.10)
        if k == "contaminate":
            s.add_argument("--target", required=True)
            s.add_argument("--source-fasta", type=Path, required=True)
            s.add_argument("--source-label", default="SRC")
            s.add_argument("--n", type=int, default=20)
            s.add_argument("--ids", type=Path, default=None, help="explicit source protein IDs")
        if k == "alttree":
            s.add_argument("--tree", type=Path, default=None)
        if k == "remove":
            s.add_argument("--taxa", required=True, help="comma-separated species codes")
    args = ap.parse_args(argv)
    panel = C.load_panel(args.panel, "full", args.panels_root)
    name = args.name or {"fragment": "frag", "contaminate": "contam", "alttree": "alttree", "remove": "taxrm"}[args.kind]
    out = _new_panel(panel, name, args.panels_root)
    man = _base_manifest(panel, args.kind, args.seed, args)
    rows = panel.species
    if args.kind == "fragment":
        fragment(panel, _read_ids(args.genes), args.fraction, args.seed, out, man)
        man["genes_file_sha256"] = C.sha256_file(args.genes)
    elif args.kind == "contaminate":
        contaminate(panel, args.target, args.source_fasta, args.n, args.source_label, args.seed,
                    _read_ids(args.ids) if args.ids else None, out, man)
    elif args.kind == "alttree":
        alttree(panel, args.tree or panel.dir / "species_tree_alt.nwk", out, man)
    else:
        rows = remove(panel, [t.strip() for t in args.taxa.split(",") if t.strip()], out, man)
    if not (out / "species_tree.nwk").exists():
        shutil.copyfile(panel.tree_path, out / "species_tree.nwk")
    _write_species(out, rows)
    _finish(out, man)
    return 0


if __name__ == "__main__":
    sys.exit(main())
