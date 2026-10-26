"""Run a HOG-inference tool on a HOBench panel and archive raw outputs + a run manifest.

    python -m wrapper.run --tool orthofinder|fastoma|of_flat --panel Y --mode full|subset \
        --variant default [--species-tree path] [--threads 6]

Work directories live on local NVMe (``$HOBENCH_WORK``, default /tmp/hobench_work); only final
outputs are copied to ``runs/<panel>/<mode>/<tool>_<variant>/``. The wrapper never reads
``families/`` (blind to the references).
"""
from __future__ import annotations

import argparse
import gzip
import io
import json
import os
import platform
import shutil
import socket
import sys
import tarfile
import time
from pathlib import Path
from typing import Optional

from . import common as C
from hobench.newick import read_newick

# Variants (SYNTHESIS E7). Each maps to extra command-line arguments.
VARIANTS = {
    "orthofinder": {"default": [], "dendroblast": ["-M", "dendroblast"]},
    "fastoma": {"default": []},
    "of_flat": {"default": [], "dendroblast": []},  # derived from orthofinder_<variant>
}
# FastOMA 0.5.0 cannot run without a species tree (fastoma-check-input: --species-tree
# required=True), so the "FastOMA without our tree" variant of E7 is not available. A
# different tree can be given with --species-tree (e.g. OrthoFinder's STAG tree) + --variant.

BIG_BYTES = 20 * 1024 * 1024  # files above this go to raw/ (gitignored)
DEFAULT_OMAMER_DB = C.ROOT / "data" / "raw" / "omamer" / "LUCA.h5"
DEFAULT_HOGMAP_DIR = C.ROOT / "data" / "hogmaps"


# ----------------------------------------------------------------------------- helpers
def run_dir(panel: str, mode: str, tool: str, variant: str, runs_root: Optional[Path] = None) -> Path:
    return (runs_root or C.ROOT / "runs") / panel / mode / f"{tool}_{variant}"


def stage_inputs(panel: C.Panel, dest: Path) -> dict:
    """Copy proteomes to dest/<CODE>.fa; return per-species input records."""
    dest.mkdir(parents=True, exist_ok=True)
    recs = {}
    for code, src in panel.fastas.items():
        tgt = dest / f"{code}.fa"
        shutil.copyfile(src.resolve(), tgt)
        n = sum(1 for _ in C.read_fasta(tgt))
        recs[code] = {"source": str(src), "sha256": C.sha256_file(tgt), "n_seqs": n}
    return recs


def check_tree(panel: C.Panel) -> tuple:
    root = read_newick(panel.tree_path)
    probs = C.validate_species_tree(root, panel.codes)
    if probs:
        raise C.WrapperError(f"species tree {panel.tree_path} invalid:\n  " + "\n  ".join(probs))
    return root


def place_output(src: Path, out_dir: Path, name: Optional[str] = None) -> dict:
    """Copy a file into out_dir (or out_dir/raw if big); return its record."""
    name = name or src.name
    size = src.stat().st_size
    tgt = (out_dir / "raw" / name) if size > BIG_BYTES else out_dir / name
    tgt.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(src, tgt)
    return {"path": str(tgt.relative_to(out_dir)), "bytes": size, "sha256": C.sha256_file(tgt)}


def gzip_to(src: Path, tmp_dir: Path) -> Path:
    out = tmp_dir / (src.name + ".gz")
    with open(src, "rb") as fi, gzip.open(out, "wb", compresslevel=6) as fo:
        shutil.copyfileobj(fi, fo)
    return out


def tar_dirs(base: Path, members: list[str], out: Path) -> Path:
    with tarfile.open(out, "w:gz") as tf:
        for m in members:
            if (base / m).exists():
                tf.add(base / m, arcname=m)
    return out


def tool_versions(tool: str) -> dict:
    v = {"python": platform.python_version()}
    if tool in ("orthofinder", "of_flat"):
        v["orthofinder"] = C.capture(["orthofinder", "-v"]).splitlines()[-1:] or ""
        v["diamond"] = C.capture(["diamond", "version"])
    if tool == "fastoma":
        v["nextflow"] = C.capture(["nextflow", "-v"])
        v["fastoma_python"] = C.capture([sys.executable, "-c",
                                         "import importlib.metadata as m;print(m.version('FastOMA'))"])
        v["fastoma_workflow_git"] = C.capture(["git", "-C", str(C.ROOT / "tools" / "FastOMA"),
                                               "describe", "--tags", "--always"])
        v["omamer"] = C.capture(["omamer", "--version"]).splitlines()[-1:]
    return v


# ----------------------------------------------------------------------------- hogmaps
def assemble_hogmaps(fasta_dir: Path, codes: list[str], hogmap_dir: Path, out_dir: Path,
                     omamer_db: Path, work: Path, threads: int = 4) -> dict:
    """Build out_dir/<CODE>.fa.hogmap restricted to the IDs in the input FASTA. IDs absent
    from the precomputed hogmap (subset/perturbed panels: fragments, contaminants) are
    searched with omamer on just those sequences and appended. If no precomputed hogmap
    exists for a species, nothing is written and FastOMA runs omamer itself (maxForks=2)."""
    out_dir.mkdir(parents=True, exist_ok=True)
    recs = {}
    for code in codes:
        fa = fasta_dir / f"{code}.fa"
        src = hogmap_dir / f"{code}.fa.hogmap"
        if not src.exists():
            recs[code] = {"source": None, "note": "no precomputed hogmap; FastOMA computes it"}
            continue
        records = {C.fasta_id(h): (h, s) for h, s in C.read_fasta(fa)}
        ids = set(records)
        head, body, seen = [], [], set()
        with open(src) as fh:
            for line in fh:
                if line.startswith("!") or line.startswith("qseqid\t"):
                    head.append(line)
                    continue
                q = line.split("\t", 1)[0]
                if q in ids:
                    body.append(line)
                    seen.add(q)
        missing = sorted(ids - seen)
        extra = []
        if missing:
            mfa = work / f"{code}.missing.fa"
            C.write_fasta(mfa, (records[i] for i in missing))
            mout = work / f"{code}.missing.hogmap"
            res = C.run_timed(["omamer", "search", "-n", "10", "-t", str(threads), "--db", str(omamer_db),
                               "--query", str(mfa), "--out", str(mout)], work / f"{code}.omamer.log")
            if res["returncode"] != 0:
                raise C.WrapperError(f"omamer on missing {code} sequences failed: {res}")
            extra = [ln for ln in open(mout) if not ln.startswith("!") and not ln.startswith("qseqid\t")]
        tgt = out_dir / f"{code}.fa.hogmap"
        with open(tgt, "w") as fo:
            fo.writelines(head)
            fo.writelines(body)
            fo.writelines(extra)
        recs[code] = {"source": str(src), "source_sha256": C.sha256_file(src), "sha256": C.sha256_file(tgt),
                      "n_ids": len(ids), "n_from_precomputed": len(seen),
                      "n_searched_now": len(missing)}
    return recs


# ----------------------------------------------------------------------------- tools
def run_orthofinder(panel, args, work: Path, out_dir: Path, man: dict) -> None:
    root = check_tree(panel)
    tree = work / "species_tree.nwk"
    shutil.copyfile(panel.tree_path, tree)
    man["inputs"] = stage_inputs(panel, work / "in")
    of_out = work / "of_out"
    if of_out.exists():
        shutil.rmtree(of_out)
    cmd = ["orthofinder", "-f", work / "in", "-s", tree, "-t", args.threads, "-a", args.threads,
           "-o", of_out, "-n", args.variant] + VARIANTS["orthofinder"][args.variant] + args.extra
    man["commands"].append([str(c) for c in cmd])
    C.write_json(out_dir / "manifest.json", man)
    if args.dry_run:
        return
    C.log(f"orthofinder start: {' '.join(map(str, cmd))}")
    res = C.run_timed(cmd, work / "orthofinder.log")
    man["timing"] = res
    if res["returncode"] != 0:
        raise C.WrapperError(f"orthofinder failed rc={res['returncode']} (log {work / 'orthofinder.log'})")
    rdirs = sorted(of_out.glob("Results_*"))
    if len(rdirs) != 1:
        raise C.WrapperError(f"expected one Results_* dir in {of_out}, got {rdirs}")
    rd = rdirs[0]
    # species-tree check: OF relabels internal nodes N0..; map them to our levels by clade
    of_tree = read_newick(rd / "Species_Tree" / "SpeciesTree_rooted_node_labels.txt")
    ours = {v: k for k, v in C.tree_clades(root).items()}
    level_map, unmatched = {}, []
    for n_name, clade in C.tree_clades(of_tree).items():
        if clade in ours:
            level_map[n_name] = ours[clade]
        else:
            unmatched.append(n_name)
    man["checks"]["of_species_tree_matches_input"] = (not unmatched and len(level_map) == len(ours))
    man["checks"]["of_unmatched_nodes"] = unmatched
    lm = work / "of_level_map.tsv"
    lm.write_text("of_node\tlevel\n" + "".join(f"{k}\t{v}\n" for k, v in sorted(level_map.items())))
    outs = man["outputs"]
    outs["of_level_map.tsv"] = place_output(lm, out_dir)
    outs["SpeciesTree_rooted_node_labels.txt"] = place_output(
        rd / "Species_Tree" / "SpeciesTree_rooted_node_labels.txt", out_dir)
    core = tar_dirs(rd, ["Phylogenetic_Hierarchical_Orthogroups", "Orthogroups", "Species_Tree",
                         "Gene_Duplication_Events", "Comparative_Genomics_Statistics", "Log.txt"],
                    work / "of_core.tar.gz")
    outs["of_core.tar.gz"] = place_output(core, out_dir)
    orth = tar_dirs(rd, ["Orthologues"], work / "of_orthologues.tar.gz")
    outs["of_orthologues.tar.gz"] = place_output(orth, out_dir)
    outs["orthofinder.log"] = place_output(work / "orthofinder.log", out_dir)


def run_fastoma(panel, args, work: Path, out_dir: Path, man: dict) -> None:
    root = check_tree(panel)
    inp = work / "in"
    if inp.exists():
        shutil.rmtree(inp)
    man["inputs"] = stage_inputs(panel, inp / "proteome")
    shutil.copyfile(panel.tree_path, inp / "species_tree.nwk")
    omamer_db = Path(args.omamer_db).resolve()
    man["omamer_db"] = {"path": str(omamer_db), "sha256": C.sha256_cached(omamer_db)}
    man["hogmaps"] = assemble_hogmaps(inp / "proteome", panel.codes, Path(args.hogmap_dir), inp / "hogmap_in",
                                      omamer_db, work, threads=min(4, int(args.threads)))
    cfg = work / "fastoma_run.config"
    cfg.write_text(
        f"includeConfig '{C.ROOT / 'config' / 'fastoma_local.config'}'\n"
        f"executor {{ name = 'local'; cpus = {int(args.threads)}; memory = '{int(args.mem_gb)} GB' }}\n"
        "process {\n"
        f"    resourceLimits = [cpus: {int(args.threads)}, memory: {int(args.mem_gb)}.GB]\n"
        f"    withName: \"hog_big\"  {{ cpus = {min(4, int(args.threads))}; memory = '12 GB' }}\n"
        "    withName: \"omamer_run\" { maxForks = 2; cpus = 4; memory = '16 GB' }\n"
        "}\n")
    out = work / "fo_out"
    if out.exists():
        shutil.rmtree(out)
    cmd = ["nextflow", "run", C.ROOT / "tools" / "FastOMA" / "FastOMA.nf", "-c", cfg,
           "--input", inp, "--output_folder", out, "--omamer_db", omamer_db,
           "-work-dir", work / "nf_work", "-with-trace", work / "trace.txt"] + \
        VARIANTS["fastoma"][args.variant] + args.extra
    man["commands"].append([str(c) for c in cmd])
    man["fastoma_config"] = cfg.read_text()
    C.write_json(out_dir / "manifest.json", man)
    if args.dry_run:
        return
    C.log(f"fastoma start: {' '.join(map(str, cmd))}")
    res = C.run_timed(cmd, work / "fastoma.log", cwd=work)
    man["timing"] = res
    if res["returncode"] != 0:
        raise C.WrapperError(f"FastOMA failed rc={res['returncode']} (log {work / 'fastoma.log'})")
    checked = read_newick(out / "species_tree_checked.nwk")
    ours, theirs = C.tree_clades(root), C.tree_clades(checked)
    renamed = sorted(k for k in ours if theirs.get(k) != ours[k])
    man["checks"]["fastoma_levels_unchanged"] = not renamed
    man["checks"]["fastoma_renamed_or_missing_levels"] = renamed
    outs = man["outputs"]
    tmp = work / "gz"
    tmp.mkdir(exist_ok=True)
    outs["FastOMA_HOGs.orthoxml.gz"] = place_output(gzip_to(out / "FastOMA_HOGs.orthoxml", tmp), out_dir)
    for f in ("RootHOGs.tsv", "OrthologousGroups.tsv"):
        if (out / f).exists():
            outs[f + ".gz"] = place_output(gzip_to(out / f, tmp), out_dir)
    for f in ("orthologs.tsv.gz", "species_tree_checked.nwk"):
        if (out / f).exists():
            outs[f] = place_output(out / f, out_dir)
    for f in ("trace.txt", "fastoma.log"):
        if (work / f).exists():
            outs[f] = place_output(work / f, out_dir)


def run_of_flat(panel, args, work: Path, out_dir: Path, man: dict) -> None:
    """Flat baseline: OrthoFinder orthogroups (+ unassigned genes as singletons) as root-level
    HOGs, no hierarchy. Derived from an existing orthofinder_<variant> run, no new inference."""
    root = check_tree(panel)
    src_dir = run_dir(panel.name, panel.mode, "orthofinder", args.variant, args.runs_root)
    src_man = src_dir / "manifest.json"
    if not src_man.exists():
        raise C.WrapperError(f"{src_man} missing: run orthofinder --variant {args.variant} first")
    sm = json.loads(src_man.read_text())
    tar_rec = sm["outputs"]["of_core.tar.gz"]
    tar_path = src_dir / tar_rec["path"]
    man["derived_from"] = {"run": str(src_dir), "manifest_sha256": C.sha256_file(src_man),
                           "of_core_sha256": C.sha256_file(tar_path)}
    man["inputs"] = sm.get("inputs", {})
    man["species_tree"] = sm.get("species_tree", man.get("species_tree"))
    rows = []
    with tarfile.open(tar_path) as tf:
        def table(name):
            m = tf.getmember(name)
            return io.TextIOWrapper(tf.extractfile(m), encoding="utf-8").read().splitlines()
        og = table("Orthogroups/Orthogroups.tsv")
        sp = og[0].split("\t")[1:]
        for line in og[1:]:
            parts = line.split("\t")
            for s, cell in zip(sp, parts[1:]):
                for g in filter(None, (x.strip() for x in cell.split(","))):
                    rows.append((parts[0], g, s))
        try:
            un = table("Orthogroups/Orthogroups_UnassignedGenes.tsv")
        except KeyError:
            un = []
        if un:
            sp = un[0].split("\t")[1:]
            for line in un[1:]:
                parts = line.split("\t")
                for s, cell in zip(sp, parts[1:]):
                    for g in filter(None, (x.strip() for x in cell.split(","))):
                        rows.append((parts[0], g, s))
    lvl = root.name
    out = work / "hogs.tsv.gz"
    with gzip.open(out, "wt") as fh:
        fh.write("level\thog_id\tgene\tspecies\n")
        for og_id, g, s in rows:
            fh.write(f"{lvl}\t{og_id}\t{g}\t{s}\n")
    man["checks"]["n_genes"] = len(rows)
    man["checks"]["n_groups"] = len({r[0] for r in rows})
    man["outputs"]["hogs.tsv.gz"] = place_output(out, out_dir)
    man["timing"] = {"wall_s": 0, "note": "derived, no inference"}


# ----------------------------------------------------------------------------- main
def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="python -m wrapper.run", description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tool", required=True, choices=sorted(VARIANTS))
    ap.add_argument("--panel", required=True)
    ap.add_argument("--mode", default="full", choices=["full", "subset"])
    ap.add_argument("--variant", default="default")
    ap.add_argument("--species-tree", type=Path, default=None,
                    help="override panels/<P>/species_tree.nwk (use a distinct --variant name)")
    ap.add_argument("--threads", type=int, default=6)
    ap.add_argument("--mem-gb", type=int, default=30, help="FastOMA executor memory cap")
    ap.add_argument("--omamer-db", default=str(DEFAULT_OMAMER_DB))
    ap.add_argument("--hogmap-dir", default=str(DEFAULT_HOGMAP_DIR))
    ap.add_argument("--panels-root", type=Path, default=None)
    ap.add_argument("--runs-root", type=Path, default=None)
    ap.add_argument("--work-root", type=Path, default=C.WORK_ROOT)
    ap.add_argument("--keep-work", action="store_true", help="do not delete the NVMe work dir on success")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--extra", nargs=argparse.REMAINDER, default=[],
                    help="extra tool args (recorded in the manifest)")
    args = ap.parse_args(argv)

    if args.variant not in VARIANTS[args.tool] and args.species_tree is None and not args.extra:
        ap.error(f"unknown variant {args.variant!r} for {args.tool}; known: {sorted(VARIANTS[args.tool])}")
    VARIANTS[args.tool].setdefault(args.variant, [])
    panel = C.load_panel(args.panel, args.mode, args.panels_root, args.species_tree)
    out_dir = run_dir(args.panel, args.mode, args.tool, args.variant, args.runs_root)
    work = args.work_root / args.panel / args.mode / f"{args.tool}_{args.variant}"
    work.mkdir(parents=True, exist_ok=True)
    out_dir.mkdir(parents=True, exist_ok=True)
    man = {
        "tool": args.tool, "variant": args.variant, "panel": args.panel, "mode": args.mode,
        "status": "running", "started": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
        "host": socket.gethostname(), "cpus_available": os.cpu_count(), "threads": args.threads,
        "git_commit": C.git_commit(), "work_dir": str(work), "argv": sys.argv,
        "species_tree": {"path": str(panel.tree_path), "sha256": C.sha256_file(panel.tree_path),
                         "newick": panel.tree_path.read_text().strip()},
        "species": panel.species, "commands": [], "checks": {}, "outputs": {},
        "tool_versions": tool_versions(args.tool),
    }
    fn = {"orthofinder": run_orthofinder, "fastoma": run_fastoma, "of_flat": run_of_flat}[args.tool]
    try:
        fn(panel, args, work, out_dir, man)
        man["status"] = "dry_run" if args.dry_run else "done"
    except Exception as e:
        man["status"] = "failed"
        man["error"] = f"{type(e).__name__}: {e}"
        C.write_json(out_dir / "manifest.json", man)
        C.log(f"FAILED {args.tool} {args.panel}/{args.mode}/{args.variant}: {e}")
        return 1
    man["finished"] = time.strftime("%Y-%m-%dT%H:%M:%S%z")
    C.write_json(out_dir / "manifest.json", man)
    C.log(f"{man['status']} {args.tool} {args.panel}/{args.mode}/{args.variant} -> {out_dir}")
    if man["status"] == "done" and not args.keep_work and args.tool != "of_flat":
        shutil.rmtree(work, ignore_errors=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
