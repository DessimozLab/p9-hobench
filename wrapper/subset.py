"""Family-neighbourhood subset builder (DESIGN §7 data mode; SYNTHESIS D5).

    python -m wrapper.subset --panel V --seeds seeds.fa|accessions.txt [--threads 8]

Seeds are *raw literature member sequences* (FASTA) or an accession list (accessions are looked
up in the full panel proteomes). Nothing is read from ``families/``.

Subset = union of
  1. neighbourhood: panel proteins with a DIAMOND hit (blastp --more-sensitive, e <= 1e-5,
     fixed --dbsize) to any seed with coverage >= --min-cov on the query OR subject side;
  2. domain-only decoys: proteins hit only with lower coverage, a seeded capped sample
     (promiscuous-domain decoys);
  3. background: a fixed seeded fraction of every proteome, independent of the seeds
     (selection key = sha256(seed|protein id), so smaller fractions give nested samples).
Total is capped (default 70k); if over the cap the background fraction is reduced, and if the
neighbourhood alone exceeds the cap the build stops (sizes are reported; use --sizes-only).
Writes panels/<P>/subset/<CODE>.fa, members.tsv (id, species, category) and manifest.json.
"""
from __future__ import annotations

import argparse
import gzip
import math
import shutil
import sys
import time
from collections import defaultdict
from pathlib import Path

from . import common as C

DEFAULT_DBSIZE = 1_000_000_000  # letters; fixed so e-values are comparable across panels/runs


def load_seeds(path: Path, full: dict[str, dict[str, tuple[str, str]]]) -> tuple[list[tuple[str, str]], dict]:
    """Return seed records (id, seq) and a report. FASTA → as is; otherwise accession list."""
    txt = path.read_text()
    report = {"source": str(path), "sha256": C.sha256_file(path)}
    if txt.lstrip().startswith(">"):
        recs = [(f"seed{i:05d}|{C.fasta_id(h)}", s) for i, (h, s) in enumerate(C.read_fasta(path))]
        report.update(kind="fasta", n=len(recs))
        return recs, report
    accs = [ln.split()[0] for ln in txt.splitlines() if ln.strip() and not ln.startswith("#")]
    index = {}
    for code, recs in full.items():
        for rid, (h, s) in recs.items():
            index[C.accession(rid)] = (rid, s)
            index[rid] = (rid, s)
    out, missing = [], []
    for i, a in enumerate(accs):
        if a in index:
            out.append((f"seed{i:05d}|{index[a][0]}", index[a][1]))
        else:
            missing.append(a)
    report.update(kind="accessions", n=len(accs), n_found=len(out), missing=missing)
    return out, report


def run_diamond(seeds: list[tuple[str, str]], full, work: Path, threads: int, dbsize: int,
                evalue: float) -> tuple[Path, list[str]]:
    work.mkdir(parents=True, exist_ok=True)
    db_fa = work / "panel.fa"
    C.write_fasta(db_fa, ((rid, s) for code in sorted(full) for rid, (h, s) in full[code].items()))
    q = work / "seeds.fa"
    C.write_fasta(q, seeds)
    env = C.tool_env()
    cmds = [["diamond", "makedb", "--in", str(db_fa), "-d", str(work / "panel"), "--threads", str(threads)],
            ["diamond", "blastp", "-q", str(q), "-d", str(work / "panel"), "-o", str(work / "hits.tsv"),
             "--more-sensitive", "--evalue", str(evalue), "--dbsize", str(dbsize), "-k", "0",
             "--threads", str(threads), "--outfmt", "6", "qseqid", "sseqid", "pident", "length", "qlen",
             "slen", "evalue", "bitscore", "qcovhsp", "scovhsp"]]
    for cmd in cmds:
        res = C.run_timed(cmd, work / f"{cmd[1]}.log", env=env)
        if res["returncode"] != 0:
            raise C.WrapperError(f"diamond failed: {res}")
    return work / "hits.tsv", [" ".join(c) for c in cmds]


def classify_hits(hits: Path, min_cov: float) -> tuple[set[str], set[str], dict]:
    """neighbourhood = subjects with any HSP covering >= min_cov of query or subject;
    domain-only = subjects with hits but none reaching min_cov."""
    best = defaultdict(float)
    with open(hits) as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            s, qcov, scov = p[1], float(p[8]) / 100.0, float(p[9]) / 100.0
            best[s] = max(best[s], qcov, scov)
    nb = {s for s, c in best.items() if c >= min_cov}
    dom = set(best) - nb
    return nb, dom, dict(best)


def select(full, nb: set[str], dom: set[str], seed: int, bg_fraction: float, dom_cap: int, cap: int):
    sp_of = {rid: code for code, recs in full.items() for rid in recs}
    dom_sorted = sorted(dom, key=lambda r: C.stable_unit(seed, "domain", r))
    dom_sel = set(dom_sorted[:dom_cap])
    core = nb | dom_sel
    if len(core) > cap:
        raise C.WrapperError(f"neighbourhood ({len(nb)}) + domain decoys ({len(dom_sel)}) exceed cap {cap}")
    keys = {rid: C.stable_unit(seed, "background", rid) for rid in sp_of}
    frac = bg_fraction
    bg = {r for r, k in keys.items() if k < frac and r not in core}
    if len(core) + len(bg) > cap:
        # reduce fraction to fit: take the (cap - |core|) non-core proteins with the smallest keys
        room = cap - len(core)
        ordered = sorted((k, r) for r, k in keys.items() if r not in core and k < frac)
        frac = ordered[room][0] if room < len(ordered) else frac
        bg = {r for k, r in ordered[:room]}
    cat = {}
    for r in bg:
        cat[r] = "background"
    for r in dom_sel:
        cat[r] = "domain_only"
    for r in nb:
        cat[r] = "neighbourhood"
    return cat, frac, sp_of


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="python -m wrapper.subset", description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--panel", required=True)
    ap.add_argument("--seeds", required=True, type=Path, help="FASTA of raw member sequences or accession list")
    ap.add_argument("--seed", type=int, default=C.GLOBAL_SEED)
    ap.add_argument("--evalue", type=float, default=1e-5)
    ap.add_argument("--dbsize", type=int, default=DEFAULT_DBSIZE)
    ap.add_argument("--min-cov", type=float, default=0.5, help="coverage threshold for a full (non-domain) hit")
    ap.add_argument("--domain-cap", type=int, default=2000, help="max domain-only decoys")
    ap.add_argument("--bg-fraction", type=float, default=0.05, help="background fraction per proteome")
    ap.add_argument("--cap", type=int, default=70000, help="max proteins in the subset panel")
    ap.add_argument("--threads", type=int, default=6)
    ap.add_argument("--panels-root", type=Path, default=None)
    ap.add_argument("--out", type=Path, default=None, help="default panels/<P>/subset")
    ap.add_argument("--work-root", type=Path, default=C.WORK_ROOT)
    ap.add_argument("--sizes-only", action="store_true", help="report sizes, write nothing to the panel")
    args = ap.parse_args(argv)

    panel = C.load_panel(args.panel, "full", args.panels_root)
    out = args.out or panel.dir / "subset"
    full = {code: {C.fasta_id(h): (h, s) for h, s in C.read_fasta(p)} for code, p in panel.fastas.items()}
    n_total = sum(len(v) for v in full.values())
    seeds, seed_report = load_seeds(args.seeds, full)
    if not seeds:
        raise C.WrapperError("no seed sequences")
    work = args.work_root / "subset" / args.panel
    if work.exists():
        shutil.rmtree(work)
    C.log(f"subset {args.panel}: {len(seeds)} seeds vs {n_total} proteins")
    hits, cmds = run_diamond(seeds, full, work, args.threads, args.dbsize, args.evalue)
    nb, dom, best = classify_hits(hits, args.min_cov)
    sizes = {"panel_total": n_total, "neighbourhood": len(nb), "domain_only_hits": len(dom)}
    C.log(f"sizes: {sizes}")
    if args.sizes_only:
        print(sizes)
        return 0
    cat, frac, sp_of = select(full, nb, dom, args.seed, args.bg_fraction, args.domain_cap, args.cap)
    if out.exists():
        shutil.rmtree(out)
    out.mkdir(parents=True)
    per_species = {}
    for code, p in panel.fastas.items():
        keep = [(h, s) for rid, (h, s) in full[code].items() if rid in cat]
        n = C.write_fasta(out / f"{code}.fa", keep)
        per_species[code] = {"n": n, "n_full": len(full[code]), "sha256": C.sha256_file(out / f"{code}.fa"),
                             "source": str(p), "source_sha256": C.sha256_file(p)}
    with open(out / "members.tsv", "w") as fh:
        fh.write("id\tspecies\tcategory\tbest_cov\n")
        for rid in sorted(cat, key=lambda r: (sp_of[r], r)):
            bc = best.get(rid)
            fh.write(f"{rid}\t{sp_of[rid]}\t{cat[rid]}\t{'' if bc is None else round(bc, 3)}\n")
    with open(hits, "rb") as fi, gzip.open(out / "diamond_hits.tsv.gz", "wb") as fo:
        shutil.copyfileobj(fi, fo)
    counts = defaultdict(int)
    for c in cat.values():
        counts[c] += 1
    man = {
        "panel": args.panel, "created": time.strftime("%Y-%m-%dT%H:%M:%S%z"), "git_commit": C.git_commit(),
        "seed": args.seed, "seeds": seed_report,
        "params": {k: getattr(args, k) for k in ("evalue", "dbsize", "min_cov", "domain_cap", "bg_fraction", "cap")},
        "background_fraction_used": frac, "diamond_version": C.capture(["diamond", "version"]),
        "commands": cmds, "sizes": sizes, "counts": dict(counts), "n_subset": len(cat),
        "species": per_species, "members": "members.tsv",
        "species_tree_sha256": C.sha256_file(panel.tree_path),
    }
    C.write_json(out / "manifest.json", man)
    C.log(f"subset written to {out}: {dict(counts)} (bg fraction {frac:.4f})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
