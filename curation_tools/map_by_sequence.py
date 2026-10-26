"""Map paper sequences onto a same-species target proteome (a cross-reference, not inference).

CLI:
  python -m curation_tools.map_by_sequence --query paper_seqs.fasta --target proteome.fasta -o map.tsv
  python -m curation_tools.map_by_sequence --ids ids.tsv --target proteome.fasta -o map.tsv \
         [--save-query fetched.fasta] [--tool diamond|mmseqs] [--min-qcov 0.5] [--threads 8]

--ids accepts a TSV with an 'id' column (e.g. extract_ids output), or one id per line.
Sequences are fetched per id: UniProt accessions from UniProt (obsolete ones via UniParc), RefSeq/GenBank
protein ids from NCBI efetch, anything else is first mapped to UniProt with map_ids.

Algorithm
  1. exact match: sha256 of the normalised sequence (upper-case, no gaps/stop) vs. every target.
     Several identical targets -> status 'exact_ambiguous'.
  2. otherwise best hit by bit score from DIAMOND (default, --very-sensitive) or MMseqs2.
     Reports mapping_identity (pident over the alignment), identity_over_query (identities/query length),
     query_cov, target_cov.
  3. status: exact | exact_ambiguous | unique | ambiguous (another target has pident within
     --ambiguity-delta (1.0) points of the best and query_cov >= --min-qcov) | partial (best hit has
     query_cov < --min-qcov or target_cov < --min-tcov; see 'flags': fragment/domain-level match) | unmapped.
Nothing is dropped: every query appears in the output.
"""
from __future__ import annotations

import argparse
import csv
import os
import re
import shutil
import subprocess
import sys
import tempfile
from collections import defaultdict
from pathlib import Path

from . import common
from .common import HttpError, get, log, norm_seq, read_fasta, seq_sha256, write_fasta

OUT_COLS = ["query_id", "query_len", "target_id", "target_len", "mapping_identity", "identity_over_query",
            "query_cov", "target_cov", "evalue", "bitscore", "method", "status", "flags", "n_exact_targets",
            "ambiguous_with"]


# ----------------------------------------------------------------------------- fetch sequences
def _parse_fasta_text(txt: str) -> list[tuple[str, str]]:
    out, h, s = [], None, []
    for ln in txt.splitlines():
        if ln.startswith(">"):
            if h is not None:
                out.append((h, "".join(s)))
            h, s = ln[1:], []
        elif ln.strip():
            s.append(ln.strip())
    if h is not None:
        out.append((h, "".join(s)))
    return out


def _acc_from_header(h: str) -> str:
    tok = h.split()[0]
    m = re.match(r"^(?:sp|tr)\|([^|]+)\|", tok)
    return m.group(1) if m else tok


def fetch_sequences(ids: list[str], id_types: dict[str, str] | None = None) -> tuple[dict[str, tuple[str, str]], list[str]]:
    """Return ({id: (source, sequence)}, [failed ids]). See module docstring for the routing."""
    from .extract_ids import find_ids

    id_types = dict(id_types or {})
    for i in ids:
        if i not in id_types:
            hs = [h for h in find_ids(i) if h.id == i or h.span == (0, len(i))]
            id_types[i] = hs[0].id_type if hs else "other"
    got: dict[str, tuple[str, str]] = {}
    ncbi = [i for i in ids if id_types[i] in ("refseq_protein", "genbank_protein")]
    uni = [i for i in ids if id_types[i] == "uniprot" or re.fullmatch(r"UPI[0-9A-F]{10}", i)]
    rest = [i for i in ids if i not in ncbi and i not in uni]

    for k in range(0, len(ncbi), 150):
        chunk = ncbi[k:k + 150]
        try:
            r = get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi",
                    params={"db": "protein", "id": ",".join(chunk), "rettype": "fasta", "retmode": "text"})
        except HttpError as e:
            log.error("NCBI efetch failed for %d ids: %s", len(chunk), e)
            continue
        by_acc = {}
        for h, s in _parse_fasta_text(r.text):
            a = _acc_from_header(h)
            by_acc[a] = s
            by_acc[re.sub(r"\.\d+$", "", a)] = s
        for i in chunk:
            s = by_acc.get(i) or by_acc.get(re.sub(r"\.\d+$", "", i))
            if s:
                got[i] = ("NCBI efetch protein", s)

    def uniprot_batch(accs: list[str]):
        for k in range(0, len(accs), 200):
            chunk = accs[k:k + 200]
            try:
                r = get("https://rest.uniprot.org/uniprotkb/accessions",
                        params={"accessions": ",".join(chunk), "format": "fasta"})
            except HttpError as e:
                log.error("UniProt fetch failed: %s", e)
                continue
            for h, s in _parse_fasta_text(r.text):
                got_acc = _acc_from_header(h)
                got[got_acc] = ("UniProtKB", s)
                # secondary accession requests come back under the primary: keep via lookup below

    plain = [i for i in uni if not i.startswith("UPI")]
    uniprot_batch([re.sub(r"-\d+$", "", i) for i in plain if "-" not in i])
    for i in plain:  # isoforms and the ones that were not returned
        if i in got:
            continue
        try:
            r = get(f"https://rest.uniprot.org/uniprotkb/{i}.fasta", ok_status=(200, 404))
            if r.status_code == 200 and r.text.startswith(">"):
                got[i] = ("UniProtKB", _parse_fasta_text(r.text)[0][1])
        except HttpError as e:
            log.error("UniProt fetch %s: %s", i, e)
    # obsolete accessions / UPIs -> UniParc
    for i in uni:
        if i in got:
            continue
        upi = i
        if not i.startswith("UPI"):
            from .map_ids import run_job
            try:
                rows, _ = run_job("UniProtKB_AC-ID", "UniParc", [i])
                upi = rows[0]["To"] if rows else None
            except Exception as e:  # noqa: BLE001
                log.error("UniParc lookup %s: %s", i, e)
                upi = None
        if upi:
            try:
                r = get(f"https://rest.uniprot.org/uniparc/{upi}.fasta", ok_status=(200, 404))
                if r.status_code == 200:
                    got[i] = (f"UniParc {upi}", _parse_fasta_text(r.text)[0][1])
            except HttpError as e:
                log.error("UniParc %s: %s", upi, e)
    if rest:   # map to UniProt first
        from .map_ids import Mapper
        rows = Mapper().map([(i, id_types[i], None) for i in rest])
        accs = {r["original_id"]: r["uniprot_acc"] for r in rows if r["uniprot_acc"]}
        uniprot_batch(sorted(set(accs.values())))
        for i, a in accs.items():
            if a in got:
                got[i] = (f"UniProtKB via {a}", got[a][1])
    failed = [i for i in ids if i not in got]
    return got, failed


# ----------------------------------------------------------------------------- search backends
def _threads(n):
    return n or min(8, os.cpu_count() or 1)


def run_diamond(query: Path, target: Path, tmp: Path, threads: int, sensitivity: str, evalue: float, k: int) -> list[dict]:
    exe = common.have_exe("diamond")
    if not exe:
        raise RuntimeError("diamond not found on PATH")
    db = tmp / "target"
    subprocess.run([exe, "makedb", "--in", str(target), "-d", str(db), "--threads", str(threads), "--quiet"],
                   check=True, capture_output=True)
    out = tmp / "diamond.tsv"
    cmd = [exe, "blastp", "-q", str(query), "-d", str(db), "-o", str(out), "--threads", str(threads),
           f"--{sensitivity}", "-e", str(evalue), "-k", str(k), "--max-hsps", "1", "--quiet",
           "--outfmt", "6", "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend",
           "sstart", "send", "evalue", "bitscore", "qlen", "slen"]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(f"diamond failed: {r.stderr[-500:]}")
    return _parse_hits(out, pident_scale=1.0)


def run_mmseqs(query: Path, target: Path, tmp: Path, threads: int, evalue: float, k: int) -> list[dict]:
    exe = common.have_exe("mmseqs")
    if not exe:
        raise RuntimeError("mmseqs not found on PATH")
    out = tmp / "mmseqs.tsv"
    cmd = [exe, "easy-search", str(query), str(target), str(out), str(tmp / "mm_tmp"), "--threads", str(threads),
           "-s", "7.5", "-e", str(evalue), "--max-seqs", str(max(k, 50)), "-v", "1", "--format-output",
           "query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,qlen,tlen"]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(f"mmseqs failed: {r.stderr[-500:]}")
    return _parse_hits(out, pident_scale=1.0)   # mmseqs 'pident' is already a percentage


def _parse_hits(path: Path, pident_scale: float) -> list[dict]:
    hits = []
    with open(path) as f:
        for ln in f:
            p = ln.rstrip("\n").split("\t")
            if len(p) < 14:
                continue
            q, t = p[0], p[1]
            pid = float(p[2]) * pident_scale
            aln, qs, qe, ts, te = int(p[3]), int(p[6]), int(p[7]), int(p[8]), int(p[9])
            qlen, tlen = int(p[12]), int(p[13])
            nid = pid / 100.0 * aln
            hits.append(dict(query_id=q, target_id=t, mapping_identity=pid, evalue=float(p[10]),
                             bitscore=float(p[11]), query_len=qlen, target_len=tlen,
                             query_cov=(abs(qe - qs) + 1) / qlen, target_cov=(abs(te - ts) + 1) / tlen,
                             identity_over_query=nid / qlen))
    return hits


# ----------------------------------------------------------------------------- main logic
def map_sequences(queries: list[tuple[str, str]], target_fasta: str | Path, tool: str = "diamond",
                  threads: int | None = None, min_qcov: float = 0.5, min_tcov: float = 0.5, ambiguity_delta: float = 1.0,
                  sensitivity: str = "very-sensitive", evalue: float = 1e-5, top_k: int = 10,
                  keep_tmp: str | None = None) -> list[dict]:
    threads = _threads(threads)
    targets = list(read_fasta(target_fasta))
    tid = [h.split()[0] for h, _ in targets]
    tseq = {h.split()[0]: s for h, s in targets}
    if len(set(tid)) != len(tid):
        raise ValueError("duplicate target ids in proteome FASTA; use fetch_proteome's canonical headers")
    exact_index: dict[str, list[str]] = defaultdict(list)
    for t, s in tseq.items():
        exact_index[seq_sha256(s)].append(t)

    results: dict[str, dict] = {}
    todo: list[tuple[str, str]] = []
    for qid, qs in queries:
        ns = norm_seq(qs)
        if not ns:
            results[qid] = _empty(qid, 0, "empty_query_sequence")
            continue
        ex = exact_index.get(seq_sha256(qs))
        if ex:
            t = ex[0]
            results[qid] = dict(query_id=qid, query_len=len(ns), target_id=",".join(ex) if len(ex) > 1 else t,
                                target_len=len(norm_seq(tseq[t])), mapping_identity=100.0, identity_over_query=1.0,
                                query_cov=1.0, target_cov=1.0, evalue=0, bitscore="", method="exact_sha256",
                                status="exact" if len(ex) == 1 else "exact_ambiguous", flags="", n_exact_targets=len(ex),
                                ambiguous_with="")
        else:
            todo.append((qid, ns))
    if todo:
        tmp_ctx = tempfile.TemporaryDirectory(dir=keep_tmp) if not keep_tmp else None
        tmp = Path(tmp_ctx.name if tmp_ctx else keep_tmp)
        try:
            qf = tmp / "q.fasta"
            write_fasta(todo, qf)
            tf = tmp / "t.fasta"
            write_fasta(((t, norm_seq(s)) for t, s in tseq.items()), tf)  # ids only, no description
            hits = (run_diamond(qf, tf, tmp, threads, sensitivity, evalue, top_k) if tool == "diamond"
                    else run_mmseqs(qf, tf, tmp, threads, evalue, top_k))
        finally:
            if tmp_ctx:
                tmp_ctx.cleanup()
        byq: dict[str, list[dict]] = defaultdict(list)
        for h in hits:
            byq[h["query_id"]].append(h)
        for qid, ns in todo:
            hs = sorted(byq.get(qid, []), key=lambda h: (-h["bitscore"], -h["mapping_identity"]))
            if not hs:
                results[qid] = _empty(qid, len(ns), "unmapped")
                continue
            b = hs[0]
            close = [h for h in hs[1:] if h["target_id"] != b["target_id"]
                     and h["mapping_identity"] >= b["mapping_identity"] - ambiguity_delta
                     and h["query_cov"] >= min_qcov]
            flags = [f for f, bad in (("low_query_cov", b["query_cov"] < min_qcov),
                                      ("low_target_cov", b["target_cov"] < min_tcov)) if bad]
            status = "ambiguous" if close else ("partial" if flags else "unique")
            results[qid] = dict(b, method=tool, status=status, flags=";".join(flags), n_exact_targets=0,
                                ambiguous_with=";".join(f"{h['target_id']}({h['mapping_identity']:.1f}%)" for h in close))
    return [results[q] for q, _ in queries if q in results]


def _empty(qid, qlen, status):
    return dict(query_id=qid, query_len=qlen, target_id="", target_len="", mapping_identity="",
                identity_over_query="", query_cov="", target_cov="", evalue="", bitscore="", method="",
                status=status, flags="", n_exact_targets=0, ambiguous_with="")


def write_results(rows: list[dict], out) -> None:
    w = csv.writer(out, delimiter="\t", lineterminator="\n")
    w.writerow(OUT_COLS)
    for r in rows:
        w.writerow([f"{r[c]:.4g}" if isinstance(r.get(c), float) and c in ("evalue",) else
                    (f"{r[c]:.3f}" if isinstance(r.get(c), float) else r.get(c, "")) for c in OUT_COLS])


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--query", help="FASTA of paper sequences")
    g.add_argument("--ids", help="TSV (id column) / list of ids whose sequences are fetched")
    ap.add_argument("--target", required=True, help="proteome FASTA of the same species")
    ap.add_argument("-o", "--out")
    ap.add_argument("--save-query", help="write fetched query sequences here (with --ids)")
    ap.add_argument("--tool", choices=["diamond", "mmseqs"], default="diamond")
    ap.add_argument("--threads", type=int)
    ap.add_argument("--min-qcov", type=float, default=0.5, help="query coverage below this -> flag low_query_cov")
    ap.add_argument("--min-tcov", type=float, default=0.5, help="target coverage below this -> flag low_target_cov")
    ap.add_argument("--ambiguity-delta", type=float, default=1.0, help="pident points (default 1.0)")
    ap.add_argument("-v", "--verbose", action="store_true")
    a = ap.parse_args(argv)
    common.setup_logging(a.verbose)
    if a.query:
        queries = [(h.split()[0], s) for h, s in read_fasta(a.query)]
    else:
        from .map_ids import read_items
        items = read_items(a.ids)
        types = {i: t for i, t, _ in items}
        got, failed = fetch_sequences([i for i, _, _ in items], types)
        queries = [(i, got[i][1]) for i, _, _ in items if i in got]
        for i in failed:
            log.error("could not fetch a sequence for %s", i)
        if a.save_query:
            write_fasta(((f"{i} source={got[i][0]}", got[i][1]) for i, _, _ in items if i in got), a.save_query)
    rows = map_sequences(queries, a.target, a.tool, a.threads, a.min_qcov, a.min_tcov, a.ambiguity_delta)
    if not a.query:
        rows += [_empty(i, 0, "sequence_not_fetched") for i in failed]
    out = open(a.out, "w") if a.out else sys.stdout
    write_results(rows, out)
    if a.out:
        out.close()
    from collections import Counter
    log.info("status counts: %s", dict(Counter(r["status"] for r in rows)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
