"""Same-species sequence cross-reference helper (CURATION_PROTOCOL section 0.5).

run_diamond(queries, code, cache_prefix): DIAMOND (--very-sensitive) of paper sequences of one species against that
species' panel proteome (panels/E/proteomes/<CODE>.fa), plus the reverse search of every hit target against all paper
sequences of the same species (reciprocal-best-hit check). Results are cached as TSV; build.py reads only the cache.
best_hits(...) applies the thresholds: identity >= ID_MIN, coverage >= 80 % of both sequences, reciprocal best hit,
margin >= 2 identity points over the runner-up target (runner-up = best other target with query coverage >= 50 %).
"""
import subprocess, tempfile
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
COLS = "qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore".split()
ID_MIN, COV_MIN, MARGIN = 98.0, 0.80, 2.0


def read_fasta(path):
    seqs, name, buf = {}, None, []
    for l in open(path):
        if l.startswith(">"):
            if name is not None:
                seqs[name] = "".join(buf)
            name, buf = l[1:].split()[0], []
        else:
            buf.append(l.strip())
    if name is not None:
        seqs[name] = "".join(buf)
    return seqs


def run_diamond(queries: dict, code: str, all_species_seqs: dict, panel="E"):
    """queries / all_species_seqs: {id: seq}. Returns (forward_df, reverse_df)."""
    target = ROOT / f"panels/{panel}/proteomes/{code}.fa"
    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        q = td / "q.fa"
        q.write_text("".join(f">{k}\n{v}\n" for k, v in queries.items() if v))
        db = td / "t.dmnd"
        subprocess.run(["diamond", "makedb", "--in", str(target), "-d", str(db), "--quiet"], check=True)
        o = td / "fw.tsv"
        subprocess.run(["diamond", "blastp", "-q", str(q), "-d", str(db), "-o", str(o), "--very-sensitive", "-k", "10",
                        "--quiet", "--tmpdir", str(td), "-f", "6", *COLS], check=True)
        fw = pd.read_csv(o, sep="\t", names=COLS) if o.stat().st_size else pd.DataFrame(columns=COLS)
        tg = set(fw.sseqid)
        tseq = {k: v for k, v in read_fasta(target).items() if k in tg}
        rv = pd.DataFrame(columns=COLS)
        if tseq:
            tf = td / "tg.fa"
            tf.write_text("".join(f">{k}\n{v}\n" for k, v in tseq.items()))
            qa = td / "qa.fa"
            qa.write_text("".join(f">{k}\n{v}\n" for k, v in all_species_seqs.items() if v))
            dbq = td / "qa.dmnd"
            subprocess.run(["diamond", "makedb", "--in", str(qa), "-d", str(dbq), "--quiet"], check=True)
            o2 = td / "rv.tsv"
            subprocess.run(["diamond", "blastp", "-q", str(tf), "-d", str(dbq), "-o", str(o2), "--very-sensitive",
                            "-k", "10", "--quiet", "--tmpdir", str(td), "-f", "6", *COLS], check=True)
            if o2.stat().st_size:
                rv = pd.read_csv(o2, sep="\t", names=COLS)
    fw["species"] = code
    rv["species"] = code
    return fw, rv


def best_hits(fw: pd.DataFrame, rv: pd.DataFrame, id_min=ID_MIN):
    """-> {query: dict(t, pid, qcov, tcov, margin, rbh, ok)} using the best (highest bitscore) HSP per pair."""
    out = {}
    if fw.empty:
        return out
    fw = fw.sort_values("bitscore", ascending=False).drop_duplicates(["qseqid", "sseqid"])
    rvb = {}
    if not rv.empty:
        r = rv.sort_values(["qseqid", "bitscore"], ascending=[True, False])
        for t, g in r.groupby("qseqid"):
            top = g.bitscore.max()
            rvb[t] = set(g[g.bitscore >= top - 1e-9].sseqid)
    for q, g in fw.groupby("qseqid"):
        g = g.sort_values("bitscore", ascending=False)
        b = g.iloc[0]
        qcov, tcov = b.length / b.qlen, b.length / b.slen
        oth = g[(g.sseqid != b.sseqid) & (g.length / g.qlen >= 0.5)]
        margin = float(b.pident - (oth.pident.max() if len(oth) else 0.0))
        rbh = q in rvb.get(b.sseqid, set())
        ok = b.pident >= id_min and qcov >= COV_MIN and tcov >= COV_MIN and rbh and margin >= MARGIN
        out[q] = dict(t=b.sseqid, pid=float(b.pident), qcov=float(qcov), tcov=float(tcov), margin=margin, rbh=bool(rbh),
                      ok=bool(ok), second=(oth.iloc[0].sseqid if len(oth) else ""))
    return out
