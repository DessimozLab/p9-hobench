#!/usr/bin/env python
"""E05_MCM - mapping preparation (network + DIAMOND); writes cache/query.fasta, cache/query_meta.tsv, cache/fwd.tsv,
cache/rev.tsv. build.py only reads these caches.

Sequences: NCBI protein accessions printed in Liu et al. 2009 Additional files 1-4 (RefSeq/GenBank/EMBL/Swiss-Prot)
are fetched with NCBI efetch. Chlamydomonas 'Chlre3:<n>' JGI ids are looked up as NCBI old locus tag
'CHLREDRAFT_<n>' (esearch) and the returned record's sequence is fetched. Other JGI/CmerDB/GiardiaDB ids are mapped
by identifier in build.py (no sequence needed). Same-species DIAMOND only.
Run from project root: python families/E05_MCM/prep_mapping.py
"""
import re, subprocess, sys, tempfile, time, urllib.parse, urllib.request
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
FAM = Path(__file__).resolve().parent
sys.path.insert(0, str(FAM))
from extract import all_leaves, header_index, id_pattern  # noqa: E402

EU = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
COLS = "qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore".split()
NCBI_RE = re.compile(r"^(?:[NX]P_?\d+|[A-Z]{3}\d{5}|[A-Z]\d{5})$")


def get(url):
    for k in range(5):
        try:
            return urllib.request.urlopen(url, timeout=60).read().decode()
        except Exception:
            time.sleep(2 + 3 * k)
    raise RuntimeError(url)


def ncbi_acc(a):
    a = a.split("(")[0].strip()
    a = re.sub(r"^([NX]P)(\d)", r"\1_\2", a)  # 'XP827264' typo -> XP_827264
    return a if NCBI_RE.match(a) else None


def fasta_records(txt):
    out, n = {}, None
    for l in txt.splitlines():
        if l.startswith(">"):
            n = l[1:].split()[0]; out[n] = [l[1:]]
        elif n:
            out[n].append(l.strip())
    return {k: (v[0], "".join(v[1:])) for k, v in out.items()}


def run(q, t, o, td):
    db = Path(td) / (Path(t).stem + ".dmnd")
    subprocess.run(["diamond", "makedb", "--in", str(t), "-d", str(db), "--quiet"], check=True)
    subprocess.run(["diamond", "blastp", "-q", str(q), "-d", str(db), "-o", str(o), "--very-sensitive", "-k", "5",
                    "--quiet", "--tmpdir", str(td), "-f", "6", *COLS], check=True)
    return pd.read_csv(o, sep="\t", names=COLS) if Path(o).stat().st_size else pd.DataFrame(columns=COLS)


def main():
    (FAM / "cache").mkdir(exist_ok=True)
    rows = [d for d in all_leaves() if d["code"]]
    meta = []
    for d in rows:
        a = ncbi_acc(d["acc"])
        m = re.match(r"Chlre3:(\d+)", d["acc"])
        if a:
            meta.append(dict(orig=d["acc"], code=d["code"], ncbi=a, route="ncbi_accession"))
        elif m:
            js = get(EU + "esearch.fcgi?db=protein&retmode=json&term=" + urllib.parse.quote(f"CHLREDRAFT_{m.group(1)}"))
            ids = re.findall(r'"idlist":\[([^\]]*)\]', js)[0].replace('"', "").split(",")
            ids = [i for i in ids if i]
            meta.append(dict(orig=d["acc"], code=d["code"], ncbi=ids[0] if len(ids) == 1 else "",
                             route=f"ncbi_old_locus_tag CHLREDRAFT_{m.group(1)} ({len(ids)} hits)"))
            time.sleep(0.4)
    md = pd.DataFrame(meta).drop_duplicates()
    ids = sorted({x for x in md.ncbi if x and not x.isdigit()})
    recs = {}
    for i in range(0, len(ids), 100):
        recs.update(fasta_records(get(EU + "efetch.fcgi?db=protein&rettype=fasta&id=" + ",".join(ids[i:i + 100]))))
        time.sleep(0.4)
    byacc = {k.split(".")[0]: (k, h, s) for k, (h, s) in recs.items()}
    for gi in sorted({x for x in md.ncbi if x.isdigit()}):  # esearch gi -> fetched one at a time
        (k, (h, s)), = fasta_records(get(EU + "efetch.fcgi?db=protein&rettype=fasta&id=" + gi)).items()
        byacc[gi] = (k, h, s)
        time.sleep(0.4)
    fetched, qrows = [], []
    for r in md.itertuples():
        hit = byacc.get(r.ncbi)
        if hit:
            qid = f"{r.code}|{r.orig}"
            fetched.append(f">{qid}\n{hit[2]}\n")
            qrows.append(dict(qid=qid, orig=r.orig, code=r.code, route=r.route, ncbi_record=hit[0], ncbi_header=hit[1],
                              qlen=len(hit[2])))
        else:
            qrows.append(dict(qid="", orig=r.orig, code=r.code, route=r.route, ncbi_record="", ncbi_header="not fetched",
                              qlen=0))
    (FAM / "cache/query.fasta").write_text("".join(sorted(set(fetched))))
    pd.DataFrame(qrows).to_csv(FAM / "cache/query_meta.tsv", sep="\t", index=False)
    sp = pd.read_csv(ROOT / "panels/E/species.tsv", sep="\t", dtype=str).set_index("code")
    allq = fasta_records((FAM / "cache/query.fasta").read_text())
    fw, rv = [], []
    with tempfile.TemporaryDirectory() as td:
        for code in sorted({k.split("|")[0] for k in allq}):
            q = Path(td) / f"{code}_q.fa"
            q.write_text("".join(f">{k}\n{s}\n" for k, (h, s) in allq.items() if k.startswith(code + "|")))
            fa = ROOT / sp.loc[code, "fasta_path"]
            d = run(q, fa, Path(td) / f"{code}.tsv", td); d["species"] = code; fw.append(d)
            tg = set(d.sseqid)
            tseq = {k: s for k, (h, s) in fasta_records(fa.read_text()).items() if k in tg}
            if tseq:
                tf = Path(td) / f"{code}_t.fa"
                tf.write_text("".join(f">{k}\n{s}\n" for k, s in tseq.items()))
                r = run(tf, q, Path(td) / f"{code}_rev.tsv", td); r["species"] = code; rv.append(r)
    pd.concat(fw).to_csv(FAM / "cache/fwd.tsv", sep="\t", index=False)
    pd.concat(rv).to_csv(FAM / "cache/rev.tsv", sep="\t", index=False)
    qc()
    idmap()


def idmap():
    """UniProt ID mapping (curation_tools.map_ids) of every NCBI accession -> cache/idmap.tsv (gene names)."""
    qm = pd.read_csv(FAM / "cache/query_meta.tsv", sep="\t", dtype=str).fillna("")
    acc = sorted({r.ncbi_record.split(".")[0] for r in qm.itertuples() if r.route == "ncbi_accession" and r.ncbi_record})
    with tempfile.TemporaryDirectory() as td:
        f = Path(td) / "ids.tsv"
        f.write_text("id\tid_type\n" + "".join(f"{a}\t{'refseq_protein' if '_' in a else 'genbank_protein'}\n" for a in acc))
        subprocess.run([sys.executable, "-m", "curation_tools.map_ids", str(f), "-o", str(FAM / "cache/idmap.tsv")],
                       check=True, cwd=ROOT)


def qc():
    """Homology QC (family membership, not paralogue label): every candidate panel target (identifier-mapped or best
    DIAMOND hit) is searched against the paper's human MCM sequences (HUMAN|NP_* in cache/query.fasta); build.py
    requires e-value <= 1e-20, the paper's own BLASTp cut-off for MCM homologues."""
    sp = pd.read_csv(ROOT / "panels/E/species.tsv", sep="\t", dtype=str).set_index("code")
    fw = pd.read_csv(FAM / "cache/fwd.tsv", sep="\t")
    cand = {(r.species, r.sseqid) for r in fw.itertuples()}
    idx = {}
    for d in all_leaves():
        tok, _ = id_pattern(d["code"], d["acc"]) if d["code"] else (None, None)
        if tok:
            if d["code"] not in idx:
                idx[d["code"]] = header_index(d["code"])[0]
            cand |= {(d["code"], a) for a in idx[d["code"]].get(tok, ())}
    allq = fasta_records((FAM / "cache/query.fasta").read_text())
    with tempfile.TemporaryDirectory() as td:
        hq = Path(td) / "human.fa"
        hq.write_text("".join(f">{k}\n{s}\n" for k, (h, s) in allq.items() if k.startswith("HUMAN|")))
        tf = Path(td) / "cand.fa"
        out = []
        for code in sorted({c for c, _ in cand}):
            want = {a for c, a in cand if c == code}
            out += [f">{code}|{k}\n{s}\n" for k, (h, s) in fasta_records((ROOT / sp.loc[code, "fasta_path"]).read_text()).items()
                    if k in want]
        tf.write_text("".join(out))
        d = run(tf, hq, Path(td) / "qc.tsv", td)
    d = d.sort_values("evalue").groupby("qseqid").head(1)
    d.to_csv(FAM / "cache/qc.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
