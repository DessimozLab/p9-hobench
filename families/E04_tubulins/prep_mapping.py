#!/usr/bin/env python
"""E04 tubulins - sequence-mapping preparation (DIAMOND, offline); writes cache/fwd.tsv and cache/rev.tsv.

Source sequences: Findeisen et al. 2014 GBE (PMID 25169981) supplementary tubulin_alignment.fas (full-length,
ungapped here); species abbreviations from DatasetS2.xls. Same-species mapping only (panel E proteome of that species).
Run from project root: python families/E04_tubulins/prep_mapping.py
"""
import io, re, subprocess, tempfile, zipfile
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
FAM = Path(__file__).resolve().parent
ZIP = ROOT / "literature/papers/PMID25169981/raw/supp/supp_evu187_suppl_data.zip"
# source abbreviation -> panel code (strain-matched where the source has several strains)
ABMAP = {"Hs": "HUMAN", "Dm": "DROME", "Nv": "NEMVE", "Amq": "AMPQE", "Tia": "TRIAD", "Mb": "MONBE", "Co": "CAPO3",
         "Sc_c": "YEAST", "Nc": "NEUCR", "Bad_b": "BATDJ", "Ec": "ENCCU", "Dd": "DICDI", "Eh": "ENTH1", "At": "ARATH",
         "Cr": "CHLRE", "Cm": "CYAME", "Thp": "THAPS", "Phr": "PHYRM", "Pf_a": "PLAF7", "Tb": "TRYB2", "Ng": "NAEGR",
         "Gii": "GIAIC", "Tct": "THETR", "Im": "ICHMU"}
COLS = "qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore".split()


def load_seqs():
    z = zipfile.ZipFile(ZIP)
    seqs, n = {}, None
    for l in io.TextIOWrapper(z.open("tubulin_alignment.fas")):
        l = l.strip()
        if l.startswith(">"):
            n = l[1:].strip(); seqs[n] = []
        elif n:
            seqs[n].append(l)
    return {k: "".join(v).replace("-", "").replace(".", "").upper() for k, v in seqs.items()}


def parse(name):
    """-> (abbrev, class digit, suffix) or None. Abbreviation = text before 'Tub<digit>'."""
    m = re.match(r"^(.+?)Tub(\d)(.*)$", name)
    return (m.group(1), m.group(2), m.group(3)) if m else None


def run(q, t, o, td):
    db = Path(td) / (Path(t).stem + ".dmnd")
    subprocess.run(["diamond", "makedb", "--in", str(t), "-d", str(db), "--quiet"], check=True)
    subprocess.run(["diamond", "blastp", "-q", str(q), "-d", str(db), "-o", str(o), "--very-sensitive", "-k", "5",
                    "--quiet", "--tmpdir", str(td), "-f", "6", *COLS], check=True)
    return pd.read_csv(o, sep="\t", names=COLS) if Path(o).stat().st_size else pd.DataFrame(columns=COLS)


def main():
    sp = pd.read_csv(ROOT / "panels/E/species.tsv", sep="\t", dtype=str).set_index("code")
    seqs = load_seqs()
    fw, rv = [], []
    with tempfile.TemporaryDirectory() as td:
        for ab, code in ABMAP.items():
            mine = {k: s for k, s in seqs.items() if parse(k) and parse(k)[0] == ab}
            if not mine:
                continue
            q = Path(td) / f"{code}_q.fa"
            q.write_text("".join(f">{k}\n{s}\n" for k, s in sorted(mine.items())))
            fa = ROOT / sp.loc[code, "fasta_path"]
            d = run(q, fa, Path(td) / f"{code}.tsv", td); d["species"] = code; fw.append(d)
            tg = set(d.sseqid)
            tf = Path(td) / f"{code}_t.fa"
            out, name, buf = [], None, []
            for l in list(open(fa)) + [">END\n"]:
                if l.startswith(">"):
                    if name in tg:
                        out.append(f">{name}\n{''.join(buf)}\n")
                    name, buf = l[1:].split()[0], []
                else:
                    buf.append(l.strip())
            if out:
                tf.write_text("".join(out))
                r = run(tf, q, Path(td) / f"{code}_rev.tsv", td); r["species"] = code; rv.append(r)
    (FAM / "cache").mkdir(exist_ok=True)
    pd.concat(fw).to_csv(FAM / "cache/fwd.tsv", sep="\t", index=False)
    pd.concat(rv).to_csv(FAM / "cache/rev.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
