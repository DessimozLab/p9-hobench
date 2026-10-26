#!/usr/bin/env python
"""E06 myosins - motor-domain mapping preparation (DIAMOND, offline); writes cache/fwd.tsv and cache/rev.tsv.

Source sequences: Kollmar & Muehlhausen 2017 (PMID 28870165) figshare 4565155 Myosin_motordomain_alignment.fas
(class-encoded CyMoBase names, motor domains only, ungapped here); species abbreviations from Additional file 2 (xls).
Same-species mapping only (panel E proteome of that species). Run from project root:
python families/E06_myosins/prep_mapping.py
"""
import io, re, subprocess, tempfile, zipfile
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
FAM = Path(__file__).resolve().parent
ALN = ROOT / "literature/papers/PMID28870165/raw/figshare_4565155/Myosin_motordomain_alignment.fas"
XLS = ROOT / "literature/papers/PMID28870165/raw/supp/12862_2017_1056_MOESM2_ESM.xls"
# source abbreviation -> panel code (strain-matched where the source has several strains)
ABMAP = {"Hs": "HUMAN", "Dm": "DROME", "Nv": "NEMVE", "Amq": "AMPQE", "Mnl": "MNELE", "Tia": "TRIAD", "Mb": "MONBE",
         "Co": "CAPO3", "Sc_c": "YEAST", "Nc": "NEUCR", "Bad_b": "BATDJ", "Ec": "ENCCU", "Dd": "DICDI", "Eh": "ENTH1",
         "At": "ARATH", "Cr": "CHLRE", "Thp": "THAPS", "Phr": "PHYRM", "Pf_a": "PLAF7", "Tb": "TRYB2", "Ng": "NAEGR",
         "Tct": "THETR", "Im": "ICHMU"}
COLS = "qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore".split()


def abbrevs():
    x = pd.read_excel(XLS, header=0)
    return set(x.iloc[:, 1].astype(str).str.strip())


def load_seqs():
    seqs, n = {}, None
    for l in open(ALN):
        l = l.strip()
        if l.startswith(">"):
            n = l[1:].strip(); seqs[n] = []
        elif n:
            seqs[n].append(l)
    return {k: "".join(v).replace("-", "").replace(".", "").upper() for k, v in seqs.items()}


_AB = None


def parse(name):
    """-> (abbrev, class, suffix) or None. Name = <abbrev>(Myo|Mhc)<class><variant>; Mhc = class 2;
    Myo_<X> = orphan (class 'orphan')."""
    global _AB
    if _AB is None:
        _AB = abbrevs()
    m = re.match(r"^(.*?)(Myo|Mhc)(.*)$", name)
    if not m or m.group(1) not in _AB:
        return None
    if m.group(2) == "Mhc":
        return (m.group(1), "2", m.group(3))
    c = re.match(r"(\d*)(.*)$", m.group(3))
    return (m.group(1), c.group(1) or "orphan", c.group(2))


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
            mine = {k: s for k, s in seqs.items() if s and parse(k) and parse(k)[0] == ab}
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
