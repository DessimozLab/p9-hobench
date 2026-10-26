#!/usr/bin/env python
"""E01 kinesins - mapping preparation (network + DIAMOND); writes cache/ files read by build.py.

Run from project root:  python families/E01_kinesins/prep_mapping.py
1. reads Wickstead et al. 2010 Additional file 3 (S3.XLS) and Additional file 8 (DatasetS1.fas, motor-domain alignment)
2. GenBank/RefSeq accessions of panel-species rows -> UniProt via curation_tools.map_ids  (cache/idmap.tsv)
3. motor-domain sequences of panel-species rows -> DIAMOND blastp vs the same species' panel proteome
   (cache/motor_vs_panel.tsv, top hits) and reverse (panel hits -> all motor domains of that species; cache/reverse.tsv)
"""
import io, subprocess, sys, zipfile, tempfile
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
FAM = ROOT / "families/E01_kinesins"
CACHE = FAM / "cache"
PAP = ROOT / "literature/papers/PMID20423470/raw/supp"
sys.path.insert(0, str(ROOT))

SPMAP = {"Homsa": "HUMAN", "Drome": "DROME", "Nemve": "NEMVE", "Triad": "TRIAD", "Monbr": "MONBE",
         "Sacce": "YEAST", "Neucr": "NEUCR", "Batde": "BATDJ", "Enccu": "ENCCU", "Dicdi": "DICDI", "Enthi": "ENTH1",
         "Arath": "ARATH", "Chlre": "CHLRE", "Cyame": "CYAME", "Thaps": "THAPS", "Plafa": "PLAF7", "Trybr": "TRYB2",
         "Naegr": "NAEGR", "Giala": "GIAIC"}


def load_table():
    x = pd.read_excel(PAP / "1471-2148-10-110-S3.XLS", header=1, dtype=str)
    x["xls_row"] = x.index + 3  # 1-based sheet row (2 header rows)
    x["sp5"] = x["Analysis ID"].str.split("_").str[0]
    return x


def load_motor():
    seqs, n = {}, None
    with zipfile.ZipFile(PAP / "1471-2148-10-110-S8.ZIP") as z:
        name = [f for f in z.namelist() if f.endswith(".fas")][0]
        for l in io.TextIOWrapper(z.open(name)):
            l = l.strip()
            if l.startswith(">"):
                n = l[1:].split()[0]; seqs[n] = []
            elif n:
                seqs[n].append(l.replace("-", "").replace(".", "").upper())
    return {k: "".join(v) for k, v in seqs.items()}


def motor_key(r, motor):
    k1 = r["Analysis ID"]
    k2 = r["sp5"] + "_" + str(r["Genome project ID"]).replace("-", "_")
    return k1 if k1 in motor else (k2 if k2 in motor else "")


def main():
    CACHE.mkdir(exist_ok=True)
    sp = pd.read_csv(ROOT / "panels/E/species.tsv", sep="\t", dtype=str).set_index("code")
    x = load_table()
    x = x[x.sp5.isin(SPMAP)].copy()
    x["species"] = x.sp5.map(SPMAP)
    # 1. ID mapping of GenBank/RefSeq accessions
    ids = []
    for r in x.to_dict("records"):
        acc = r["GenBank Acc.No."] if isinstance(r["GenBank Acc.No."], str) else ""
        acc = acc.strip()
        if not acc:
            continue
        t = "refseq_protein" if acc[:3] in ("NP_", "XP_", "YP_") else ("uniprot" if acc[0] in "OPQ" and "." in acc and acc[1].isdigit() else "genbank_protein")
        if t == "uniprot":
            acc = acc.split(".")[0]
        ids.append((acc, t, sp.loc[r["species"], "taxid"]))
    inp = CACHE / "idmap_input.tsv"
    pd.DataFrame(ids, columns=["id", "id_type", "taxon"]).drop_duplicates().to_csv(inp, sep="\t", index=False)
    if not (CACHE / "idmap.tsv").exists(): subprocess.run([sys.executable, "-m", "curation_tools.map_ids", str(inp), "-o", str(CACHE / "idmap.tsv")],
                   check=True, cwd=ROOT)
    # 2. motor-domain DIAMOND vs same-species panel proteome
    motor = load_motor()
    x["motor_key"] = [motor_key(r, motor) for r in x.to_dict("records")]
    fw, rv = [], []
    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        for code, g in x.groupby("species"):
            q = td / f"{code}_q.fa"
            with open(q, "w") as fh:
                for aid in sorted(set(g["motor_key"]) - {""}):
                    fh.write(f">{aid}\n{motor[aid]}\n")
            # all motor domains of this source species (incl. non-panel-row duplicates) for reverse search
            sp5 = g.sp5.iloc[0]
            qa = td / f"{code}_all.fa"
            with open(qa, "w") as fh:
                for k, s in motor.items():
                    if k.split("_")[0] == sp5:
                        fh.write(f">{k}\n{s}\n")
            db = td / f"{code}.dmnd"
            subprocess.run(["diamond", "makedb", "--in", str(ROOT / sp.loc[code, "fasta_path"]), "-d", str(db), "--quiet"], check=True)
            cols = "qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore".split()
            o = td / f"{code}.tsv"
            subprocess.run(["diamond", "blastp", "-q", str(q), "-d", str(db), "-o", str(o), "--very-sensitive",
                            "-k", "5", "--quiet", "--tmpdir", str(td), "-f", "6", *cols], check=True)
            if o.stat().st_size:
                d = pd.read_csv(o, sep="\t", names=cols); d["species"] = code; fw.append(d)
            # reverse: hit targets vs all source motor domains of that species
            if o.stat().st_size:
                tg = set(pd.read_csv(o, sep="\t", names=cols).sseqid)
            else:
                tg = set()
            tf = td / f"{code}_t.fa"
            n = 0
            with open(tf, "w") as fh:
                name, buf = None, []
                for l in open(ROOT / sp.loc[code, "fasta_path"]):
                    if l.startswith(">"):
                        if name in tg:
                            fh.write(f">{name}\n{''.join(buf)}\n"); n += 1
                        name, buf = l[1:].split()[0], []
                    else:
                        buf.append(l.strip())
                if name in tg:
                    fh.write(f">{name}\n{''.join(buf)}\n"); n += 1
            if n:
                dbq = td / f"{code}_all.dmnd"
                subprocess.run(["diamond", "makedb", "--in", str(qa), "-d", str(dbq), "--quiet"], check=True)
                o2 = td / f"{code}_rev.tsv"
                subprocess.run(["diamond", "blastp", "-q", str(tf), "-d", str(dbq), "-o", str(o2), "--very-sensitive",
                                "-k", "5", "--quiet", "--tmpdir", str(td), "-f", "6", *cols], check=True)
                if o2.stat().st_size:
                    d = pd.read_csv(o2, sep="\t", names=cols); d["species"] = code; rv.append(d)
    pd.concat(fw).to_csv(CACHE / "motor_vs_panel.tsv", sep="\t", index=False)
    pd.concat(rv).to_csv(CACHE / "reverse.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
