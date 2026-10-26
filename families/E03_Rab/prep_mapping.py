#!/usr/bin/env python
"""E03 Rab GTPases - mapping preparation (DIAMOND); writes cache/ files read by build.py.

Run from project root:  python families/E03_Rab/prep_mapping.py
Source: Kloepper TH, Kienle N, Fasshauer D, Munro S (2012) BMC Biol 10:71 (PMID 22873208), Additional files 4-9 =
alignments of Rab supergroups I-VI (FASTA, CR line ends; names <SpeciesCode>_<RabLabel>[_n], no accessions).
Paper sequences (degapped) of panel species -> DIAMOND vs the same species' panel proteome + reverse search.
"""
import sys
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
FAM = Path(__file__).resolve().parent
CACHE = FAM / "cache"
sys.path.insert(0, str(FAM))
from seqmap import run_diamond  # noqa: E402

SUPP = "literature/papers/PMID22873208/raw/supp"
FILES = {"I": "1741-7007-10-71-S4.TXT", "II": "1741-7007-10-71-S5.TXT", "III": "1741-7007-10-71-S6.TXT",
         "IV": "1741-7007-10-71-S7.TXT", "V": "1741-7007-10-71-S8.TXT", "VI": "1741-7007-10-71-S9.TXT"}
# Additional file 3 short names -> panel E codes
SPMAP = {"HoSa": "HUMAN", "DrMe": "DROME", "NeVe": "NEMVE", "AmQu": "AMPQE", "TrAd": "TRIAD", "MoBr": "MONBE",
         "CaOw": "CAPO3", "SaCe": "YEAST", "ChRe": "CHLRE", "ArTh": "ARATH", "PlFa": "PLAF7", "TrBr": "TRYB2",
         "NaGr": "NAEGR", "ThPs": "THAPS", "PhRa": "PHYRM"}


def load_alignments():
    """-> DataFrame name, sp, group, file, order (1-based record index in file), seq (degapped)."""
    rows = []
    for grp, fn in FILES.items():
        txt = (ROOT / SUPP / fn).read_text().replace("\r\n", "\n").replace("\r", "\n")
        name, buf, k = None, [], 0
        for l in txt.split("\n") + [">"]:
            if l.startswith(">"):
                if name:
                    k += 1
                    s = "".join(buf).replace("-", "").replace(".", "").upper()
                    rows.append(dict(name=name, sp=name.split("_")[0], group=grp, file=fn, order=k, seq=s))
                name, buf = l[1:].strip(), []
            else:
                buf.append(l.strip())
    return pd.DataFrame(rows)


def main():
    CACHE.mkdir(exist_ok=True)
    a = load_alignments()
    fws, rvs = [], []
    for sp5, code in SPMAP.items():
        g = a[a.sp == sp5]
        if g.empty:
            continue
        q = dict(zip(g.name, g.seq))
        fw, rv = run_diamond(q, code, q)
        fws.append(fw); rvs.append(rv)
    pd.concat(fws).to_csv(CACHE / "fw.tsv", sep="\t", index=False)
    pd.concat(rvs).to_csv(CACHE / "rv.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
