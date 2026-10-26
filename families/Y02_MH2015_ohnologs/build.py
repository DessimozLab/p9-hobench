#!/usr/bin/env python
"""Y02_MH2015_ohnologs: cross-check of Marcet-Houben & Gabaldon 2015 (PMID26252497) S1 Data sheet 'Fig1C'
(S. cerevisiae ohnolog pairs + inferred topology class) against YGOB v1 pillars and the Y01 units.
Writes crosscheck.tsv. No reference family is produced (see curation_log.md)."""
import csv, hashlib
from pathlib import Path
import xlrd
ROOT = Path(__file__).resolve().parents[2]
S1 = ROOT / "literature/papers/PMID26252497/raw/supp/pbio.1002220.s001.xls"
V1 = ROOT / "literature/raw/B1/ygob/wayback/v1_Pillars.tab"
Y01 = ROOT / "families/Y01_YGOB_pillars/members.tsv"
sh = xlrd.open_workbook(str(S1)).sheet_by_name("Fig1C")
mh = {}
for r in range(1, sh.nrows):
    a, b, t = sh.row_values(r)[:3]
    mh.setdefault(frozenset((a, b)), []).append((r + 1, a, b, t))
v1 = set()
for ln, l in enumerate(V1.read_text().splitlines(), 1):
    f = l.split("\t")
    if f[2] != "---" and f[7] != "---":
        v1.add(frozenset((f[2], f[7])))
unit, yg = {}, set()
for r in csv.DictReader(Y01.open(), delimiter="\t"):
    if r["species"] == "YEAST":
        unit.setdefault(r["unit"], []).append(r["original_id_in_paper"]); yg.add(r["original_id_in_paper"])
ypair = {frozenset(v): u for u, v in unit.items() if len(v) == 2}
with (Path(__file__).parent / "crosscheck.tsv").open("w") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["gene1", "gene2", "mh_rows", "mh_topology", "in_ygob_v1_pillar_pair", "y01_status", "y01_unit"])
    for p, rows in sorted(mh.items(), key=lambda x: sorted(x[0])):
        g1, g2 = sorted(p)
        st = "paired" if p in ypair else ("both_genes_present_not_paired" if p <= yg else "not_in_Y01")
        w.writerow([g1, g2, ";".join(f"Fig1C:{r[0]}" for r in rows), ";".join(sorted({r[3] for r in rows})),
                    p in v1, st, ypair.get(p, "")])
print("sha256 S1", hashlib.sha256(S1.read_bytes()).hexdigest())
