#!/usr/bin/env python
"""Extract COG 2020 rows + protein sequences for the panel-K genomes (curation helper for K01/K02/K03).

Reads raw/cog-20.cog.csv, raw/cog-20.org.csv, raw/cog-20.fa.gz (NCBI COG 2020, Galperin et al. 2021 NAR,
ftp.ncbi.nih.gov/pub/COG/COG2020/data/, md5 verified against raw/checksums.md5.txt) and writes
derived/cog20_panelK.csv (all COG rows of the panel-K assemblies; same columns as cog-20.cog.csv, + panel code)
derived/cog20_panelK.fa (sequences of those proteins; header = COG protein ID as in cog-20.cog.csv).
Panel species -> COG assembly: same organism/strain only (taxid identical or strain-level equivalent, listed below).
"""
import csv, gzip, re
from pathlib import Path
D = Path(__file__).resolve().parent
# panel code -> COG 2020 assembly (cog-20.org.csv). BRADU (B. diazoefficiens USDA 110) and GEOSL (G. sulfurreducens PCA)
# are not in COG 2020 (only B. japonicum USDA 6 and G. metallireducens GS-15 = different species) -> not extracted.
ASM = {
 'ECOLI': 'GCF_000005845.2',  # E. coli K-12 MG1655 (511145) = panel K-12 (83333)
 'PSEAE': 'GCF_000006765.1', 'HELPY': 'GCF_000008525.1', 'BACSU': 'GCF_000009045.1', 'MYCTU': 'GCF_000195955.2',
 'CHLTR': 'GCF_000008725.1', 'SYNY3': 'GCF_000009725.1',  # Synechocystis PCC 6803 (1148) = panel substr. Kazusa (1111708)
 'DEIRA': 'GCF_000008565.1', 'THEMA': 'GCF_000008545.1', 'AQUAE': 'GCF_000008625.1', 'BACTN': 'GCF_000011065.1',
 'LEPIN': 'GCF_000092565.1', 'SACS2': 'GCF_000007005.1', 'KORCO': 'GCF_000019605.1', 'NITMS': 'GCA_000018465.1',
 'THEKO': 'GCF_000009965.1', 'METJA': 'GCF_000091665.1', 'METAC': 'GCA_000007345.1', 'HALSA': 'GCF_000006805.1',
}
inv = {v: k for k, v in ASM.items()}
org = {r[0]: r for r in csv.reader(open(D / 'raw/cog-20.org.csv'))}
assert all(a in org for a in inv), [a for a in inv if a not in org]
prots = set()
(D / 'derived').mkdir(exist_ok=True)
with open(D / 'raw/cog-20.cog.csv') as f, open(D / 'derived/cog20_panelK.csv', 'w', newline='') as o:
    w = csv.writer(o)
    for row in csv.reader(f):
        if row[1] in inv:
            w.writerow(row + [inv[row[1]]]); prots.add(row[2])
keep = False; n = 0
with gzip.open(D / 'raw/cog-20.fa.gz', 'rt') as f, open(D / 'derived/cog20_panelK.fa', 'w') as o:
    for line in f:
        if line.startswith('>'):
            pid = line[1:].split()[0]
            pid = re.sub(r'_(\d+)$', r'.\1', pid)  # fasta header AAM01497_1 == protein id AAM01497.1
            keep = pid in prots
            if keep: o.write(f'>{pid} {line[1:].strip()}\n'); n += 1
        elif keep: o.write(line)
print(len(prots), 'proteins in COG rows;', n, 'sequences written')
