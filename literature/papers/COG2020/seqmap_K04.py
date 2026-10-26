#!/usr/bin/env python
"""Same-species sequence cross-reference of COG 2020 panel-K proteins of the K04 isofunctional-decoy COGs onto the
panel proteomes (same procedure as seqmap_panelK.py, separate output dir derived/seqmap_K04/ so K01/K02 files are
untouched). Tool: curation_tools.map_by_sequence (exact sha256 first, else DIAMOND --very-sensitive),
--ambiguity-delta 2 --min-qcov 0.8 --min-tcov 0.8."""
import csv, glob, subprocess, sys
from pathlib import Path
D = Path(__file__).resolve().parent; ROOT = D.parents[2]
COGS = {'COG1384', 'COG1190', 'COG0605', 'COG2032', 'COG0588', 'COG0696', 'COG0703', 'COG1685', 'COG0710', 'COG0757'}
fa = {}; k = None
for l in open(D / 'derived/cog20_panelK.fa'):
    if l[0] == '>': k = l[1:].split()[0]; fa[k] = []
    else: fa[k].append(l.strip())
want = {}
for r in csv.reader(open(D / 'derived/cog20_panelK.csv')):
    if r[6] in COGS: want.setdefault(r[-1], set()).add(r[2])
out = D / 'derived/seqmap_K04'; out.mkdir(exist_ok=True)
for r in csv.DictReader(open(ROOT / 'literature/raw/C/panels_verified.tsv'), delimiter='\t'):
    c = r['mnemonic']
    if r['panel'] != 'K' or c not in want: continue
    tgt = glob.glob(str(ROOT / f"data/proteomes/{r['upid']}_*/proteome.fasta"))[0]
    q = out / f'{c}.query.fa'
    with open(q, 'w') as o:
        for p in sorted(want[c]): o.write(f'>{p}\n' + ''.join(fa[p]) + '\n')
    subprocess.run([sys.executable, '-m', 'curation_tools.map_by_sequence', '--query', str(q), '--target', tgt,
                    '-o', str(out / f'{c}.tsv'), '--ambiguity-delta', '2', '--min-qcov', '0.8', '--min-tcov', '0.8',
                    '--threads', '8'], check=True, cwd=ROOT, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
