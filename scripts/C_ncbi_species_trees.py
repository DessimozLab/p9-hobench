#!/usr/bin/env python3
"""Build NCBI-taxonomy species trees (ete3 get_topology) for each panel in panels_verified.tsv.
Leaves are named by UniProt mnemonic. Reports polytomies. Uses a private taxa.sqlite under literature/raw/C/trees."""
import csv, sys
from ete3 import NCBITaxa
D="/scratch/nromashc/agentic-ai/workspaces/mywork/literature/raw/C/"
ncbi=NCBITaxa(dbfile=D+"trees/taxa.sqlite")
rows=list(csv.DictReader(open(D+"panels_verified.tsv"),delimiter="\t"))
for p in sorted({r["panel"] for r in rows}):
    sp={}
    for r in rows:
        if r["panel"]==p: sp[int(r["taxid"])]=r["mnemonic"]
    # NOTE: ete3 get_topology collapses by taxid; two mnemonics sharing a taxid would collide (none expected)
    # map strain-level taxids to themselves; ete3 handles any valid taxid
    t=ncbi.get_topology(list(sp), intermediate_nodes=False)
    for l in t.iter_leaves():
        l.name=sp.get(int(l.name), l.name)
    npoly=[(n.name,len(n.children)) for n in t.traverse() if not n.is_leaf() and len(n.children)>2]
    t.write(format=9,outfile=D+f"trees/panel_{p}_ncbi_topology.nwk")
    print(p,len(sp),"leaves; polytomies:",len(npoly), [(n.get_leaf_names()[:3],c) for n,c in [(n,len(n.children)) for n in t.traverse() if not n.is_leaf() and len(n.children)>2]][:12])
