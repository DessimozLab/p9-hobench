#!/usr/bin/env python3
"""Render markdown tables for each panel from literature/raw/C/panels_verified.tsv (core tier then optional tier) + totals."""
import csv, sys
D="/scratch/nromashc/agentic-ai/workspaces/mywork/literature/raw/C/"
rows=list(csv.DictReader(open(D+"panels_verified.tsv"),delimiter="\t"))
def canon(r):
    if r["upid"]=="UP000005640": return 20652   # RP-table value 147503 is wrong; FASTA has 20652 (verified by streaming count)
    return int(r["qfo_canonical"] or r["rp_canonical"])
def fmt(r):
    q="QfO" if r["QfO2026_02"]=="yes" else "RP"
    oma=r["OMA_codes_by_taxid"]; ens="Y" if r["Ensembl_by_taxid"]!="-" else "-"
    return f"| {r['mnemonic']} | {r['name']} | {r['taxid']} | {r['upid']} | {q} | {canon(r):,} | {oma} | {ens} | {r['why']} |"
names={"V":"Vertebrata","Y":"Saccharomycetaceae / budding yeasts","P":"Viridiplantae","E":"Eukaryota-wide (+ archaea/bacteria outgroups)","K":"Prokaryotes"}
out=[]
for p in "VYPEK":
    pr=[r for r in rows if r["panel"]==p]
    out.append(f"\n#### Panel {p}: {names[p]}\n")
    for tier in ("core","optional"):
        sub=[r for r in pr if r["tier"]==tier]
        if not sub: continue
        tot=sum(canon(r) for r in sub)
        out.append(f"**{tier.capitalize()}: {len(sub)} proteomes, {tot:,} canonical proteins.**\n")
        out.append("| Code | Species | TaxID | UniProt proteome | Set | Canonical proteins | OMA code (exact taxid) | Ensembl (exact taxid) | Why / note |")
        out.append("|---|---|---|---|---|---|---|---|---|")
        out += [fmt(r) for r in sub]
        out.append("")
print("\n".join(out))
