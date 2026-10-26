#!/usr/bin/env python3
"""Verify panel candidates: UniProt REST (type, taxon, assembly), RP 2026_03 FTP-table canonical counts,
QfO 2026_02 membership/counts, OMA All.May2026 presence, Ensembl division presence. Writes panels_verified.tsv"""
import csv, json, urllib.request, time, sys
D="/scratch/nromashc/agentic-ai/workspaces/mywork/literature/raw/C/"
rp={r[0]:r for r in csv.reader(open(D+"UniProt_RP_2026_03_species.tsv"),delimiter="\t")}
qfo={r[0]:r for r in csv.reader(open(D+"QfO_2026_02_species.tsv"),delimiter="\t")}
oma=json.load(open(D+"oma_genomes_All.May2026.json")); omat={}
for g in oma: omat.setdefault(g["taxon_id"],[]).append(g["code"])
ens={}
for div in ["Vertebrates","Plants","Fungi","Protists","Metazoa"]:
    for s in json.load(open(D+f"ensembl_species_Ensembl{div}.json"))["species"]:
        ens.setdefault(int(s["taxon_id"]),[]).append(div[:4]+":"+s["name"])
def rest(upid):
    u=f"https://rest.uniprot.org/proteomes/{upid}?format=json"
    for i in range(3):
        try:
            with urllib.request.urlopen(u,timeout=60) as r: return json.load(r)
        except Exception as e: time.sleep(2)
    return {}
cache={}
out=csv.writer(open(D+"panels_verified.tsv","w"),delimiter="\t")
out.writerow(["panel","mnemonic","name","taxid","upid","rest_type","rest_taxid","rest_assembly","rest_total_entries","rp_canonical","rp_additional","QfO2026_02","qfo_canonical","OMA_codes_by_taxid","Ensembl_by_taxid","tier","why"])
for r in csv.DictReader(open(D+"panel_candidates.tsv"),delimiter="\t"):
    u=r["upid"]
    if u not in cache: cache[u]=rest(u)
    j=cache[u]
    row=rp.get(u); q=qfo.get(u)
    out.writerow([r["panel"],r["mnemonic"],r["name"],r["taxid"],u,j.get("proteomeType","?"),j.get("taxonomy",{}).get("taxonId","?"),
      j.get("genomeAssembly",{}).get("assemblyId",""),j.get("proteinCount","?"),row[4] if row else "NOT_IN_RP",row[5] if row else "",
      "yes" if q else "no",q[4] if q else "",",".join(omat.get(int(r["taxid"]),[])) or "-",";".join(ens.get(int(r["taxid"]),[]))[:60] or "-",r["tier"],r["why"]])
