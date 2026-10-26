#!/usr/bin/env python3
"""Minimal UniProt ID-mapping REST client (verified pattern). usage: C_uniprot_idmap.py FROM TO id1,id2,...
e.g. RefSeq_Protein UniProtKB NP_001230.1 ; Ensembl_Protein UniProtKB ENSP00000269305 ; EMBL-GenBank-DDBJ_CDS UniProtKB AAA12345.1"""
import sys, json, time, urllib.request, urllib.parse
frm,to,ids=sys.argv[1:4]
r=urllib.request.urlopen(urllib.request.Request("https://rest.uniprot.org/idmapping/run",data=urllib.parse.urlencode({"from":frm,"to":to,"ids":ids}).encode()),timeout=60)
job=json.load(r)["jobId"]
for _ in range(60):
    s=json.load(urllib.request.urlopen(f"https://rest.uniprot.org/idmapping/status/{job}",timeout=60))
    if "results" in s or "failedIds" in s or s.get("jobStatus")=="FINISHED": break
    time.sleep(2)
url=f"https://rest.uniprot.org/idmapping/uniprotkb/results/stream/{job}?format=tsv&fields=accession,id,gene_names,reviewed,length,organism_id" if to.startswith("UniProtKB") else f"https://rest.uniprot.org/idmapping/results/stream/{job}"
print(urllib.request.urlopen(url,timeout=120).read().decode())
