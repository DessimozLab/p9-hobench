#!/usr/bin/env python3
"""Fetch assembly metadata (name, level, BUSCO, annotation PCG count) from NCBI Datasets v2 for accessions in panels_verified.tsv"""
import csv, json, urllib.request, sys
D="/scratch/nromashc/agentic-ai/workspaces/mywork/literature/raw/C/"
rows=list(csv.DictReader(open(D+"panels_verified.tsv"),delimiter="\t"))
accs=sorted({r["rest_assembly"] for r in rows if r["rest_assembly"]})
out=csv.writer(open(D+"assembly_meta.tsv","w"),delimiter="\t")
out.writerow(["accession","organism","assembly_name","level","refseq_cat","busco_complete","busco_lineage","annot_name","annot_pcg","annot_total_genes","release"])
for i in range(0,len(accs),50):
    chunk=accs[i:i+50]
    u="https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/"+",".join(chunk)+"/dataset_report"
    d=json.load(urllib.request.urlopen(u,timeout=120))
    for r in d.get("reports",[]):
        ai=r.get("assembly_info",{}); an=r.get("annotation_info",{}); bu=ai.get("busco",{}) or r.get("assembly_stats",{}) 
        b=ai.get("busco",{})
        out.writerow([r["accession"],r["organism"]["organism_name"],ai.get("assembly_name"),ai.get("assembly_level"),ai.get("refseq_category",""),
            b.get("complete",""),b.get("busco_lineage",""),an.get("name",""),an.get("stats",{}).get("gene_counts",{}).get("protein_coding",""),an.get("stats",{}).get("gene_counts",{}).get("total",""),ai.get("release_date","")])
