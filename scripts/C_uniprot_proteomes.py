#!/usr/bin/env python3
"""Query UniProt proteomes REST API for a list of NCBI taxids; print all proteomes (ref/non-ref)
with upid, type, protein_count (all-proteome entries, incl. non-canonical? see note), BUSCO, CPD, assembly.
Usage: C_uniprot_proteomes.py taxid [taxid ...]   (or --file taxids.txt)"""
import sys, json, urllib.request, urllib.parse, time
def q(taxid):
    url = "https://rest.uniprot.org/proteomes/search?" + urllib.parse.urlencode(
        {"query": f"taxonomy_id:{taxid}", "format": "json", "size": 50})
    for i in range(3):
        try:
            with urllib.request.urlopen(url, timeout=60) as r:
                return json.load(r)["results"]
        except Exception as e:
            time.sleep(2)
    return []
def main(ids):
    print("\t".join(["query_taxid","upid","organism","taxid","type","protein_count","busco_complete","cpd","assembly","modified","isRedundant","excluded"]))
    for t in ids:
        for r in q(t):
            b = r.get("proteomeCompletenessReport", {}).get("busco", {})
            bc = ""
            if b: bc = f"C:{b.get('complete')}%[{b.get('lineageDbVersion','')}]" if 'complete' in b else str(b)
            cp = r.get("proteomeCompletenessReport", {}).get("cpdReport", {}).get("status", "")
            ga = r.get("genomeAssembly", {}).get("assemblyId", "")
            print("\t".join(str(x) for x in [t, r["id"], r["taxonomy"]["scientificName"], r["taxonomy"]["taxonId"],
                  r.get("proteomeType"), r.get("proteinCount"), bc, cp, ga, r.get("modified"), r.get("redundantTo",""), r.get("excludedReason","")]))
if __name__ == "__main__":
    a = sys.argv[1:]
    if a and a[0] == "--file": a = [l.split()[0] for l in open(a[1]) if l.strip() and not l.startswith("#")]
    main(a)
