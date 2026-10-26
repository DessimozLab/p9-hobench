#!/usr/bin/env python3
"""Verify DOIs via Europe PMC REST. Input: file with lines 'label<TAB>doi'. Output TSV: label doi status title journal year firstauthor"""
import sys, json, urllib.request, urllib.parse, time
def q(doi):
    url="https://www.ebi.ac.uk/europepmc/webservices/rest/search?"+urllib.parse.urlencode({"query":f'DOI:"{doi}"',"format":"json","resultType":"lite","pageSize":3})
    for i in range(3):
        try:
            with urllib.request.urlopen(url,timeout=60) as r: return json.load(r)["resultList"]["result"]
        except Exception: time.sleep(2)
    return []
for line in open(sys.argv[1]):
    if not line.strip() or line.startswith("#"): continue
    label,doi=line.rstrip("\n").split("\t")[:2]
    r=q(doi)
    if r:
        x=r[0]; print("\t".join([label,doi,"VERIFIED",x.get("title","")[:110],x.get("journalTitle",x.get("bookOrReportDetails",{}).get("publisher","preprint")),str(x.get("pubYear")),x.get("authorString","").split(",")[0]]))
    else: print("\t".join([label,doi,"NOT_FOUND_IN_EPMC","","","",""]))
