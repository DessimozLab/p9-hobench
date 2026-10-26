#!/usr/bin/env python3
import sys, json, urllib.request, urllib.parse, re
for doi in sys.argv[1:]:
    u="https://www.ebi.ac.uk/europepmc/webservices/rest/search?"+urllib.parse.urlencode({"query":f'DOI:"{doi}"',"format":"json","resultType":"core","pageSize":1})
    r=json.load(urllib.request.urlopen(u,timeout=60))["resultList"]["result"]
    if not r: print(doi,"NOT FOUND"); continue
    x=r[0]; print("##",doi,x.get("pubYear"),x.get("journalTitle"),x.get("pmcid",""),"\n",x.get("title"),"\n",re.sub(r'<[^>]+>','',x.get("abstractText","(no abstract)")),"\n")
