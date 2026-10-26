#!/usr/bin/env python3
"""Search Europe PMC by title/keywords; print top hits (doi, year, journal, title). usage: C_epmc_title.py 'query' [n]"""
import sys, json, urllib.request, urllib.parse
q=sys.argv[1]; n=int(sys.argv[2]) if len(sys.argv)>2 else 4
u="https://www.ebi.ac.uk/europepmc/webservices/rest/search?"+urllib.parse.urlencode({"query":q,"format":"json","resultType":"lite","pageSize":n})
for x in json.load(urllib.request.urlopen(u,timeout=60))["resultList"]["result"]:
    print(x.get("doi"),x.get("pubYear"),x.get("journalTitle"),"|",x.get("title","")[:100],"|",x.get("authorString","")[:30])
