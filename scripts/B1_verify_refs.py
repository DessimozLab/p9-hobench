#!/usr/bin/env python3
"""Verify literature references against Europe PMC (B1 survey: animal/fungal families).

Input: text file, one DOI per line (lines starting with '#' ignored; optional
'<tag>\t<doi>').  Output: TSV with the Europe PMC record (title, first author,
year, journal, PMID, PMCID, isOpenAccess, inPMC) or NOT_FOUND.

Usage: python scripts/B1_verify_refs.py dois.txt out.tsv
"""
import sys
import time
import json
import urllib.parse
import requests

API = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"


def query(doi):
    q = f'DOI:"{doi}"'
    url = f"{API}?query={urllib.parse.quote(q)}&format=json&resultType=core&pageSize=3"
    for attempt in range(3):
        try:
            r = requests.get(url, timeout=30)
            r.raise_for_status()
            return r.json().get("resultList", {}).get("result", [])
        except Exception as e:  # network hiccup -> retry
            err = e
            time.sleep(2)
    raise RuntimeError(f"{doi}: {err}")


def main(inp, out):
    rows = []
    for line in open(inp):
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        tag, doi = (line.split("\t", 1) if "\t" in line else ("", line))
        res = query(doi)
        if not res:
            rows.append([tag, doi, "NOT_FOUND"] + [""] * 8)
            continue
        h = res[0]
        first = (h.get("authorList", {}).get("author") or [{}])[0].get("lastName", "")
        rows.append([
            tag, doi, "OK", h.get("pmid", ""), h.get("pmcid", ""),
            first, h.get("pubYear", ""),
            h.get("journalInfo", {}).get("journal", {}).get("medlineAbbreviation", "")
            or h.get("journalInfo", {}).get("journal", {}).get("title", ""),
            h.get("isOpenAccess", ""), h.get("inPMC", ""),
            h.get("title", "").replace("\t", " "),
        ])
        time.sleep(0.15)
    with open(out, "w") as fh:
        fh.write("\t".join(["tag", "doi", "status", "pmid", "pmcid", "first_author",
                            "year", "journal", "isOpenAccess", "inPMC", "title"]) + "\n")
        for r in rows:
            fh.write("\t".join(map(str, r)) + "\n")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
