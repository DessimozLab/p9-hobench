"""Map identifiers to UniProtKB accessions with the UniProt ID-mapping REST API.

CLI:
  python -m curation_tools.map_ids ids.tsv -o mapped.tsv [--taxon 559292] [--id-col id] [--type-col id_type]
  python -m curation_tools.map_ids --id NP_001234.1 --id ENSDARP00000012345 --id YAL001C --taxon 559292

Input: a TSV with a header containing an id column (default 'id', e.g. the output of
`extract_ids`; optional columns id_type, taxon) or a plain one-id-per-line file, plus any number of --id.
Identifiers of the same kind are batched into one job (default 1000/job); jobs are polled until done.
Each id type has a cascade of from-databases; ids still unmapped after the cascade are written
with mapping_method=unmapped (never dropped).

Output TSV: original_id, id_type, uniprot_acc, taxid, gene_name, reviewed, length, entry_name,
            mapping_method, n_hits, note, uniparc
Rows for ids that map to several UniProtKB entries are all kept (n_hits>1; reviewed first).
A sidecar <out>.meta.json records the UniProt release and date.
"""
from __future__ import annotations

import argparse
import csv
import io
import json
import re
import sys
import time
from collections import OrderedDict, defaultdict

from . import common
from .common import HttpError, get, log, post

API = "https://rest.uniprot.org"
FIELDS = "accession,id,reviewed,gene_primary,organism_id,length"
TSV_COLS = ["From", "Entry", "Entry Name", "Reviewed", "Gene Names (primary)", "Organism (ID)", "Length"]
BATCH = 1000

# id_type -> cascade of (from_db, transform) tried in order on the still-unmapped ids.
# transform: 'asis' | 'unversioned' | 'upper'
CASCADE: dict[str, list[tuple[str, str]]] = {
    "refseq_protein": [("RefSeq_Protein", "asis"), ("RefSeq_Protein", "unversioned"),
                       ("UniProtKB_AC-ID", "unversioned")],
    "genbank_protein": [("EMBL-GenBank-DDBJ_CDS", "asis"), ("EMBL-GenBank-DDBJ_CDS", "unversioned"),
                        ("RefSeq_Protein", "unversioned")],
    "ensembl_protein": [("Ensembl_Protein", "unversioned"), ("Ensembl_Genomes_Protein", "unversioned")],
    "ensembl_transcript": [("Ensembl_Transcript", "unversioned"), ("Ensembl_Genomes_Transcript", "unversioned")],
    "ensembl_gene": [("Ensembl", "unversioned"), ("Ensembl_Genomes", "unversioned")],
    "refseq_nucleotide": [("RefSeq_Nucleotide", "asis"), ("RefSeq_Nucleotide", "unversioned")],
    "sgd_systematic": [("SGD", "asis"), ("Gene_Name", "asis")],
    "sgd_id": [("SGD", "asis")],
    "arabidopsis_agi": [("Araport", "agi"), ("Gene_Name", "agi")],
    "flybase_gene": [("FlyBase", "asis")],
    "flybase_protein": [("Ensembl_Genomes_Protein", "asis"), ("Ensembl_Genomes_Transcript", "asis")],
    "flybase_transcript": [("Ensembl_Genomes_Transcript", "asis")],
    "flybase_annotation": [("FlyBase", "asis"), ("Gene_Name", "cg_strip")],
    "wormbase_gene": [("WormBase", "asis")],
    "wormbase_protein": [("WormBase_Protein", "asis")],
    "wormbase_seqname": [("WormBase", "asis"), ("Gene_Name", "asis")],
    "uniprot": [("UniProtKB_AC-ID", "asis")],
    "uniprot_entry_name": [("UniProtKB_AC-ID", "asis")],
    "gene_name": [("Gene_Name", "asis")],
    "other": [("Gene_Name", "asis")],
}
# default taxon (UniProt Gene_Name lookups need one) when the user gives none
DEFAULT_TAXON = {"sgd_systematic": "559292", "arabidopsis_agi": "3702", "flybase_annotation": "7227",
                 "wormbase_seqname": "6239"}
# from-dbs that require a taxon
NEEDS_TAXON = {"Gene_Name"}


class MappingError(RuntimeError):
    pass


def _transform(i: str, how: str) -> str:
    if how == "unversioned":
        return re.sub(r"\.\d+$", "", i)
    if how == "upper":
        return i.upper()
    if how == "agi":
        return re.sub(r"\.\d+$", "", i.upper())
    if how == "cg_strip":
        return re.sub(r"-P[A-Z]$", "", i)
    return i


def guess_type(i: str) -> str:
    """Classify a bare identifier (used when the input has no id_type column)."""
    from .extract_ids import find_ids
    hs = find_ids(i)
    for h in hs:
        if h.span == (0, len(i)) or h.id == i:
            return h.id_type
    return "gene_name"


def run_job(from_db: str, to_db: str, ids: list[str], taxon: str | None = None, poll: float = 2.0,
            timeout: float = 1800) -> tuple[list[dict], str | None]:
    """Run one ID-mapping job; return (rows as dicts with 'From' + TSV_COLS fields, release)."""
    data = {"from": from_db, "to": to_db, "ids": ",".join(ids)}
    if taxon:
        data["taxId"] = str(taxon)
    r = post(f"{API}/idmapping/run", data=data)
    job = r.json()["jobId"]
    t0, delay = time.time(), poll
    while True:
        st = get(f"{API}/idmapping/status/{job}", ok_status=(200, 303), params=None, retries=5)
        if st.status_code == 303 or st.url != f"{API}/idmapping/status/{job}":
            break
        js = {}
        try:
            js = st.json()
        except ValueError:
            pass
        status = js.get("jobStatus")
        if status in ("FINISHED",) or "results" in js or "failedIds" in js:
            break
        if status in ("ERROR", "FAILED") or js.get("errors"):
            raise MappingError(f"UniProt job {job} failed: {js}")
        if time.time() - t0 > timeout:
            raise MappingError(f"UniProt job {job} timed out after {timeout}s")
        time.sleep(delay)
        delay = min(delay * 1.3, 15)
    if to_db == "UniProtKB":
        url = f"{API}/idmapping/uniprotkb/results/stream/{job}"
        res = get(url, params={"format": "tsv", "fields": FIELDS}, timeout=600)
    else:
        url = f"{API}/idmapping/stream/{job}"
        res = get(url, params={"format": "tsv"}, timeout=600)
    release = res.headers.get("X-UniProt-Release")
    rows = list(csv.DictReader(io.StringIO(res.text), delimiter="\t"))
    return rows, release


class Mapper:
    def __init__(self, taxon: str | None = None, batch: int = BATCH):
        self.taxon, self.batch = taxon, batch
        self.release: str | None = None
        self.jobs: list[dict] = []

    def _map_db(self, from_db: str, ids: list[str], taxon: str | None):
        out = []
        for i in range(0, len(ids), self.batch):
            chunk = ids[i:i + self.batch]
            try:
                rows, rel = run_job(from_db, "UniProtKB", chunk, taxon)
            except (HttpError, MappingError) as e:
                log.error("mapping job %s (%d ids) failed: %s", from_db, len(chunk), e)
                self.jobs.append({"from": from_db, "n_ids": len(chunk), "error": str(e)})
                continue
            self.release = rel or self.release
            self.jobs.append({"from": from_db, "n_ids": len(chunk), "n_rows": len(rows), "taxon": taxon})
            out.extend(rows)
        return out

    def map(self, items: list[tuple[str, str, str | None]]) -> list[dict]:
        """items: (original_id, id_type, taxon). Returns output rows (one or more per item)."""
        results: dict[tuple, list[dict]] = {}
        by_type: dict[str, list[tuple[str, str | None]]] = defaultdict(list)
        for oid, t, tax in items:
            by_type[t].append((oid, tax or self.taxon or DEFAULT_TAXON.get(t)))
        for t, lst in by_type.items():
            pending = list(OrderedDict.fromkeys(lst))  # unique (id, taxon)
            for from_db, how in CASCADE.get(t, CASCADE["other"]):
                if not pending:
                    break
                groups: dict[str | None, list[str]] = defaultdict(list)
                for oid, tax in pending:
                    if from_db in NEEDS_TAXON and not tax:
                        continue
                    groups[tax if from_db in NEEDS_TAXON else None].append(oid)
                still = []
                done_keys = set()
                for tax, oids in groups.items():
                    tmap: dict[str, str] = {}
                    for o in oids:
                        tmap.setdefault(_transform(o, how), o)
                    rows = self._map_db(from_db, list(tmap), tax)
                    for rw in rows:
                        o = tmap.get(rw["From"])
                        if o is None:
                            continue
                        acc = rw.get("Entry", "")
                        if from_db == "UniProtKB_AC-ID" and not rw.get("Length"):
                            continue  # inactive (deleted) entry: handled by the UniParc fallback below
                        method = f"{from_db}" + ("" if how == "asis" else f"[{how}]") + (f"+taxon={tax}" if tax else "")
                        note = "" if rw["From"] != acc else ""
                        if t in ("uniprot", "uniprot_entry_name") and rw["From"] not in (acc, rw.get("Entry Name")):
                            note = f"secondary_or_merged:{rw['From']}->{acc}"
                        results.setdefault((o, tax if from_db in NEEDS_TAXON else None, t), []).append(
                            dict(uniprot_acc=acc, taxid=rw.get("Organism (ID)", ""),
                                 gene_name=rw.get("Gene Names (primary)", ""),
                                 reviewed={"reviewed": "yes", "unreviewed": "no"}.get(rw.get("Reviewed", ""), rw.get("Reviewed", "")),
                                 length=rw.get("Length", ""), entry_name=rw.get("Entry Name", ""),
                                 mapping_method=method, note=note))
                        done_keys.add((o, tax if from_db in NEEDS_TAXON else None))
                done_ids = {k[0] for k in done_keys}
                pending = [(o, tx) for o, tx in pending if o not in done_ids]
        # UniParc fallback for uniprot-type ids that failed (obsolete/deleted accessions)
        unparc = {}
        missing = [(o, t) for o, t, _ in items if t == "uniprot" and not any(k[0] == o and k[2] == t for k in results)]
        if missing:
            try:
                rows, rel = run_job("UniProtKB_AC-ID", "UniParc", sorted({m[0] for m in missing}))
                for rw in rows:
                    unparc[rw["From"]] = rw.get("To") or rw.get("UniParc") or list(rw.values())[-1]
            except (HttpError, MappingError) as e:
                log.error("UniParc fallback failed: %s", e)
        out = []
        for oid, t, tax in items:
            tax = tax or self.taxon or DEFAULT_TAXON.get(t)
            hits = results.get((oid, tax, t)) or results.get((oid, None, t)) or []
            if hits:
                hits = sorted(hits, key=lambda h: (h["reviewed"] != "yes", -int(h["length"] or 0)))
                seen = set()
                hits = [h for h in hits if not (h["uniprot_acc"] in seen or seen.add(h["uniprot_acc"]))]
                for h in hits:
                    if tax and h["taxid"] and h["taxid"] != str(tax) and "taxon=" not in h["mapping_method"]:
                        h["note"] = (h["note"] + "; " if h["note"] else "") + f"taxon_differs_from_requested({tax})"
                    out.append(dict(original_id=oid, id_type=t, n_hits=len(hits), uniparc="", **h))
            elif oid in unparc:
                out.append(dict(original_id=oid, id_type=t, uniprot_acc="", taxid="", gene_name="", reviewed="",
                                length="", entry_name="", mapping_method="UniParc_fallback", n_hits=0,
                                note="accession obsolete/deleted in UniProtKB; sequence archived in UniParc",
                                uniparc=unparc[oid]))
            else:
                out.append(dict(original_id=oid, id_type=t, uniprot_acc="", taxid="", gene_name="", reviewed="",
                                length="", entry_name="", mapping_method="unmapped", n_hits=0,
                                note="no UniProtKB match via " + ",".join(f for f, _ in CASCADE.get(t, CASCADE["other"])),
                                uniparc=""))
        return out


OUT_COLS = ["original_id", "id_type", "uniprot_acc", "taxid", "gene_name", "reviewed", "length", "entry_name",
            "mapping_method", "n_hits", "note", "uniparc"]


def read_items(path: str, id_col="id", type_col="id_type", taxon_col="taxon") -> list[tuple[str, str, str | None]]:
    items, seen = [], set()
    with open(path, newline="") as f:
        head = f.readline()
        f.seek(0)
        if id_col in head.rstrip("\n").split("\t"):
            for r in csv.DictReader(f, delimiter="\t"):
                i = (r.get(id_col) or "").strip()
                t = (r.get(type_col) or "").strip() or guess_type(i)
                tax = (r.get(taxon_col) or "").strip() or None
                if i and (i, t, tax) not in seen:
                    seen.add((i, t, tax))
                    items.append((i, t, tax))
        else:
            for ln in f:
                i = ln.strip().split("\t")[0]
                if i and not i.startswith("#") and (i, "x", None) not in seen:
                    seen.add((i, "x", None))
                    items.append((i, guess_type(i), None))
    return items


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("input", nargs="?", help="TSV with an id column, or one id per line")
    ap.add_argument("--id", action="append", default=[], help="identifier (repeatable)")
    ap.add_argument("-o", "--out", help="output TSV (default stdout)")
    ap.add_argument("--taxon", help="NCBI taxon id: used for Gene_Name / SGD-style lookups (default for all ids)")
    ap.add_argument("--id-col", default="id")
    ap.add_argument("--type-col", default="id_type")
    ap.add_argument("--batch", type=int, default=BATCH)
    ap.add_argument("-v", "--verbose", action="store_true")
    a = ap.parse_args(argv)
    common.setup_logging(a.verbose)
    items = read_items(a.input, a.id_col, a.type_col) if a.input else []
    items += [(i, guess_type(i), None) for i in a.id]
    if not items:
        ap.error("no identifiers given")
    m = Mapper(a.taxon, a.batch)
    rows = m.map(items)
    out = open(a.out, "w") if a.out else sys.stdout
    w = csv.writer(out, delimiter="\t", lineterminator="\n")
    w.writerow(OUT_COLS)
    for r in rows:
        w.writerow([r.get(c, "") for c in OUT_COLS])
    if a.out:
        out.close()
        json.dump({"uniprot_release": m.release, "date": common.now_iso(), "jobs": m.jobs, "n_input": len(items),
                   "n_unmapped": sum(1 for r in rows if r["mapping_method"] == "unmapped")},
                  open(a.out + ".meta.json", "w"), indent=2)
    n_un = sum(1 for r in rows if r["mapping_method"] == "unmapped")
    log.info("%d ids -> %d rows; %d unmapped; UniProt release %s", len(items), len(rows), n_un, m.release)
    if n_un:
        log.warning("UNMAPPED: %s", ", ".join(r["original_id"] for r in rows if r["mapping_method"] == "unmapped")[:500])
    return 0


if __name__ == "__main__":
    sys.exit(main())
