"""Fetch a paper (DOI / PMID / PMCID) and its supplements, convert full text to plain text.

CLI:
  python -m curation_tools.fetch_paper resolve 10.1101/gr.3672305
  python -m curation_tools.fetch_paper fetch  16169922 [--out literature/papers] [--no-supp] [--force]
  python -m curation_tools.fetch_paper search "orthologous groups gene family phylogeny" -n 10 [--oa]

Retrieval cascade (every step is logged in meta.json -> "events"; nothing fails silently):
  full text : Europe PMC fullTextXML -> PMC OA open-data S3 -> legacy PMC OA service (oa.fcgi)
              -> NCBI efetch db=pmc (only if it has a body) -> PMC HTML article page
              -> PDF (PMC / open-access link) -> pdftotext/pypdf
  supplements: Europe PMC supplementaryFiles zip -> legacy OA tgz -> files linked from the PMC page
Layout:  <out>/<safe_id>/meta.json, fulltext.txt, tables/*.tsv, raw/{fulltext.xml,*.pdf,supp/...}
Files of articles without a CC/PD licence (redistributable=false) are written under raw/ only.
"""
from __future__ import annotations

import argparse
import io
import json
import re
import sys
import tarfile
import zipfile
from pathlib import Path
from urllib.parse import quote, urljoin

from lxml import etree

from . import common
from .common import HttpError, get, log, now_iso, sha256_file
from .textconv import Converted, jats_to_text, pmc_html_to_text

EPMC = "https://www.ebi.ac.uk/europepmc/webservices/rest"
S3 = "https://pmc-oa-opendata.s3.amazonaws.com"
OA_FCGI = ["https://www.ncbi.nlm.nih.gov/pmc/utils/oa/oa.fcgi", "https://pmc.ncbi.nlm.nih.gov/utils/oa/oa.fcgi"]
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
PMC_PAGE = "https://pmc.ncbi.nlm.nih.gov/articles"
DEFAULT_OUT = "literature/papers"
IMAGE_EXT = (".jpg", ".jpeg", ".gif", ".png", ".tif", ".tiff")
MIN_TEXT = 1500  # chars; below this a "full text" is considered missing


# ----------------------------------------------------------------------------- ids
def classify_id(s: str) -> tuple[str, str]:
    s = s.strip()
    s = re.sub(r"^(https?://(dx\.)?doi\.org/)", "", s, flags=re.I)
    m = re.match(r"^(pmid|pmcid|doi|pmc)\s*:\s*(.+)$", s, re.I)
    if m:
        k, v = m.group(1).lower(), m.group(2).strip()
        s = v if k != "pmc" else "PMC" + v
        if k in ("pmid", "pmcid", "doi"):
            return k, (v.upper() if k == "pmcid" else v)
    if re.fullmatch(r"PMC\d+", s, re.I):
        return "pmcid", s.upper()
    if re.fullmatch(r"\d{1,9}", s):
        return "pmid", s
    if s.startswith("10."):
        return "doi", s
    raise ValueError(f"cannot classify identifier {s!r} (expected DOI, PMID or PMCID)")


def _query_for(kind: str, val: str) -> str:
    if kind == "pmid":
        return f"EXT_ID:{val} AND SRC:MED"
    if kind == "pmcid":
        return f"PMCID:{val}"
    return f'DOI:"{val}"'


def parse_epmc_record(r: dict) -> dict:
    jr = (r.get("journalInfo") or {}).get("journal") or {}
    authors = [a.get("fullName") or f"{a.get('firstName', '')} {a.get('lastName', '')}".strip()
               for a in (r.get("authorList") or {}).get("author", [])]
    if not authors and r.get("authorString"):
        authors = [a.strip() for a in r["authorString"].rstrip(".").split(",")]
    lic = r.get("license")
    oa = r.get("isOpenAccess") == "Y"
    redistributable = bool(oa and lic and re.match(r"(cc|pd|public)", lic.strip(), re.I))
    return {
        "pmid": r.get("pmid"), "pmcid": r.get("pmcid"), "doi": r.get("doi"),
        "source": r.get("source"), "title": (r.get("title") or "").rstrip("."),
        "authors": authors, "year": int(r["pubYear"]) if r.get("pubYear", "").isdigit() else None,
        "journal": jr.get("title") or r.get("bookOrReportDetails", {}).get("publisher"),
        "is_open_access": oa, "license": lic, "redistributable": redistributable,
        "cited_by_count": r.get("citedByCount"), "has_supplements_flag": r.get("hasSuppl") == "Y",
        "in_epmc": r.get("inEPMC") == "Y", "in_pmc": r.get("inPMC") == "Y",
        "abstract": re.sub(r"<[^>]+>", "", r.get("abstractText") or "").strip(),
        "fulltext_urls": [{"style": u.get("documentStyle"), "availability": u.get("availability"),
                           "site": u.get("site"), "url": u.get("url")}
                          for u in (r.get("fullTextUrlList") or {}).get("fullTextUrl", [])],
    }


def resolve(identifier: str) -> dict:
    """Resolve DOI/PMID/PMCID to metadata through Europe PMC. Raises LookupError if not found."""
    kind, val = classify_id(identifier)
    r = get(f"{EPMC}/search", params={"query": _query_for(kind, val), "resultType": "core",
                                      "format": "json", "pageSize": 5})
    res = r.json().get("resultList", {}).get("result", [])
    if not res:
        raise LookupError(f"Europe PMC has no record for {kind}={val}")
    if len(res) > 1:  # e.g. a DOI shared by a preprint and a journal version: prefer MED/PMC
        res.sort(key=lambda x: {"MED": 0, "PMC": 1}.get(x.get("source"), 2))
    meta = parse_epmc_record(res[0])
    meta["query"] = {"kind": kind, "value": val}
    meta["europepmc_url"] = r.url
    return meta


def search(query: str, n: int = 10, oa_only: bool = False, sort: str | None = None,
           with_suppl: bool = False) -> list[dict]:
    """Top-n Europe PMC hits. sort: 'cited' (citation count desc), 'date', or None (relevance)."""
    q = query
    if oa_only:
        q += " AND OPEN_ACCESS:y"
    if with_suppl:
        q += " AND HAS_SUPPL:y"
    if sort == "cited":
        q += " sort_cited:y"
    elif sort == "date":
        q += " sort_date:y"
    rows, cursor = [], "*"
    while len(rows) < n:
        r = get(f"{EPMC}/search", params={"query": q, "resultType": "lite", "format": "json",
                                          "pageSize": min(100, n - len(rows)), "cursorMark": cursor})
        js = r.json()
        res = js.get("resultList", {}).get("result", [])
        for x in res:
            rows.append({"pmid": x.get("pmid") or "", "pmcid": x.get("pmcid") or "", "doi": x.get("doi") or "",
                         "title": (x.get("title") or "").rstrip("."), "year": x.get("pubYear") or "",
                         "journal": x.get("journalTitle") or "", "open_access": x.get("isOpenAccess") or "N",
                         "has_suppl": x.get("hasSuppl") or "N", "cited_by": x.get("citedByCount", 0),
                         "source": x.get("source")})
        nxt = js.get("nextCursorMark")
        if not res or not nxt or nxt == cursor:
            break
        cursor = nxt
    return rows[:n]


# ----------------------------------------------------------------------------- fetch
class Fetcher:
    def __init__(self, meta: dict, outdir: Path, force: bool = False):
        self.meta, self.dir, self.force = meta, outdir, force
        self.raw = outdir / "raw"
        self.events: list[dict] = []
        self.files: list[dict] = []
        self.pmcid = meta.get("pmcid")
        self.html_page: bytes | None = None
        self.html_checked = False

    # logging ---------------------------------------------------------
    def ev(self, step: str, ok: bool, url: str = "", detail: str = ""):
        self.events.append({"step": step, "ok": ok, "url": url, "detail": detail})
        (log.info if ok else log.warning)("%s: %s %s %s", step, "OK" if ok else "FAILED", url, detail)

    def save(self, rel: str, data: bytes, role: str, url: str = "") -> Path:
        p = self.dir / rel
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_bytes(data)
        self.files.append({"path": rel, "role": role, "bytes": len(data), "sha256": common.sha256_bytes(data),
                           "url": url})
        return p

    def txt_rel(self, name: str) -> str:
        return name if self.meta["redistributable"] else f"raw/{name}"

    # full text ---------------------------------------------------------
    def try_epmc_xml(self) -> Converted | None:
        url = f"{EPMC}/{self.pmcid}/fullTextXML"
        try:
            r = get(url, retries=2)
        except HttpError as e:
            self.ev("epmc_fulltext_xml", False, url, str(e))
            return None
        return self._take_xml("epmc_fulltext_xml", r.content, url)

    def _take_xml(self, step: str, content: bytes, url: str) -> Converted | None:
        if b"<errorBean" in content[:600] or b"errMsg" in content[:600]:
            self.ev(step, False, url, " ".join(re.sub(r"<[^>]+>", " ", content[:600].decode("utf8", "replace")).split())[:160])
            return None
        try:
            conv = jats_to_text(content)
        except Exception as e:  # noqa: BLE001
            self.ev(step, False, url, f"XML parse/convert error: {e}")
            return None
        if len(conv.text) < MIN_TEXT or "no <body>" in " ".join(conv.warnings):
            self.ev(step, False, url, f"no usable body ({len(conv.text)} chars); {conv.warnings}")
            return None
        self.save("raw/fulltext.xml", content, "fulltext_xml", url)
        self.ev(step, True, url, f"{len(conv.text)} chars, {conv.n_figures} figures, {conv.n_tables} tables")
        return conv

    def try_s3(self) -> Converted | None:
        url = f"{S3}/?list-type=2&prefix={self.pmcid}."
        try:
            r = get(url, retries=2)
            keys = re.findall(r"<Key>([^<]+)</Key>", r.text)
        except HttpError as e:
            self.ev("pmc_s3_opendata", False, url, str(e))
            return None
        xmls = sorted(k for k in keys if k.endswith(".xml"))
        if not xmls:
            self.ev("pmc_s3_opendata", False, url, "article not in PMC OA open-data bucket (not in OA subset)")
            return None
        k = xmls[-1]
        try:
            x = get(f"{S3}/{k}")
        except HttpError as e:
            self.ev("pmc_s3_opendata", False, f"{S3}/{k}", str(e))
            return None
        return self._take_xml("pmc_s3_opendata", x.content, f"{S3}/{k}")

    def try_oa_fcgi(self) -> tuple[Converted | None, list[str]]:
        """Legacy PMC OA web service. Returns (converted_text, supplement files extracted)."""
        for base in OA_FCGI:
            url = f"{base}?id={self.pmcid}"
            try:
                r = get(url, retries=1)
            except HttpError as e:
                self.ev("pmc_oa_fcgi", False, url, str(e))
                continue
            try:
                root = etree.fromstring(r.content)
            except etree.XMLSyntaxError:
                self.ev("pmc_oa_fcgi", False, url, "non-XML response")
                continue
            err = root.find(".//error")
            if err is not None:
                self.ev("pmc_oa_fcgi", False, url, (err.get("code") or "") + " " + (err.text or ""))
                continue
            tgz = next((l.get("href") for l in root.iter("link") if l.get("format") == "tgz"), None)
            if not tgz:
                self.ev("pmc_oa_fcgi", False, url, "no tgz link")
                continue
            tgz = re.sub(r"^ftp://", "https://", tgz)
            try:
                blob = get(tgz, timeout=600).content
            except HttpError as e:
                self.ev("pmc_oa_fcgi", False, tgz, str(e))
                continue
            conv, supp = self._unpack_tgz(blob, tgz)
            self.ev("pmc_oa_fcgi", conv is not None, tgz, f"{len(supp)} supplementary files")
            return conv, supp
        return None, []

    def _unpack_tgz(self, blob: bytes, url: str):
        conv, supp = None, []
        with tarfile.open(fileobj=io.BytesIO(blob), mode="r:gz") as tf:
            for m in tf.getmembers():
                if not m.isfile():
                    continue
                data = tf.extractfile(m).read()
                name = Path(m.name).name
                if name.endswith((".nxml", ".xml")) and conv is None:
                    conv = self._take_xml("pmc_oa_fcgi_xml", data, url)
                else:
                    self.save(f"raw/supp/{common.safe_name(name)}", data, "supplement", url)
                    supp.append(name)
        return conv, supp

    def try_efetch(self) -> Converted | None:
        url = f"{EUTILS}/efetch.fcgi?db=pmc&id={self.pmcid}&rettype=xml"
        try:
            r = get(url, retries=2)
        except HttpError as e:
            self.ev("ncbi_efetch_pmc", False, url, str(e))
            return None
        return self._take_xml("ncbi_efetch_pmc", r.content, url)

    def pmc_page(self) -> bytes | None:
        if self.html_checked:
            return self.html_page
        self.html_checked = True
        url = f"{PMC_PAGE}/{self.pmcid}/"
        try:
            self.html_page = get(url, retries=2).content
        except HttpError as e:
            self.ev("pmc_html_page", False, url, str(e))
        return self.html_page

    def try_html(self) -> Converted | None:
        page = self.pmc_page()
        url = f"{PMC_PAGE}/{self.pmcid}/"
        if page is None:
            return None
        conv = pmc_html_to_text(page)
        if len(conv.text) < MIN_TEXT or conv.warnings and "no article body" in " ".join(conv.warnings):
            self.ev("pmc_html_page", False, url, f"no usable body ({len(conv.text)} chars)")
            self.save("raw/article.html", page, "html", url)
            return None
        self.save("raw/article.html", page, "html", url)
        self.ev("pmc_html_page", True, url, f"{len(conv.text)} chars (free-to-read page; not necessarily re-usable)")
        return conv

    def try_pdf(self) -> str | None:
        cands = []
        page = self.pmc_page() if self.pmcid else None
        if page:
            from lxml import html as lhtml
            root = lhtml.fromstring(page)
            for h in root.xpath("//a[contains(@href,'.pdf')]/@href"):
                if "/bin/" in h:       # supplementary PDFs are not the paper
                    continue
                cands.append(urljoin(f"{PMC_PAGE}/{self.pmcid}/", h))
        cands += [u["url"] for u in self.meta.get("fulltext_urls", [])
                  if u.get("style") == "pdf" and u.get("availability") in ("Open access", "Free")]
        seen = set()
        for u in cands:
            if u in seen:
                continue
            seen.add(u)
            try:
                r = get(u, retries=1, timeout=180)
            except HttpError as e:
                self.ev("pdf_download", False, u, str(e))
                continue
            if not r.content.startswith(b"%PDF"):
                self.ev("pdf_download", False, u, f"not a PDF (content-type {r.headers.get('content-type')})")
                continue
            p = self.save("raw/fulltext.pdf", r.content, "fulltext_pdf", u)
            try:
                text, method = common.pdf_to_text(p)
            except Exception as e:  # noqa: BLE001
                self.ev("pdf_to_text", False, str(p), str(e))
                return None
            if len(text) < MIN_TEXT:
                self.ev("pdf_to_text", False, str(p), f"only {len(text)} chars (scanned PDF? needs OCR)")
                return None
            self.ev("pdf_to_text", True, u, f"{method}, {len(text)} chars")
            return text
        return None

    def get_fulltext(self):
        m = self.meta
        conv, how, text = None, None, None
        if self.pmcid:
            for name, fn in (("europepmc_xml", self.try_epmc_xml), ("pmc_s3_xml", self.try_s3)):
                conv = fn()
                if conv:
                    how = name
                    break
            if not conv:
                conv, supp = self.try_oa_fcgi()
                self._oa_supp = supp
                how = "pmc_oa_fcgi" if conv else None
            if not conv:
                conv = self.try_efetch()
                how = "ncbi_efetch" if conv else None
            if not conv:
                conv = self.try_html()
                how = "pmc_html" if conv else None
        else:
            self.ev("pmcid", False, "", "no PMCID: article is not in PMC/Europe PMC full-text collection")
        if not conv:
            text = self.try_pdf()
            how = "pdf" if text else None
        if conv:
            text = conv.text
            for tid, tsv in conv.tables.items():
                self.save(self.txt_rel(f"tables/{common.safe_name(tid)}.tsv"), (tsv + "\n").encode(), "table_tsv")
        if text:
            self.save(self.txt_rel("fulltext.txt"), text.encode(), "fulltext_text")
            m["fulltext"] = {"source": how, "text_file": self.txt_rel("fulltext.txt"), "chars": len(text),
                             "n_figures": conv.n_figures if conv else None,
                             "n_tables": conv.n_tables if conv else None,
                             "warnings": conv.warnings if conv else []}
            self.supp_hrefs = conv.supp_hrefs if conv else []
        else:
            m["fulltext"] = None
            self.supp_hrefs = []
            others = [u for u in m.get("fulltext_urls", []) if u.get("url")]
            m.setdefault("not_retrievable", []).append(
                "full text: all automatic routes failed; publisher/other URLs for manual retrieval: "
                + "; ".join(u["url"] for u in others[:5]))

    # supplements ---------------------------------------------------------
    def _save_supp(self, name: str, data: bytes, url: str) -> str:
        rel = f"raw/supp/{common.safe_name(name) or 'file'}"
        self.save(rel, data, "supplement", url)
        return rel

    def get_supplements(self):
        m = self.meta
        got: list[str] = list(getattr(self, "_oa_supp", []) or [])
        if self.pmcid:
            url = f"{EPMC}/{self.pmcid}/supplementaryFiles"
            try:
                r = get(url, retries=1, timeout=600)
                if r.content[:2] == b"PK":
                    self.save("raw/supplementary.zip", r.content, "supplement_zip", url)
                    n_img = 0
                    with zipfile.ZipFile(io.BytesIO(r.content)) as z:
                        for info in z.infolist():
                            if info.is_dir():
                                continue
                            if info.filename.lower().endswith(IMAGE_EXT):
                                n_img += 1      # figure images stay inside supplementary.zip only
                                continue
                            rel = self._save_supp(info.filename.replace("/", "__"), z.read(info), url)
                            got.append(rel)
                    self.ev("epmc_supplementary_zip", True, url,
                            f"{len(got)} files extracted to raw/supp/; {n_img} image files left inside raw/supplementary.zip")
                else:
                    msg = re.sub(r"<[^>]+>", " ", r.text[:500])
                    self.ev("epmc_supplementary_zip", False, url, " ".join(msg.split())[:200])
            except HttpError as e:
                self.ev("epmc_supplementary_zip", False, url, str(e))
            except zipfile.BadZipFile as e:
                self.ev("epmc_supplementary_zip", False, url, f"bad zip: {e}")
        if not got and self.pmcid:
            page = self.pmc_page()
            if page:
                from lxml import html as lhtml
                root = lhtml.fromstring(page)
                base = f"{PMC_PAGE}/{self.pmcid}/"
                links = []
                for h in root.xpath("//section[contains(@class,'supplementary-materials')]//a/@href"):
                    u = urljoin("https://pmc.ncbi.nlm.nih.gov/", h)
                    if "/bin/" in u and not u.endswith("_index.html") and u not in links:
                        links.append(u)
                if not links:
                    self.ev("pmc_page_supplements", False, base, "no supplementary links on PMC page")
                for u in links:
                    try:
                        r = get(u, retries=2, timeout=600)
                        if r.headers.get("content-type", "").startswith("text/html") and b"<html" in r.content[:2000].lower() \
                                and not u.lower().endswith((".html", ".htm")):
                            pow_ = b"POW_CHALLENGE" in r.content
                            self.ev("pmc_page_supplements", False, u,
                                    "PMC served an interactive proof-of-work page instead of the file; automated "
                                    "bypass deliberately NOT attempted - download manually (requested/)" if pow_
                                    else "got HTML instead of the file")
                            if pow_:
                                break
                            continue
                        got.append(self._save_supp(u.rsplit("/", 1)[-1], r.content, u))
                        self.ev("pmc_page_supplements", True, u, f"{len(r.content)} bytes")
                    except HttpError as e:
                        self.ev("pmc_page_supplements", False, u, str(e))
        m["supplements"] = got
        if not got:
            if m.get("has_supplements_flag") or getattr(self, "supp_hrefs", []):
                m.setdefault("not_retrievable", []).append(
                    "supplementary files: Europe PMC flags supplements / XML lists "
                    f"{getattr(self, 'supp_hrefs', [])} but none could be downloaded; get them from the publisher site")
            else:
                self.ev("supplements", True, "", "article has no supplementary files listed")


def fetch(identifier: str, out_root: str | Path = DEFAULT_OUT, supp: bool = True, force: bool = False) -> dict:
    meta = resolve(identifier)
    pid = (f"PMID{meta['pmid']}" if meta.get("pmid") else meta.get("pmcid")
           or "doi_" + common.safe_name(meta.get("doi") or identifier))
    outdir = Path(out_root) / pid
    outdir.mkdir(parents=True, exist_ok=True)
    f = Fetcher(meta, outdir, force)
    meta["safe_id"] = pid
    meta["retrieval_date"] = now_iso()
    meta["toolkit_version"] = __import__("curation_tools").__version__
    meta["not_retrievable"] = []
    if not meta["is_open_access"]:
        f.ev("license", True, "", "not in the Europe PMC open-access subset: derived text stored under raw/ "
                                  "(not for redistribution); check publisher terms")
    f.get_fulltext()
    if supp:
        f.get_supplements()
    for x in f.files:  # ensure hashes are of the final file content
        x["sha256"] = sha256_file(outdir / x["path"])
    meta["files"] = f.files
    meta["events"] = f.events
    (outdir / "meta.json").write_text(json.dumps(meta, indent=2, ensure_ascii=False) + "\n")
    return meta


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-v", "--verbose", action="store_true")
    sub = ap.add_subparsers(dest="cmd", required=True)
    a = sub.add_parser("resolve", help="print metadata for DOI/PMID/PMCID")
    a.add_argument("id")
    a = sub.add_parser("fetch", help="download full text + supplements")
    a.add_argument("ids", nargs="+")
    a.add_argument("--out", default=DEFAULT_OUT)
    a.add_argument("--no-supp", action="store_true")
    a.add_argument("--force", action="store_true")
    a = sub.add_parser("search", help="search Europe PMC, print TSV")
    a.add_argument("query")
    a.add_argument("-n", type=int, default=10)
    a.add_argument("--oa", action="store_true", help="open access only")
    a.add_argument("--suppl", action="store_true", help="only papers with supplementary files")
    a.add_argument("--sort", choices=["cited", "date"], default=None)
    args = ap.parse_args(argv)
    common.setup_logging(args.verbose)
    if args.cmd == "resolve":
        meta = resolve(args.id)
        meta.pop("abstract", None)
        print(json.dumps(meta, indent=2))
        return 0
    if args.cmd == "search":
        rows = search(args.query, args.n, args.oa, args.sort, args.suppl)
        cols = ["pmid", "pmcid", "doi", "title", "year", "journal", "open_access", "has_suppl", "cited_by"]
        print("\t".join(cols))
        for r in rows:
            print("\t".join(str(r[c]).replace("\t", " ") for c in cols))
        if not rows:
            log.warning("no hits for %r", args.query)
            return 2
        return 0
    rc = 0
    for i in args.ids:
        try:
            m = fetch(i, args.out, supp=not args.no_supp, force=args.force)
        except (LookupError, ValueError, HttpError) as e:
            log.error("%s: cannot resolve/fetch: %s", i, e)
            rc = max(rc, 1)
            continue
        ft = m["fulltext"]
        print(f"{m['safe_id']}\t{m['title'][:70]}\tOA={m['is_open_access']}\tlicense={m['license']}\t"
              f"fulltext={ft['source'] if ft else 'NONE'}\tsupplements={len(m['supplements'])}")
        for msg in m["not_retrievable"]:
            log.error("NOT RETRIEVABLE [%s]: %s", m["safe_id"], msg)
            rc = max(rc, 2)
    return rc


if __name__ == "__main__":
    sys.exit(main())
