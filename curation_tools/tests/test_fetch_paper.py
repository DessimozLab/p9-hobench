import json

import pytest

from curation_tools import fetch_paper as fp
from curation_tools.common import HttpError


def test_classify_id():
    assert fp.classify_id("16169922") == ("pmid", "16169922")
    assert fp.classify_id("PMC1240090") == ("pmcid", "PMC1240090")
    assert fp.classify_id("pmc1240090") == ("pmcid", "PMC1240090")
    assert fp.classify_id("10.1101/gr.3672305") == ("doi", "10.1101/gr.3672305")
    assert fp.classify_id("https://doi.org/10.1101/gr.3672305") == ("doi", "10.1101/gr.3672305")
    assert fp.classify_id("doi:10.1/x") == ("doi", "10.1/x")
    assert fp.classify_id("PMID: 123") == ("pmid", "123")
    with pytest.raises(ValueError):
        fp.classify_id("not-an-id")


def test_parse_epmc_record(fix):
    r = json.loads((fix / "epmc_core.json").read_text())["resultList"]["result"][0]
    m = fp.parse_epmc_record(r)
    assert m["pmid"] == "23193287" and m["pmcid"] == "PMC3531190" and m["year"] == 2013
    assert m["is_open_access"] and m["license"] == "cc by-nc" and m["redistributable"]
    assert m["authors"] == ["Benson DA", "Cavanaugh M"]
    assert m["abstract"].endswith("GenBank is a database.")
    assert m["fulltext_urls"][0]["style"] == "pdf"
    r2 = dict(r, isOpenAccess="N", license=None)
    assert not fp.parse_epmc_record(r2)["redistributable"]


class FakeResp:
    def __init__(self, content=b"", status=200, headers=None, url=""):
        self.content, self.status_code, self.headers, self.url = content, status, headers or {}, url
        self.text = content.decode("utf8", "replace")

    def json(self):
        return json.loads(self.text)


def test_fetch_cascade_offline(tmp_path, fix, monkeypatch):
    """Non-OA article: EPMC XML 500 -> S3 absent -> oa.fcgi 404 -> efetch has no body -> PMC HTML works;
    supplements blocked by PoW page. Everything must be logged; nothing silent."""
    core = (fix / "epmc_core.json").read_text()
    core = json.loads(core)
    rec = core["resultList"]["result"][0]
    rec.update(isOpenAccess="N", pmcid="PMC1", pmid="1")
    rec.pop("license")
    html = (fix / "pmc_sample.html").read_bytes()
    big = html.replace(b"<p>Results paragraph.</p>", b"<p>" + b"Long results text. " * 200 + b"</p>")

    def fake_get(url, **kw):
        if "europepmc/webservices/rest/search" in url:
            return FakeResp(json.dumps(core).encode(), url=url)
        if url.endswith("/fullTextXML"):
            raise HttpError(url, 500, "server error")
        if "pmc-oa-opendata" in url:
            return FakeResp(b"<ListBucketResult></ListBucketResult>")
        if "oa.fcgi" in url:
            raise HttpError(url, 404)
        if "efetch" in url:
            return FakeResp(b"<pmc-articleset><article><front><article-meta><article-title>T</article-title></article-meta></front></article></pmc-articleset>")
        if url.endswith("/PMC1/"):
            return FakeResp(big)
        if url.endswith("supplementaryFiles"):
            return FakeResp(b"<errorBean><errMsg>not open access</errMsg></errorBean>")
        if "/bin/" in url:
            return FakeResp(b"<html>POW_CHALLENGE</html>", headers={"content-type": "text/html"})
        raise AssertionError("unexpected url " + url)

    monkeypatch.setattr(fp, "get", fake_get)
    meta = fp.fetch("PMC1", tmp_path)
    d = tmp_path / "PMID1"
    assert meta["fulltext"]["source"] == "pmc_html"
    assert (d / "raw" / "fulltext.txt").exists() and not (d / "fulltext.txt").exists()   # non-redistributable -> raw/
    assert (d / "raw" / "tables" / "tbl1.tsv").exists()
    steps = {(e["step"], e["ok"]) for e in meta["events"]}
    assert ("epmc_fulltext_xml", False) in steps and ("pmc_s3_opendata", False) in steps
    assert ("pmc_oa_fcgi", False) in steps and ("ncbi_efetch_pmc", False) in steps and ("pmc_html_page", True) in steps
    assert any("proof-of-work" in e["detail"] for e in meta["events"])
    assert any("supplementary" in m for m in meta["not_retrievable"]) or meta["supplements"] == []
    saved = json.loads((d / "meta.json").read_text())
    assert all(len(f["sha256"]) == 64 for f in saved["files"])
    assert saved["redistributable"] is False and saved["retrieval_date"]


def test_fetch_no_pmcid_is_reported(tmp_path, fix, monkeypatch):
    core = json.loads((fix / "epmc_core.json").read_text())
    rec = core["resultList"]["result"][0]
    rec.update(pmcid=None, isOpenAccess="N", fullTextUrlList={"fullTextUrl": []})
    monkeypatch.setattr(fp, "get", lambda url, **kw: FakeResp(json.dumps(core).encode(), url=url))
    meta = fp.fetch("23193287", tmp_path)
    assert meta["fulltext"] is None
    assert any("full text" in m for m in meta["not_retrievable"])


@pytest.mark.network
def test_resolve_and_search_live():
    m = fp.resolve("16169922")
    assert m["pmcid"] == "PMC1240090" and m["doi"] == "10.1101/gr.3672305" and m["year"] == 2005
    assert fp.resolve("10.1101/gr.3672305")["pmid"] == "16169922"
    assert fp.resolve("PMC1240090")["pmid"] == "16169922"
    hits = fp.search("orthology benchmark quest for orthologs", n=5)
    assert 1 <= len(hits) <= 5 and {"pmid", "doi", "title", "cited_by", "open_access"} <= set(hits[0])


@pytest.mark.network
def test_fetch_oa_live(tmp_path):
    m = fp.fetch("PMC3531190", tmp_path, supp=False)   # GenBank NAR 2013, OA
    assert m["fulltext"]["source"] in ("europepmc_xml", "pmc_s3_xml")
    assert (tmp_path / "PMID23193287" / "fulltext.txt").stat().st_size > 10000
