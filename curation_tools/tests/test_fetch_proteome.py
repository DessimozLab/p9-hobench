import gzip
import json

import pytest

from curation_tools import fetch_proteome as fpm
from curation_tools.common import read_fasta

SP = (">sp|P00498|HIS1_YEAST ATP phosphoribosyltransferase OS=Saccharomyces cerevisiae (strain ATCC 204508 / S288c) "
      "OX=559292 GN=HIS1 PE=1 SV=1\nMDLVNHLTDRLLFAIPKKGRLYSKS\nVSILNG\n")
TR = ">tr|Q12345|Q12345_YEAST Putative protein OS=Saccharomyces cerevisiae OX=559292 PE=4 SV=2\nMKV\n"


def test_parse_header():
    p = fpm.parse_header(SP.splitlines()[0][1:])
    assert (p["db"], p["acc"], p["entry"], p["gn"], p["pe"], p["sv"], p["ox"]) == \
           ("sp", "P00498", "HIS1_YEAST", "HIS1", "1", "1", "559292")
    assert p["name"] == "ATP phosphoribosyltransferase"
    q = fpm.parse_header(TR.splitlines()[0][1:])
    assert q["gn"] == "" and q["acc"] == "Q12345" and q["pe"] == "4"
    odd = fpm.parse_header("weird header without structure")
    assert odd["acc"] == "weird"


def test_rewrite_headers_and_map(tmp_path):
    src = tmp_path / "in.fasta"
    src.write_text(SP + TR)
    n, nres = fpm.rewrite(src, tmp_path / "o.fa", tmp_path / "m.tsv", {"P00498": ("SGD:S000000857", "SGD:S000000857")}, "canonical")
    assert (n, nres) == (2, 25 + 6 + 3)
    recs = list(read_fasta(tmp_path / "o.fa"))
    assert [h for h, _ in recs] == ["P00498", "Q12345"] and recs[0][1] == "MDLVNHLTDRLLFAIPKKGRLYSKSVSILNG"
    rows = [l.split("\t") for l in (tmp_path / "m.tsv").read_text().splitlines()]
    assert rows[0][:2] == ["accession", "original_header"] and rows[1][1].startswith("sp|P00498|HIS1_YEAST")
    assert rows[1][12] == "SGD:S000000857" and rows[1][10] == "31"


def test_duplicate_accession_rejected(tmp_path):
    src = tmp_path / "in.fasta"
    src.write_text(SP + SP)
    with pytest.raises(ValueError):
        fpm.rewrite(src, tmp_path / "o.fa", tmp_path / "m.tsv", None, "canonical")


def test_load_gene2acc(tmp_path):
    p = tmp_path / "g.gz"
    with gzip.open(p, "wt") as f:
        f.write("SGD:S1\tP1\tSGD:S1\nSGD:S2\tP2\tSGD:S2\nSGD:S2\tP3\tSGD:S2\n")
    d = fpm.load_gene2acc(p)
    assert d["P3"] == ("SGD:S2", "SGD:S2") and len(d) == 3


@pytest.mark.network
def test_live_yeast_reference_proteome(tmp_path):
    m = fpm.fetch("UP000002311", tmp_path)
    assert m["source"] == "ftp" and m["is_reference_or_representative"] and m["taxid"] == 559292
    assert m["n_proteins"] == m["gene_count_reported"] == m["n_canonical"]
    assert m["uniprot_release"] and len(m["outputs"]["proteome_fasta"]["sha256"]) == 64
    d = tmp_path / "UP000002311_559292"
    assert (d / "proteome.fasta").exists() and (d / "header_map.tsv").exists()
    ids = [h for h, _ in read_fasta(d / "proteome.fasta")]
    assert "P34111" in ids and len(set(ids)) == len(ids) and all("|" not in i for i in ids)
    assert json.loads((d / "manifest.json").read_text())["proteome_id"] == "UP000002311"
