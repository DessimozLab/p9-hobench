import io
import zipfile

import pytest

from curation_tools import extract_ids as ex


def ids(text):
    return {(h.id, h.id_type) for h in ex.find_ids(text)}


def test_patterns_positive():
    got = ids("NP_001234.1 XP_123456789.2 YP_009724390.1 WP_000000001.1 ENSDARP00000012345 ENSP00000269305.4 "
              "ENSMUSG00000000001 FBpp0070001 FBgn0000008 WBGene00000001 YAL001C YDR001W-A At1g01010.1 AT5G10000 "
              "P12345 Q9H0H5-2 A0A024R161 AAA12345.1 CAA1234567.1 RAD51_HUMAN NM_000546.6")
    for x in [("NP_001234.1", "refseq_protein"), ("XP_123456789.2", "refseq_protein"), ("YP_009724390.1", "refseq_protein"),
              ("WP_000000001.1", "refseq_protein"), ("ENSDARP00000012345", "ensembl_protein"),
              ("ENSP00000269305.4", "ensembl_protein"), ("ENSMUSG00000000001", "ensembl_gene"),
              ("FBpp0070001", "flybase_protein"), ("FBgn0000008", "flybase_gene"), ("WBGene00000001", "wormbase_gene"),
              ("YAL001C", "sgd_systematic"), ("YDR001W-A", "sgd_systematic"), ("AT1G01010.1", "arabidopsis_agi"),
              ("AT5G10000", "arabidopsis_agi"), ("P12345", "uniprot"), ("Q9H0H5-2", "uniprot"),
              ("A0A024R161", "uniprot"), ("AAA12345.1", "genbank_protein"), ("CAA1234567.1", "genbank_protein"),
              ("RAD51_HUMAN", "uniprot_entry_name"), ("NM_000546.6", "refseq_nucleotide")]:
        assert x in got, x


def test_patterns_other_family():
    got = ids("jgi|Aspni7|12345 Potri.001G123400.1 LOC_Os01g01010.1 TGME49_123456 PF3D7_0100100 DDB_G0267170 AGAP001234 SPAC1002.01")
    assert {t for _, t in got} == {"other"}
    assert len(got) == 8


def test_patterns_negative():
    # non-ID look-alikes must not be reported
    got = ids("PMC1240090 SRR1234567 GSE12345 Figure 1.2 Table S3 p<0.05 10.1101/gr.3672305 ABCD12 H2B A1B2 2005 NP_12")
    assert got == set()


def test_specific_beats_generic():
    # an Ensembl/AGI/SGD id must not additionally be reported as UniProt/GenBank-shaped
    got = ids("AT1G01010 YAL001C ENSDARP00000012345")
    assert {t for _, t in got} == {"arabidopsis_agi", "sgd_systematic", "ensembl_protein"}


def test_scan_text_and_tsv_context(tmp_path):
    (tmp_path / "a.txt").write_text("line one\nOrtholog is NP_014368.1 here\nline three\n")
    (tmp_path / "b.tsv").write_text("gene\tspecies\tnote\nYAL001C\tyeast\tTFC3\nENSDARP00000012345\tzebrafish\tx\n")
    rows, err = ex.scan([tmp_path])
    assert not err
    by = {r["id"]: r for r in rows}
    assert by["NP_014368.1"]["location"] == "line 2"
    assert "line one" in by["NP_014368.1"]["context"] and "line three" in by["NP_014368.1"]["context"]
    assert by["YAL001C"]["location"].endswith("row 2") and by["YAL001C"]["col_header"] == "gene"
    assert by["YAL001C"]["column"] == "A" and "yeast" in by["YAL001C"]["context"] and "TFC3" in by["YAL001C"]["context"]
    assert by["ENSDARP00000012345"]["column"] == "A"


def test_scan_xlsx_docx_pdf_zip(tmp_path):
    openpyxl = pytest.importorskip("openpyxl")
    docx = pytest.importorskip("docx")
    wb = openpyxl.Workbook()
    ws = wb.active
    ws.title = "Orthologs"
    ws.append(["Gene", "Protein"])
    ws.append(["TFC3", "NP_014368.1; XP_001234567.1"])
    wb.save(tmp_path / "t.xlsx")
    d = docx.Document()
    d.add_paragraph("Paragraph mentioning FBgn0000008.")
    t = d.add_table(rows=2, cols=2)
    t.rows[0].cells[0].text, t.rows[0].cells[1].text = "sp", "id"
    t.rows[1].cells[0].text, t.rows[1].cells[1].text = "worm", "WBGene00000001"
    d.save(tmp_path / "t.docx")
    (tmp_path / "t.pdf").write_bytes(_mini_pdf("Fly protein FBpp0070001 on page one"))
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w") as z:
        z.writestr("inner/list.txt", "AT1G01010.1\n")
        z.writestr("img.png", b"\x89PNG")
    (tmp_path / "s.zip").write_bytes(buf.getvalue())
    rows, err = ex.scan([tmp_path])
    assert not err
    s = {(r["id"], r["file"].replace(str(tmp_path) + "/", "")): r for r in rows}
    assert s[("NP_014368.1", "t.xlsx")]["location"] == "sheet 'Orthologs' row 2"
    assert s[("XP_001234567.1", "t.xlsx")]["col_header"] == "Protein"
    assert ("FBgn0000008", "t.docx") in s
    assert s[("WBGene00000001", "t.docx")]["location"] == "table 1 row 2"
    assert ("FBpp0070001", "t.pdf") in s
    assert ("AT1G01010.1", "s.zip!inner/list.txt") in s


def test_unique_and_filters(tmp_path):
    (tmp_path / "a.txt").write_text("P12345 and NP_014368.1\nNP_014368.1 again\n")
    rows, _ = ex.scan([tmp_path])
    u = ex.unique(rows)
    assert {r["id"]: r["n_occurrences"] for r in u}["NP_014368.1"] == 2
    rows_hi, _ = ex.scan([tmp_path], min_conf="high")
    assert {r["id"] for r in rows_hi} == {"NP_014368.1"}
    rows_t, _ = ex.scan([tmp_path], types={"uniprot"})
    assert {r["id"] for r in rows_t} == {"P12345"}


def test_unreadable_file_is_reported(tmp_path):
    (tmp_path / "bad.xlsx").write_bytes(b"not an xlsx")
    rows, err = ex.scan([tmp_path])
    assert rows == [] and len(err) == 1 and "bad.xlsx" in err[0]


def test_fasta_only_headers(tmp_path):
    (tmp_path / "x.fasta").write_text(">NP_014368.1 histone\nMSGRGKGGKGLGKGGAKRHRKILRDNIQGITKPAIRRLARRGGVKRISGLIYEEVRAVLKSFLESVIRDS\n")
    rows, _ = ex.scan([tmp_path])
    assert [r["id"] for r in rows] == ["NP_014368.1"]


def _mini_pdf(text: str) -> bytes:
    objs = [b"<< /Type /Catalog /Pages 2 0 R >>", b"<< /Type /Pages /Kids [3 0 R] /Count 1 >>",
            b"<< /Type /Page /Parent 2 0 R /MediaBox [0 0 612 792] /Contents 4 0 R /Resources << /Font << /F1 5 0 R >> >> >>"]
    stream = f"BT /F1 12 Tf 72 700 Td ({text}) Tj ET".encode()
    objs.append(b"<< /Length %d >>\nstream\n" % len(stream) + stream + b"\nendstream")
    objs.append(b"<< /Type /Font /Subtype /Type1 /BaseFont /Helvetica >>")
    out = b"%PDF-1.4\n"
    offs = []
    for i, o in enumerate(objs, 1):
        offs.append(len(out))
        out += b"%d 0 obj\n" % i + o + b"\nendobj\n"
    x = len(out)
    out += b"xref\n0 %d\n0000000000 65535 f \n" % (len(objs) + 1)
    for o in offs:
        out += b"%010d 00000 n \n" % o
    out += b"trailer\n<< /Size %d /Root 1 0 R >>\nstartxref\n%d\n%%%%EOF\n" % (len(objs) + 1, x)
    return out
