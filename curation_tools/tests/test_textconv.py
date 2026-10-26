from curation_tools.textconv import jats_to_text, pmc_html_to_text, _grid_from_rows


def test_jats_structure(fix):
    c = jats_to_text((fix / "jats_sample.xml").read_bytes())
    t = c.text
    assert t.startswith("# Evolution of the TOY family")
    assert "Authors: Jane Smith" in t
    assert "## Introduction" in t and "### Subsection" in t       # nested heading levels
    assert "H_2O" in t and "10^3" in t                              # sub/sup preserved
    assert "a link <http://example.org/x>" in t
    assert "[FIGURE f1] Figure 1 Tree. Gene tree of the TOY family" in t
    assert c.n_figures == 1 and c.n_tables == 1
    assert "[SUPPLEMENT sd1]" in t and "media-1.xlsx" in t
    assert c.supp_hrefs == ["media-1.xlsx"]
    assert "[r1] Doe J, A paper" in t or "[r1] Doe J." in t         # reference rendered with spaces


def test_jats_table_tsv_spans(fix):
    c = jats_to_text((fix / "jats_sample.xml").read_bytes())
    rows = c.tables["t1"].split("\n")
    # colspan=2 header expands to an empty cell; rowspan values are repeated on each spanned row so that
    # every TSV row is self-contained
    assert rows[0].split("\t") == ["Species", "Genes", ""]
    assert rows[1].split("\t") == ["Yeast", "YAL001C", "YBR001C"]
    assert rows[2].split("\t") == ["Yeast", "YAL002W", "-"]
    assert "[TABLE t1] Table 1 Orthologs." in c.text and "[/TABLE t1]" in c.text
    assert "[table note] Genes are SGD systematic names." in c.text


def test_grid_rowspan_colspan():
    tsv = _grid_from_rows([[("a", 2, 1), ("b", 1, 2)], [("c", 1, 1), ("d", 1, 1)], [("e", 1, 1), ("f", 1, 1), ("g", 1, 1)]])
    assert tsv.split("\n")[0] == "a\t\tb"
    assert tsv.split("\n")[1] == "c\td\tb"


def test_jats_without_body_warns():
    c = jats_to_text(b"<article><front><article-meta><article-title>T</article-title></article-meta></front></article>")
    assert any("no <body>" in w for w in c.warnings)


def test_pmc_html(fix):
    c = pmc_html_to_text((fix / "pmc_sample.html").read_bytes())
    t = c.text
    assert t.startswith("# The TOY gene order browser")
    assert t.count("## Abstract") == 1                  # abstract not duplicated from the body
    assert "[FIGURE fig1] Figure 1. Phylogeny of yeasts." in t
    assert "Open in a new tab" not in t
    assert c.tables["tbl1"].split("\n") == ["Gene 1\tGene 2", "SPO21 (YOL091W)\tYSW1 (YBR148W)"]
    assert c.supp_hrefs == ["/articles/instance/1/bin/supp.xls"]
    assert "## References" in t and "[ref1] Doe J. 1999. A paper." in t
