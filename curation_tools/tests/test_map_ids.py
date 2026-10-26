import pytest

from curation_tools import map_ids as mi


def fake_rows(mapping, to="UniProtKB"):
    """Build a run_job replacement from {(from_db, id): [(acc, name, reviewed, taxid, length)]}."""
    calls = []

    def run_job(from_db, to_db, ids, taxon=None, **kw):
        calls.append((from_db, to_db, list(ids), taxon))
        rows = []
        for i in ids:
            if to_db == "UniParc":
                if (from_db, i) in mapping:
                    rows.append({"From": i, "To": mapping[(from_db, i)]})
                continue
            for acc, gene, rev, tax, ln in mapping.get((from_db, i), []):
                rows.append({"From": i, "Entry": acc, "Entry Name": acc + "_X", "Reviewed": rev,
                             "Gene Names (primary)": gene, "Organism (ID)": tax, "Length": ln})
        return rows, "2099_01"

    return run_job, calls


def test_cascade_unversioned_and_unmapped(monkeypatch):
    run, calls = fake_rows({
        ("RefSeq_Protein", "NP_014368"): [("P02309", "HHF1", "reviewed", "559292", "103")],
        ("Ensembl_Protein", "ENSMUSP00000000001"): [("Q9DC51", "Gnai3", "reviewed", "10090", "354")],
    })
    monkeypatch.setattr(mi, "run_job", run)
    m = mi.Mapper()
    rows = m.map([("NP_014368.1", "refseq_protein", None), ("ENSMUSP00000000001.5", "ensembl_protein", None),
                  ("NP_999999.1", "refseq_protein", None)])
    by = {r["original_id"]: r for r in rows}
    assert by["NP_014368.1"]["uniprot_acc"] == "P02309" and by["NP_014368.1"]["mapping_method"] == "RefSeq_Protein[unversioned]"
    assert by["NP_014368.1"]["reviewed"] == "yes" and by["NP_014368.1"]["length"] == "103"
    assert by["ENSMUSP00000000001.5"]["uniprot_acc"] == "Q9DC51"
    assert by["NP_999999.1"]["mapping_method"] == "unmapped" and by["NP_999999.1"]["uniprot_acc"] == ""
    assert m.release == "2099_01"
    # first stage tried the versioned id for refseq
    assert calls[0][0] == "RefSeq_Protein" and "NP_014368.1" in calls[0][2]


def test_multiple_hits_reviewed_first_and_n_hits(monkeypatch):
    run, _ = fake_rows({("FlyBase", "FBgn0000008"): [("A0A1", "a", "unreviewed", "7227", "1087"),
                                                     ("Q9W2", "a", "reviewed", "7227", "900"),
                                                     ("Q0E8", "a", "unreviewed", "7227", "1329")]})
    monkeypatch.setattr(mi, "run_job", run)
    rows = mi.Mapper().map([("FBgn0000008", "flybase_gene", None)])
    assert [r["uniprot_acc"] for r in rows] == ["Q9W2", "Q0E8", "A0A1"]
    assert all(r["n_hits"] == 3 for r in rows)


def test_gene_name_needs_taxon_and_default_taxon(monkeypatch):
    run, calls = fake_rows({("SGD", "YAL001C"): [], ("Gene_Name", "YAL001C"): [("P34111", "TFC3", "reviewed", "559292", "1160")]})
    monkeypatch.setattr(mi, "run_job", run)
    rows = mi.Mapper().map([("YAL001C", "sgd_systematic", None)])
    assert rows[0]["uniprot_acc"] == "P34111" and "taxon=559292" in rows[0]["mapping_method"]
    assert ("Gene_Name", "UniProtKB", ["YAL001C"], "559292") in calls
    # gene_name without taxon cannot be queried -> reported, not dropped
    rows = mi.Mapper().map([("TFC3", "gene_name", None)])
    assert rows[0]["mapping_method"] == "unmapped"


def test_inactive_accession_falls_back_to_uniparc(monkeypatch):
    run, _ = fake_rows({("UniProtKB_AC-ID", "Q9Y3G6"): [("Q9Y3G6", "", "", "", "")]})
    base_run = run

    def run2(from_db, to_db, ids, taxon=None, **kw):
        if to_db == "UniParc":
            return [{"From": i, "To": "UPI000006D6F7"} for i in ids], "x"
        return base_run(from_db, to_db, ids, taxon)

    monkeypatch.setattr(mi, "run_job", run2)
    rows = mi.Mapper().map([("Q9Y3G6", "uniprot", None)])
    assert rows[0]["uniprot_acc"] == "" and rows[0]["uniparc"] == "UPI000006D6F7"
    assert rows[0]["mapping_method"] == "UniParc_fallback"


def test_secondary_accession_note(monkeypatch):
    run, _ = fake_rows({("UniProtKB_AC-ID", "Q99999"): [("P11111", "G", "reviewed", "9606", "100")]})
    monkeypatch.setattr(mi, "run_job", run)
    rows = mi.Mapper().map([("Q99999", "uniprot", None)])
    assert rows[0]["uniprot_acc"] == "P11111" and "secondary_or_merged" in rows[0]["note"]


def test_batching(monkeypatch):
    run, calls = fake_rows({})
    monkeypatch.setattr(mi, "run_job", run)
    ids = [("NP_%06d.1" % i, "refseq_protein", None) for i in range(25)]
    mi.Mapper(batch=10).map(ids)
    first_stage = [c for c in calls if c[0] == "RefSeq_Protein" and "NP_000000.1" in c[2] or len(c[2]) in (10, 5)]
    assert sorted(len(c[2]) for c in calls[:3]) == [5, 10, 10]


def test_read_items(tmp_path):
    f = tmp_path / "ids.tsv"
    f.write_text("id\tid_type\ttaxon\nNP_014368.1\trefseq_protein\t\nNP_014368.1\trefseq_protein\t\nYAL001C\t\t559292\n")
    items = mi.read_items(str(f))
    assert items == [("NP_014368.1", "refseq_protein", None), ("YAL001C", "sgd_systematic", "559292")]
    g = tmp_path / "list.txt"
    g.write_text("# comment\nP12345\nENSP00000269305\n")
    assert mi.read_items(str(g)) == [("P12345", "uniprot", None), ("ENSP00000269305", "ensembl_protein", None)]


@pytest.mark.network
def test_live_mapping():
    rows = mi.Mapper().map([("NP_014368.1", "refseq_protein", None), ("YAL001C", "sgd_systematic", None),
                            ("ENSMUSP00000000001", "ensembl_protein", None), ("Q9Y3G6", "uniprot", None),
                            ("AT1G01010.1", "arabidopsis_agi", None), ("NP_000000000.1", "refseq_protein", None)])
    by = {r["original_id"]: r for r in rows}
    assert by["NP_014368.1"]["uniprot_acc"] == "P02309" and by["NP_014368.1"]["taxid"] == "559292"
    assert by["YAL001C"]["uniprot_acc"] == "P34111" and by["YAL001C"]["gene_name"] == "TFC3"
    assert by["ENSMUSP00000000001"]["uniprot_acc"] == "Q9DC51"
    assert by["Q9Y3G6"]["uniparc"].startswith("UPI")
    assert by["AT1G01010.1"]["uniprot_acc"] == "Q0WV96"
    assert by["NP_000000000.1"]["mapping_method"] == "unmapped"
