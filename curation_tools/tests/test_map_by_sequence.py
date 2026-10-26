import csv
import random
import shutil

import pytest

from curation_tools import map_by_sequence as ms
from curation_tools.common import read_fasta, write_fasta

pytestmark = pytest.mark.skipif(shutil.which("diamond") is None and not ms.common.have_exe("diamond"),
                                reason="diamond not installed")


def mutate(s, positions, aa="W"):
    s = list(s)
    for p in positions:
        s[p] = aa if s[p] != aa else "A"
    return "".join(s)


def queries(fix):
    t = dict((h, s) for h, s in read_fasta(fix / "mini_proteome.fasta"))
    rnd = random.Random(1)
    rand = "".join(rnd.choice("ACDEFGHIKLMNPQRSTVWY") for _ in range(200))
    return [("exact_T2", t["T2"]),
            ("exact_T3", t["T3_dup_of_T4"]),
            ("mut_T1", mutate(t["T1"], [10, 50, 100, 150, 200])),                  # ~98.5 % to T1 only
            ("amb_T3T4", mutate(t["T3_dup_of_T4"], [20])),                          # 1 mismatch to T3, 2 to T4
            ("partial_T2", t["T2"][:30] + "G"),                                    # short fragment w/ 1 mismatch
            ("random", rand)], t


def test_statuses(fix):
    qs, t = queries(fix)
    rows = {r["query_id"]: r for r in ms.map_sequences(qs, fix / "mini_proteome.fasta", threads=2)}
    assert len(rows) == len(qs)                           # nothing dropped
    assert rows["exact_T2"]["status"] == "exact" and rows["exact_T2"]["target_id"] == "T2"
    assert rows["exact_T2"]["mapping_identity"] == 100.0 and rows["exact_T2"]["method"] == "exact_sha256"
    assert rows["exact_T3"]["target_id"] == "T3_dup_of_T4" and rows["exact_T3"]["status"] == "exact"
    m = rows["mut_T1"]
    assert m["status"] == "unique" and m["target_id"] == "T1" and m["method"] == "diamond"
    assert 97.5 < m["mapping_identity"] < 99.0 and m["query_cov"] > 0.99 and m["target_cov"] > 0.99
    a = rows["amb_T3T4"]
    assert a["status"] == "ambiguous" and a["target_id"] == "T3_dup_of_T4" and "T4_dup_of_T3" in a["ambiguous_with"]
    assert rows["random"]["status"] == "unmapped" and rows["random"]["target_id"] == ""
    p = rows["partial_T2"]
    assert p["status"] == "partial" and p["flags"] == "low_target_cov" and p["target_id"] == "T2"   # fragment of T2


def test_exact_ambiguous_identical_targets(tmp_path, fix):
    seqs = dict(read_fasta(fix / "mini_proteome.fasta"))
    tgt = tmp_path / "t.fa"
    write_fasta([("A", seqs["T2"]), ("B", seqs["T2"]), ("C", seqs["T1"])], tgt)
    rows = ms.map_sequences([("q", seqs["T2"].lower() + "*")], tgt, threads=2)   # normalisation: case + stop codon
    assert rows[0]["status"] == "exact_ambiguous" and rows[0]["n_exact_targets"] == 2
    assert rows[0]["target_id"] == "A,B"


def test_duplicate_target_ids_rejected(tmp_path):
    tgt = tmp_path / "t.fa"
    write_fasta([("A", "MKV"), ("A", "MKW")], tgt)
    with pytest.raises(ValueError):
        ms.map_sequences([("q", "MKV")], tgt)


def test_cli_roundtrip(tmp_path, fix):
    qs, _ = queries(fix)
    qf = tmp_path / "q.fa"
    write_fasta(qs, qf)
    out = tmp_path / "o.tsv"
    assert ms.main(["--query", str(qf), "--target", str(fix / "mini_proteome.fasta"), "-o", str(out), "--threads", "2"]) == 0
    rows = list(csv.DictReader(open(out), delimiter="\t"))
    assert [r["query_id"] for r in rows] == [q for q, _ in qs] and set(ms.OUT_COLS) == set(rows[0])


@pytest.mark.skipif(ms.common.have_exe("mmseqs") is None, reason="mmseqs not installed")
def test_mmseqs_backend(fix):
    qs, _ = queries(fix)
    rows = {r["query_id"]: r for r in ms.map_sequences(qs, fix / "mini_proteome.fasta", tool="mmseqs", threads=2)}
    assert rows["mut_T1"]["target_id"] == "T1" and 97.5 < rows["mut_T1"]["mapping_identity"] < 99.0


@pytest.mark.network
def test_fetch_sequences_live():
    got, failed = ms.fetch_sequences(["NP_014368.1", "P34111", "Q9Y3G6", "YAL001C"])
    assert not failed
    assert got["NP_014368.1"][1].startswith("MSGRGKGGKGLGKGG") and got["P34111"][0] == "UniProtKB"
    assert got["Q9Y3G6"][0].startswith("UniParc")
    assert got["YAL001C"][1] == got["P34111"][1]
