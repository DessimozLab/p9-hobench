"""Tests for wrapper/: tree checks, timing parser, subset builder, perturbations, of_flat.

Everything runs in pytest tmp dirs (never in panels/ or runs/). The subset test runs real
DIAMOND on a tiny synthetic panel and is skipped if diamond is not installed.
"""
from __future__ import annotations

import gzip
import io
import json
import random
import shutil
import tarfile
from pathlib import Path

import pytest

from hobench.newick import read_newick
from wrapper import common as C
from wrapper import perturb, run, subset

TREE = "((YEAST,CANGA)PostWGD,(KLULA,ASHGO)KLE)Saccharomycetaceae;"
CODES = ["YEAST", "CANGA", "KLULA", "ASHGO"]
AA = "ACDEFGHIKLMNPQRSTVWY"


def rand_seq(rng, n):
    return "M" + "".join(rng.choice(AA) for _ in range(n - 1))


def mutate(rng, s, rate):
    return "".join(rng.choice(AA) if (i > 0 and rng.random() < rate) else c for i, c in enumerate(s))


@pytest.fixture
def toy_panel(tmp_path):
    """4 species x 60 random proteins + one 'family' (fam0..fam3 per species, homologous)."""
    rng = random.Random(1)
    root = tmp_path / "panels"
    d = root / "TOY"
    (d / "proteomes").mkdir(parents=True)
    fam = rand_seq(rng, 300)
    for i, c in enumerate(CODES):
        recs = [(f"tr|{c}FAM{i}|FAM_{c} family member", mutate(rng, fam, 0.15))]
        # a domain-only homolog: 80 aa of the family glued to 400 random aa
        recs.append((f"tr|{c}DOM|DOM_{c} domain carrier", rand_seq(rng, 200) + mutate(rng, fam[100:180], 0.1)
                     + rand_seq(rng, 250)))
        recs += [(f"tr|{c}{k:03d}|P{k}_{c} random", rand_seq(rng, rng.randint(80, 400))) for k in range(60)]
        C.write_fasta(d / "proteomes" / f"{c}.fa", recs)
    (d / "species.tsv").write_text("code\ttaxid\n" + "".join(f"{c}\t{i}\n" for i, c in enumerate(CODES)))
    (d / "species_tree.nwk").write_text(TREE + "\n")
    return root


# ----------------------------------------------------------------------------- common
def test_validate_tree_ok_and_errors():
    t = read_newick(TREE)
    assert C.validate_species_tree(t, CODES) == []
    bad = read_newick("((YEAST,CANGA)KE,(KLULA,ASHGO))Saccharomycetaceae;")
    probs = C.validate_species_tree(bad, CODES)
    assert any("'KE'" in p for p in probs) and any("''" in p for p in probs)
    multi = read_newick("(YEAST,CANGA,(KLULA,ASHGO)KLE)Root;")
    assert any("binary" in p for p in C.validate_species_tree(multi, CODES))
    assert any("only_panel" in p for p in C.validate_species_tree(t, CODES + ["XXX"]))


def test_prune_tree_removes_unary_levels():
    t = read_newick(TREE)
    p = C.prune_tree(t, {"YEAST", "KLULA", "ASHGO"})
    assert C.tree_clades(p) == {"KLE": frozenset({"KLULA", "ASHGO"}),
                                "Saccharomycetaceae": frozenset({"YEAST", "KLULA", "ASHGO"})}
    assert sorted(p.leaf_names()) == ["ASHGO", "KLULA", "YEAST"]


def test_parse_time_v():
    txt = ("\tUser time (seconds): 12.50\n\tSystem time (seconds): 1.25\n"
           "\tPercent of CPU this job got: 350%\n"
           "\tElapsed (wall clock) time (h:mm:ss or m:ss): 1:02:03\n"
           "\tMaximum resident set size (kbytes): 123456\n\tExit status: 0\n")
    d = C.parse_time_v(txt)
    assert d == {"user_s": 12.5, "sys_s": 1.25, "cpu_percent": "350%", "wall_clock": "1:02:03",
                 "max_rss_kb": 123456, "exit_status": 0}


def test_stable_unit_nested_and_deterministic():
    ids = [f"id{i}" for i in range(2000)]
    a = {i for i in ids if C.stable_unit(5, "bg", i) < 0.05}
    b = {i for i in ids if C.stable_unit(5, "bg", i) < 0.10}
    assert a <= b and 50 < len(a) < 150
    assert a == {i for i in reversed(ids) if C.stable_unit(5, "bg", i) < 0.05}
    assert a != {i for i in ids if C.stable_unit(6, "bg", i) < 0.05}


def test_load_panel_and_missing(toy_panel):
    p = C.load_panel("TOY", "full", toy_panel)
    assert p.codes == CODES and all(f.exists() for f in p.fastas.values())
    with pytest.raises(C.WrapperError):
        C.load_panel("TOY", "subset", toy_panel)


# ----------------------------------------------------------------------------- subset
needs_diamond = pytest.mark.skipif(shutil.which("diamond", path=C.tool_env()["PATH"]) is None,
                                   reason="diamond not installed")


@needs_diamond
def test_subset_categories_and_cap(toy_panel, tmp_path):
    seeds = tmp_path / "seeds.txt"
    seeds.write_text("YEASTFAM0\nNOT_IN_PANEL\n")
    out = tmp_path / "sub"
    rc = subset.main(["--panel", "TOY", "--seeds", str(seeds), "--panels-root", str(toy_panel),
                      "--out", str(out), "--work-root", str(tmp_path / "w"), "--threads", "2",
                      "--bg-fraction", "0.1", "--domain-cap", "2"])
    assert rc == 0
    man = json.loads((out / "manifest.json").read_text())
    assert man["seeds"]["missing"] == ["NOT_IN_PANEL"]
    rows = [ln.split("\t") for ln in (out / "members.tsv").read_text().splitlines()[1:]]
    cat = {r[0]: r[2] for r in rows}
    fam = {f"tr|{c}FAM{i}|FAM_{c}" for i, c in enumerate(CODES)}
    dom = {f"tr|{c}DOM|DOM_{c}" for c in CODES}
    assert all(cat.get(f) == "neighbourhood" for f in fam)          # positive control
    assert sum(cat.get(d) == "domain_only" for d in dom) == 2        # capped domain decoys
    assert not any(cat.get(d) == "neighbourhood" for d in dom)
    n_bg = sum(v == "background" for v in cat.values())
    assert 0 < n_bg < 0.25 * 62 * 4
    # FASTA contents == members, per species
    for c in CODES:
        ids = {C.fasta_id(h) for h, _ in C.read_fasta(out / f"{c}.fa")}
        assert ids == {r[0] for r in rows if r[1] == c}
    # background is independent of the seeds and deterministic: rerun with another seed list
    seeds.write_text("KLULAFAM2\n")
    out2 = tmp_path / "sub2"
    subset.main(["--panel", "TOY", "--seeds", str(seeds), "--panels-root", str(toy_panel), "--out", str(out2),
                 "--work-root", str(tmp_path / "w"), "--threads", "2", "--bg-fraction", "0.1", "--domain-cap", "2"])
    rows2 = [ln.split("\t") for ln in (out2 / "members.tsv").read_text().splitlines()[1:]]
    bg1 = {r[0] for r in rows if r[2] == "background"} | {r[0] for r in rows if r[2] != "background"}
    bg2 = {r[0] for r in rows2 if r[2] == "background"}
    assert bg2 <= bg1  # same background keys; only core membership differs
    # cap: total limited to 10 proteins -> background shrinks to fit
    out3 = tmp_path / "sub3"
    subset.main(["--panel", "TOY", "--seeds", str(seeds), "--panels-root", str(toy_panel), "--out", str(out3),
                 "--work-root", str(tmp_path / "w"), "--threads", "2", "--bg-fraction", "0.5",
                 "--domain-cap", "0", "--cap", "10"])
    assert json.loads((out3 / "manifest.json").read_text())["n_subset"] == 10
    with pytest.raises(C.WrapperError):  # neighbourhood alone exceeds cap
        subset.main(["--panel", "TOY", "--seeds", str(seeds), "--panels-root", str(toy_panel),
                     "--out", str(tmp_path / "sub4"), "--work-root", str(tmp_path / "w"), "--threads", "2",
                     "--cap", "2"])


# ----------------------------------------------------------------------------- perturbations
def _seqs(path):
    return {C.fasta_id(h): s for h, s in C.read_fasta(path)}


def test_fragment(toy_panel, tmp_path):
    genes = tmp_path / "genes.txt"
    ids = [f"{c}{k:03d}" for c in CODES for k in range(10)]  # accessions; 40 genes
    genes.write_text("\n".join(ids + ["MISSING1"]) + "\n")
    perturb.main(["fragment", "--panel", "TOY", "--genes", str(genes), "--panels-root", str(toy_panel)])
    out = toy_panel / "TOY_frag"
    man = json.loads((out / "perturbation.json").read_text())
    assert man["genes_not_found"] == ["MISSING1"]
    n_el = man["n_eligible"]
    assert man["n_fragmented"] == -(-n_el // 10)
    for ch in man["changes"]:
        orig = _seqs(toy_panel / "TOY" / "proteomes" / f"{ch['species']}.fa")[ch["orig_id"]]
        new = _seqs(out / "proteomes" / f"{ch['species']}.fa")
        f1, f2 = (f["id"] for f in ch["fragments"])
        assert ch["orig_id"] not in new and new[f1] + new[f2] == orig
        assert min(len(new[f1]), len(new[f2])) >= perturb.MIN_FRAG
        assert f1 == C.accession(ch["orig_id"]) + "_f1"
    # untouched species are symlinks; total protein count grows by n_fragmented
    tot = lambda d: sum(len(_seqs(p)) for p in (d / "proteomes").glob("*.fa"))
    assert tot(out) == tot(toy_panel / "TOY") + man["n_fragmented"]
    # deterministic
    perturb.main(["fragment", "--panel", "TOY", "--genes", str(genes), "--panels-root", str(toy_panel),
                  "--name", "frag2"])
    man2 = json.loads((toy_panel / "TOY_frag2" / "perturbation.json").read_text())
    assert man2["changes"] == man["changes"]


def test_contaminate(toy_panel, tmp_path):
    src = tmp_path / "human.fa"
    rng = random.Random(3)
    C.write_fasta(src, [(f"sp|H{i:04d}|X{i}_HUMAN", rand_seq(rng, 120)) for i in range(50)])
    perturb.main(["contaminate", "--panel", "TOY", "--target", "KLULA", "--source-fasta", str(src),
                  "--n", "7", "--source-label", "HUMAN", "--panels-root", str(toy_panel)])
    out = toy_panel / "TOY_contam"
    man = json.loads((out / "perturbation.json").read_text())
    new = _seqs(out / "proteomes" / "KLULA.fa")
    assert man["n_inserted"] == 7 and len(new) == 62 + 7
    assert all(ch["new_id"] in new and ch["new_id"].startswith("CONTAM_HUMAN_H") for ch in man["changes"])
    assert (out / "proteomes" / "YEAST.fa").is_symlink()


def test_alttree_and_remove(toy_panel, tmp_path):
    alt = toy_panel / "TOY" / "species_tree_alt.nwk"
    alt.write_text("((YEAST,KLULA)Alt1,(CANGA,ASHGO)Alt2)Saccharomycetaceae;\n")
    perturb.main(["alttree", "--panel", "TOY", "--panels-root", str(toy_panel)])
    man = json.loads((toy_panel / "TOY_alttree" / "perturbation.json").read_text())
    assert man["clades"]["shared"] == [["Saccharomycetaceae", "Saccharomycetaceae"]]
    assert {x[0] for x in man["clades"]["only_parent"]} == {"PostWGD", "KLE"}
    perturb.main(["remove", "--panel", "TOY", "--taxa", "CANGA", "--panels-root", str(toy_panel)])
    out = toy_panel / "TOY_taxrm"
    man = json.loads((out / "perturbation.json").read_text())
    assert man["removed_protein_counts"] == {"CANGA": 62}
    assert man["removed_levels"] == ["PostWGD"]
    assert man["clades"]["only_parent"] == [] and man["clades"]["only_new"] == []
    p = C.load_panel("TOY_taxrm", "full", toy_panel)
    assert p.codes == ["YEAST", "KLULA", "ASHGO"]
    assert C.validate_species_tree(read_newick(p.tree_path), p.codes) == []


# ----------------------------------------------------------------------------- of_flat
def test_of_flat_from_fake_run(toy_panel, tmp_path):
    runs = tmp_path / "runs"
    src = runs / "TOY" / "full" / "orthofinder_default"
    src.mkdir(parents=True)
    og = "Orthogroup\tYEAST\tCANGA\tKLULA\tASHGO\nOG0000000\ty1, y2\tc1\t\ta1\n"
    un = "Orthogroup\tYEAST\tCANGA\tKLULA\tASHGO\nOG0000001\t\t\tk9\t\n"
    with tarfile.open(src / "of_core.tar.gz", "w:gz") as tf:
        for name, txt in (("Orthogroups/Orthogroups.tsv", og), ("Orthogroups/Orthogroups_UnassignedGenes.tsv", un)):
            b = txt.encode()
            ti = tarfile.TarInfo(name)
            ti.size = len(b)
            tf.addfile(ti, io.BytesIO(b))
    (src / "manifest.json").write_text(json.dumps({"outputs": {"of_core.tar.gz": {"path": "of_core.tar.gz"}},
                                                   "inputs": {}}))
    rc = run.main(["--tool", "of_flat", "--panel", "TOY", "--panels-root", str(toy_panel),
                   "--runs-root", str(runs), "--work-root", str(tmp_path / "w")])
    assert rc == 0
    rows = gzip.open(runs / "TOY" / "full" / "of_flat_default" / "hogs.tsv.gz", "rt").read().splitlines()
    assert rows[0] == "level\thog_id\tgene\tspecies"
    body = sorted(r.split("\t") for r in rows[1:])
    assert body == sorted([["Saccharomycetaceae", "OG0000000", g, s] for g, s in
                           (("y1", "YEAST"), ("y2", "YEAST"), ("c1", "CANGA"), ("a1", "ASHGO"))]
                          + [["Saccharomycetaceae", "OG0000001", "k9", "KLULA"]])


def test_assemble_hogmaps_filters(tmp_path):
    fa = tmp_path / "fa"
    fa.mkdir()
    C.write_fasta(fa / "AAA.fa", [("p1 x", "MAAA"), ("p2 y", "MCCC")])
    hm = tmp_path / "hm"
    hm.mkdir()
    (hm / "AAA.fa.hogmap").write_text("!omamer-version: 2.1.0\nqseqid\thogid\np1\tH1\np2\tH2\np3\tH3\n")
    recs = run.assemble_hogmaps(fa, ["AAA", "BBB"], hm, tmp_path / "out", Path("/nonexistent"), tmp_path)
    txt = (tmp_path / "out" / "AAA.fa.hogmap").read_text()
    assert "p3" not in txt and "p1\tH1" in txt and txt.startswith("!omamer")
    assert recs["AAA"]["n_searched_now"] == 0 and recs["BBB"]["source"] is None
