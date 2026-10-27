#!/usr/bin/env python
"""Build collection V07_Xenopus_LS (collection XLAEVIS): X. laevis L/S homeolog pairs + X. tropicalis ortholog.

Source (Session et al. 2016 Nature, doi 10.1038/nature19840, PMID 27762356):
  Supplementary Table 7 (`Supplementary Table 7_ Ks.xlsx`, sheet `XT-XL.ortho`): one row per X. tropicalis locus with
  its 1 (singleton) or 2 (homeolog) X. laevis (co-)orthologues, by gene name and JGI model id. Rows classified
  `homeolog` (8,806) are the candidate units. Identification basis (Supplementary Note 8.1): BLASTP best hits filtered
  by synteny and subgenome identity against the BAC-FISH map.

Unit = one `homeolog` row whose three gene names each map to exactly one canonical protein of the pinned panel-V
proteome of the same species (exact gene-name match in UniProt 2026_03; no cross-species mapping), with
  * the L name ending in `.L`, the S name in `.S`, both sharing the X. tropicalis root name (consistency check);
  * L and S proteins not sequence-identical (identical copies excluded, protocol rule 6);
  * a homology sanity check (DIAMOND blastp, L->XT and S->XT best hit must be the unit's XT protein with >=70 %
    identity over >=50 % of the query), guarding against gene-name drift since 2016.
At most 500 passing units are kept (random.Random(20261026).sample over units sorted by id).

Tree per unit: ((L, S)[D tag=ALLO], XT)[S]  -- tier T2, src Session2016. The D carries no lvl: its default
level is its LCA mapping = the XENLA leaf (panels/V/levels.tsv lists internal levels only).
L/S pair: apart at XENLA, merge at Xenopus; each homeolog vs XT: orthologs (merge Xenopus).

Run from the project root: python families/V07_Xenopus_LS/build.py
"""
from __future__ import annotations

import csv
import gzip
import hashlib
import random
import shutil
import subprocess
import sys
import tempfile
from collections import Counter
from pathlib import Path

import openpyxl
import yaml

ROOT = Path(__file__).resolve().parents[2]
FAM = Path(__file__).resolve().parent
ID = "V07_Xenopus_LS"
PAPER = ROOT / "literature/papers/PMID27762356"
SUPP = PAPER / "raw/supp"
XLSX = SUPP / "Supplementary Table 7_ Ks.xlsx"
SHEET = "XT-XL.ortho"
XB = PAPER / "raw/xenbase"
PEP = {"XENLA": XB / "XL_9.1_v1.8.3.2.primaryTranscripts.pep.fa.gz",   # Session 2016 models (Xelaev18*m ids)
       "XENTR": XB / "Xtropicalisv9.0.Named.primaryTrs.pep.fa.gz"}       # X. tropicalis v9.0 models (by name)
SEQ_ID, SEQ_COV, SEQ_MARGIN = 99.0, 80.0, 2.0
SEED, MAX_UNITS = 20261026, 500
QC_PIDENT, QC_QCOV = 70.0, 50.0

SRC = "Session2016"
QUOTES = {
    "synteny_basis": "Finally, we only accepted X. laevis homoeologues whose synteny and subgenome identity agree "
                     "with the BAC FISH map",
    "pairs_defined": "24,419 X. laevis protein-coding genes can be placed in 2:1 or 1:1 correspondence with 15,613 X. "
                     "tropicalis genes, defining 8,806 homeologous pairs of X. laevis genes with X. tropicalis "
                     "orthologs, and 6,807 single copy orthologs.",
    "allo": "The L and S chromosome sets therefore represent the descendants of two distinct diploid progenitors, "
            "confirming the allotetraploid hypothesis even in the absence of extant progenitor species.",
    "ls_naming": "appending an “L” and “S” for the longer and shorter homeologs, respectively",
    "table7": "Supplemental Table 7: Tables calculating the amount of sequence change between Xenopus orthologs.",
}
QUOTE_LOC = {
    "synteny_basis": "08_supp1.pdf (Supplementary Information), Supplementary Note 8.1",
    "pairs_defined": "PMID27762356 fulltext.txt, Results paragraph 'We annotated 45,099 protein-coding genes'",
    "allo": "PMID27762356 fulltext.txt, Results (subgenome-specific transposable elements)",
    "ls_naming": "PMID27762356 fulltext.txt, Results (chromosome nomenclature)",
    "table7": "Supplementary Tables 1-11.docx, legend of Supplemental Table 7",
}


def sha256(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def species_rows() -> dict:
    with open(ROOT / "panels/V/species.tsv") as f:
        return {r["code"]: r for r in csv.DictReader(f, delimiter="\t")}


def read_fasta(p: Path) -> dict:
    seqs, cur = {}, None
    op = gzip.open if str(p).endswith(".gz") else open
    with op(p, "rt") as f:
        for line in f:
            if line.startswith(">"):
                cur = line[1:].split()[0].split("|")[0].removesuffix(".g")
                seqs[cur] = []
            elif cur:
                seqs[cur].append(line.strip())
    return {k: "".join(v).upper().rstrip("*") for k, v in seqs.items()}


def name_index(code: str, sp: dict) -> tuple[dict, dict]:
    """gene name -> [accessions] restricted to the pinned panel FASTA; plus sequences."""
    fasta = ROOT / "panels/V/proteomes" / f"{code}.fa"
    if not fasta.exists():
        fasta = ROOT / sp[code]["fasta_path"]
    seqs = read_fasta(fasta)
    hm = (ROOT / sp[code]["fasta_path"]).parent / "header_map.tsv"
    idx: dict[str, list[str]] = {}
    with open(hm) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if r["accession"] in seqs and r["gene_name"]:
                idx.setdefault(r["gene_name"], []).append(r["accession"])
    return idx, seqs


def _diamond(query: dict, target: dict, k: int) -> dict:
    """query id -> list of (target, pident, qcov%, tcov%, bitscore), best first (best HSP per pair)."""
    diamond = shutil.which("diamond") or str(ROOT / "conda_envs/ortho/bin/diamond")
    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        for name, d in (("q", query), ("t", target)):
            with open(td / f"{name}.fa", "w") as f:
                for a, sq in d.items():
                    f.write(f">{a}\n{sq}\n")
        subprocess.run([diamond, "makedb", "--in", str(td / "t.fa"), "-d", str(td / "t"), "--quiet", "--threads", "4"],
                       check=True)
        subprocess.run([diamond, "blastp", "-q", str(td / "q.fa"), "-d", str(td / "t"), "-o", str(td / "o.tsv"),
                        "--quiet", "--threads", "4", "--max-target-seqs", str(k), "--evalue", "1e-10",
                        "--outfmt", "6", "qseqid", "sseqid", "pident", "qstart", "qend", "qlen", "sstart", "send",
                        "slen", "bitscore"], check=True)
        out: dict[str, list] = {}
        seen = set()
        for line in open(td / "o.tsv"):
            q, t, pid, qs, qe, ql, ss, se, sl, bs = line.rstrip("\n").split("\t")
            if (q, t) in seen:
                continue
            seen.add((q, t))
            out.setdefault(q, []).append((t, float(pid), 100.0 * (int(qe) - int(qs) + 1) / int(ql),
                                          100.0 * (int(se) - int(ss) + 1) / int(sl), float(bs)))
    for v in out.values():
        v.sort(key=lambda x: -x[4])
    return out


def seq_map(query: dict, paper_all: dict, panel: dict) -> dict:
    """Same-species mapping of paper proteins onto the panel proteome.
    -> {qid: (status, acc, method, identity, coverage, margin, rbh, reason)}"""
    by_seq: dict[str, list[str]] = {}
    for a, sq in panel.items():
        by_seq.setdefault(sq, []).append(a)
    res, todo = {}, {}
    for q, sq in query.items():
        hits = by_seq.get(sq, [])
        if len(hits) == 1:
            res[q] = ("mapped", hits[0], "exact", "100.0", "1.000", "", "", "identical sequence")
        elif len(hits) > 1:
            res[q] = ("ambiguous", "", "", "", "", "", "", f"identical to {len(hits)} panel proteins")
        else:
            todo[q] = sq
    fwd = _diamond(todo, panel, 5) if todo else {}
    cand = {}
    for q in todo:
        hs = [h for h in fwd.get(q, []) if min(h[2], h[3]) >= 50.0]
        if not hs:
            res[q] = ("unmapped", "", "", "", "", "", "", "no panel hit")
            continue
        b = hs[0]
        runner = max((h[1] for h in hs[1:]), default=0.0)
        margin = b[1] - runner
        cov = min(b[2], b[3])
        why = []
        if b[1] < SEQ_ID:
            why.append(f"identity {b[1]:.1f}<{SEQ_ID}")
        if cov < SEQ_COV:
            why.append(f"coverage {cov:.0f}<{SEQ_COV}")
        if margin < SEQ_MARGIN:
            why.append(f"margin {margin:.1f}<{SEQ_MARGIN}")
        if why:
            res[q] = ("unmapped" if b[1] < SEQ_ID or cov < SEQ_COV else "ambiguous", b[0], "", f"{b[1]:.2f}",
                      f"{cov / 100:.3f}", f"{margin:.1f}", "", "best hit " + b[0] + ": " + "; ".join(why))
        else:
            cand[q] = (b, margin, cov)
    rev = _diamond({b[0]: panel[b[0]] for b, _, _ in cand.values()}, paper_all, 5) if cand else {}
    for q, (b, margin, cov) in cand.items():
        r = rev.get(b[0], [])
        top = [h for h in r if r and h[4] == r[0][4]]
        rbh = bool(r) and r[0][0] == q and len(top) == 1
        if rbh:
            res[q] = ("mapped", b[0], "sequence_rbh", f"{b[1]:.2f}", f"{cov / 100:.3f}", f"{margin:.1f}", "1",
                      "diamond reciprocal best hit")
        else:
            res[q] = ("ambiguous", "", "", f"{b[1]:.2f}", f"{cov / 100:.3f}", f"{margin:.1f}", "0",
                      f"best hit {b[0]} is not reciprocal")
    return res


def map_name(name: str, idx: dict) -> tuple[str, str, str]:
    """-> (status, accession, reason)"""
    if not name or name == "NA":
        return "unmapped", "", "no name in source"
    if name.upper().startswith("LOC"):
        return "excluded", "", "NCBI LOC placeholder name (denylisted as a name basis)"
    hits = idx.get(name, [])
    if len(hits) == 1:
        return "mapped", hits[0], "gene_name_exact"
    if len(hits) > 1:
        return "ambiguous", "", f"{len(hits)} proteins share the gene name"
    return "unmapped", "", "gene name not found in panel proteome (renamed/retired since 2016?)"


def diamond_qc(units: list[dict], seqs: dict) -> dict:
    """Best hit of each L/S protein against the units' XT proteins -> {acc: (target, pident, qcov)}."""
    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        with open(td / "q.fa", "w") as q, open(td / "t.fa", "w") as t:
            for u in units:
                q.write(f">{u['L']}\n{seqs['XENLA'][u['L']]}\n>{u['S']}\n{seqs['XENLA'][u['S']]}\n")
            for xt in sorted({u["XT"] for u in units}):
                t.write(f">{xt}\n{seqs['XENTR'][xt]}\n")
        diamond = shutil.which("diamond") or str(ROOT / "conda_envs/ortho/bin/diamond")
        subprocess.run([diamond, "makedb", "--in", str(td / "t.fa"), "-d", str(td / "t"), "--quiet", "--threads", "1"],
                       check=True)
        subprocess.run([diamond, "blastp", "-q", str(td / "q.fa"), "-d", str(td / "t"), "-o", str(td / "o.tsv"),
                        "--quiet", "--threads", "1", "--max-target-seqs", "1", "--evalue", "1e-5",
                        "--outfmt", "6", "qseqid", "sseqid", "pident", "qcovhsp", "bitscore"], check=True)
        best: dict[str, tuple] = {}
        for line in open(td / "o.tsv"):
            qid, sid, pid, qc, bs = line.split("\t")
            if qid not in best or float(bs) > best[qid][3]:
                best[qid] = (sid, float(pid), float(qc), float(bs))
    return best


SHARED_SNAPSHOT = FAM / "shared_with_other_families.tsv"


def shared_genes() -> dict:
    """gene_id -> other family id, from every families/*/members.tsv (and units/*/members.tsv) except this one.
    The result is snapshotted to shared_with_other_families.tsv, which is used when families/ is not available
    (e.g. the verifier's isolated rebuild)."""
    out: dict[str, str] = {}
    fams = ROOT / "families"
    others = [d for d in sorted(fams.iterdir()) if d.is_dir() and d.name != ID] if fams.exists() else []
    if not others and SHARED_SNAPSHOT.exists():
        with open(SHARED_SNAPSHOT) as f:
            return {r["gene_id"]: r["family"] for r in csv.DictReader(f, delimiter="\t")}
    for d in others:
        for m in [d / "members.tsv", *sorted(d.glob("units/*/members.tsv"))]:
            if not m.exists():
                continue
            with open(m) as f:
                for r in csv.DictReader(f, delimiter="\t"):
                    if r.get("gene_id"):
                        out.setdefault(r["gene_id"], d.name)
    return out


def write_shared_snapshot(used: dict) -> None:
    with open(SHARED_SNAPSHOT, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["gene_id", "family"])
        for g in sorted(used):
            w.writerow([g, used[g]])


def main() -> int:
    sp = species_rows()
    idx, seqs = {}, {}
    for c in ("XENLA", "XENTR"):
        idx[c], seqs[c] = name_index(c, sp)

    wb = openpyxl.load_workbook(XLSX, read_only=True)
    rows = list(wb[SHEET].iter_rows(values_only=True))
    header = rows[0]
    assert header[:7] == ("XT.GeneName", "XT.SeqID", "Singleton \nor Homology", "XL.GeneName1", "XL.SeqID1",
                          "XL.GeneName2", "XL.SeqID2"), header

    # paper protein sequences (Session 2016 gene models): XL by JGI model id, XT by v9.0 gene name
    paper = {c: read_fasta(PEP[c]) for c in PEP}
    parsed = []
    for i, r in enumerate(rows[1:], start=2):
        if r[2] != "homeolog":
            continue
        xt, xl1, xl2 = (str(x) if x is not None else "" for x in (r[0], r[3], r[5]))
        loc = f"{XLSX.name}:{SHEET}:row {i}"
        uid = "XL_" + "".join(ch if (ch.isalnum() or ch in "_-") else "_" for ch in (str(r[1]) if xt == "unnamed" else xt))
        # key into the paper FASTA: XT v9.0 gene name (unnamed models are keyed by their Xetrov90 id) / XL model id
        names = {"XT": (xt, "XENTR", str(r[1]), xt if xt in paper["XENTR"] else str(r[1]))}
        for nm, mid in ((xl1, r[4]), (xl2, r[6])):
            slot = "L" if nm.endswith(".L") else "S" if nm.endswith(".S") else "?"
            names.setdefault(slot if slot not in names else slot + "2", (nm, "XENLA", str(mid), str(mid)))
        parsed.append((uid, loc, names))
    smap = {}
    for c in ("XENLA", "XENTR"):
        q = {key: paper[c][key] for _, _, ns in parsed for (nm, code, mid, key) in ns.values()
             if code == c and key in paper[c]}
        smap[c] = seq_map(q, paper[c], seqs[c])

    report, cands = [], []
    for uid, loc, names in parsed:
        res = {}
        for slot, (nm, code, mid, key) in names.items():
            sm = smap[code].get(key)
            g = dict(slot=slot, name=nm, code=code, model=mid, identity="", coverage="", margin="", rbh="")
            if sm and sm[0] == "mapped":
                g.update(status="mapped", acc=sm[1], method=sm[2], identity=sm[3], coverage=sm[4], margin=sm[5],
                         rbh=sm[6], reason=sm[7])
            else:
                st, acc, why = map_name(nm, idx[code])
                seqwhy = "paper sequence not found" if sm is None else sm[7]
                g.update(status=st, acc=acc, method="gene_name_exact" if st == "mapped" else "",
                         reason=f"sequence: {seqwhy}; name: {why}")
                if st == "mapped" and sm is not None and sm[1] and sm[1] != acc:
                    g.update(status="ambiguous", acc="", method="",
                             reason=f"name maps to {acc} but best sequence hit is {sm[1]} ({sm[7]})")
                elif st == "mapped" and sm is not None and sm[0] == "ambiguous":
                    g.update(status="ambiguous", acc="", method="", reason=f"sequence ambiguous ({sm[7]}); name not used")
            res[slot] = g
        problem = ""
        if set(res) != {"XT", "L", "S"}:
            problem = "names do not form one .L and one .S homeolog"
        elif any(v["status"] != "mapped" for v in res.values()):
            problem = "not all three genes uniquely mapped"
        elif seqs["XENLA"][res["L"]["acc"]] == seqs["XENLA"][res["S"]["acc"]]:
            problem = "L and S proteins identical (excluded, protocol rule 6)"
        unit = dict(uid=uid, loc=loc, res=res, problem=problem)
        report.append(unit)
        if not problem:
            cands.append(dict(uid=uid, L=res["L"]["acc"], S=res["S"]["acc"], XT=res["XT"]["acc"], unit=unit))

    # duplicate unit ids / accessions used twice -> exclude all involved
    uc = Counter(c["uid"] for c in cands)
    ac = Counter(a for c in cands for a in (c["L"], c["S"], c["XT"]))
    for c in cands:
        if uc[c["uid"]] > 1 or any(ac[a] > 1 for a in (c["L"], c["S"], c["XT"])):
            c["unit"]["problem"] = "protein or unit id shared with another source row"
    cands = [c for c in cands if not c["unit"]["problem"]]

    best = diamond_qc(cands, seqs)
    passing = []
    for c in cands:
        bad = []
        for s in ("L", "S"):
            b = best.get(c[s])
            if not b or b[0] != c["XT"] or b[1] < QC_PIDENT or b[2] < QC_QCOV:
                bad.append(f"{s}:{'nohit' if not b else f'{b[0]} {b[1]:.1f}%id {b[2]:.0f}%cov'}")
        if bad:
            c["unit"]["problem"] = "homology QC failed (" + "; ".join(bad) + ")"
        else:
            c["qc"] = {s: best[c[s]] for s in ("L", "S")}
            passing.append(c)
    passing.sort(key=lambda c: c["uid"])
    chosen = sorted(random.Random(SEED).sample(passing, min(MAX_UNITS, len(passing))), key=lambda c: c["uid"])
    # Erratum 2026-10-27 (PI): drop sampled units sharing any gene with another family (the other curation is richer).
    shared = shared_genes()
    mine = {a for c in passing for a in (c["L"], c["S"], c["XT"])}
    shared = {g: f for g, f in shared.items() if g in mine}
    write_shared_snapshot(shared)
    dropped = {}
    for c in chosen:
        hit = sorted({f"{a}:{shared[a]}" for a in (c["L"], c["S"], c["XT"]) if a in shared})
        if hit:
            dropped[c["uid"]] = "dropped: shares gene(s) with another family (" + ", ".join(hit) + ")"
    chosen = [c for c in chosen if c["uid"] not in dropped]
    chosen_ids = {c["uid"] for c in chosen}

    # ---------------------------------------------------------------- sources
    files = [XLSX, SUPP / "Supplementary Tables 1-11.docx", SUPP / "08_supp1.pdf", PAPER / "fulltext.txt",
             PEP["XENLA"], PEP["XENTR"]]
    source = dict(id=SRC, pmid=27762356, doi="10.1038/nature19840",
                  title="Genome evolution in the allotetraploid frog Xenopus laevis", year=2016,
                  source_method="synteny",
                  locator=f"{XLSX.name} sheet {SHEET} (homeolog rows); Supplementary Note 8.1",
                  files=[dict(path=str(p.relative_to(ROOT)), sha256=sha256(p)) for p in files],
                  quote=QUOTES["synteny_basis"])
    unit_source = {k: v for k, v in source.items() if k != "files"}

    # ---------------------------------------------------------------- outputs
    units_dir = FAM / "units"
    if units_dir.exists():
        shutil.rmtree(units_dir)
    units_dir.mkdir()
    mem_rows, asr_rows = [], []
    for k in ("synteny_basis", "pairs_defined", "allo", "ls_naming", "table7"):
        asr_rows.append(dict(assertion_id=f"RULE_{k}", type="tag", what="compilation_rule", target="collection",
                             genes="", tier="T1", source_id=SRC, locator=QUOTE_LOC[k], quote_key=k, quote=QUOTES[k]))
    mcols = ["gene_id", "species", "original_id_in_paper", "label", "mapping_method", "mapping_identity",
             "mapping_coverage", "mapping_margin", "reciprocal_best_hit", "source", "locator", "unit"]
    for c in chosen:
        u, res = c["uid"], c["unit"]["res"]
        ud = units_dir / u
        ud.mkdir()
        rows_u = []
        for s in ("L", "S", "XT"):
            g = res[s]
            rows_u.append(dict(gene_id=g["acc"], species=g["code"], original_id_in_paper=f"{g['name']}|{g['model']}",
                               label=f"{u}:{s}", mapping_method=g["method"], mapping_identity=g["identity"],
                               mapping_coverage=g["coverage"], mapping_margin=g["margin"],
                               reciprocal_best_hit=g["rbh"], source=SRC, locator=c["unit"]["loc"], unit=u))
        mem_rows.extend(rows_u)
        with open(ud / "members.tsv", "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=mcols[:-1], delimiter="\t", extrasaction="ignore")
            w.writeheader()
            w.writerows(rows_u)
        L, S, XT = c["L"], c["S"], c["XT"]
        tree = (f"(({L},{S})[&&NHX:ev=D:id={u}_LS:tag=ALLO:tier=T2:src={SRC}],{XT})"
                f"[&&NHX:ev=S:id={u}_Xenopus:tier=T2:src={SRC}];\n")
        (ud / "reference_tree.nhx").write_text(tree)
        meta = dict(id=u, name=f"X. laevis homeologs {res['L']['name']}/{res['S']['name']} + X. tropicalis {res['XT']['name']}",
                    panel="V", collection="XLAEVIS", phenomena=["ALLO"], sources=[unit_source], root_level="Xenopus",
                    root_single_copy=False, complete_species=[], protocol_version=1,
                    notes=f"{c['unit']['loc']}; homology QC L->XT {c['qc']['L'][1]:.1f}% id, "
                          f"S->XT {c['qc']['S'][1]:.1f}% id")
        (ud / "family.yaml").write_text(yaml.safe_dump(meta, sort_keys=False, allow_unicode=True))
        asr_rows.append(dict(assertion_id=f"{u}_LS", type="event", what="D_ALLO_homeolog_pair", target=f"{u}_LS",
                             genes=f"{L},{S}", tier="T2", source_id=SRC, locator=c["unit"]["loc"],
                             quote_key="synteny_basis", quote=QUOTES["synteny_basis"]))
        asr_rows.append(dict(assertion_id=f"{u}_Xenopus", type="event", what="S_Xenopus_orthologs",
                             target=f"{u}_Xenopus", genes=f"{L},{S},{XT}", tier="T2", source_id=SRC,
                             locator=c["unit"]["loc"], quote_key="pairs_defined", quote=QUOTES["pairs_defined"]))

    with open(FAM / "members.tsv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=mcols, delimiter="\t")
        w.writeheader()
        w.writerows(mem_rows)
    with open(FAM / "assertions.tsv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(asr_rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(asr_rows)
    n_status = Counter()
    with open(FAM / "mapping_report.tsv", "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["unit", "unit_status", "original_id", "species", "slot", "status", "gene_id", "mapping_method",
                    "mapping_identity", "mapping_coverage", "reason", "locator"])
        for u in report:
            ust = ("unit" if u["uid"] in chosen_ids else dropped[u["uid"]] if u["uid"] in dropped else
                   f"not_a_unit: {u['problem']}" if u["problem"] else "not_sampled (passes all checks)")
            for s, g in u["res"].items():
                n_status[g["status"]] += 1
                w.writerow([u["uid"], ust, f"{g['name']}|{g['model']}", g["code"], s, g["status"], g["acc"],
                            g["method"], g["identity"], g["coverage"], g["reason"], u["loc"]])
    probs = Counter(u["problem"].split(" (")[0] for u in report if u["problem"])

    coll = dict(id=ID, name="Xenopus laevis L/S homeolog pairs with X. tropicalis orthologs (Session et al. 2016)",
                panel="V", collection="XLAEVIS",
                layout="one unit (= one X. tropicalis locus with an L/S homeolog pair) per directory under units/; "
                       "each unit is a loadable family (hobench.reference.load_family)",
                phenomena=["ALLO"], sources=[source], root_level="Xenopus", root_single_copy=False,
                complete_species=[], decoys=[], distinct_from=[],
                tier_summary=dict(T2_D_ALLO=len(chosen), T2_S_Xenopus=len(chosen)),
                sampling=dict(seed=SEED, max_units=MAX_UNITS, candidates_passing=len(passing),
                              units=len(chosen), dropped_shared_with_other_families=sorted(dropped), method="random.Random(seed).sample over passing units sorted by id"),
                counts=dict(source_homeolog_rows=len(report), gene_slots_by_status=dict(n_status),
                            unit_rejections=dict(probs)),
                protocol_version=1,
                notes="T2 rule: Session 2016 Table 7 homeolog rows (BLASTP + synteny/subgenome check vs BAC-FISH "
                      "map) -> D at the XENLA leaf tagged ALLO (L and S descend from two diploid progenitors, so the "
                      "L/S split is not a duplication inside the X. laevis lineage; per PI/SYNTHESIS A6 it is encoded "
                      "on the XENLA branch) and S between the pair and the X. tropicalis gene at Xenopus.")
    (FAM / "collection.yaml").write_text(yaml.safe_dump(coll, sort_keys=False, allow_unicode=True))
    print(f"rows={len(report)} candidates={len(cands)} passing={len(passing)} units={len(chosen)} "
          f"slots={dict(n_status)} rejections={dict(probs)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
