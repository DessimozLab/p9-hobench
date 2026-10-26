"""Shared helpers for Panel P family builds (P01-P06). Not a family (no family.yaml).

All network results are cached as TSV/JSON inside the calling family's ``cache/`` directory so that
``build.py`` reruns offline and deterministically once the cache exists.
"""
from __future__ import annotations

import csv
import hashlib
import json
import re
import sys
import time
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from hobench.newick import Node  # noqa: E402
from hobench.speciestree import SpeciesTree  # noqa: E402

PANEL = ROOT / "panels" / "P"


def sha256(path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def species_table() -> dict[str, dict]:
    with open(PANEL / "species.tsv") as fh:
        return {r["code"]: r for r in csv.DictReader(fh, delimiter="\t")}


def species_tree() -> SpeciesTree:
    return SpeciesTree.from_newick(PANEL / "species_tree.nwk")


_PROT: dict[str, dict[str, str]] = {}
_HDR: dict[str, dict[str, dict]] = {}


def proteome(code: str) -> dict[str, str]:
    """accession -> sequence of the pinned panel proteome."""
    if code not in _PROT:
        path = ROOT / species_table()[code]["fasta_path"]
        seqs, acc, buf = {}, None, []
        with open(path) as fh:
            for line in fh:
                if line.startswith(">"):
                    if acc:
                        seqs[acc] = "".join(buf)
                    acc, buf = line[1:].split()[0], []
                else:
                    buf.append(line.strip())
        if acc:
            seqs[acc] = "".join(buf)
        _PROT[code] = seqs
    return _PROT[code]


def header_map(code: str) -> dict[str, dict]:
    if code not in _HDR:
        path = ROOT / species_table()[code]["fasta_path"]
        hm = path.parent / "header_map.tsv"
        with open(hm) as fh:
            _HDR[code] = {r["accession"]: r for r in csv.DictReader(fh, delimiter="\t")}
    return _HDR[code]


# --------------------------------------------------------------------------- network (cached)
def http_json(url: str, tries: int = 4):
    for k in range(tries):
        try:
            req = urllib.request.Request(url, headers={"Content-Type": "application/json",
                                                       "Accept": "application/json"})
            with urllib.request.urlopen(req, timeout=60) as r:
                return json.loads(r.read().decode())
        except urllib.error.HTTPError as e:
            if e.code in (400, 404):
                return None
            time.sleep(2 + 3 * k)
        except Exception:
            time.sleep(2 + 3 * k)
    return None


class JsonCache:
    def __init__(self, path: Path):
        self.path = Path(path)
        self.data = json.loads(self.path.read_text()) if self.path.exists() else {}

    def get(self, key, fn):
        if key not in self.data:
            self.data[key] = fn()
            self.save()
        return self.data[key]

    def save(self):
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self.path.write_text(json.dumps(self.data, indent=1, sort_keys=True))


ENS = "https://rest.ensembl.org"


def ensembl_overlap_genes(cache: JsonCache, species: str, region: str):
    """protein-coding Ensembl genes overlapping region 'chr:start-end' (current assembly)."""
    url = f"{ENS}/overlap/region/{species}/{region}?feature=gene;content-type=application/json"
    res = cache.get("overlap|" + species + "|" + region, lambda: http_json(url))
    return [g for g in (res or []) if g.get("biotype") == "protein_coding"]


def ensembl_lookup(cache: JsonCache, eid: str):
    url = f"{ENS}/lookup/id/{eid}?content-type=application/json"
    return cache.get("lookup|" + eid, lambda: http_json(url))


def ensembl_archive(cache: JsonCache, eid: str):
    url = f"{ENS}/archive/id/{eid}?content-type=application/json"
    return cache.get("archive|" + eid, lambda: http_json(url))


def uniprot_map(cache: JsonCache, ids: list[tuple[str, str]], taxon: str | None = None):
    """ids: [(original_id, id_type)] -> {original_id: [rows]} via curation_tools.map_ids (cached)."""
    from curation_tools.map_ids import Mapper
    todo = [(i, t) for i, t in ids if "uniprot|" + i not in cache.data]
    if todo:
        rows = Mapper(taxon=taxon).map([(i, t, taxon) for i, t in todo])
        got: dict[str, list] = {}
        for r in rows:
            got.setdefault(r["original_id"], []).append(r)
        for i, _ in todo:
            cache.data["uniprot|" + i] = got.get(i, [])
        cache.save()
    return {i: cache.data["uniprot|" + i] for i, _ in ids}


def pick_in_proteome(code: str, rows: list[dict]):
    """Return (acc, note) for mapping rows that hit the panel proteome: unique hit -> acc; else ('', reason)."""
    accs = sorted({r["uniprot_acc"] for r in rows if r.get("uniprot_acc")})
    prot = proteome(code)
    inp = [a for a in accs if a in prot]
    if len(inp) == 1:
        return inp[0], f"{len(accs)} UniProt hit(s), 1 in panel proteome"
    if not accs:
        return "", "no UniProt match"
    if not inp:
        return "", f"UniProt hit(s) {accs[:4]} not in panel proteome (non-canonical/other entries)"
    return "", f"ambiguous: {len(inp)} panel accessions {inp[:4]}"


# --------------------------------------------------------------------------- trees
def nhx(node: Node) -> str:
    """Serialise with all features (ev first) — supports id/tier/lvl_lo/lvl_hi keys."""
    order = ["ev", "id", "lvl", "lvl_lo", "lvl_hi", "tag", "tier", "src"]

    def rec(n: Node) -> str:
        s = ""
        if n.children:
            s = "(" + ",".join(rec(c) for c in n.children) + ")"
        s += n.name or ""
        if n.children and n.features:
            keys = [k for k in order if k in n.features] + sorted(k for k in n.features if k not in order)
            s += "[&&NHX:" + ":".join(f"{k}={n.features[k]}" for k in keys if n.features[k] != "") + "]"
        return s
    return rec(node) + ";"


def species_s_tree(st: SpeciesTree, genes: list[tuple[str, str]], feats: dict | None = None,
                   inparalog_ev: str = "U") -> Node:
    """S-node tree following the species tree over genes [(gene, species)]; several genes of one
    species are joined under an ``inparalog_ev`` node."""
    by: dict[str, list[str]] = {}
    for g, s in genes:
        by.setdefault(s, []).append(g)
    if len(genes) == 1:
        return Node(genes[0][0])
    pruned = st.prune_to(by)

    def conv(n: Node) -> Node:
        if n.is_leaf:
            gs = by[n.name]
            if len(gs) == 1:
                return Node(gs[0])
            out = Node(features={"ev": inparalog_ev, **(feats or {})})
            for g in gs:
                out.add_child(Node(g))
            return out
        out = Node(features={"ev": "S", **(feats or {})})
        for c in n.children:
            out.add_child(conv(c))
        return out
    return conv(pruned)


def mk(ev: str, children: list[Node], **feats) -> Node:
    kids = [c for c in children if c is not None]
    if len(kids) == 1:
        return kids[0]
    if not kids:
        return None
    n = Node(features={"ev": ev, **{k: str(v) for k, v in feats.items() if v not in (None, "")}})
    for c in kids:
        n.add_child(c)
    return n


def quote_in(text: str, quote: str) -> bool:
    norm = lambda s: re.sub(r"\s+", " ", s)
    return norm(quote) in norm(text)


# --------------------------------------------------------------------------- generic ID -> panel mapping
ENS_SP = {"LEPOC": "lepisosteus_oculatus", "DANRE": "danio_rerio", "ORYLA": "oryzias_latipes",
          "HUMAN": "homo_sapiens", "MOUSE": "mus_musculus", "CHICK": "gallus_gallus",
          "XENTR": "xenopus_tropicalis", "LATCH": "latimeria_chalumnae", "CALMI": "callorhinchus_milii",
          "PETMA": "petromyzon_marinus", "CIOIN": "ciona_intestinalis", "SALSA": "salmo_salar"}


def id_type(oid: str) -> str | None:
    if re.match(r"^ENS[A-Z]*P\d{11}$", oid):
        return "ensembl_protein"
    if re.match(r"^ENS[A-Z]*T\d{11}$", oid):
        return "ensembl_transcript"
    if re.match(r"^ENS[A-Z]*G\d{11}$", oid):
        return "ensembl_gene"
    if re.match(r"^[NX]P_\d+(\.\d+)?$", oid):
        return "refseq_protein"
    if re.match(r"^[NX]M_\d+(\.\d+)?$", oid):
        return "refseq_nucleotide"
    if re.match(r"^([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})(\.\d+)?$", oid):
        return "uniprot"
    if re.match(r"^[A-Z]{3}\d{5}(\.\d+)?$", oid):
        return "genbank_protein"
    return None


def map_one(cache: JsonCache, code: str, oid: str, region: str | None = None) -> dict:
    """Map a paper identifier (or, if ``region`` 'chr:start-end' is given, a locus of the current Ensembl
    assembly of that species) to exactly one accession of the panel proteome. Returns dict(status, gene_id,
    mapping_method, reason)."""
    out = dict(status="unmapped", gene_id="", mapping_method="", reason="", mapping_identity="", mapping_coverage="")
    oid = (oid or "").strip()
    t = id_type(re.sub(r"\.\d+$", "", oid) if t_is_uniprot_versioned(oid) else oid)
    if t:
        q = re.sub(r"\.\d+$", "", oid) if t == "uniprot" else oid
        rows = uniprot_map(cache, [(q, t)])[q]
        acc, note = pick_in_proteome(code, rows)
        if not acc and t == "uniprot" and q in proteome(code):
            acc, note = q, "accession present in panel proteome"
        if not acc and t.startswith("ensembl"):
            arc = ensembl_archive(cache, oid)
            reps = [p["stable_id"] for p in (arc or {}).get("possible_replacement") or []]
            cur = None
            if arc and arc.get("is_current") == "1":
                # current ID: try its gene (lookup) -> UniProt
                lk = ensembl_lookup(cache, oid)
                if lk and lk.get("Parent") and t != "ensembl_gene":
                    cur = lk["Parent"]
                    if lk.get("object_type") == "Translation":
                        lk2 = ensembl_lookup(cache, cur)
                        cur = (lk2 or {}).get("Parent") or cur
            elif len(reps) == 1:
                cur = reps[0]
            if cur and cur != oid:
                acc2, note2 = pick_in_proteome(code, uniprot_map(cache, [(cur, id_type(cur))])[cur])
                if acc2:
                    acc, note = acc2, f"{oid} -> Ensembl {cur}; {note2}"
                else:
                    note += f"; via Ensembl {cur}: {note2}"
            elif arc and arc.get("is_current") != "1":
                note += f"; retired Ensembl ID ({len(reps)} archive replacements)"
        if acc:
            meth = rows[0]["mapping_method"] if rows else "ensembl"
            out.update(status="mapped", gene_id=acc, mapping_method="exact_id:" + meth, reason=note,
                       mapping_identity="100", mapping_coverage="100")
        else:
            out.update(status="ambiguous" if note.startswith("ambiguous") else "unmapped", reason=note)
        return out
    if region and code in ENS_SP:
        genes = ensembl_overlap_genes(cache, ENS_SP[code], region)
        if len(genes) == 1:
            g = genes[0]["id"]
            acc, note = pick_in_proteome(code, uniprot_map(cache, [(g, "ensembl_gene")])[g])
            if acc:
                out.update(status="mapped", gene_id=acc, mapping_method=f"ensembl_locus:{region}->{g}",
                           reason=f"Ensembl gene {g} {genes[0].get('external_name', '')}; {note}",
                           mapping_identity="100", mapping_coverage="100")
            else:
                out.update(reason=f"locus {region} -> {g}; {note}")
        else:
            out.update(status="ambiguous" if genes else "unmapped",
                       reason=f"{len(genes)} protein-coding Ensembl genes at {region}")
        return out
    out["reason"] = "no gene-model identifier (genomic contig/scaffold/EST/manual prediction)"
    return out


def t_is_uniprot_versioned(oid: str) -> bool:
    return bool(re.match(r"^([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9][A-Z][A-Z0-9]{2}[0-9])\.\d+$", oid))


# --------------------------------------------------------------------------- sequence mapping (protocol §0.5)
import hashlib as _hl
import subprocess as _sp
import os as _os

DIAMOND = str(ROOT / "conda_envs" / "ortho" / "bin" / "diamond")


def norm_seq(s: str) -> str:
    return re.sub(r"[^A-Za-z]", "", str(s)).upper().rstrip("*")


def seqmap(code: str, queries: dict[str, str], cache_tsv: Path, min_id: float = 98.0, min_cov: float = 0.8,
           margin: float = 2.0) -> dict[str, dict]:
    """Same-species sequence cross-reference of {qid: seq} onto panel proteome ``code``.

    exact (identical sequence, unique) | unique (best hit >= min_id identity, query and target coverage
    >= min_cov, reciprocal best hit, >= ``margin`` identity points over the runner-up target) | otherwise a
    failure status with the reason. DIAMOND hits are cached in ``cache_tsv`` (deterministic reruns)."""
    prot = proteome(code)
    by_hash: dict[str, list[str]] = {}
    for a, s in prot.items():
        by_hash.setdefault(_hl.sha256(norm_seq(s).encode()).hexdigest(), []).append(a)
    out: dict[str, dict] = {}
    todo = {}
    for q, s in queries.items():
        s = norm_seq(s)
        hits = by_hash.get(_hl.sha256(s.encode()).hexdigest(), [])
        if len(hits) == 1:
            out[q] = dict(status="exact", acc=hits[0], ident=100.0, qcov=1.0, tcov=1.0, margin="", note="identical sequence")
        elif len(hits) > 1:
            out[q] = dict(status="ambiguous", acc="", ident=100.0, qcov=1.0, tcov=1.0, margin=0,
                          note=f"identical to {len(hits)} panel proteins {hits[:4]}")
        else:
            todo[q] = s
    if not todo:
        return out
    cache_tsv = Path(cache_tsv)
    if not cache_tsv.exists():
        import tempfile
        with tempfile.TemporaryDirectory() as td:
            qf = Path(td) / "q.fa"
            qf.write_text("".join(f">{q}\n{s}\n" for q, s in todo.items()))
            tf = ROOT / species_table()[code]["fasta_path"]
            db = Path(td) / "t"
            _sp.run([DIAMOND, "makedb", "--in", str(tf), "-d", str(db), "--quiet"], check=True)
            fmt = ["6", "qseqid", "sseqid", "pident", "length", "qlen", "slen", "qstart", "qend", "sstart", "send",
                   "evalue", "bitscore"]
            _sp.run([DIAMOND, "blastp", "-q", str(qf), "-d", str(db), "-o", str(cache_tsv) + ".tmp", "--very-sensitive",
                     "-k", "10", "--quiet", "--outfmt", *fmt], check=True)
            # reciprocal: best target -> query set
            _os.replace(str(cache_tsv) + ".tmp", cache_tsv)
    hits: dict[str, list] = {}
    with open(cache_tsv) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            q, t = f[0], f[1]
            pid, qlen, slen = float(f[2]), int(f[4]), int(f[5])
            qcov = (int(f[7]) - int(f[6]) + 1) / qlen
            tcov = (int(f[9]) - int(f[8]) + 1) / slen
            hits.setdefault(q, []).append((float(f[11]), pid, t, qcov, tcov))
    # best query per target (reciprocal-best-hit check, among this query set)
    best_q_for_t: dict[str, tuple] = {}
    for q, hs in hits.items():
        if q not in todo:
            continue
        for bs, pid, t, qc, tc in hs:
            if t not in best_q_for_t or bs > best_q_for_t[t][0]:
                best_q_for_t[t] = (bs, q)
    for q in todo:
        hs = sorted({h[2]: h for h in sorted(hits.get(q, []))}.values(), reverse=True)  # best hsp per target
        if not hs:
            out[q] = dict(status="unmapped", acc="", ident="", qcov="", tcov="", margin="", note="no DIAMOND hit")
            continue
        bs, pid, t, qc, tc = hs[0]
        runner = max([h[1] for h in hs[1:] if h[3] >= 0.5] or [0.0])
        mg = round(pid - runner, 1)
        rec = dict(acc=t, ident=pid, qcov=round(qc, 3), tcov=round(tc, 3), margin=mg)
        problems = []
        if pid < min_id:
            problems.append(f"identity {pid}<{min_id}")
        if qc < min_cov or tc < min_cov:
            problems.append(f"coverage q={qc:.2f}/t={tc:.2f}<{min_cov}")
        if mg < margin:
            problems.append(f"margin {mg}<{margin} over runner-up")
        if best_q_for_t.get(t, (0, q))[1] != q:
            problems.append(f"not reciprocal best (target best-matched by {best_q_for_t[t][1]})")
        if problems:
            st = "ambiguous" if (mg < margin or "reciprocal" in " ".join(problems)) and pid >= min_id and qc >= min_cov and tc >= min_cov else "unmapped"
            rec.update(status=st, acc_candidate=t, note="; ".join(problems))
            rec["acc"] = ""
        else:
            rec.update(status="unique", note="RBH, thresholds met")
        out[q] = rec
    return out


_IDX: dict[str, dict] = {}


def id_index(code: str) -> dict[str, set]:
    """Gene-name / locus index of a panel proteome: NAME (upper case) -> set of accessions. Uses UniProt gene
    names (GN=), ordered-locus / ORF names when present in the header, and the gene2acc symbol columns."""
    if code in _IDX:
        return _IDX[code]
    idx: dict[str, set] = {}
    for acc, h in header_map(code).items():
        names = set()
        for k in ("gene_name", "gene_symbol_gene2acc", "canonical_of_gene_group"):
            for v in re.split(r"[;,\s]+", h.get(k) or ""):
                if v:
                    names.add(v.split(":", 1)[-1].upper())
        m = re.search(r"\bGN=(\S+)", h.get("original_header", ""))
        if m:
            names.add(m.group(1).upper())
        for v in names:
            idx.setdefault(v, set()).add(acc)
    _IDX[code] = idx
    return idx


def _norm_text(s: str) -> str:
    import unicodedata
    s = unicodedata.normalize("NFKC", s)
    s = s.replace("–", "-").replace("—", "-").replace("−", "-").replace("’", "'").replace("‘", "'")
    s = s.replace("“", '"').replace("”", '"')
    return re.sub(r"\s+", " ", s).casefold()


def write_assertions(path: Path, rows: list[dict], texts: dict[str, str]) -> None:
    """rows: assertion_id, type, target, tier, source_id, locator, quote. Every quote must occur (normalised) in
    texts[source_id]; raises otherwise."""
    norm = {k: _norm_text(v) for k, v in texts.items()}
    for r in rows:
        q = r["quote"]
        if len(q) > 300:
            raise SystemExit(f"quote too long ({len(q)}): {r['assertion_id']}")
        if _norm_text(q) not in norm[r["source_id"]]:
            raise SystemExit(f"quote not found in {r['source_id']}: {r['assertion_id']}: {q[:80]}")
    cols = ["assertion_id", "type", "target", "tier", "source_id", "locator", "quote"]
    with open(path, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[c]).replace("\t", " ") for c in cols) + "\n")


# --------------------------------------------------------------------------- rice MSU7 -> RAP (IRGSP-1.0) loci
MSU_INFO = ROOT / "literature/raw/P_shared/all.locus_brief_info.7.0"
_MSU: dict[str, tuple] = {}


def msu_locus(locus: str):
    """(chrom number, start, stop, strand) of an MSU7 locus (pseudomolecules = IRGSP-1.0 sequences)."""
    if not _MSU:
        with open(MSU_INFO) as fh:
            next(fh)
            for line in fh:
                f = line.split("\t")
                c, loc, s, e, o = f[0], f[1], int(f[3]), int(f[4]), f[5]
                old = _MSU.get(loc)
                _MSU[loc] = (c.replace("Chr", ""), min(s, old[1]) if old else s, max(e, old[2]) if old else e, o)
    return _MSU.get(locus.split(".")[0])


def msu_to_rap(cache: "JsonCache", locus: str):
    """Unique RAP-DB protein-coding gene overlapping the MSU7 locus on the same strand by >=50 % of the MSU span
    (Ensembl Plants IRGSP-1.0 annotation). Returns (rap_id, note) or ('', reason)."""
    pos = msu_locus(locus)
    if not pos:
        return "", f"{locus} not in MSU7 locus table"
    c, s, e, o = pos
    genes = ensembl_overlap_genes(cache, "oryza_sativa", f"{c}:{s}-{e}")
    strand = 1 if o == "+" else -1
    good = []
    for g in genes:
        ov = min(e, g["end"]) - max(s, g["start"]) + 1
        if g.get("strand") == strand and ov >= 0.5 * (e - s + 1):
            good.append(g["id"])
    if len(good) == 1:
        return good[0], f"{locus} Chr{c}:{s}-{e}{o} -> RAP {good[0]} (unique same-strand overlap >=50 %)"
    return "", f"{locus} Chr{c}:{s}-{e}{o}: {len(good)} RAP genes overlap ({good[:3]})"
