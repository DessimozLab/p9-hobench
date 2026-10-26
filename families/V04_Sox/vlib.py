"""Shared helpers for Panel V family builds (V01-V04). Not a family (no family.yaml).

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
import urllib.parse
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from hobench.newick import Node  # noqa: E402
from hobench.speciestree import SpeciesTree  # noqa: E402

PANEL = ROOT / "panels" / "V"


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
        if not acc and "not in panel proteome" in note:
            cands = re.findall(r"'([A-Z0-9]+)'", note.split("not in panel")[0])
            out["candidates"] = cands
        if not acc and t.startswith("ensembl"):
            sq = ensembl_or_uniparc_seq(cache, oid, t)
            if sq:
                nv = cache.data.get("nver|ens:" + oid)
                keys = [f"ens:{oid}#{j}" for j in range(nv)] if nv else ["ens:" + oid]
                out["candidates"] = out.get("candidates", []) + keys
                note += f"; sequence for {oid} from {sq}"
        if acc:
            meth = rows[0]["mapping_method"] if rows else "ensembl"
            out.update(status="mapped", gene_id=acc, mapping_method="exact_id",
                       mapping_detail=f"UniProt ID mapping ({meth})", reason=note,
                       mapping_identity="100", mapping_coverage="1.0")
        else:
            out.update(status="ambiguous" if note.startswith("ambiguous") else "unmapped", reason=note)
        return out
    if region and code in ENS_SP:
        genes = ensembl_overlap_genes(cache, ENS_SP[code], region)
        if len(genes) == 1:
            g = genes[0]["id"]
            acc, note = pick_in_proteome(code, uniprot_map(cache, [(g, "ensembl_gene")])[g])
            if not acc and "not in panel proteome" in note:
                out["candidates"] = re.findall(r"'([A-Z0-9]+)'", note.split("not in panel")[0])
            if not acc and ensembl_or_uniparc_seq(cache, g, "ensembl_gene"):
                out["candidates"] = out.get("candidates", []) + ["ens:" + g]
                out["mapping_detail"] = f"Ensembl locus {region} -> {g}"
            if acc:
                out.update(status="mapped", gene_id=acc, mapping_method="exact_id",
                           mapping_detail=f"Ensembl locus {region} -> {g} -> UniProt ID mapping",
                           reason=f"Ensembl gene {g}; {note}",
                           mapping_identity="100", mapping_coverage="1.0")
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


def fetch_uniprot_seqs(cache: JsonCache, accs: list[str]) -> dict[str, str]:
    todo = [a for a in accs if "seq|" + a not in cache.data]
    for i in range(0, len(todo), 100):
        chunk = todo[i:i + 100]
        url = ("https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=" +
               urllib.parse.quote(" OR ".join(f"accession:{a}" for a in chunk)))
        txt = None
        for k in range(4):
            try:
                with urllib.request.urlopen(url, timeout=120) as r:
                    txt = r.read().decode()
                break
            except Exception:
                time.sleep(3 + 3 * k)
        got = {}
        for blk in (txt or "").split(">")[1:]:
            h, *seq = blk.split("\n")
            got[h.split("|")[1]] = "".join(seq)
        for a in chunk:
            cache.data["seq|" + a] = got.get(a, "")
        cache.save()
    return {a: cache.data["seq|" + a] for a in accs}


def resolve_by_sequence(cache: JsonCache, code: str, recs: list[dict], strict_family: bool = False) -> None:
    """Same-species sequence fallback for records whose UniProt hits are outside the panel proteome
    (``rec['candidates']``): every candidate sequence is mapped to the panel proteome; the record is mapped
    when all candidates agree on one target with identity >= 98% (99% for strict species), coverage >= 80%
    of query and target, and no other target within 2 identity points."""
    strict = species_table()[code].get("strict_mapping") == "1"
    min_id = 99.0 if (strict or strict_family) else 98.0
    todo = [r for r in recs if r.get("status") != "mapped" and r.get("candidates")]
    if not todo:
        return
    accs = sorted({a for r in todo for a in r["candidates"]})
    seqs = fetch_uniprot_seqs(cache, [a for a in accs if not a.startswith("ens:")])
    seqs.update({a: cache.data.get("seq|" + a, "") for a in accs if a.startswith("ens:")})
    key = "seqmap|" + code + "|" + ",".join(accs)
    def run():
        from curation_tools.map_by_sequence import map_sequences
        q = [(a, s) for a, s in seqs.items() if s]
        path = ROOT / species_table()[code]["fasta_path"]
        return map_sequences(q, path, ambiguity_delta=2.0, min_qcov=0.8, min_tcov=0.8)
    res = {r["query_id"]: r for r in cache.get(key, run)}
    for r in todo:
        hits = [res.get(a) for a in r["candidates"]]
        ok = [h for h in hits if h and h["status"] in ("exact", "unique")
              and float(h["mapping_identity"] or 0) >= min_id and float(h["query_cov"] or 0) >= 0.8
              and float(h["target_cov"] or 0) >= 0.8]
        tg = {h["target_id"] for h in ok}
        if len(tg) == 1:
            h = ok[0]
            r.update(status="mapped", gene_id=h["target_id"], mapping_method="seq_same_species",
                     mapping_detail=(r.get("mapping_detail", "") + "; " if r.get("mapping_detail") else "")
                     + f"same-species {h['method']} {h['status']} (candidate {h['query_id']}, min_id {min_id:g})",
                     mapping_identity=f"{float(h['mapping_identity']):.1f}",
                     mapping_coverage=f"{min(float(h['query_cov']), float(h['target_cov'])):.3f}",
                     reason=r["reason"] + f"; sequence of {h['query_id']} -> {h['target_id']} ({h['status']})")
        else:
            r["reason"] += f"; sequence fallback failed ({len(tg)} distinct qualifying targets)"


def ensembl_or_uniparc_seq(cache: JsonCache, oid: str, t: str) -> str:
    """Protein sequence for an Ensembl ID: current IDs via Ensembl REST (gene -> canonical translation),
    retired protein/transcript IDs via UniParc cross-references (UniParc keeps inactive Ensembl xrefs).
    Stores it under cache['seq|ens:<oid>']; returns the provenance string ('' if none)."""
    k = "seq|ens:" + oid
    if k in cache.data and cache.data[k]:
        return cache.data.get("seqsrc|ens:" + oid, "cache")
    seq, src = "", ""
    lk = ensembl_lookup(cache, oid)
    if lk:
        tid = None
        if lk.get("object_type") == "Gene":
            lk2 = cache.get("lookupx|" + oid, lambda: http_json(f"{ENS}/lookup/id/{oid}?expand=1;content-type=application/json"))
            tid = (lk2 or {}).get("canonical_transcript", "").split(".")[0] or None
            if not tid:
                tr = [x for x in (lk2 or {}).get("Transcript", []) if x.get("is_canonical")]
                tid = tr[0]["id"] if tr else None
        else:
            tid = oid
        if tid:
            r = cache.get("pseq|" + tid, lambda: http_json(f"{ENS}/sequence/id/{tid}?type=protein;content-type=application/json"))
            if isinstance(r, dict) and r.get("seq"):
                seq, src = r["seq"], f"Ensembl REST {tid} (release current)"
    if not seq and t in ("ensembl_protein", "ensembl_transcript"):
        def up():
            url = "https://rest.uniprot.org/uniparc/stream?format=fasta&query=" + urllib.parse.quote(f"dbid:{oid}")
            try:
                with urllib.request.urlopen(url, timeout=120) as r:
                    return r.read().decode()
            except Exception:
                return ""
        txt = cache.get("uniparc|" + oid, up)
        blks = (txt or "").split(">")[1:]
        if blks:
            # several UniParc sequences = successive versions of the retired ID; each becomes a candidate
            for j, b in enumerate(blks):
                h, *sq = b.split("\n")
                cache.data[f"seq|ens:{oid}#{j}"] = "".join(sq)
            cache.data["nver|ens:" + oid] = len(blks)
            seq, src = "".join(blks[0].split("\n")[1:]), "UniParc " + ",".join(b.split()[0] for b in blks)
    cache.data[k] = seq
    cache.data["seqsrc|ens:" + oid] = src
    cache.save()
    return src
