#!/usr/bin/env python
"""Build the YGOB pillar collection (Y01_YGOB_pillars, collection YGOB) for panel Y.

Sources (all stored locally; sha256 recorded in collection.yaml):
  * YGOB v7 (Aug 2012) Pillars.tab, Wayback copy truncated at 1 MiB
    (literature/raw/B1/ygob/Pillars.TRUNCATED_wayback_20210728.tab): 5,173 complete pillar lines are used,
    the partial 5,174th line is dropped.  Columns: 12 post-WGD species x 2 slots + ancestor + 8 pre-WGD species.
    Slots of a post-WGD species are NOT track-assigned (YGOB README + Byrne & Wolfe 2005).
  * YGOB v1 (Sep 2005, the release described in Byrne & Wolfe 2005) Synteny_<outgroup>.tab files (complete Wayback
    copies): per pre-WGD outgroup locus, the S. cerevisiae / C. glabrata / S. castellii genes on Track A and Track B,
    with a six-character synteny score.  README: same track => orthologs if the score is robust (1 or 0),
    Track A vs Track B => paralogs.  Used to resolve YEAST x CANGA orthology/paralogy inside v7 pillars.
  * UniProt 2026_03 ordered-locus / ORF names per panel proteome (families/Y01_YGOB_pillars/source_cache/).

Encoding (protocol v1 §2-3; SYNTHESIS A2 / ME-2, allopolyploid WGD):
  unit tree = S@Saccharomycetaceae( KLE clade S-polytomy , U( ZT clade S-polytomy , post-WGD part ) )
  post-WGD part: for YEAST/CANGA with robust v1 tracks: D[WGD](trackA S(YEAST,CANGA), trackB S(...)) or
  S(YEAST,CANGA) (same track, single copies); other species with two copies -> per-species D[WGD];
  remaining single copies -> leaves under the U node (pseudo-orthology risk; not asserted).
  WGD D nodes: lvl_lo = PostWGD, lvl_hi = ZT+WGD ancestor  =>  merge in [ZT+WGD, Saccharomycetaceae].
Run from the project root:  python families/Y01_YGOB_pillars/build.py [--species-tree PATH]
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import re
import shutil
import sys
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
import yaml  # noqa: E402

from hobench.speciestree import SpeciesTree  # noqa: E402

FAM = Path(__file__).resolve().parent
FID = "Y01_YGOB_pillars"
COLLECTION = "YGOB"
CACHE = FAM / "source_cache"
YG = ROOT / "literature/raw/B1/ygob"
V7 = YG / "Pillars.TRUNCATED_wayback_20210728.tab"
V7_README = YG / "README"
V1_README = YG / "wayback/v1_README"
V1_SYN = {o: YG / f"wayback/v1_Synteny_{o}.tab" for o in ("Klactis", "Agossypii", "Kwaltii", "Skluyveri")}
PAPER_YGOB = ROOT / "literature/papers/PMID16169922/raw/fulltext.txt"
PAPER_SC07 = ROOT / "literature/papers/PMID17494770/raw/fulltext.txt"
PAPER_MH15 = ROOT / "literature/papers/PMID26252497/fulltext.txt"

# v7 Pillars.tab columns -> (panel code or None, slot); README order
V7_COLS = [
    ("VANPO", 1), ("TETPH", 1), ("HENBL", 1), ("NAUDC", 1), ("NAUCC", 1), ("KAZNA", 1), ("KAZAF", 1),
    ("CANGA", 1), ("SACUV", 1), ("Skud", 1), ("Smik", 1), ("YEAST", 1), ("ANC", 0),
    ("ZYGRC", 0), ("TORDC", 0), ("KLULA", 0), ("ASHGO", 0), ("ERECY", 0), ("Lkluyveri", 0), ("LACTH", 0),
    ("Lwaltii", 0),
    ("YEAST", 2), ("Smik", 2), ("Skud", 2), ("SACUV", 2), ("CANGA", 2), ("KAZAF", 2), ("KAZNA", 2),
    ("NAUCC", 2), ("NAUDC", 2), ("HENBL", 2), ("TETPH", 2), ("VANPO", 2),
]
POST_WGD = ["YEAST", "SACUV", "CANGA", "KAZAF", "KAZNA", "NAUCC", "NAUDC", "HENBL", "TETPH", "VANPO"]
ZT = ["ZYGRC", "TORDC"]
KLE = ["KLULA", "LACTH", "ASHGO", "ERECY"]
UPID = {"YEAST": "UP000002311_559292", "CANGA": "UP000002428_284593", "NAUCC": "UP000001640_27288",
        "NAUDC": "UP000000689_1071378", "KAZAF": "UP000005220_1071382", "VANPO": "UP000000267_436907",
        "TETPH": "UP000005666_1071381", "ZYGRC": "UP000008536_559307", "TORDC": "UP000005627_1076872",
        "LACTH": "UP000002036_559295", "KLULA": "UP000000598_284590", "ASHGO": "UP000000591_284811",
        "ERECY": "UP000006790_931890", "HENBL": "UP000002866_1071380", "KAZNA": "UP000006310_1071383",
        "SACUV": "UP001162085_230603"}

# ---------------------------------------------------------------- quotes (verbatim; whitespace-normalised)
QUOTES = {
    "pillar_homology": ("YGOB2005", "PMID16169922 fulltext, Results 'pillars and tracks' paragraph",
                        "Each slot either is vacant or contains a gene from that genome homologous to the other "
                        "genes in the pillar."),
    "slots": ("YGOB2005", "PMID16169922 fulltext, Results 'pillars and tracks' paragraph",
              "A pillar has a slot for each gene that can be present, i.e., two slots per post-WGD species and one "
              "per pre-WGD species."),
    "slots_untracked": ("YGOB2005", "PMID16169922 fulltext, Results 'pillars and tracks' paragraph",
                        "Thus, in a pillar, the two slots for genes from a post-WGD species are not preassigned to "
                        "particular tracks."),
    "ohnolog_def": ("YGOB2005", "PMID16169922 fulltext, Discussion",
                    "We have suggested (Wolfe 2000) that paralogs arising from a WGD should be called ohnologs"),
    "pre_wgd_tracks": ("YGOB2005", "PMID16169922 fulltext, Results 'pillars and tracks' paragraph",
                       "There is one track for each pre-WGD species and two tracks (track A and track B) for each "
                       "post-WGD species."),
    "track_orthologs": ("YGOBv1README", "literature/raw/B1/ygob/wayback/v1_README, section (4) Synteny_*.tab",
                        "Genes from two post-WGD species that are both on the same Track in this file can be "
                        "considered orthologs if the synteny score is robust (1 or 0) at all the positions of "
                        "interest."),
    "track_paralogs": ("YGOBv1README", "literature/raw/B1/ygob/wayback/v1_README, section (4) Synteny_*.tab",
                       "Similarly genes on Track A in one species and Track B in another are paralogs."),
    "pseudo_orthologs": ("Scannell2007", "PMID17494770 fulltext, Results (single-copy orthologs and paralogs)",
                         "the remaining 44% are paralogs (these genes became duplicated in the WGD, and after "
                         "speciation the two species reciprocally lost different copies)"),
    "single_copy_risk": ("Scannell2007", "PMID17494770 fulltext, Abstract",
                         "Almost half of their surviving single-copy genes are not orthologs but paralogs formed "
                         "by WGD, as would be expected if most gene pairs were resolved independently."),
    "hybrid": ("MH2015", "PMID26252497 fulltext, Abstract",
               "We used phylogenomics to study the ancient genome duplication occurred in the yeast Saccharomyces "
               "cerevisiae lineage and found compelling evidence for the existence of a contemporaneous "
               "interspecies hybridization."),
    "parents": ("MH2015", "PMID26252497 fulltext, Results (parental lineages)",
                "From now on, we consider the ZT cluster as the extant clade closest to one of the parents (parent "
                "A), while the KLE cluster will be considered as the closest to the other parent (parent B)."),
    "kle_split": ("MH2015", "PMID26252497 fulltext, Results",
                  "at the branch preceding the divergence between Saccharomyces and a clade containing the genera "
                  "Kluyveromyces, Lachancea, and Eremothecium"),
}
QUOTE_FILES = {"YGOB2005": PAPER_YGOB, "YGOBv1README": V1_README, "Scannell2007": PAPER_SC07, "MH2015": PAPER_MH15,
               "YGOBv7README": V7_README}


def norm_ws(s: str) -> str:
    return re.sub(r"\s+", " ", s).strip()


def sha256(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def check_quotes() -> None:
    for k, (sid, _loc, q) in QUOTES.items():
        txt = norm_ws(QUOTE_FILES[sid].read_text(errors="replace"))
        if norm_ws(q) not in txt:
            raise SystemExit(f"quote {k!r} not found verbatim in {QUOTE_FILES[sid]}")


# ---------------------------------------------------------------- ID mapping
def ygob_key(code: str, name: str) -> str:
    n = name.strip()
    if code == "VANPO":            # Kpol_1049.17 (YGOB) == Kpol_1049p17 (GenBank locus tag / UniProt ORF name)
        n = n.replace(".", "p")
    return uni_key(n)


def uni_key(tok: str) -> str:
    # A. gossypii: YGOB 'ADR419C' == UniProt OLN 'ADR419C' or locus tag 'AGOS_ADR419C'
    t = tok.strip().upper()
    if t.startswith("AGOS_"):
        t = t[5:]
    return t.replace("_", "")


def load_index(code: str):
    """key -> set(accessions) for ordered-locus/ORF names (primary index) and gene names/synonyms (fallback),
    restricted to accessions present in the panel proteome FASTA."""
    up = UPID[code].split("_")[0]
    fasta = proteome_path(code)
    accs = {l[1:].split()[0] for l in fasta.open() if l.startswith(">")}
    prim, fb = defaultdict(set), defaultdict(set)
    with gzip.open(CACHE / f"uniprot_oln_{up}.tsv.gz", "rt") as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        for r in rd:
            a = r["Entry"]
            if a not in accs:
                continue
            for col, idx in (("Gene Names (ordered locus)", prim), ("Gene Names (ORF)", prim),
                             ("Gene Names (primary)", fb), ("Gene Names (synonym)", fb)):
                for tok in re.split(r"[ ;/,]+", r.get(col) or ""):
                    if tok:
                        idx[uni_key(tok)].add(a)
    return prim, fb, accs


PANEL_FASTA_DIR = None


def proteome_path(code: str) -> Path:
    if PANEL_FASTA_DIR is not None and (PANEL_FASTA_DIR / f"{code}.fa").exists():
        return PANEL_FASTA_DIR / f"{code}.fa"
    return ROOT / "data/proteomes" / UPID[code] / "proteome.fasta"


# ---------------------------------------------------------------- v1 tracks
def load_v1_tracks():
    """Return {frozenset(YEAST+CANGA genes): {'rel': {(g1,g2): 'same'|'diff'}, 'loc': [...]}}, keyed by the
    pre-WGD outgroup gene, plus conflicts.  Only robust positions (score 1 or 0 at Cgl A/B and Scer A/B)."""
    by_out = {}
    for org, p in V1_SYN.items():
        for ln, line in enumerate(p.read_text().splitlines(), 1):
            f = line.split("\t")
            if len(f) < 8:
                continue
            score, out = f[0], f[1]
            cglA, scasA, scerA, scerB, scasB, cglB = f[2:8]
            # score positions: CglA ScasA ScerA ScerB ScasB CglB
            pos = {"CANGA_A": score[0], "YEAST_A": score[2], "YEAST_B": score[3], "CANGA_B": score[5]}
            robust = all(pos[k] in "10" for k in pos)
            genes = {"YEAST_A": scerA, "YEAST_B": scerB, "CANGA_A": cglA, "CANGA_B": cglB}
            genes = {k: v for k, v in genes.items() if v and v != "undef"}
            by_out[(org, out)] = dict(robust=robust, genes=genes, score=score, loc=f"{p.name}:{ln}")
    return by_out


# ---------------------------------------------------------------- main
def main():
    global PANEL_FASTA_DIR
    ap = argparse.ArgumentParser()
    ap.add_argument("--species-tree", default=str(ROOT / "panels/Y/species_tree.nwk"))
    ap.add_argument("--panel-dir", default=str(ROOT / "panels/Y"))
    ap.add_argument("--min-species", type=int, default=4)
    a = ap.parse_args()
    if Path(a.panel_dir, "proteomes").exists():
        PANEL_FASTA_DIR = Path(a.panel_dir, "proteomes")
    st = SpeciesTree.from_newick(a.species_tree)
    check_quotes()
    global SOURCES
    SOURCES = sources_meta()  # computed once; per-unit family.yaml omits the file hashes (see collection.yaml)
    panel = set(st.species)
    L_POST = st.lca(*[s for s in POST_WGD if s in panel])
    L_ZTW = st.lca("ZYGRC", "YEAST")
    L_SAC = st.lca(*[s for s in POST_WGD + ZT + KLE if s in panel])
    L_KLE = st.lca(*[s for s in KLE if s in panel])
    L_ZT = st.lca(*ZT)
    assert st.parent(L_POST) == L_ZTW and st.is_ancestor(L_SAC, L_ZTW, strict=True), (L_POST, L_ZTW, L_SAC)
    # upper bound of the WGD interval: the lineage just below Saccharomycetaceae that contains ZT+WGD
    L_HI = L_ZTW
    while st.parent(L_HI) != L_SAC:
        L_HI = st.parent(L_HI)
    print(f"levels: PostWGD={L_POST} ZT+WGD={L_ZTW} WGD_lvl_hi={L_HI} Saccharomycetaceae={L_SAC} KLE={L_KLE} ZT={L_ZT}")

    # ---- mapping indices
    idx = {c: load_index(c) for c in UPID if c in panel}

    def map_gene(code, name):
        if code == "SACUV":
            return None, "unmapped", ("YGOB Suva_ names refer to the Scannell et al. 2011 S. uvarum assembly; UniProt "
                                      "UP001162085 (GCA_947243795) has no cross-reference and the YGOB AA.fsa "
                                      "archive copy is truncated (no sequence available for same-species mapping)")
        prim, fb, _ = idx[code]
        k = ygob_key(code, name)
        hits = prim.get(k, set())
        meth = "uniprot_locus_name"
        if not hits:
            hits = fb.get(k, set())
            meth = "uniprot_gene_name"
        if len(hits) == 1:
            return next(iter(hits)), meth, ""
        if not hits:
            return None, "unmapped", "locus name not found in UniProt 2026_03 entries of the panel proteome"
        return None, "ambiguous", "locus name matches several accessions: " + ",".join(sorted(hits))

    # ---- v7 pillars
    raw = V7.read_bytes()
    lines = raw.split(b"\n")
    partial = lines[-1]
    lines = lines[:-1] if not raw.endswith(b"\n") else lines
    assert all(len(l.split(b"\t")) == 33 for l in lines), "unexpected column count"
    print(f"v7 pillars: {len(lines)} complete lines; dropped partial line ({len(partial)} bytes)")

    v1 = load_v1_tracks()
    v1_by_outgene = defaultdict(list)
    for (org, out), rec in v1.items():
        v1_by_outgene[out].append(rec)

    report, members_all, assertions, unit_rows = [], [], [], []
    gene_seen = Counter()
    pillars = []
    for ln, l in enumerate(lines, 1):
        cells = [c.decode() for c in l.split(b"\t")]
        entries = []
        for (code, slot), g in zip(V7_COLS, cells):
            if g == "---" or code == "ANC":
                continue
            entries.append((code, slot, g))
        anc = cells[12] if cells[12] != "---" else ""
        pillars.append((ln, anc, entries))
        for code, slot, g in entries:
            gene_seen[(code, g)] += 1

    # pre-pass: accessions hit by >1 YGOB name (identical copies merged into one UniProt entry) -> excluded
    acc_from = defaultdict(set)
    for ln, anc, entries in pillars:
        for code, slot, g in entries:
            if code in panel and code != "SACUV":
                acc0, _, _ = map_gene(code, g)
                if acc0:
                    acc_from[acc0].add((code, g))
    multi = {acc0: names for acc0, names in acc_from.items() if len(names) > 1}

    unit_dir = FAM / "units"
    if unit_dir.exists():
        shutil.rmtree(unit_dir)
    unit_dir.mkdir()
    acc_used = {}
    stats = Counter()
    for ln, anc, entries in pillars:
        uid = f"YGOB_{anc.replace('.', '_')}" if anc else f"YGOB_L{ln:05d}"
        loc = f"{V7.name}:line {ln}"
        mem = []   # (code, slot, ygob_name, acc)
        for code, slot, g in entries:
            row = dict(unit=uid, original_id=g, species=code, slot=slot, locator=loc)
            if code not in panel:
                row.update(status="not_in_panel", gene_id="", reason="species not in panel Y")
            elif gene_seen[(code, g)] > 1:
                row.update(status="excluded", gene_id="", reason="gene occurs in >1 pillar of the file")
            else:
                acc, meth, why = map_gene(code, g)
                if acc and acc in multi:
                    row.update(status="excluded", gene_id="", reason=f"identical/merged copies: UniProt {acc} carries "
                               f"several YGOB loci ({','.join(sorted(n for _, n in multi[acc]))})")
                elif acc:
                    acc_used[acc] = (code, g)
                    row.update(status="mapped", gene_id=acc, reason=meth)
                    mem.append((code, slot, g, acc, meth))
                else:
                    row.update(status=meth, gene_id="", reason=why)
            report.append(row)
        in_panel = [r for r in report if r["unit"] == uid and r["status"] != "not_in_panel"]
        n_map = sum(r["status"] == "mapped" for r in in_panel)
        sp = {m[0] for m in mem}
        if not in_panel:
            continue
        frac = n_map / len(in_panel)
        reason = None
        if len(sp) < a.min_species:
            reason = f"<{a.min_species} panel species mapped"
        elif frac < 0.8:
            reason = f"mapping {n_map}/{len(in_panel)} < 80%"
        if reason:
            stats["unit_skipped:" + reason.split(" ")[0]] += 1
            for r in report:
                if r["unit"] == uid:
                    r["unit_status"] = "not_a_unit: " + reason
            continue
        res = build_unit(uid, ln, anc, mem, st, L_POST, L_HI, L_SAC, v1_by_outgene, entries)
        if res is None:
            stats["unit_skipped:no_asserted_pair"] += 1
            for r in report:
                if r["unit"] == uid:
                    r["unit_status"] = "not_a_unit: no asserted pair"
            continue
        tree, nodes, tags, track_info = res
        for r in report:
            if r["unit"] == uid:
                r["unit_status"] = "unit"
        stats["units"] += 1
        for t in tags:
            stats["tag:" + t] += 1
        stats["track_resolved"] += bool(track_info)
        write_unit(unit_dir / uid, uid, ln, anc, mem, tree, tags, L_SAC, sp, track_info)
        for m in mem:
            members_all.append(dict(gene_id=m[3], species=m[0], original_id_in_paper=m[2], label=uid,
                                    mapping_method=m[4], mapping_identity="", mapping_coverage="", source="YGOBv7",
                                    locator=loc, unit=uid))
        for nd in nodes:
            assertions.append(nd)
            stats[f"node:{nd['type']}:{nd['what']}"] += 1

    # ---- outputs (collection level)
    with (FAM / "mapping_report.tsv").open("w") as fh:
        cols = ["unit", "unit_status", "original_id", "species", "slot", "status", "gene_id", "reason", "locator"]
        w = csv.DictWriter(fh, cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in report:
            r.setdefault("unit_status", "")
            w.writerow(r)
    with (FAM / "members.tsv").open("w") as fh:
        cols = ["gene_id", "species", "original_id_in_paper", "label", "mapping_method", "mapping_identity",
                "mapping_coverage", "source", "locator", "unit"]
        w = csv.DictWriter(fh, cols, delimiter="\t")
        w.writeheader()
        w.writerows(members_all)
    with (FAM / "assertions.tsv").open("w") as fh:
        cols = ["assertion_id", "type", "what", "target", "genes", "tier", "source_id", "locator", "quote_key", "quote"]
        w = csv.DictWriter(fh, cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for k, (sid, lc, q) in QUOTES.items():
            w.writerow(dict(assertion_id=f"RULE_{k}", type="rule", target="collection", genes="", tier="T1",
                            source_id=sid, locator=lc, quote_key=k, quote=q))
        for nd in assertions:
            w.writerow(nd)
    mstat = Counter((r["species"], r["status"]) for r in report)
    write_collection_yaml(stats, L_POST, L_ZTW, L_SAC, mstat, len(lines), partial)
    print(json.dumps(dict(sorted(stats.items())), indent=1))
    per_sp = defaultdict(Counter)
    for r in report:
        if r.get("unit_status") == "unit":
            per_sp[r["species"]][r["status"]] += 1
    print({k: dict(v) for k, v in sorted(per_sp.items())})


def track_resolution(mem, v1_by_outgene):
    """YEAST/CANGA track relations from v1 synteny tables, if consistent with this v7 pillar."""
    yc = {(m[0], m[2]) for m in mem if m[0] in ("YEAST", "CANGA")}
    if not {c for c, _ in yc} >= {"YEAST", "CANGA"}:
        return None
    outs = [m[2] for m in mem if m[0] in ("KLULA", "ASHGO")]
    recs = [r for o in outs for r in v1_by_outgene.get(o, [])]
    good = []
    for r in recs:
        if not r["robust"]:
            continue
        gset = {(k.split("_")[0], v) for k, v in r["genes"].items()}
        if gset != yc:
            continue   # v1 and v7 disagree on YEAST/CANGA membership: do not transfer tracks
        good.append(r)
    if not good:
        return None
    # relation signature: which genes share a track
    sigs = set()
    for r in good:
        A = frozenset(v for k, v in r["genes"].items() if k.endswith("_A"))
        B = frozenset(v for k, v in r["genes"].items() if k.endswith("_B"))
        sigs.add(frozenset([A, B]) - {frozenset()})
    if len(sigs) != 1:
        return {"conflict": True, "locs": [r["loc"] for r in good]}
    tracks = sorted(sorted(t) for t in next(iter(sigs)))
    return {"conflict": False, "tracks": tracks, "locs": [r["loc"] + " score=" + r["score"] for r in good]}


def build_unit(uid, ln, anc, mem, st, L_POST, L_ZTW, L_SAC, v1_by_outgene, entries):
    """Return (newick-ish nested structure, assertion rows, tags, track_info) or None if nothing is asserted."""
    by_sp = defaultdict(list)
    for m in mem:
        by_sp[m[0]].append(m)
    acc = {m[2]: m[3] for m in mem}
    kle = [m[3] for s in KLE for m in by_sp.get(s, [])]
    zt = [m[3] for s in ZT for m in by_sp.get(s, [])]
    nodes, tags = [], set()
    n_assert_pairs = 0
    # WGD interval upper bound: just below the youngest common ancestor of post-WGD with a KLE species present in
    # this unit (= Saccharomycetaceae when KLE is monophyletic); without KLE genes: the panel-wide bound L_ZTW (L_HI)
    hi = L_ZTW
    kle_sp = [s for s in KLE if by_sp.get(s)]
    if kle_sp:
        b = max((st.lca(L_POST, s) for s in kle_sp), key=st.depth_of)
        hi = L_POST
        while st.parent(hi) != b:
            hi = st.parent(hi)
    WGD_ATTR = dict(lvl_lo=L_POST, lvl_hi=hi, tag="WGD", tier="T2", src="YGOBv7")

    def node(ev, kids, nid, **attr):
        return {"ev": ev, "kids": kids, "id": nid, **attr}

    post_children = []
    tr = track_resolution(mem, v1_by_outgene)
    track_info = None
    used = set()
    if tr and not tr.get("conflict"):
        tracks = tr["tracks"]
        track_nodes = []
        for i, t in enumerate(tracks):
            ids = [acc[g] for g in t if g in acc]
            used.update(t)
            if len(ids) == 1:
                track_nodes.append(ids[0])
            elif len(ids) > 1:
                track_nodes.append(node("S", ids, f"{uid}_trk{i+1}", tier="T2", src="YGOBv1"))
                nodes.append(dict(assertion_id=f"{uid}_trk{i+1}", type="event", target=f"{uid}_trk{i+1}",
                                  genes=",".join(ids), tier="T2", source_id="YGOBv1", what="S_same_track",
                                  locator="; ".join(tr["locs"]), quote_key="track_orthologs",
                                  quote=QUOTES["track_orthologs"][2]))
                n_assert_pairs += 1
        if len(track_nodes) == 2:
            single = all(len(t) == 1 for t in tracks)
            gsp = {m[2]: m[0] for m in mem}
            cross_sp = single and len({gsp[g] for t in tracks for g in t}) == 2
            tg = "WGD,LOSS" if cross_sp else "WGD"
            if cross_sp:
                tags.add("LOSS")
            d = node("D", track_nodes, f"{uid}_wgdYC", **{**WGD_ATTR, "tag": tg, "src": "YGOBv1"})
            post_children.append(d)
            genes = [x for t in track_nodes for x in (t["kids"] if isinstance(t, dict) else [t])]
            nodes.append(dict(assertion_id=f"{uid}_wgdYC", type="event", target=f"{uid}_wgdYC",
                              genes=",".join(genes), tier="T2", source_id="YGOBv1",
                              what="D_WGD_pseudoortholog" if cross_sp else "D_WGD_tracks",
                              locator="; ".join(tr["locs"]) + f"; {V7.name}:line {ln}",
                              quote_key="pseudo_orthologs" if cross_sp else "track_paralogs",
                              quote=QUOTES["pseudo_orthologs" if cross_sp else "track_paralogs"][2]))
            n_assert_pairs += 1
            tags.add("WGD")
        elif len(track_nodes) == 1:
            post_children.append(track_nodes[0])
        track_info = tr
    # per-species ohnolog pairs (no cross-species track information)
    for s in POST_WGD:
        ms = [m for m in by_sp.get(s, []) if m[2] not in used]
        if len(ms) == 2:
            nid = f"{uid}_wgd{s}"
            post_children.append(node("D", [m[3] for m in ms], nid, **WGD_ATTR))
            nodes.append(dict(assertion_id=nid, type="event", target=nid, genes=",".join(m[3] for m in ms),
                              tier="T2", source_id="YGOBv7", what="D_WGD_ohnolog_pair",
                              locator=f"{V7.name}:line {ln} (both {s} slots filled)", quote_key="slots",
                              quote=QUOTES["slots"][2]))
            n_assert_pairs += 1
            tags.add("WGD")
        else:
            post_children.extend(m[3] for m in ms)
    if any(isinstance(c, dict) and c["ev"] == "D" for c in post_children) and \
            any(len(by_sp.get(s, [])) == 1 for s in POST_WGD):
        tags.add("LOSS")
    zt_node = None
    if len(zt) > 1:
        zt_node = node("S", zt, f"{uid}_ZT", tier="T2", src="YGOBv7")
        nodes.append(dict(assertion_id=f"{uid}_ZT", type="event", target=f"{uid}_ZT", genes=",".join(zt),
                          tier="T2", source_id="YGOBv7", what="S_preWGD", locator=f"{V7.name}:line {ln}",
                          quote_key="pillar_homology", quote=QUOTES["pillar_homology"][2]))
        n_assert_pairs += 1
    elif zt:
        zt_node = zt[0]
    ztw_kids = ([zt_node] if zt_node else []) + post_children
    if len(ztw_kids) > 1:
        ztw = node("U", ztw_kids, f"{uid}_ZTWGD", tier="T2", src="YGOBv7")
    elif ztw_kids:
        ztw = ztw_kids[0]
    else:
        ztw = None
    kle_node = None
    if len(kle) > 1:
        kle_node = node("S", kle, f"{uid}_KLE", tier="T2", src="YGOBv7")
        nodes.append(dict(assertion_id=f"{uid}_KLE", type="event", target=f"{uid}_KLE", genes=",".join(kle),
                          tier="T2", source_id="YGOBv7", what="S_preWGD", locator=f"{V7.name}:line {ln}",
                          quote_key="pillar_homology", quote=QUOTES["pillar_homology"][2]))
        n_assert_pairs += 1
    elif kle:
        kle_node = kle[0]
    if isinstance(kle_node, dict) and st.is_ancestor(st.lca(*{m[0] for s in KLE for m in by_sp.get(s, [])}), L_ZTW):
        # KLE not monophyletic in the panel tree: KLE genes become direct children of the Saccharomycetaceae S node
        kle_node = None
        nodes.pop()  # drop the KLE S-polytomy assertion; the pairs are asserted by the root S node instead
        kle_kids = kle
    else:
        kle_kids = [kle_node] if kle_node else []
    if kle_kids and ztw:
        root = node("S", kle_kids + [ztw], f"{uid}_root", tier="T2", src="YGOBv7")
        nodes.append(dict(assertion_id=f"{uid}_root", type="event", target=f"{uid}_root",
                          genes="KLE x (ZT,post-WGD)", tier="T2", source_id="YGOBv7", what="S_Saccharomycetaceae",
                          locator=f"{V7.name}:line {ln}", quote_key="pre_wgd_tracks",
                          quote=QUOTES["pre_wgd_tracks"][2]))
        n_assert_pairs += 1
    else:
        root = kle_node or ztw or (node("S", kle_kids, f"{uid}_root", tier="T2", src="YGOBv7") if len(kle_kids) > 1 else None)
    if n_assert_pairs == 0 or not isinstance(root, dict):
        return None
    return root, nodes, sorted(tags), track_info


def to_nhx(n) -> str:
    if isinstance(n, str):
        return n
    kids = ",".join(to_nhx(k) for k in n["kids"])
    feats = [f"ev={n['ev']}", f"id={n['id']}"]
    for k in ("lvl_lo", "lvl_hi", "tag", "tier", "src"):
        if n.get(k):
            feats.append(f"{k}={n[k]}")
    return f"({kids})[&&NHX:{':'.join(feats)}]"


SOURCES = None


def sources_meta():
    return [
        dict(id="YGOBv7", pmid=16169922, doi="10.1101/gr.3672305",
             title="YGOB v7 (Aug 2012) Pillars.tab (Wayback 2021-07-28, truncated at 1 MiB)", year=2012,
             source_method="synteny", locator="Pillars.tab line = pillar",
             files=[{"path": str(V7.relative_to(ROOT)), "sha256": sha256(V7)},
                    {"path": str(V7_README.relative_to(ROOT)), "sha256": sha256(V7_README)}],
             quote=QUOTES["pillar_homology"][2]),
        dict(id="YGOBv1", pmid=16169922, doi="10.1101/gr.3672305",
             title="YGOB v1 (Sep 2005) Synteny_*.tab track tables (Wayback)", year=2005, source_method="synteny",
             locator="Synteny_<outgroup>.tab line = outgroup locus",
             files=[{"path": str(p.relative_to(ROOT)), "sha256": sha256(p)} for p in [V1_README, *V1_SYN.values()]],
             quote=QUOTES["track_orthologs"][2]),
        dict(id="YGOBv1README", pmid=16169922, doi="10.1101/gr.3672305", title="YGOB v1 README", year=2005,
             source_method="synteny", locator="README section (4)",
             files=[{"path": str(V1_README.relative_to(ROOT)), "sha256": sha256(V1_README)}],
             quote=QUOTES["track_paralogs"][2]),
        dict(id="YGOB2005", pmid=16169922, doi="10.1101/gr.3672305",
             title="The Yeast Gene Order Browser: combining curated homology and syntenic context reveals gene fate "
                   "in polyploid species", year=2005, source_method="synteny", locator="Results; Discussion",
             files=[{"path": str(PAPER_YGOB.relative_to(ROOT)), "sha256": sha256(PAPER_YGOB)}],
             quote=QUOTES["slots_untracked"][2]),
        dict(id="Scannell2007", pmid=17494770, doi="10.1073/pnas.0608218104",
             title="Independent sorting-out of thousands of duplicated gene pairs in two yeast species descended "
                   "from a whole-genome duplication", year=2007, source_method="synteny", locator="Abstract; Results",
             files=[{"path": str(PAPER_SC07.relative_to(ROOT)), "sha256": sha256(PAPER_SC07)}],
             quote=QUOTES["pseudo_orthologs"][2]),
        dict(id="MH2015", pmid=26252497, doi="10.1371/journal.pbio.1002220",
             title="Beyond the Whole-Genome Duplication: Phylogenetic Evidence for an Ancient Interspecies "
                   "Hybridization in the Baker's Yeast Lineage", year=2015, source_method="phylogeny",
             locator="Abstract; Results", files=[{"path": str(PAPER_MH15.relative_to(ROOT)),
                                                  "sha256": sha256(PAPER_MH15)}],
             quote=QUOTES["hybrid"][2]),
    ]


def write_unit(d: Path, uid, ln, anc, mem, tree, tags, L_SAC, sp, track_info):
    d.mkdir(parents=True)
    has_kle = any(m[0] in KLE for m in mem)
    meta = dict(id=uid, name=f"YGOB pillar {anc or 'line ' + str(ln)}", panel="Y", collection=COLLECTION,
                phenomena=sorted(set(tags)) + (["ALLO"] if "WGD" in tags else []),
                sources=[{k: v for k, v in s.items() if k != "files"} for s in SOURCES],
                root_level=L_SAC if has_kle else None, root_single_copy=False, complete_species=[],
                protocol_version=1,
                notes=(f"Pillar {V7.name} line {ln} (ancestral locus {anc or 'none'}). "
                       + ("YEAST/CANGA tracks from YGOB v1 synteny tables: " + "; ".join(track_info["locs"])
                          if track_info and not track_info.get("conflict") else
                          "No robust consistent YGOB v1 track data; cross-species post-WGD pairs unasserted (U).")))
    if meta["root_level"] is None:
        del meta["root_level"]
    (d / "family.yaml").write_text(yaml.safe_dump(meta, sort_keys=False, allow_unicode=True))
    with (d / "members.tsv").open("w") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["gene_id", "species", "original_id_in_paper", "label", "mapping_method", "mapping_identity",
                    "source", "locator"])
        for m in mem:
            w.writerow([m[3], m[0], m[2], uid, m[4], "", "YGOBv7", f"{V7.name}:line {ln}"])
    (d / "reference_tree.nhx").write_text(to_nhx(tree) + ";\n")


def write_collection_yaml(stats, L_POST, L_ZTW, L_SAC, mstat, nlines, partial):
    meta = dict(
        id=FID, name="YGOB pillars (yeast WGD ohnologs and pre-/post-WGD orthologs)", panel="Y",
        collection=COLLECTION, layout="one unit (= one YGOB pillar) per directory under units/; each unit is a "
                                      "loadable family (hobench.reference.load_family)",
        phenomena=["WGD", "LOSS", "ALLO"],
        levels=dict(PostWGD=L_POST, ZT_WGD=L_ZTW, Saccharomycetaceae=L_SAC),
        sources=sources_meta(),
        truncation=dict(file=str(V7.relative_to(ROOT)), sha256=sha256(V7), complete_lines=nlines,
                        dropped_partial_line=partial.decode(errors="replace")),
        stats=dict(sorted(stats.items())),
        mapping=sorted([f"{s}:{st}={n}" for (s, st), n in mstat.items()]),
        complete_species=[], protocol_version=1,
        notes="Encoding per SYNTHESIS A2/ME-2 (allopolyploid WGD): ohnologs D[lvl_lo=PostWGD,lvl_hi=ZT+WGD]; "
              "KLE vs post-WGD orthologs at Saccharomycetaceae (root S); ZT vs post-WGD unasserted (U); "
              "cross-species post-WGD pairs asserted only where YGOB v1 tracks resolve YEAST/CANGA.")
    (FAM / "collection.yaml").write_text(yaml.safe_dump(meta, sort_keys=False, allow_unicode=True))


if __name__ == "__main__":
    main()
