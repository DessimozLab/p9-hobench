#!/usr/bin/env python
"""P04_MIPs: plant major intrinsic proteins / aquaporins (Danielson & Johanson 2008 BMC Plant Biol 8:45;
Anderberg, Kjellbom & Johanson 2012 Front Plant Sci 3:33).

Sources:
  DJ2008 Additional file 1 (1471-2229-8-45-S1.pdf): full-length alignment of all 91 MIPs of A. thaliana (35),
         Z. mays (33) and P. patens (23); names carry the subfamily/group (e.g. AtPIP1.1). Parsed with pypdf.
  DJ2008 Table 1: P. patens MIPs with JGI v1.1 protein ids (-> PHYPADRAFT_<id>), pseudogenes.
  AKJ2012 Table A2: GenBank accessions of the 33 rice MIPs (OsPIP1;1 ...), sequences fetched from NCBI (cached).
Mapping: same-species sequence (exact / >=98 % id, >=80 % cov, RBH, margin 2), else UniProt gene name
(<SUBFAMILY><g>-<n>, e.g. PIP1-1, single hit in ARATH/MAIZE) or PHYPADRAFT_<id>, or GenBank -> UniProt ID mapping.

Tree:
  root U (tag HGT) [ native MIPs, PpGIP1;1 (xenolog) ]
  native = D_MIP_subfamilies (lo Embryophyta, hi ABOVE_RANGE)[PIP, TIP, NIP, SIP, XIP, HIP]
    PIP = U[ PpPIP3;1, D_PIP1_PIP2 (lo Embryophyta, hi ABOVE_RANGE)[PIP1 (U), PIP2 (U)] ]
    TIP = U[ PpTIP6s, D_TIP_groups (lvl Magnoliopsida: groups evolved in the seed-plant ancestor)[TIP1..TIP5] ]
    SIP = S_SIP (Embryophyta)[ PpSIP1s (U), D_SIP1_SIP2 (lo Mesangiospermae, hi Magnoliopsida)[SIP1, SIP2] ]
    NIP: U (At and monocot NIP group nomenclatures differ); XIP, HIP: P. patens only.
"""
from __future__ import annotations

import json
import re
import sys
import urllib.request
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))  # plib.py: copy of families/_shared_P/plib.py (self-contained rebuild)
sys.path.insert(0, str(ROOT))
import plib  # noqa: E402
from plib import Node, mk, nhx, sha256  # noqa: E402

FAM = "P04_MIPs"
DJ, AKJ = "DanielsonJohanson2008", "Anderberg2012"
SUPP = ROOT / "literature/raw/B2/supp/PMC2386804"
ALN = SUPP / "1471-2229-8-45-S1.pdf"
DJ_TXT = ROOT / "literature/papers/PMID18430224/fulltext.txt"
DJ_T1 = ROOT / "literature/papers/PMID18430224/tables/T1.tsv"
AKJ_TXT = ROOT / "literature/papers/PMID22639644/fulltext.txt"
AKJ_TA2 = ROOT / "literature/papers/PMID22639644/tables/TA2.tsv"
CACHE = HERE / "cache"
PREFIX = {"At": "ARATH", "Zm": "MAIZE", "Pp": "PHYPA", "Os": "ORYSJ"}

ASSERTIONS = [
    (AKJ, "Introduction", "event", "D_MIP_subfamilies", "T1",
     "While the moss Physcomitrella patens contains representatives from all seven groups the dicots at most have five "
     "of these subfamilies, the PIPs, TIPs, NIPs, SIPs, and XIPs and monocots only four, having lost also the XIPs"),
    (DJ, "Abstract", "event", "D_MIP_subfamilies", "T1",
     "Our results suggest an early diversification of MIPs resulting in a large number of subfamilies already in "
     "primitive terrestrial plants."),
    (DJ, "Discussion, PIPs - the most conserved MIPs in plants", "event", "D_PIP1_PIP2", "T1",
     "Both PIP1s and PIP2s are highly conserved in P. patens indicating that these groups must have formed early on in "
     "the evolution of land plants"),
    (AKJ, "Abstract", "event", "D_TIP_groups", "T1",
     "Based on ESTs from non-angiosperms we conclude that the specialized groups of TIPs present in higher plants are "
     "not found in primitive vascular plants but evolved later in a common ancestor of seed plants."),
    (AKJ, "Discussion, SIP2s evolved later", "event", "D_SIP1_SIP2", "T1",
     "suggesting that SIP2s have only evolved in the angiosperms whereas SIP1s go deeper and are present in all "
     "terrestrial plants sharing a common ancestor with the mosses"),
    (AKJ, "Discussion, SIP2s evolved later", "event", "S_SIP", "T2",
     "SIP1s go deeper and are present in all terrestrial plants sharing a common ancestor with the mosses"),
    (DJ, "Results, Phylogeny and classification", "tag", "MIP_root_GIP_xenolog", "T1",
     "One of these was the MIP xenolog (homolog resulting from horizontal gene transfer) PpGIP1;1 previously "
     "identified as a GlpF-like MIP and named accordingly"),
    (AKJ, "Introduction", "tag", "D_MIP_subfamilies", "T1",
     "This indicates that the seventh subfamily, the GIPs, was lost early on in the lineage leading to vascular plants, "
     "already before the appearance of the lycophytes."),
    (AKJ, "Discussion, XIPs", "tag", "D_MIP_subfamilies", "T1",
     "Interestingly, the XIPs have recently been lost in both monocots and in Arabidopsis."),
    (DJ, "Abstract", "completeness", "PHYPA", "T1",
     "Analyses of the whole genome of Physcomitrella patens resulted in the identification of 23 MIPs, belonging to "
     "seven different subfamilies"),
]


def alignment() -> dict[str, str]:
    from pypdf import PdfReader
    pat = re.compile(r"(At|Zm|Pp)[A-Z]{3}\d\.\d")
    seqs: dict[str, str] = {}
    for p in PdfReader(str(ALN)).pages:
        lines = [l.strip() for l in p.extract_text().split("\n") if l.strip()]
        names = [l for l in lines if pat.fullmatch(l)]
        res = [l for l in lines if re.fullmatch(r"[-A-Z ]+", l) and not pat.fullmatch(l)]
        if len(names) != len(res):
            raise SystemExit(f"alignment page parse error: {len(names)} names vs {len(res)} rows")
        for n, x in zip(names, res):
            seqs[n] = seqs.get(n, "") + x.replace(" ", "")
    return {n: s.replace("-", "") for n, s in seqs.items()}


def ncbi_fetch(accs: list[str]) -> dict[str, str]:
    cf = CACHE / "ncbi_rice_mips.fa"
    if not cf.exists():
        url = ("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&rettype=fasta&retmode=text&id="
               + ",".join(accs))
        cf.write_text(urllib.request.urlopen(url, timeout=120).read().decode())
    out, cur = {}, None
    for line in cf.read_text().splitlines():
        if line.startswith(">"):
            cur = line[1:].split()[0].split(".")[0]
            out[cur] = ""
        elif cur:
            out[cur] += line.strip()
    return out


def group_of(name: str) -> tuple[str, str]:
    """'AtPIP1.1' / 'OsPIP1;1' -> ('PIP', 'PIP1')."""
    m = re.match(r"(?:At|Zm|Pp|Os)([A-Z]{3})(\d)", name)
    return m.group(1), m.group(1) + m.group(2)


def main() -> None:
    st = plib.species_tree()
    aln = alignment()
    assert len(aln) == 91, len(aln)
    t1 = pd.read_csv(DJ_T1, sep="\t")
    pp_ids = {}
    for _, r in t1.iterrows():
        nm = str(r.iloc[0]).strip()
        if re.fullmatch(r"[A-Z]{3}\d;\d", nm):
            pp_ids["Pp" + nm.replace(";", ".")] = str(r.iloc[4]).strip()
    rice = {}
    for line in AKJ_TA2.read_text().splitlines():
        f = line.split("\t")
        if f and re.fullmatch(r"OsN?[A-Z]{3}\d;\d", f[0]) and len(f) > 2:
            m = re.match(r"([A-Z]{3}\d{5})", f[2].strip())
            if m:
                rice[f[0]] = m.group(1)
    rice_seq = ncbi_fetch(sorted(rice.values()))
    entries = []  # (name, code, seq, idlist, locator, source)
    for n, s in aln.items():
        code = PREFIX[n[:2]]
        ids = []
        sub, grp = group_of(n)
        if code in ("ARATH", "MAIZE"):
            ids = [f"{grp}-{n.split('.')[-1]}"]
        elif code == "PHYPA" and pp_ids.get(n, "").isdigit():
            ids = [f"PHYPADRAFT_{pp_ids[n]}"]
        entries.append((n, code, s, ids, f"Additional file 1 alignment (S1.pdf), row {n}", DJ))
    for n, acc in rice.items():
        entries.append((n, "ORYSJ", rice_seq.get(acc, ""), [acc], f"Table A2 row {n} (GenBank {acc})", AKJ))
    queries: dict[str, dict[str, str]] = {}
    for n, code, s, *_ in entries:
        if s:
            queries.setdefault(code, {})[n] = s
    seqres = {}
    for code, q in sorted(queries.items()):
        seqres.update(plib.seqmap(code, q, CACHE / f"diamond_{code}.tsv"))
    cache = plib.JsonCache(CACHE / "uniprot_idmap.json")
    gb = sorted({(i, "genbank_protein") for n, c, s, ids, *_ in entries if c == "ORYSJ" for i in ids})
    upmap = plib.uniprot_map(cache, gb, taxon="39947") if gb else {}
    report, members, info, by_acc = [], [], {}, {}
    for n, code, s, ids, loc, src in entries:
        sub, grp = group_of(n)
        base = dict(original_id=n, species=code, label=f"{sub}|{grp}", gene_id="", locator=loc, source=src)
        r = seqres.get(n)
        d = dict(base=base, code=code)
        if r and r["status"] in ("exact", "unique"):
            d.update(acc=r["acc"], method=f"sequence_{r['status']}", ident=r["ident"], cov=f"{r['qcov']}/{r['tcov']}",
                     note=r["note"] + (f"; margin {r['margin']}" if r.get("margin") != "" else ""))
        else:
            hits = set()
            for i in ids:
                if code == "ORYSJ":
                    hits |= {x["uniprot_acc"] for x in upmap.get(i, []) if x.get("uniprot_acc") in plib.proteome(code)}
                else:
                    hits |= plib.id_index(code).get(i.upper(), set())
            seqnote = f"sequence: {r['status']} ({r.get('note', '')})" if r else "sequence: none"
            if len(hits) == 1:
                acc = next(iter(hits))
                d.update(acc=acc, method="exact_id", ident="", cov="",
                         note=f"exact identifier mapping {ids} -> {acc} (single panel hit); {seqnote}")
            elif hits:
                d.update(status="ambiguous", reason=f"{ids} -> several panel proteins {sorted(hits)[:4]}; {seqnote}")
            else:
                d.update(status="unmapped", reason=f"{seqnote}; identifiers {ids} give no panel hit")
        info[n] = d
        if d.get("acc"):
            by_acc.setdefault(d["acc"], []).append(n)
    for acc, ns in by_acc.items():
        if len(ns) > 1:
            for x in ns:
                info[x].pop("acc")
                info[x].update(status="ambiguous", reason=f"{len(ns)} paper genes ({', '.join(ns)}) map to one panel protein {acc}")
    for n, d in info.items():
        b = d["base"]
        if d.get("acc"):
            report.append({**b, "gene_id": d["acc"], "status": "mapped", "reason": f"{d['method']}: {d['note']}"})
            members.append(dict(gene_id=d["acc"], species=d["code"], original_id_in_paper=n, label=b["label"],
                                mapping_method=d["method"], mapping_identity=d["ident"], mapping_coverage=d["cov"],
                                source=b["source"], locator=b["locator"]))
        else:
            report.append({**b, "status": d["status"], "reason": d["reason"]})
    # Table 1 pseudogenes / partial genes (not in the alignment)
    for _, r in t1.iterrows():
        nm = str(r.iloc[0]).strip()
        if nm.startswith(("Pseudo", "Partial")) or nm == "-":
            report.append(dict(original_id=nm if nm != "-" else f"ProteinID {r.iloc[4]}", species="PHYPA",
                               label="MIP", gene_id="", status="excluded",
                               reason=f"Table 1: {str(r.iloc[5])[:150]}", locator="DJ2008 Table 1", source=DJ))
    mem = pd.DataFrame(members)

    # ---------------------------------------------------------------- tree
    T = lambda src=DJ: dict(tier="T2", src=src)  # noqa: E731
    by: dict[str, list] = {}
    for m in members:
        sub, grp = m["label"].split("|")
        by.setdefault(grp, []).append((m["gene_id"], m["species"], m["original_id_in_paper"]))

    def u(genes, nid, src=DJ):
        if not genes:
            return None
        if len(genes) == 1:
            return Node(genes[0][0])
        nd = Node(features=dict(ev="U", id=nid, **T(src)))
        for g in genes:
            nd.add_child(Node(g[0]))
        return nd

    pip = mk("U", [u(by.get("PIP3", []), "PIP3"),
                   mk("D", [u(by.get("PIP1", []), "PIP1"), u(by.get("PIP2", []), "PIP2")], id="D_PIP1_PIP2",
                      lvl_lo="Embryophyta", lvl_hi="ABOVE_RANGE", tier="T1", src=DJ)], id="PIP", **T())
    tip = mk("U", [u(by.get("TIP6", []), "TIP6_Pp"),
                   mk("D", [u(by.get(f"TIP{i}", []), f"TIP{i}") for i in range(1, 6)], id="D_TIP_groups",
                      lvl="Magnoliopsida", tier="T1", src=AKJ)], id="TIP", **T())
    sip_pp = [g for g in by.get("SIP1", []) if g[1] == "PHYPA"]
    sip1_ang = [g for g in by.get("SIP1", []) if g[1] != "PHYPA"]
    d_sip = mk("D", [u(sip1_ang, "SIP1_angiosperm"), u(by.get("SIP2", []), "SIP2")], id="D_SIP1_SIP2",
               lvl_lo="Mesangiospermae", lvl_hi="Magnoliopsida", tier="T1", src=AKJ)
    sip = mk("S", [u(sip_pp, "SIP1_Pp"), d_sip], id="S_SIP", tier="T2", src=AKJ)
    nip = u([g for k, v in by.items() if k.startswith("NIP") for g in v], "NIP")
    native = mk("D", [pip, tip, nip, sip, u(by.get("XIP1", []), "XIP"), u(by.get("HIP1", []), "HIP")],
                id="D_MIP_subfamilies", lvl_lo="Embryophyta", lvl_hi="ABOVE_RANGE", tier="T1", src=AKJ)
    gip = by.get("GIP1", [])
    # GIP is a xenolog (HGT, donor outside the panel). An H root would bound the native subfamily
    # duplications at Embryophyta (hobench rule), which the sources do not support -> U root; the HGT is
    # recorded as a tag assertion on the U node and in notes (xenolog pairs are not scored as orthologs/paralogs).
    root = mk("U", [native] + [Node(g[0]) for g in gip], id="MIP_root_GIP_xenolog", tag="HGT", tier="T2", src=DJ)
    tree = nhx(root)
    (HERE / "reference_tree.nhx").write_text(tree + "\n")
    ids_in = set(re.findall(r"id=([^:\]]+)", tree))
    texts = {DJ: DJ_TXT.read_text(), AKJ: AKJ_TXT.read_text()}
    plib.write_assertions(HERE / "assertions.tsv", [
        dict(assertion_id=f"P04_A{i + 1}", type=a[2], target=a[3], tier=a[4], source_id=a[0], locator=a[1], quote=a[5])
        for i, a in enumerate([a for a in ASSERTIONS if a[3] in ids_in or a[2] == "completeness"])], texts)
    mem.to_csv(HERE / "members.tsv", sep="\t", index=False)
    rep = pd.DataFrame(report)[["original_id", "species", "label", "gene_id", "status", "reason", "locator", "source"]]
    rep.to_csv(HERE / "mapping_report.tsv", sep="\t", index=False)
    inp = rep[rep.status.isin(["mapped", "unmapped", "ambiguous"])]
    frac = inp.groupby("species").apply(lambda x: (x.status == "mapped").mean(), include_groups=False)
    n_in, n_map = len(inp), int((inp.status == "mapped").sum())
    complete = sorted(c for c in ("PHYPA",) if frac.get(c, 0) >= 0.9)
    ts = {}
    for ev, tier in re.findall(r"ev=([SDUH])[^\]]*?tier=(T\w+)", tree):
        ts[f"{ev}_{tier}"] = ts.get(f"{ev}_{tier}", 0) + 1
    import yaml
    meta = {
        "id": FAM, "name": "Plant major intrinsic proteins (aquaporins): PIP, TIP, NIP, SIP, XIP, HIP, GIP",
        "panel": "P", "phenomena": ["HGT", "LOSS", "DEEP"],
        "descriptors": {"n_members": int(len(mem)), "n_species": int(mem.species.nunique()),
                        "n_in_panel_paper_members": n_in, "mapped_fraction": round(n_map / n_in, 3),
                        "mapped_fraction_by_species": {k: round(float(v), 3) for k, v in frac.items()},
                        "groups": mem.label.str.split("|").str[1].value_counts().to_dict()},
        "sources": [
            {"id": DJ, "doi": "10.1186/1471-2229-8-45", "pmid": "18430224", "pmcid": "PMC2386804",
             "title": "Unexpected complexity of the aquaporin gene family in the moss Physcomitrella patens",
             "year": 2008, "source_method": "phylogeny",
             "locator": "Additional file 1 (alignment), Table 1, Results 'Phylogeny and classification', Discussion",
             "quote": "Analyses of the whole genome of Physcomitrella patens resulted in the identification of 23 MIPs",
             "files": [{"path": str(p.relative_to(ROOT)), "sha256": sha256(p)} for p in (ALN, DJ_TXT, DJ_T1)]},
            {"id": AKJ, "doi": "10.3389/fpls.2012.00033", "pmid": "22639644", "pmcid": "PMC3355642",
             "title": "Annotation of Selaginella moellendorffii major intrinsic proteins and the evolution of the "
                      "protein family in terrestrial plants",
             "year": 2012, "source_method": "phylogeny", "locator": "Table A2 (rice accessions), Introduction, Discussion",
             "quote": "the specialized groups of TIPs present in higher plants are not found in primitive vascular plants",
             "files": [{"path": str(p.relative_to(ROOT)), "sha256": sha256(p)} for p in (AKJ_TXT, AKJ_TA2, CACHE / "ncbi_rice_mips.fa")]},
        ],
        "root_level": "Embryophyta", "root_single_copy": False, "complete_species": complete,
        "decoys": [], "distinct_from": [], "tier_summary": ts, "protocol_version": 1,
        "notes": ("Root is U: PpGIP1;1 is a GlpF-like MIP xenolog (HGT, donor outside panel P); an H node would force "
                  "the native subfamily D to Embryophyta, so GIP vs native pairs are left unasserted (tag HGT). The six native subfamilies were present in the land-plant ancestor (GIP lost "
                  "before lycophytes, HIP lost in seed plants, XIP lost in monocots and Arabidopsis). TIP groups 1-5 "
                  "evolved in the seed-plant ancestor (panel branch into Magnoliopsida). SIP2 evolved in angiosperms "
                  "(interval Mesangiospermae..Magnoliopsida); P. patens SIP1s vs angiosperm SIP1/SIP2 = S at Embryophyta. "
                  "NIPs: U (A. thaliana NIP1-7 vs monocot NIP1-4 nomenclatures are not comparable). Selaginella MIPs "
                  "(Anderberg 2012) not included: only scaffold coordinates are published, no sequences/ids. Rice "
                  "names follow the monocot scheme; maize and Arabidopsis follow DJ2008."),
    }
    (HERE / "family.yaml").write_text(yaml.safe_dump(meta, sort_keys=False, allow_unicode=True, width=120))
    print(f"{FAM}: {n_map}/{n_in} mapped ({n_map / n_in:.1%}); {frac.round(2).to_dict()}")


if __name__ == "__main__":
    main()
