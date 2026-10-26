#!/usr/bin/env python
"""P03_GID1_CXE: gibberellin receptor GID1 nested inside the plant carboxylesterase (CXE) family
(Yoshida et al. 2018 PNAS 115:E7844).

Sources: Dataset S2 (sd02.xlsx: species, name, locus, source; 260 CXE/GID1 genes), Dataset S1 (sd01: alignment of
CXEs/GID1s of At, Os, Sm, Pp) and Dataset S3 (sd03: alignment of 169 GID1s); ungapped alignment rows are used as
query sequences for same-species sequence mapping. Identifiers: AGI, LOC_Os (MSU7 -> RAP-DB locus overlap),
SELMODRAFT_<JGI id>, Phvul.* -> PHAVU_*g, XP_ accessions via UniProt ID mapping.

Tree:
  root U[ Physcomitrella CXEs/GID1-likes,
          D_GID1_origin (lo=Tracheophyta, hi=ABOVE_RANGE; GID1 arose after the moss/vascular split)
            [ U(tracheophyte CXEs incl. GID1-like CXEs), GID1 ] ]
  GID1 = U[ lycophyte (L-type) GID1s, S_Magnoliopsida[ Amborella, S_Mesangiospermae[ monocot M-type (U),
             D_GID1_AB (lvl=Rosids: A/B split just before/after core eudicots)[ A-type (U), B-type (U) ] ] ] ]
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))  # plib.py: copy of families/_shared_P/plib.py (self-contained rebuild)
sys.path.insert(0, str(ROOT))
import plib  # noqa: E402
from plib import Node, mk, nhx, sha256  # noqa: E402

FAM = "P03_GID1_CXE"
SRC = "Yoshida2018"
SUPP = ROOT / "literature/raw/B2/supp/PMC6099883"
S1, S2, S3 = (SUPP / f"pnas.1806040115.sd0{i}.xlsx" for i in (1, 2, 3))
FULLTEXT = ROOT / "literature/papers/PMID30068603/fulltext.txt"
CACHE = HERE / "cache"
SPECIES = {"Arabidopsis thaliana": "ARATH", "Oryza sativa": "ORYSJ", "Selaginella moellendorff": "SELMO",
           "Physcomitrella patens": "PHYPA", "Glycine max": "SOYBN", "Phaseolus vulgaris": "PHAVU",
           "Vitis vinifera": "VITVI", "Amborella trichopoda": "AMBTC", "Zea mays": "MAIZE", "Brassica rapa": "BRARP"}

ASSERTIONS = [
    ("Introduction", "event", "D_GID1_origin", "T1",
     "Previous studies have indicated that GID1 appeared after the divergence of vascular plants from the moss lineage "
     "as no GID1 homologs are found in Physcomitrella mosses or the liverwort Marchantia polymorpha"),
    ("Results, Establishment of GID1 from CXEs", "event", "D_GID1_origin", "T1",
     "confirming that GID1 was derived from one specific CXE group"),
    ("Discussion", "event", "D_GID1_AB", "T1",
     "Core eudicots, with a few exceptions, generally contain a set of A- and B-type GIDs, indicating that "
     "diversification of GID1 into A- and B-types occurred just before or after the establishment of core eudicots."),
    ("Results, Diversification of GID1 in Angiosperms", "event", "S_GID1_Magnoliopsida", "T2",
     "all gymnosperms studied and Amborella trichopoda, the earliest angiosperm, have one gene copy"),
    ("Results, Diversification of GID1 in Angiosperms", "event", "S_GID1_Mesangiospermae", "T2",
     "Thus, the default copy number in the early stage of seed plants might have been one."),
]


def alignment_seqs(path) -> dict[str, str]:
    x = pd.read_excel(path, header=None)
    out: dict[str, str] = {}
    for row in x.itertuples(index=False):
        name = row[0]
        if not isinstance(name, str) or name.startswith("Dataset"):
            continue
        res = "".join(str(v) for v in row[1:] if isinstance(v, str) and re.fullmatch(r"[A-Za-z]", v))
        out[name.strip()] = out.get(name.strip(), "") + res.upper()
    return out


def gtype(name: str, code: str) -> str:
    m = re.search(r"_(GID1[a-z]?)", name)
    if not m and re.search(r"GID1L", name):
        return "CXE_GID1-like"
    if not m:
        return "CXE"
    t = m.group(1)
    if code == "SELMO":
        return "GID1_L"
    if code == "AMBTC":
        return "GID1_BA"
    if code in ("ORYSJ", "MAIZE"):
        return "GID1_M"
    if t.startswith("GID1a") or t.startswith("GID1c"):
        return "GID1_A"
    if t.startswith("GID1b"):
        return "GID1_B"
    return "GID1_unclassified"


def main() -> None:
    st = plib.species_tree()
    t = pd.read_excel(S2, header=1)
    t["xls_row"] = range(3, len(t) + 3)
    t = t[t["Name in this study"].notna()]
    t["species"] = t["species"].ffill().astype(str).str.strip()
    seqs = alignment_seqs(S1)
    seqs.update({k.replace("_", " ", 1) if False else k: v for k, v in alignment_seqs(S3).items()})
    report, members, info, queries = [], [], {}, {}
    for _, r in t.iterrows():
        name = str(r["Name in this study"]).strip()
        code = next((c for k, c in SPECIES.items() if r["species"].startswith(k)), None)
        locus = str(r["Locus name"]).strip() if pd.notna(r["Locus name"]) else ""
        loc = f"Dataset S2 row {r['xls_row']}"
        base = dict(original_id=name, species=code or r["species"], label="", gene_id="", locator=loc, source=SRC)
        if code is None:
            report.append({**base, "label": "GID1" if "GID1" in name else "CXE", "status": "not_in_panel",
                           "reason": "species not in panel P"})
            continue
        lab = gtype(name, code)
        base["label"] = lab
        info[name] = dict(code=code, locus=locus, lab=lab, base=base)
        key = name.replace(" ", "_")
        s = seqs.get(name) or seqs.get(key)
        if s:
            queries.setdefault(code, {})[name] = s
    seqres = {}
    for code, q in sorted(queries.items()):
        seqres.update(plib.seqmap(code, q, CACHE / f"diamond_{code}.tsv"))
    cache = plib.JsonCache(CACHE / "uniprot_idmap.json")
    ens = plib.JsonCache(CACHE / "ensembl_overlap.json")
    xp = sorted({(d["locus"], "refseq_protein") for d in info.values() if re.match(r"[XN]P_\d+", d["locus"])})
    upmap = plib.uniprot_map(cache, xp) if xp else {}
    by_acc: dict[str, list[str]] = {}
    for name, d in info.items():
        code, s, locus = d["code"], seqres.get(name), d["locus"]
        if s and s["status"] in ("exact", "unique"):
            d.update(acc=s["acc"], method=f"sequence_{s['status']}", ident=s["ident"], cov=f"{s['qcov']}/{s['tcov']}",
                     note=s["note"] + (f"; margin {s['margin']}" if s.get("margin") != "" else ""))
        else:
            ids, extra = [], ""
            if code == "ARATH":
                ids = [locus.upper()]
            elif code == "SELMO" and re.fullmatch(r"\d+", locus):
                ids = [f"SELMODRAFT_{locus}"]
            elif code == "PHAVU" and locus.startswith("Phvul."):
                ids = [locus.replace("Phvul.", "PHAVU_") + "g"]
            elif code == "ORYSJ" and locus.startswith("LOC_Os"):
                rap, extra = plib.msu_to_rap(ens, locus.split(".")[0])
                ids = [rap] if rap else []
            hits = set()
            for i in ids:
                hits |= plib.id_index(code).get(i.upper(), set())
            if re.match(r"[XN]P_\d+", locus):
                hits |= {x["uniprot_acc"] for x in upmap.get(locus, []) if x.get("uniprot_acc") in plib.proteome(code)}
                ids.append(locus)
            seqnote = f"sequence: {s['status']} ({s.get('note', '')})" if s else "sequence: not in Datasets S1/S3"
            if len(hits) == 1:
                acc = next(iter(hits))
                d.update(acc=acc, method="exact_id", ident="", cov="",
                         note=f"exact identifier mapping {locus} via {ids} -> {acc} (single panel hit) {extra}; {seqnote}")
            elif hits:
                d.update(status="ambiguous", reason=f"{locus} -> several panel proteins {sorted(hits)[:4]}; {seqnote}")
            else:
                d.update(status="unmapped", reason=f"{seqnote}; locus '{locus}' ({ids}) gives no panel hit {extra}")
        if d.get("acc"):
            by_acc.setdefault(d["acc"], []).append(name)
    for acc, ns in by_acc.items():
        if len(ns) > 1:
            for n in ns:
                info[n].pop("acc")
                info[n].update(status="ambiguous", reason=f"{len(ns)} paper genes ({', '.join(ns)}) map to one panel protein {acc}")
    for name, d in info.items():
        b = d["base"]
        if d.get("acc"):
            report.append({**b, "gene_id": d["acc"], "status": "mapped", "reason": f"{d['method']}: {d['note']}"})
            members.append(dict(gene_id=d["acc"], species=d["code"], original_id_in_paper=name, label=d["lab"],
                                mapping_method=d["method"], mapping_identity=d["ident"], mapping_coverage=d["cov"],
                                source=SRC, locator=b["locator"]))
        else:
            report.append({**b, "status": d["status"], "reason": d["reason"]})
    mem = pd.DataFrame(members)

    T = dict(tier="T2", src=SRC)
    g = {}
    for m in members:
        g.setdefault(m["label"], []).append((m["gene_id"], m["species"]))

    def u(genes, nid):
        if not genes:
            return None
        if len(genes) == 1:
            return Node(genes[0][0])
        n = Node(features=dict(ev="U", id=nid, **T))
        for x in genes:
            n.add_child(Node(x[0]))
        return n

    d_ab = mk("D", [u(g.get("GID1_A", []), "GID1_A"), u(g.get("GID1_B", []), "GID1_B")], id="D_GID1_AB",
              lvl="Rosids", tier="T1", src=SRC)
    s_mes = mk("S", [u(g.get("GID1_M", []), "GID1_M"), d_ab], id="S_GID1_Mesangiospermae", **T)
    s_mag = mk("S", [u(g.get("GID1_BA", []), "GID1_BA"), s_mes], id="S_GID1_Magnoliopsida", **T)
    gid1 = mk("U", [Node(x[0]) for x in g.get("GID1_L", [])] + [s_mag], id="GID1", **T)
    cxe_all = g.get("CXE", []) + g.get("CXE_GID1-like", [])
    trach = [x for x in cxe_all if x[1] != "PHYPA"]
    moss = [x for x in cxe_all if x[1] == "PHYPA"]
    d_orig = mk("D", [u(trach, "CXE_Tracheophyta"), gid1], id="D_GID1_origin", lvl_lo="Tracheophyta",
                lvl_hi="ABOVE_RANGE", tier="T1", src=SRC)
    root = mk("U", [Node(x[0]) for x in moss] + [d_orig], id="CXE_root", **T)
    (HERE / "reference_tree.nhx").write_text(nhx(root) + "\n")
    ids_in_tree = set(re.findall(r"id=([^:\]]+)", nhx(root)))
    plib.write_assertions(HERE / "assertions.tsv", [
        dict(assertion_id=f"P03_A{i + 1}", type=a[1], target=a[2], tier=a[3], source_id=SRC, locator=a[0], quote=a[4])
        for i, a in enumerate([a for a in ASSERTIONS if a[2] in ids_in_tree])], {SRC: FULLTEXT.read_text()})
    mem.to_csv(HERE / "members.tsv", sep="\t", index=False)
    rep = pd.DataFrame(report)[["original_id", "species", "label", "gene_id", "status", "reason", "locator", "source"]]
    rep.to_csv(HERE / "mapping_report.tsv", sep="\t", index=False)
    inp = rep[rep.status.isin(["mapped", "unmapped", "ambiguous"])]
    frac = inp.groupby("species").apply(lambda x: (x.status == "mapped").mean(), include_groups=False)
    n_in, n_map = len(inp), int((inp.status == "mapped").sum())
    import yaml
    txt = (HERE / "reference_tree.nhx").read_text()
    ts = {}
    for ev, tier in re.findall(r"ev=([SDUH])[^\]]*?tier=(T\w+)", txt):
        ts[f"{ev}_{tier}"] = ts.get(f"{ev}_{tier}", 0) + 1
    meta = {
        "id": FAM, "name": "Gibberellin receptor GID1 nested within the plant carboxylesterase (CXE) family",
        "panel": "P", "phenomena": ["NESTED", "NEOF"],
        "descriptors": {"n_members": int(len(mem)), "n_species": int(mem.species.nunique()),
                        "n_in_panel_paper_members": n_in, "mapped_fraction": round(n_map / n_in, 3),
                        "mapped_fraction_by_species": {k: round(float(v), 3) for k, v in frac.items()},
                        "labels": mem.label.value_counts().to_dict()},
        "sources": [{
            "id": SRC, "doi": "10.1073/pnas.1806040115", "pmid": "30068603", "pmcid": "PMC6099883",
            "title": "Evolution and diversification of the plant gibberellin receptor GID1", "year": 2018,
            "source_method": "phylogeny",
            "locator": "Dataset S2 (gene list), Datasets S1/S3 (alignments), Fig. 3, Results, Discussion",
            "quote": "confirming that GID1 was derived from one specific CXE group",
            "files": [{"path": str(p.relative_to(ROOT)), "sha256": sha256(p)} for p in (S1, S2, S3, FULLTEXT, plib.MSU_INFO)]}],
        "root_level": "Embryophyta", "root_single_copy": False, "complete_species": [],
        "decoys": [], "distinct_from": [], "tier_summary": ts, "protocol_version": 1,
        "notes": ("CXE members are those listed in Dataset S2 (Marshall et al. 2003 sets for At/Os/Sm/Pp, partial; "
                  "not a complete CXE inventory -> no complete_species). GID1 types (A/B eudicot, M monocot, BA "
                  "Amborella, L lycophyte) follow the gene names/Fig. 3. D_GID1_origin separates GID1 from the "
                  "tracheophyte CXEs listed (lo=Tracheophyta, hi=ABOVE_RANGE: the progenitor CXE lineage is not "
                  "identified, so only 'separate at Tracheophyta and below' is asserted); moss CXEs and the lycophyte "
                  "GID1s vs seed-plant GID1s are U. D_GID1_AB at the panel branch into Rosids (no basal eudicot in "
                  "panel P; the gamma hexaploidy is on the same branch but the paper does not attribute the split to it, "
                  "so no WGD tag). Within types: U (Glycine/Phaseolus/Brassicaceae copies; WGD lineages)."),
    }
    (HERE / "family.yaml").write_text(yaml.safe_dump(meta, sort_keys=False, allow_unicode=True, width=120))
    print(f"{FAM}: {n_map}/{n_in} mapped ({n_map / n_in:.1%}); {frac.round(2).to_dict()}; labels {meta['descriptors']['labels']}")


if __name__ == "__main__":
    main()
