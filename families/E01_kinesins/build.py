#!/usr/bin/env python
"""E01_kinesins: build members.tsv, mapping_report.tsv, reference_tree.nhx, assertions.tsv, family.yaml.

Source: Wickstead, Gull & Richards 2010 BMC Evol Biol 10:110 (PMID 20423470), Additional file 3 (S3.XLS: 1,624
kinesins with Family/Subfamily labels) and Additional file 8 (DatasetS1.fas: motor-domain alignment, used only to
cross-check ID mappings). Network/DIAMOND results are cached in cache/ by prep_mapping.py; this script is offline.
Run from the project root:  python families/E01_kinesins/build.py
"""
import csv, hashlib, json, re, sys
from collections import defaultdict
from pathlib import Path
import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[2]
FAM = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(FAM))
from prep_mapping import SPMAP, load_table, load_motor, motor_key  # noqa: E402
from hobench.newick import Node, write_newick  # noqa: E402

SRC = "Wickstead2010"
PAPER = ROOT / "literature/papers/PMID20423470"
XLS_REL = "literature/papers/PMID20423470/raw/supp/1471-2148-10-110-S3.XLS"
ALN_REL = "literature/papers/PMID20423470/raw/supp/1471-2148-10-110-S8.ZIP"
TXT_REL = "literature/papers/PMID20423470/fulltext.txt"
PANEL = "E"
# Minimal ancestral repertoire (MAR): families present in the LCEA under all 4 tested eukaryote topologies.
MAR = ["Kinesin-1", "Kinesin-2", "Kinesin-3", "Kinesin-4/10", "Kinesin-5", "Kinesin-8", "Kinesin-9", "Kinesin-13",
       "Kinesin-14", "Kinesin-17"]
Q_MAR = ("These include 11 kinesin paralogs - namely, Kinesin-1, 2, 3, 4/10, 5, 8, 9A, 9B, 13, 14, 17 - and 3 protein "
         "domain architectures - KISc, KISc-FHA, SAM-KISc.")
Q_MAR2 = ("our analysis identifies a complex core set of ancestral characters that were present in the LCEA under any "
          "of these 4 tree topologies.")
Q_SUBFAM = ("we have identified well-supported subfamilies by appending a letter to the family name (e.g. Kinesin-9A "
            "and 9B)")
Q_X = ("We designated these additional tentative paralog families X1-X14")
Q_COMPLETE = ("From these datasets, we extracted complete kinesin repertoires using HMMERv2.3.2 [22] to find all "
              "predicted proteins with a match to the Pfam 'kinesin motor domain' profile (PF00225; [69]).")
ID_MIN, QCOV_MIN, MARGIN = 98.0, 0.80, 2.0


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def clean(v):
    return v.strip() if isinstance(v, str) and v.strip() and v.strip().lower() != "nan" else ""


def panel_sets():
    sp = pd.read_csv(ROOT / f"panels/{PANEL}/species.tsv", sep="\t", dtype=str).set_index("code")
    accs, names = {}, {}
    for code in set(SPMAP.values()):
        fa = ROOT / sp.loc[code, "fasta_path"]
        hm = fa.parent / "header_map.tsv"
        accs[code] = {l[1:].split()[0] for l in open(fa) if l.startswith(">")}
        d = defaultdict(set)
        h = pd.read_csv(hm, sep="\t", dtype=str, keep_default_na=False, usecols=["accession", "gene_name", "gene_symbol_gene2acc"])
        for a, g1, g2 in h.itertuples(index=False):
            if a not in accs[code]:
                continue
            for g in (g1, g2):
                for t in re.split(r"[;,/ ]+", g.split(":")[-1] if ":" in g else g):
                    if t:
                        d[t.upper()].add(a)
        names[code] = d
    return sp, accs, names


def name_candidates(sp5, gid):
    gid = clean(gid)
    if not gid:
        return []
    c = {"Nemve": [f"NEMVEDRAFT_v1g{gid}"], "Triad": [f"TRIADDRAFT_{gid}"], "Monbr": [f"MONBRDRAFT_{gid}"],
         "Naegr": [f"NAEGRDRAFT_{gid}"], "Thaps": [f"THAPSDRAFT_{gid}", f"THAPS_{gid}"], "Cyame": [f"CYME_{gid}"],
         "Neucr": [gid], "Arath": [gid.split(".")[0]], "Trybr": [gid], "Sacce": [gid], "Enthi": [],
         "Giala": [gid] + ([f"GL50803_{gid.split('_')[1].zfill(n)}" for n in (7, 8)] if gid.startswith("GL50803_") else []),
         "Drome": [gid.split("-")[0]], "Plafa": [gid], "Dicdi": [], "Capsp": [], "Chlre": [], "Batde": [],
         "Homsa": [], "Enccu": []}.get(sp5, [])
    return [x.upper() for x in c]


def main():
    sp, accs, names = panel_sets()
    x = load_table()
    motor = load_motor()
    x = x[x.sp5.isin(SPMAP)].copy()
    x["species"] = x.sp5.map(SPMAP)
    x["motor_key"] = [motor_key(r, motor) for r in x.to_dict("records")]
    idm = pd.read_csv(FAM / "cache/idmap.tsv", sep="\t", dtype=str, keep_default_na=False)
    idmap = defaultdict(set)
    for r in idm.itertuples(index=False):
        if r.uniprot_acc:
            idmap[r.original_id].add(r.uniprot_acc)
    fw = pd.read_csv(FAM / "cache/motor_vs_panel.tsv", sep="\t")
    rv = pd.read_csv(FAM / "cache/reverse.tsv", sep="\t")
    rv_best = rv.sort_values(["qseqid", "bitscore"], ascending=[True, False]).groupby("qseqid").head(1)
    rv_best = dict(zip(rv_best.qseqid, rv_best.sseqid))
    fwg = {k: g.sort_values("bitscore", ascending=False) for k, g in fw.groupby("qseqid")}

    rep, mem = [], []
    for r in x.to_dict("records"):
        code, fam, sub = r["species"], clean(r["Family"]), clean(r["Subfamily"])
        label = fam + (f"|{sub}" if sub else "")
        loc = f"S3.XLS sheet Suppl_info row {r['xls_row']}: {r['Analysis ID']}; Family={fam}; Subfamily={sub or '-'}"
        oid = r["Analysis ID"]
        base = dict(original_id=oid, species=code, label=label, locator=loc)
        if clean(r["Excluded (divergent)"]):
            rep.append({**base, "gene_id": "", "status": "excluded",
                        "reason": "divergent sequence excluded from the phylogeny by the authors (Family=ND)"}); continue
        if clean(r["Excluded (near ID included)"]):
            rep.append({**base, "gene_id": "", "status": "excluded",
                        "reason": "near-identical copy (paper: " + clean(r["Excluded (near ID included)"]) + ")"}); continue
        cands, how = set(), []
        acc = clean(r["GenBank Acc.No."])
        if acc:
            a2 = acc.split(".")[0] if (acc[0] in "OPQ" and acc[1].isdigit()) else acc
            hit = idmap.get(a2, set()) & accs[code]
            if hit:
                cands |= hit; how.append(f"UniProt ID mapping of {acc}")
        suffix = oid.split("_", 1)[1].upper() if "_" in oid else ""
        for nm in dict.fromkeys(name_candidates(r["sp5"], r["Genome project ID"]) + [suffix]):
            hit = names[code].get(nm, set())
            if hit:
                cands |= hit; how.append(f"gene name {nm}")
        # motor-domain check
        mk = r["motor_key"]
        best = None
        if mk and mk in fwg:
            g = fwg[mk]
            b = g.iloc[0]
            others = g[g.sseqid != b.sseqid]
            margin = b.pident - (others.iloc[0].pident if len(others) else 0.0)
            best = dict(t=b.sseqid, pid=float(b.pident), qcov=b.length / b.qlen, tcov=b.length / b.slen,
                        margin=float(margin), rbh=rv_best.get(b.sseqid) == mk)
        good_motor = best and best["pid"] >= ID_MIN and best["qcov"] >= QCOV_MIN and best["rbh"] and best["margin"] >= MARGIN
        if len(cands) == 1:
            gene = next(iter(cands))
            if best and best["pid"] >= ID_MIN and best["t"] != gene and best["qcov"] >= QCOV_MIN:
                rep.append({**base, "gene_id": "", "status": "ambiguous",
                            "reason": f"ID mapping -> {gene} but motor domain best hit {best['t']} ({best['pid']:.1f}%)"}); continue
            chk = (f"motor-domain check {best['pid']:.1f}% id to {best['t']}" if best else "no motor-domain sequence in DatasetS1")
            rep.append({**base, "gene_id": gene, "status": "mapped", "reason": "; ".join(how) + "; " + chk,
                        "_m": dict(mapping_method="id_mapping", mapping_identity=f"{best['pid']:.1f}" if best and best["t"] == gene else "",
                                   mapping_coverage="")})
        elif len(cands) > 1:
            if good_motor and best["t"] in cands:
                gene = best["t"]
                rep.append({**base, "gene_id": gene, "status": "mapped",
                            "reason": "; ".join(how) + f" -> {len(cands)} entries; resolved by motor domain {best['pid']:.1f}% RBH",
                            "_m": dict(mapping_method="id_mapping", mapping_identity=f"{best['pid']:.1f}", mapping_coverage="")})
            else:
                rep.append({**base, "gene_id": "", "status": "ambiguous",
                            "reason": "ID maps to several panel entries: " + ",".join(sorted(cands))})
        else:
            if good_motor:
                rep.append({**base, "gene_id": "", "status": "unmapped",
                            "reason": (f"no ID route; motor domain only: {best['pid']:.1f}% id to {best['t']} over "
                                       f"{best['qcov']:.0%} of query but {best['tcov']:.0%} of target (<80% coverage, protocol 0.5)")})
            elif best:
                rep.append({**base, "gene_id": "", "status": "unmapped",
                            "reason": f"no ID route; best motor-domain hit {best['t']} {best['pid']:.1f}% (qcov {best['qcov']:.0%}, margin {best['margin']:.1f}, rbh {best['rbh']})"})
            else:
                rep.append({**base, "gene_id": "", "status": "unmapped", "reason": "no ID route and no motor-domain sequence/hit"})
    # same UniProt entry from >1 paper row -> keep none if labels conflict, else keep first, others excluded
    byg = defaultdict(list)
    for i, r in enumerate(rep):
        if r["status"] == "mapped":
            byg[r["gene_id"]].append(i)
    for g, ii in byg.items():
        if len(ii) > 1:
            labs = {rep[i]["label"] for i in ii}
            if len(labs) > 1:
                for i in ii:
                    rep[i].update(status="ambiguous", reason=rep[i]["reason"] + f"; {g} also mapped by another row with a different label")
                    rep[i].pop("_m", None)
            else:
                for i in ii[1:]:
                    rep[i].update(status="excluded", reason=f"same UniProt entry {g} as {rep[ii[0]]['original_id']} (isoform/duplicate model)")
                    rep[i].pop("_m", None)
    for r in rep:
        if r["status"] == "mapped":
            mem.append(dict(gene_id=r["gene_id"], species=r["species"], original_id_in_paper=r["original_id"], label=r["label"],
                            **r["_m"], source=SRC, locator=r["locator"]))
    members = pd.DataFrame(mem).sort_values(["label", "species", "gene_id"]).reset_index(drop=True)
    members.to_csv(FAM / "members.tsv", sep="\t", index=False, quoting=csv.QUOTE_MINIMAL)
    pd.DataFrame([{k: r.get(k, "") for k in ("original_id", "species", "label", "gene_id", "status", "reason", "locator")}
                  for r in rep]).to_csv(FAM / "mapping_report.tsv", sep="\t", index=False)

    # ---------------------------------------------------------------- tree
    ids = defaultdict(int)

    def nid(base):
        base = re.sub(r"[^A-Za-z0-9_]", "_", base)
        ids[base] += 1
        return base if ids[base] == 1 else f"{base}_{ids[base]}"

    def terminal(name, genes):  # genes: list of (gene, species)
        if len(genes) == 1:
            return Node(genes[0][0])
        spc = defaultdict(int)
        for _, s in genes:
            spc[s] += 1
        multi = any(v > 1 for v in spc.values())
        i = nid(name)
        nd = Node(i, features={"id": i, "ev": "U" if multi else "S", "tier": "T2", "src": SRC})
        for g, _ in sorted(genes):
            nd.add_child(Node(g))
        return nd

    def group(name, children, ev="U", **feat):
        kids = [c for c in children if c is not None]
        if not kids:
            return None
        if len(kids) == 1:
            return kids[0]
        i = nid(name)
        nd = Node(i, features={"id": i, "ev": ev, "tier": feat.pop("tier", "T2"), "src": SRC, **feat})
        for c in kids:
            nd.add_child(c)
        return nd

    fams = defaultdict(lambda: defaultdict(list))  # fam -> sub ('' = family-level only) -> genes
    loose = []
    for m in members.itertuples(index=False):
        fam, _, sub = m.label.partition("|")
        if fam.startswith("(") or fam in ("ungrouped", ""):
            loose.append(m.gene_id); continue
        if sub.startswith("(") or not sub:
            sub = ""
        fams[fam][sub].append((m.gene_id, m.species))

    def family_node(fam):
        subs = fams[fam]
        named = sorted(s for s in subs if s)
        if fam == "Kinesin-9":
            ab = group("K9A_9B_LECA", [terminal("K9A", subs.get("9A", [])) if subs.get("9A") else None,
                                        terminal("K9B", subs.get("9B", [])) if subs.get("9B") else None],
                       ev="D", tier="T1", lvl_lo="Eukaryota", lvl_hi="ABOVE_RANGE")
            rest = [Node(g) for g, _ in subs.get("", [])]
            return group("K9_family", [ab] + rest)
        if not named:
            return terminal(fam.replace("Kinesin-", "K"), subs[""])
        kids = [terminal(f"K{s}", subs[s]) for s in named] + [Node(g) for g, _ in sorted(subs.get("", []))]
        return group(fam.replace("Kinesin-", "K") + "_family", kids)

    mar_nodes = [family_node(f) for f in MAR if f in fams]
    leca = group("KIN_MAR_LECA", mar_nodes, ev="D", tier="T1", lvl_lo="Eukaryota", lvl_hi="ABOVE_RANGE")
    others = [family_node(f) for f in sorted(fams, key=lambda s: (s.startswith("X"), int(re.sub(r"\D", "", s.split("/")[0]) or 0))) if f not in MAR]
    root = group("KIN_root", [leca] + others + [Node(g) for g in sorted(loose)])
    (FAM / "reference_tree.nhx").write_text(write_newick(root, features=["ev", "lvl", "lvl_lo", "lvl_hi", "id", "tier", "src", "tag"], lengths=False) + "\n")

    # ---------------------------------------------------------------- assertions
    present = {f: sorted({s for g in fams[f].values() for _, s in g}) for f in fams}
    rows = [dict(assertion_id="E01_A1", type="event", target="KIN_MAR_LECA", tier="T1", source_id=SRC,
                 locator=f"{TXT_REL} Results, 'The kinesin repertoire of the last common eukaryotic ancestor', para 4",
                 quote=Q_MAR),
            dict(assertion_id="E01_A2", type="event", target="KIN_MAR_LECA", tier="T1", source_id=SRC,
                 locator=f"{TXT_REL} Results, same paragraph", quote=Q_MAR2),
            dict(assertion_id="E01_A3", type="event", target="K9A_9B_LECA", tier="T1", source_id=SRC,
                 locator=f"{TXT_REL} Results, same paragraph (9A and 9B listed separately in the MAR)", quote=Q_MAR),
            dict(assertion_id="E01_A4", type="tag", target="subfamily labels (T2 S/U groups)", tier="T2", source_id=SRC,
                 locator=f"{TXT_REL} Results, 'Kinesin families and subfamilies'", quote=Q_SUBFAM),
            dict(assertion_id="E01_A5", type="tag", target="X1-X14 groups (outside LECA node; U)", tier="T2", source_id=SRC,
                 locator=f"{TXT_REL} Results", quote=Q_X)]
    rp0 = pd.read_csv(FAM / "mapping_report.tsv", sep="\t", dtype=str, keep_default_na=False)
    inp0 = rp0[rp0.status.isin(["mapped", "unmapped", "ambiguous"])]
    fr0 = inp0.groupby("species").status.apply(lambda s: (s == "mapped").mean())
    for k, c in enumerate(sorted(fr0[fr0 >= 0.9].index), 1):
        rows.append(dict(assertion_id=f"E01_C{k}", type="completeness", target=c, tier="T1", source_id=SRC,
                         locator=f"{TXT_REL} Methods (kinesin identification); Additional file 7 lists the proteome used",
                         quote=Q_COMPLETE))
    with open(FAM / "assertions.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader(); w.writerows(rows)

    # ---------------------------------------------------------------- family.yaml
    rp = pd.read_csv(FAM / "mapping_report.tsv", sep="\t", dtype=str, keep_default_na=False)
    inp = rp[rp.status.isin(["mapped", "unmapped", "ambiguous"])]
    frac_sp = inp.groupby("species").status.apply(lambda s: (s == "mapped").mean())
    complete = sorted(frac_sp[frac_sp >= 0.9].index)
    panel_euk = ["HUMAN", "DROME", "NEMVE", "AMPQE", "MNELE", "TRIAD", "MONBE", "CAPO3", "YEAST", "NEUCR", "BATDJ",
                 "ENCCU", "DICDI", "ENTH1", "ARATH", "CHLRE", "CYAME", "THAPS", "PHYRM", "PLAF7", "TRYB2", "NAEGR",
                 "GIAIC", "THETR", "ICHMU"]
    absent = [c for c in panel_euk if c not in set(SPMAP.values())]
    meta = dict(
        id="E01_kinesins", name="Kinesin superfamily (eukaryote-wide; LCEA minimal ancestral repertoire)", panel=PANEL,
        collection="E_cytoskeleton_motors", phenomena=["DEEP", "LOSS", "PROM", "LSE"],
        descriptors=dict(n_members=int(len(members)), n_species=int(members.species.nunique()),
                         mapped_fraction=round(float((inp.status == "mapped").mean()), 3),
                         n_rows_panel_species=int(len(rp)),
                         status_counts={k: int(v) for k, v in rp.status.value_counts().items()},
                         per_species_mapped_fraction={k: round(float(v), 3) for k, v in frac_sp.items()},
                         families_present={f: present[f] for f in sorted(present)}),
        sources=[dict(id=SRC, doi="10.1186/1471-2148-10-110", pmid="20423470",
                      title="Patterns of kinesin evolution reveal a complex ancestral eukaryote with a multifunctional cytoskeleton",
                      year=2010, source_method="phylogeny", table="Additional file 3 (S3.XLS)", locator="Results; Figure 1, 3; Additional files 3, 8",
                      quote=Q_MAR,
                      files=[dict(path=p, sha256=sha(ROOT / p)) for p in (XLS_REL, ALN_REL, TXT_REL)])],
        root_level="Eukaryota", root_single_copy=False,
        complete_species=complete, decoys=[], distinct_from=[],
        tier_summary=dict(T1_D_events=2, T2_S_nodes=None, T2_U_nodes=None),
        protocol_version=1,
        notes=("Kinesin-myosin homology: kinesin and myosin motor domains share a common P-loop NTPase fold/ancestry, so no "
               "distinct_from is declared between E01_kinesins and E06_myosins (homologous at fold level, but no "
               "orthology/paralogy level asserted). D nodes: (1) KIN_MAR_LECA = split among the 10 MAR families "
               "(Kinesin-9 counted once; 9A/9B split is the separate D node K9A_9B_LECA), both lvl=Eukaryota with "
               "interval [Eukaryota, ABOVE_RANGE] (paper infers presence in the LCEA; the duplications may be older). "
               "Non-MAR families (6, 7, 15, 16, 18, 19, 20), tentative X1-X14 groups, parenthesised (uncertain) family "
               "labels and 'ungrouped' sequences hang from the U root (no assertion against MAR families). Within a "
               "family, subfamily splits are undated in the paper -> U. Terminal label groups: <=1 gene per species -> "
               "S polytomy (T2 rule b; no two panel-E species share a WGD, so rule-b pseudo-orthology exception does "
               "not apply); >1 gene in a species -> U (rule c). Mapping: GenBank/RefSeq -> UniProt ID mapping or exact "
               "gene-model names (JGI DRAFT names) in the panel proteome, cross-checked with the paper's motor-domain "
               "sequence (DatasetS1); rows with only a motor-domain sequence (no ID route) fail the 80% target-coverage "
               "rule and are reported unmapped. Source used old (2006-2007) proteome versions; Batde = strain JEL423 "
               "(panel JAM81). Panel-E eukaryotes absent from the source: " + ", ".join(absent) +
               " (source has Phytophthora sojae, not P. ramorum; source code 'Capsp' is the annelid Capitella sp. I, not Capsaspora). complete_species = species with a genome-wide HMMER "
               "search in the source and >=90% of in-panel rows uniquely mapped; caveat: current gene models may "
               "contain kinesins absent from the 2007 datasets."),
    )
    cnt = defaultdict(int)

    def walk(n):
        if n.children:
            cnt[(n.features.get("ev"), n.features.get("tier"))] += 1
            for c in n.children:
                walk(c)
    walk(root)
    meta["tier_summary"] = {f"{ev}_{t}": int(v) for (ev, t), v in sorted(cnt.items())}
    meta["completeness_quote"] = Q_COMPLETE
    (FAM / "family.yaml").write_text(yaml.safe_dump(meta, sort_keys=False, allow_unicode=True, width=120))


if __name__ == "__main__":
    main()
