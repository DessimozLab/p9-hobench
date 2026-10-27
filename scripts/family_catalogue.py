#!/usr/bin/env python
"""Catalogue of frozen reference families -> results/catalogue/{families.tsv,families.md,coverage.md,abandoned.md}.
Usage: python scripts/family_catalogue.py   (run from repo root)"""
from __future__ import annotations
import re, sys, tarfile, collections
from pathlib import Path
import pandas as pd, yaml
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from hobench.reference import load_families, reconcile
from hobench.speciestree import SpeciesTree

FAM, OUT = ROOT / "families", ROOT / "results/catalogue"
OUT.mkdir(parents=True, exist_ok=True)
CONTROL_PREFIX = ("K01_", "K02_", "E05_")
COLLECTIONS = {"Y01_YGOB_pillars", "V07_Xenopus_LS"}
def md(df, index=False):
    if index: df = df.reset_index()
    h = [str(c) for c in df.columns]
    return "\n".join(["| " + " | ".join(h) + " |", "|" + "---|" * len(h)] +
                     ["| " + " | ".join(str(x).replace("|", "/") for x in r) + " |" for r in df.itertuples(index=False)])
PANEL_N = {p: len(pd.read_csv(ROOT / f"panels/{p}/species.tsv", sep="\t")) for p in "EKPVY"}
STS = {p: SpeciesTree.from_newick(ROOT / f"panels/{p}/species_tree.nwk") for p in "EKPVY"}

def track_row_id(fam):
    c = fam.meta.get("collection")
    if fam.id.startswith("K01_"): return "K01_RPROT", "control"
    if fam.id.startswith(("K02_", "E05_")): return fam.id, "control"
    if fam.path.parent.name == "units":
        return fam.path.parent.parent.name, "collection"
    return fam.id, "family"

# ---- NCBI ranks for depth bins ------------------------------------------------------
def load_ranks():
    need = {}
    for p in "EKPVY":
        lv = pd.read_csv(ROOT / f"panels/{p}/levels.tsv", sep="\t", dtype=str)
        for t in lv.ncbi_taxid.dropna(): need[int(t)] = 1
    parent, rank = {}, {}
    with tarfile.open(ROOT / "taxdump.tar.gz") as tf:
        for line in tf.extractfile("nodes.dmp"):
            f = line.decode().split("\t|\t")
            parent[int(f[0])] = int(f[1]); rank[int(f[0])] = f[2]
    return parent, rank
PARENT, RANK = load_ranks()
RECENT = {"species", "subspecies", "strain", "genus", "subgenus", "species group", "species subgroup", "family",
          "subfamily", "tribe", "subtribe", "section", "series", "varietas", "forma", "isolate", "serotype", "biotype"}
INTER = {"order", "suborder", "infraorder", "parvorder", "superfamily", "class", "subclass", "infraclass", "superclass",
         "cohort", "subcohort", "superorder", "division", "subdivision", "subphylum", "infraphylum", "microphylum"}
def rank_bin(r):
    return "recent" if r in RECENT else "intermediate" if r in INTER else "deep"
def level_bins(panel):
    lv = pd.read_csv(ROOT / f"panels/{panel}/levels.tsv", sep="\t", dtype=str, keep_default_na=False)
    out = {}
    for _, r in lv.iterrows():
        if r.level_id in ("LUCA", "Eukaryota"): out[r.level_id] = "LECA/LUCA"; continue
        t = int(r.ncbi_taxid) if r.ncbi_taxid else 1
        rk = RANK.get(t, "no rank")
        while rk in ("no rank", "clade") and t in PARENT and PARENT[t] != t:
            t = PARENT[t]; rk = RANK.get(t, "no rank")
        out[r.level_id] = rank_bin(rk)
    return out
BINS = {p: level_bins(p) for p in "EKPVY"}
def bin_of(panel, name, st):
    if st.is_species(name): return "recent"
    return BINS[panel].get(name, "deep")

# ---- text sources ------------------------------------------------------------------
def verification_map():
    rows = []
    for l in (FAM / "VERIFICATION.md").read_text().splitlines():
        if l.startswith("|") and not l.startswith("|---"):
            c = [x.strip() for x in l.strip("|").split("|")]
            if len(c) >= 5 and c[0] != "family": rows.append((c[0], c[4]))
    return rows
VER = verification_map()
def verdict(fid):
    for k, v in VER:
        key = k.split()[0]
        if fid == key or fid.startswith(key + "_") or fid.startswith(key + " "): return re.sub(r"\*", "", v)
    if fid.startswith("K01_"): return next(v for k, v in VER if k.startswith("K01"))
    if re.match(r"K04[a-j]_", fid): return next(v for k, v in VER if k.startswith("K04a"))
    if fid.startswith("YGOB_"): return next(v for k, v in VER if k.startswith("Y01"))
    if fid.startswith("XL_"): return next(v for k, v in VER if k.startswith("V07"))
    return "not requested/unlisted"
ERR = (FAM / "ERRATA.md").read_text()

def log_status(d):
    p = d / "curation_log.md"
    if not p.exists(): return ""
    for l in p.read_text().splitlines()[:6]:
        if re.match(r"(?i)\s*status", l): return l
    return ""

def mapping_stats(dirs, collection=False):
    m = u = 0
    for d in dirs:
        p = d / "mapping_report.tsv"
        if not p.exists(): continue
        df = pd.read_csv(p, sep="\t", dtype=str, keep_default_na=False)
        if "unit_status" in df: df = df[df.unit_status == "unit"]
        s = df.status.value_counts()
        m += s.get("mapped", 0); u += s.get("mapped", 0) + s.get("unmapped", 0) + s.get("ambiguous", 0)
    return (m / u) if u else float("nan")

def author_year(meta):
    s = (meta.get("sources") or [{}])[0]
    return f"{s.get('id', '?')} {s.get('year', '')} doi:{s.get('doi', s.get('pmid', '?'))}".replace("  ", " ")

TAGMAP = {"WGD": "WGD", "ALLO": "ALLO", "LOSS": "LOSS", "TAND": "TAND/LSE", "LSE": "TAND/LSE", "HGT": "HGT", "EGT": "HGT",
          "DEEP": "DEEP/LECA/LUCA", "LECA": "DEEP/LECA/LUCA", "LUCA": "DEEP/LECA/LUCA", "PROM": "PROM", "FUS": "FUS",
          "FAST": "FAST", "ANALOG": "ANALOG", "CONTROL": "CONTROL"}

def main():
    fams = load_families(FAM)
    groups = collections.OrderedDict()
    for f in fams:
        rid, track = track_row_id(f)
        groups.setdefault(rid, {"track": track, "fams": []})["fams"].append(f)
    rows, ev_rows = [], []
    for rid, g in groups.items():
        fl = g["fams"]; panel = fl[0].panel; st = STS[panel]
        genes, spp = set(), set(); tiers = collections.Counter(); nint = nleaf = nS = nH = 0
        levels = []; ndec = ncs = 0; tags = set(); ntree_events = 0
        for f in fl:
            genes |= set(f.members.gene_id); spp |= set(f.members.species)
            ndec += len(f.decoys) + len(f.distinct_from); ncs += len(f.complete_species)
            tags |= {TAGMAP.get(t, t) for t in f.phenomena}
            rf = reconcile(f, st, check=False)
            for e in rf.events:
                if e.ev == "S": nS += 1
                elif e.ev in ("D", "H"):
                    if e.ev == "H": nH += 1
                    else:
                        tiers[e.tier or "NA"] += 1
                        if st.is_leaf[e.lvl_lo]: nleaf += 1
                        else: nint += 1
                    nm = st.names[e.lvl_lo]
                    levels.append(nm)
                    ev_rows.append((rid, panel, e.ev, e.tier or "NA", bin_of(panel, nm, st), nm))
        if g["track"] == "collection": d0 = FAM / rid; meta = yaml.safe_load((d0 / "collection.yaml").read_text()); dirs = [d0]
        else: dirs = [f.path for f in fl]; meta = fl[0].meta
        pct = mapping_stats(dirs)
        if pct != pct:  # fall back to descriptors
            vals = [f.meta.get("descriptors", {}).get("mapped_fraction") for f in fl]
            vals = [v for v in vals if v is not None]; pct = sum(vals) / len(vals) if vals else float("nan")
        # headline
        why = []
        ls = " ".join(log_status(d) for d in dirs[:1])
        ver = verdict(rid if g["track"] != "control" or rid != "K01_RPROT" else "K01_RPROT_x")
        if g["track"] == "control" and rid == "K01_RPROT":
            ok = [mapping_stats([f.path]) >= 0.8 for f in fl]
            headline = f"y ({sum(ok)}/{len(ok)} units >=80% mapped; rest non-headline)"
        else:
            if pct == pct and pct < 0.8: why.append(f"mapped {pct:.0%} < 80%")
            if "FAIL" in ver and "PASS" not in ver.split("FAIL")[0][-4:]: why.append("verification FAIL")
            if not tiers and not nS and not nH: why.append("no events")
            headline = "y" if not why else "n: " + "; ".join(why)
            if g["track"] == "collection": headline += " (per-unit)"
        dm = (levels and sorted(levels, key=lambda n: int(st.depth[st.idx(n)]))) or []
        rows.append(dict(id=rid, panel=panel, track=g["track"], headline=headline, n_units=len(fl), n_members=len(genes),
                         n_species=f"{len(spp)}/{PANEL_N[panel]}", pct_unique_mapping=("%.1f" % (100 * pct)) if pct == pct else "NA",
                         phenomena=",".join(sorted(tags)),
                         D_T1=tiers.get("T1", 0), D_T2=tiers.get("T2", 0), D_T3=tiers.get("T3", 0) + tiers.get("T3fig", 0),
                         D_other=sum(v for k, v in tiers.items() if k not in ("T1", "T2", "T3", "T3fig")),
                         D_internal=nint, D_leaf=nleaf, n_S=nS, n_H=nH,
                         deepest_event=(dm[0] if dm else ""), shallowest_event=(dm[-1] if dm else ""),
                         n_decoys_distinct=ndec, n_complete_species=ncs, verification=ver,
                         source=author_year(meta) + (" (+%d more)" % (len(meta.get("sources") or []) - 1) if len(meta.get("sources") or []) > 1 else "")))
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "families.tsv", sep="\t", index=False)
    (OUT / "families.md").write_text("# HOBench family catalogue\n\nGenerated by scripts/family_catalogue.py. T3 column includes T3fig. "
        "D leaf = duplication whose (youngest-bound) level is a species leaf; D internal = any ancestral level.\n\n" + md(df) + "\n")
    ev = pd.DataFrame(ev_rows, columns=["row", "panel", "ev", "tier", "bin", "level"])
    ev.to_csv(OUT / "events_by_level.tsv", sep="\t", index=False)
    coverage(df, ev)
    abandoned()
    print(df.groupby("track")[["n_units", "n_members", "D_T1", "D_T2", "D_T3", "n_H"]].sum())
    return df, ev

def coverage(df, ev):
    out = ["# Coverage\n"]
    cols = ["WGD", "ALLO", "LOSS", "TAND/LSE", "HGT", "DEEP/LECA/LUCA", "PROM", "FUS", "FAST", "ANALOG", "CONTROL"]
    out.append("## Panel x phenomenon (number of catalogue rows tagged; K01 = 1 row, Y01/V07 = 1 row each; tag on any unit counts)\n")
    mat = pd.DataFrame({c: [sum(1 for _, r in df[df.panel == p].iterrows() if c in r.phenomena.split(",")) for p in "EKPVY"] for c in cols}, index=list("EKPVY"))
    mat["ANALOG (rows with decoys/distinct_from)"] = [int((df[df.panel == p].n_decoys_distinct > 0).sum()) for p in "EKPVY"]
    mat["rows"] = [int((df.panel == p).sum()) for p in "EKPVY"]
    out.append(md(mat, True) + "\n")
    allt = sorted({t for s in df.phenomena for t in s.split(",") if t})
    out.append("Tag vocabulary seen (normalised): " + ", ".join(allt) + "\n")
    out.append("## Depth coverage (D and H events)\n")
    out.append("Binning rule (no dated species trees are stored, so bins are an NCBI-rank proxy for age): each event is placed at its youngest-bound level "
               "`lvl_lo` (the species-tree node of the LCA mapping/lower interval bound). Level -> NCBI rank via `ncbi_taxid` in panels/<P>/levels.tsv "
               "(no-rank/clade nodes inherit the nearest ranked ancestor). **recent** (~<100 Mya): species leaf, or rank species/genus/family (incl. sub-/tribe); "
               "**intermediate**: order, class, cohort, subphylum and their super-/infra- ranks; **deep**: phylum, kingdom, domain and above; "
               "**LECA/LUCA**: level `LUCA` or `Eukaryota` (panel nodes exactly at LUCA / LECA). Intervals are binned by lower bound, so deep-hi intervals are counted young. "
               "Ranks are only a proxy: e.g. Chordata-level (phylum) events are ~550 Mya, family-level may exceed 100 Mya in bacteria.\n")
    b = ev.groupby(["panel", "bin"]).size().unstack(fill_value=0).reindex(columns=["recent", "intermediate", "deep", "LECA/LUCA"], fill_value=0)
    b.loc["all"] = b.sum(); out.append("All D+H events (units counted individually):\n\n" + md(b, True) + "\n")
    hl = set(df[df.headline.str.startswith("y")].id)
    e2 = ev[ev.row.isin(hl) & ev.tier.isin(["T1", "T2"])]
    b = e2.groupby(["panel", "bin"]).size().unstack(fill_value=0).reindex(columns=["recent", "intermediate", "deep", "LECA/LUCA"], fill_value=0)
    b.loc["all"] = b.sum(); out.append("Headline rows only, T1+T2 D events (H events have no tier and are excluded):\n\n" + md(b, True) + "\n")
    b = ev.groupby(["bin", "tier"]).size().unstack(fill_value=0); out.append("By tier:\n\n" + md(b, True) + "\n")
    out.append("""## Comparison with existing benchmarks
Facts for the other benchmarks are taken only from LITERATURE.md and literature/A_benchmarks_metrics.md; "not stated" = not in those files.

| | HOBench (this work) | OrthoBench (Trachana 2011; Emms & Kelly 2020) | SwissTree (Boeckmann 2011) | QfO benchmark service (Altenhoff 2016/2020; Nevers 2022; Altenhoff 2024) |
|---|---|---|---|---|
| Taxonomic scope | 5 panels (E eukaryotes+archaea+bacteria 34 spp; K prokaryote-centred 21; P plants 14; V vertebrates 17; Y yeasts 16) | 12 metazoans; Bilateria only, animal-centred | 19 expert-curated reference gene trees spanning all domains (also described as few "tractable" families) | Reference proteome set of ~78-81 species (48 eukaryotes, 23 bacteria, 7 archaea in 2019-22); test sets vary (VGNC: 8 mammals) |
| Levels scored | Every species-tree node (panel levels.tsv; see depth table) | One (Bilateria); described as a limitation | Pairwise orthology only; no level concept stated | Pairwise; no test scores nesting; OrthoBench group test prototyped, one level |
| Events asserted | Per-family reconciled trees: S, D (exact level or interval), H, with tiers T1/T2/T3 | 70 RefOGs (groups only; no events) | Reconciled/curated reference trees (events implicit in trees; pairwise scored) | None scored: duplication timing, paralogy and nesting are not scored |
| HGT / WGD content | Explicit H nodes (Y03, Y04, K03 etc.) and WGD/ALLO D tags (Y01, V03, V07, ...) | Not stated; challenge tags are duplication/loss, rate, domain architecture, alignment quality | Not stated | YGOB-style WGD "traps" appear only in the separate Salichos 2011 yeast benchmark; no HGT test stated |
| Hierarchical scoring | Per-level HOG P/R/F, B-cubed, three-class relation, stratified by depth (design) | Single-level pair P/R/F per RefOG, (n-1) normalised | Pairwise precision/recall | HOGs flattened to ortholog pairs, nesting discarded |
""")
    (OUT / "coverage.md").write_text("\n".join(out))

def abandoned():
    out = ["# Abandoned / non-built families and units\n", "Source: families/*/curation_log.md (directories lacking family.yaml/collection.yaml), plus the K01 unit table in the shared RPROT log.\n",
           "| id | reason (from curation log) |", "|---|---|"]
    k01 = None
    for d in sorted(FAM.iterdir()):
        if not d.is_dir() or d.name.startswith("_") or (d / "family.yaml").exists() or (d / "collection.yaml").exists(): continue
        if d.name.startswith("K01_RPROT_"):
            k01 = k01 or next(FAM.glob("K01_RPROT_*/curation_log.md"))
            unit = d.name.split("_")[-1]
            for l in Path(k01).read_text().splitlines():
                if l.startswith(f"| {unit} |"):
                    c = [x.strip() for x in l.strip("|").split("|")]; out.append(f"| {d.name} | {c[1]}: {c[2]} |"); break
            else: out.append(f"| {d.name} | no unit-table row; see shared K01 log |")
            continue
        txt = (d / "curation_log.md").read_text() if (d / "curation_log.md").exists() else ""
        why = re.search(r"(?s)## Why abandoned\s*\n(.*?)(\n## |\Z)", txt)
        reason = (why.group(1).strip() if why else log_status(d) or "no curation log")
        out.append(f"| {d.name} | {' '.join(reason.split())[:600]} |")
    # built but below headline / excluded within built families
    out.append("\n## Built but not headline-eligible\n")
    df = pd.read_csv(OUT / "families.tsv", sep="\t", dtype=str, keep_default_na=False)
    for _, r in df[df.headline.str.startswith("n")].iterrows(): out.append(f"- {r.id}: {r.headline[3:]}")
    out.append("- K01_RPROT units L21, L29, L32, S16: <80 % uniquely mapped (built, non-headline; per shared log)")
    for c in sorted(COLLECTIONS):
        mr = pd.read_csv(FAM / c / "mapping_report.tsv", sep="\t", dtype=str, keep_default_na=False)
        u = mr.drop_duplicates("unit")
        out.append(f"\n## {c}: candidate units not retained\n")
        for k, v in u[u.unit_status != "unit"].unit_status.str.replace(r"^not_a_unit: ", "", regex=True).value_counts().items(): out.append(f"- {v} candidates: {k}")
    (OUT / "abandoned.md").write_text("\n".join(out) + "\n")

if __name__ == "__main__":
    main()
