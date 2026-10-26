#!/usr/bin/env python
"""E07 adaptor-complex subunit families: build members.tsv, mapping_report.tsv, reference_tree.nhx, assertions.tsv,
family.yaml for the family whose directory holds this file (the same build.py is used by E07a_adaptins_large,
E07b_adaptins_medium_small and E07c_adaptins_scaffold; FAMILY is chosen from the directory name).

Source: Hirst J, Schlacht A, Norcott JP, Traynor D, Bloomfield G, Antrobus R, Kay RR, Dacks JB, Robinson MS (2014)
Characterization of TSET, an ancient and widespread membrane trafficking complex. eLife 3:e02866 (PMID 24867644).
Figure 3-source data 1 (every TSET / COPI / AP subunit sequence used in the phylogenies: section = subunit,
annotation, accession). Figure 1-source data 1/2 (reverse-HHpred hits; place beta-/gamma-COP among the large and
delta-/zeta-COP among the medium/small subunit homologues). Network results are cached in cache/ by
families/E07a_adaptins_large/prep_mapping.py; this script is offline.
Run from the project root:  python families/<E07 family>/build.py
"""
import csv, hashlib, re, sys
from collections import defaultdict
from pathlib import Path
import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[2]
FAM = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(FAM))
from prep_mapping import load_table, id_tokens, XLSX_REL  # noqa: E402
from hobench.newick import Node, write_newick  # noqa: E402
from hobench.speciestree import SpeciesTree  # noqa: E402

FID = FAM.name
SRC, PANEL = "Hirst2014", "E"
PAP = "literature/papers/PMID24867644"
TXT_REL = f"{PAP}/fulltext.txt"
F1D1 = f"{PAP}/raw/supp/elife-02866-fig1-data1.xlsx"
F1D2 = f"{PAP}/raw/supp/elife-02866-fig1-data2.xlsx"
Q = dict(
    APLECA=("our most recent analyses indicate that, contrary to previous reports (Hirst et al., 2011), AP-5 diverged "
            "basally within the adaptin clade, followed by AP-3, AP-4, and APs 1 and 2, all prior to the LECA."),
    APCOPI=("The AP and COPI complexes share a similar architecture, due to common ancestry predating the last "
            "eukaryotic common ancestor (LECA)."),
    TSETLECA="Although TSET is deduced to have been present in LECA",
    TSETDIST=("Thus, TSET is clearly an ancient component of the eukaryotic membrane-trafficking system, distinct from "
              "the known heterotetramers."),
    LARGE2=("The gene encoding the large subunit then duplicated and mutated to generate the two distinct types of "
            "large subunits (red and magenta)"),
    SMALL2=("the gene encoding the small subunit also duplicated and mutated (yellow and orange), with one of the two "
            "proteins (orange) acquiring a μ homology domain (MHD) to form the ancestral heterotetramer"),
    SCAF=("Phylogenetic analysis of the TTRAYs and their closest relatives, β′-COP and α-COP (Figure 4—figure "
          "supplement 9), showed that the paralogues are due to ancient duplications in the TSET and COPI families "
          "respectively, which occurred prior to the divergence of the LECA."),
    HOMOL="The large adaptor subunits share sequence and structure homology, as do the medium and small subunits.",
    TPLATE="an Arabidopsis protein related to the AP β subunits and β-COP",
    TSAUCER="(the γαδεζ-like subunit)",
    TCUP="(the μ-like subunit)",
    TSPOON="the higher ratios for the small (σ-like/TSPOON) subunit",
    F1D1="Figure 1—source data 1. Large subunit homologues found by reverse HHpred in different organisms.",
    F1D2="Figure 1—source data 2. Medium and small subunit homologues found by reverse HHpred in different organisms.",
    BCOP="coatomer subunit beta-1 [Arabidopsis thaliana]",
    GCOP="coatomer subunit gamma [Arabidopsis thaliana]",
    LOSS=("Although TSET is deduced to have been present in LECA, the complex appears to have been entirely or "
          "partially lost in various lineages (Figure 3B)."),
    COMPLETE="HMMs were used to search the protein databases for the organisms in Figure 3A",
)
FAMS = {
    "E07a_adaptins_large": dict(
        name="Adaptor-complex large subunits (AP-1..5 beta and gamma/alpha/delta/epsilon/zeta, beta-/gamma-COP, TPLATE, TSAUCER)",
        sections=["TPLATE", "TSAUCER", "COPI BETA", "COPI GAMMA", "AP Beta", "AP GADEZ"],
        leca_q=["APLECA", "APCOPI", "TSETLECA", "LARGE2"],
        member_q=[("tag", "TPLATE = beta-like large subunit", "TPLATE"), ("tag", "TSAUCER = gamma-like large subunit", "TSAUCER"),
                  ("tag", "beta-COP / gamma-COP among large-subunit homologues", "F1D1"),
                  ("tag", "beta-COP in Fig1-source data 1", "BCOP"), ("tag", "gamma-COP in Fig1-source data 1", "GCOP"),
                  ("tag", "large subunits homologous", "HOMOL")]),
    "E07b_adaptins_medium_small": dict(
        name="Adaptor-complex medium and small subunits (AP-1..5 mu and sigma, delta-/zeta-COP, TCUP, TSPOON)",
        sections=["TCUP", "TSPOON", "COPI DELTA", "COPI ZETA", "AP Mu", "AP Sigma"],
        leca_q=["APLECA", "APCOPI", "TSETLECA", "SMALL2"],
        member_q=[("tag", "TCUP = mu-like subunit", "TCUP"), ("tag", "TSPOON = sigma-like subunit", "TSPOON"),
                  ("tag", "delta-COP / zeta-COP among medium/small homologues", "F1D2"),
                  ("tag", "medium and small subunits homologous", "HOMOL")]),
    "E07c_adaptins_scaffold": dict(
        name="COPI / TSET scaffolding subunits (alpha-COP, beta'-COP, TTRAY1, TTRAY2)",
        sections=["TTRAY1", "TTRAY2", "COPI ALPHA", "COPI BETA'"],
        leca_q=["SCAF"],
        member_q=[]),
}
CFG = FAMS[FID]
AP_GROUP = {"AP Beta": {"B": None}, "AP Mu": {"M": None}, "AP Sigma": {"S": None},
            "AP GADEZ": {"G": "1", "A": "2", "D": "3", "E": "4", "Z": "5"}}
APNAME = {"B": "beta", "M": "mu", "S": "sigma", "G": "gamma", "A": "alpha", "D": "delta", "E": "epsilon", "Z": "zeta"}


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def group_of(section, ann):
    """subunit group label, or '' if the annotation does not parse consistently with its section."""
    if not section.startswith("AP "):
        return section.replace("COPI ", "COPI-").replace("BETA'", "BETAprime")
    m = re.search(r"([1-5])([BMSGADEZ])\d*$", ann)
    if not m:
        return ""
    cx, let = m.groups()
    allowed = AP_GROUP[section]
    if let not in allowed or (allowed[let] and allowed[let] != cx):
        return ""
    if let == "B" and cx in ("1", "2"):
        return "AP1/2-beta"
    return f"AP{cx}-{APNAME[let]}"


DRAFT = {"BATDJ": ["BATDEDRAFT_{}"], "MONBE": ["MONBRDRAFT_{}"], "NAEGR": ["NAEGRDRAFT_{}"],
         "THAPS": ["THAPSDRAFT_{}", "THAPS_{}"], "NEMVE": ["NEMVEDRAFT_v1g{}"], "TRIAD": ["TRIADDRAFT_{}"]}
JGI_ORG = {"Batde5": "BATDJ", "Monbr1": "MONBE", "Naegr1": "NAEGR", "Thaps3": "THAPS", "Nemve1": "NEMVE"}


def name_candidates(tok, code, acc_cell):
    """gene-model names to look up in the species' panel header map (exact, single hit required)."""
    t = tok.strip()
    c = []
    for pat, f in ((r"^(AT\dG\d+)(\.\d+)?$", r"\1"), (r"^(EHI_\d+)$", r"\1"), (r"^(PF3D7_\d+)$", r"\1"),
                   (r"^(Tb927\.\d+\.\d+)$", r"\1"), (r"^(AMSG_\d+)T\d+$", r"\1"), (r"^(NCU\d+)T\d+$", r"\1"),
                   (r"^(DDB_G\d+)$", r"\1"), (r"^(CM[A-Z]\d+C)$", r"CYME_\1")):
        if re.match(pat, t):
            c.append(re.sub(pat, f, t))
    m = re.match(r"^Cre(\d+)\.g(\d+)\.t\d", t)
    if m:
        c.append(f"CHLRE_{m.group(1)}g{m.group(2)}v5")
    m = re.match(r"^CAOG_(\d+)T\d+$", t)
    if m:
        c.append(f"CAOG_{m.group(1).zfill(6)}")
    m = re.match(r"^GL50803_(\d+)$", t)
    if m:
        c += [f"GL50803_{m.group(1).zfill(n)}" for n in (6, 7, 8)]
    # JGI protein ids: 'jgi|Org|<id>|...' or a bare number in a JGI-sourced species -> <ORG>DRAFT_<id>
    m = re.match(r"^jgi\|(\w+)\|(\d+)\|", acc_cell)
    num = m.group(2) if (m and JGI_ORG.get(m.group(1)) == code and t == m.group(2)) else (t if re.match(r"^\d+$", acc_cell.strip()) else "")
    if num and code in DRAFT:
        c += [f.format(num) for f in DRAFT[code]]
    return [x.upper() for x in c]


def main():
    st = SpeciesTree.from_newick(ROOT / f"panels/{PANEL}/species_tree.nwk")
    t = load_table()
    t = t[t.section.isin(CFG["sections"])]
    p = t[t.species != ""].copy()
    sp = pd.read_csv(ROOT / f"panels/{PANEL}/species.tsv", sep="\t", dtype=str).set_index("code")
    accs, names, lens = {}, {}, {}
    for code in p.species.unique():
        fa = ROOT / f"panels/{PANEL}/proteomes/{code}.fa"
        accs[code] = {l[1:].split()[0] for l in open(fa) if l.startswith(">")}
        d = defaultdict(set)
        hm = ROOT / Path(sp.loc[code, "fasta_path"]).parent / "header_map.tsv"
        h = pd.read_csv(hm, sep="\t", dtype=str, keep_default_na=False, usecols=["accession", "gene_name", "gene_symbol_gene2acc", "length"])
        lens[code] = {a: int(L) for a, L in zip(h.accession, h.length) if L.isdigit()}
        for a, g1, g2 in h[["accession", "gene_name", "gene_symbol_gene2acc"]].itertuples(index=False):
            if a in accs[code]:
                for g in (g1, g2):
                    for tok in re.split(r"[;,/ ]+", g):
                        tok = tok.split(":")[-1]
                        if tok:
                            d[tok.upper()].add(a)
        names[code] = d
    idm = pd.read_csv(FAM / "cache/idmap.tsv", sep="\t", dtype=str, keep_default_na=False)
    idmap = defaultdict(set)
    for r in idm.itertuples(index=False):
        if r.uniprot_acc:
            idmap[r.original_id].add(r.uniprot_acc)
    seqm = {}
    for f in (FAM / "cache").glob("seqmap_*.tsv"):
        for r in pd.read_csv(f, sep="\t", dtype=str, keep_default_na=False).to_dict("records"):
            seqm[r["query_id"]] = r

    rep = []
    for r in p.to_dict("records"):
        code, grp = r["species"], group_of(r["section"], r["annotation"])
        loc = (f"{XLSX_REL} (Figure 3-source data 1) row {r['xlsx_row']}: section {r['section']}; {r['annotation']}; "
               f"{r['accession']}")
        base = dict(original_id=r["annotation"], species=code, label=grp or f"({r['section']}:{r['annotation']})", locator=loc)
        if code == "GIAIC" and "GL50581" in r["accession"]:
            rep.append({**base, "gene_id": "", "status": "not_in_panel", "reason": "Giardia assemblage B (GS, GL50581) - panel GIAIC is WB (GL50803)"}); continue
        cands, how, seqinfo = set(), [], None
        for tok, typ in id_tokens(r["accession"]):
            if typ:
                h = (idmap.get(tok, set()) | idmap.get(tok.split(".")[0], set())) & accs[code]
                if h:
                    cands |= h; how.append(f"UniProt ID mapping of {tok}")
                elif tok in seqm:
                    s = seqm[tok]
                    ok = s["status"] == "exact" or (s["status"] == "unique" and s["mapping_identity"] and
                                                    float(s["mapping_identity"]) >= 99.5 and float(s["query_cov"] or 0) >= 0.8
                                                    and float(s["target_cov"] or 0) >= 0.8)
                    if ok and s["target_id"] in accs[code]:
                        cands.add(s["target_id"]); seqinfo = s
                        how.append(f"same-species sequence of {tok}: {s['status']} {s['target_id']} {s['mapping_identity']}% "
                                   f"(cov q{s['query_cov']}/t{s['target_cov']})")
                    else:
                        how.append(f"sequence of {tok}: {s['status'] or 'n/a'} {s['target_id']} {s['mapping_identity']}")
            else:
                for nm in name_candidates(tok, code, r['accession']):
                    h = names[code].get(nm, set())
                    if h:
                        L = r.get("length", "")
                        bad = [a for a in h if L.isdigit() and abs(lens[code].get(a, 0) - int(L)) > 0.2 * int(L)]
                        if bad:
                            how.append(f"gene name {nm} -> {','.join(sorted(h))} rejected: length {lens[code].get(bad[0])} vs table {L} (>20% diff)")
                        else:
                            cands |= h; how.append(f"gene name {nm}" + (f" (length {lens[code].get(next(iter(h)))} vs table {L})" if L else ""))
        if len(cands) == 1:
            g = next(iter(cands))
            meth = "sequence" if seqinfo and not any(x.startswith(("UniProt", "gene")) for x in how) else "id_mapping"
            mi = seqinfo["mapping_identity"] if meth == "sequence" else ""
            mc = (f"{min(float(seqinfo['query_cov'] or 1), float(seqinfo['target_cov'] or 1)):.2f}" if meth == "sequence" else "")
            rep.append({**base, "gene_id": g, "status": "mapped", "reason": "; ".join(how),
                        "_m": dict(mapping_method=meth, mapping_identity=mi if mi else ("100.0" if meth == "sequence" else ""),
                                   mapping_coverage=mc)})
        elif len(cands) > 1:
            rep.append({**base, "gene_id": "", "status": "ambiguous", "reason": "; ".join(how) + " -> " + ",".join(sorted(cands))})
        else:
            rep.append({**base, "gene_id": "", "status": "unmapped", "reason": ("; ".join(how) or "no mappable identifier (JGI/EST id)")})
    byg = defaultdict(list)
    for i, r in enumerate(rep):
        if r["status"] == "mapped":
            byg[r["gene_id"]].append(i)
    for g, ii in byg.items():
        if len(ii) > 1:
            for i in ii:
                rep[i].update(status="ambiguous", reason=rep[i]["reason"] + f"; {g} also mapped by another row")
                rep[i].pop("_m", None)
    mem = [dict(gene_id=r["gene_id"], species=r["species"], original_id_in_paper=r["original_id"], label=r["label"],
                **r["_m"], source=SRC, locator=r["locator"]) for r in rep if r["status"] == "mapped"]
    members = pd.DataFrame(mem).sort_values(["label", "species", "gene_id"]).reset_index(drop=True)
    members.to_csv(FAM / "members.tsv", sep="\t", index=False, quoting=csv.QUOTE_MINIMAL)
    pd.DataFrame([{k: r.get(k, "") for k in ("original_id", "species", "label", "gene_id", "status", "reason", "locator")}
                  for r in rep]).to_csv(FAM / "mapping_report.tsv", sep="\t", index=False)

    # ------------------------------------------------------------------- tree
    ids = defaultdict(int)

    def mk(name, ev, kids, tier="T2", **feat):
        kids = [k for k in kids if k is not None]
        if not kids:
            return None
        if len(kids) == 1:
            return kids[0]
        base = re.sub(r"[^A-Za-z0-9_]", "_", name)
        ids[base] += 1
        i = base if ids[base] == 1 else f"{base}_{ids[base]}"
        nd = Node(i, features={"ev": ev, "id": i, "tier": tier, "src": SRC, **feat})
        for k in kids:
            nd.add_child(k)
        return nd

    def label_group(name, genes, force_u=False):
        by = defaultdict(list)
        for g, s in genes:
            by[s].append(g)
        if force_u:
            return mk(f"{name}_grp", "U", [Node(g) for g, _ in sorted(genes)])
        s_kids = [Node(v[0]) for s, v in sorted(by.items()) if len(v) == 1]
        u_kids = [Node(g) for s, v in sorted(by.items()) if len(v) > 1 for g in sorted(v)]
        s_node = mk(f"{name}_S", "S", s_kids) if len(s_kids) > 1 else (s_kids[0] if s_kids else None)
        return mk(f"{name}_grp", "U", [s_node] + u_kids) if u_kids else s_node

    groups, loose = defaultdict(list), []
    for m in members.itertuples(index=False):
        if m.label.startswith("("):
            loose.append(m.gene_id)
        else:
            groups[m.label].append((m.gene_id, m.species))
    kids = [label_group(g, v, force_u=(g == "AP1/2-beta")) for g, v in sorted(groups.items())]
    leca = mk(f"{FID.split('_')[0]}_LECA", "D", kids, tier="T1", lvl_lo="Eukaryota", lvl_hi="ABOVE_RANGE")
    root = mk(f"{FID.split('_')[0]}_root", "U", [leca] + [Node(g) for g in sorted(loose)]) if loose else leca
    (FAM / "reference_tree.nhx").write_text(
        write_newick(root, features=["ev", "lvl", "lvl_lo", "lvl_hi", "id", "tier", "src", "tag"], lengths=False) + "\n")

    # ------------------------------------------------------------------- assertions
    rows = []
    for k, qk in enumerate(CFG["leca_q"], 1):
        rows.append(dict(assertion_id=f"{FID}_A{k}", type="event", target=leca.name, tier="T1", source_id=SRC,
                         locator=f"{TXT_REL} (Introduction / Results / Figure 4 legend)", quote=Q[qk]))
    extra = list(CFG["member_q"]) + [("tag", "DEEP (pre-LECA subunit paralogs)", CFG["leca_q"][0]),
                                     ("tag", "LOSS (TSET lost in many lineages)", "LOSS"),
                                     ("tag", "AP1/2-beta pooled as U group (beta1/beta2 labels; no AP-1/AP-2 beta split asserted)", "APLECA")]
    for typ, tgt, qk in extra:
        loc = {"F1D1": F1D1, "F1D2": F1D2, "BCOP": F1D1, "GCOP": F1D1}.get(qk, TXT_REL)
        rows.append(dict(assertion_id=f"{FID}_A{len(rows) + 1}", type=typ, target=tgt, tier="T2" if qk not in ("LOSS",) else "T1",
                         source_id=SRC, locator=loc, quote=Q[qk]))
    with open(FAM / "assertions.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader(); w.writerows(rows)

    # ------------------------------------------------------------------- family.yaml
    rp = pd.read_csv(FAM / "mapping_report.tsv", sep="\t", dtype=str, keep_default_na=False)
    inp = rp[rp.status.isin(["mapped", "unmapped", "ambiguous"])]
    frac = inp.groupby("species").status.apply(lambda s: (s == "mapped").mean())
    cnt = defaultdict(int)

    def walk(nd):
        if nd.children:
            cnt[f"{nd.features.get('ev')}_{nd.features.get('tier')}"] += 1
            for c in nd.children:
                walk(c)
    walk(root)
    meta = dict(
        id=FID, name=CFG["name"], panel=PANEL, collection="E_adaptor_complexes",
        phenomena=["DEEP", "LOSS"],
        descriptors=dict(n_members=int(len(members)), n_species=int(members.species.nunique()),
                         mapped_fraction=round(float((inp.status == "mapped").mean()), 3) if len(inp) else 0.0,
                         n_rows_panel_species=int(len(rp)),
                         status_counts={k: int(v) for k, v in rp.status.value_counts().items()},
                         per_species_mapped_fraction={k: round(float(v), 3) for k, v in frac.items()},
                         groups={k: len(v) for k, v in sorted(groups.items())}),
        sources=[dict(id=SRC, doi="10.7554/eLife.02866", pmid="24867644",
                      title="Characterization of TSET, an ancient and widespread membrane trafficking complex",
                      year=2014, journal="eLife", source_method="phylogeny",
                      table="Figure 3-source data 1", locator="Results; Figure 1, 3, 4; Figure 1-source data 1-2; Figure 3-source data 1",
                      quote=Q[CFG["leca_q"][0]],
                      files=[dict(path=x, sha256=sha(ROOT / x)) for x in (XLSX_REL, F1D1, F1D2, TXT_REL)])],
        root_level="Eukaryota", root_single_copy=False,
        complete_species=[], decoys=[], distinct_from=[],
        tier_summary=dict(sorted(cnt.items())),
        protocol_version=1,
        notes=("One of three E07 families split from Hirst et al. 2014 by subunit homology class (large; medium+small; "
               "COPI/TSET scaffold), because the classes are not asserted homologous to each other. Tree: root = T1 D "
               "(lvl_lo=Eukaryota, lvl_hi=ABOVE_RANGE) over subunit groups (section of Figure 3-source data 1; AP "
               "subunits from the annotation code <complex><B|M|S|G|A|D|E|Z>), all pre-LECA paralogs per the quoted "
               "AP-1..5 / COPI / TSET statements. AP-1 and AP-2 beta labels (1B/2B) pooled into one U group: the paper "
               "does not resolve a beta1/beta2 split. Groups: <=1 gene/species -> S, else U. Rows whose annotation "
               "does not parse consistently with its section are loose under a U root. Table lists sequences used for "
               "phylogenies (not exhaustive for AP/COPI; opisthokont AP/COPI sampling sparse) -> complete_species "
               "empty. Mapping: RefSeq/GenBank/UniProt/SGD -> UniProt ID mapping; else exact gene-model name in the "
               "panel header map (AGI, DDB_G, BDEG_, CAOG_, Cre, NCU, CYME_, GL50803_); else NCBI sequence of the "
               "RefSeq id with an exact (or unique >=99.5%, cov >=0.8 both) same-species match. JGI numeric/EST ids "
               "are unmapped. DICDI TSET genes (XP_639969, XP_640471, XP_629998) are absent from the panel proteome."),
    )
    (FAM / "family.yaml").write_text(yaml.safe_dump(meta, sort_keys=False, allow_unicode=True, width=120))


if __name__ == "__main__":
    main()
