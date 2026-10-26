#!/usr/bin/env python
"""E03_Rab: build members.tsv, mapping_report.tsv, reference_tree.nhx, assertions.tsv, family.yaml.

Source: Kloepper TH, Kienle N, Fasshauer D, Munro S (2012) Untangling the evolution of Rab G proteins: implications of a
comprehensive genomic analysis. BMC Biol 10:71 (PMID 22873208). Additional files 4-9: alignments of Rab supergroups
I-VI with sequence names <Species>_<Rab label>; Additional file 3: species short names. No accessions -> same-species
sequence mapping (cache/ from prep_mapping.py; this script is offline).
Run from the project root:  python families/E03_Rab/build.py
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
from prep_mapping import load_alignments, SPMAP, SUPP, FILES  # noqa: E402
from seqmap import best_hits, ID_MIN, COV_MIN  # noqa: E402
from hobench.newick import Node, write_newick  # noqa: E402
from hobench.speciestree import SpeciesTree  # noqa: E402

FID, SRC, PANEL = "E03_Rab", "Kloepper2012", "E"
TXT_REL = "literature/papers/PMID22873208/fulltext.txt"
Q = dict(
    LECA=("Of these proteins, Rab1, Rab2, Rab4, Rab5, Rab6, Rab7, Rab8, Rab11, Rab18, Rab21, Rab23, and Rab28 are known "
          "likely candidates, and Rab14, Rab32, and RabL4 have just recently been added to this list [16]. To these "
          "proteins we can now add Rab22, Rab24, Rab29, RabX1, and Rab7L1."),
    LECA2=("Our analysis identifies a set of 20 Rab proteins that are likely to have been present in the LECA based on "
           "the arguments above"),
    SG=("it became apparent that the twenty LECA Rabs fell into a set of six larger supergroups, suggesting a "
        "primitive pre-LECA eukaryote with just six Rabs."),
    G1="Group I: Rab1, Rab8, Rab18",
    G2="Group II: Rab5, Rab21, Rab22, Rab24 and RabX1",
    G3="Group III: Rab7, Rab23, Rab29, Rab32 and Rab7L1",
    G4="Group IV: Rab2, Rab4, Rab11 and Rab14",
    G5="Group V: Rab6",
    G6="Group IV: Rab28 and RabL4",
    FILE1="Additional file 4 Alignment of group I Rab proteins in FASTA format.",
    RAB1=("Rab1 is the 'indispensible' member of the supergroup, and seems to have been duplicated a number of times "
          "in metazoans (Rab19, Rab30, Rab33, Rab35, RabX6)."),
    RAB8=("More strikingly, it has a large set of duplications in metazoans (Rab3, Rab10, Rab15, Rab26, Rab27, Rab34, "
          "Rab44, Rab45, RabX4)."),
    RAB8V=("In addition, Rab8 seems to have undergone another duplication in deuterostomes (Rab12) and two in "
           "vertebrates (Rab8b, Rab13)."),
    RAB18=("It shows a duplication in angiosperms and also duplication in bilateria (Rab40), which expanded further in "
           "vertebrates (Rab40b, Rab40c) and then in primates (Rab40a, Rab40aL)."),
    UNLINKED=("To resolve the five metazoan groups that we could not link unambiguously to one LECA Rab (Rab19, Rab33, "
              "Rab34, Rab44 and Rab45), we used a pBLAST-based similarity search"),
    RAB5=("Rab5 is apparently indispensible, and also has the most complex evolutionary history of this supergroup. It "
          "seems to have been duplicated independently in basal fungi after the loss of their flagellum (Ypt52)"),
    RAB5V="In addition, it has quadrupled in vertebrates (Rab5a, Rab5b, Rab5c, Rab17).",
    RAB5Y="Similarly, there are independent duplications in Saccharomycotina (Ypt10)",
    RAB22V="Rab21 and Rab22 have no major expansions (Rab22 has a duplication in vertebrates)",
    RAB24="Rab24 is rather unusual as it is well conserved in metazoans, where it has also been duplicated (Rab20)",
    RAB9=("However, we identified a conserved set of fungal proteins that must either be an independent duplication of "
          "fungal Rab7 or must belong to the Rab9 subfamily"),
    RAB9b="and hence Rab9 appeared with the opisthokonts",
    RAB7="Rab7 seems to be indispensible, and is also the member of the family that shows the most expansions.",
    RAB32V="Rab32 has been duplicated in vertebrates (Rab38)",
    RAB11Y=("Rab11 seems to be the indispensible member of this supergroup, and is certainly the only Rab from this "
            "family that is still present in Saccharomycotina, where it has duplicated (Ypt31/32)."),
    RAB2="Rab2 is present in all major eukaryotic phyla and has been independently duplicated in metazoa (Rab39)",
    RAB4V="By contrast, Rab4 and Rab14 have no major duplications (Rab4 has been duplicated in vertebrates)",
    RAB6V=("an expansion in vertebrates (Rab6b). In the latter case, there was a later expansion to generate Rab41, "
           "which is present in primates and dolphins"),
    RAB6C="which seems to be specific to Hominidae (Rab6c, [53])",
    METHOD=("We defined a set of 20 different subgroups, each subgroup representing a Rab present in the LECA. Because "
            "these subgroups have undergone extensive duplications within the metazoans we supplemented these 20 "
            "subgroups with 15 further subgroups representing those metazoan duplications."),
    LOSS=("Tracing the fate of the LECA Rabs revealed extensive losses with many extant eukaryotes having fewer Rabs, "
          "and none having the full complement."),
    LSE=("Metazoans show expansions in twelve different LECA Rabs, plants in eight, heterokonts in five, apicomplexans "
         "in three, fungi in two and kinetoplastids in just one"),
    RAS=("This may indicate that these groups are not classic Rab proteins, but would be better placed somewhere "
         "between the Ran and Rab supergroups of the Ras family."),
)
# supergroup -> LECA Rab -> labels assigned to it (quoted in Q), with the quote keys supporting the assignment
LECA = {
    "I": {"Rab1": (["Rab1", "Rab30", "Rab35", "RabX6"], ["RAB1"]),
          "Rab8": (["Rab8", "Rab3", "Rab10", "Rab15", "Rab26", "Rab27", "RabX4", "Rab12", "Rab13"], ["RAB8", "RAB8V"]),
          "Rab18": (["Rab18", "Rab40"], ["RAB18"])},
    "II": {"Rab5": (["Rab5", "Rab17", "Ypt52", "Ypt10"], ["RAB5", "RAB5V", "RAB5Y"]), "Rab21": (["Rab21"], []),
           "Rab22": (["Rab22"], []), "Rab24": (["Rab24", "Rab20"], ["RAB24"]), "RabX1": (["RabX1"], [])},
    "III": {"Rab7": (["Rab7", "Rab9"], ["RAB9b"]), "Rab23": (["Rab23"], []), "Rab29": (["Rab29"], []),
            "Rab32": (["Rab32", "Rab38"], ["RAB32V"]), "Rab7L1": (["Rab7L1"], [])},
    "IV": {"Rab2": (["Rab2", "Rab39"], ["RAB2"]), "Rab4": (["Rab4"], []), "Rab11": (["Rab11", "Ypt31"], ["RAB11Y"]),
           "Rab14": (["Rab14"], [])},
    "V": {"Rab6": (["Rab6", "Rab41"], ["RAB6V"])},
    "VI": {"Rab28": (["Rab28"], []), "RabL4": (["RabL4"], [])},
}
GQ = {"I": "G1", "II": "G2", "III": "G3", "IV": "G4", "V": "G5", "VI": "G6"}
# quoted vertebrate (HUMAN-branch) duplications: LECA Rab -> (HUMAN gene labels, quote keys)
HUMAN_D = {"Rab5": ({"Rab5a", "Rab5b", "Rab5c", "Rab17"}, ["RAB5V"]),
           "Rab8": ({"Rab8a", "Rab8b", "Rab13", "Rab12"}, ["RAB8V"]),
           "Rab22": ({"Rab22a", "Rab22b"}, ["RAB22V"]),
           "Rab32": ({"Rab32", "Rab38"}, ["RAB32V"]),
           "Rab4": ({"Rab4a", "Rab4b"}, ["RAB4V"]),
           "Rab6": ({"Rab6a", "Rab6b", "Rab6c", "Rab41"}, ["RAB6V", "RAB6C"])}
LBL = re.compile(r"^(Rab7L1|RabX\d+|RabL\d+|Rab\d+|Ypt\d+|Sec4)(.*)$")


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def parse_label(name):
    """'AuAn_Rab8B_2' -> ('Rab8B', 'Rab8'); composite/unknown -> (rest, '')."""
    rest = name.split("_", 1)[1]
    m = LBL.match(rest)
    if not m:
        return rest, ""
    base, suf = m.groups()
    if re.match(r"^(_\d+)+(_\d+)*$", suf) and suf.count("_") > 1:   # e.g. Rab3_8_10_13
        return rest, ""
    if re.match(r"^_\d+(_\d+)+", suf):
        return rest, ""
    lab = re.sub(r"_\d+$", "", rest)
    return lab, base


def main():
    st = SpeciesTree.from_newick(ROOT / f"panels/{PANEL}/species_tree.nwk")
    a = load_alignments()
    a = a[a.sp.isin(SPMAP)].copy()
    a["species"] = a.sp.map(SPMAP)
    fw = pd.read_csv(FAM / "cache/fw.tsv", sep="\t")
    rv = pd.read_csv(FAM / "cache/rv.tsv", sep="\t")
    hits = {c: best_hits(fw[fw.species == c], rv[rv.species == c]) for c in a.species.unique()}
    rep = []
    for r in a.to_dict("records"):
        lab, base = parse_label(r["name"])
        loc = f"{SUPP}/{r['file']} (Additional file, supergroup {r['group']}) record {r['order']}: >{r['name']}"
        b = hits[r["species"]].get(r["name"])
        d = dict(original_id=r["name"], species=r["species"], label=lab, supergroup=r["group"], base=base, group=r["group"],
                 locator=loc)
        if b is None:
            rep.append({**d, "gene_id": "", "status": "unmapped", "reason": "no DIAMOND hit in panel proteome"}); continue
        note = (f"same-species sequence: best {b['t']} {b['pid']:.1f}% id over aligned residues, cov q{b['qcov']:.0%}/"
                f"t{b['tcov']:.0%}, margin {b['margin']:.1f}, rbh {b['rbh']}")
        if b["ok"]:
            rep.append({**d, "gene_id": b["t"], "status": "mapped", "reason": note,
                        "_m": dict(mapping_method="sequence", mapping_identity=f"{b['pid']:.1f}",
                                   mapping_coverage=f"{min(b['qcov'], b['tcov']):.2f}")})
        elif b["pid"] >= ID_MIN and b["qcov"] >= COV_MIN and b["tcov"] >= COV_MIN:
            rep.append({**d, "gene_id": "", "status": "ambiguous", "reason": note + " (fails RBH/margin)"})
        else:
            rep.append({**d, "gene_id": "", "status": "unmapped", "reason": note + " (below thresholds)"})
    byg = defaultdict(list)
    for i, r in enumerate(rep):
        if r["status"] == "mapped":
            byg[r["gene_id"]].append(i)
    for g, ii in byg.items():
        if len(ii) > 1:
            for i in ii:
                rep[i].update(status="ambiguous", reason=rep[i]["reason"] + f"; {g} also hit by another paper sequence")
                rep[i].pop("_m", None)
    mem = [dict(gene_id=r["gene_id"], species=r["species"], original_id_in_paper=r["original_id"], label=r["label"],
                supergroup=r["supergroup"],
                **r["_m"], source=SRC, locator=r["locator"]) for r in rep if r["status"] == "mapped"]
    members = pd.DataFrame(mem).sort_values(["label", "species", "gene_id"]).reset_index(drop=True)
    members.to_csv(FAM / "members.tsv", sep="\t", index=False, quoting=csv.QUOTE_MINIMAL)
    pd.DataFrame([{k: r.get(k, "") for k in ("original_id", "species", "supergroup", "label", "gene_id", "status", "reason", "locator")}
                  for r in rep]).to_csv(FAM / "mapping_report.tsv", sep="\t", index=False)

    # ------------------------------------------------------------------- tree
    ids = defaultdict(int)
    asserted = []

    def nid(base):
        base = re.sub(r"[^A-Za-z0-9_]", "_", base)
        ids[base] += 1
        return base if ids[base] == 1 else f"{base}_{ids[base]}"

    def mk(name, ev, kids, tier="T2", **feat):
        kids = [k for k in kids if k is not None]
        if not kids:
            return None
        if len(kids) == 1:
            return kids[0]
        i = nid(name)
        nd = Node(i, features={"ev": ev, "id": i, "tier": tier, "src": SRC, **feat})
        for k in kids:
            nd.add_child(k)
        return nd

    def label_group(name, genes, human_d=None):
        """genes: [(gene, species, human_label)]; <=1/species -> S, else U; HUMAN copies under a quoted D if all are
        named in human_d."""
        if not genes:
            return None
        by = defaultdict(list)
        for g, s, hl in genes:
            by[s].append((g, hl))
        s_kids, u_kids = [], []
        for s in sorted(by):
            gs = sorted(by[s])
            if len(gs) == 1:
                s_kids.append(Node(gs[0][0]))
            elif s == "HUMAN" and human_d and all(hl in human_d[0] for _, hl in gs):
                d = mk(f"{name}_HUMAN_dup", "D", [Node(g) for g, _ in gs], tier="T1")
                asserted.extend((d.name, k) for k in human_d[1])
                s_kids.append(d)
            else:
                u_kids += [Node(g) for g, _ in gs]
        s_node = mk(f"{name}_S", "S", s_kids) if len(s_kids) > 1 else (s_kids[0] if s_kids else None)
        return mk(f"{name}_grp", "U", [s_node] + u_kids) if u_kids else s_node

    # assign members
    lab2leca = {}
    for sg, d in LECA.items():
        for leca, (labs, _) in d.items():
            for l in labs:
                lab2leca[l] = (sg, leca)
    pool = defaultdict(lambda: defaultdict(list))  # (sg, leca) -> core/derived label -> genes
    loose = defaultdict(list)                       # sg -> genes
    for m in members.itertuples(index=False):
        sg, lab = m.supergroup, m.label
        _, base = parse_label(m.original_id_in_paper)
        hl = lab if m.species == "HUMAN" else ""
        if base in lab2leca and lab2leca[base][0] == sg and not re.search(r"(partial|miss)", lab):
            key = lab2leca[base]
            hd = HUMAN_D.get(key[1])
            if m.species == "HUMAN" and hd and lab in hd[0]:
                pool[key][key[1]].append((m.gene_id, m.species, hl))   # quoted vertebrate paralog -> core group
            else:
                pool[key][base].append((m.gene_id, m.species, hl))
        else:
            loose[sg].append(m.gene_id)

    opis = set(st.clade_species("Opisthokonta"))
    bila = set(st.clade_species("Bilateria"))

    def leca_node(sg, leca):
        groups = pool.get((sg, leca), {})
        core = groups.get(leca, [])
        hd = HUMAN_D.get(leca)
        if leca == "Rab7" and groups.get("Rab9"):
            non = [x for x in core if x[1] not in opis]
            op = [x for x in core if x[1] in opis]
            d = mk("Rab7_Rab9_Opisthokonta", "D", [label_group("Rab7_opisthokonta", op), label_group("Rab9", groups["Rab9"])],
                   tier="T1", lvl="Opisthokonta")
            if d is not None and d.features.get("ev") == "D":
                asserted.extend([(d.name, "RAB9b"), (d.name, "RAB9")])
            return mk("Rab7_complex", "U", [label_group("Rab7_nonopisthokonta", non), d])
        if leca == "Rab18" and groups.get("Rab40"):
            non = [x for x in core if x[1] not in bila]
            bi = [x for x in core if x[1] in bila]
            d = mk("Rab18_Rab40_Bilateria", "D", [label_group("Rab18_bilateria", bi), label_group("Rab40", groups["Rab40"])],
                   tier="T1", lvl="Bilateria")
            if d is not None and d.features.get("ev") == "D":
                asserted.append((d.name, "RAB18"))
            return mk("Rab18_complex", "U", [label_group("Rab18_nonbilateria", non), d])
        kids = [label_group(leca, core, hd)] + [label_group(l, g) for l, g in sorted(groups.items()) if l != leca]
        return mk(f"{leca}_LECA", "U", kids)

    sg_nodes = []
    for sg, d in LECA.items():
        kids = [leca_node(sg, l) for l in d]
        dd = mk(f"SG{sg}_LECA", "D", kids, tier="T1", lvl_lo="Eukaryota", lvl_hi="ABOVE_RANGE")
        if dd is not None and dd.features.get("ev") == "D":
            asserted.extend([(dd.name, "LECA"), (dd.name, GQ[sg])])
        sg_nodes.append(mk(f"SG{sg}", "U", [dd] + [Node(g) for g in sorted(loose.get(sg, []))]))
    root = mk("RAB_LECA_supergroups", "D", sg_nodes, tier="T1", lvl_lo="Eukaryota", lvl_hi="ABOVE_RANGE")
    asserted.extend([(root.name, "SG"), (root.name, "LECA2")])
    (FAM / "reference_tree.nhx").write_text(
        write_newick(root, features=["ev", "lvl", "lvl_lo", "lvl_hi", "id", "tier", "src", "tag"], lengths=False) + "\n")

    # ------------------------------------------------------------------- assertions
    rows = [dict(assertion_id=f"{FID}_A{k}", type="event", target=n, tier="T1", source_id=SRC,
                 locator=f"{TXT_REL} (Results/Discussion)", quote=Q[q]) for k, (n, q) in enumerate(asserted, 1)]
    n = len(rows)
    extra = [("tag", "DEEP (20 LECA Rabs in 6 supergroups)", "LECA2", "T1"), ("tag", "LOSS", "LOSS", "T1"),
             ("tag", "LSE", "LSE", "T1"),
             ("tag", "T2 rule: subgroup classification (HMM + phylogeny); label = sequence name suffix", "METHOD", "T2"),
             ("tag", "Rab19, Rab33, Rab34, Rab44, Rab45 (and other unassigned labels) left loose in their supergroup U", "UNLINKED", "T2"),
             ("tag", "supergroup of each sequence = alignment file (Additional files 4-9)", "FILE1", "T2"),
             ("tag", "homology with ARF family (E02): Rab and Arf are Ras-superfamily families; no distinct_from", "RAS", "T2")]
    for sg, d in LECA.items():
        for leca, (labs, qs) in d.items():
            for qk in qs:
                extra.append(("tag", f"labels {','.join(labs)} assigned to LECA {leca}", qk, "T2"))
    for typ, tgt, qk, tier in extra:
        n += 1
        rows.append(dict(assertion_id=f"{FID}_A{n}", type=typ, target=tgt, tier=tier, source_id=SRC, locator=TXT_REL, quote=Q[qk]))
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
    present = sorted(set(SPMAP.values()))
    panel_euk = sorted(st.clade_species("Eukaryota"))
    meta = dict(
        id=FID, name="Rab GTPases (eukaryote-wide; 20 LECA Rabs in six supergroups)", panel=PANEL,
        collection="E_small_GTPases", phenomena=["DEEP", "LOSS", "LSE", "PROM"],
        descriptors=dict(n_members=int(len(members)), n_species=int(members.species.nunique()),
                         mapped_fraction=round(float((inp.status == "mapped").mean()), 3),
                         n_rows_panel_species=int(len(rp)),
                         status_counts={k: int(v) for k, v in rp.status.value_counts().items()},
                         per_species_mapped_fraction={k: round(float(v), 3) for k, v in frac.items()}),
        sources=[dict(id=SRC, doi="10.1186/1741-7007-10-71", pmid="22873208",
                      title="Untangling the evolution of Rab G proteins: implications of a comprehensive genomic analysis",
                      year=2012, journal="BMC Biology", source_method="HMM_classification",
                      table="Additional files 4-9 (supergroup alignments)", locator="Results; Figures 1-5; Additional files 3-9",
                      quote=Q["LECA"],
                      files=[dict(path=p, sha256=sha(ROOT / p)) for p in
                             [f"{SUPP}/{f}" for f in FILES.values()] + [f"{SUPP}/1741-7007-10-71-S3.XLS", TXT_REL]])],
        root_level="Eukaryota", root_single_copy=False,
        complete_species=[], decoys=[], distinct_from=[],
        homologous_families=[dict(family="E02_ARF", relation="Ras superfamily (homologous small GTPases); deep merges not penalised", src=SRC)],
        tier_summary=dict(sorted(cnt.items())),
        protocol_version=1,
        notes=("Root = T1 D over the six supergroups (pre-LECA, interval [Eukaryota, ABOVE_RANGE]); each supergroup = U "
               "wrapper over a T1 D splitting its LECA Rabs, plus loose genes whose label is not quote-assigned to one "
               "LECA Rab (Ypt1, Sec4, Ypt7, Ypt51/53, Ypt3, Ypt6, Rab19/33/34/43/44/45, Rab25, Rab36/37, Rab42, RabX2, "
               "composite/partial labels; supergroup known from the alignment file). LECA Rab node = U over label groups "
               "(core label + quoted metazoan-derived labels); undated metazoan splits are not asserted. T1 nested D: "
               "Rab9 vs opisthokont Rab7 at Opisthokonta; Rab40 vs bilaterian Rab18 at Bilateria; vertebrate "
               "duplications (Rab5a/b/c/17, Rab8a/8b/13/12, Rab22a/b, Rab32/38, Rab4a/b, Rab6a/b/c/41) as D over the "
               "HUMAN copies (HUMAN = only vertebrate). Label groups: <=1 gene/species -> S; else U (T2 rules b/c). "
               "Mapping: no accessions; degapped alignment sequences (authors removed gap-rich columns) -> same-species "
               "DIAMOND; identity computed over aligned residue pairs (trimmed columns appear as gaps), coverage of both "
               "sequences >= 80% by aligned span, RBH, margin >= 2. complete_species empty: alignments are HMM training "
               "subsets, no statement that they list every Rab of a species. Source_method HMM_classification "
               "(iterative HMM + tree; not a denylisted tool). Elias et al. 2012 (J Cell Sci, PMID 22366452) fetched "
               "(PMC HTML, supplements not retrievable) but not used for assertions. Panel-E eukaryotes absent from the "
               "source alignments: " + ", ".join(c for c in panel_euk if c not in present) +
               ". Homology: Rab and Arf/Sar are homologous Ras-superfamily GTPases; no distinct_from with E02_ARF."),
    )
    (FAM / "family.yaml").write_text(yaml.safe_dump(meta, sort_keys=False, allow_unicode=True, width=120))


if __name__ == "__main__":
    main()
