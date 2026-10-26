#!/usr/bin/env python
"""E07 adaptor-complex families (E07a large, E07b medium/small, E07c scaffold) - mapping preparation (network).

Run from project root:  python families/E07a_adaptins_large/prep_mapping.py
Source: Hirst J, Schlacht A, Norcott JP, et al. (2014) Characterization of TSET, an ancient and widespread membrane
trafficking complex. eLife 3:e02866 (PMID 24867644). Figure 3-source data 1 (elife-02866-fig3-data1.xlsx: species,
annotation, accession, length of every TSET / COPI / adaptin sequence used in the phylogenies).
RefSeq / GenBank / UniProt / SGD accessions of panel-species rows -> UniProt (curation_tools.map_ids). The result is
written to cache/idmap.tsv of each of the three E07 family directories (their build.py files are offline).
"""
import io, re, shutil, subprocess, sys
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
FAM = Path(__file__).resolve().parent
XLSX_REL = "literature/papers/PMID24867644/raw/supp/elife-02866-fig3-data1.xlsx"
SPMAP = {"homo sapiens": "HUMAN", "drosophila melanogaster": "DROME", "nematostella vectensis": "NEMVE",
         "monosiga brevicollis": "MONBE", "capsaspora owczarzaki": "CAPO3", "capsaspora woczarzaki": "CAPO3",
         "saccharomyces cerevisiae": "YEAST", "sacchromyces cerevisiae": "YEAST", "neurospora crassa": "NEUCR",
         "batrachochytrium dendrobatidis": "BATDJ", "batrachocytrium dendrobabtidis": "BATDJ",
         "encephalitozoon cuniculi": "ENCCU", "dictyostelium discoideum": "DICDI", "entamoeba histolytica": "ENTH1",
         "enthamoeba histolytica": "ENTH1", "arabidopsis thaliana": "ARATH", "chlamydomonas reinhardtii": "CHLRE",
         "cyanidioschyzon merolae": "CYAME", "thalassiosira pseudonana": "THAPS", "plasmodium falciparum": "PLAF7",
         "trypanosoma brucei": "TRYB2", "naegleria gruberi": "NAEGR", "nagleria gruberi": "NAEGR",
         "giardia intestinalis": "GIAIC", "thecamonas trahens": "THETR"}
SECTIONS = {"TPLATE", "TSAUCER", "TCUP", "TSPOON", "TTRAY1", "TTRAY2", "COPI GAMMA", "COPI ZETA", "COPI ALPHA",
            "COPI BETA'", "COPI BETA", "COPI DELTA", "AP Beta", "AP Mu", "AP GADEZ", "AP Sigma"}


def load_table():
    x = pd.read_excel(ROOT / XLSX_REL, header=None, dtype=str)
    rows, sec = [], None
    for i, r in x.iterrows():
        a = r[0].strip() if isinstance(r[0], str) else ""
        if isinstance(r[1], str) and r[1].strip() and a != "Species":
            rows.append(dict(section=sec, species_name=a, annotation=r[1].strip(),
                             accession=(r[2].strip() if isinstance(r[2], str) else ""),
                             length=(r[3].strip() if isinstance(r[3], str) else ""), xlsx_row=i + 1))
        elif a and not a.startswith("http") and a in SECTIONS | {"T-SET", "COPI", "Adaptins"}:
            sec = a if a in SECTIONS else sec
    t = pd.DataFrame(rows)
    t["species"] = t.species_name.str.lower().str.strip().map(SPMAP).fillna("")
    return t


def id_tokens(acc):
    """split an accession cell into candidate IDs with a map_ids type (or None = offline gene-name lookup)."""
    out = []
    for tok in re.split(r"[/|,;\s]+", acc):
        tok = tok.strip()
        if not tok:
            continue
        if re.match(r"^[NX]P_\d+(\.\d+)?$", tok):
            out.append((tok, "refseq_protein"))
        elif re.match(r"^[OPQ]\d[A-Z0-9]{3}\d$|^[A-NR-Z]\d[A-Z][A-Z0-9]{2}\d$", tok):
            out.append((tok, "uniprot"))
        elif re.match(r"^S\d{9}$", tok.replace("SGDID:", "")) and "SGDID" in acc:
            out.append((tok.replace("SGDID:", ""), "sgd_id"))
        elif re.match(r"^[A-Z]{3}\d{5}(\.\d+)?$", tok):
            out.append((tok, "genbank_protein"))
        else:
            out.append((tok, None))
    return out


def main():
    t = load_table()
    sp = pd.read_csv(ROOT / "panels/E/species.tsv", sep="\t", dtype=str).set_index("code")
    ids = []
    for r in t[t.species != ""].to_dict("records"):
        for tok, typ in id_tokens(r["accession"]):
            if typ:
                ids.append((tok, typ, sp.loc[r["species"], "taxid"]))
    cache = FAM / "cache"
    cache.mkdir(exist_ok=True)
    inp = cache / "idmap_input.tsv"
    pd.DataFrame(ids, columns=["id", "id_type", "taxon"]).drop_duplicates().to_csv(inp, sep="\t", index=False)
    if not (cache / "idmap.tsv").exists():
        subprocess.run([sys.executable, "-m", "curation_tools.map_ids", str(inp), "-o", str(cache / "idmap.tsv")],
                       check=True, cwd=ROOT)
    # RefSeq/GenBank ids whose UniProt mapping misses the panel proteome -> fetch sequence, exact same-species match
    idm = pd.read_csv(cache / "idmap.tsv", sep="\t", dtype=str, keep_default_na=False)
    for code, g in t[t.species != ""].groupby("species"):
        accs = {l[1:].split()[0] for l in open(ROOT / f"panels/E/proteomes/{code}.fa") if l.startswith(">")}
        need = []
        for r in g.to_dict("records"):
            for tok, typ in id_tokens(r["accession"]):
                if typ in ("refseq_protein", "genbank_protein"):
                    hit = set(idm[idm.original_id.isin([tok, tok.split(".")[0]])].uniprot_acc) & accs
                    if not hit:
                        need.append(tok)
        out = cache / f"seqmap_{code}.tsv"
        if need and not out.exists():
            q = cache / f"seqmap_{code}_ids.txt"
            q.write_text("\n".join(sorted(set(need))) + "\n")
            subprocess.run([sys.executable, "-m", "curation_tools.map_by_sequence", "--ids", str(q), "--target",
                            str(ROOT / f"panels/E/proteomes/{code}.fa"), "-o", str(out), "--ambiguity-delta", "2.0",
                            "--min-qcov", "0.8", "--min-tcov", "0.8"], check=True, cwd=ROOT)
    for other in ("E07b_adaptins_medium_small", "E07c_adaptins_scaffold"):
        d = ROOT / "families" / other / "cache"
        d.mkdir(parents=True, exist_ok=True)
        for f in cache.iterdir():
            shutil.copy(f, d / f.name)


if __name__ == "__main__":
    main()
