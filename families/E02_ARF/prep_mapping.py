#!/usr/bin/env python
"""E02 ARF family - mapping preparation (network + DIAMOND); writes cache/ files read by build.py.

Run from project root:  python families/E02_ARF/prep_mapping.py
Source: Vargova et al. 2021 GBE 13:evab157 (PMID 34247240), Supplementary table 1 (ARF family sequences, annotation,
accessions) and Supplementary dataset 1 (FASTA of analysed protein sequences).
1. RefSeq/GenBank accessions of panel-species rows -> UniProt (curation_tools.map_ids)        -> cache/idmap.tsv
2. paper protein sequences (dataset 1, else corrected sequence in table 1) -> DIAMOND vs the same species' panel
   proteome, plus reverse search for the reciprocal-best-hit check                              -> cache/fw.tsv, rv.tsv
"""
import subprocess, sys, zipfile, io, re
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
FAM = Path(__file__).resolve().parent
CACHE = FAM / "cache"
sys.path.insert(0, str(FAM))
from seqmap import run_diamond  # noqa: E402

ZIP_REL = "literature/papers/PMID34247240/raw/supp/evab157_supplementary_data.zip"
# species header in table 1 (prefix) -> panel code
SPMAP = {"Homo sapiens": "HUMAN", "Drosophila melanogaster": "DROME", "Nematostella vectensis": "NEMVE",
         "Trichoplax adhaerens": "TRIAD", "Mnemiopsis leidyi": "MNELE", "Amphimedon queenslandica": "AMPQE",
         "Monosiga brevicollis": "MONBE", "Capsaspora owczarzaki": "CAPO3", "Saccharomyces cerevisiae": "YEAST",
         "Neurospora crassa": "NEUCR", "Batrachochytrium dendrobatidis": "BATDJ", "Dictyostelium discoideum": "DICDI",
         "Entamoeba histolytica": "ENTH1", "Thecamonas trahens": "THETR", "Giardia intestinalis": "GIAIC",
         "Trypanosoma brucei": "TRYB2", "Naegleria gruberi": "NAEGR", "Arabidopsis thaliana": "ARATH",
         "Chlamydomonas reinhardtii": "CHLRE", "Cyanidioschyzon merolae": "CYAME", "Plasmodium falciparum": "PLAF7",
         "Thalassiosira pseudonana": "THAPS"}
COLS = ["species_header", "Abbreviation", "Annotation", "Scrollsaw", "Annotation basis", "Conclusively annotated",
        "S4", "S5", "S6", "S11", "Notes", "Number of genes", "Accession", "Corr_model", "No_model", "Corr_genreads",
        "Corr_transcr", "Corrected protein sequence", "Additional domain", "Myr", "TM", "Palm"]


def _zip_bytes(name):
    with zipfile.ZipFile(ROOT / ZIP_REL) as z:
        return z.read(name)


def load_table():
    x = pd.read_excel(io.BytesIO(_zip_bytes("Supplementary_table_1.xlsx")), header=None, dtype=str)
    x.columns = COLS
    rows, sp = [], None
    for i, r in x.iterrows():
        if i < 3:
            continue
        if isinstance(r.species_header, str) and r.species_header.strip() and not isinstance(r.Abbreviation, str):
            sp = r.species_header.strip()
            continue
        if isinstance(r.Abbreviation, str) and r.Abbreviation.strip():
            d = {k: (v.strip() if isinstance(v, str) else "") for k, v in r.items()}
            d["species_name"] = sp
            d["xlsx_row"] = i + 1  # 1-based spreadsheet row
            rows.append(d)
    t = pd.DataFrame(rows)
    t["species"] = [next((c for k, c in SPMAP.items() if s.startswith(k)), "") for s in t.species_name]
    return t


def load_dataset():
    seqs, name = {}, None
    for l in _zip_bytes("Supplementary_dataset_1.txt").decode("utf-8", "replace").splitlines():
        l = l.strip()
        if l.startswith(">"):
            name = l[1:].split()[0]
            seqs[name] = []
        elif name:
            seqs[name].append(re.sub(r"[^A-Za-z]", "", l).upper())
    return {k: "".join(v) for k, v in seqs.items()}


def row_seq(r, ds):
    s = ds.get(r["Abbreviation"], "")
    if not s and r["Corrected protein sequence"]:
        s = re.sub(r"[^A-Za-z]", "", r["Corrected protein sequence"]).upper()
    return s.replace("*", "")


def accessions(r):
    return [a.strip() for a in re.split(r"[;,\s]+", r["Accession"]) if re.match(r"^(NP_|XP_|[A-Z]{3}\d{5})[\d.]*", a.strip())]


def main():
    CACHE.mkdir(exist_ok=True)
    t = load_table()
    ds = load_dataset()
    sp = pd.read_csv(ROOT / "panels/E/species.tsv", sep="\t", dtype=str).set_index("code")
    p = t[t.species != ""]
    ids = []
    for r in p.to_dict("records"):
        for a in accessions(r):
            ids.append((a, "refseq_protein" if a[:3] in ("NP_", "XP_") else "genbank_protein", sp.loc[r["species"], "taxid"]))
    inp = CACHE / "idmap_input.tsv"
    pd.DataFrame(ids, columns=["id", "id_type", "taxon"]).drop_duplicates().to_csv(inp, sep="\t", index=False)
    if not (CACHE / "idmap.tsv").exists():
        subprocess.run([sys.executable, "-m", "curation_tools.map_ids", str(inp), "-o", str(CACHE / "idmap.tsv")],
                       check=True, cwd=ROOT)
    fws, rvs = [], []
    for code, g in p.groupby("species"):
        allsp = {r["Abbreviation"]: row_seq(r, ds) for r in t[t.species_name == g.species_name.iloc[0]].to_dict("records")}
        q = {r["Abbreviation"]: row_seq(r, ds) for r in g.to_dict("records")}
        fw, rv = run_diamond(q, code, allsp)
        fws.append(fw); rvs.append(rv)
    pd.concat(fws).to_csv(CACHE / "fw.tsv", sep="\t", index=False)
    pd.concat(rvs).to_csv(CACHE / "rv.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
