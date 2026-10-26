#!/usr/bin/env python3
"""Build the benchmark species panels: panels/<P>/{species.tsv, proteomes/<CODE>.fa, species_tree.nwk,
species_tree_alt.nwk, levels.tsv} for P in Y K V E P (panel D dropped, reviews/round1/SYNTHESIS.md).

Sources (no silent assumptions):
  * species lists: literature/raw/C/panels_verified.tsv rows with tier == core (survey C, section 3), plus the
    explicit EXTRA species below (SYNTHESIS A5/A6, PI decision 2026-10-26);
  * proteomes: data/proteomes/<UPID>_<taxid>/ written by `python -m curation_tools.fetch_proteome` (manifest.json
    gives release, counts, sha256; the sha256 is re-computed here and must match);
  * species trees: hand-written Newick strings below, each internal node named by a level id
    (^[A-Za-z0-9_-]{3,}$), topology and contested nodes cited in LEVELS (DOIs checked in Europe PMC,
    cache panels/citations.tsv; refresh with --verify-dois);
  * NCBI taxid/name per level: computed with ete3 from data/raw/ncbi/taxa.sqlite (built from taxdump.tar.gz)
    as the NCBI LCA of the level's species; `ncbi_exact` = y iff that NCBI taxon contains exactly the level's
    panel species (i.e. the level corresponds to an NCBI taxon restricted to the panel).

Usage: python scripts/build_panels.py [--panels Y K] [--verify-dois] [--no-sha]
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import re
import sys
import time
import urllib.parse
import urllib.request
from collections import Counter
from pathlib import Path

from ete3 import Tree

ROOT = Path(__file__).resolve().parents[1]
PROTEOMES = ROOT / "data" / "proteomes"
PANELS_DIR = ROOT / "panels"
VERIFIED = ROOT / "literature" / "raw" / "C" / "panels_verified.tsv"
TAXDB = ROOT / "data" / "raw" / "ncbi" / "taxa.sqlite"
CITATIONS = PANELS_DIR / "citations.tsv"
PINNED_RELEASE = "2026_03"
LEVEL_RE = re.compile(r"^[A-Za-z0-9_-]{3,}$")

# ---------------------------------------------------------------------------------------------------------------
# Citations (key -> DOI). Every key used in LEVELS must be here; --verify-dois checks each DOI in Europe PMC.
CIT = {
    "Shen2018": "10.1016/j.cell.2018.10.023",
    "Shen2016": "10.1534/g3.116.034744",
    "MarcetHouben2015": "10.1371/journal.pbio.1002220",
    "Wolfe2015": "10.1371/journal.pbio.1002221",
    "WolfeShields1997": "10.1038/42711",
    "Coleman2021": "10.1126/science.abe0511",
    "Hug2016": "10.1038/nmicrobiol.2016.48",
    "Parks2018": "10.1038/nbt.4229",
    "Parks2022GTDB": "10.1093/nar/gkab776",
    "Williams2020": "10.1038/s41559-019-1040-x",
    "Eme2023": "10.1038/s41586-023-06186-2",
    # vertebrates
    "Delsuc2006": "10.1038/nature04336",
    "Irisarri2017": "10.1038/s41559-017-0240-5",
    "Betancur2017": "10.1186/s12862-017-0958-3",
    "Hughes2018": "10.1073/pnas.1719358115",
    "Meyer2021lungfish": "10.1038/s41586-021-03198-8",
    "Simakov2020": "10.1038/s41559-020-1156-z",
    "Nakatani2021": "10.1038/s41467-021-24573-z",
    "Braasch2016": "10.1038/ng.3526",
    "Lien2016": "10.1038/nature17164",
    "Robertson2017": "10.1186/s13059-017-1241-z",
    "Session2016": "10.1038/nature19840",
    "Putnam2008": "10.1038/nature06967",
    # eukaryotes
    "Burki2020": "10.1016/j.tree.2019.08.008",
    "Derelle2015": "10.1073/pnas.1420657112",
    "Whelan2015": "10.1073/pnas.1503453112",
    "Schultz2023": "10.1038/s41586-023-05936-6",
    "Simion2017": "10.1016/j.cub.2017.02.031",
    "Redmond2021": "10.1038/s41467-021-22074-7",
    "Kapli2020": "10.1126/sciadv.abc5162",
    "Laumer2019": "10.1098/rspb.2019.0831",
    "Brown2013": "10.1098/rspb.2013.1755",
    "James2006": "10.1038/nature05110",
    "CapellaGutierrez2012": "10.1186/1741-7007-10-47",
    "Zaremba2017": "10.1038/nature21031",
    # plants
    "OneKP2019": "10.1038/s41586-019-1693-2",
    "Puttick2018": "10.1016/j.cub.2018.01.063",
    "Harris2020": "10.1016/j.cub.2020.03.048",
    "Qiu2006": "10.1073/pnas.0603335103",
    "Amborella2013": "10.1126/science.1241089",
    "Zuntini2024": "10.1038/s41586-024-07324-0",
    "Schmutz2010": "10.1038/nature08670",
    "Wang2011Brassica": "10.1038/ng.919",
    "Jiao2011": "10.1038/nature09916",
}

# ---------------------------------------------------------------------------------------------------------------
# Extra species beyond the survey-C core lists: code -> (upid, why)
EXTRA = {
    "V": {"ONCMY": ("UP000694395", "rainbow trout (Ss4R salmonid; Salmoninae); reference proteome"),
          "ESOLU": ("UP000265140", "northern pike (Esociformes; unduplicated sister lineage of salmonids for Ss4R)"),
          "CORBA": ("UP000502089", "Coregonus sp. 'balchen' (Coregoninae; splits LORe/AORe Ss4R regions, "
                                   "SYNTHESIS A5); NON-reference proteome, UniParc sequences")},
    "P": {"PHAVU": ("UP000000226", "common bean (Phaseoleae; splits the soybean-specific WGD from the "
                                   "older legume WGD, SYNTHESIS A6)")},
}
# species needing the stricter sequence-mapping identity threshold (hobench.panel STRICT_SPECIES + Coregonus)
STRICT = {"XENLA", "SALSA", "ONCMY", "CORBA"}

# Species-specific notes added to the `why` column
NOTES = {
    "CORBA": "UniProt mnemonic 9TELE is not species-specific, code CORBA chosen here; UniParc UPIs as IDs "
             "(41,647 distinct sequences for 41,750 proteins); no gene-level canonical collapsing (QC: OMArk)",
    "ONCMY": "PI request named it OMYKI; UniProt/OMA mnemonic ONCMY used (already in hobench STRICT_SPECIES)",
    "ZYGRC": "FTP canonical file not found because REST omits superkingdom; fetched via REST stream (4,987 = geneCount)",
}

# ---------------------------------------------------------------------------------------------------------------
# Species trees. Leaves = panel codes; internal labels = level ids. Branch lengths deliberately omitted
# (topology-only reference; tools that need lengths estimate them).
#
# Panel Y: pruned from Shen et al. 2018 Cell (doi 10.1016/j.cell.2018.10.023), figshare 10.6084/m9.figshare.5854692
# file data_in_Figure2/343_2408OGs_unpartitioned_ML_species-names_updated.newick (CC BY 4.0), rooted on
# S. pombe, pruned with ete3 to the 16 panel species (Lachancea_thermotolerans, not Phaffomyces_thermotolerans).
# The 332-taxon time-calibrated tree in the same archive gives the identical topology for the 15 species it has.
# Bootstrap support after pruning: 100 at every node except Saccharomyces+Nakaseomyces vs
# Naumovozyma+Kazachstania (52) -> PostWGD_crown is flagged contested.
TREES = {
    "Y": "(SCHPO,(YARLI,(CANAL,((((ERECY,ASHGO)Eremothecium,KLULA)Kluyveromyces_Eremothecium,LACTH)KLE_clade,"
         "((TORDC,ZYGRC)ZT_clade,((((NAUDC,NAUCC)Naumovozyma,KAZAF)Naumovozyma_Kazachstania,"
         "(YEAST,CANGA)Saccharomyces_Nakaseomyces)PostWGD_crown,(VANPO,TETPH)Vanderwaltozyma_Tetrapisispora)PostWGD)"
         "ZT_PostWGD)Saccharomycetaceae)CUG_Saccharomycetaceae)Saccharomycotina)Ascomycota;",
    # Panel K: two-domain rooting between Bacteria and Archaea (Williams 2020). Bacteria rooted between
    # Terrabacteria and Gracilicutes (Coleman 2021). Phylum-level relationships inside the two bacterial
    # superphyla and inside archaea follow GTDB-based phylogenomics (Parks 2018/2022) and Hug 2016 and are
    # LOW CONFIDENCE (HGT-rich, short internal branches): every non-trivial internal bacterial node is flagged
    # contested in levels.tsv. Placement of Thermotogota (with Deinococcota as "DST" in Terrabacteria) and of
    # Aquificota (with Campylobacterota in Gracilicutes) is the most uncertain choice.
    "K": "((((DEIRA,THEMA)DST_clade,(SYNY3,(BACSU,MYCTU)Bacillota_Actinomycetota)Terrabacteria_core)Terrabacteria,"
         "((LEPIN,(BACTN,CHLTR)FCB_PVC)Spirochaetota_FCB_PVC,(((ECOLI,PSEAE)Gammaproteobacteria,BRADU)Pseudomonadota,"
         "(GEOSL,(HELPY,AQUAE)Campylobacterota_Aquificota)Desulfobacterota_Campylobacterota)Proteobacteria_sl)"
         "Gracilicutes)Bacteria,"
         "((THEKO,(METJA,(METAC,HALSA)Methanomicrobia_Halobacteria)Euryarchaeota_core)Euryarchaeota,"
         "(NITMS,(KORCO,SACS2)Crenarchaeota_Korarchaeota)TACK)Archaea)LUCA;",
}

# One plausible alternative topology per panel at a contested node (species-tree-error robustness experiment).
ALT_TREES = {
    # Alternative at ZT_PostWGD: post-WGD lineage sister to KLE instead of ZT, i.e. the signal of the KLE-like
    # parent of the hybrid WGD ancestor (Marcet-Houben & Gabaldon 2015; Wolfe 2015).
    "Y": "(SCHPO,(YARLI,(CANAL,((TORDC,ZYGRC)ZT_clade,((((ERECY,ASHGO)Eremothecium,KLULA)Kluyveromyces_Eremothecium,"
         "LACTH)KLE_clade,((((NAUDC,NAUCC)Naumovozyma,KAZAF)Naumovozyma_Kazachstania,"
         "(YEAST,CANGA)Saccharomyces_Nakaseomyces)PostWGD_crown,(VANPO,TETPH)Vanderwaltozyma_Tetrapisispora)PostWGD)"
         "KLE_PostWGD)Saccharomycetaceae)CUG_Saccharomycetaceae)Saccharomycotina)Ascomycota;",
    # Alternative: Aquificota sister to Thermotogota (classic hyperthermophile grouping of rRNA trees) inside
    # Terrabacteria instead of with Campylobacterota in Gracilicutes.
    "K": "((((DEIRA,(THEMA,AQUAE)Thermotogota_Aquificota)DST_clade,(SYNY3,(BACSU,MYCTU)Bacillota_Actinomycetota)"
         "Terrabacteria_core)Terrabacteria,"
         "((LEPIN,(BACTN,CHLTR)FCB_PVC)Spirochaetota_FCB_PVC,(((ECOLI,PSEAE)Gammaproteobacteria,BRADU)Pseudomonadota,"
         "(GEOSL,HELPY)Desulfobacterota_Campylobacterota)Proteobacteria_sl)Gracilicutes)Bacteria,"
         "((THEKO,(METJA,(METAC,HALSA)Methanomicrobia_Halobacteria)Euryarchaeota_core)Euryarchaeota,"
         "(NITMS,(KORCO,SACS2)Crenarchaeota_Korarchaeota)TACK)Archaea)LUCA;",
}

# Level annotations: level_id -> (literature_name, contested y/n, alternative topology, citation keys, note)
_ = None
LEVELS = {
    "Y": {
        "Ascomycota": ("Ascomycota (Taphrinomycotina + Saccharomycotina)", "n", "", ["Shen2018"], "root; SCHPO outgroup"),
        "Saccharomycotina": ("Saccharomycotina", "n", "", ["Shen2018", "Shen2016"], ""),
        "CUG_Saccharomycetaceae": ("CUG-Ser1 clade + Saccharomycetaceae (Shen 2018 major clades)", "n", "",
                                   ["Shen2018"], "Yarrowia (Dipodascaceae/Trigonopsidaceae) branches outside"),
        "Saccharomycetaceae": ("Saccharomycetaceae", "n", "", ["Shen2018"], ""),
        "KLE_clade": ("KLE clade (Kluyveromyces, Lachancea, Eremothecium)", "y",
                      "Lachancea outside KLE / KLE paraphyletic in some analyses; KLE-like parent of the WGD hybrid",
                      ["Shen2018", "MarcetHouben2015"], "BS 100 in Shen 2018; treat as soft for WGD-gene scoring"),
        "Kluyveromyces_Eremothecium": ("Kluyveromyces + Eremothecium", "y",
                                       "(Lachancea,Kluyveromyces) sister to Eremothecium",
                                       ["Shen2018"], "Lachancea is the unstable member (survey C 4.2)"),
        "Eremothecium": ("Eremothecium", "n", "", ["Shen2018"], ""),
        "ZT_PostWGD": ("ZT clade + post-WGD clade", "y",
                       "post-WGD sister to KLE (species_tree_alt.nwk); ~30% of post-WGD genes cluster with KLE "
                       "because the WGD followed an interspecies hybridization",
                       ["Shen2018", "MarcetHouben2015", "Wolfe2015"], "species tree = ZT parent signal"),
        "ZT_clade": ("ZT clade (Zygosaccharomyces, Torulaspora)", "n", "", ["Shen2018"], ""),
        "PostWGD": ("post-WGD clade", "n", "", ["Shen2018", "WolfeShields1997"], "WGD on the branch into this level"),
        "PostWGD_crown": ("Saccharomyces/Nakaseomyces + Naumovozyma/Kazachstania", "y",
                          "(Vanderwaltozyma,Tetrapisispora) sister to Naumovozyma/Kazachstania",
                          ["Shen2018"], "bootstrap 52 in the pruned Shen 2018 ML tree"),
        "Naumovozyma_Kazachstania": ("Naumovozyma + Kazachstania", "n", "", ["Shen2018"], ""),
        "Naumovozyma": ("Naumovozyma", "n", "", ["Shen2018"], ""),
        "Saccharomyces_Nakaseomyces": ("Saccharomyces + Nakaseomyces", "n", "", ["Shen2018"], ""),
        "Vanderwaltozyma_Tetrapisispora": ("Vanderwaltozyma + Tetrapisispora", "n", "", ["Shen2018"], ""),
        "KLE_PostWGD": ("KLE + post-WGD (alternative tree only)", "y", "alt only", ["MarcetHouben2015"], ""),
    },
    "K": {
        "LUCA": ("Last universal common ancestor (Bacteria/Archaea split)", "y",
                 "root position inside Bacteria (not assessed)", ["Williams2020", "Coleman2021"],
                 "two-domain rooting between Bacteria and Archaea"),
        "Bacteria": ("Bacteria", "y", "root of Bacteria elsewhere (e.g. CPR- or Thermotogota-first)",
                     ["Coleman2021"], "rooted between Terrabacteria and Gracilicutes"),
        "Terrabacteria": ("Terrabacteria", "y", "membership of Thermotogota disputed", ["Coleman2021", "Hug2016"], ""),
        "DST_clade": ("Deinococcota + Thermotogota (DST: Deinococcota, Synergistota, Thermotogota)", "y",
                      "Thermotogota near the bacterial root or with Aquificota (species_tree_alt.nwk)",
                      ["Coleman2021", "Hug2016"], "low confidence"),
        "Terrabacteria_core": ("Cyanobacteriota + Bacillota + Actinomycetota", "y",
                               "Cyanobacteriota sister to Actinomycetota or to Chloroflexota-group",
                               ["Coleman2021", "Parks2022GTDB"], "low confidence"),
        "Bacillota_Actinomycetota": ("Bacillota (Firmicutes) + Actinomycetota", "y",
                                     "Bacillota sister to the rest of Terrabacteria",
                                     ["Coleman2021", "Parks2022GTDB"], "low confidence"),
        "Gracilicutes": ("Gracilicutes", "n", "", ["Coleman2021"], ""),
        "Spirochaetota_FCB_PVC": ("Spirochaetota + FCB + PVC", "y", "Spirochaetota sister to Proteobacteria s.l.",
                                  ["Coleman2021", "Hug2016"], "low confidence"),
        "FCB_PVC": ("FCB (Bacteroidota) + PVC (Chlamydiota) superphyla", "y", "FCB and PVC not sisters",
                    ["Hug2016", "Coleman2021"], "low confidence"),
        "Proteobacteria_sl": ("Proteobacteria sensu lato (incl. Desulfobacterota, Campylobacterota)", "y",
                              "Proteobacteria s.l. paraphyletic (GTDB splits it)", ["Parks2018", "Parks2022GTDB"], ""),
        "Pseudomonadota": ("Pseudomonadota (Alpha + Gammaproteobacteria)", "n", "", ["Parks2018", "Parks2022GTDB"], ""),
        "Gammaproteobacteria": ("Gammaproteobacteria", "n", "", ["Parks2022GTDB"], ""),
        "Desulfobacterota_Campylobacterota": ("Desulfobacterota + Campylobacterota (+Aquificota)", "y",
                                              "Campylobacterota sister to Pseudomonadota",
                                              ["Parks2018", "Coleman2021"], "low confidence"),
        "Campylobacterota_Aquificota": ("Campylobacterota + Aquificota", "y",
                                        "Aquificota with Thermotogota (species_tree_alt.nwk) or near the root",
                                        ["Coleman2021", "Hug2016"], "low confidence; extensive HGT in Aquificota"),
        "Thermotogota_Aquificota": ("Thermotogota + Aquificota (alternative tree only)", "y", "alt only",
                                    ["Hug2016"], ""),
        "Archaea": ("Archaea", "n", "", ["Williams2020", "Eme2023"], ""),
        "Euryarchaeota": ("Euryarchaeota (Thermococci, Methanococci, Methanomicrobia, Halobacteria)", "y",
                          "Euryarchaeota paraphyletic w.r.t. archaeal root", ["Williams2020", "Parks2022GTDB"], ""),
        "Euryarchaeota_core": ("Methanococci + Methanomicrobia/Halobacteria", "y",
                               "Thermococci sister to Methanococci", ["Parks2022GTDB"], "low confidence"),
        "Methanomicrobia_Halobacteria": ("Methanomicrobia + Halobacteria (Stenosarchaea)", "n", "",
                                         ["Parks2022GTDB"], ""),
        "TACK": ("TACK superphylum (Thaumarchaeota, Crenarchaeota, Korarchaeota)", "n", "",
                 ["Williams2020", "Eme2023"], ""),
        "Crenarchaeota_Korarchaeota": ("Crenarchaeota + Korarchaeota", "y", "Korarchaeota sister to Thaum+Cren",
                                       ["Williams2020"], "low confidence"),
    },
}

# ---------------------------------------------------------------------------------------------------------------
# Panel V. Chordate backbone: tunicates sister to vertebrates (Delsuc 2006); jawed-vertebrate timetree
# (Irisarri 2017); ray-finned fishes Betancur-R 2017 / Hughes 2018 (gar = Holostei sister to Teleostei;
# Otomorpha (zebrafish) sister to Euteleosteomorpha; Esociformes sister to Salmoniformes; Coregoninae sister to
# Salmoninae). Only topology; no polytomy needed resolving except the trivial ones of the NCBI tree.
TREES["V"] = ("(BRAFL,(CIOIN,(PETMA,(CALMI,((LEPOC,(DANRE,((ESOLU,(CORBA,(SALSA,ONCMY)Salmoninae)Salmonidae)"
              "Protacanthopterygii,ORYLA)Euteleosteomorpha)Teleostei)Neopterygii,(LATCH,((XENTR,XENLA)Xenopus,"
              "(CHICK,(HUMAN,MOUSE)Euarchontoglires)Amniota)Tetrapoda)Sarcopterygii)Euteleostomi)Gnathostomata)"
              "Vertebrata)Olfactores)Chordata;")
# Alternative: "Euchordata" (cephalochordates sister to vertebrates), the pre-2006 consensus refuted by Delsuc 2006.
ALT_TREES["V"] = TREES["V"].replace("(BRAFL,(CIOIN,(PETMA", "(CIOIN,(BRAFL,(PETMA").replace(
    ")Vertebrata)Olfactores)Chordata;", ")Vertebrata)Euchordata)Chordata;")
LEVELS["V"] = {
    "Chordata": ("Chordata", "n", "", ["Putnam2008", "Delsuc2006"], "root; BRAFL = pre-1R outgroup"),
    "Olfactores": ("Olfactores (tunicates + vertebrates)", "y",
                   "Euchordata: cephalochordates sister to vertebrates (species_tree_alt.nwk)", ["Delsuc2006"],
                   "Delsuc 2006 phylogenomics; tunicates fast-evolving"),
    "Vertebrata": ("Vertebrata", "n", "", ["Simakov2020", "Nakatani2021"],
                   "1R on branch into Vertebrata (Simakov 2020); timing relative to cyclostome split contested "
                   "(Nakatani 2021: cyclostomes after 1R, plus cyclostome-specific hexaploidization) -> interval event"),
    "Gnathostomata": ("Gnathostomata", "n", "", ["Simakov2020", "Nakatani2021", "Irisarri2017"],
                      "2R (allotetraploid) on branch into Gnathostomata per Simakov 2020 / Nakatani 2021"),
    "Euteleostomi": ("Osteichthyes (bony vertebrates)", "n", "", ["Irisarri2017"], ""),
    "Neopterygii": ("Neopterygii (Holostei + Teleostei)", "n", "", ["Betancur2017", "Hughes2018", "Braasch2016"],
                    "gar = unduplicated pre-3R outgroup"),
    "Teleostei": ("Teleostei (sampled crown = Clupeocephala; no Elopomorpha/Osteoglossomorpha)", "n", "", ["Hughes2018", "Betancur2017"], "teleost 3R (TGD) on branch into Teleostei"),
    "Euteleosteomorpha": ("Euteleosteomorpha", "n", "", ["Betancur2017", "Hughes2018"], ""),
    "Protacanthopterygii": ("Esociformes + Salmoniformes", "n", "", ["Betancur2017", "Lien2016"],
                            "pike = unduplicated outgroup of Ss4R"),
    "Salmonidae": ("Salmonidae", "n", "", ["Lien2016", "Robertson2017"],
                   "Ss4R on branch into Salmonidae; LORe regions rediploidized independently in Coregoninae and "
                   "Salmoninae (Robertson 2017) -> per-lineage D events"),
    "Salmoninae": ("Salmoninae (Salmo + Oncorhynchus)", "n", "", ["Robertson2017"], ""),
    "Sarcopterygii": ("Sarcopterygii", "n", "", ["Meyer2021lungfish", "Irisarri2017"],
                      "coelacanth outside lungfish+tetrapods; no lungfish in panel"),
    "Tetrapoda": ("Tetrapoda", "n", "", ["Irisarri2017"], ""),
    "Xenopus": ("Xenopus", "n", "", ["Session2016"], "X. laevis allotetraploidy (L/S) = D on the XENLA leaf"),
    "Amniota": ("Amniota", "n", "", ["Irisarri2017"], ""),
    "Euarchontoglires": ("Euarchontoglires", "n", "", ["Irisarri2017"], ""),
    "Euchordata": ("Cephalochordata + Vertebrata (alternative tree only)", "y", "alt only", ["Delsuc2006"], ""),
}

# Panel E. Two-domain tree: eukaryotes nested in Archaea as sister to Asgard (Williams 2020; Eme 2023), TACK
# sister to Asgard+Eukaryota; bacterial part as panel K. Eukaryotic root between Opimoda (here = Amorphea) and
# Diphoda (Derelle 2015; Excavata not monophyletic: Discoba + Metamonada on the Diphoda side); supergroups after
# Burki 2020; Obazoa (apusomonads sister to opisthokonts, Brown 2013); Microsporidia earliest-diverging sequenced
# fungi (Capella-Gutierrez 2012; James 2006). Animal root: ctenophores-first (Whelan 2015; Schultz 2023 synteny),
# Placozoa sister to Cnidaria+Bilateria (Laumer 2019 = ParaHoxozoa); alternative sponges-first.
TREES["E"] = ("(((BACSU,SYNY3)Terrabacteria,(ECOLI,BRADU)Pseudomonadota)Bacteria,((THEKO,HALSA)Euryarchaeota,"
              "((NITMS,SACS2)TACK,(LOKOS,(((DICDI,ENTH1)Amoebozoa,(THETR,((ENCCU,(BATDJ,(YEAST,NEUCR)Ascomycota)"
              "Chytridiomycota_Dikarya)Fungi,(CAPO3,(MONBE,(MNELE,(AMPQE,(TRIAD,(NEMVE,(HUMAN,DROME)Bilateria)"
              "Planulozoa)ParaHoxozoa)Myriazoa)Metazoa)Choanozoa)Holozoa)Opisthokonta)Obazoa)Amorphea,"
              "(((CYAME,(ARATH,CHLRE)Viridiplantae)Archaeplastida,((THAPS,PHYRM)Stramenopiles,(PLAF7,ICHMU)Alveolata)"
              "Stramenopiles_Alveolata)Diaphoretickes,((TRYB2,NAEGR)Discoba,GIAIC)Discoba_Metamonada)Diphoda)"
              "Eukaryota)Asgard_Eukaryota)TACK_Asgard_Eukaryota)Archaea_Eukaryota)LUCA;")
ALT_TREES["E"] = TREES["E"].replace("(MNELE,(AMPQE,(TRIAD,(NEMVE,(HUMAN,DROME)Bilateria)Planulozoa)ParaHoxozoa)Myriazoa)",
                                    "(AMPQE,(MNELE,(TRIAD,(NEMVE,(HUMAN,DROME)Bilateria)Planulozoa)ParaHoxozoa)"
                                    "Ctenophora_ParaHoxozoa)")
LEVELS["E"] = {
    "LUCA": ("LUCA (Bacteria / Archaea+Eukaryota split)", "y", "three-domain tree (eukaryotes sister to Archaea)",
             ["Williams2020", "Coleman2021"], ""),
    "Bacteria": ("Bacteria", "n", "", ["Coleman2021"], "rooted between Terrabacteria and Gracilicutes"),
    "Terrabacteria": ("Terrabacteria", "n", "", ["Coleman2021"], ""),
    "Pseudomonadota": ("Pseudomonadota (Alpha + Gamma)", "n", "", ["Parks2022GTDB"], "mitochondrial ancestry (BRADU)"),
    "Archaea_Eukaryota": ("Archaea + Eukaryota (two-domain)", "y", "Eukaryota sister to monophyletic Archaea",
                          ["Williams2020", "Eme2023"], ""),
    "Euryarchaeota": ("Euryarchaeota", "n", "", ["Williams2020"], ""),
    "TACK_Asgard_Eukaryota": ("TACK + Asgard + Eukaryota", "y", "TACK sister to Euryarchaeota+Asgard",
                              ["Williams2020", "Eme2023"], ""),
    "TACK": ("TACK", "n", "", ["Williams2020"], ""),
    "Asgard_Eukaryota": ("Asgard archaea + Eukaryota", "y", "eukaryotes nested within Asgard (Heimdallarchaeia; "
                         "Eme 2023) - same panel topology since only Lokiarchaeia is sampled",
                         ["Eme2023", "Zaremba2017", "Williams2020"], ""),
    "Eukaryota": ("Eukaryota (LECA)", "y", "root between Discoba and the rest, or within Excavata",
                  ["Derelle2015", "Burki2020"], "root between Opimoda and Diphoda"),
    "Amorphea": ("Amorphea (= Opimoda in this panel: no malawimonads/CRuMs)", "n", "", ["Derelle2015", "Burki2020"], ""),
    "Amoebozoa": ("Amoebozoa", "n", "", ["Burki2020"], ""),
    "Obazoa": ("Obazoa (Apusomonadida + Opisthokonta)", "n", "", ["Brown2013", "Burki2020"], ""),
    "Opisthokonta": ("Opisthokonta", "n", "", ["Burki2020"], ""),
    "Fungi": ("Fungi (incl. Microsporidia)", "y", "Microsporidia long-branch placement (attraction to other fast taxa)",
              ["CapellaGutierrez2012", "James2006"], ""),
    "Chytridiomycota_Dikarya": ("Chytridiomycota + Dikarya", "n", "", ["James2006"], ""),
    "Ascomycota": ("Ascomycota", "n", "", ["James2006"], ""),
    "Holozoa": ("Holozoa", "n", "", ["Burki2020"], ""),
    "Choanozoa": ("Choanoflagellata + Metazoa", "n", "", ["Burki2020"], ""),
    "Metazoa": ("Metazoa", "y", "sponges-first (species_tree_alt.nwk; Simion 2017, Redmond 2021) vs "
                "ctenophores-first (main; Whelan 2015, Schultz 2023); may reflect systematic error (Kapli 2020)",
                ["Schultz2023", "Whelan2015", "Simion2017", "Redmond2021", "Kapli2020"], "score root-agnostically"),
    "Myriazoa": ("Porifera + ParaHoxozoa (ctenophores-first)", "y", "alt: Ctenophora_ParaHoxozoa",
                 ["Schultz2023", "Whelan2015"], ""),
    "ParaHoxozoa": ("ParaHoxozoa (Placozoa + Cnidaria + Bilateria)", "y", "Placozoa sister to Cnidaria",
                    ["Laumer2019"], ""),
    "Planulozoa": ("Cnidaria + Bilateria", "y", "Placozoa+Cnidaria clade", ["Laumer2019"], ""),
    "Bilateria": ("Bilateria", "n", "", ["Laumer2019"], ""),
    "Diphoda": ("Diphoda (Diaphoretickes + Discoba + Metamonada)", "y", "Excavata monophyletic or root within it",
                ["Derelle2015", "Burki2020"], ""),
    "Diaphoretickes": ("Diaphoretickes", "n", "", ["Burki2020"], ""),
    "Archaeplastida": ("Archaeplastida", "y", "Archaeplastida non-monophyletic (cryptists inside)", ["Burki2020"], ""),
    "Viridiplantae": ("Viridiplantae", "n", "", ["OneKP2019"], ""),
    "Stramenopiles_Alveolata": ("SAR without Rhizaria (TSAR)", "n", "", ["Burki2020"], "no Rhizaria in panel"),
    "Stramenopiles": ("Stramenopiles", "n", "", ["Burki2020"], ""),
    "Alveolata": ("Alveolata", "n", "", ["Burki2020"], ""),
    "Discoba_Metamonada": ("Discoba + Metamonada (excavates)", "y", "Discoba and Metamonada not sisters",
                           ["Derelle2015", "Burki2020"], ""),
    "Discoba": ("Discoba", "n", "", ["Burki2020"], ""),
    "Ctenophora_ParaHoxozoa": ("Ctenophora + ParaHoxozoa (sponges-first; alternative tree only)", "y", "alt only",
                               ["Simion2017", "Redmond2021"], ""),
}

# Panel P. OneKP 2019 backbone; bryophyte monophyly (Puttick 2018; Harris 2020; OneKP 2019) - with only a
# liverwort and a moss this is the Setaphyta node; Amborella sister to all other angiosperms (Amborella 2013;
# Zuntini 2024); Vitales sister to the remaining rosids; Brassicaceae (malvids) and Fabaceae (fabids) as eurosids.
TREES["P"] = ("((OSTLU,CHLRE)Chlorophyta,(KLENI,((MARPO,PHYPA)Setaphyta,(SELMO,(AMBTC,((ORYSJ,MAIZE)Poaceae,(VITVI,"
              "((ARATH,BRARP)Brassicaceae,(SOYBN,PHAVU)Phaseoleae)Eurosids)Rosids)Mesangiospermae)Magnoliopsida)"
              "Tracheophyta)Embryophyta)Streptophyta)Viridiplantae;")
# Alternative: liverworts sister to all other land plants (bryophyte paraphyly; Qiu 2006).
ALT_TREES["P"] = TREES["P"].replace("((MARPO,PHYPA)Setaphyta,(SELMO,", "(MARPO,(PHYPA,(SELMO,").replace(
    ")Tracheophyta)Embryophyta)", ")Tracheophyta)Mosses_Tracheophyta)Embryophyta)")
LEVELS["P"] = {
    "Viridiplantae": ("Viridiplantae", "n", "", ["OneKP2019"], "root"),
    "Chlorophyta": ("Chlorophyta", "n", "", ["OneKP2019"], ""),
    "Streptophyta": ("Streptophyta", "n", "", ["OneKP2019", "Puttick2018"],
                     "KLENI = Klebsormidiophyceae; Zygnematophyceae (land-plant sister) not sampled"),
    "Embryophyta": ("Embryophyta (land plants)", "n", "", ["Puttick2018", "OneKP2019"], ""),
    "Setaphyta": ("Setaphyta (liverworts + mosses; bryophytes)", "y",
                  "bryophyte paraphyly, liverworts-first (species_tree_alt.nwk; Qiu 2006)",
                  ["Puttick2018", "Harris2020", "OneKP2019", "Qiu2006"], ""),
    "Tracheophyta": ("Tracheophyta", "n", "", ["OneKP2019"], ""),
    "Magnoliopsida": ("Angiospermae", "y", "Amborella+Nymphaeales sister to other angiosperms (same panel topology; "
                      "no Nymphaeales sampled)", ["Amborella2013", "Zuntini2024", "OneKP2019"],
                      "epsilon (angiosperm) WGD on branch into Magnoliopsida; zeta (seed plant) on the same panel "
                      "branch since no gymnosperm is sampled (Jiao 2011)"),
    "Mesangiospermae": ("Mesangiospermae (monocots + eudicots)", "n", "", ["Zuntini2024"], ""),
    "Poaceae": ("Poaceae", "n", "", ["Zuntini2024"], "rho WGD before the rice/maize split; maize tetraploidy on MAIZE leaf"),
    "Rosids": ("Rosids incl. Vitales (core eudicot representatives)", "n", "", ["Zuntini2024"],
               "gamma hexaploidy before the VITVI split"),
    "Eurosids": ("Eurosids (malvids + fabids)", "n", "", ["Zuntini2024"], ""),
    "Brassicaceae": ("Brassicaceae", "n", "", ["Wang2011Brassica"],
                     "alpha/beta WGDs before ARATH/BRARP split; Brassica triplication on BRARP leaf"),
    "Phaseoleae": ("Phaseoleae (Glycine + Phaseolus)", "n", "", ["Schmutz2010"],
                   "legume WGD (~59 Ma) before the split; Glycine WGD (~13 Ma) on SOYBN leaf"),
    "Mosses_Tracheophyta": ("mosses + vascular plants (alternative tree only)", "y", "alt only", ["Qiu2006"], ""),
}


# ---------------------------------------------------------------------------------------------------------------
def sha256_file(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def proteome_dir(upid: str) -> Path:
    c = sorted(PROTEOMES.glob(f"{upid}_*"))
    if len(c) != 1:
        raise SystemExit(f"ERROR: expected exactly one proteome dir for {upid}, found {c}")
    return c[0]


def uniprot_mnemonic(hmap: Path) -> str:
    cnt = Counter()
    with open(hmap) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if "_" in r["entry_name"]:
                cnt[r["entry_name"].rsplit("_", 1)[1]] += 1
    return cnt.most_common(1)[0][0] if cnt else ""


def panel_species(panel: str) -> list[dict]:
    rows = [r for r in csv.DictReader(open(VERIFIED), delimiter="\t") if r["panel"] == panel and r["tier"] == "core"]
    out = [dict(code=r["mnemonic"], upid=r["upid"], why=r["why"],
                source_set="QfO" if r["QfO2026_02"] == "yes" else "RP") for r in rows]
    for code, (upid, why) in EXTRA.get(panel, {}).items():
        out.append(dict(code=code, upid=upid, why=why, source_set="RP"))
    return out


def build_species(panel: str, check_sha: bool) -> list[dict]:
    pdir = PANELS_DIR / panel
    (pdir / "proteomes").mkdir(parents=True, exist_ok=True)
    rows = []
    for s in panel_species(panel):
        d = proteome_dir(s["upid"])
        m = json.loads((d / "manifest.json").read_text())
        if m["uniprot_release"] != PINNED_RELEASE:
            raise SystemExit(f"ERROR: {s['code']} {s['upid']} release {m['uniprot_release']} != {PINNED_RELEASE}")
        fa = d / "proteome.fasta"
        sha = m["outputs"]["proteome_fasta"]["sha256"]
        if check_sha and sha256_file(fa) != sha:
            raise SystemExit(f"ERROR: sha256 of {fa} differs from manifest")
        n = sum(1 for ln in open(fa) if ln.startswith(">"))
        if n != m["n_proteins"] or n == 0:
            raise SystemExit(f"ERROR: {fa}: {n} sequences vs manifest {m['n_proteins']}")
        link = pdir / "proteomes" / f"{s['code']}.fa"
        target = os.path.relpath(fa, link.parent)
        if link.is_symlink() or link.exists():
            link.unlink()
        link.symlink_to(target)
        why = s["why"] + (f"; {NOTES[s['code']]}" if s["code"] in NOTES else "")
        rows.append(dict(code=s["code"], species=m["organism"] or "", taxid=m["taxid"], upid=s["upid"],
                         mnemonic=uniprot_mnemonic(d / "header_map.tsv"), uniprot_release=m["uniprot_release"],
                         proteome_type=m["proteome_type"], fetch_source=m["source"],
                         fasta_path=os.path.relpath(fa, ROOT), n_proteins=n, sha256=sha,
                         source_set=s["source_set"], strict_mapping=int(s["code"] in STRICT), why=why))
    # remove stale links
    codes = {r["code"] for r in rows}
    for f in (pdir / "proteomes").glob("*.fa"):
        if f.stem not in codes:
            f.unlink()
    if len(codes) != len(rows):
        raise SystemExit(f"ERROR: duplicate codes in panel {panel}")
    with open(pdir / "species.tsv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    return rows


def validate_tree(nwk: str, codes: set[str], what: str) -> Tree:
    t = Tree(nwk, format=8)
    errs = []
    if len(t.children) != 2:
        errs.append(f"root has {len(t.children)} children (must be rooted binary)")
    for n in t.traverse():
        if not n.is_leaf() and len(n.children) != 2:
            errs.append(f"node {n.name} has {len(n.children)} children")
        if not n.is_leaf() and not LEVEL_RE.match(n.name or ""):
            errs.append(f"bad level id {n.name!r}")
    names = [n.name for n in t.traverse()]
    dup = [k for k, v in Counter(names).items() if v > 1]
    if dup:
        errs.append(f"duplicate node names {dup}")
    leaves = set(t.get_leaf_names())
    if leaves != codes:
        errs.append(f"leaves != species codes: missing {sorted(codes - leaves)}, extra {sorted(leaves - codes)}")
    if errs:
        raise SystemExit(f"ERROR in {what}:\n  " + "\n  ".join(errs))
    return t


def load_citations() -> dict:
    if not CITATIONS.exists():
        return {}
    return {r["key"]: r for r in csv.DictReader(open(CITATIONS), delimiter="\t")}


def verify_dois(keys) -> None:
    have = load_citations()
    for k in sorted(keys):
        doi = CIT[k]
        if k in have and have[k]["doi"] == doi and have[k]["status"] == "VERIFIED":
            continue
        url = ("https://www.ebi.ac.uk/europepmc/webservices/rest/search?" +
               urllib.parse.urlencode({"query": f'DOI:"{doi}"', "format": "json", "resultType": "lite"}))
        res = json.load(urllib.request.urlopen(url, timeout=60))["resultList"]["result"]
        r = res[0] if res else {}
        ok = bool(r) and (r.get("doi") or "").lower() == doi.lower()
        have[k] = dict(key=k, doi=doi, status="VERIFIED" if ok else "NOT_FOUND", title=r.get("title", ""),
                       journal=r.get("journalTitle", ""), year=r.get("pubYear", ""),
                       checked=time.strftime("%Y-%m-%d"), source="Europe PMC REST search DOI:")
        time.sleep(0.2)
    with open(CITATIONS, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["key", "doi", "status", "title", "journal", "year", "checked", "source"],
                           delimiter="\t", lineterminator="\n")
        w.writeheader()
        for k in sorted(have):
            w.writerow(have[k])


def ncbi_annot(t: Tree, taxid_of: dict, ncbi) -> dict:
    """level -> (ncbi_taxid, ncbi_name, exact y/n)."""
    lineages = {c: ncbi.get_lineage(int(x)) for c, x in taxid_of.items()}
    out = {}
    for n in t.traverse():
        if n.is_leaf():
            continue
        sp = n.get_leaf_names()
        common = set(lineages[sp[0]])
        for s in sp[1:]:
            common &= set(lineages[s])
        lca = [x for x in lineages[sp[0]] if x in common][-1]
        inside = {c for c, lin in lineages.items() if lca in lin}
        out[n.name] = (lca, ncbi.get_taxid_translator([lca])[lca], "y" if inside == set(sp) else "n")
    return out


def harmonise_alt(t: Tree, alt: Tree) -> None:
    """A level id shared by main and alt trees must denote the same species set; otherwise the alt node is
    renamed <id>_alt (so scoring never confuses two different clades under one name)."""
    main = {n.name: frozenset(n.get_leaf_names()) for n in t.traverse() if not n.is_leaf()}
    for n in alt.traverse():
        if not n.is_leaf() and n.name in main and main[n.name] != frozenset(n.get_leaf_names()):
            n.name = n.name + "_alt"


def write_levels(panel: str, t: Tree, alt: Tree, taxid_of: dict, ncbi, cits: dict) -> None:
    meta = LEVELS[panel]
    ann = ncbi_annot(t, taxid_of, ncbi)
    alt_ann = ncbi_annot(alt, taxid_of, ncbi)
    rows = []
    for tree, tag, a in ((t, "main", ann), (alt, "alt_only", alt_ann)):
        for n in tree.traverse("preorder"):
            if n.is_leaf() or (tag == "alt_only" and n.name in ann):
                continue
            base = n.name[:-4] if n.name.endswith("_alt") and n.name not in meta else n.name
            if base not in meta:
                raise SystemExit(f"ERROR: level {n.name} of panel {panel} has no LEVELS annotation")
            lit, contested, alt_txt, keys, note = meta[base]
            if base != n.name:
                contested, alt_txt, note = "y", "alt only", f"alternative-tree clade differing from main-tree {base}"
            for k in keys:
                if k not in CIT:
                    raise SystemExit(f"ERROR: citation key {k} not in CIT")
                if cits.get(k, {}).get("status") != "VERIFIED":
                    raise SystemExit(f"ERROR: citation {k} not verified (run --verify-dois)")
            rows.append(dict(level_id=n.name, tree=tag, ncbi_taxid=a[n.name][0], ncbi_name=a[n.name][1],
                             ncbi_exact=a[n.name][2], literature_name=lit,
                             children=",".join(c.name for c in n.children), n_species=len(n.get_leaf_names()),
                             contested=contested, alternative=alt_txt,
                             citation="; ".join(f"{k} doi:{CIT[k]}" for k in keys), note=note))
    unused = set(meta) - {r["level_id"] for r in rows} - {r["level_id"][:-4] for r in rows}
    if unused:
        raise SystemExit(f"ERROR: LEVELS[{panel}] has annotations for absent levels {sorted(unused)}")
    with open(PANELS_DIR / panel / "levels.tsv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--panels", nargs="+", default=list(TREES))
    ap.add_argument("--verify-dois", action="store_true")
    ap.add_argument("--no-sha", action="store_true", help="skip re-hashing FASTAs (still uses manifest sha256)")
    a = ap.parse_args(argv)
    from ete3 import NCBITaxa
    if not TAXDB.exists():
        raise SystemExit(f"ERROR: {TAXDB} missing; build with NCBITaxa(dbfile=..., taxdump_file='taxdump.tar.gz')")
    ncbi = NCBITaxa(dbfile=str(TAXDB))
    PANELS_DIR.mkdir(exist_ok=True)
    if a.verify_dois:
        verify_dois({k for p in a.panels for v in LEVELS[p].values() for k in v[3]})
    cits = load_citations()
    sys.path.insert(0, str(ROOT))
    from hobench.speciestree import SpeciesTree
    for p in a.panels:
        rows = build_species(p, not a.no_sha)
        codes = {r["code"] for r in rows}
        fa_codes = {f.stem for f in (PANELS_DIR / p / "proteomes").glob("*.fa")}
        assert fa_codes == codes, (p, fa_codes ^ codes)
        t = validate_tree(TREES[p], codes, f"{p} species_tree")
        alt = validate_tree(ALT_TREES[p], codes, f"{p} species_tree_alt")
        harmonise_alt(t, alt)
        if t.robinson_foulds(alt, unrooted_trees=False)[0] == 0:
            raise SystemExit(f"ERROR: {p} alt tree identical to main tree")
        for tr, fn in ((t, "species_tree.nwk"), (alt, "species_tree_alt.nwk")):
            (PANELS_DIR / p / fn).write_text(tr.write(format=8, format_root_node=True) + "\n")
            SpeciesTree.from_newick(PANELS_DIR / p / fn)  # project parser must accept it
        write_levels(p, t, alt, {r["code"]: r["taxid"] for r in rows}, ncbi, cits)
        print(f"panel {p}: {len(rows)} species, {sum(r['n_proteins'] for r in rows):,} proteins, "
              f"{len(t) - 1} internal levels; RF(main, alt) = {t.robinson_foulds(alt, unrooted_trees=False)[0]}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
