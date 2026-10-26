# B2 — Candidate reference gene families: plants, protists / eukaryote-wide (LECA), prokaryotes (incl. archaeal–eukaryotic)

Survey date: 2026-10-26 (T+~0h25 → T+~1h30 relative to START.txt). Author: B2 literature agent.

## 0. Method and verification standard

* **Citations.** Every paper below was resolved through the Europe PMC REST API (`search?query=TITLE:"…"` / `EXT_ID:` with `resultType=core`), which returns PMID, PMCID, DOI, `isOpenAccess` and `inPMC`. All raw hits are logged in `literature/raw/B2/epmc_hits.jsonl` (script: `literature/raw/B2/epmc_search.py`). The **OA** column below uses Europe PMC's flags: `OA` = isOpenAccess Y; `PMC` = in PMC (free to read, but not an OA licence, so the full-text XML/supplements are *not* served by Europe PMC); `PW` = neither (paywalled, or free only at the publisher). Numbers quoted from abstracts are in `literature/raw/B2/abstracts_numbers.txt`.
* **Data checks.** For OA papers, `literature/raw/B2/supp_scan.py` parsed the Europe PMC full-text XML to list `<supplementary-material>` items and repository links (Dryad, Zenodo, figshare, GitHub, TreeBASE). For 17 candidates I **downloaded and opened** the supplements (Europe PMC `/{PMCID}/supplementaryFiles`, GitHub raw, NCBI FTP) to confirm that they hold accessions and subfamily labels. Checksums are in `literature/raw/B2/MANIFEST.tsv`. Large files are git-ignored (`literature/raw/B2/.gitignore`), so they can be fetched again from the URLs in the manifest.
* **Status tags used in the text.** **[VERIFIED-DATA]**: I opened the file and saw accessions/labels. **[LISTED]**: the supplement or repository exists and its caption/metadata describes the content, but I did not open it. **UNVERIFIED**: a statement from memory that I did not check here (mostly approximate family sizes). Every bibliographic record is verified unless marked UNVERIFIED.
* **Circularity policy.** Reference memberships should come from literature-curated, phylogeny-based classifications. Some sources classify sequences with HMMs built from curated seeds (Rab, SNARE, HistoneDB, ARF "HMM" rows), use orthology tools (pan-NLRome orthogroups, COG/arCOG) or rely on synteny (WGD pair lists). Those are flagged *tool-assisted*. Use them only as secondary labels, or restrict them to the rows the authors marked as conclusively assigned.

Practical lessons (for LESSONS.md):
1. The Europe PMC `/{PMCID}/supplementaryFiles` endpoint returns a zip of all supplements for OA articles. It works for BMC, PLOS, eLife, GBE/MBE and Frontiers. One call per paper is fast, apart from the 79 MB homeobox bundle.
2. Dryad's v2 API returns metadata and file lists, but **file downloads return 401/403 to scripts** (bearer token needed, or bot-blocked). Plan a human/browser download, or a Dryad API token, for the phytochrome and SMC datasets.
3. Two curated web databases from the Fasshauer/Klöpper lab, the Rab DB (`bioinformatics.mpibpc.mpg.de/rab/`) and the SNARE DB (`…/snare/`), **no longer resolve**. Use the paper supplements instead. CyMoBase (cymobase.org), diArk, BLDB (bldb.eu), HistoneDB 2.0 and Nelson's P450 site (drnelson.uthsc.edu) are up (HTTP 200, 2026-10-26).
4. `pandas`/`xlrd` are not installed in base python. I installed them to the scratchpad (`pip --target`). The `conda_envs/ortho` env has pandas 3.0.6 but no xlrd for old `.xls` files.

---

## 1. Master table (42 candidates)

Phenomena tags: **WGD** (plant polyploidy events), **TD** (tandem/lineage-specific expansion), **LOSS** (differential loss), **HGT**, **EGT** (endosymbiotic gene transfer), **LECA** (ancient eukaryotic paralogs), **LUCA** (pre-LUCA/LUCA paralogs), **DOM** (promiscuous/repeat domains), **FUS** (fusion/fission), **FAST** (long branches / fast evolution), **CONV** (convergence / non-homologous analogues), **SC** (near-single-copy control), **ORG** (organelle-encoded members: absent from nuclear proteomes).

| ID | Family | Scope | Phenomena | Best machine-readable source | OA | Prio |
|---|---|---|---|---|---|---|
| E1 | Kinesin superfamily | 45 eukaryotes (all supergroups) | LECA, LOSS, TD, DOM | Wickstead 2010 Add. file 3 (1,626 rows: IDs, GenBank acc., family/subfamily) **[VERIFIED-DATA]** | OA | **5** |
| E2 | ARF/Arl/Sar GTPases | 114 eukaryotes | LECA, LOSS, FAST, TD | Vargová 2021 Suppl. Table 1 (2,286 rows, 31 OG labels, accessions) **[VERIFIED-DATA]** | OA | **5** |
| E3 | Rab GTPases | 81 spp. in alignments (LECA ~20–23 Rabs) | LECA, TD, LOSS | Klöpper 2012 FASTA alignments with Rab-labelled names (1,669 seqs) **[VERIFIED-DATA]**; TreeBASE S12811 | OA | **5** |
| E4 | Tubulin superfamily (α β γ δ ε ζ + bacterial BtubA/B; FtsZ outgroup) | 504 spp. | LECA, TD (recent dups), LOSS, SC-like core | Findeisen 2014 alignment (3,286 seqs, subgroup in name) + species table **[VERIFIED-DATA]** | OA | **5** |
| E5 | Myosins | 300+ eukaryotes | LECA, DOM, FUS, LOSS | Sebé-Pedrós 2014 DataS1.fasta (725) **[VERIFIED-DATA]**; Kollmar 2017 figshare 4565155 (CC BY; alignments + trees) **[LISTED]**; CyMoBase | OA | **5** |
| E6 | Dynein heavy chains (+IC, LC) | 24 → ~hundreds of eukaryotes | LECA, LOSS (land plants), long genes (annotation errors) | Kollmar 2016 Newick trees (~3,000 DHC leaves) **[VERIFIED-DATA]**; Wickstead 2007 File S2 IDs **[LISTED]** | OA | 4 |
| E7 | MCM2-7 (+MCM8/9, archaeal MCM) | eukaryotes + archaea | LECA, LOSS (MCM8/9), SC per paralog | Liu 2009 Newick treefiles (labelled, e.g. `TbMcm9`) **[VERIFIED-DATA]** | OA | **5** |
| E8 | SMC1-6 (+Rad50/SbcC) | eukaryotes; prokaryotic SMC | LECA, LUCA-ish, SC | Yoshinaga 2021 Dryad (accession lists xlsx, trees; CC0) **[LISTED]** | OA | 4 |
| E9 | Adaptins/AP1-5, COPI, TSET | eukaryotes | LECA (heterotetramer paralogy), LOSS | Hirst 2014 eLife Fig3-data1 (1,413 rows; 874 with accession) **[VERIFIED-DATA]**; Hirst 2011 Table S1 **[LISTED]** | OA | 4 |
| E10 | SNAREs (Qa/Qb/Qc/R; 20 subclasses) | 145 spp. | LECA, TD | Kienle 2009 Nexus trees zip **[LISTED]**; Kloepper 2007 (PMC) | OA/PMC | 3 |
| E11 | Histones + variants | eukaryotes | FAST(slow), TD (identical clusters), CONV (variants) | HistoneDB 2.0 (curated seeds + HMM; tool-assisted) | OA | 3 |
| E12 | Aminoacyl-tRNA synthetases (cyto vs mito vs plastid) | all domains | LUCA, HGT, EGT, dual targeting | Chaliotis 2017 Suppl. zip **[LISTED]**; Brindefalk 2007 (PW) | OA | **5** |
| E13 | EF-1α / EFL | eukaryotes (protists, green algae) | HGT, LOSS (mutually exclusive) | Keeling & Inagaki 2004 (PMC); Gile 2009 (PW) | PMC | 4 |
| E14 | Photolyase/cryptochrome superfamily | bacteria + eukaryotes | HGT, LUCA-ish, TD, LOSS | Mei & Dvornyk 2015 S1 Table (GenBank nt + protein acc.) **[VERIFIED-DATA]** | OA | 4 |
| E15 | Proteasome regulators (19S, PA28, PA200, PI31…) | eukaryotes | LECA, LOSS | Fort 2015 Supplemental_Tables.xls **[LISTED]** | OA | 3 |
| E16 | Actin-related proteins Arp1–10 | eukaryotes | LECA, FAST | Muller 2005 (PW) | PW | 3 |
| E17 | ESCRT (incl. Asgard homologs) | eukaryotes + Asgard | LECA, archaeal–euk. | Leung 2008 (PW); Hatano 2022 (OA, MOESM6 xlsx **[LISTED]**) | mixed | 3 |
| P1 | MADS-box MIKC^c (AP1/FUL, AP3/PI, AG, SEP, FLC, SOC1…) | streptophytes, angiosperm focus | **WGD** (γ, ρ, α, Brassica ×3), TD, LOSS | Parenicová 2003 (107 At genes; PW); Genetics papers (PMC) with accession tables in figures/text | PW/PMC | **5** |
| P2 | Phytochromes (+ neochrome HGT) | green plants + streptophyte algae | WGD (seed-plant, angiosperm), LOSS, HGT (hornwort→fern neochrome) | Li 2015 Dryad `10.5061/dryad.5rs50` (alignments + trees, CC0) **[LISTED]**, download blocked to scripts | OA | **5** |
| P3 | Homeobox TFs (14 classes: KNOX, BEL, HD-ZIP I–IV, WOX…) | 10 genomes (red alga→angiosperm) | WGD, LOSS, DOM | Mukherjee 2009 msp201_5.html table (name, accession, class, species, sequence) **[VERIFIED-DATA]** | OA | **5** |
| P4 | Class III HD-ZIP | streptophytes | WGD, miR165/166 conservation | Floyd 2006 (PMC); Prigge & Clark 2006 (PW) | PMC | 3 |
| P5 | Nuclear auxin pathway: ARF A/B/C, Aux/IAA, TIR1/AFB | streptophytes (1KP transcriptomes) | WGD, LOSS, origin in charophytes | Mutte 2018 eLife (species list xlsx; trees only on iTOL `itol.embl.de/shared/dolfweijers`) | OA | 4 |
| P6 | GID1 receptor within the CXE carboxylesterase family (+DELLA) | bacteria→land plants | neofunctionalisation; family nested in larger family; WGD | Yoshida 2018 Dataset S2 (260 CXE/GID1 genes from 70 spp. with locus/accession) **[VERIFIED-DATA]** | OA | 4 |
| P7 | bHLH (26 subfamilies) | green plants (+ outgroups) | TD, WGD, DOM (short domain) | Pires & Dolan 2010 File009 (566 rows: species, protein, locus) + FASTA **[VERIFIED-DATA]** | OA | 4 |
| P8 | PEBP (MFT / FT / TFL1) | land plants | WGD (seed-plant), TD | Karlgren 2011 (PW) | PW | 3 |
| P9 | LEAFY | land plants (+ streptophyte algae) | SC (control), WGD retention in polyploids | Sayou 2014; Maizel 2005 (Science, PW) | PW | 3 |
| P10 | NLR / NB-LRR (TNL, CNL, RNL) | angiosperms (RefPlantNLR: 78 host spp.) | **DOM**, TD, FUS (integrated domains), LOSS (TNL loss in monocots), FAST | RefPlantNLR S1 Table (481 entries, GenPept, locus, subclass) + NB-ARC Newick **[VERIFIED-DATA]**; Zenodo 10.5281/zenodo.3936022; pan-NLRome GitHub | OA | **5** |
| P11 | Receptor-like kinases (RLK/Pelle) and plant kinome | land plants | DOM (kinase), TD, WGD | Shiu & Bleecker 2001 (≥610 At RLKs; PMC); Shiu 2004 (PW); Lehti-Shiu 2009 (PW) | PMC/PW | 4 |
| P12 | PPR proteins (P, PLS, E, DYW) | land plants | DOM (repeat), TD, LOSS (DYW/editing loss) | Lurin 2004 (441 At; PW); Cheng 2016 (PW); Gutmann 2020 (PW) | PW | 3 |
| P13 | F-box superfamily | Arabidopsis/poplar/rice | DOM, extreme TD (>100× expansion) | Xu 2009 PNAS (PMC; 692/337/779 genes, 42 families) | PMC | 3 |
| P14 | Terpene synthases (TPS-a…h) + microbial-type TPS | land plants | TD, HGT (MTPSL in non-seed plants) | Chen 2011 (PW); Jia 2016 PNAS (PMC) | PW/PMC | 3 |
| P15 | Expansins (EXPA, EXPB, EXLA, EXLB) | land plants (+ bacterial EXLX) | WGD, TD, HGT (bacterial expansins, UNVERIFIED citation) | Sampedro & Cosgrove 2005 (OA; alignment as PDF only) | OA | 3 |
| P16 | Major intrinsic proteins / aquaporins (PIP, TIP, NIP, SIP, XIP, GIP, HIP) | land plants (+ bacterial AqpN) | HGT (NIP origin), LOSS (XIP), WGD | Danielson & Johanson 2008 Nexus trees with TRANSLATE blocks (gene names) **[VERIFIED-DATA]**; Zardoya 2002 (PMC) | OA | 4 |
| P17 | LHC antenna superfamily (LHCA/B, LHCSR, ELIP, OHP, SEP, PsbS; FCP/LHCR) | cyanobacteria (HLIPs) + all plastid lineages | EGT, secondary-plastid transfer, TD | Engelken 2010 (OA; tables in PDF only) | OA | 4 |
| P18 | Plant cytochrome P450s | land plants (11 ancestral angiosperm clans) | TD, WGD; nomenclature at 40/55 % identity (not strictly phylogenetic) | Nelson P450 site (curated names) `drnelson.uthsc.edu/plants/` | site | 3 |
| P19 | AP2/ERF | Arabidopsis/rice (122/139 ERF genes) | TD, WGD, DOM | Nakano 2006 (PW) | PW | 2 |
| P20 | R2R3-MYB | land plants | TD, WGD, DOM | Dubos 2010 (PW) | PW | 2 |
| P21 | GAPDH (GapA/B plastid vs GapC cytosolic) + Calvin-cycle EGT | bacteria, plants, algae | EGT, HGT, replacement | Figge 1999 (PW); Martin 2002 (PMC) | PW/PMC | 4 |
| X1 | RecA / RadA / Rad51 / Dmc1 (+ Rad51 paralogs) | all three domains | LUCA, LECA dup (Rad51/Dmc1), EGT (plant organellar RecA) | Lin 2006 PNAS (PMC; SI at PNAS, not checked) | PMC | **5** |
| X2 | Ribosomal proteins | 995 prokaryotic genomes | **SC (negative control)**, rare paralogs (Zn-ribbon variants), LOSS | Yutin 2012 Tables S1/S3 (r-protein list w/ COG; copy-number matrix) **[VERIFIED-DATA]** | OA | **5** |
| X3 | Translation factors EF-Tu/EF-1α vs EF-G/EF-2 | all domains | LUCA paralogs (rooting), EGT (tufA→nucleus) | Baldauf 1996 (PMC); Iwabe 1989 (PMC); Baldauf 1990 (PMC) | PMC | 4 |
| X4 | F/V/A-ATPase catalytic & non-catalytic subunits | all domains | LUCA paralogs, HGT (V-type into bacteria), ORG | Gogarten 1989 (PMC); Mulkidjanian 2007 (PW) | PMC/PW | 3 |
| X5 | Nitrogenase NifHDK/VnfDGK/AnfDGK + Bch/Chl paralogs | bacteria + archaea | HGT, paralog vs ortholog, LOSS | Garcia 2020 GitHub alignment (319 seqs labelled Nif/Anf/Vnf/BchChl + species) **[VERIFIED-DATA]**; Raymond 2004 (PW); Boyd 2015 (PMC) | OA | 4 |
| X6 | β-lactamases (classes A, C, D serine; B1/B2/B3 metallo) | bacteria | **HGT** (plasmids), **CONV** (class B unrelated fold), allelic families | NCBI Reference Gene Catalog (6,128 β-lactamase entries with class, family, RefSeq/GenBank accession, PMIDs) **[VERIFIED-DATA]**; BLDB (Naas 2017, OA) | OA | 4 |
| X7 | B-family DNA polymerases (Pol α, δ, ε, ζ; archaeal PolB1–3) | archaea + eukaryotes | **FUS** (chimeric origin, inactivated Pol ε module), archaeal–eukaryotic, LECA | Tahirov 2009 alignment + Newick trees (GI\|taxon labels) **[VERIFIED-DATA]**; Makarova 2014 (OA) | OA | 4 |
| X8 | Two-component systems (HPK / RR) | 200+ bacteria | DOM, HGT vs lineage-specific expansion | Alm 2006 Dataset S4 (FASTA of all HPK domains), S1 per-genome table **[LISTED]** | OA | 3 |
| X9 | Chaperonins GroEL / thermosome / CCT(α–θ) | all domains | LECA paralogs (CCT), recurrent archaeal paralogy | Archibald 1999 (PW); Lund 2009 (PW) | PW | 3 |
| X10 | Rubisco forms I–IV + RLPs | all domains | HGT (form II to dinoflagellates), ORG (rbcL plastid-encoded) | Tabita 2008 (PW); Banda 2020 (PW); Schulz 2022 (PW) | PW | 2 |
| X11 | CRISPR-Cas (Cas1, effectors) | bacteria + archaea | extreme HGT, modularity | Makarova 2020 (PMC); Makarova 2025 (OA) | PMC/OA | 2 |
| X12 | Asgard ESPs / archaeal–eukaryotic families (ESCRT, small GTPases…) | Asgard + eukaryotes | archaeal–eukaryotic, FAST | Eme 2023 (OA); Zaremba-Niedzwiedzka 2017 (PW) | OA/PW | 3 |

(E16–E17, P19–P20, X9–X12 are covered briefly in §3. Rubisco and CRISPR are rated low because plastid encoding or HGT modularity makes a "reference orthology" poorly defined for nuclear-proteome benchmarks.)

---

## 2. Detailed entries — top-priority candidates

### E1. Kinesin superfamily — **Priority 5**
* **Papers.** Wickstead B, Gull K, Richards TA (2010) *Patterns of kinesin evolution reveal a complex ancestral eukaryote with a multifunctional cytoskeleton.* BMC Evol Biol 10:110. PMID 20423470, PMC2867816, doi:10.1186/1471-2148-10-110 (**OA**). Wickstead B, Gull K (2006) *A "holistic" kinesin phylogeny…* Mol Biol Cell. PMID 16481395, PMC1415282, doi:10.1091/mbc.e05-11-1090 (PMC).
* **Data [VERIFIED-DATA].** Add. file 3 (`1471-2148-10-110-S3.XLS`) has 1,626 rows with Analysis ID, genome-project ID, GenBank accession (595 rows; the rest use genome-project IDs), Alternative source, **Family**, **Subfamily**, exclusion flags and domain architecture. Add. files 8/9 are full and trimmed motor-domain alignments (1,458 / 1,263 seqs; headers like `Galga_ENSGALP00000010647`). Add. file 7 lists proteome sources and versions. 45 species. Per-species counts: Paramecium 133, human 99, poplar 86, Tetrahymena 82, Arabidopsis 69, Physcomitrella 68 … Theileria 2.
* **Phenomena.** 51 anciently derived paralogs; ≥11 families in LECA (abstract). Massive lineage-specific expansions (ciliates, land-plant Kinesin-14 and -7). Losses correlated with cilia (Kinesin-2/-9/-16/-17). Promiscuous domain architectures (105).
* **Species.** The 45 eukaryotes (old proteome versions). Use current proteomes of the same species and map by sequence.
* **Difficulty.** High. A shared motor domain links ~20 families, so graph methods merge families or split by architecture. Hierarchical levels span LECA → clade-specific. Some labels are "ND"/"ungrouped"/parenthesised (uncertain), so keep only firm labels.
* **Why 5.** Best-in-class machine-readable, phylogeny-based labels across all supergroups, with explicit LECA inference.

### E2. ARF GTPase family (Arf, Arl, Sar, SR-β, TRIM23) — **Priority 5**
* **Paper.** Vargová R, Wideman JG, Derelle R, Klimeš V, Kahn RA, Dacks JB, Eliáš M (2021) GBE 13:evab157. PMID 34247240, PMC8358228, doi:10.1093/gbe/evab157 (**OA**).
* **Data [VERIFIED-DATA].** Suppl. Table 1 has 2,286 rows from 115 species. Columns include *Annotation* (31 orthologous-group labels: Arf1 442, Sar1 160, Arl1 116, SRβ 116, Arl2 115, Arfrp1 112, Arl3 107, Arl8 106 … 153 "unannotated"), *Annotation basis* (HMM / phylogeny / "IQ"), *Conclusively annotated* flag, *Accession number(s)*, and gene-model corrections. Suppl. dataset 1 is a FASTA of 2,054 sequences. The authors report ~16 ARF members in LECA.
* **Difficulty.** High. Short (~180 aa), very similar paralogs; many protist-specific long-branch paralogs; many corrected gene models (an annotation-error phenomenon in itself). Use the "conclusively annotated" subset as gold. Rows annotated only by HMM are tool-assisted.
* **Why 5.** Curated, eukaryote-wide, explicit LECA orthogroups, accessions present. Pairs naturally with Rab (E3) as the "small-GTPase" test.

### E3. Rab GTPases — **Priority 5**
* **Papers.** Klöpper TH, Kienle N, Fasshauer D, Munro S (2012) BMC Biol 10:71. PMID 22873208, PMC3425129, doi:10.1186/1741-7007-10-71 (**OA**). Elias M, Brighouse A, Gabernet-Castello C, Field MC, Dacks JB (2012) J Cell Sci 125:2500. PMID 22366452, PMC3383260, doi:10.1242/jcs.101378 (PMC; LECA ≈23 Rabs).
* **Data [VERIFIED-DATA].** Add. files 4–9 are FASTA alignments of Rab groups I–VI: 1,669 sequences from 81 species. Names encode species and Rab subfamily (e.g. `AuAn_Rab8B_2`; most frequent: Rab18 59, Rab6 52, Rab8 52…). Add. file 3 lists genome sources. The trees are on TreeBASE S12811 (the PURL redirects; I did not confirm that the TreeBASE host is online). The authors' Rab DB is offline. Note: there are **no accessions** in the alignment names; map by sequence identity.
* **Difficulty.** Very high. Short paralogs, extreme lineage expansions (e.g. Trichomonas, Entamoeba; UNVERIFIED counts), frequent losses. Classification was HMM-assisted from curated seeds, so it is partly tool-derived.
* **Why 5.** It is the canonical eukaryote-wide "many LECA paralogs + expansions" family and is strongly complementary to E2.

### E4. Tubulin superfamily — **Priority 5**
* **Papers.** Findeisen P, Mühlhausen S, Dempewolf S, Hertzog J, Zietlow A, Carlomagno T, Kollmar M (2014) GBE 6:2274. PMID 25169981, PMC4202323, doi:10.1093/gbe/evu187 (**OA**; 3,524 tubulins, 504 species). Outgroup: Santana-Molina C, Del Saz-Navarro D, Devos DP (2022) *Early origin and evolution of the FtsZ/tubulin protein family.* Front Microbiol 13:1100249. PMID 36704558, PMC9871819 (**OA**; Data_Sheet_2.zip **[LISTED]**).
* **Data [VERIFIED-DATA].** `tubulin_alignment.fas` has 3,286 sequences whose names encode the subgroup (aTub 330 tags, bTub 216, gTub 210, dTub 122, eTub 104, iTub 91, zTub 8; many plain `Tub` names need the paper's mapping). `DatasetS2.xls` gives the 504-row species table with abbreviation and full NCBI lineage (includes *Prosthecobacter* BtubA/B). CyMoBase holds the sequences (site up).
* **Difficulty.** Medium. α/β/γ are slow-evolving and near-universal (good control), while δ/ε/ζ are lost in many lineages together with centriole loss. Extensive *recent* duplications (paper title) test fine-grained in-paralog handling.
* **Why 5.** Cleanly labelled, ultra-broad, mixes a slow core with patchy paralogs, and offers a prokaryotic outgroup (deep homology without orthology).

### E5. Myosins — **Priority 5**
* **Papers.** Sebé-Pedrós A, Grau-Bové X, Richards TA, Ruiz-Trillo I (2014) GBE 6:290. PMID 24443438, PMC3942036, doi:10.1093/gbe/evu013 (**OA**). Kollmar M, Mühlhausen S (2017) BMC Evol Biol 17:211. PMID 28870165, PMC5583752, doi:10.1186/s12862-017-1056-2 (**OA**). Odronitz F, Kollmar M (2007) Genome Biol 8:R196. PMID 17877792, PMC2375034, doi:10.1186/gb-2007-8-9-r196 (**OA**; 2,269 manually annotated myosins, 328 species; Add. file 5 is FASTA of all sequences **[LISTED]**).
* **Data.** Sebé-Pedrós `DataS1.fasta` has 725 sequences **[VERIFIED-DATA]**, but its headers do **not** carry class labels; the class assignments are in the supplementary PDF. Kollmar 2017 figshare `10.6084/m9.figshare.4565155` (CC BY 4.0) contains `Myosin_motordomain_alignment.fas` (28 MB), `Phylogenetic_trees.zip` and an intron-pattern alignment (listed via the figshare API). CyMoBase is the curated source for names (e.g. class-based naming).
* **Difficulty.** High. Multi-domain tails (MyTH4-FERM, SH3, PH, kinase…) with domain shuffling. 30+ classes, many lineage-specific. Long genes, so gene-model errors.
* **Why 5.** Manually annotated, eukaryote-wide, classic domain-architecture test.

### E7. MCM2-9 helicases — **Priority 5**
* **Paper.** Liu Y, Richards TA, Aves SJ (2009) *Ancient diversification of eukaryotic MCM DNA replication proteins.* BMC Evol Biol 9:60. PMID 19292915, PMC2667178, doi:10.1186/1471-2148-9-60 (**OA**).
* **Data [VERIFIED-DATA].** Add. files 5 and 7 are Newick treefiles (275 labelled leaves in file 7, names like `TbMcm9`, `PtMcm9`, archaeal `ThvoMcm`), rooted with archaeal MCMs. Add. file 6 lists genomes and proteomes (.doc). Labels encode species and paralog but **no accessions**, so map by sequence (via the PDF trees) or re-identify.
* **Phenomena.** Six LECA paralogs (MCM2-7) arose by duplications from an archaeal MCM, so one copy per paralog per genome (near-single-copy). MCM8/9 show repeated independent losses (e.g. absent from many fungi; UNVERIFIED detail).
* **Difficulty.** Low to moderate, which makes it an ideal "easy, deep" calibrator. Tests whether methods resolve 6 deep paralogs as separate HOGs at the LECA level instead of one HOG.

### E12. Aminoacyl-tRNA synthetases — **Priority 5**
* **Papers.** Woese CR, Olsen GJ, Ibba M, Söll D (2000) MMBR 64:202. PMID 10704480, PMC98992 (PMC). Brindefalk B et al. (2007) *Origin and evolution of the mitochondrial aminoacyl-tRNA synthetases.* MBE 24:743–756. PMID 17182897, doi:10.1093/molbev/msl202 (**PW**). Duchêne AM et al. (2005) PNAS *Dual targeting is the rule for organellar aaRS in Arabidopsis.* PMID 16251277, doi:10.1073/pnas.0504682102 (PW per EPMC). Chaliotis A et al. (2017) *The complex evolutionary history of aminoacyl-tRNA synthetases.* NAR 45:1059–1068. PMID 28180287, PMC5388404, doi:10.1093/nar/gkw1182 (**OA**; `gkw1182_Supp.zip` **[LISTED]**).
* **Phenomena.** LUCA-level class I/II paralogy; rampant inter-domain HGT; mitochondrial paralogs of α-proteobacterial or other origin; plastid/mitochondrial dual targeting in plants (organellar aaRS that are cytosolic in some species). This is the archetypal **EGT/HGT** benchmark family.
* **Difficulty.** Very high. Within one species, the cytosolic and organellar enzymes for the same amino acid are *not* orthologous to their namesakes elsewhere (xenologs).
* **Data gap.** No single labelled supplement confirmed yet. Chaliotis Supp should be opened next. The Arabidopsis organellar aaRS list from Duchêne 2005 is paywalled, so request it.

### P1. MADS-box (MIKC^c) — **Priority 5** (central WGD test)
* **Papers.** Parenicová L et al. (2003) Plant Cell 15:1538. PMID 12837945, doi:10.1105/tpc.011544 (PW in EPMC; 107 Arabidopsis MADS genes). Gramzow L, Theissen G (2010) *A hitchhiker's guide to the MADS world of plants.* Genome Biol 11:214. PMID 20587009, PMC2911102 (**OA**). Event-resolving subfamily papers (all Genetics, in PMC):
  * Litt & Irish 2003, AP1/FUL (γ-triplication origin of euAP1/euFUL): PMID 14573491, PMC1462802.
  * Kramer, Dorit & Irish 1998, AP3/PI (paleoAP3 → euAP3/TM6): PMID 9611190, PMC1460198.
  * Kramer, Jaramillo & Di Stilio 2004, AGAMOUS subfamily (AG/PLE, SHP via core-eudicot duplication): PMID 15020484, PMC1470751.
  * Zahn LM et al. 2005, SEPALLATA: PMID 15687268, PMC1449606.
  * Preston & Kellogg 2006, grass AP1/FUL paralogs (ρ WGD): PMID 16816429, PMC1569798.
  * Schranz ME et al. 2002, FLC copies in *Brassica rapa* (Brassica triplication): PMID 12454088, PMC1462321.
  * Ruelens P et al. 2013, FLC in monocots (tandem origin): PMID 23955420 (PW).
  * Thangavel G, Nayar S 2018, MIKC in non-seed plants: PMID 29720991, PMC5915566 (**OA**; no supplement found by the scan).
  * Becker A, Theissen G 2003 MPE (PW); Nam J et al. 2003 MBE (PW); Gramzow & Theißen 2015 J Exp Zool B (PW).
* **Data.** No single machine-readable labelled dataset was confirmed. Subfamily membership for Arabidopsis, rice and other genomes is defined in tables of these papers (Parenicová Table 1 lists the 107 At genes and their groups). Expect manual curation from PDFs: ~10–15 subfamilies × ~15 genomes.
* **Phenomena.** The textbook WGD family. γ (core eudicots) generated euAP1/euFUL, euAP3/TM6 and AG/SHP-like pairs. ρ (grasses) produced FUL1/FUL2-type pairs. α (Brassicaceae) produced SHP1/SHP2-type pairs (UNVERIFIED assignment; check against Kramer 2004). The Brassica triplication produced the BrFLC copies. Type I (Mα/β/γ) genes evolve fast with tandem duplication, a good contrast.
* **Species.** Arabidopsis, *B. rapa*, papaya or grape (γ-only reference), tomato, rice, maize/sorghum (ρ), *Amborella*, a gymnosperm, *Physcomitrium*, *Marchantia*, *Selaginella*, a charophyte.
* **Difficulty.** High. Short conserved MADS domain plus a fast K/C region; WGD-derived paralogs close in age to the speciation nodes (γ predates the asterid/rosid split).
* **Request list.** Parenicová 2003 full text and supplements; Becker & Theissen 2003.

### P2. Phytochromes (+ neochrome) — **Priority 5**
* **Papers.** Li F-W, Melkonian M, Rothfels CJ, Villarreal JC, Stevenson DW, Graham SW, Wong GK-S, Pryer KM, Mathews S (2015) *Phytochrome diversity in green plants and the origin of canonical plant phytochromes.* Nat Commun 6:7852. PMID 26215968, PMC4525182, doi:10.1038/ncomms8852 (**OA**). Li F-W et al. (2014) *Horizontal transfer of an adaptive chimeric photoreceptor from bryophytes to ferns.* PNAS 111:6672. PMID 24733898, PMC4020063 (PMC). Mathews S et al. (2000) Am J Bot (PW).
* **Data [LISTED].** Dryad `doi:10.5061/dryad.5rs50` (CC0) contains `alignments_and_trees.zip` (2.0 MB) and `baits-120-60.fas` (file list confirmed via the Dryad API). **The binary download is blocked to scripts (401/403)**, so a browser download or API token is needed. The supplementary PDF has Tables 1–2.
* **Phenomena.** Seed-plant duplications (PHYA/PHYC vs PHYB/E, angiosperm PHYB/PHYE) linked to ζ/ε WGD-era events (interpretation to be checked against the paper); independent duplications in ferns, mosses and zygnematophytes; phytochrome loss in some algae; neochrome (PHY–PHOT chimera, i.e. **fusion**) transferred **horizontally** hornwort → fern.
* **Species.** Many are 1KP transcriptomes, so for a genome panel use the subset with genomes (Arabidopsis, rice, *Amborella*, *Picea*/*Gnetum*, *Azolla*/*Salvinia*, *Ceratopteris*, *Physcomitrium*, *Marchantia*, *Anthoceros*, *Mesotaenium*/*Spirogloea*, *Klebsormidium*, *Chlamydomonas*).
* **Difficulty.** Medium-high. Large multi-domain protein, and the neochrome fusion should be split or attached correctly by methods.

### P3. Plant homeobox (14 classes) — **Priority 5**
* **Paper.** Mukherjee K, Brocchieri L, Bürglin TR (2009) MBE 26:2775. PMID 19734295, PMC2775110, doi:10.1093/molbev/msp201 (**OA**).
* **Data [VERIFIED-DATA].** `msp201_5.html` is a table with Gene number, Gene name, **Accession No.**, alternative names (incl. AGI codes), species code, chromosome, **Gene class** (HD-ZIP I–IV, BEL, KNOX, WOX, PLINC, PHD, DDT, NDX, LD, SAWADEE, PINTOX), notes and protein sequence. 10 genomes: *Cyanidioschyzon* (red alga), *Chlamydomonas*, *Ostreococcus*, *Physcomitrella*, *Selaginella*, Arabidopsis, poplar, rice, sorghum, maize (genome list per abstract: red alga, green algae, moss, *Selaginella*, flowering plants). PDFs hold the trees (msp201_1.pdf is 70 MB).
* **Phenomena.** All 14 classes were present in the moss–vascular ancestor. KNOX/BEL (TALE) heterodimer partners are ancient (algae). Grass WGD (ρ) and maize allotetraploidy give recent paralogs. Homeodomain promiscuity forces separation of 14 classes that share one domain.
* **Difficulty.** High at class boundaries (the shared ~60-aa homeodomain), moderate within classes.
* **Why 5.** Machine-readable, accession-bearing, explicit class labels, and a ready plant panel that includes maize (recent allopolyploid).

### P10. NLR immune receptors — **Priority 5** (hardest plant family)
* **Papers.** Kourelis J, Sakai T, Adachi H, Kamoun S (2021) *RefPlantNLR…* PLoS Biol 19:e3001124. PMID 34669691, PMC8559963, doi:10.1371/journal.pbio.3001124 (**OA**). Van de Weyer A-L et al. (2019) *A species-wide inventory of NLR genes and alleles in Arabidopsis thaliana.* Cell 178:1260. PMID 31442410, PMC6709784, doi:10.1016/j.cell.2019.07.038 (**OA**). Shao Z-Q et al. (2016) Plant Physiol 170:2095. PMID 26839128, doi:10.1104/pp.15.01487 (PW; >6,000 NLRs in 22 angiosperms; 7 TNL + 14 CNL + 2 RNL lineages in the ancestral angiosperm).
* **Data [VERIFIED-DATA].** RefPlantNLR S1 Table has 481 experimentally validated NLRs from 78 host species with GenPept, CDS, Locus, **Subclass**, Structure (CNL 294, TNL 90, RNL 7, …) and PMID of the functional evidence. S15 Dataset is a Newick NB-ARC ML tree. S13 holds animal/bacterial/archaeal NB-ARC outgroups. Zenodo `10.5281/zenodo.3936022` and GitHub `JKourelis/refplantnlR` mirror it. The pan-NLRome GitHub (`weigelworld/pan-nlrome/results/orthogroups/nlrome.refined.groups`) gives orthogroups across 64 accessions (tool-derived, then manually refined; use as secondary).
* **Phenomena.** Promiscuous NB-ARC and LRR domains; integrated decoy domains (fusion); massive tandem clusters with gene conversion; birth–death; TNL loss in monocots (and in some eudicot lineages; UNVERIFIED specifics); presence/absence polymorphism.
* **Difficulty.** Extreme. Shao 2016's ancestral-lineage framework (23 ancestral angiosperm lineages) is the right level for HOG evaluation, but its data are paywalled, so **request**. RefPlantNLR gives functional anchors (named genes such as RPS4/RRS1, ZAR1, NRC; UNVERIFIED that all are present) to test whether known orthologs (e.g. ZAR1 across eudicots) are recovered.

### X1. RecA / RadA / Rad51 / Dmc1 — **Priority 5**
* **Paper.** Lin Z, Kong H, Nei M, Ma H (2006) *Origins and evolution of the recA/RAD51 gene family: evidence for ancient gene duplication and endosymbiotic gene transfer.* PNAS 103:10328. PMID 16798872, PMC1502457, doi:10.1073/pnas.0604232103 (PMC; Europe PMC has no full text, SI at PNAS not checked).
* **Phenomena.** LUCA-depth family (bacterial RecA vs archaeal RadA). Rad51/Dmc1 duplication in LECA; Rad51 paralogs (Rad51B/C/D, XRCC2/3) are ancient and fast. Plant nuclear genes for plastid- and mitochondrion-targeted RecA come from **EGT**. Dmc1 is lost in several lineages (e.g. *Drosophila*, *C. elegans*; UNVERIFIED here).
* **Difficulty.** Medium (core) to very high (Rad51 paralogs, fast). Small enough (~5–15 genes per eukaryote) to curate by hand.

### X2. Ribosomal proteins — **Priority 5** (negative control)
* **Paper.** Yutin N, Puigbò P, Koonin EV, Wolf YI (2012) *Phylogenomics of prokaryotic ribosomal proteins.* PLoS One 7:e36972. PMID 22615861, PMC3353972, doi:10.1371/journal.pone.0036972 (**OA**).
* **Data [VERIFIED-DATA].** Table S1 lists 102 r-protein families with COG IDs and names (universal 34; bacterial/archaeal-specific). Table S3 is a copy-number matrix (995 genomes × r-proteins; most cells = 1). Files S1/S4/S6 are Newick trees. S3 also exposes the known exceptions: paralogous zinc-ribbon/zinc-less variants and duplicates in some genomes (File S3 "Distribution of paralogous r-proteins in bacteria").
* **Use.** Near-single-copy orthogroups with known species membership. Any split or merge is a method error (a hierarchical "speciation-only" control). Eukaryotic r-proteins (cytosolic plus mitochondrial/plastid paralogs of different origin) add an EGT layer: the nuclear-encoded organellar r-proteins of plants are bacterial-derived and not orthologous to cytosolic ones.
* **Difficulty.** Low by design, except short proteins (L34, L36 ~40 aa), which are often missed by annotation.

---

## 3. Brief entries — remaining candidates

**E6 Dynein heavy chains — P4.** Wickstead & Gull 2007 *Traffic* 8:1708 (PMID 17897317, PMC2239267, **OA**; File S2 `.xls` = IDs of DHC/IC/Tctex repertoires in 24 eukaryotes **[LISTED]**). Kollmar 2016 MBE 33:3249 (PMID 27880711, PMC5100056, **OA**; Supplementary_Data1/2 = Newick trees with ~3,000 labelled DHC leaves such as `CrDHC8`, `PhpDHC8A` **[VERIFIED-DATA]**, no accessions; sequences in CyMoBase/diArk). Phenomena: LECA paralogs (9 DHC classes; UNVERIFIED count); complete loss of axonemal dyneins in angiosperms and other non-ciliated lineages, a clean **loss** test; ~4,500-aa proteins so annotation splits are common (fission artefacts).

**E8 SMC — P4.** Yoshinaga M, Inagaki Y (2021) GBE 13:evab256 (PMID 34894224, PMC8665677, **OA**). Dryad `10.5061/dryad.ghx3ffbq3` (CC0, v4) holds `Sequence_files.zip`, `Treefiles.zip`, `Yoshinaga_512_arcSMC_list.xlsx`, `Yoshinaga_SMC_Rad50SbcC_accessary_list.xlsx` (accession lists, per the README description) **[LISTED; download blocked to scripts]**. Cobbe & Heck 2004 MBE (PMID 14660695, PW). SMC1-6 are near-single-copy LECA paralogs, with prokaryotic SMC/Rad50 as deep outgroups.

**E9 Adaptor complexes / COPI / TSET — P4.** Hirst J et al. (2014) eLife 3:e02866 (PMID 24867644, PMC4031984, **OA**): Fig3-source-data-1 has 1,413 rows (species, annotation, accession for 874) **[VERIFIED-DATA]**, plus alignments (supp1.zip). Hirst J et al. (2011) *The fifth adaptor protein complex.* PLoS Biol 9:e1001170 (PMID 22022230, PMC3191125, **OA**; Table S1 = all sequences with taxon/genome URLs **[LISTED]**). Phenomena: subunits of AP1–5/COPI/TSET are paralogous (β/γ/α large subunits, μ/σ), a "complex-level paralogy" test with many independent losses (AP-5, TSET).

**E10 SNAREs — P3.** Kloepper TH, Kienle CN, Fasshauer D (2007) MBoC 18:3463 (PMID 17596510, PMC1951749, PMC-only; 145 species, 20 ancestral subclasses). Kienle N, Kloepper TH, Fasshauer D (2009) BMC Evol Biol 9:19 (PMID 19166604, PMC2639358, **OA**; Add. file 3 = Nexus trees and alignments **[LISTED]**). SNARE DB is offline. Short coiled-coil SNARE motif means weak signal.

**E11 Histones — P3.** Talbert PB et al. (2012) Epigenetics Chromatin 5:7 (PMID 22650316, PMC3380720, **OA**; nomenclature, no data supplement). Draizen EJ et al. (2016) HistoneDB 2.0, Database baw014 (PMID 26989147, PMC4795928, **OA**; web DB live at ncbi.nlm.nih.gov/research/histonedb; curated seeds plus HMM-classified NR sequences, so tool-assisted). Phenomena: ultra-conserved H3/H4 (no phylogenetic signal); identical replication-dependent copies (in-paralogs); recurrent independent origins of variants (H3.3-like and CenH3 per lineage), which is **convergence**. A good "methods should not over-resolve" test.

**E13 EF-1α / EFL — P4.** Keeling PJ, Inagaki Y (2004) PNAS 101:15380 (PMID 15492217, PMC524445, PMC). Gile GH et al. (2009) J Eukaryot Microbiol (PMID 19602082, PW). Mutually exclusive EF-1α vs EFL across protists and green algae, with multiple HGT of EFL. Small dataset; test of xenolog detection.

**E14 Photolyase/cryptochrome — P4.** Mei Q, Dvornyk V (2015) PLoS One 10:e0135940 (PMID 26352435, PMC4564169, **OA**). S1 Table (docx) lists CRY/PHR sequences with GenBank nucleotide and protein accessions, taxon and figure designation **[VERIFIED-DATA]**. Subfamilies: CPD I/II/III, (6-4)/animal CRY, CRY-DASH, plant CRY. Spans bacteria to eukaryotes, with HGT and many losses. Plant photoreceptor angle: plant CRY1/CRY2 come from the Brassicaceae/seed-plant duplication (check).

**E15 Proteasome regulators — P3.** Fort P, Kajava AV, Delsuc F, Coux O (2015) GBE 7:1363 (PMID 25943340, PMC4453063, **OA**; Supplemental_Tables.xls **[LISTED]**).

**E16 ARPs — P3.** Muller J et al. (2005) MBoC 16:5736 (PMID 16195354, PW in EPMC). Arp1–10 are LECA paralogs of actin, with fast-evolving lineage-specific ARPs. **Request.**

**E17 ESCRT — P3.** Leung KF, Dacks JB, Field MC (2008) Traffic 9:1698 (PMID 18637903, PW). Hatano T et al. (2022) Nat Commun 13:3398 (PMID 35697693, PMC9192718, **OA**; MOESM6 xlsx **[LISTED]**). Archaeal–eukaryotic (Asgard) origin.

**P4 Class III HD-ZIP — P3.** Floyd SK, Zalewski CS, Bowman JL (2006) Genetics 173:373 (PMID 16489224, PMC1461458, PMC). Prigge MJ, Clark SE (2006) Evol Dev (PMID 16805899, PW). Small (5 in Arabidopsis: REV, PHB, PHV, CNA, ATHB8). Seed-plant and angiosperm duplications define clean subclades. Also covered by P3.

**P5 Auxin ARF/Aux-IAA/TIR1-AFB — P4.** Mutte SK et al. (2018) eLife 7:e33399 (PMID 29580381, PMC5873896, **OA**). Supp file 1 = species xlsx; Supp file 2 = alignments (in **pptx**, so not machine-friendly); full trees only on iTOL (`itol.embl.de/shared/dolfweijers`, not tested for export). Parry G et al. (2009) PNAS TIR1/AFB (PMID 20018756, PMC2799741). Finet C et al. (2010) Plant J ARFs (PMID 20626651, PW). Phenomena: A/B/C ARF classes originate in charophytes; proto-ARF; WGD duplicates in angiosperms (e.g. TIR1/AFB1, AFB2/AFB3; Brassicaceae α; UNVERIFIED pairing).

**P6 GID1 within CXE — P4.** Yoshida H et al. (2018) PNAS 115:E7844 (PMID 30068603, PMC6099883, **OA**). Dataset S2 lists 260 CXE and GID1 genes from 70 species, from *E. coli* to angiosperms, with "Locus name" (NCBI WP_/XP_ or genome locus) **[VERIFIED-DATA]**. Datasets S1/S3/S4 are alignments. Phenomenon: a receptor nested inside a large enzyme family, so the family boundary must be set at the right HOG level. Pair with DELLA: Hernández-García J et al. (2019) MBE 36:908 (PMID 30668817, PW; **request**).

**P7 bHLH — P4.** Pires N, Dolan L (2010) MBE 27:862 (PMID 19942615, PMC2839125, **OA**). File009 = 566 rows (Arabidopsis 158, rice 173, *Selaginella* 103, *Physcomitrella* 98, algae and outgroups); File011 = FASTA **[VERIFIED-DATA]**. 26 subfamilies. The short domain gives weak signal, so subfamily assignment is hard, which makes this a good stress test.

**P8 PEBP — P3.** Karlgren A et al. (2011) Plant Physiol 156:1967 (PMID 21642442, PW). MFT is ancestral; FT/TFL1-like came from a seed-plant duplication; many WGD copies (e.g. soybean FT homologs; UNVERIFIED). **Request.**

**P9 LEAFY — P3.** Sayou C et al. (2014) Science 343:645 (PMID 24436181, PW). Maizel A et al. (2005) Science 308:260 (PMID 15821093, PW). Single copy in most diploid angiosperms, with retained WGD copies in polyploids and gymnosperm LFY/NEEDLY paralogs. A single-copy positive control with known exceptions.

**P11 RLK/Pelle — P4.** Shiu S-H, Bleecker AB (2001) PNAS 98:10763 (PMID 11526204, PMC58549; ≥610 At members). Shiu S-H et al. (2004) Plant Cell 16:1220 (PMID 15105442, PW). Lehti-Shiu MD et al. (2009) Plant Physiol 150:12 (PMID 19321712, PW). Lehti-Shiu & Shiu (2012) Phil Trans B 367:2619 (PMID 22889912, PMC3415837). Subfamily lists (e.g. LRR-I…XV) are in the PW papers' supplements, so **request** Shiu 2004 and Lehti-Shiu 2009 supplements.

**P12 PPR — P3.** Lurin C et al. (2004) Plant Cell 16:2089 (PMID 15269332; 441 At PPRs). Cheng S et al. (2016) Plant J 85:532 (PMID 26764122). Gutmann B et al. (2020) Mol Plant 13:215 (PMID 31760160). All PW. Extreme repeat family. Loss of PLS/DYW correlates with loss of organellar editing (e.g. *Marchantia*).

**P13 F-box — P3.** Xu G, Ma H, Nei M, Kong H (2009) PNAS 106:835 (PMID 19126682, PMC2630105). 692/337/779 genes in Arabidopsis/poplar/rice; 42 families; >100-fold expansions.

**P14 TPS — P3.** Chen F, Tholl D, Bohlmann J, Pichersky E (2011) Plant J 66:212 (PMID 21443633, PW; *Selaginella* 18, other genomes 40–152 TPS). Jia Q et al. (2016) PNAS 113:12328 (PMID 27791023, PMC5087002): microbial-type TPS in non-seed land plants (HGT).

**P15 Expansins — P3.** Sampedro J, Cosgrove DJ (2005) Genome Biol 6:242 (PMID 16356276, PMC1414085, **OA**; alignment only as PDF). Sampedro J, Carey RE, Cosgrove DJ (2006) J Plant Res (PMID 16411016, PW). Bacterial EXLX expansins (HGT) need a citation; that citation is UNVERIFIED and not yet found.

**P16 MIPs/aquaporins — P4.** Danielson JA, Johanson U (2008) BMC Plant Biol 8:45 (PMID 18430224, PMC2386804, **OA**): 23 moss MIPs in 7 subfamilies; 13 Nexus tree files with TRANSLATE blocks naming At, Zm and Pp genes (`ATPIP1.1`…) **[VERIFIED-DATA]**. Zardoya R et al. (2002) PNAS 99:14893 (PMID 12397183, PMC137515): NIPs originate by HGT from bacteria. Abascal F, Irisarri I, Zardoya R (2014) BBA (PMID 24355433, PW). Moderate size (~35 in angiosperms; UNVERIFIED), cleanly named subfamilies, and an HGT origin.

**P17 LHC superfamily — P4.** Engelken J, Brinkmann H, Adamska I (2010) BMC Evol Biol 10:233 (PMID 20673336, PMC3020630, **OA**; all tables inside one PDF, so curation by hand). Neilson JA, Durnford DG (2010) Photosynth Res (PMID 20596891, PW). Unique across photosynthetic eukaryotes (green, red, secondary plastids: diatom FCPs). Cyanobacterial one-helix HLIPs make it an EGT + secondary-endosymbiosis test.

**P18 Plant P450s — P3.** Nelson D, Werck-Reichhart D (2011) Plant J 66:194 (PMID 21443632, PW; 11 ancestral angiosperm genes). Nelson DR (2009) *The cytochrome P450 homepage.* Hum Genomics 4:59 (PMID 19951895, PMC3500189, **OA**). Curated nomenclature at drnelson.uthsc.edu/plants/ (site up). Names follow % identity thresholds, so they are not strict phylogenetic labels. Good for clan-level (A-type/non-A) HOGs.

**P19 AP2/ERF — P2.** Nakano T et al. (2006) Plant Physiol 140:411 (PMID 16407444, PW; 122 At / 139 rice ERFs, 12/15 groups).

**P20 R2R3-MYB — P2.** Dubos C et al. (2010) Trends Plant Sci (PMID 20674465, PW).

**P21 GAPDH / Calvin-cycle EGT — P4.** Figge RM, Schubert M, Brinkmann H, Cerff R (1999) MBE 16:429 (PMID 10331270, PW). Martin W et al. (2002) PNAS 99:12246 (PMID 12218172, PMC129430; thousands of cyanobacterial genes in the Arabidopsis nucleus). GapA/B (plastid, cyanobacterial origin) vs GapC (cytosolic) plus GapC1 replacements in algae make a compact EGT/replacement test. Related: Baldauf SL, Palmer JD (1990) *tufA transfer to the nucleus* (PMID 2314461, PW) and Baldauf, Manhart & Palmer 1990 (PMID 2371274, PMC54314).

**X3 EF-Tu/EF-G — P4.** Baldauf SL, Palmer JD, Doolittle WF (1996) PNAS 93:7749 (PMID 8755547, PMC38819). Iwabe N et al. (1989) PNAS 86:9355 (PMID 2531898, PMC298494). Ancient paralog pairs used to root the tree of life. Near-single-copy, but with mitochondrial/plastid paralogs in eukaryotes.

**X4 ATP synthase / V-ATPase — P3.** Gogarten JP et al. (1989) PNAS 86:6661 (PMID 2528146, PMC297905). Mulkidjanian AY et al. (2007) Nat Rev Microbiol 5:892 (PMID 17938630, PW). Some subunits are organelle-encoded (atpA/B in plastid; atp6/8/9 in mito), so nuclear-proteome panels miss them. Choose V-ATPase A/B (nuclear in eukaryotes) for the LUCA-paralogy test.

**X5 Nitrogenase — P4.** Garcia AK, McShea H, Kolaczkowski B, Kaçar B (2020) Geobiology 18:394 (PMID 32065506, PMC7216921, **OA**). GitHub `kacarlab/AncientNitrogenase` `Extant-MLAnc_Align.fasta` has 319 sequences: extant labelled `Nif_<species>`, `Anf_…`, `Vnf_…`, `BchChl_…`, plus reconstructed ancestors `AncA-1…` **[VERIFIED-DATA]**. Accessions are in the docx supplements (not opened). Raymond J et al. (2004) MBE 21:541 (PMID 14694078, PW; 110 genomes). Boyd ES et al. (2015) J Bacteriol (PMID 25733617, PMC4403663). Prokaryotic HGT-rich, with paralogous Bch/Chl outgroups.

**X6 β-lactamases — P4.** The NCBI AMR Reference Gene Catalog (`ReferenceGeneCatalog.txt`, db_version 2026-08-07.1, md5 20efd3f3…) has 6,128 β-lactamase entries, each with allele, gene family (blaOXA 1,431, blaPDC 702, blaADC 383, blaKPC 295, blaCTX-M 268, blaTEM 233, …), product name with Ambler class (A 2,247; B1 472; B2 26; B3 151; C 1,784; D 1,443), RefSeq/GenBank protein + nucleotide accessions and PubMed references **[VERIFIED-DATA]**. Curated by NCBI with literature references (not phylogeny-based). Naas T et al. (2017) BLDB, J Enzyme Inhib Med Chem 32:917 (PMID 28719998, PMC6445328, **OA**; bldb.eu up). Hall BG, Barlow M (2004) Drug Resist Updat (PMID 15158767, PW). **Benchmark use:** class B (metallo, MBL fold) is non-homologous to A/C/D, so it is a **negative control against false merging**. Families are allele clusters on mobile elements, so orthology vs xenology tests are possible with host genomes (chromosomal intrinsic blaPDC in *Pseudomonas*, blaADC/blaOXA-51-like in *Acinetobacter* are vertical). Requires bacterial genome panel design (strain-level).

**X7 B-family DNA polymerases — P4.** Tahirov TH, Makarova KS, Rogozin IB, Pavlov YI, Koonin EV (2009) Biol Direct 4:11 (PMID 19296856, PMC2669801, **OA**). Add. file 1 = Exo-Pol alignment; Add. file 2 = Newick trees with labels `GI|taxon-abbrev|species` **[VERIFIED-DATA]** (GI numbers are deprecated but still resolvable via NCBI E-utilities). Makarova KS, Krupovic M, Koonin EV (2014) Front Microbiol 5:354 (PMID 25101062, PMC4104785, **OA**). Phenomena: Pol ε contains an inactivated second Exo-Pol module (**fusion**); eukaryotic polymerases have a chimeric archaeal origin. Strong archaeal–eukaryotic test.

**X8 Two-component systems — P3.** Alm E, Huang K, Arkin A (2006) PLoS Comput Biol 2:e143 (PMID 17083272, PMC1630713, **OA**). Dataset S4 = FASTA of all HPK domains (3.0 MB), S1 = per-genome HGT vs LSE counts, S3 = species tree **[LISTED]**. Capra EJ, Laub MT (2012) Annu Rev Microbiol (PMID 22746333, PMC4097194). The HPK subfamily classification (11 classes, per memory; UNVERIFIED) is tool-derived, so use it as a stress test, not a gold standard.

**X9 Chaperonins — P3.** Archibald JM, Logsdon JM, Doolittle WF (1999) Curr Biol 9:1053 (PMID 10508614, PW). Lund PA (2009) FEMS Microbiol Rev 33:785 (PMID 19416363, PW). The eight CCT subunits are LECA paralogs, so this is a near-single-copy-per-paralog test at LECA.

**X10 Rubisco — P2.** Tabita FR et al. (2008) J Exp Bot 59:1515 (PMID 18281717, PW). Banda DM et al. (2020) Nat Plants (PMID 32868887, PW). Schulz L et al. (2022) Science (PMID 36227987, PW). rbcL is plastid-encoded, which is unsuitable for nuclear-proteome panels apart from form II in dinoflagellates and the RLPs.

**X11 CRISPR-Cas — P2.** Makarova KS et al. (2020) Nat Rev Microbiol 18:67 (PMID 31857715, PMC8905525). Makarova KS et al. (2015) (PMID 26411297, PMC5426118). Makarova KS et al. (2025) Nat Microbiol (PMID 41198952, PMC12669027, **OA**). Extreme HGT/modularity. Optional stress test only.

**X12 Asgard ESPs — P3.** Eme L et al. (2023) Nature 618:992 (PMID 37316666, PMC10307638, **OA**; there is a 2026 Author Correction, PMID 41673166). Zaremba-Niedzwiedzka K et al. (2017) Nature 541:353 (PMID 28077874, PW). Use specific families (ESCRT, small GTPases, actin/gelsolin) rather than the ESP concept.

**Tool-assisted references (do NOT use as ground truth; secondary only).** COG (Galperin MY et al. 2025, NAR, PMID 39494517, PMC11701660; 2021, PMID 33167031) and arCOG (Makarova KS, Wolf YI, Koonin EV 2015, Life 5:818, PMID 25764277, PMC4390880, **OA**; FTP live at ftp.ncbi.nlm.nih.gov/pub/wolf/COGs/arCOG/). These are manually curated on top of tool-generated clusters, which helps check prokaryotic membership, but they are flat (non-hierarchical) and partly circular with BLAST-based methods. The same applies to PLAZA, PlantTFDB, iTAK, OrthoDB, eggNOG and 1KP gene-family trees (One Thousand Plant Transcriptomes Initiative 2019, Nature 574:679; PMID 31645766, PMC6872490, **OA**).

**Plant WGD reference events (for event-level assertions, not families).** Bowers JE et al. (2003) Nature 422:433 (α/β/γ; PMID 12660784, PW). Thomas BC, Pedersen B, Freeling M (2006) Genome Res 16:934 (α homeolog fractionation; PMID 16760422, PMC1484460). Jaillon O et al. (2007) Nature 449:463 (γ hexaploidy; PMID 17721507, PW). Jiao Y et al. (2011) Nature 473:97 (ζ seed plant / ε angiosperm; PMID 21478875, PW). Paterson AH, Bowers JE, Chapman BA (2004) PNAS 101:9903 (ρ cereals; PMID 15161969, PMC470771). Wang X et al. (2011) Nat Genet 43:1035 (*B. rapa* triplication; PMID 21873998, PW). Schmutz J et al. (2010) Nature 463:178 (soybean; PMID 20075913, PW). Synteny-derived homeolog lists (e.g. Thomas 2006, BRAD for *B. rapa*) are tool-derived but independent of sequence-similarity orthology. They can supply **event-level** assertions (this paralog pair = α), which are ideal for hierarchical metrics on any family above.

---

## 4. Ranked shortlist (top 15, chosen to maximise diversity)

| Rank | ID | Family | Why it is in the top 15 (diversity axis covered) |
|---|---|---|---|
| 1 | E1 | Kinesins | Eukaryote-wide (45 spp., all supergroups), LECA depth, family/subfamily labels with accessions in one XLS. Domain-sharing merge/split test. |
| 2 | P3 | Plant homeobox (14 classes) | Plant superfamily with accession+class table over red alga→maize. WGD (ρ, maize) plus a shared-domain class boundary. |
| 3 | E2 | ARF/Arl/Sar | 2,286 curated rows, 115 spp., "conclusively annotated" flags. Short fast paralogs and gene-model corrections. |
| 4 | P1 | MADS-box MIKC | The canonical plant WGD family (γ, ρ, α, Brassica ×3) with explicit duplication timing in the literature. Needs PDF curation and a full-text request. |
| 5 | P10 | NLRs (RefPlantNLR + pan-NLRome) | Hardest plant case: promiscuous domains, tandem birth–death, integrated-domain fusions. Functionally validated anchors with GenPept accessions. |
| 6 | X2 | Ribosomal proteins | Prokaryotic near-single-copy negative control (995 genomes, copy-number matrix). Calibrates false splits and merges. |
| 7 | E7 | MCM2-9 | Clean, slow, deep (archaea→LECA) paralog set, near single copy. "Easy-deep" calibrator for hierarchical levels. |
| 8 | E12 | Aminoacyl-tRNA synthetases | LUCA paralogy plus HGT plus EGT/dual targeting (xenologs in one genome). |
| 9 | X1 | RecA/RadA/Rad51/Dmc1 | Three-domain LUCA family, LECA duplication, EGT of plant organellar RecA, fast Rad51 paralogs. |
| 10 | E4 | Tubulin superfamily | 3,286 labelled sequences from 504 spp. Slow core vs patchy δ/ε/ζ losses; recent duplications; prokaryotic BtubA/B and FtsZ outgroups. |
| 11 | P2 | Phytochromes + neochrome | Seed-plant/angiosperm-era duplications, fern/moss/algal independent duplications, a fusion protein transferred by HGT. Dryad data (CC0, manual download). |
| 12 | X6 | β-lactamases (NCBI catalog) | Prokaryotic HGT on plasmids; non-homologous class B as negative merge control; fully machine-readable with accessions. |
| 13 | X7 | B-family DNA polymerases | Archaeal–eukaryotic chimera; Pol ε inactivated-module fusion; Newick trees with GI accessions. |
| 14 | E5 | Myosins | Multi-domain shuffling, 30+ classes, manually annotated (CyMoBase, figshare CC BY). |
| 15 | P16 | Aquaporins/MIPs (+ P6 GID1/CXE as alternate) | Mid-size plant family with HGT origin (NIPs) and labelled Nexus trees. GID1-in-CXE is the alternate for "family nested in a family". |

Coverage check:
* **Depth.** LUCA (X1, X2, E12, X7). LECA (E1, E2, E4, E5, E7). Land plant / angiosperm (P1, P2, P3, P16). Recent: Brassica ×3 and maize (P1, P3), plus allelic β-lactamases (X6).
* **Rate.** Slow (E4 α/β, E7, X2) vs fast (E2, P10, Rad51 paralogs).
* **Phenomena.** WGD (P1, P2, P3). HGT (X6, E12, P2, P16). EGT (E12, X1). Fusion (X7, P2 neochrome, P10). Promiscuous domains (E1, P10, E5). Loss (E4, E1). Single-copy control (X2, E7). Non-homology control (X6 class B).
* **Taxonomy.** Viridiplantae, all eukaryotic supergroups, bacteria and archaea.

Alternates if any top-15 fails curation: E9 (TSET/AP, 874 accessions), E14 (cryptochromes, accession table), P7 (bHLH table), X5 (nitrogenase GitHub), E6 (dynein Newick), P17 (LHC).

---

## 5. Paywalled or non-OA papers worth a human full-text/supplement request (ranked)

| # | Paper | Why |
|---|---|---|
| 1 | Shao Z-Q et al. 2016 Plant Physiol 170:2095, doi:10.1104/pp.15.01487 (PMID 26839128) + supplements | Ancestral-lineage framework for >6,000 angiosperm NLRs (23 ancestral lineages). This is the hierarchy for P10. |
| 2 | Parenicová L et al. 2003 Plant Cell 15:1538, doi:10.1105/tpc.011544 (PMID 12837945) + supplements | Complete Arabidopsis MADS list with phylogenetic groups (P1). |
| 3 | Shiu S-H et al. 2004 Plant Cell 16:1220, doi:10.1105/tpc.020834 (PMID 15105442) + supplements | Arabidopsis + rice RLK subfamily membership lists (P11). |
| 4 | Lehti-Shiu MD et al. 2009 Plant Physiol 150:12, doi:10.1104/pp.108.134353 (PMID 19321712) | RLK/Pelle across land plants (P11). |
| 5 | Brindefalk B et al. 2007 MBE, doi:10.1093/molbev/msl202 (PMID 17182897) | Mitochondrial aaRS origins (E12). |
| 6 | Duchêne A-M et al. 2005 PNAS, doi:10.1073/pnas.0504682102 (PMID 16251277) | Arabidopsis organellar aaRS dual targeting list (E12). |
| 7 | Karlgren A et al. 2011 Plant Physiol, doi:10.1104/pp.111.176206 (PMID 21642442) | PEBP across land plants (P8). |
| 8 | Hernández-García J et al. 2019 MBE, doi:10.1093/molbev/msz009 (PMID 30668817) | DELLA origin (pairs with P6). |
| 9 | Becker A, Theissen G 2003 MPE, doi:10.1016/s1055-7903(03)00207-0 (PMID 14615187) | MADS clade definitions (P1). |
| 10 | Nam J et al. 2003 MBE, doi:10.1093/molbev/msg152 (PMID 12777513) | MADS antiquity and duplication timing (P1). |
| 11 | Muller J et al. 2005 MBoC, doi:10.1091/mbc.e05-06-0508 (PMID 16195354) | ARP classification (E16). |
| 12 | Cobbe N, Heck MM 2004 MBE, doi:10.1093/molbev/msh023 (PMID 14660695) | SMC phylogeny (E8). |
| 13 | Leung KF, Dacks JB, Field MC 2008 Traffic, doi:10.1111/j.1600-0854.2008.00797.x (PMID 18637903) | ESCRT across eukaryotes (E17). |
| 14 | Raymond J et al. 2004 MBE, doi:10.1093/molbev/msh047 (PMID 14694078) | Nitrogenase natural history (X5). |
| 15 | Mulkidjanian AY et al. 2007 Nat Rev Microbiol, doi:10.1038/nrmicro1767 (PMID 17938630) | F/V-ATPase evolution (X4). |
| 16 | Figge RM et al. 1999 MBE, doi:10.1093/oxfordjournals.molbev.a026125 (PMID 10331270) | GAPDH EGT/HGT (P21). |
| 17 | Lurin C et al. 2004 Plant Cell, doi:10.1105/tpc.104.022236 (PMID 15269332) | PPR list (P12). |
| 18 | Chen F et al. 2011 Plant J, doi:10.1111/j.1365-313x.2011.04520.x (PMID 21443633) | TPS subfamilies (P14). |
| 19 | Nelson D, Werck-Reichhart D 2011 Plant J, doi:10.1111/j.1365-313x.2011.04529.x (PMID 21443632) | Plant P450 clans (P18). |
| 20 | Engelken-type LHC: Neilson JA, Durnford DG 2010 Photosynth Res, doi:10.1007/s11120-010-9576-2 (PMID 20596891) | LHC diversification (P17). |
| 21 | Sayou C et al. 2014 Science, doi:10.1126/science.1248229 (PMID 24436181) | LEAFY across land plants (P9). |
| 22 | Hall BG, Barlow M 2004 Drug Resist Updat, doi:10.1016/j.drup.2004.02.003 (PMID 15158767) | Serine β-lactamase phylogeny (X6). |
| 23 | Archibald JM et al. 1999 Curr Biol, doi:10.1016/s0960-9822(99)80457-6 (PMID 10508614) | Archaeal chaperonin paralogy (X9). |
| 24 | Nakano T et al. 2006 Plant Physiol, doi:10.1104/pp.105.073783 (PMID 16407444) | ERF groups (P19). |
| 25 | Abascal F, Irisarri I, Zardoya R 2014 BBA, doi:10.1016/j.bbagen.2013.12.001 (PMID 24355433) | MIP diversity (P16). |

Also needed: the **PNAS SI** of Lin et al. 2006 (RecA/Rad51, PMC-only, so supplements are not served by Europe PMC; try the pnas.org SI link) and the **JCS supplement** of Elias et al. 2012 (Rab). Both are free at the publisher but not reachable through the EPMC API. Dryad bundles `10.5061/dryad.5rs50` (phytochrome) and `10.5061/dryad.ghx3ffbq3` (SMC) are CC0 but need a manual browser download.

---

## 6. Caveats for the design phase

1. **Label provenance varies.** Kinesin, ARF, homeobox and RefPlantNLR are phylogeny- or experiment-curated with accessions. Rab, SNARE and HistoneDB are HMM-classified from curated seeds. MADS, RLK, PEBP and DELLA need PDF curation. Treat each assertion with a `basis` field (phylogeny / HMM / experimental / synteny).
2. **Old proteome versions.** Most eukaryote-wide datasets (2007–2014) use outdated gene models. Mapping to current proteomes must go by sequence (≥95–100 % identity with coverage checks). Unmapped members should be recorded explicitly, not silently dropped.
3. **Transcriptome-based datasets** (phytochrome, auxin, 1KP) include species without genomes. Restrict to genome-backed species, or use the authors' transcripts as the "proteome" (this changes the benchmark's error model).
4. **Organelle-encoded members** (rbcL, atp, some r-proteins, tufA in some algae) are absent from nuclear proteomes. Either include organellar proteomes or exclude those members.
5. **Prokaryotic panels** need strain-level genome selection (β-lactamase, nitrogenase, TCS). Plasmid content varies within species.
6. **Repository brittleness.** Dryad needs manual download, two Fasshauer-lab DBs are offline, and TreeBASE availability is uncertain. Snapshot everything we use into `data/raw` with checksums.

## 7. Files produced

* `literature/B2_families_plants_protists_prokaryotes.md` (this file)
* `literature/raw/B2/epmc_search.py`: Europe PMC citation verifier. `epmc_hits.jsonl` holds the raw verification records.
* `literature/raw/B2/supp_scan.py`: full-text supplementary/data-link scanner. `scan1.txt`–`scan3.txt` hold its output.
* `literature/raw/B2/abstracts_numbers.txt`: numeric claims extracted from abstracts.
* `literature/raw/B2/cited_pmid_metadata.tsv`: Europe PMC journal/volume/issue/pages/OA flags for all 124 cited PMIDs (every volume:page in this file was checked against it).
* `literature/raw/B2/MANIFEST.tsv`: md5/size/source of all downloaded files.
* `literature/raw/B2/supp/PMC*/`: opened supplements (git-ignored, re-downloadable).
* `literature/raw/B2/github/nitrogenase_Extant-MLAnc_Align.fasta`, `literature/raw/B2/amr_ReferenceGeneCatalog.txt`
