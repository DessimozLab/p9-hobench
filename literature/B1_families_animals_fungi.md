# B1 — Candidate reference gene families: animals (incl. early-branching Metazoa & unicellular holozoans) and fungi

Survey date: 2026-10-26 (session start 13:17Z; this report written ~13:50Z). Author: B1 literature agent.
Scope: families whose evolutionary history (membership, subfamily structure, timing of duplications/losses, HGT, WGD fate) is resolved in publications, assessed for use as literature-based reference sets for **hierarchical** orthology benchmarking.

## 0. How this was produced (reproducibility)

* Every citation below was checked against Europe PMC REST by DOI with `scripts/B1_verify_refs.py`
  → `literature/raw/B1/refs_verified.tsv` (131 references, **131/131 status OK**; columns: PMID, PMCID, first author, year, journal, `isOpenAccess`, `inPMC`, title). Input DOI list: `literature/raw/B1/dois_round1.txt`. Abstracts of all 131 are cached in `literature/raw/B1/abstracts.txt` (counts quoted from abstracts below come from there).
  13 DOIs I first wrote from memory were wrong and were corrected by title search before inclusion (lesson: never trust remembered DOIs).
* OA flags: "OA" = Europe PMC `isOpenAccess=Y`; "PMC" = `inPMC=Y` but not in the OA subset (full text readable, bulk supplement download not possible by API); "closed" = neither (candidate for human full-text request).
* Supplements of 50+ papers were fetched with the Europe PMC `.../{PMCID}/supplementaryFiles` endpoint (works only for the OA subset); non-image files were unpacked to `literature/raw/B1/supp/<tag>__<PMCID>/` (266 MB — **should be git-ignored or LFS'd; not committed by me**). Direct PMC "bin/" downloads for non-OA papers are blocked by a JavaScript/recaptcha challenge — they must be requested from a human.
* Databases were probed live (2026-10-26). Several flagship literature-curated resources are **currently offline**; Wayback Machine copies were used where possible (see §1).
* Nothing here is ground truth yet: "data verified" means I opened the file and confirmed it contains identifiers + labels, not that labels were cross-checked against current proteomes.

## 1. Infrastructure check: curated databases (live status today)

| Resource | Content relevant to us | Status 2026-10-26 | What I could get |
|---|---|---|---|
| **YGOB** (Byrne & Wolfe 2005, PMID 16169922) – yeast WGD pillars | manually curated, synteny-aware homology "pillars" for 20 Saccharomycetaceae + reconstructed ancestor; ohnolog pairs on 2 tracks | **ygob.ucd.ie refuses connections (port 80, "Connection refused")**; same host serves CGOB → also down | Wayback copy of `Pillars.tab` (v7 Aug2012) is **truncated at 1 MiB** (original 2,778,211 B): 5,173 pillars, 297 *S. cerevisiae* ohnolog pairs, 33 columns verified. Saved as `raw/B1/ygob/Pillars.TRUNCATED_wayback_20210728.tab` + README. **Full file needs human/alternate network or e-mail to K. Wolfe lab.** |
| **KinBase** (Manning et al.) – kinomes with Group/Family/Subfamily | kinase.com port 80 connects but returns nothing; https refused | Wayback: `Kincat_Hsap.08.02.xls` (620 human kinase entries: 516 non-pseudogene, 10 groups, 133 families; RefSeq/Entrez IDs) and `kinase_subfamily_evolution.06.10.txt` (294 subfamilies × 7 species copy numbers + "Predicted History", e.g. Metazoan 69, Primordial 39, Unikont-lost-in-yeast 12). Totals: human 518, urchin 353, fly 239, worm 443, yeast 131, *Dictyostelium* 295, *Tetrahymena* 1063. Saved in `raw/B1/kinbase/`. |
| **HomeoDB / HomeoDB2** (Zhong & Holland 2008, 2011) | homeobox genes by class/family for ~10 animals | **homeodb.zoo.ox.ac.uk unreachable**; Wayback has per-species HTML family pages only, FASTA export was form-driven (not archived) | fallback: Holland et al. 2007 (OA) human homeodomain FASTA (S7) + HGNC "Homeoboxes" group (342 genes; HOXL subclass 52). |
| **CyMoBase** (Kollmar) – myosins, kinesins, dyneins, tubulins, actins | 27,253 manually annotated sequences, 1,529 species | up (https://www.cymobase.org) | interactive search only; no bulk download link found → use paper supplements (below) or ask authors. |
| **OHNOLOGS v2** (Singh & Isambert 2020) | synteny-based 2R/3R ohnolog families, 27 vertebrates | up | downloaded human 2R (2,142 families) and zebrafish 3R (2,033 families), criterion set "[A]" (exact stringency label to be checked in the help page) → `raw/B1/ohnologs/`. |
| **HGNC gene groups** | expert-curated human gene-family membership with hierarchy | up | `family.csv`, `gene_has_family.csv`, `hierarchy.csv` downloaded (2,569 groups) → `raw/B1/hgnc/`. Counts include pseudogenes (e.g. OR 873, CYP 113). |
| Xenbase GenePageReports | *X. laevis* .L/.S homeolog naming, *X. tropicalis* orthology | up (download.xenbase.org directory listing OK; web UI behind recaptcha) | not downloaded yet. |
| Nelson Cytochrome P450 homepage | manually curated CYP nomenclature (family/subfamily = hierarchy) | root up (drnelson.uthsc.edu 200) | not downloaded. |
| HistoneDB 2.0 (NCBI) | curated histone variant classification | up | not downloaded. |
| RabDB | Rab subfamily classification | up (rabdb.org 200) | not downloaded. |
| FCPD fungal P450 (Moktali 2012) | 22,940 CYPs, 213 species, clans | p450.riceblast.snu.ac.kr redirects (302) – status unclear | — |

**Implication for design:** YGOB, KinBase and HomeoDB — the three most "literature-curated" databases on the wish list — are not reliably reachable. The benchmark must snapshot whatever we obtain with checksums (DATA.md) and should not depend on live servers.

## 2. Summary table of candidates (43)

Tags: WGD2R / WGD3R / WGD4R / WGDyeast / ALLO (allopolyploid) · TAND (tandem arrays) · LSE (lineage-specific expansion) · HGT · ILS · FAST (fast evolution / long branches) · PROM (promiscuous / multi-domain) · FUS (fusion/fission) · DEEP (pre-metazoan / LECA-age events) · LOSS (differential loss / hidden paralogy) · CTRL (simple/near 1:1 control).
"Data" column: ✔ = machine-readable IDs + labels verified by me; ◐ = labels exist but in PDF / not yet retrieved / needs mapping; ✘ = text/figures only.

| # | Family / resource | Scope | Tags | Key refs (OA status) | Data | Prio |
|---|---|---|---|---|---|---|
| F1 | Yeast WGD ohnologs & orthologs (YGOB pillars) | Saccharomycetaceae (20 spp. + ancestor) | WGDyeast, LOSS (pseudo-orthology), ALLO | Byrne&Wolfe 2005 (PMC); Scannell 2006 (closed); Scannell 2007 (PMC); Marcet-Houben&Gabaldón 2015 (OA) | ◐ (truncated archive ✔; MH2015 365 ohnolog pairs + topology ✔) | **5** |
| F2 | GAL pathway (GAL1/GAL3 ohnologs; GAL cluster origins incl. HGT to *Schizosaccharomyces*; parallel losses) | Fungi (80 spp. in Slot&Rokas) | WGDyeast, HGT, LOSS, cluster | Hittinger 2004 (PMC); Hittinger&Carroll 2007 (closed); Slot&Rokas 2010 (PMC) | ✘/◐ | 3 |
| F3 | Fungal HGT set: URA1 DHOD (bacterial→Saccharomycetaceae), W/S-clade bacterial genes (ADH reacquisition, sucrose enzyme), enterobactin operon | Saccharomycotina | HGT | Hall 2005 (PMC); Gonçalves 2018 (OA); Kominek 2019 (PMC); Shen 2018 (PMC); Marcet-Houben&Gabaldón 2010 (closed); Fitzpatrick 2012 (closed) | ◐ | **4** |
| F4 | Sterigmatocystin cluster HGT *Aspergillus*→*Podospora* | Pezizomycotina | HGT, cluster (23 genes) | Slot&Rokas 2011 (closed) | ✘ | 3 |
| F5 | Subtelomeric MAL/IMA/SUC families | *Saccharomyces* & relatives | LSE, TAND-like, FAST | Brown 2010 (PMC; 4 xls in PMC, not scriptable) | ◐ | 3 |
| F6 | Fungal kinome | 30 fungi | PROM | Kosti 2010 (OA; RAR supplement) | ◐ | 2 |
| F7 | Fungal CYPome (FCPD) | 213 fungi/oomycetes; 1→>300 CYPs/genome | LSE, PROM-like (huge superfamily) | Moktali 2012 (OA) | ◐ (DB status unclear) | 3 |
| F8 | Gene fusions/fissions across Fungi | 115 fungal genomes; 63 fusions, 110 fissions | FUS, HGT | Leonard&Richards 2012 (PMC; SD01.xls not scriptable); Durrens 2008 (OA) | ◐ | **4** (phenomenon) |
| F9 | CTG-clade (*Candida*) families; CGOB | *Candida*/*Lodderomyces* | LSE, gene-finding errors | Butler 2009 (OA); Maguire 2013 (OA, Table S1 ✔) | ◐ (CGOB down) | 2 |
| F10 | Microsporidia as fast-evolving overlay | Microsporidia | FAST, LOSS (reduced genomes) | Nakjang 2013 (OA; Tables S1–S8 ✔); Inagaki 2004 (closed) | ◐ | **4** (as overlay) |
| F11 | Fungal hydrophobins | Fungi | FAST, short proteins | Lovett 2022 (bioRxiv preprint only) | ✘ | 1 |
| C1 | **Kinesin superfamily** | 45 eukaryotes incl. animals, fungi, choanoflagellate | DEEP, PROM (motor + diverse tails), LOSS | Wickstead 2010 (OA); Lawrence 2004 nomenclature (OA) | **✔** 1,624 seqs, GenBank acc + Family + Subfamily (S3.XLS); full alignments (S8/S9) | **5** |
| C2 | **Myosin superfamily** | pan-eukaryote (2,269 seqs/328 spp.; 7,852/929 spp.) | DEEP, PROM, LSE, LOSS | Odronitz&Kollmar 2007 (OA); Sebé-Pedrós 2014 (OA); Kollmar&Mühlhausen 2017 (OA) | **✔** class-encoded FASTA (2,267 names, S5) + alignment; 725-seq FASTA (SP2014) | **5** |
| C3 | Dynein heavy/intermediate/light chains | 24 eukaryotes; Kollmar 2016 broader | DEEP, PROM (huge multi-domain DHC), LOSS | Wickstead&Gull 2007 (OA); Kollmar 2016 (OA) | **✔** per-subunit sheets with RefSeq IDs (158 DHC rows); Newick trees (Kollmar 2016 Data1/2) | 4 |
| C4 | Tubulin superfamily (α–ζ) | 504 species, 3,524 tubulins | DEEP, recent dups, CTRL-ish core | Findeisen 2014 (OA) | ✔ alignment of 3,286 named seqs; species table; subfamily labels only partly encoded in names (mapping needed) | 4 |
| C5 | Rab GTPases | pan-eukaryote; LECA ≥20 Rabs; human 66 | DEEP, PROM (small-GTPase domain → merging with Ras/Rho/Arf), LSE | Klöpper 2012 (OA); Elias 2012 (PMC); Diekmann 2011 RabDB (OA) | ✔ subfamily-labelled alignments (578 seqs in S4 alone) | 4 |
| C6 | MCM2–9 | all eukaryotes | DEEP (7 pre-LECA duplications), CTRL (no Mcm2-7 loss) | Liu 2009 (OA) | ✔ Newick trees (S5/S7) | **4** (control) |
| C7 | SMC proteins | eukaryotes | DEEP, CTRL | Cobbe&Heck 2004 (closed) | ✘ | 3 |
| C8 | **Protein kinome** (ePK groups/families/subfamilies) | human, fly, worm, yeast, urchin, *Dictyostelium*, *Tetrahymena*, *Monosiga* | PROM, DEEP, LSE (nematode, ciliate), LOSS | Manning 2002 Science (closed); Manning 2002 TIBS (closed); Manning 2008 Monosiga (PMC); Pincus 2008 (PMC); Suga 2012 (closed) | **✔** (Wayback KinBase tables) | **5** |
| C9 | SH2-domain proteins (decoy for domain-based merging) | human/mouse; premetazoan pTyr toolkit | PROM | Liu 2006 Mol Cell (closed); Pincus 2008 (PMC); HGNC "SH2 domain containing" (101 genes) | ◐ | 3 |
| A1 | **Homeobox incl. Hox/ParaHox** | Metazoa (ctenophore, sponge, cnidarian → vertebrates) | TAND, WGD2R/3R, DEEP, LOSS, short domain | Holland 2007 (OA); Ryan 2006 (OA); Ryan 2010 (EvoDevo; EPMC flags not OA); Amores 1998 (closed); Brooke 1998 (closed); Kuraku&Meyer 2009 (closed); HomeoDB2 (OA) | ◐ (HomeoDB down; human FASTA ✔; HGNC ✔) | **5** |
| A2 | Nuclear receptors | Bilateria; *C. elegans* >270 NRs (HNF4 burst) | LSE (nematode), WGD2R, FAST (supnrs), LOSS | Bertrand 2004 (closed); NR nomenclature 1999 (closed); Robinson-Rechavi 2005 (closed) | ◐ (HGNC 49 human) | 4 |
| A3 | bHLH | 12–13 metazoans incl. sponge, cnidarians (695 bHLHs) | DEEP, LOSS | Simionato 2007 (OA); Ledent 2001 (PMC) | ✔ family × species count tables (S1…S18 xls) | 4 |
| A4 | Sox | chordates; 11 teleosts | WGD3R (11/19 sox duplicated in teleosts), WGD2R | Voldoire 2017 (OA); Bowles 2000 (closed); Heenan 2016 (closed) | ✔ accession table per subgroup × species (s006.xls) | 4 |
| A5 | Fox | Metazoa | DEEP, clusters | Shimeld 2010 (closed); Larroux 2008 (closed) | ✘ | 3 |
| A6 | T-box | Holozoa incl. *Capsaspora*, fungi(!) | DEEP, LOSS | Sebé-Pedrós 2013 (PMC) | ◐ | 3 |
| A7 | Wnt | Metazoa (13 subfamilies; protostome losses) | DEEP, LOSS, TAND (wnt1-6-10 cluster) | Kusserow 2005 (closed); Cho 2010 (PMC); Janssen 2010 (OA) | ✔ FASTA alignments (S3/S4) | 4 |
| A8 | TGF-β pathway (ligands, receptors, Smads) | Metazoa + *Trichoplax*, sponge | DEEP | Huminiecki 2009 (OA) | ✔ Smad FASTA + MrBayes trees (S5.zip) | 3 |
| A9 | **Opsins** | Metazoa (14 phyla), vertebrates, teleosts | DEEP, WGD2R, WGD3R, LSE (teleost RH1/RH2), TAND, gene conversion | Ramirez 2016 (OA); Feuda 2012 (PMC); Lagman 2013 (OA); Musilova 2019 (PMC); Musilova 2023 (closed); Davies 2012 (closed) | **✔** 768 opsins with GenBank IDs+names (Ramirez T1); Ensembl IDs + subtype per species (Lagman S1) | **5** |
| A10 | Neuropeptide GPCRs | Bilateria | DEEP, LOSS, PROM (7TM) | Mirabeau&Joly 2013 (PMC); Jékely 2013 (PMC) | ◐ | 3 |
| V1 | **Globins** | vertebrates (+ amphioxus synteny); placental β-cluster | WGD2R, TAND, LOSS (embryonic genes), convergent cooption (lamprey Hb) | Hoffmann 2012 (PMC); Hoffmann 2010 (PMC); Opazo 2008 (PMC); Storz 2011 (PMC) | ◐ (Table S1 PDF in PMC; 1006756107 sd01.txt exists) | **5** |
| V2 | OHNOLOGS v2 genome-wide 2R/3R families | 27 vertebrates | WGD2R, WGD3R, LOSS | Singh 2015 (OA); Singh&Isambert 2020 (OA) | **✔** downloaded | 4 (overlay, computational) |
| V3 | Teleost 3R ohnolog pairs anchored on gar | zebrafish vs spotted gar | WGD3R, LOSS | Braasch 2016 (OA) | **✔** ST23: 1,606 zebrafish and 1,315 medaka TGD ohnolog pairs + 10,415 / 9,265 singletons, all anchored to gar gene IDs | 4 |
| V4 | **Salmonid Ss4R ohnologs: AORe vs LORe** | salmonids + pike outgroup | WGD4R; ohnolog divergence *after* speciation (LORe) | Lien 2016 (OA); Berthelot 2014 (OA); Robertson 2017 (OA); Gundappa 2022 (OA) | **✔** 394 genes with AORe/LORe + per-gene NEXUS tree + alignment (Robertson MOESM2); genome-wide AORe/LORe gene lists (MOESM4); Figshare 10.6084/m9.figshare.14724684 (121,864 trees) | **5** |
| V5 | *Xenopus laevis* L/S homeologs | *X. laevis* vs *X. tropicalis* | ALLO (subgenomes diverged ~34 Ma, merged ~17–18 Ma; >56% retained) | Session 2016 (OA) | ◐ (Xenbase .L/.S naming; Nature SI not in EPMC zip) | 4 |
| V6 | Aquaporins | Metazoa→vertebrates; zebrafish 18 | WGD2R, WGD3R, TAND clusters (AQP2/5/6) | Finn 2014 (OA); Tingaud-Sequeira 2010 (OA) | ◐ (accessions in PDFs) | 4 |
| V7 | Voltage-gated Na channels | choanoflagellate→vertebrates | PROM (4×6TM), WGD2R/3R, LOSS | Zakon 2012 (PMC); Widmark 2011 (closed) | ✘ | 3 |
| V8 | Keratins | vertebrates | TAND, gene conversion, LSE | Vandebergh&Bossuyt 2012 (closed) | ✘ (HGNC 64 human) | 3 |
| V9 | Mammalian olfactory receptors (781 OGGs, 13 placentals) | Placentalia | LSE (elephant ~2000), TAND, pseudogenes, LOSS | Niimura 2014 (OA); Niimura&Nei 2007 (OA) | ◐ supplement on genome.cshlp.org (login redirect for scripts) | 4 |
| V10 | **KRAB-ZNFs** | tetrapods; human ~378–423 | PROM (C2H2 arrays), LSE, TAND, FAST, young (primate-specific) | de Tribolet-Hardy 2023 (OA); Imbeault 2017 (closed); Huntley 2006 (PMC); Emerson&Thomas 2009 (OA) | **✔** 467 human KZFP loci incl. pseudogenes with `age_MA`, cluster, classification (S1); ZNF FASTA for 12+ species (Emerson s009) | **4** |
| V11 | **APOBEC3 (Z1/Z2/Z3 domains)** | Mammals | FUS (domain fusions), LSE, FAST, LOSS | Münk 2012 (OA); LaRue 2009 nomenclature (PMC) | ✔ species × Z-domain gene table (S1.xls, 165 rows) + S7 txt | **4** |
| V12 | TLRs | vertebrates | FAST, LOSS, LSE (teleost) | Roach 2005 (closed) | ✘ | 3 |
| V13 | Chemokines & receptors | 18 vertebrates (553 genes) | FAST, TAND, WGD2R/3R | Nomiyama 2013 (OA) | ✔ ID-to-systematic-name tables (534 + 234 rows, SD8–SD11) | 4 |
| V14 | Interferons | vertebrates | FAST, intronless/retroposed, TAND | Secombes&Zou 2017 (OA, review) | ✘ | 2 |
| V15 | MHC | vertebrates | FAST, polymorphism, assembly gaps | Kelley 2005 (closed) | ✘ | 1 |
| V16 | Ependymins | Metazoa; teleost-specific groups | FAST, WGD3R | Suárez-Castillo 2007 (OA) | ✔ xls tables of new sequences | 3 |
| I1 | **Drosophila P450s (12 genomes)** | *Drosophila* | LSE, TAND, LOSS (D. sechellia), close species (ILS) | Good 2014 (OA); Pollard 2006 ILS (OA) | **✔** Newick of all CYPs labelled `Cyp…_SPECIES`, all-CYP FASTA, gain/loss-on-tree table | **5** |
| I2 | Arthropod P450 | Arthropoda | LSE | Dermauw 2020 (closed); Nelson 2004 nomenclature (closed) | ✘ | 3 |
| I3 | Insect chemoreceptors / OBPs / CSPs | *Drosophila*; Arthropoda | LSE, FAST, PROM (7TM) | McBride 2007 (PMC); Vieira&Rozas 2011 (OA) | ◐ | 3 |
| I4 | Drosophila immune genes | 12 *Drosophila* | FAST, LSE | Sackton 2007 (closed) | ✘ | 2 |
| I5 | Aphid carotenoid genes (fungus→aphid) / whitefly BtPMaT1 (plant→insect) | Hemiptera | HGT (+duplication) | Moran&Jarvik 2010 (closed); Xia 2021 (closed) | ✘ | 3 |
| N1 | **Caenorhabditis chemoreceptors (19 families)** | *C. elegans/briggsae/remanei* (+2 spp. counts) | LSE, TAND, pseudogenes, FAST | Thomas&Robertson 2008 (OA) | **✔** 4 FASTA files, 6,769 seqs, headers `ce|cosmid|gene|locus`; per-family counts, strict-ortholog fractions (S1/S2) | **4** |
| N2 | **Nematode HGT: plant cell-wall-degrading enzymes, GH5 cellulases** | Tylenchida, *Pristionchus*/diplogastrids | HGT + subsequent duplication/loss | Danchin 2010 (PMC); Mayer 2011 (OA); Haegeman 2011 (closed) | ◐ | **4** |
| N3 | Nematode MATH-BTB / F-box expansions | *Caenorhabditis* | LSE, FAST | Thomas 2006 (PMC) | ✘ | 2 |
| R1 | Bdelloid rotifer HGT | *Adineta*, *Rotaria* | HGT, degenerate tetraploidy (homeologs) | Gladyshev 2008 (closed); Boschetti 2012 (OA); Flot 2013 (closed); Eyres 2015 (OA) | ◐ per-transcript alien-index tables (Boschetti s007/s008) | 3–4 |
| D1 | TF repertoires across eukaryotes / holozoans | Holozoa + other eukaryotes | DEEP | de Mendoza 2013 (PMC); Larroux 2008 (closed); Suga 2013 Capsaspora (OA) | ◐ | 3 |
| X1 | Histones (replication-dependent arrays + variants) | eukaryotes | TAND (concerted evolution), CTRL-extreme conservation, variants | Draizen 2016 HistoneDB 2.0 (OA) | ◐ (DB up) | 3 |

(Actins/ARPs, ABC transporters (Dean & Annilo 2005, closed), connexins/innexins (Abascal & Zardoya 2013, closed), TRIMs etc. considered but not pursued; no OA data source found in the time box.)

## 3. Detailed assessments (priority ≥4; brief for others)

### F1. Yeast WGD: YGOB pillars (+ Scannell 2006/2007, Marcet-Houben & Gabaldón 2015) — PRIORITY 5
* **Phenomena:** single allo-/auto-polyploidisation event (~100 Ma) shared by *Saccharomyces*, *Nakaseomyces* (*C. glabrata*), *Naumovozyma*, *Kazachstania*, *Tetrapisispora*, *Vanderwaltozyma*; rapid reciprocal loss of ohnologs → **pseudo-orthology / hidden paralogy**: "4–7% of single-copy genes compared between any two species are not orthologues" (Scannell 2006 abstract). Marcet-Houben & Gabaldón 2015 add evidence for **interspecies hybridisation** at the WGD (allopolyploid), giving ohnolog-pair topologies inconsistent with a simple autopolyploid model (their Fig 1C table: 365 *S. cerevisiae* ohnolog pairs with topology class A/B/C, divergence age).
* **Data:** YGOB `Pillars.tab` (33 columns: 12 post-WGD species × 2 tracks + ancestor + 8 non-WGD species; manually curated using synteny). Byrne & Wolfe 2005 report 551 *S. cerevisiae*, 599 *N. castellii*, 404 *C. glabrata* ohnolog pairs. Full file currently not reachable (§1). MH2015 `s001.xls` (Fig1C pairs, Fig4 topology support) and `s019.xls` (4,739 rows of pre-WGD species gene-order table) downloaded ✔.
* **Hierarchical value:** the gold-standard case for hierarchical orthology — every pillar defines (i) an ancestral gene at the pre-WGD node, (ii) the WGD duplication node, (iii) per-species track assignment that resolves which single-copy genes are true orthologs vs. pseudo-orthologs. Unique in giving *synteny-adjudicated* labels independent of sequence trees.
* **Size:** genome-wide (~5,000 ancestral loci); per family 1–2 genes per post-WGD species.
* **Species needed:** *S. cerevisiae*, *S. kudriavzevii/uvarum*, *C. glabrata*, *N. castellii*, *K. naganishii*, *T. phaffii*, *V. polyspora* (post-WGD) + *Z. rouxii*, *T. delbrueckii*, *K. lactis*, *E. gossypii*, *L. thermotolerans*, *L. kluyveri* (non-WGD). Old YGOB annotations ≠ current UniProt proteomes → ID mapping step needed (gene names like KLLA0F00484g, ADR419C, NCAS0A07650 are GenBank locus tags → mappable via UniProt `gene_orf` / ENA).
* **Expected difficulty:** high for graph/score-based methods (pseudo-orthologs look like 1:1 BBH); moderate for tree-based. Ancient (≈100 Myr) WGD + fast post-WGD loss.
* **Action:** obtain full Pillars.tab (human request / Wolfe lab mirror; Wolfe lab now at UCD). Until then the truncated file (5,173 pillars) is usable as a *subset* if pillars are sampled without bias (they are ordered by ancestral position, so truncation biases to some ancestral chromosomes — note).

### C1. Kinesin superfamily (Wickstead, Gull & Richards 2010) — PRIORITY 5
* 1,624 kinesins from 45 eukaryotes (incl. human, *Drosophila*, *C. elegans*, *S. cerevisiae*, *S. pombe*, *Monosiga*…), 51 anciently derived paralogs, 105 domain architectures. **S3.XLS verified:** columns Analysis ID, Genome project ID, GenBank accession, **Family (Kinesin-1…-20, orphan), Subfamily (e.g. 8A, 3B)**, exclusion flags, Pfam/CDD architecture. Alignments S8 (DatasetS1, 7.7 MB) / S9.
* **Phenomena:** DEEP (51 anciently derived paralogs inferred, mapped onto LCEA), extensive secondary loss (fungi lack several families), PROM (motor domain shared; tails differ — kinesin-3 with FHA/PX/PH), fast-evolving orphans.
* Size: human ~45 (HGNC 46), yeast 6 (approx., from general knowledge — verify from S3).
* Difficulty: high — methods may merge all kinesins into one family (shared motor domain) or split divergent orphans; correct hierarchical labels at LECA/opisthokont/metazoan levels are available.
* Species: any subset of the 45; recommend human, mouse, zebrafish, *Drosophila*, *C. elegans*, *Nematostella*, *Monosiga*, *S. cerevisiae*, *S. pombe*, *Dictyostelium* (outgroup).

### C2. Myosin superfamily (Odronitz & Kollmar 2007; Sebé-Pedrós et al. 2014; Kollmar & Mühlhausen 2017) — PRIORITY 5
* 2,269 manually annotated motor domains / 328 species, 35 classes (2007); 7,852 myosins / 929 eukaryotes, 79 classes (2017). Sebé-Pedrós 2014 focuses on holozoans/fungi (class structure at the animal–fungal origin).
* **Verified files:** Odronitz S5 (FASTA, 2,267 names encoding species + class, e.g. `Ajc_aMyo1_fl`), S6 (alignment); Sebé-Pedrós DataS1.fasta (725 seqs, names `Species_class_index`); Kollmar 2017 xls is species list only — class assignments for 7,852 sequences are in CyMoBase (no bulk export found) and the 31 MB PDF.
* Phenomena: PROM (motor + IQ + MyTH4-FERM/SH3/PH tails), DEEP, LSE (class-1, class-5 expansions), losses; mammalian myosin-II heavy-chain tandem cluster (MYH1/2/3/4/8/13) gives TAND.
* Difficulty: high (domain-sharing + long branches). Excellent for split/merge metrics.

### C8. Protein kinome (KinBase; Manning et al. 2002, 2008) — PRIORITY 5
* Hierarchical classification (Group → Family → Subfamily) is literally a nested labelling. **Verified:** human catalog (620 entries incl. 102 pseudogenes; RefSeq + Entrez IDs) and subfamily-evolution matrix for 7 species with "Predicted History" (origin node of each subfamily: Primordial, Unikont, Opisthokont, Metazoan, Coelomate, Deuterostome, Chordate, Nematode-specific …, "lost in Yeast").
* Size: human 518, worm 443, fly 239, yeast 131 (from the matrix).
* Phenomena: PROM (kinase domain is the archetypal promiscuous domain), LSE (nematode CK1/FER etc., ciliates), DEEP, LOSS.
* Difficulty: high; the subfamily-origin labels allow node-level (hierarchical) scoring directly.
* Caveat: KinBase IDs are 2002–2008 gene models; mapping to current proteomes required. *Monosiga* kinome (Manning 2008, PMC) adds the choanoflagellate tyrosine-kinase expansion (pTyr toolkit, Pincus 2008).

### A1. Homeobox superclass incl. Hox/ParaHox — PRIORITY 5
* Holland, Booth & Bruce 2007 (OA): classification of all human homeobox genes into 11 classes/102 families; S7 = homeodomain FASTA. HGNC groups (Homeoboxes 342 genes, HOXL 52). Ryan 2006 (OA): 130 *Nematostella* homeodomains placed against 228 human / 97 *Drosophila* (S5 NEXUS alignment ✔). Ryan 2010: ctenophore complement. Teleost Hox cluster duplication/loss: Amores 1998, Kuraku & Meyer 2009 (closed).
* Phenomena: TAND (Hox clusters), WGD2R (4 clusters HoxA–D), WGD3R (7 zebrafish clusters with differential loss), DEEP (ParaHox/Hox origin), FAST (*Drosophila* Hox3/zen/bcd), short 60-aa domain → low phylogenetic signal → methods merge/split classes arbitrarily.
* Size: human 300 homeobox loci = 235 probable functional genes + 65 pseudogenes, 11 classes / 102 families (Holland 2007 abstract); HGNC "Homeoboxes" group 342 members; 39 HOX genes in human (general knowledge).
* Data gap: HomeoDB offline → membership for non-human species must be rebuilt from papers (Ryan 2006 alignments; Holland lab papers). Hox clusters alone (≈40 genes/vertebrate) are a clean, well-labelled sub-benchmark (gene symbols HOXA1…HOXD13 are themselves orthology labels across vertebrates; teleost a/b suffix encodes 3R).

### V1. Globin superfamily — PRIORITY 5
* Hoffmann, Opazo & Storz 2012: Hb/Mb/Cygb as 2R paralogs, with synteny to amphioxus; Hoffmann 2010: lamprey vs gnathostome hemoglobins evolved convergently from different ancestors (cooption) — a classic case where *functional names mislead orthology*; Opazo 2008: placental β-cluster embryonic gene differential loss (ε, γ, η…); Storz 2011 review.
* Data: Table S1 (PDF) in PMC; Hoffmann 2010 `sd01.txt` listed in PMC (not scriptable) → **human request for supplements**. Family size ≈ 10–15 per vertebrate (HGNC "Hemoglobin subunits" 13).
* Phenomena: WGD2R, TAND (α- and β-clusters), LOSS, convergent naming. Small and tractable → ideal for exact node-level scoring.

### A9. Opsins — PRIORITY 5
* Ramirez 2016 (OA): reconciled opsin phylogeny across 14 phyla, ≥9 bilaterian ancestral paralogs, ≥4 in cnidarian–bilaterian ancestor. **T1 verified:** 768 rows (GenBank/UniProt ID, organism, opsin name, aligned sequence).
* Lagman 2013 (OA): vertebrate visual opsins + GNAT/GNAI/OT-VP receptors/CACNA1-L paralogon from 2R — **S1 verified**: per-species tables with approved symbols, subtype (LWS/SWS1/RH2…), Ensembl gene/transcript IDs (GRCh37 era).
* Musilova 2019 (PMC): up to 38 RH1 in *Diretmus*; Musilova 2023 (closed): RH2 0–8 copies in 132 teleosts, ≥4 ancestral duplications + gene conversion.
* Phenomena: DEEP, 2R, 3R, LSE, gene conversion (homogenises paralogs → apparent recent duplications), 7TM GPCR → risk of merging with melatonin receptors / other GPCRs (Feuda 2012).

### V4. Salmonid Ss4R: AORe vs LORe ohnologs — PRIORITY 5
* Robertson 2017 (OA): ~25% of each salmonid genome (≥4,550 ohnologs) evolved under **lineage-specific ohnologue resolution (LORe)**: rediploidisation happened *after* speciation, so ohnolog pairs in different salmonid lineages are not 1:1 across species — the duplication node is "smeared" across the species tree. **MOESM2 verified:** 394 genes with ohnolog locations, pike ortholog, mode (AORe/LORe), "independent rediploidization" call, per-gene **NEXUS tree** and alignment. MOESM4: genome-wide gene lists by mode with `hom_id` (homeolog-pair IDs) and NCBI protein IDs (XP_…).
* Gundappa 2022 (OA): 121,864 ohnolog trees across salmonids incl. huchen; Figshare 10.6084/m9.figshare.14724684 (files listed via API ✔; alignments/trees ~500 MB).
* Lien 2016 (OA) Atlantic salmon genome; Berthelot 2014 (OA) rainbow trout.
* **Hierarchical value:** a direct test of whether a method places a duplication at the Ss4R node (AORe regions) *vs.* recognises lineage-specific resolution (LORe) — no other dataset tests this. Pike (*Esox lucius*) = pre-Ss4R outgroup; also contains teleost 3R.
* Species: *Salmo salar*, *Oncorhynchus mykiss*/*kisutch*/*tshawytscha*, *Salvelinus*, *Coregonus*/*Thymallus*, *Hucho hucho*, *Esox lucius*, + non-salmonid teleost.
* Difficulty: very high (young, highly similar ohnologs; LORe ambiguity; collapse of ohnolog regions in assemblies — Gundappa 2022 flags collapsed regions in *T. thymallus*).

### I1. Drosophila P450s across 12 genomes (Good et al. 2014) — PRIORITY 5
* 114 gains, 74 losses mapped on the 12-species phylogeny; 74–94 CYPs per species; ≥77 ancestral. **Verified:** `Supp data1.nwk` — one tree with every leaf named `Cyp<name>_<SPECIES>` (e.g. `Cyp311a1_ERE`), i.e. a **literature-labelled gene tree**; `Supp data2_allCYPAA.fasta`; `Supp data4 gain and loss on tree.xlsx`.
* Phenomena: LSE/TAND (Cyp6/Cyp4 clusters), LOSS (*D. sechellia* elevated loss), interparalog chimeras (FUS-like), very close species (*mel–sim–sec*) where ILS is documented genome-wide (Pollard 2006, OA).
* Insects + recent → complements vertebrate/yeast WGD cases. CYP nomenclature (Nelson) gives an independent family/subfamily hierarchy.
* Caveat: FlyBase gene models of 2007 assemblies; map via FlyBase IDs/CYP names.

### N1. Caenorhabditis chemoreceptor families (Thomas & Robertson 2008) — PRIORITY 4
* ~1,300 genes + ~400 pseudogenes in 19 families in *C. elegans*; *C. briggsae*, *C. remanei* covered; per-family "fraction strict orthologs" and "fraction paralogs" reported. **Verified:** S33–S36 FASTA (6,769 sequences; headers carry WormBase cosmid ID + gene name + locus) and S1/S2 count tables.
* Phenomena: massive LSE in a nematode, TAND, high pseudogene rate, low identity (~20–26% mean pairwise within families) → hard family delimitation; 7TM GPCR merging. Size per genome: hundreds per family.
* Use: subsample 2–3 families (e.g. *srh*, *str*, *srd*) to keep runtime manageable.

### V10. KRAB zinc-finger proteins — PRIORITY 4
* de Tribolet-Hardy 2023 (OA): census of human KZFPs (378 in human per abstract); **S1 verified**: 467 loci incl. pseudogenes with classification, `age_MA` (estimated age), cluster ID — gives a node-of-origin label per human gene. Imbeault 2017 (closed): phylogenetic dating of KZFPs, root reassigned to coelacanth–tetrapod ancestor. Huntley 2006 (PMC): 423 human KRAB-ZNF loci, 103 conserved in mammals. Emerson & Thomas 2009 (OA): ZNF protein/CDS FASTA for ≥12 species (human, mouse, cow, dog, chicken, fugu, *Ciona*, fly, *Anopheles*, *C. elegans* …).
* Phenomena: PROM (C2H2 zinc-finger arrays → over-merging), LSE, TAND clusters (chr19), very young primate-specific genes (ILS among great apes possible — Scally 2012), gene conversion, fast evolution of DNA-contact residues.
* Difficulty: extreme; recommend a primate/mammal panel only.

### V11. APOBEC3 — PRIORITY 4 (fusion/fission representative)
* Münk, Willemsen & Bravo 2012 (OA): A3Z1–A3Z2–A3Z3 arrangement in the placental ancestor; independent duplications/fusions in primates, felids, horses; human 7 A3 genes, some double-domain (fusion of Z-domains). LaRue 2009 (PMC): nomenclature by Z-domain type (labels are hierarchical).
* **Verified:** S1.xls (species × Z1a/Z2a/…/Z3, AICDA, APOBEC1 table, 165 rows), S5.xls, S7.txt.
* Phenomena: FUS (single- vs double-domain genes → which domain is orthologous?), FAST (positive selection), LSE, LOSS. Small family (HGNC "Apolipoprotein B mRNA editing enzyme catalytic subunits": 12 human genes) → tractable.

### C6. MCM2–9 (+ Microsporidia overlay, F10) — PRIORITY 4 (positive control)
* Liu, Richards & Aves 2009 (OA): Mcm2-9 formed by 7 duplications before LECA; **Mcm2-7 present in all genomes studied (no loss)**; Mcm8/9 lost repeatedly. Newick trees ✔.
* Use as an "easy" control (1:1 within each paralog across animals and fungi; ancient duplication nodes well separated). Adding Microsporidia (Nakjang 2013: dramatic reduction, fast evolution) and other long-branch taxa creates a controlled *fast-evolution* stress test on an otherwise easy family. Same overlay applies to tubulins (C4) and kinesins (C1).

### F3/N2/R1/I5. HGT into fungi and animals — PRIORITY 4 (F3, N2), 3–4 (R1), 3 (I5)
* **F3 (fungi):** URA1 dihydroorotate dehydrogenase acquired from Lactobacillales by Saccharomycetaceae after divergence from the *C. albicans* lineage (Hall 2005); W/S clade (*Wickerhamiella/Starmerella*) lost alcoholic fermentation and reinstated it partly through HGT of bacterial genes, plus a bacterial sucrose-assimilation enzyme (Gonçalves 2018, OA, data on Figshare); horizontal operon transfer of an enterobactin-type siderophore pathway from *E. coli* relatives (Kominek 2019, PMC); GAL cluster of *Schizosaccharomyces* acquired from a *Candida* yeast (Slot & Rokas 2010). Shen 2018 (332 yeasts, PMC) quantifies HGT as rare in budding yeasts — useful to calibrate.
* **N2 (nematodes):** plant-parasitic nematode cell-wall-degrading enzymes (GH5 cellulases, pectate lyases, etc.) acquired by several independent transfers from different bacteria, then duplicated (Danchin 2010, PMC); *Pristionchus* cellulases (donor lineage to be extracted from the paper — UNVERIFIED here) with gene turnover across 7 *Pristionchus* + 3 diplogastrids (Mayer 2011, OA; supplement is PDF/DOC only).
* **R1 (bdelloids):** foreign genes concentrated in telomeric regions (Gladyshev 2008); per-transcript alien index tables (Boschetti 2012 s007/s008 ✔), four *Rotaria* species (Eyres 2015). Note *A. vaga* is degenerate tetraploid (Flot 2013) → HGT + homeolog combination.
* **I5:** aphid carotenoid desaturases/cyclases of fungal origin, duplicated after transfer (Moran & Jarvik 2010); whitefly BtPMaT1 from plants (Xia 2021) — both closed.
* **Benchmark logic:** reference assertion = "gene X in recipient clade is orthologous to donor-lineage genes, not to recipient paralogs; acquisition node = …". Methods restricted to eukaryote-only species panels will place these as singletons or as deep paralogs — the metric should test that HGT genes are *not* merged into unrelated native families, and that post-HGT duplications are recovered. Needs donor taxa (bacteria/fungi/plants) in the panel.

### Other priority-4 entries (brief)
* **F8 fusion/fission in Fungi** — 63 fusions, 110 fissions across 115 genomes (Leonard & Richards 2012); supplementary SD01.xls exists in PMC but needs human download. Unique fungal FUS resource.
* **F10 Microsporidia** — see C6; Nakjang 2013 Tables S1–S8 ✔ (gene family gains/losses in microsporidia).
* **C3 dyneins** — Wickstead & Gull 2007: 9 DHC families (2 cytoplasmic, 7 axonemal), 6 IC families; xls with RefSeq IDs per subunit ✔; Kollmar 2016 trees ✔. Loss of all dyneins in some lineages (higher plants, *Entamoeba*; Wickstead & Gull 2007) and of axonemal families in non-ciliated fungi (general expectation — verify in the xls) → LOSS test.
* **C4 tubulins** — 3,524 tubulins/504 species; α/β recent lineage-specific duplications → "many in-paralogs vs. few out-paralogs" test; subfamily labels require parsing (names like `PfTub1A`, `Pf_aTub1A`).
* **C5 Rab GTPases** — ≥20 LECA Rabs, six supergroups (Klöpper 2012); 66 in human; labelled alignments ✔ (`AuAn_Rab8B_2` style). Ras-superfamily domain → merging with other small GTPases.
* **A2 nuclear receptors** — ~25 in Urbilateria (Bertrand 2004), >270 in *C. elegans* via HNF4 burst (Robinson-Rechavi 2005); NR nomenclature (NR1A1…) is hierarchical. All key papers closed → human request.
* **A3 bHLH** — 695 bHLHs in 12 metazoans; family × species count tables ✔ (no accession list seen yet in S1; check S2–S18).
* **A4 Sox** — 11/19 sox genes duplicated in teleosts, all but *sox17* via 3R (Voldoire 2017); s006.xls ✔ = per-subgroup accession table (elephant shark, mouse, human, teleosts) = ready-made orthology table.
* **A7 Wnt** — 13 subfamilies; losses in protostomes; Janssen 2010 alignments ✔.
* **V2/V3 ohnolog overlays** — OHNOLOGS (computational, synteny-based) and gar-anchored zebrafish/medaka 3R pairs (Braasch 2016 ST23 ✔) are *not* literature-curated family histories but independent, synteny-derived assertions of WGD paralogy; useful as genome-wide "duplication-at-2R/3R-node" checks. Must be flagged as a different evidence class from curated families.
* **V5 Xenopus laevis** — L/S homeologs (>56% retained); Xenbase symbols with `.L/.S` encode the duplication node; pair with *X. tropicalis*. Easy to extract, important for ALLO coverage.
* **V6 aquaporins** — 17 vertebrate classes (Aqp0–16), zebrafish 18 paralogs with "dual paralogy"; accessions in PDF supplements → parse or request.
* **V9 mammalian ORs** — 781 OGGs over >10,000 intact ORs in 13 placentals, only 3 OGGs with 1:1 orthology across all 13 (Niimura 2014). Supplement retrieval blocked for scripts → human request. Very high value for LSE; large.
* **V13 chemokines** — 553 genes in 18 vertebrates mapped to systematic names (Nomiyama 2013) ✔ — fast-evolving small proteins where sequence-only methods fail; classification is synteny-based.

### Priority ≤3 (one-liners)
F2 GAL (excellent biology but small; text-level extraction) · F4 STC cluster HGT (closed) · F5 subtelomeric MAL (PMC xls not scriptable) · F6 fungal kinome (RAR, dated) · F7 fungal CYPs (DB status unclear) · F9 CTG/CGOB (CGOB down) · F11 hydrophobins (preprint only) · C7 SMC (closed, small) · C9 SH2 (closed; use as decoy set from HGNC) · A5 Fox, A6 T-box, A8 TGF-β, A10 neuropeptide GPCRs (deep but extraction-heavy) · V7 Nav (closed key paper) · V8 keratins (closed) · V12 TLR (closed) · V14 IFN (review only) · V15 MHC (avoid: polymorphism + assembly) · V16 ependymins (small, FAST, OA tables ✔ — good optional fast-evolving case) · I2–I4 · N3 · D1 · X1 histones (useful negative control for over-merging of variants; DB up).

## 4. Ranked shortlist (top 15) — chosen jointly for diversity

| Rank | Family | Taxon / depth | Main phenomena | Data readiness | Why in |
|---|---|---|---|---|---|
| 1 | Yeast WGD pillars (YGOB + MH2015) | Fungi, ~100 Ma | WGDyeast, pseudo-orthology, ALLO | ◐ (full file needed) | Only synteny-curated, genome-wide hidden-paralogy truth set |
| 2 | Salmonid Ss4R AORe/LORe | Teleosts, 80–95 Ma | WGD4R (+3R), duplication after speciation | ✔ | Unique test of duplication-node placement under delayed rediploidisation |
| 3 | Kinesins | Pan-eukaryote (LECA) | DEEP, PROM, LOSS | ✔ | 1,624 accessions with family/subfamily labels — best-structured deep set |
| 4 | Myosins | Pan-eukaryote | PROM, DEEP, LSE, TAND (MYH) | ✔ | Domain-architecture complexity, curated classes |
| 5 | Protein kinome (KinBase) | Eukaryotes, 7 spp. | PROM, LSE (nematode), DEEP | ✔ | 3-level labels + subfamily-origin node per subfamily |
| 6 | Homeobox / Hox–ParaHox | Metazoa | TAND, 2R/3R, short domain | ◐ | Canonical tandem-array + WGD case; deep-animal coverage |
| 7 | Globins | Vertebrates (+ lamprey) | 2R, TAND, LOSS, convergent naming | ◐ | Small, exactly resolved, synteny-backed |
| 8 | Opsins | Metazoa → teleosts | DEEP, 2R, 3R, LSE, gene conversion | ✔ | Spans cnidarian→teleost-specific depths in one family |
| 9 | Drosophila P450s | 12 *Drosophila* | LSE, TAND, LOSS, ILS-prone | ✔ (labelled Newick) | Insect, recent, gain/loss mapped per branch |
| 10 | Caenorhabditis chemoreceptors | 3–5 *Caenorhabditis* | extreme LSE, pseudogenes, low identity | ✔ | Nematode + family-delimitation stress |
| 11 | KRAB-ZNFs | Tetrapods/primates | PROM (C2H2), LSE, young genes | ✔ | Worst-case merging; per-gene age labels |
| 12 | APOBEC3 | Mammals | FUS (domain fusion), FAST | ✔ | Fusion/fission with explicit domain-level labels |
| 13 | Nematode HGT cell-wall enzymes | Tylenchida, *Pristionchus* | HGT + post-HGT duplication | ◐ | HGT into animals with resolved donor phylogeny |
| 14 | Yeast HGT set (URA1, W/S clade, enterobactin, GAL in *S. pombe*) | Saccharomycotina | HGT | ◐ | Fungal HGT, several independent events |
| 15 | MCM2–9 (+ Microsporidia overlay) | All eukaryotes | DEEP, CTRL; FAST overlay | ✔ | Positive control; controlled long-branch stress |

**Alternates (swap in for coverage or if data fail):** mammalian ORs (V9, needs human-fetched supplement), *Xenopus laevis* L/S homeologs (V5 — recommended as a species-panel overlay on ranks 6–8 rather than a separate family), Sox teleost 3R table (A4), chemokines (V13), Rab GTPases (C5), tubulins (C4), bdelloid HGT (R1), nuclear receptors (A2), fungal fusions/fissions (F8).

**Coverage check of the shortlist:** taxa — fungi/yeasts (1, 3–5, 14, 15), vertebrates (2, 6–8, 11, 12), insects (9), nematodes (10, 13, 5), deep metazoa/holozoa (3–6, 8), Microsporidia (15). Depth — LECA (3, 4, 15), pre-metazoan (5), metazoan/bilaterian (6, 8), 2R ~500 Ma (6–8), yeast WGD ~100 Ma, Ss4R 80–95 Ma, *Drosophila* <60 Ma (9), primate-young KZFPs (11). Rate — slow (15, 3), fast (11, 12, 10). Phenomena — all requested tags covered; ILS only indirectly (rank 9 species set; great apes for 11) — **no family in the literature is resolved specifically for ILS**; recommend ILS be exercised through species-panel choice (Drosophila *mel* subgroup, great apes; see Pollard 2006, Scally 2012, Rasmussen & Kellis 2012 DLCoal) rather than a dedicated family.

## 5. Paywalled / not retrievable by script → candidates for REQUESTED.md

Full text closed (not OA, not in PMC), ordered by importance to the shortlist:
1. Scannell et al. 2006 Nature, 10.1038/nature04562 (yeast reciprocal ohnolog loss) — **plus full YGOB Pillars.tab (server down)**
2. Manning et al. 2002 Science, 10.1126/science.1075762 (human kinome; supplementary tables)
3. Manning et al. 2002 TIBS, 10.1016/S0968-0004(02)02179-5 (kinomes yeast→man)
4. Amores et al. 1998 Science, 10.1126/science.282.5394.1711 (zebrafish Hox clusters)
5. Kuraku & Meyer 2009 Int J Dev Biol, 10.1387/ijdb.072533km (teleost Hox)
6. Imbeault et al. 2017 Nature, 10.1038/nature21683 (KZFP phylogeny/ages)
7. Bertrand et al. 2004 MBE, 10.1093/molbev/msh200 (NR evolution)
8. Robinson-Rechavi et al. 2005 J Mol Evol, 10.1007/s00239-004-0175-8 (nematode HNF4 burst)
9. Gladyshev et al. 2008 Science, 10.1126/science.1156407; Flot et al. 2013 Nature, 10.1038/nature12326 (bdelloids)
10. Moran & Jarvik 2010 Science, 10.1126/science.1187113; Xia et al. 2021 Cell, 10.1016/j.cell.2021.02.014 (insect HGT)
11. Slot & Rokas 2011 Curr Biol, 10.1016/j.cub.2010.12.020; Marcet-Houben & Gabaldón 2010 TIG, 10.1016/j.tig.2009.11.007; Fitzpatrick 2012 FEMS Microbiol Lett, 10.1111/j.1574-6968.2011.02465.x (fungal HGT)
12. Hittinger & Carroll 2007 Nature, 10.1038/nature06151 (GAL1/GAL3)
13. Musilova et al. 2023 Vision Res, 10.1016/j.visres.2023.108204 (teleost RH2)
14. Kusserow et al. 2005 Nature, 10.1038/nature03158 (Wnt in *Nematostella*)
15. Larroux et al. 2008 MBE, 10.1093/molbev/msn047; Shimeld et al. 2010 Genomics, 10.1016/j.ygeno.2009.08.002 (TF classes, Fox)
16. Brooke et al. 1998 Nature, 10.1038/31933 (ParaHox); Zhong et al. 2008 Evol Dev, 10.1111/j.1525-142X.2008.00266.x (HomeoDB)
17. Vandebergh & Bossuyt 2012 MBE, 10.1093/molbev/msr269 (keratins); Widmark et al. 2011 MBE, 10.1093/molbev/msq257 (Nav)
18. Dermauw et al. 2020 IBMB, 10.1016/j.ibmb.2020.103490; Nelson et al. 2004 Pharmacogenetics, 10.1097/00008571-200401000-00001 (P450)
19. Sackton et al. 2007 Nat Genet, 10.1038/ng.2007.60 (Drosophila immunity)
20. Roach et al. 2005 PNAS, 10.1073/pnas.0502272102 (TLRs); Suga et al. 2012 Sci Signal, 10.1126/scisignal.2002733; Liu et al. 2006 Mol Cell, 10.1016/j.molcel.2006.06.001 (SH2)
21. Cobbe & Heck 2004 MBE, 10.1093/molbev/msh023 (SMC); Bowles et al. 2000 Dev Biol, 10.1006/dbio.2000.9883; Heenan et al. 2016 Gene, 10.1016/j.gene.2015.09.013 (Sox)
22. Haegeman et al. 2011 MPMI, 10.1094/MPMI-03-11-0055; Davies et al. 2012 Mol Ecol, 10.1111/j.1365-294X.2012.05617.x; NR nomenclature 1999 Cell, 10.1016/S0092-8674(00)80726-6

Supplements of PMC (non-OA-subset) or publisher-hosted papers that scripts cannot fetch (JS challenge / login redirect) — need a human browser download:
* Niimura et al. 2014 Genome Res, 10.1101/gr.169532.113 — OGG membership tables (genome.cshlp.org suppl)
* Hoffmann et al. 2012 MBE — `supp_msr207_Supplementary_Table1.pdf` (PMC3245541); Hoffmann 2010 PNAS `1006756107_sd01.txt` (PMC2922537)
* Brown et al. 2010 Curr Biol — 4 xls (PMC2877759); Leonard & Richards 2012 PNAS — `1210909110_sd01.xls` (PMC3535628)
* Kominek et al. 2019 Cell — Supplemental Tables 1–4 (PMC7295392); Shen et al. 2018 Cell — supplementary xlsx 8–13 (PMC6291210)
* Session et al. 2016 Nature — Supplementary Tables (homeolog lists); Sebé-Pedrós 2013 PNAS SI; de Mendoza 2013 PNAS SI; Danchin 2010 PNAS SI; Elias 2012 JCS SI
* Kollmar & Mühlhausen 2017 — per-sequence class assignments (CyMoBase export; contact author)

## 6. Notes for DESIGN.md (risks and recommendations)

1. **Evidence classes differ** and should be stored as such: (a) curated phylogeny with labelled trees (kinesins, myosins, Drosophila CYPs, Robertson 2017 trees); (b) curated nomenclature hierarchies (KinBase, HGNC groups, CYP/NR/kinesin/APOBEC3 names, Xenbase .L/.S); (c) synteny-curated orthology/paralogy (YGOB, Lagman 2013, Nomiyama 2013, gar ST23); (d) computational synteny (OHNOLOGS). Metrics should report per class.
2. **Identifier drift** is the main practical cost: most curated sets use 2005–2014 gene models (KinBase, YGOB, Wickstead, Thomas & Robertson, Good 2014). Budget a mapping step (UniProt ID mapping / sequence-identity mapping to current proteomes at ≥95% identity with coverage threshold), and record unmapped fraction as a data-quality metric — never silently drop.
3. **Avoid circularity:** Lagman 2013 lists Ensembl Family IDs and some later resources (OHNOLOGS uses Ensembl gene sets/paralogs) — usable as membership, but not as orthology truth where they rely on Compara trees.
4. **Pseudogenes:** KinBase (102), KZFP S1, Thomas & Robertson (~400 in *C. elegans*), HGNC OR counts include pseudogenes — filter to protein-coding before scoring.
5. **Long-branch overlay design:** add Microsporidia (e.g. *Encephalitozoon cuniculi*) and a fast-evolving animal (e.g. *Oikopleura* or a parasitic nematode) to the panels of conserved families (MCM, tubulins, kinesins) to create controlled FAST cases with known labels.
6. **Server fragility:** snapshot every external file with sha256 (already done for downloaded ones; see `sha256sum` of `raw/B1/{kinbase,ygob,ohnologs,hgnc}`) — KinBase/YGOB/HomeoDB are down now and only partially archived.
7. `literature/raw/B1/supp/` (266 MB) should be git-ignored or tracked via a manifest + checksums rather than committed.

## 7. Files produced

* `literature/B1_families_animals_fungi.md` (this file)
* `scripts/B1_verify_refs.py` — DOI → Europe PMC verification
* `literature/raw/B1/refs_verified.tsv`, `dois_round1.txt`, `abstracts.txt`, `pmc_supp_list.txt` (`refs_round1.tsv` = first pass incl. the 13 wrong DOIs, kept for audit)
* `literature/raw/B1/ygob/` (README + truncated Pillars.tab), `kinbase/`, `ohnologs/`, `hgnc/`, `gundappa2022_figshare/`
* `literature/raw/B1/supp/<tag>__<PMCID>/` — unpacked OA supplements for ~50 papers (non-image files only)
