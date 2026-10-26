# A. Orthology benchmarks & metrics — survey focused on hierarchical orthology (HOGs / reconciled gene trees)

Compiled 2026-10-26. All references were checked against Europe PMC (title/DOI/PMID) or Crossref/doi.org; see the reference table at the end. Anything stated without a verified source is marked **UNVERIFIED**. Facts marked "(full text)" were read in the open-access full text, not only in the abstract.

---

## 1. Quest for Orthologs (QfO) benchmark service & reference proteomes

**Service.** [Altenhoff-2016] set up a web service that runs standardized tests on predictions made over a shared proteome set. In the first round, 15 methods were run on 20 benchmark variants. The common denominator is **ortholog pairs**. HOGs and groups uploaded as OrthoXML are flattened into pairs (full text). Results are plotted as precision proxy vs recall proxy with a Pareto frontier. The authors deliberately avoid a combined score because the right trade-off depends on the application. Since 2020 the service runs on ELIXIR OpenEBench (Nextflow plus containers) [Altenhoff-2020].

**Reference proteomes.** These are one representative protein per gene, taken from UniProt reference proteomes and released yearly.
- 2011 release: 66 species (40 eukaryotes, 26 prokaryotes), 754,149 proteins [Altenhoff-2016].
- 2019, 2020 and 2022 releases: 78 species (48 eukaryotes, 23 bacteria, 7 archaea), about 1.0–1.4 M sequences including isoforms [Altenhoff-2020; Nevers-2022; Altenhoff-2024].
- A 2025 review cites 81 proteomes [Sarton-Lohéac-2025].
- The species tree is a literature consensus built from 6 tree-of-life sources, with confidence levels L10–L90 per node. Only 87 of 146 splits were congruent across sources [Boeckmann-2015].

**Tests.**

| Test | What it measures | Notes / limitations |
|---|---|---|
| STD → **GSTD** (generalized species-tree discordance) | Samples a chain of orthologs covering 10 species. Builds an MSA (MAFFT) and a least-squares distance tree, then takes the RF distance to the reference species subtree. Precision proxy = mean RF; recall proxy = number of completed samplings (50,000 species-set samples). Branches shorter than 10 Myr are collapsed to limit ILS [Altenhoff-2016; Altenhoff-2020]. | Confounds orthology error with **gene-tree reconstruction error**. Fast-evolving but correct members raise RF, so more complete HOGs are penalized (simulation by Train, shown in [Sarton-Lohéac-2025], Fig. 9). Only sampled pairs count, and in-paralogs are not rewarded [Altenhoff-2020]. |
| **SwissTree** | Pairwise orthology vs. 19 expert-curated reference gene trees spanning all domains [Boeckmann-2011; FastOMA-2025]. | Few families, chosen to be "tractable". Showed almost no precision–recall trade-off, possibly because ambiguous sequences are excluded from curated trees [Altenhoff-2016]. |
| **TreeFam-A** | Pairs from curated TreeFam trees (1,314 trees, metazoan-centred) [Ruan-2008; Sarton-Lohéac-2025]. | Taxonomically narrow. MetaPhOrs was partly trained on it, so the test is circular [Altenhoff-2016]. TreeFam's last release is v9 [Schreiber-2014]; that it is no longer maintained is **UNVERIFIED**. |
| **GO / EC conservation** | Schlicker semantic similarity of experimentally backed annotations; recall = number of pairs [Altenhoff-2016]. | Relies on the ortholog conjecture, which is only weakly true [Altenhoff-2012]. Function-based tests often miss false assignments [Trachana-2014]. Experiments are guided by sequence similarity, which adds circularity [Altenhoff-2016]. |
| **VGNC** (2020/2022) | 18,282 human-curated ortholog sets in 8 mammals, 4 of them in QfO. True positive rate (TPR) and positive predictive value (PPV). A false positive is counted only when both species already have a gene in a different set, which is conservative because there is no negative set [Nevers-2022]. | Covers only recent divergences and is near saturation: PPV > 0.95 and TPR ≥ 0.9 for almost all methods. HGNC staff co-curate it. **We found no separate "HGNC test"**; HGNC appears through HCOP/consensus calls [Altenhoff-2020]. |
| **FAS** (2024) | Feature-architecture similarity (domains, TM regions, disorder) of predicted pairs [Altenhoff-2024]. | OMA HOGs score lowest because they contain many in-paralogs. The test penalizes hierarchical groups by design. |
| OrthoBench in QfO | A group-level test was prototyped in 2022 [Nevers-2022]. | Covers only one taxonomic level. |

**Stated limitations relevant to HOGs.**
- Flattening HOGs into pairs discards nesting. Pairs between in-paralogs under the root may be scored as false positives [Sarton-Lohéac-2025].
- No QfO test scores paralogy, **duplication timing**, or the hierarchical nesting itself [Sarton-Lohéac-2025].
- Coverage and recall are hard to score fairly. "Inparalog coverage" was proposed but never implemented [Altenhoff-2020].
- The tests assume whole-protein orthology, so domain-level orthology is penalized [Altenhoff-2020; Altenhoff-2024].
- 78 proteomes give poor coverage of many clades [Altenhoff-2024].
- Missing predictions count as false negatives, so subset submissions are discouraged [Altenhoff-2016].
- Knowing the species tree gave methods only a "negligible" advantage on STD [Altenhoff-2016].

## 2. Other benchmarks & resources

- **OrthoBench** [Trachana-2011]: 70 manually curated Bilateria-level reference orthogroups (RefOGs) from 12 metazoans, tagged by challenge type (duplication/loss, rate, domain architecture, alignment quality). The paper found genome annotation to be the largest single factor, affecting up to about 30% of performance.
  - **Revisited** in [Emms-2020]: 31 of 70 RefOGs changed (24 extensively). Up-to-date proteomes, the HMMER-to-tree evidence and a commentary per RefOG are on GitHub and Zenodo.
  - Metric: gene-pair precision/recall/F computed within RefOGs. TP, FP and FN are divided by (n−1) per RefOG to avoid bias towards large families.
  - Findings: best methods reached about 77% F; the best method got only 23 of 70 RefOGs exactly right. The main error was **splitting** (low recall) rather than missing genes. Recall fell for large RefOGs, poor MSAs (norMD) and single-domain families; evolutionary rate mattered less (full text).
  - Limitation: a single taxonomic level, which invites overfitting to that level [Sarton-Lohéac-2025].
- **Bacterial RefOGs** ("orthobench2") [Trachana-2014]: phylogeny-curated bacterial groups. The study shows that function-based validation misjudges methods.
- **LEMOrtho** (lemortho.ezlab.org; reported in the OrthoDB v11 paper [Kuznetsov-2023]): an interactive benchmark on the revised RefOGs.
  - Per-RefOG F1, precision and recall, exact matches, "akin" matches, and split events. Variation of Information is mentioned by [Sarton-Lohéac-2025].
  - Split / fusion / "complex mismatch" categories follow [Kriventseva-2015].
  - The benchmark is run by a method developer (OrthoDB); see bias in §5.
- **YGOB-based yeast benchmark** [Salichos-2011]: 2,723 synteny-curated ortholog groups across 6 yeasts, 3 of them post-WGD. **Differential loss after WGD creates "traps"** (the retained genes are paralogs), and the false discovery rate (FDR) of all methods rises sharply there. This is the clearest real-data benchmark of hidden paralogy [Byrne-2005; Scannell-2006].
- **SwissTree**: concept and pilot [Boeckmann-2011]. A 2017 preprint on how taxon sampling affects nodes exists (bioRxiv 10.1101/181966; only on bioRxiv in EPMC).
- **PANTHER** [Thomas-2022], **Ensembl Compara** [Vilella-2009; Herrero-2016], **eggNOG/OrthoDB/OMA** [Altenhoff-2024b OMA; Kuznetsov-2023]: these are tree- or hierarchy-based *resources*, not benchmarks. They often serve as references, which is a circularity risk. Vilella-2009 introduced "two non-sequence-based metrics of gene tree correctness". The exact Compara duplication-consistency (species-overlap) score formula is **UNVERIFIED** here because the Ensembl docs were unreachable.
- **Method papers with their own HOG or group evaluation:**
  - **GETHOGs** [Altenhoff-2013]: ALF simulations with known orthologs and paralogs.
  - **OMA 2.0** [Train-2017]: robustness to asymmetric rates and fragments.
  - **OrthoFinder 3** [Emms-2026]:
    - Evaluated on OrthoBench with 7 measures: % missing RefOGs, % missing genes, % RefOGs fused, % RefOGs split (fission), recall, precision, and **entropy** of RefOG dispersion. Methods are ranked on each measure and the ranks averaged.
    - Reports that FastOMA fused 10% of RefOGs vs about 2.85% for the others.
    - QfO GSTD.
    - **Only root-level orthogroups were benchmarked, not the per-node HOGs** (full text).
  - **FastOMA** [Majidian-2025]:
    - QfO results: SwissTree precision 0.955 at recall 0.69; GSTD Eukaryota normalized RF 0.225.
    - A **reference-free "implied losses" measure** per rootHOG, used to show that a better-resolved species tree (TimeTree vs NCBI) gives more parsimonious histories.
    - Robustness to species-tree errors injected by the authors (Supp. Figs 18–20).
    - **No per-level HOG accuracy test** (full text).
  - **SonicParanoid2** [Cosentino-2024]: domain-aware language model; QfO and empirical sets.
  - **Broccoli** [Derelle-2020]: detects chimeric/fused proteins; "two benchmark data sets".
  - **Possvm** [Grau-Bové-2021]: species-overlap plus MCL on gene trees; high precision and recall vs manually curated orthogroups when given accurate trees.
  - **Hieranoid 2** [Kaduk-2017]; **Domainoid** [Persson-2019].
  - **REvolutionH-tl** [Ramírez-Rafael-2026]: best-match-graph method benchmarked on synthetic data against OrthoFinder, GeneRax and RANGER-DTL, including reconciliation output.
- **Reconciliation / gene-tree accuracy (simulation):**
  - **GeneRax** [Morel-2020]: relative RF to the true tree (best in 90% of simulations).
  - **ALE** [Szöllősi-2013]; **AleRax** [Morel-2024]; **Treerecs** [Comte-2020].
  - **SimPhy** [Mallo-2016]: shows a systematic **overestimation of duplication time** by standard (LCA) reconciliation under ILS.
  - **DLCoal** [Rasmussen-2012]: ILS-aware reconciliation improves orthologs, paralogs, duplications and losses.
  - **Hahn-2007**: gene-tree error biases reconciliation, pushing **duplications towards the root and losses towards the tips**.
- **Simulators:**
  - **ALF** [Dalquen-2012]: genome-level evolution, including LGT, fusion/fission and indels.
  - **ZOMBI** [Davín-2020]: genomes with extinct lineages and transfers.
  - **SaGePhy** [Kundu-2019]: gene and sub-gene (domain) evolution.
  - **SimPhy** [Mallo-2016]: duplication-loss + ILS + HGT + gene conversion.
- **Tooling:**
  - OrthoXML [Schmitt-2011]; **OrthoXML-Tools** [Kharrazi-2025].
  - **pyHam/iHam** [Train-2019]: per-branch duplication, loss and gain from HOGs; usable for event extraction.
  - **OMAmer** [Rossier-2021]: per-level subfamily placement. It shows 18–62% of closest-sequence assignments are wrong, usually *over-specific*. This is a useful precedent for per-level evaluation.
- **Reviews and perspectives:** [Gabaldón-2009; Glover-2019; Linard-2021; Majidian-2025b; Glover-2025]. **The 2025 HOG review [Sarton-Lohéac-2025] states explicitly that "the field still lacks robust benchmarking tailored to HOGs"**. It calls for metrics on consistency across levels, duplication placement and ancestral gene content, and for RefOGs at multiple levels.

## 3. Metrics usable for HOG / reconciled-tree accuracy

| Metric | Definition | Strengths | Weaknesses | Ref |
|---|---|---|---|---|
| Pairwise ortholog precision/recall | Compare the predicted ortholog relation with the reference relation over gene pairs. | Works for any method; QfO standard. | Loses nesting; rewards neither in-paralogs nor timing; needs a negative set (VGNC workaround). | Altenhoff-2016; Nevers-2022 |
| Orthology/paralogy relation accuracy (3-class) | For each reference pair, label it ortholog / in-paralog / out-paralog relative to a level, or by its LCA event; give a confusion matrix. | Uses the HOG structure. The relation graph is equivalent to an event-labelled tree when it forms a cograph. | Many pairs in big families; needs weighting. | Hellmuth-2013; Lafond-2014 |
| Per-level group P/R/F (OrthoBench-style) | At each species-tree node L, compare predicted HOGs at L with reference HOGs at L, using pair counts with (n−1) normalization. | Directly measures "HOG at level L"; extends OrthoBench to many levels. | Needs a level mapping between species trees; reference must state membership at L. | Emms-2020 (single level) |
| Split / fusion / complex mismatch; % exact | Count RefOGs spread over ≥2 predicted groups (fission) or sharing a group with another RefOG (fusion). | Diagnostic and intuitive; separates over-splitting from over-merging. | Counts depend on what non-reference genes are allowed. | Kriventseva-2015; Emms-2026 |
| Entropy / Variation of Information | Conditional entropy of the predicted partition given the reference (and the reverse); VI = H(P\|R) + H(R\|P). | Information-theoretic; separates homogeneity from completeness. | Less intuitive; sensitive to singletons. | Meilă-2003; Emms-2026 |
| V-measure | Harmonic mean of homogeneity and completeness (entropy-based). | Bounded [0,1]. | Biased with many clusters. | Rosenberg-2007 (ACL D07-1043) |
| **B-cubed P/R** | Per-item precision and recall of its cluster, then averaged. | Only metric satisfying all 4 formal constraints in Amigó et al.; per-gene, so it is naturally size-aware. | Treats every item equally; must be restricted to reference genes. | Amigó-2009 |
| ARI | Chance-corrected pair agreement. | Standard way to measure tool–tool agreement. | Dominated by large clusters; hard to interpret per family. | Hubert-1985; used in Emms-2020 |
| Duplication-node recovery & **timing** | For each reference duplication (with its species-tree branch), check whether a predicted duplication separates the same paralogous clades, and on which branch. Report signed offset (too deep / too shallow). | Directly tests the hierarchical content. Captures the known root-ward bias. | Needs literature duplications located on the species tree; ILS confounds it. | Hahn-2007; Mallo-2016 |
| RF / normalized RF | Symmetric difference of bipartitions. | Linear time; well known. | Ignores event labels; sensitive to single misplaced taxa; GSTD conflates it with gene-tree error. | Robinson-1981 |
| Labeled RF (LRF) / ELRF | RF extended with node-label flips and edits (speciation/duplication). LRF is a metric computed in linear time. | Compares reconciled trees directly. Found that denser taxon sampling gives better topology but **worse labelling**. | Needs fully reconciled reference trees; not normalized per level. | Briand-2020; Briand-2022 (pylabeledrf) |
| PLR (path-label reconciliation) | Semi-metric combining topology, the gene→species map and event labels; tunable weight; linear time. | Designed for reconciled gene trees. | Recent (WABI 2024, arXiv 2407.06367); not peer-reviewed in a journal (**not in EPMC**). | López-Sánchez-2024 |
| Triplet / path distances | Fraction of rooted triplets resolved differently; path-length difference. | Triplets are robust to a single rogue leaf; work on rooted trees. | Ignore labels unless extended. | Critchlow-1996; Steel-1993 |
| Reconciliation distances (event/path) | Distances between DTL reconciliations. | Score the event mapping itself. | Designed for alternative reconciliations of one tree. | Grueter-2021 (TCBB, Crossref-verified) |
| Reference-free: implied losses, completeness, HOG counts per level, tool–tool VI | Parsimony and plausibility checks. Mammalian-level HOG counts range from 16,214 to 56,136 across databases, and converge to about 13–16k after an ≥80% completeness filter. | No reference needed; flags over-splitting. | Not accuracy; optimizing them can be gamed (e.g., by over-merging). | Majidian-2025; Sarton-Lohéac-2025; Altenhoff-2024b |

## 4. Documented sources of error (and quantification)

- **Annotation and gene models.**
  - OrthoBench: annotation was the largest single factor, affecting up to about 30% of performance [Trachana-2011].
  - The choice of NCBI, Ensembl, UniProt or Augustus changes HOG completeness and QfO scores [Prieto-Baños-2025]. Isoform choice matters too (same paper).
  - Mixing annotation methods inflates the apparent number of lineage-specific genes [Weisman-2022].
  - Fragments are handled specifically by OMA 2.0 and FastOMA [Train-2017; Majidian-2025].
- **Contamination.**
  - Over 2 M contaminated GenBank entries [Steinegger-2020].
  - OMArk found strong contamination in 73 of 1,805 UniProt eukaryotic proteomes, and error propagation in bird annotations from a fragmented zebra finch reference [Nevers-2025].
  - Case study: the tardigrade "extensive HGT" claim [Boothby-2015] was largely contamination [Koutsovoulos-2016].
  - We found **no dedicated study quantifying how contamination affects orthology/HOG accuracy**. This is a gap our perturbation tests could fill.
- **Evolutionary rate and long branches.** In simulations with no gains or losses, orthology errors grow with evolutionary rate and mimic real gain/loss and phylostratigraphy patterns [Natsidis-2021]. Asymmetric rates after duplication are a further problem [Train-2017].
- **Duplication-rich clades.** Bidirectional best hits (BBH) misses up to 60% of orthologs in plants and animals [Dalquen-Dessimoz-2013].
- **HGT/LGT.** In simulation, LGT disrupts *all* methods. Ambiguous sequences hurt alignment-score-based methods more than distance-based ones [Dalquen-2013].
- **WGD, hidden paralogy and pseudo-orthology.**
  - Differential loss after yeast WGD [Scannell-2006] creates traps that inflate FDR for every method [Salichos-2011].
  - Hidden paralogy in cyclostomes after the vertebrate 2R duplications [Kuraku-2010].
  - Pseudo-orthologs in general [Koonin-2005]; paralogs entering phylogenomic matrices drive artefactual topologies [Siu-Ting-2019].
- **ILS.** Standard reconciliation misreads deep coalescence as duplication plus loss, and overestimates duplication age [Rasmussen-2012; Mallo-2016]. QfO collapses branches shorter than 10 Myr [Altenhoff-2016].
- **Gene-tree and alignment error.**
  - Duplications are pushed to the root and losses to the tips [Hahn-2007].
  - Alignment uncertainty changes downstream inference [Wong-2008].
  - Systematic tree errors [Kapli-2021]; incongruence review [Steenwyk-2023].
- **Species-tree errors.**
  - Tree-of-life sources disagree on 40% of splits [Boeckmann-2015].
  - A better-resolved tree reduces implied losses in FastOMA, which is also robust to injected taxonomy errors [Majidian-2025].
  - Misplaced branches split or merge HOGs [Sarton-Lohéac-2025].
- **Domains and promiscuous domains.**
  - Domain architecture is mostly conserved in orthologs [Forslund-2011].
  - Single-domain RefOGs were *hardest* for recall in OrthoBench [Emms-2020].
  - Domain-aware methods: Domainoid, SonicParanoid2, Broccoli (chimeras).
  - QfO acknowledges partial orthology as unbenchmarked [Altenhoff-2024].

## 5. Simulation-based vs literature-based benchmarks; known biases

- **Simulation** (ALF, ZOMBI, SimPhy, SaGePhy) gives full ground truth, including duplication timing, and isolates individual factors [Dalquen-2013]. Its weaknesses:
  - It is model-bound: sequence models are too clean, and there are no annotation or contamination artefacts unless injected.
  - Methods that share the generative model (DL/DTL likelihood) may be favoured.
- **Literature / curated** benchmarks (OrthoBench, SwissTree, TreeFam-A, YGOB, VGNC) contain real artefacts. Their known biases:
  - **Ascertainment:** they focus on well-studied, "tractable" families and model organisms, and ambiguous sequences are dropped. This removes the precision–recall trade-off [Altenhoff-2016].
  - **Taxonomic narrowness:** mostly metazoan/vertebrate, and a single level [Sarton-Lohéac-2025].
  - **Circularity:** reference trees built with the same tools or databases; parameters tuned on the benchmark (MetaPhOrs on TreeFam-A; InParanoid and PhylomeDB on QfO) [Altenhoff-2016]. FastOMA places sequences using the OMA database via OMAmer [Majidian-2025], so references derived from OMA would leak into it.
  - **Tree-based tests** (GSTD, LRF to the species tree) penalize correct but divergent members and conflate gene-tree error with orthology error [Sarton-Lohéac-2025].
  - **Function-based tests** misjudge accuracy [Trachana-2014].
  - **Self-benchmarking:** OF3 vs FastOMA claims conflict, and each group runs its own benchmark (LEMOrtho by OrthoDB, OrthoBench by the OrthoFinder authors).
  - **Absence of a negative set** makes precision ill-defined [Nevers-2022].
  - **Species-tree dependence:** two estimates are compared with each other [Sarton-Lohéac-2025].

## 6. Practical: OrthoFinder 3 and FastOMA

**OrthoFinder 3.1.5** (installed in `conda_envs/ortho`; GitHub now lists v3.2).
- Input: a directory with one FASTA per species, one representative isoform per gene (`primary_transcript.py` helps).
- Species tree: optional via `-s` (rooted Newick). Otherwise it is inferred with STAG and rooted.
- Defaults in the installed code (`process_args.py`): `-M msa` trees with **FAMSA + FastTree**, DIAMOND search. `--core/--assign` is the scalable SHOOT workflow, recommended above about 64–100 species.
- Outputs (`Results_*/`):
  - `Phylogenetic_Hierarchical_Orthogroups/N<k>.tsv`: one file per species-tree node. Columns: `HOG, OG, Gene Tree Parent Clade, <species…>`, with comma-separated genes per species.
  - Node names map to clades via `Species_Tree/SpeciesTree_rooted_node_labels.txt`.
  - `Orthogroups/Orthogroups.tsv`, `Gene_Duplication_Events/Duplications.tsv` (one row per duplication, with species-tree node and support), `Resolved_Gene_Trees/`, `Orthologues/`.
  - The v3.2 README says N0.tsv was removed and merged into `Orthogroups/Orthogroups.tsv`. The 3.1.5 code still writes N0.tsv and deletes it only with `rm_legacy`. **Verify on a toy run.**
- HOGs exist only for species-tree nodes, so supplying the *same* species tree to both tools is required for comparable levels.

**FastOMA 0.5.0** (bioconda; GitHub latest v0.5.1). It is a Nextflow pipeline (`nextflow run dessimozlab/FastOMA`).
- Inputs:
  - `in_folder/proteome/*.fa` (file name = species name).
  - **Rooted species tree required**. It may be non-binary; leaf names must equal file names; internal node names must be unique.
  - **OMAmer database**: the default `LUCA.h5` is about 7.7 GB and downloads automatically.
  - Optional: `splice/` isoform files and precomputed `hogmap_in/`.
- Outputs:
  - `FastOMA_HOGs.orthoxml`: nested HOGs with taxonomic levels; `paralogGroup` elements mark duplications.
  - `RootHOGs.tsv`, `OrthologousGroups.tsv` (strict, ≤1 gene per species), `orthologs.tsv.gz`, `species_tree_checked.nwk`, `hogmap/`, `RootHOGsFasta/`, `OrthologousGroupsFasta/`, `stats/`, `report.html`, `phylostratigraphy.html`.
  - Parse with pyHam or orthoxml-tools.

**Runtime (verified figures only):**
- FastOMA: 2,086 eukaryotic proteomes in under 24 h on 300 cores [Majidian-2025]. The 3-species test data runs in minutes (README).
- OF3: 1,024 species in 128 h with 32 threads; FastOMA could not finish 1,024 within the 7-day limit. On a 64-proteome set, runtime gains plateau above about 32 threads [Emms-2026].
- **No verified wall-clock figure for 20–80 proteomes on 12 cores**; the OF3 Fig. 4 / Supp. Data 2 numbers were not extracted. Expect hours (**UNVERIFIED**) and time a pilot run.

---

## Implications for our design

1. **Encode references as level-specific, provenance-tagged assertions**, not single flat groups. Each assertion should name its species-tree node: "genes X are one HOG at node L"; "a duplication on branch b separates clades A and B"; "pair (g1, g2) is ortholog/paralog". Unresolved items are recorded and excluded from scoring, because there is no negative set.
2. **Primary metric: per-level HOG P/R/F (pair-based, normalized by n−1 per family) plus B-cubed.** Report it at every reference level and stratified by depth. This directly fills the gap identified by QfO and Sarton-Lohéac-2025.
3. **Report fission/fusion/exact counts and VI/entropy** alongside the primary metric, so that over-splitting (typical of OMA-like methods) and over-merging (reported for FastOMA in OF3) are visible separately.
4. **Duplication recovery and timing metric:** for each literature duplication, report recovered (yes/no) and signed branch offset (too deep / correct / too shallow). Expect a root-ward bias [Hahn-2007; Mallo-2016].
5. **Three-class relation accuracy** (ortholog / in-paralog / out-paralog at the reference level), so in-paralog coverage is rewarded, which QfO cannot do.
6. **Use labelled tree distances (LRF via pylabeledrf; PLR optional) only for families with fully reconciled reference trees.** Keep them secondary.
7. **Do not use GSTD-style RF to the species tree or GO/EC similarity as accuracy measures.** They conflate gene-tree error and the ortholog conjecture.
8. **Score only reference genes.** Extra non-reference genes inside a HOG count as false positives only if the reference asserts family completeness for those species. Document this rule (VGNC-style conservative precision).
9. **Give both tools the identical rooted species tree** (OF `-s`; FastOMA requires one) and map levels by leaf sets. Add a robustness arm with perturbed or collapsed species trees and with the NCBI taxonomy.
10. **Control circularity:**
    - Record how each reference tree was built and exclude references derived from OMA, OrthoFinder, Ensembl Compara or PANTHER outputs, or flag them.
    - Log the OMAmer DB version, and test FastOMA sensitivity to the reference DB (for example, a smaller or clade DB).
11. **Proteome hygiene and phenomenon tags:**
    - One isoform per gene; record the annotation source.
    - Run OMArk/BUSCO to tag contamination and fragmentation.
    - Tag families by WGD, ILS (internodes shorter than about 10 Myr), HGT, fast rate, multi-domain and size, and report metrics per tag.
12. **Perturbation tests as controls:** fragment genes, inject a contaminant proteome, swap the annotation source, drop or add species. Report metric deltas; these quantify robustness, which the literature has not measured for contamination.
13. **Validate metric code on simulations with known truth** (small ALF/ZOMBI/SimPhy runs) and on toy trees: positive control = truth scored against itself; negative controls = random and over-merged partitions.
14. **No single combined score.** Show precision–recall per level with bootstrap CIs over families (resampling families, not genes) and Pareto fronts. Report method–method agreement (ARI/VI) as a difficulty signal, not as accuracy.
15. **Taxonomic and level breadth matters more than family count.** Prefer families with literature-dated duplications at several depths and outside metazoa (plants, fungi, protists, prokaryotes), because existing sets are single-level and animal-centred.

---

## References (verified via Europe PMC unless noted; DOI | PMID)

- Altenhoff-2016 Standardized benchmarking in the QfO. Nat Methods | 10.1038/nmeth.3830 | 27043882
- Altenhoff-2020 QfO benchmark service and consensus calls in 2020. NAR | 10.1093/nar/gkaa308 | 32374845
- Nevers-2022 QfO orthology benchmark service in 2022. NAR | 10.1093/nar/gkac330 | 35552456
- Altenhoff-2024 New developments for the QfO benchmark service. NAR Genom Bioinform | 10.1093/nargab/lqae167 | 39664814
- Altenhoff-2024b OMA orthology in 2024. NAR | 10.1093/nar/gkad1020 | 37962356
- Altenhoff-2009 Phylogenetic and functional assessment of orthologs inference projects. PLoS Comput Biol | 10.1371/journal.pcbi.1000262 | 19148271
- Altenhoff-2012 Resolving the ortholog conjecture. PLoS Comput Biol | 10.1371/journal.pcbi.1002514 | 22615551
- Altenhoff-2013 Inferring HOGs from orthologous gene pairs. PLoS One | 10.1371/journal.pone.0053786 | 23342000
- Altenhoff-2019 OMA standalone. Genome Res | 10.1101/gr.243212.118 | 31235654
- Boeckmann-2011 Benchmark phylogenomic DBs on reference gene trees. Brief Bioinform | 10.1093/bib/bbr034 | 21737420
- Boeckmann-2015 QfO entails quest for tree of life. GBE | 10.1093/gbe/evv121 | 26133389
- Boothby-2015 Tardigrade HGT. PNAS | 10.1073/pnas.1510461112 | 26598659
- Briand-2020 Generalized RF for labeled trees. BMC Genomics | 10.1186/s12864-020-07011-0 | 33208096
- Briand-2022 Linear-time labeled RF. Syst Biol | 10.1093/sysbio/syac028 | 35426933
- Byrne-2005 Yeast Gene Order Browser. Genome Res | 10.1101/gr.3672305 | 16169922
- Chen-2007 Assessing orthology detection strategies. PLoS One | 10.1371/journal.pone.0000383 | 17440619
- Comte-2020 Treerecs. Bioinformatics | 10.1093/bioinformatics/btaa615 | 33085745
- Cosentino-2024 SonicParanoid2. Genome Biol | 10.1186/s13059-024-03298-4 | 39054525
- Dalquen-2012 ALF. Mol Biol Evol | 10.1093/molbev/msr268 | 22160766
- Dalquen-2013 Impact of duplication, indels, LGT, sequencing error: simulation. PLoS One | 10.1371/journal.pone.0056925 | 23451112
- Dalquen-Dessimoz-2013 BBH miss many orthologs. GBE | 10.1093/gbe/evt132 | 24013106
- Davín-2020 Zombi. Bioinformatics | 10.1093/bioinformatics/btz710 | 31566657
- Derelle-2020 Broccoli. Mol Biol Evol | 10.1093/molbev/msaa159 | 32602888
- Emms-2015 OrthoFinder. Genome Biol | 10.1186/s13059-015-0721-2 | 26243257
- Emms-2019 OrthoFinder (v2). Genome Biol | 10.1186/s13059-019-1832-y | 31727128
- Emms-2020 Revisiting Orthobench. GBE | 10.1093/gbe/evaa211 | 33022036
- Emms-2026 OrthoFinder: improved... (v3). Nat Methods | 10.1038/s41592-026-03126-6 | 42265210 (preprint 10.1101/2025.07.15.664860; correction 10.1038/s41592-026-03238-z)
- Forslund-2011 Domain architecture conservation in orthologs. BMC Bioinformatics | 10.1186/1471-2105-12-326 | 21819573
- Gabaldón-2009 Joining forces in the QfO. Genome Biol | 10.1186/gb-2009-10-9-403 | 19785718
- Glover-2019 Advances and applications in the QfO. Mol Biol Evol | 10.1093/molbev/msz150 | 31241141
- Glover-2025 Perspectives on orthology during the QfO. J Mol Evol | 10.1007/s00239-025-10290-4 | 41263969
- Grau-Bové-2021 Possvm. Mol Biol Evol | 10.1093/molbev/msab234 | 34352080
- Hahn-2007 Bias in tree reconciliation methods. Genome Biol | 10.1186/gb-2007-8-7-r141 | 17634151
- Hellmuth-2013 Orthology relations, symbolic ultrametrics, cographs. J Math Biol | 10.1007/s00285-012-0525-x | 22456957
- Herrero-2016 Ensembl comparative genomics resources. Database | 10.1093/database/baw053 | 27141089
- Hulsen-2006 Benchmarking ortholog identification with functional genomics data. Genome Biol | 10.1186/gb-2006-7-4-r31 | 16613613
- Kaduk-2017 Hieranoid 2. Bioinformatics | 10.1093/bioinformatics/btw774 | 28096085
- Kapli-2021 Systematic errors in phylogenetic trees. Curr Biol | 10.1016/j.cub.2020.11.043 | 33497629
- Kharrazi-2025 OrthoXML-Tools. J Mol Evol | 10.1007/s00239-025-10271-7 | 41003686
- Koonin-2005 Orthologs, paralogs, and evolutionary genomics. Annu Rev Genet | 10.1146/annurev.genet.39.073003.114725 | 16285863
- Koutsovoulos-2016 No evidence for extensive HGT in tardigrade. PNAS | 10.1073/pnas.1600338113 | 27035985
- Kriventseva-2015 OrthoDB v8. NAR | 10.1093/nar/gku1220 | 25428351
- Kundu-2019 SaGePhy. Bioinformatics | 10.1093/bioinformatics/btz081 | 30715213
- Kuraku-2010 Hidden paralogy, hagfish/lamprey. Integr Comp Biol | 10.1093/icb/icq044 | 21558193
- Kuznetsov-2023 OrthoDB v11 (LEMOrtho described therein). NAR | 10.1093/nar/gkac998 | 36350662
- Lafond-2014 Orthology and paralogy constraints. BMC Genomics | 10.1186/1471-2164-15-S6-S12 | 25572629
- Linard-2021 Ten years of the QfO. Mol Biol Evol | 10.1093/molbev/msab098 | 33822172
- Majidian-2025 FastOMA. Nat Methods | 10.1038/s41592-024-02552-8 | 39753922
- Majidian-2025b QfO in the era of data deluge and AI. J Mol Evol | 10.1007/s00239-025-10272-6 | 41085653
- Mallo-2016 SimPhy. Syst Biol | 10.1093/sysbio/syv082 | 26526427
- Menet-2022 Phylogenetic reconciliation (review). PLoS Comput Biol | 10.1371/journal.pcbi.1010621 | 36327227
- Morel-2020 GeneRax. Mol Biol Evol | 10.1093/molbev/msaa141 | 32502238
- Morel-2024 AleRax. Bioinformatics | 10.1093/bioinformatics/btae162 | 38514421
- Natsidis-2021 Systematic errors in orthology inference. iScience | 10.1016/j.isci.2021.102110 | 33659875
- Nevers-2025 OMArk. Nat Biotechnol | 10.1038/s41587-024-02147-w | 38383603
- Persson-2019 Domainoid. BMC Bioinformatics | 10.1186/s12859-019-3137-2 | 31660857
- Prieto-Baños-2025 Annotation matters. Bioinformatics | 10.1093/bioinformatics/btaf365 | 40560601
- Ramírez-Rafael-2026 REvolutionH-tl 2.0. PLoS Comput Biol | 10.1371/journal.pcbi.1013017 | 42585261
- Rasmussen-2012 DLCoal. Genome Res | 10.1101/gr.123901.111 | 22271778
- Rossier-2021 OMAmer. Bioinformatics | 10.1093/bioinformatics/btab219 | 33787851
- Ruan-2008 TreeFam 2008. NAR | 10.1093/nar/gkm1005 | 18056084
- Salichos-2011 Ortholog prediction in a yeast model clade. PLoS One | 10.1371/journal.pone.0018755 | 21533202
- Sarton-Lohéac-2025 Reconstructing evolutionary histories with HOGs. J Mol Evol | 10.1007/s00239-025-10277-1 | 41269329
- Scannell-2006 Reciprocal gene loss in polyploid yeasts. Nature | 10.1038/nature04562 | 16541074
- Schmitt-2011 SeqXML and OrthoXML. Brief Bioinform | 10.1093/bib/bbr025 | 21666252
- Schreiber-2014 TreeFam v9. NAR | 10.1093/nar/gkt1055 | 24194607
- Siu-Ting-2019 Inadvertent paralog inclusion. Mol Biol Evol | 10.1093/molbev/msz067 | 30903171
- Steenwyk-2023 Incongruence in the phylogenomics era. Nat Rev Genet | 10.1038/s41576-023-00620-x | 37369847
- Steinegger-2020 Terminating contamination. Genome Biol | 10.1186/s13059-020-02023-1 | 32398145
- Szöllősi-2013 ALE. Syst Biol | 10.1093/sysbio/syt054 | 23925510
- Thomas-2022 PANTHER. Protein Sci | 10.1002/pro.4218 | 34717010
- Trachana-2011 OrthoBench. BioEssays | 10.1002/bies.201100062 | 21853451
- Trachana-2014 Phylogeny-based benchmark, limits of function-based validation. PLoS One | 10.1371/journal.pone.0111122 | 25369365
- Train-2017 OMA 2.0. Bioinformatics | 10.1093/bioinformatics/btx229 | 28881964
- Train-2019 iHam & pyHam. Bioinformatics | 10.1093/bioinformatics/bty994 | 30508066
- Vilella-2009 EnsemblCompara GeneTrees. Genome Res | 10.1101/gr.073585.107 | 19029536
- Weisman-2022 Mixing annotation methods. Curr Biol | 10.1016/j.cub.2022.04.085 | 35588743
- Wong-2008 Alignment uncertainty. Science | 10.1126/science.1151532 | 18218900
- Crossref/doi.org-verified (not in PubMed): Robinson-1981 Math Biosci 10.1016/0025-5564(81)90043-2; Critchlow-1996 Syst Biol 10.1093/sysbio/45.3.323; Steel-1993 Syst Biol 10.1093/sysbio/42.2.126; Hubert-1985 J Classif 10.1007/BF01908075; Amigó-2009 Inf Retr 10.1007/s10791-008-9066-8; Meilă-2003 LNCS 10.1007/978-3-540-45167-9_14; Grueter-2021 IEEE/ACM TCBB 10.1109/tcbb.2019.2942015; Rosenberg-2007 V-measure, EMNLP-CoNLL, ACL Anthology D07-1043 (URL verified); López-Sánchez-2024 PLR, arXiv 10.48550/arXiv.2407.06367 (WABI 2024).
- Web resources checked 2026-10-26: github.com/DessimozLab/FastOMA (README), github.com/OrthoFinder/OrthoFinder (README), lemortho.ezlab.org (cited in OrthoDB v11 text; site itself not fetched).
