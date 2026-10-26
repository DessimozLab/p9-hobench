# LITERATURE — survey synthesis

Detailed, citation-verified surveys (every DOI/PMID checked against Europe PMC; unverifiable items flagged UNVERIFIED):
- `literature/A_benchmarks_metrics.md` — existing orthology benchmarks, HOG/reconciliation metrics, documented error sources, tool I/O.
- `literature/B1_families_animals_fungi.md` — 43 candidate families (Metazoa, Fungi), data readiness, ranked top-15.
- `literature/B2_families_plants_protists_prokaryotes.md` — 42 candidates (plants, LECA-level, prokaryotes), ranked top-15.
- `literature/C_genomes_species.md` — proteome sources, licences, QfO 2026_02 set, 5 proposed species panels, species-tree sources.
- `requested/REQUESTED.md` — 50 paywalled papers / unreachable data files ranked for human retrieval.

## 1. State of the art in orthology benchmarking
- **Quest for Orthologs benchmark service** (Altenhoff 2016 Nat Methods 10.1038/nmeth.3830; Altenhoff 2020 NAR 10.1093/nar/gkaa308;
  Nevers 2022 NAR 10.1093/nar/gkac330; Altenhoff 2024 NARGAB 10.1093/nargab/lqae167) scores *pairwise* orthology: species-tree
  discordance (STD/GSTD), reference gene trees (SwissTree, TreeFam-A), function conservation (GO/EC), and VGNC/HGNC tests.
  HOGs are flattened to ortholog pairs; duplication timing, paralogy and nesting are not scored. A 2025 HOG review states
  that the field "still lacks robust benchmarking tailored to HOGs".
- **OrthoBench** (Trachana 2011; Emms & Kelly 2020 GBE 10.1093/gbe/evaa211): 70 curated RefOGs at one level (Bilateria),
  animal-centred; dominant error is family splitting (best method exactly right in only 23/70).
- **SwissTree** (Boeckmann 2011/2015): manually curated reconciled trees — closest in spirit to us, but few families and
  pairwise scoring only.
- **Tool papers** (OrthoFinder 3: Emms 2026 Nat Methods 10.1038/s41592-026-03126-6; FastOMA 2024/25) evaluate root-level
  groups only. OF3 reports FastOMA over-merges ~10% of reference groups. FastOMA places sequences with OMAmer against the
  OMA database ⇒ references derived from OMA would leak information into FastOMA.
- **Simulation benchmarks** (ALF: Dalquen 2012; Zombi: Davín 2020; SimPhy) know the truth but miss real-data error modes
  (annotation, contamination, real rate heterogeneity); literature benchmarks are realistic but incomplete and biased toward
  well-studied families — no negative set ("absence of evidence" problem).

## 2. Metrics available for hierarchical orthology
Pairwise ortholog P/R; cluster metrics (ARI, VI, B-cubed); RefOG-style per-family P/R/F with fission/fusion counts
(OrthoBench); labelled tree distances for reconciled trees (labeled RF: Briand 2020, 2022; Treerecs; path/triplet measures);
event-based scoring is essentially absent. Known biases: reconciliation pushes duplications toward the root when gene trees
are wrong (Hahn 2007 Genome Biol 10.1186/gb-2007-8-7-r141); ILS mimics duplication+loss; GSTD conflates gene-tree and
orthology error; GO/EC tests rely on the ortholog conjecture. **Gap we fill: per-level, event-aware scoring of HOGs against
literature-asserted histories, stratified by the phenomenon that makes a family hard.**

## 3. Documented error sources (what the benchmark must contain)
Annotation errors (up to ~30% effect on performance; fragmented/fused gene models), WGD "traps" with reciprocal ohnolog loss
(pseudo-orthology; Scannell 2006), delayed rediploidisation (salmonid LORe; Robertson 2017), allopolyploid homeologs
(*X. laevis*), HGT (into fungi: URA1; into nematodes: cell-wall enzymes; prokaryote families), long branches (Microsporidia),
promiscuous domains (kinases, KRAB-ZNF, NLR, small GTPases), tandem/lineage-specific expansions (P450s, chemoreceptors),
species-tree uncertainty (ctenophore/sponge, Amborella, Lachancea), contamination (effect on orthology never quantified).

## 4. Usable literature-based family data (machine-readable, verified by surveyors)
| Family | Source | Data | Scope / phenomena |
|---|---|---|---|
| Yeast WGD pillars | YGOB (Byrne & Wolfe 2005); Marcet-Houben & Gabaldón 2015 | truncated Pillars.tab (5,173 pillars; full file requested) + 365 ohnolog pairs | WGD, pseudo-orthology |
| Teleost 3R ohnologs | Braasch 2016 ST23 | 1,606 zebrafish / 1,315 medaka TGD pairs anchored to gar | WGD3R, loss |
| Salmonid 4R | Robertson 2017 (394 genes AORe/LORe + trees) | trees + labels | WGD4R, delayed rediploidisation |
| Kinesins | Wickstead 2010 S3 | 1,624 accessions, family/subfamily | LECA, promiscuous motor domain |
| Myosins | Sebé-Pedrós 2014; Kollmar 2017 | class-labelled FASTA/trees | LECA, domain shuffling |
| ARF/Arl/Sar | Vargová 2021 | 2,286 rows, 31 OG labels | LECA, short fast GTPases |
| Rab | Klöpper 2012 | 1,669 labelled seqs | LECA, GTPase-domain merging |
| Tubulins | Findeisen 2014 | 3,286 named seqs | LECA, conserved |
| MCM2–9 | Liu 2009 | labelled Newick | LECA control |
| Dyneins | Kollmar 2016; Wickstead & Gull 2007 | Newick, RefSeq IDs | LECA, long genes |
| Adaptins/COPI/TSET | Hirst 2014 | 1,413 rows | LECA heterotetramer paralogy |
| DNA pol B | Tahirov 2009 | labelled trees | fusion origin, archaea–eukaryote |
| Kinome | KinBase (Manning 2002/2008, archived) | 3-level labels | promiscuous domain, LSE |
| Opsins | Ramirez 2016 T1; Lagman 2013 | 768 GenBank IDs + names | deep → 2R/3R |
| Sox | Voldoire 2017 | accession table per subgroup×species | 2R/3R |
| KRAB-ZNF | de Tribolet-Hardy 2023 S1 | 467 human loci with ages | merging stress, young |
| APOBEC3 | Münk 2012 | species × Z-domain table | fusion/fission |
| Chemokines | Nomiyama 2013 | ID→name tables | fast, tandem |
| Drosophila P450 | Good 2014 | labelled Newick over 12 spp. | LSE, tandem, ILS-prone |
| Caenorhabditis chemoreceptors | Thomas & Robertson 2008 | 6,769 seqs | extreme LSE |
| Plant homeobox | Mukherjee 2009 | name/acc/class over 10 genomes | plant WGD |
| NLR | RefPlantNLR (Kourelis 2021) | 481 entries + tree | promiscuous, LSE |
| plant bHLH | Pires & Dolan 2010 | 566 rows | short domain, WGD |
| GID1/CXE | Yoshida 2018 | 260 genes | nested family |
| MIPs/aquaporins | Danielson & Johanson 2008 | Nexus trees | HGT (NIP), WGD |
| Ribosomal proteins | Yutin 2012 | copy-number matrix | single-copy control |
| Nitrogenase/Bch | Garcia 2020 | 319 labelled seqs | prokaryote paralogy/HGT |
| β-lactamases | NCBI Reference Gene Catalog | 6,128 entries with class | non-homologous classes (merge decoys), HGT |
| Photolyase/CRY | Mei & Dvornyk 2015 | accession table | HGT, deep |

Practical constraints found: YGOB, KinBase, HomeoDB, RabDB, SNARE DB are offline (partial Wayback snapshots);
Dryad and some publisher supplements block scripts; most curated sets use 2005–2014 gene models ⇒ identifier drift,
requiring sequence-based mapping to current proteomes. No family has a literature-resolved ILS history; ILS must be
represented through species choice (short internodes) and treated as an untagged nuisance.

## 5. Genomes
UniProt Reference Proteomes (canonical, one protein per gene, CC BY 4.0; release 2026_03) incl. the QfO 2026_02 set
(81 species). QfO is thin for plants, non-mammalian vertebrates and lacks WGD yeasts ⇒ we add UniProt reference proteomes.
Phytozome excluded (redistribution restrictions). Five panels proposed (V 14 spp, Y 16, P 13, E 34, K 21); species trees from
NCBI taxonomy with polytomies resolved from published phylogenomic trees.
