# E03_Rab curation log
status: curated, protocol v1 (headline-eligible: 87.4% of in-panel sequences uniquely mapped; 388 members, 15 species)

## Read
- Klöpper, Kienle, Fasshauer & Munro 2012 BMC Biol 10:71 (PMID 22873208, PMC3425129, OA CC BY) via fetch_paper ->
  literature/papers/PMID22873208/ (fulltext.txt; raw/supp Additional files 3-9). Additional files 4-9 = FASTA alignments
  of Rab supergroups I-VI (1,669 seqs, CR line ends, names <Species>_<RabLabel>[_n], no accessions); file 3 = species
  short names (HoSa, DrMe, NeVe, AmQu, TrAd, MoBr, CaOw, SaCe, ChRe, ArTh, PlFa, TrBr, NaGr, ThPs, PhRa = panel species).
- Elias et al. 2012 J Cell Sci (PMID 22366452): fetched (PMC HTML only, not OA subset, supplements not retrievable);
  not needed for assertions (Klöpper has explicit LECA list + supergroup expansions). TreeBASE S12811 not used (trees
  are not needed for T1/T2 encoding).

## Decisions
- Mapping: no accessions -> same-species DIAMOND of the degapped alignment sequence vs panel proteome. The authors removed
  columns with >25% gaps, so internal segments are missing and appear as alignment gaps; identity is therefore computed
  over aligned residue pairs (nident / (length - gaps)); coverage of BOTH sequences by aligned span >= 80%; RBH;
  margin >= 2 points. With gap-inclusive identity only 32/444 would pass although HUMAN hits are 100% over residues.
  Two sequences hitting the same entry -> both ambiguous.
- Tree (T1 quotes in assertions.tsv):
  - root D (lvl_lo=Eukaryota, lvl_hi=ABOVE_RANGE) over 6 supergroups ("twenty LECA Rabs fell into a set of six larger
    supergroups, suggesting a primitive pre-LECA eukaryote with just six Rabs"); supergroup of each sequence = its file.
  - each supergroup: U wrapper over a T1 D splitting its LECA Rabs (list quote) + loose genes whose label is not
    quote-assigned to a single LECA Rab (yeast Ypt1/Sec4/Ypt7/Ypt51/Ypt53/Ypt3/Ypt6; Rab19/33/34/44/45 which the paper
    says it "could not link unambiguously"; Rab43, Rab25, Rab36, Rab37, Rab42, RabX2; composite/partial names).
  - LECA Rab node = U over label groups (core + quoted metazoan-derived labels, e.g. Rab8 <- Rab3/10/15/26/27/X4/12/13);
    undated metazoan splits are not asserted.
  - nested T1 D: Rab9 vs opisthokont Rab7 at Opisthokonta ("Rab9 appeared with the opisthokonts"); Rab40 vs bilaterian
    Rab18 at Bilateria ("duplication in bilateria (Rab40)"); vertebrate duplications as D over the HUMAN copies (no lvl;
    level = HUMAN branch): Rab5a/b/c+Rab17, Rab8a/8b+Rab13+Rab12 (deuterostome+vertebrate), Rab22a/b, Rab32+Rab38,
    Rab4a/b, Rab6a/b/c+Rab41. Parent S nodes get the expected verifier warning (>1 HUMAN member, under the D).
- T2 groups: <=1 gene/species -> S; else U (rules b/c; no two panel-E species share a WGD). Label = name suffix with the
  copy index removed; lineage-specific letter suffixes (Rab8A/B, Rab11A1...) are pooled into the base label because the
  same letters are reused independently in different lineages.
- complete_species: none (alignments are HMM-training subsets; no per-species completeness statement).
- source_method HMM_classification (iterative HMM + phylogeny, own curated models; not denylisted).
- Homology: Rab and Arf/Sar homologous (Ras superfamily; quote) -> no distinct_from with E02_ARF.

## Problems
- AMPQE 73%, PHYRM 76%, DROME 81% mapped (old gene models, partial/EST sequences, trimmed long proteins <80% target cov).
- Several verifier 'quote does not name leaf label' warnings: quotes name the family (e.g. 'Rab22') while human leaves
  are 'Rab22a/b'.

Time: ~1 h.
