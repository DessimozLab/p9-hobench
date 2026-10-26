# Review round 1 — Prior-work expert (orthology benchmarking)

Overall: the merge-level formalism is sound and does fill a real gap. QfO flattens HOGs to pairs, and OrthoBench, LEMOrtho and the OF3/FastOMA papers score one level only. The main risks are that the design ignores a directly reusable curated resource (SwissTree), that metric definitions drift from established ones without justification, and that the nomenclature rules hide tool-derived orthology.

1. **PW-1 — MUST. Reuse SwissTree as a reference source and cross-check.**
   - Change: add the SwissTree gold-standard trees to §4 as a source tier "curated reconciled tree" (≈T1).
     - What it is: 19 families with NHX consensus trees, explicit duplication nodes, evidence trees, `gene_relationships.txt`, sequences (`swisstree.fst`) and a species tree. All are at https://swisstree.vital-it.ch/ftp/ (checked today).
     - Coverage: ST011 GH14 (Viridiplantae, 39 dups), ST014 CASP (Mesangiospermae, 43 dups), ST006 HOX9–14 (Bilateria, 18 dups), ST007 RPS10/S20 (cellular organisms), ST002 NOX and ST004 SERINC (Eukaryota).
     - Mapping: through the sequences, since the IDs follow 2011 proteomes.
   - Rationale: these are the only published, expert-reconciled S/D trees. They fill your weakest cells (plant tandem/WGD depth, Hox), cost little to compile, and let you show that metric 4 restricted to SwissTree reproduces what the QfO SwissTree test would measure (Boeckmann 2011; Altenhoff 2016).
   - If you exclude them, DESIGN must say why.

2. **PW-2 — MUST. Fix factual errors in §4 controls.**
   - β-lactamases: classes A, C and D are homologous serine enzymes (PBP-like fold). Only class B (metallo-β-lactamases) is unrelated (Hall & Barlow 2004). Use only B-vs-{A,C,D} as non-homologous decoys, and mark A/C/D relations above class level as U.
   - Ribosomal proteins: Yutin 2012 documents paralogs (e.g., zinc-ribbon L31/L33/L36 variants). Restrict the "single-copy control" to proteins the copy-number matrix shows as single-copy in every panel species.

3. **PW-3 — MUST. Close the nomenclature circularity loophole (§2.1 T2, §4 V).**
   - Problem: "gene names in curated nomenclatures" is not tool-free. VGNC sets are seeded from HCOP consensus of Ensembl Compara, OMA, PANTHER and others (Nevers 2022). ZFIN and NCBI names often inherit from them too. That contradicts §4's exclusion of OMA, OrthoFinder, Compara and PANTHER sources, and it favours FastOMA through OMAmer.
   - Change: record the evidence basis per nomenclature source and accept only synteny- or phylogeny-based names. Examples that qualify: Xenbase .L/.S from subgenome synteny (Session 2016), Hox cluster position, YGOB synteny pillars. Drop or demote prediction-seeded names.

4. **PW-4 — MUST. Fix within-family pair weighting.**
   - Problem: `metrics.py` counts raw pairs, which grow as n². Inside a kinesin or KRAB-ZNF family, the largest subfamily dominates every per-level F1. OrthoBench divides TP/FP/FN by (n−1) per RefOG (Emms 2020) for exactly this reason. The survey's own recommendation (A, Implications #2: per-level P/R plus B-cubed) was dropped without explanation.
   - Change: make per-level **B-cubed P/R** (Amigó 2009, per gene, size-aware) primary next to pair F1, or apply the (n−1) normalisation. Document the choice.

5. **PW-5 — SHOULD. Use established root-level definitions.**
   - Change: replace the ad-hoc "#HOGs covering ≥90%" in metric 5 with the OrthoBench/OF3 measures computed on root-level HOGs per family: precision, recall, F, % fission, % fusion, % exact, and entropy of RefOG dispersion (Emms 2020; Emms 2026). Keep `split_excess` as a supplement.
   - Validation: unit-test equivalence against `benchmark.py` from github.com/davidemms/Open_Orthobench on a toy partition.
   - COULD: run OF3 on the 12 OrthoBench proteomes (~1 h) and reproduce the published figures as an end-to-end check of parsing and metrics.

6. **PW-6 — SHOULD. Anchor the formalism in theory and state the relation to tree distances.**
   - Ground it: the reference merge-level matrix is a symbolic ultrametric whose orthology relation is a cograph (Hellmuth 2013; Lafond 2014). That is why the validator's transitivity checks are needed. Predicted merge levels from OrthoXML follow pyHam/QfO flattening (Train 2019; Altenhoff 2016). Cite these in §2.
   - Relate: merge-level agreement is a pair-based reconciled-tree comparison, close to path-label reconciliation (López-Sánchez 2024).
   - Survey recommendation #6 (LRF via pylabeledrf, Briand 2022) was dropped. Either add LRF as secondary on fully resolved references (SwissTree, MCM), or justify dropping it (HOG-derived trees are multifurcating, and LRF penalises polytomies).

7. **PW-7 — SHOULD. Make "QfO-comparable" (metric 4) precise.**
   - QfO scores OrthoFinder from `Orthologues/` (gene-tree-based pairs), not from `N*.tsv` HOGs. Compute OF orthology both ways and report the difference.
   - Also adopt VGNC's conservative false-positive rule explicitly: an FP counts only when both genes are in declared-complete reference sets (Nevers 2022).

8. **PW-8 — SHOULD. Event metric: add spurious events and an explicit-event cross-check.**
   - Problem: Hahn 2007 and SimPhy (Mallo 2016) predict extra, root-ward duplications. "Recovered ≥0.5 of cross-copy pairs" measures only the missed side.
   - Change:
     - Report a spurious-duplication count: predicted duplications inside reference S clades.
     - Cross-check placements against tool-native events (pyHam per-branch duplications; OF `Gene_Duplication_Events/Duplications.tsv`).
     - Report sensitivity to the 0.5 threshold.

9. **PW-9 — SHOULD. Make level error comparable across panels and ILS-aware.**
   - Report signed error also in Myr (TimeTree) or as a fraction of family depth, since an edge count depends on panel density.
   - Add a tolerance variant that ignores ±1-level errors across internodes shorter than 10 Myr, mirroring QfO's GSTD collapse (Altenhoff 2016). This is the only ILS handling you can justify, since no family has a literature-resolved ILS history.

10. **PW-10 — SHOULD. Keep genome-scale collections in their own stratum.**
    - YGOB pillars, Braasch 3R pairs and Marcet-Houben ohnologs are thousands of small units. Report them separately from the macro-average over curated families, as Salichos & Rokas 2011 and VGNC report.
    - Reproduce Salichos 2011's "trap" analysis (pillars with reciprocal ohnolog loss) so the yeast results can be compared directly.

11. **PW-11 — COULD.** Add FastOMA's reference-free "implied losses" (Majidian 2025) as a difficulty covariate only, never as accuracy. GSTD and GO/EC tests stay excluded, which is correct as is.

12. **PW-12 — MUST (documentation). Correct LITERATURE.md.**
    - §1: SwissTree is Boeckmann 2011 (Brief Bioinform 12:423). Boeckmann 2015 is the QfO species-tree paper.
    - §1 lists "VGNC/HGNC tests", but survey A found no separate HGNC test.
    - §2 says event-based scoring is "essentially absent". SwissTree annotates duplication nodes but QfO scores them only pairwise; say that instead.
    - The gap claim should read "no **real-data** benchmark scores HOGs per level". GETHOGs and OMA 2.0 evaluated hierarchical groups on ALF simulations (Altenhoff 2013; Train 2017).

No changes needed in the curation blinding, the frozen-reference tag, the OMAmer leakage experiment (E3) or the no-composite-score policy. These match QfO practice (Altenhoff 2016).
