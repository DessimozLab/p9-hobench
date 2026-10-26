# Round 1 review: External PI (comparative genomics)

**Overall.** Sound and timely: merge levels, tiers and event recovery fill a real gap (QfO flattens HOGs to pairs). Risks: breadth vs ~22 h, decoy/contamination errors, no adoptable release form.

1. **PI-1 (MUST): build one pipeline end to end before scaling.** In DESIGN.md §11, add a gate at T+6:
   - 2 frozen families each in Y and K (e.g. a YGOB pillar subset + URA1; ribosomal proteins + EF-Tu/EF-G);
   - each scored end to end with both tools, the baselines and the controls.

   Fan out curation only after this works. Add a second gate at T+12: freeze what has passed verification. Nothing curated after T+14 enters the headline set. *Why:* avoids 30 half-curated, unscored families.

2. **PI-2 (MUST): fix decoy and contamination semantics.**
   - (a) β-lactamase classes A, C and D are distant homologs (serine-β-lactamase/PBP fold). Only class B (MBL fold) is non-homologous to them. In DESIGN.md §4 K, declare only B vs {A,C,D} as decoys, and set A/C/D to `U` at the root.
   - (b) Metric 5 (ENGINEERING_NOTES, "Metric 5") counts any other reference family as contamination. But Rab/ARF and kinesin/myosin are homologous superfamilies. Merging them is wrong only at levels where the literature says they were distinct genes. Add explicit assertions in family.yaml (`distinct_from: {family, at_level, src}`). Score contamination only against these assertions plus true decoys.

   *Why:* as written, the metric penalises correct deep merges.

3. **PI-3 (MUST): name three headline metrics in §6 and keep the rest as secondary.**
   - (i) Per-level HOG F1: depth profile plus macro mean over levels.
   - (ii) Event recovery per named literature event (yeast WGD, TGD, 2R, salmonid 4R, LECA paralogs).
   - (iii) Family integrity (split_excess, contamination).

   Present them QfO-style: precision–recall scatter with a Pareto front, no composite. The per-event table is the most novel and citable output.

4. **PI-4 (MUST): define the unit of evaluation for genome-scale collections.** YGOB pillars (~5k) and Braasch ST23 (~1.6k pairs) are not families. Report each as a separate *genome-scale track*, macro-averaged over its own units (pillar or ohnolog pair). Never average them with kinesins.
   - Add `source_method` (manual phylogeny / synteny / classification / tool-assisted) to `sources` in family.yaml.
   - Braasch and OHNOLOGS are synteny- or Compara-assisted (B1 §V2/V3).

   *Why:* the aggregates stay valid and circularity becomes auditable (§10).

5. **PI-5 (MUST): release the benchmark input itself.** Release:
   - the exact input FASTAs, in both full-panel and subset mode, with md5s (CC BY 4.0 permits this);
   - the panel species trees and a manifest.

   Call subset mode "HOBench-core" and make it the official input, not a shortcut. *Why:* external methods need identical inputs, since QfO counts missing predictions as false negatives (A §1).

6. **PI-6 (SHOULD): interoperability formats.**
   - (a) Input: accept OrthoXML 0.3/0.5 and a documented TSV (`level, hog_id, gene`). The reader exists as `predictions.read_tsv`; document it.
   - (b) Reference export: recPhyloXML (the standard for reconciled trees) alongside NHX. Also export per-level reference HOGs in OrthoXML, with U regions marked via `<property>`. Also export flattened ortholog/paralog pairs, which QfO could use as a SwissTree-like test.
   - (c) Name each level with NCBI taxid + name where an exact clade exists, so FastOMA, OMA and pyHam match.

7. **PI-7 (SHOULD): add a QfO projection and more methods.**
   - Make reference trees prunable to the QfO 2026_02 species and tree. D levels are re-mapped by clade species sets.
   - Two methods are too few for the inter-method agreement in §6. Add precomputed HOGs as methods, never as truth: OMA All.May2026, or Ensembl Compara on V. Add cheap variants: OrthoFinder `-M msa`, and FastOMA with vs without a user species tree.

   Together these give 4–6 methods at little compute.

8. **PI-8 (SHOULD): QC references independently of the tools before the freeze (§5).** Run an all-vs-all DIAMOND search within each family. Flag members with no hit to any other member, and members whose best hit lies in another family. Fix or document each flagged member. Curators stay blind to the tools.

9. **PI-9 (SHOULD): narrow E2 and E3.**
   - Keep three perturbations, on Y plus one subset panel, 3 seeds each:
     - (a) fragmentation;
     - (c) contamination, whose effect has never been quantified;
     - (d) species-tree error, which is specific to hierarchical methods.
   - Downgrade (b) fusion and (e) taxon sampling to COULD.
   - Drop E3: an OMAmer database without the panel species cannot be rebuilt in this budget. State this as a limitation instead, and report the overlap between FastOMA root HOGs and OMA HOGs descriptively.

10. **PI-10 (SHOULD): triage the panels.**
    - Drop D.
    - Keep β-lactamases only if chromosomal copies (blaPDC/ADC/OXA-51-like) occur in K-panel proteomes. Plasmid alleles in strain-specific proteomes will mostly be absent.

11. **PI-11 (SHOULD): add a coverage table.** Compare clade × depth × phenomenon for HOBench, OrthoBench (70 RefOGs, Bilateria) and SwissTree.

12. **PI-12 (SHOULD): keep the statistics honest.**
    - Report n per stratum.
    - Give bootstrap CIs only when n ≥ 5 families; otherwise show individual points.
    - Make paired comparisons on shared families only.

13. **PI-13 (COULD): packaging.** Reserve 1.5 h at T+20 for:
    - a README quickstart (`hobench evaluate --orthoxml x --panel Y`);
    - per-family cards (sources, tiers, mapping stats, caveats);
    - LICENSE (code MIT, data CC BY 4.0), an environment file and CITATION.cff;
    - a Zenodo-ready layout.

    OpenEBench container: follow-up.

**No change needed:**
merge-level formalism, T1–T3 tiers, blind frozen curation, 10% re-derivation, controls, family bootstrap.
