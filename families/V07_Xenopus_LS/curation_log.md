# V07_Xenopus_LS — curation log (collection XLAEVIS)

Status: **done** (panel V), protocol v1. Curator agent, 2026-10-26/27. Time ~1.5 h.

## Source
Session et al. 2016 Nature (PMID 27762356, doi 10.1038/nature19840, CC BY). Main text via fetch_paper
(`literature/papers/PMID27762356/fulltext.txt`). Supplements delivered by the human (`requested/08_*`, sha256 in
`requested/RECEIVED.tsv`); copied/extracted to `literature/papers/PMID27762356/raw/supp/` (Supplementary Table 7 xlsx and
Tables docx extracted from 08_supp2.zip; 08_supp1.pdf copied). Paper gene-model peptides (needed for sequence mapping)
from Xenbase: `raw/xenbase/XL_9.1_v1.8.3.2.primaryTranscripts.pep.fa.gz` (45,099 models = paper count) and
`Xtropicalisv9.0.Named.primaryTrs.pep.fa.gz`. sha256 of all files in collection.yaml.
Member table: Suppl. Table 7, sheet `XT-XL.ortho`, 8,806 rows `homeolog` (= the paper's 8,806 homeologous pairs).
Basis (Suppl. Note 8.1, quoted): BLASTP best hits, then "only accepted X. laevis homoeologues whose synteny and
subgenome identity agree with the BAC FISH map" → source_method synteny. Not tool-derived orthology DBs (no denylist
source). UniProt/Xenbase gene names are used only as a fallback mapping key, not as evidence.

## Encoding (T2)
Unit = XT locus + L + S: `((L,S)[D tag=ALLO], XT)[S]`, root_level Xenopus. D has no `lvl` (defaults to the XENLA
leaf; verifier only accepts internal levels for lvl) → L/S apart at XENLA, merge at Xenopus; L–XT, S–XT orthologs.
ALLO quote: "The L and S chromosome sets therefore represent the descendants of two distinct diploid progenitors…".
complete_species: [] (pairs are a filtered subset, not a genome-wide complete list).

## Mapping (same species only)
Per paper protein: exact sequence (unique) → DIAMOND RBH ≥99 % id, ≥80 % cov both, margin ≥2 → else exact unique gene
name (rejected if it disagrees with the best sequence hit; NCBI LOC names never used as a name key).
Overall source slots: 20,904 mapped / 26,418 (unique_map 0.83 per verifier); 3,489 unmapped, 801 ambiguous, 1,224 excluded.
Unit filters: both .L and .S named (819 rows have an unassigned `.1`-type copy → rejected), all 3 mapped, L≠S sequence,
homology QC (L→XT and S→XT DIAMOND best hit = unit XT, ≥70 % id, ≥50 % qcov; 159 rows failed — name drift guard).
4,671 passing units → seeded sample of 500 (random.Random(20261026), units sorted by id); others `not_sampled`.
Members: 1,500 (XENLA 1,000, XENTR 500); methods exact 989, sequence_rbh 310, gene_name_exact 201.

## Verification
hobench validate_family: 500/500 units 0 errors/0 warnings. verify_provenance.py on an assembled single-family view
(collection.yaml as family.yaml + unit trees under a root U; `--root` = project): PASS, 0 errors, 0 warnings (a–f);
1,005/1,005 quotes found. build.py deterministic (members.tsv identical across reruns). g (--rebuild) not run.

## Caveats
- The 8,806 pairs are the paper's conservative 2:1 set; genes in larger families / tandem copies are absent.
- Session's XT orthologue assignment is BLAST + synteny, not a phylogeny; tier T2.

## Erratum 2026-10-27 (post-freeze, PI request)
- mapping_coverage now written as a fraction (min of query/target coverage, e.g. 1.000, 0.964), as hobench/panel.py
  requires; identity stays in percent. No mapping decision changed.
- Units sharing any member gene with another family are dropped (other curation is richer): XL_hoxa4, XL_hoxb7
  (genes in V03_Hox). Recorded in mapping_report.tsv unit_status `dropped: shares gene(s) with another family (...)`
  and collection.yaml sampling.dropped_shared_with_other_families. The sample itself is unchanged (no re-draw):
  500 -> 498 units, 1,500 -> 1,494 members. build.py reads families/*/members.tsv and snapshots the overlapping genes
  to shared_with_other_families.tsv (used when families/ is absent, e.g. isolated rebuild; verified byte-identical).
- Re-validated: hobench validate_family + panel validate_members 0 issues on all 498 units; verify_provenance (assembled
  view) PASS 0/0.
