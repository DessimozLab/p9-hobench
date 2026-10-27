# CLAIM_CHECK: hobench_manuscript.tex, S1_Text.tex, S2_Text.tex, README.md

Independent check against repo files. Verdicts: OK / WRONG / UNSUPPORTED / IMPRECISE. Documents were not edited.
Not verifiable here: the test count of 95 (pytest not installed; 82 `def test_` plus 2 parametrized blocks, REPORT_CLAIM_CHECK.md says 95 passed), the 469 MB archive size, "no identifiers typed by hand", reviewer/verifier blindness, and the claim that curators could not open tool outputs.

## A. Claims that need correction (WRONG / UNSUPPORTED / IMPRECISE)

| # | Claim | Location | Source | Verdict | Correction |
|---|---|---|---|---|---|
| 1 | "HGT has one headline family" | Discussion, limits | results/analysis/phenomena.md (HGT: PY, 2 items); Table 2 (Y03, P04) | WRONG | Two headline items carry HGT (Y03 URA1, P04 MIPs/NIPs). |
| 2 | "Five species panels (90 species...)" | Abstract | panels/*/species.tsv | IMPRECISE | Panels sum to 102 species slots (16+21+34+17+14); 90 are unique taxa, because K and E share prokaryotes. Say "90 unique species". |
| 3 | "In both tools about half of the contaminants went into HOGs with yeast reference genes" | Robustness | results/Y_contam200_full/summary.md | IMPRECISE | OrthoFinder 103/200 (52%), FastOMA 133/200 (67%). |
| 4 | "Kept each copy's orthologs together in 79% (OF) / 57% (FastOMA)" of events | WGD events | phenomena.md 2: "min cohesion" 0.786 / 0.572 | IMPRECISE | These are mean minimum fractions of each copy's orthologs that share a HOG, not percentages of events. The 26%/64% separation values are also means of a per-event fraction (Fig 3 caption says it correctly). |
| 5 | "FastOMA ... split 56% of [X. laevis units]" | Vertebrate duplications | error_anatomy.md (b): over-split 0.571; REPORT 55.6% | IMPRECISE | Two sources differ (57% vs 55.6%). Use one, and say which node type (D nodes, 333 units). |
| 6 | "Running trivial baselines on real data exposed five metric defects" | Metrics; S1 s3; Discussion | results/ERRATA.md E5; S2 text | IMPRECISE | E5 (membership F1, unresolved root) came from an error analysis of P02, and the collapsed-items defect from the evaluator. S2 says so; the main text and S1 do not. Say "controls and error analysis exposed five". |
| 7 | "For all tools and testable primary endpoints, subset minus full was within +-0.05" | Subset section | results/equivalence.md | IMPRECISE | K event recovery (n=1 event, delta +1.0) is "no" for OrthoFinder and of_flat; the notes call it inconclusive. State the exception. "Pre-registered" means committed in the repo before runs, not registered externally. |
| 8 | "Of 1,592 YGOB WGD events, 202 are testable"; "201 testable YGOB WGD events" | Metrics; Phenomena; Table 4; S1 | results/ERRATA.md E1; Y01 curation_log | IMPRECISE | 1,592 = 1,591 YGOB WGD events + 1 URA1 HGT event; 202 = 201 WGD + 1 HGT. Use 1,591 / 201 consistently, or call 1,592/202 "Y-panel events". |
| 9 | "about 31,800 genes, 85 T1, 2,138 T2" | Reference families | results/catalogue/families.tsv; families/ERRATA.md | IMPRECISE | Catalogue is stale for V07 (500 units, 1,500 members). After the documented correction to 498 units, T2 = 2,136 and genes = 31,761. Regenerate the catalogue or note it. Table 2 already uses 498/1,494. |
| 10 | "The agents ... ran autonomously for 24 hours" | Methods, Agent setup | RESEARCH_TREE.md M7: closed at T+20h35; "about 21 hours" in Discussion | IMPRECISE | 24 h was the budget; the run lasted about 21 h. Fix the inconsistency. |
| 11 | "49 of 50 ... supplied"; S1: "A human delivered 49 items within about 1.5 hours" | Methods; S1 s1 | RECEIVED.tsv (49 ranks); LESSONS.md ("~90% in 1.5 h ... 45/50"); RESEARCH_TREE.md ("~45/50") | IMPRECISE | RECEIVED.tsv supports 49 delivered in total; the 1.5 h figure applies to about 45 (90%). Repo docs disagree with each other (45 vs 49). |
| 12 | "22 of 23 verification groups passed" | Curation | VERIFICATION.md (23 groups, 22 pass) | OK, but RESEARCH_TREE.md says "21/22" (older) | Update RESEARCH_TREE.md or mention it was superseded. |
| 13 | K01 "40 units" (Table 2) | Table 2 | VERIFICATION.md ("30 of 60 units"); families.tsv n_units 40; S2 "30 units" | IMPRECISE | Verification file says 60 units; catalogue says 40. Reconcile. |
| 14 | Hahn 2007 cited as agreeing with the direction of OrthoFinder's WGD error | WGD events | hahn2007bias (bias towards inferring spurious, often deeper, duplications) | UNSUPPORTED | OrthoFinder here merges ohnologs (too few/too young duplications), the opposite direction of the reconciliation bias. Drop or reword ("a different, documented bias"). |
| 15 | "such tests separate methods in ways that pairwise benchmarks do not show"; "A pairwise benchmark that pools these pairs would show a single high number" | Discussion | no pairwise benchmark was run | UNSUPPORTED | Hypothesis only. Say "would probably". Pooled correctness here is 0.88-0.91 for either tool; compute it if the claim stays. |
| 16 | "specific, repeatable failure modes" | Author summary | single runs; FastOMA replicate not bit-identical (2% units differ); one species pair for hidden paralogy (49 pairs) | UNSUPPORTED | Remove "repeatable" or say "consistent across the tested runs". |
| 17 | "Existing benchmarks do not test whether these hierarchies are correct at each level" | Author summary | Intro hedges: "We found no real-data benchmark..." | IMPRECISE | Match the Intro hedge. Same for Abstract "Existing benchmarks score HOGs only through flattened views". |
| 18 | "They failed on WGDs in opposite ways" | Abstract | phenomena.md 2c: V Hox WGD, FastOMA separation 0.73-0.82 below OrthoFinder (0.82-0.92), and cohesion lower | IMPRECISE | Holds for the 201 yeast events only. In vertebrate WGD events FastOMA is worse on both. Add "in yeast". |
| 19 | "Gene fragmentation lowered accuracy more than contamination, species-tree changes or taxon removal" | Abstract | robustness.md: one seed, panel Y, 10% fragments, one tree node, one clade | IMPRECISE | Depends on credit rule (OF -0.050 strict vs -0.003 lenient). Add "in yeast, at the tested strength". |
| 20 | "FastOMA's small deltas are therefore within its run-to-run noise" (S1: "no effect beyond FastOMA's run-to-run noise") | Robustness; S1 s5 | noise_floor_fastoma.tsv: mean -0.001, per-unit SD 0.038 (n=1,977) | IMPRECISE | The SD is per unit; with n=1,977 the mean's standard error is about 0.001, so alt-tree merge delta -0.0035 is not clearly inside noise. Say "of the same order as". |
| 21 | "the relations among the yeast genes did not change" (200 contaminants); S1 "still did not change reference relations" | Robustness; S1 | Y_contam200_full: FastOMA 50 units changed in the 20-protein test; deltas -0.001 | IMPRECISE | Mean deltas were negligible, but relations did change in some units. Say "changed by less than 0.001 on average". |
| 22 | "A recent review notes that HOG-specific benchmarking is still missing" | Intro | sartonloheac2025hogs | Could not verify | Check the review states this; otherwise cite as general context. |
| 23 | Lafond & El-Mabrouk 2014 cited for "merge-level matrix is a symbolic ultrametric" | Metrics | lafond2014orthology (orthology constraints, satisfiability) | IMPRECISE | Hellmuth 2013 supports the symbolic ultrametric statement; Lafond supports orthology-relation consistency only. |
| 24 | S2 path `design/MUTATION_REPORT.md`; results/ERRATA.md same | S2 header | actual: docs/MUTATION_REPORT.md | WRONG | Fix path. |
| 25 | `\cite` missing in S1 ("[Marcet-Houben and Gabaldón 2015]", "[Omelchenko et al. 2010]") | S1 | S1 uses bracketed text | IMPRECISE | Not linked to the reference list. |
| 26 | UniProt 2025 paper cited for release 2026_03 | Table 1, Methods | uniprot2025 | IMPRECISE | Release postdates the paper; add a release note or cite the release page. |
| 27 | Stated "31 subagents" | Methods | list sums to 31 without the 7 reviewers launched by the coordinator | IMPRECISE | Say "31 subagents plus 7 reviewers launched by the coordinator" (or 38 total). |

## B. Claims checked and OK

| Claim | Source |
|---|---|
| Table 3 values (13/13/11/5 items; +0.054/+0.134/+0.285/+0.315; Holm 0.035/0.058/0.12/0.12) | results/confirmatory.tsv |
| Pooled merge accuracy 0.886/0.846/0.712 (n=13); event recovery 0.617/0.479/0.194 (n=11); "not on every panel" (V) | main_table.md |
| Panel E: one big HOG 0.838, OF 0.729, FastOMA 0.531; S-pair over-split 66.7%/20.4%; FastOMA over-merge 24.4% (OF 0.1%), mostly Eukaryota level | main_table.md, error_anatomy.md |
| bHLH FastOMA precision 0.959, recall 0.025, 186 groups | REPORT.md, error_anatomy.md (a) |
| Hidden paralogy: 49 pairs, 4%/12% correct; 1,189 true pairs, 91%/87% | phenomena.md s4 |
| YGOB D-S accuracy gap 0.68 (OF) and 0.21 (FastOMA), n=532 | phenomena.md s3 |
| Hox 2R 8% / 3R 9% (FastOMA), 77% / 82% (OF); 13 and 11 events | phenomena.md 2c |
| Fragmentation -0.050/-0.114 (both) and -0.003/-0.024 (either); contamination/alt-tree/taxon removal at most 0.005; contam200 +0.0004/-0.0010 | robustness.md, Y_contam200_full |
| OF replicate identical; FastOMA replicate differs in 2% of units, mean -0.001, SD 0.038 | noise_floor_*.tsv |
| Table 1 species and protein counts; subset sizes 25,855 / 70,000 / 25,429; Y subset 54%, K 6% | panels/*/species.tsv, subset/ids |
| Pair identity 96-99.6% | equivalence.md |
| Catalogue split 13 headline / 3 controls / 10 decoys / 6 secondary = 32 | families.tsv |
| Table 2 members and species counts for all rows | families.tsv |
| 87% of YGOB events untestable (1,390/1,592) | results/ERRATA.md E1 |
| Verification table in S2 (all 23 rows, bounds, V02 resample, V01 fail) | VERIFICATION.md |
| 95 tests / 10 of 10 mutants killed | MUTATION_REPORT.md, REPORT.md |
| Tool versions OrthoFinder 3.1.5, FastOMA 0.5.0, DIAMOND 2.2.8, OMAmer 2.1.0 | runs/*/manifest.json |
| 13 subagents, 3.8 h lost; outages 2h50 + 4h20 (about 7 h); about 11 h lost in total; about 85 candidates | RESEARCH_TREE.md |
| Difficulty analysis: rho -0.41, p=0.19, n=12; earlier -0.58 with n=20 | difficulty.md, RESEARCH_TREE.md |
| V07 500 to 498 correction; Y03 quote demoted; LOC audit note | families/ERRATA.md |
| E1, E2, E4 precede E/P/V scoring; E5 from P02 error analysis | results/ERRATA.md |
| Kinesin domain rule rescued 38 (361 to 399) | E01 curation_log |
| README file paths and commands all exist (score_method.py, build_subsets.py, run_all.sh, DATA.md, AUDIT.md, MANIFEST.sha256) | filesystem |
| All \cite keys exist in references.bib (7 bib entries unused: braasch2016gar, hoffmann2012globins, kollmar2017myosin, mukherjee2009homeobox, nevers2025omark, robertson2017salmonid, voldoire2017sox) | references.bib |

Citation support: all other citations match their sentences by title (QfO, OrthoBench, SwissTree, ALF, GETHOGs, pyHam, OMAmer, DIAMOND, Holm, TOST, Clopper-Pearson, B-cubed, AMI, family sources in Table 2). Doubtful ones are rows 14, 22, 23, 26.

## C. Over-claims and tensions with the stated limitations

- Author summary and Discussion generalise from 13 items and one pooled test (p=0.035, items not independent, E and V drive the gain). The Abstract is hedged; the Author summary ("found specific, repeatable failure modes") is not.
- "Only the pooled test is inferential" (Discussion) agrees with the text, but Phenomena subsection percentages (49 pairs, 201 events, 13 Hox events) have no tests and rest on one species pair or one family; say so.
- Robustness: Abstract drops the "mild perturbation, single seed, Y only" qualifier that the Results include.
- "Pre-registered" is used for choices committed in git by agents, with no external registry; two items (S/D-balanced accuracy, E3 reporting changes) were post-results.
- "HOBench shows ... separate methods in ways that pairwise benchmarks do not show": untested (row 15).
- Causal wording: "OrthoFinder merged ohnolog copies, FastOMA split each copy's orthologs" is descriptive and fine; "mistaking paralogs for orthologs after WGD" implies mechanism but is a classification outcome. The reconciliation-bias sentence is already hedged ("we did not test the cause") but cites a work that predicts the opposite direction.

## D. Clarity list (edit-for-clarity rules 1 to 6; 22 instances)

| # | Rule | Location | Sentence (shortened) | Suggested rewrite |
|---|---|---|---|---|
| 1 | 2 | Abstract | "Existing benchmarks score HOGs only through these flattened views." | OK as is; but "Contested event timing gives an interval" (Fig 1) becomes "When event timing is contested, the reference allows an interval." |
| 2 | 2 | Intro/Concept | "This design addresses the main weakness of literature benchmarks" | "We designed the reference to avoid the main weakness of literature benchmarks." |
| 3 | 2 | Metrics | "Requiring both sides prevents trivial predictions from scoring." | "Because both sides must hold, trivial predictions cannot score." |
| 4 | 2 | Metrics | "Running trivial baselines on real data exposed five metric defects" | "When we ran trivial baselines and error analyses on real data, we found five metric defects." |
| 5 | 2 | Tool comparison | "Differences in the other endpoints pointed the same way but were not significant." | "The other endpoints differed in the same direction, but not significantly." |
| 6 | 2 | Robustness | "These perturbations were mild, so they show that the tools are not sensitive to them" | "The perturbations were mild. Neither tool changed much, but we cannot say how larger errors would affect them." |
| 7 | 2 | Phenomena | "YGOB synteny tracks show which single-copy genes ... are paralogs" | "Using YGOB synteny tracks, we identified single-copy genes in two post-WGD species that are paralogs after reciprocal loss." |
| 8 | 6 | Discussion | "the interval for the yeast WGD, together with copies that consist of a single gene, makes 87% of YGOB events untestable" (subject 17 words from verb) | "Fourth, some choices are ours and debatable. The yeast WGD interval and single-gene copies make 87% of YGOB events untestable." |
| 9 | 6 | Methods/Curation | "The rule that genes with the same label and at most one gene per species are orthologs (an S node) was not applied across WGD lineages" | "We did not apply the rule 'same label and at most one gene per species means orthologs (an S node)' across WGD lineages." |
| 10 | 6 | Metrics | "The merge-level matrix of a reconciled tree is a symbolic ultrametric [cites], This is why..." | Fix the comma+capital: "...symbolic ultrametric [cites]. Hence per-level HOG partitions, pairwise relations and duplication timing can all be read from one matrix." |
| 11 | 5 | Metrics | "Per-level B-cubed F1 ... Family-membership F1: precision and recall of family membership" | "Family-membership F1: precision and recall for which genes belong to the family." |
| 12 | 5 | Tool comparison | "the pooled confirmatory test (paired permutation, OrthoFinder default versus FastOMA, Holm correction over the four endpoints)" | Split: "Our confirmatory test pools all headline items. It is a paired permutation test of OrthoFinder default against FastOMA, with Holm correction over the four endpoints." |
| 13 | 5 | Abstract | "contamination, species-tree changes or taxon removal" | "adding foreign proteins, changing the species tree, or removing species." |
| 14 | 5 | Subsets | "a capped sample of domain-only hits" | "a limited number of proteins that share only a domain with reference members." |
| 15 | 5 | Table 4 | "Allopolyploid interval for the yeast WGD" | "Interval for the timing of the yeast WGD, which was an allopolyploidy." |
| 16 | 5 | S1 | "OMArk completeness and contamination scores" | "OMArk scores for completeness and contamination." |
| 17 | 3 | Methods | "the blind verifier re-derived min(all, 30) seeded assertions" | "the blind verifier re-derived a random sample of up to 30 assertions per family (fixed seed)." |
| 18 | 3 | Phenomena | "ohnolog copies may coalesce at the ancestor of the post-WGD and ZT clades" | "ohnolog copies may trace back to a single gene in the ancestor of ..." ("coalesce" is not used this way for paralogs). |
| 19 | 3 | Results | "FastOMA ... spread the 435 members over 186 groups" | "FastOMA placed the 435 members in 186 groups." |
| 20 | 3 | Discussion | "the account reached its usage limit three times" | "we hit the account's usage limit three times." |
| 21 | 1 | Intro | "weakness of literature benchmarks" (could mean benchmarks about literature) | "weakness of benchmarks built from published sources." |
| 22 | 1 | Discussion | "ideas taken from the literature" paragraph: "level induction in OrthoXML from pyHam" | "we took the method for inducing levels in OrthoXML from pyHam." |

Rule 4 (missing "that"): no clear instance found in the three texts; the closest case is Methods "Subagents ... shared files through the project directory and git", which reads well.
