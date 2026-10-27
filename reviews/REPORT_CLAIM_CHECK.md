# REPORT.md claim check (independent)

Tests: `pytest -q tests` gives 95 passed (matches REPORT). Verdicts: OK / WRONG / UNSUPPORTED / IMPRECISE.

## Problems

| Claim | Location | Source | Verdict | Correct value / note |
|---|---|---|---|---|
| "5 species panels (98 species)" | Summary Data | panels/*/species.tsv | WRONG | 16+21+34+17+14 = 102 panel-species slots; 90 unique taxa (prokaryotes shared by E and K). Neither equals 98. |
| "22 headline-eligible" of 32 families | Summary References | results/catalogue/families.md | WRONG | 16 rows are headline "y" (13 scored items plus 3 controls K01, K02, E05); 16 are "n". 22 not derivable. |
| "82 text-asserted (T1)" events | Summary References | families.tsv D_T1 | WRONG | Sum of D_T1 = 85. T2 = 2,138, HGT = 1 and ~31.8k genes (31,767) are OK. |
| "21/22 families passed blind re-derivation (<=2/30 discrepancies)" | Summary; Sec 1.3 | families/VERIFICATION.md | WRONG | 23 verification groups: 22 PASS, 1 FAIL (V01, 2 discrepancies). V02 passed only after resample (1/45). The kill rule fails at >=2, so passes have <=1 discrepancy, not <=2. |
| Section 2.2 lists K04a-j under "Headline" | 2.2 | families.md | WRONG | K04a-j are headline "n: no events". Sec 5.1 correctly says they only test contamination. |
| "V subset ... V07 sampled to 300 units" | Limitation 6 | results/V_subset/per_family.tsv | UNSUPPORTED | 333 XLAEVIS units are scored in V subset, as the REPORT itself says in 5.3. No source gives 300. |
| P02 bHLH "precision 0.95 but recall ~0.01" | 5.4 | P_subset/per_family.tsv | WRONG | membership precision 0.959, recall 0.025 (F1 0.049; error_anatomy header says "~0.02"). |
| OrthoFinder S-pair over-split "22.7%" | 5.1 note | error_anatomy.md (b) | WRONG | 20.4% (FastOMA 66.7% is OK). |
| Subset-vs-full pair merge-level identity "0.97-0.996" | 5.5 | equivalence.md | WRONG | Range is 0.9613 (K OrthoFinder) to 0.9956 (Y of_flat). |
| "Hardest headline families: P02, E02, E01" | 5.4 | difficulty.md | IMPRECISE | Rank order is E02 (1.77), E01 (1.51), P02 (1.46). P02's difficulty rests on 3 skill values from 1 endpoint. |
| "decisions D0-D10" | header | RESEARCH_TREE.md | IMPRECISE | D11 exists (T+18h). |
| V proteins "~400k" | 2.1 | panels/V/species.tsv | IMPRECISE | 447,758. |
| V07 "498 units" | Summary; 2.2 | catalogue says 500; curation_log says 500->498 | IMPRECISE | The frozen catalogue says 500 units and 1,500 members; 498 comes from a later curation log. State which is frozen. |
| "E/P/V results are therefore reported in subset mode" | 5.5 | DESIGN 7 | IMPRECISE | DESIGN says that if equivalence fails, results are reported conditional on subset mode. K event recovery was not equivalent (n=1), and K OF merge-acc CI reaches +0.044 (margin 0.05). Equivalence on Y/K does not license E/V/P; the REPORT only half-says this. |

## Over-claims

| Claim | Verdict | Note |
|---|---|---|
| "OrthoFinder beats FastOMA" (Summary i) | IMPRECISE | Holm p=0.035 is one pooled test over 13 non-independent items (Y has 2, E 5, V 3, P 3). The gain is driven by E and V; on Y and P FastOMA is equal or better in several endpoints. It is not a general statement about the tools. |
| "(iii) ... opposite ways" | IMPRECISE | Rests on the Y/V WGD pools only (separation 0.36 vs 0.67, cohesion 0.80 vs 0.61). The V Hox WGD events show FastOMA failing on separation (0.73-0.82), not only cohesion. |
| "(vi) DendroBLAST scores best" | IMPRECISE | True only for the pooled merge accuracy and event recovery. On V it is worse than OF default (merge 0.928 vs 0.960, events 0.759 vs 0.907), and there are no CIs or tests. |
| "(vii) Harder families tend to show lower agreement" | IMPRECISE | rho=-0.41, p=0.19, n=12. The caveat is given, but "tend" suggests more than a non-significant trend supports. The earlier n=20 rho=-0.58 is in RESEARCH_TREE only. |
| "nonzero deltas are attributable to the perturbation" | OK with caveat | Holds for OrthoFinder only (bit-identical replicate). FastOMA was not replicated, as the REPORT says. |
| "FastOMA over-merges LECA paralogs (24.4%)" | IMPRECISE | 24.4% is over-merged share for all E D-node pairs (headline subset), 24.0% at Eukaryota. It is not LECA-specific, and "OF 0.1%" is the same pool. |

## Verified OK (spot list)

- Main table (Y/E/P/V rows, baselines incl. one_big_hog 0.838, membership 0.180): all match main_table.md.
- Confirmatory: 13/13/11/5 items; +0.054/+0.134/+0.285/+0.315; Holm 0.035/0.058/0.12/0.12 (confirmatory.tsv).
- Pooled dendroblast/default/FastOMA: merge 0.886/0.846/0.712 (n=13); events 0.617/0.479/0.194 (n=11).
- Hidden-paralogy trap: 0.041 [0.000, 0.105] and 0.122 [0.038, 0.224], n=49; same-track 0.910/0.865, n=1,189. WGD events: 243; OF 0.222, FastOMA 0.235; separation 0.362; cohesion 0.606. YGOB D-S: -0.683/-0.206, n=532.
- Hox: 0.77/0.82 vs 0.08/0.09 vs flat 0/0 (n=13/11).
- X. laevis: 0.988 / 0.429, 55.6% over-split, 333 units; FastOMA -0.26 [-0.29, -0.23]; OF -0.005.
- Zebrafish 3R 91% over-split; V FastOMA over-merge <=1.7%.
- Robustness table and the 0.005 bound; 2,615 fragmented genes; 4,156 unit-metrics; 8-13 units touched.
- Difficulty rho=-0.41, p=0.19, n=12.
- Equivalence: Y -0.0003 [-0.0008, +0.0002], n=1,976; Y subset 54%, K 6%.
- Catalogue phenomenon counts (WGD 5 ... controls 3); event depth 2,136/38/5/13.
- YGOB 202/1,592 testable (ERRATA E1); mutation 10/10 killed (mutant c killed only incidentally in the original suite, later given a direct test); 95 tests; RECEIVED: 49 distinct ranks delivered; ~11 h outages (3 in RESEARCH_TREE); AUDIT 0 violations; species counts per panel (16/21/34/17/14) and protein counts (85k/68k/340k/418k); subset sizes 25.9k/25.4k/70k.
- Unverifiable here: "T+19h", "~25 primary sources", "<=98% identity" thresholds (protocol-level).
