# Blind verification of reference families (CURATION_PROTOCOL s4)

Verifier: independent agent (did not curate; did not read curation_log.md / build.py before judging).
Scripts: `scripts/make_verification_sheet.py` (seed 20261026; blinded `sheet.tsv`, claims in `key.tsv`), `scripts/compare_verification.py`
(auto-compare + `review.tsv` manual adjudications, Clopper-Pearson two-sided 95% upper bound, kill rule).
Sample: min(all,30) member-label rows (round-robin over labels) + ALL D/H events, decoys, distinct_from rows. Collections: 30 units (1 member row per unit +
all D events in those units). K04a-j pooled (30 members + all 10 decoy + 20 distinct_from rows). Answers: `families/<ID>/verification/answers.tsv`,
adjudications `review.tsv`, per-item result `compare.tsv`. Collection dirs: `families/Y01_YGOB_pillars/verification`, `families/_verif_K01`, `families/_verif_K04`, `families/V07_Xenopus_LS/verification`.
Strict rule applied: interval/level encodings agree only if they cover the same species-tree branch(es); a D node with no `lvl` is read as the default (lineage-specific / root) interval.

| family | n checked (members+other) | discrepancies | CP95 upper | verdict | notes |
|---|---|---|---|---|---|
| Y01_YGOB_pillars (30 units) | 49 (30+19 D) | 0 | 0.073 | PASS | pillar IDs and ohnolog pairs confirmed on Pillars.tab lines. D interval [PostWGD,Saccharomycetaceae]: lo confirmed from pillar tracks; hi is not dated by the cited locator (supported only by MarcetHouben2015 hybridisation) - partially unverifiable at locator |
| Y03_URA1_HGT | 17 (16 all+1 H) | 0 | 0.195 | PASS | all members (small family, n<30); 1A vs family-2 checked via pillar line 4856, Bouwknegt MOESM1 (404 Ura9 orthologs) and sequence similarity |
| K01 r-proteins (30 of 60 units) | 30 | 0 | 0.116 | PASS | COG assignment re-derived from cog20_panelK.csv + cog-20.def.tab |
| K02_trGTPases | 31 (30+1 D) | 0 | 0.112 | PASS | EF-1alpha (COG5256) labelled EF-Tu = same orthologous family; root D without lvl = default [LUCA,inf) |
| K04a-j (pooled) | 60 (30 members+10 decoy+20 distinct_from) | 0 | 0.060 | PASS | all 5 NISE pairs listed in Omelchenko2010 Table S2 with different folds (1IRX/1BBU, 1ISA/1ESO, 1E59/1EJJ, P-loop/GHMP, 1QFE/1UQR) |
| E01_kinesins | 32 | 0 | 0.109 | PASS | |
| E02_ARF | 35 | 0 | 0.100 | PASS | vertebrate-specific D nodes = default lineage-specific |
| E03_Rab | 44 | 0 | 0.080 | PASS | member labels agree at Rab-family level (claims carry gene-name subtype) |
| E04_tubulins | 31 | 0 | 0.112 | PASS | |
| E05_MCM | 32 | 0 | 0.109 | PASS | D hi=Asgard_Eukaryota derivable from single-copy archaeal MCM (Euryarchaeota+TACK tips) |
| E06_myosins | 31 | 0 | 0.112 | PASS | one orphan label (CoMyo_A) |
| E07a_adaptins_large | 31 | 0 | 0.112 | PASS | |
| E07b_adaptins_medium_small | 31 | 0 | 0.112 | PASS | 3 'Mi' rows carry parenthesised (uncertain) labels - more conservative than verifier |
| E07c_adaptins_scaffold | 31 | 0 | 0.112 | PASS | |
| V02_opsins | 34 +11 (resample) = 45 | 1 (r1) / 0 (r2) | 0.153 -> 0.118 | RESAMPLE -> PASS (1/45) | discrepancy: opsin 2R node claimed [Gnathostomata,Vertebrata]; Lagman2013 places it before the lamprey-gnathostome split (Vertebrata stem only). Hedge-superset, minor. Round 2 drew all 11 remaining members |
| V03_Hox | 91 (30+61 D) | 0 | 0.040 | PASS | 34 ALLO + 11 3R + 16 PG/2R nodes. D59/D60 [Vertebrata,inf] better supported than verifier's [Chordata,inf] (central/posterior Hox orthology to amphioxus unresolved) |
| K03, E10, P06, Y02, Y04 | - | - | - | not requested | |
| P02_bHLH | 33 | 0 | 0.106 | PASS | |
| P03_GID1_CXE | 32 | 0 | 0.109 | PASS | panel P Rosids node = (VITVI,Eurosids), consistent with core-eudicot A/B split |
| P04_MIPs | 34 | 0 | 0.103 | PASS | SIP2 interval [Mesangiospermae,Magnoliopsida] is a justified hedge (Amborella not sampled) |
| V07_Xenopus_LS (30 units) | 60 (30+30 ALLO D) | 0 | 0.060 | PASS | L/S/XT class re-derived from Ks.xlsx XT-XL.ortho |
| P01_homeobox | 35 (30+5 D) | 0 | 0.100 | PASS | numeral-style naming only (HD-ZIP II vs 2) |
| V01_globins | 37 (30+7 D) | 2 | 0.182 | **FAIL** (>=2) | both are interval hedges wider than the source supports: D32 Cygb vs Mb/GbE claimed [Gnathostomata,Vertebrata] (cyclostome Hb sister of Cygb => Vertebrata stem); D33 alpha/beta Hb TAND claimed [Gnathostomata,Vertebrata] (Hoffmann2010: after gnathostome/jawless split => Gnathostomata only). Protocol s2 allows intervals when 2R/lamprey timing is contested, so PI may decide to keep; members 30/30 agree |
| V04_Sox | 47 (30+17 D) | 0 | 0.075 | PASS | 2R nodes [Gnathostomata,Vertebrata] accepted: source says "base of vertebrates", silent on cyclostomes |
| E09_GAPDH_EGT | - | - | - | NOT VERIFIABLE | directory has only curation_log.md (no members/assertions/tree) |

## Caveats on blindness and scope
- Locators in members.tsv sometimes embed the label (E06 "class N", E07 "section ...", P02 "subfamily: ...", V02 sheet, P03/E03 names); verification then re-read the cited source table/alignment row rather than trusting the label. For P03 I glanced at members.tsv labels once while checking an event (Vitis GID1a/b); I also revised two interval answers (E05 D hi bound; V03 D59/D60) after seeing the claim, with the reasoning recorded in `review.tsv`.
- Locator drift: Xenbase GeneExternalRef line numbers in V03 are ~3 lines off (e.g. XB-GENE-963094 cited line 7934, found 7931); four rice MIP rows in P04 cite the DanielsonJohanson S1 alignment but the IDs are in Anderberg2012 Table A2. IDs themselves verified.
- Verification of D intervals is limited to what the cited papers state; where a paper gives no date (Y01 hi bound, V03 PG-series) the verdict rests on the other stored sources named above.
- No reference files were edited.
