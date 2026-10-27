# Abandoned / non-built families and units

Source: families/*/curation_log.md (directories lacking family.yaml/collection.yaml), plus the K01 unit table in the shared RPROT log.

| id | reason (from curation log) |
|---|---|
| E09_GAPDH_EGT | - The quotable EGT (cyanobacterial Gap2 -> Archaeplastida GapA/B, H node donor=SYNY3) has source-listed members in only 3 panel species (SYNY3, ARATH, CHLRE). Red algal GapAB is stated (Figge 1999) but no source lists a C. merolae sequence; THAPS uses GapC1 (plastid-retargeted cytosolic copy), not GapA. - Extending to cytosolic GapC (HUMAN, YEAST, DICDI, ...) would add only unasserted (U) members: the GapC mitochondrial-donor hypothesis has explicit exceptions (Euglena/trypanosome GapC), no panel donor member is source-listed (no BRADU/ECOLI gap1 in any source), and eukaryotic GapC copy number |
| E10_kinome | - Table S1 (sequences + group/family/subfamily/pseudogene) and S4 (pseudogenes) list HUMAN kinases only. - Table S6 gives per-subfamily domain COUNTS for yeast, worm, fly, human - no yeast/fly gene identifiers or sequences. - TIBS paper has no per-gene lists (it points to www.kinase.com, offline; the Wayback KinBase files in literature/raw/B1/kinbase/ are again human-only + the same count table). - => members mappable in panel E: HUMAN only (1 species). Cross-species orthology cannot be asserted; protocol requires >=4 panel species. |
| K01_RPROT_L28 | excluded_multicopy: MYCTU: 3 COG 2020 members; MYCTU: Yutin Table S3 count 3 (2nd paralog) |
| K01_RPROT_L31 | excluded_multicopy: ECOLI: 2 COG 2020 members; ECOLI: Yutin Table S3 count 2 (2nd paralog); PSEAE: 2 COG 2020 members; PSEAE: Yutin Table S3 count 2 (2nd paralog); BACSU: 2 COG 2020 members; BACSU: Yutin Table S3 count 2 (2nd paralog) |
| K01_RPROT_L33 | excluded_multicopy: BACSU: 2 COG 2020 members; BACSU: Yutin Table S3 count 3 (2nd paralog); MYCTU: 2 COG 2020 members; MYCTU: Yutin Table S3 count 2 (2nd paralog) |
| K01_RPROT_L36 | excluded_multicopy: ECOLI: 2 COG 2020 members; ECOLI: Yutin Table S3 count 2 (2nd paralog); PSEAE: 2 COG 2020 members; PSEAE: Yutin Table S3 count 2 (2nd paralog) |
| K01_RPROT_L7L12 | excluded_by_curator: Table S1 assigns COG2058 (archaeal L12E/P1/P2 family) to bacterial L7/L12; in COG 2020 COG2058 has no bacterial member in the panel genomes (Yutin Table S3 lists L7L12 in every bacterium), so the COG cannot serve as the membership source for this universal unit |
| K01_RPROT_L7ae | excluded_multicopy: BACSU: 2 COG 2020 members; BACSU: Yutin Table S3 count 2 (2nd paralog) |
| K01_RPROT_S1 | excluded_by_curator: Yutin 2012 excluded S1p from the study (varied domain architecture, ubiquitous S1-like domains) |
| K01_RPROT_S10 | excluded_multicopy: HALSA: 2 COG 2020 members; HALSA: Yutin Table S3 count 2 (2nd paralog) |
| K01_RPROT_S11 | excluded_multicopy: METAC: 2 COG 2020 members |
| K01_RPROT_S12 | excluded_multicopy: AQUAE: 2 COG 2020 members; AQUAE: Yutin Table S3 count 2 (2nd paralog) |
| K01_RPROT_S14 | excluded_multicopy: BACSU: 2 COG 2020 members; BACSU: Yutin Table S3 count 2 (2nd paralog); MYCTU: 2 COG 2020 members; MYCTU: Yutin Table S3 count 2 (2nd paralog); DEIRA: 2 COG 2020 members |
| K01_RPROT_S17 | excluded_multicopy: BACTN: 2 COG 2020 members; BACTN: Yutin Table S3 count 2 (2nd paralog); METAC: 2 COG 2020 members; METAC: Yutin Table S3 count 2 (2nd paralog) |
| K01_RPROT_S18 | excluded_multicopy: MYCTU: 2 COG 2020 members; MYCTU: Yutin Table S3 count 2 (2nd paralog) |
| K01_RPROT_S21 | excluded_multicopy: BRADU: Yutin Table S3 count 2 (2nd paralog); GEOSL: Yutin Table S3 count 2 (2nd paralog) |
| K01_RPROT_S7 | excluded_multicopy: AQUAE: 2 COG 2020 members; AQUAE: Yutin Table S3 count 2 (2nd paralog) |
| K03_nitrogenase | Status: **abandoned** (fewer than 4 panel species have members). Curator: K-panel curation agent, 2026-10-26, ~5 min. |
| P06_MADS | Status: **abandoned** (protocol/brief feasibility rule: < 4 panel species with labelled, mappable IDs). ~15 min. |
| V05_teleost3R | Status: **not done** |
| V06_salmonid4R | Status: **not done** |
| V08_chemokines | Status: **not done** |
| V09_KRAB_ZNF | Status: **not done** |
| Y02_MH2015_ohnologs | Status: **abandoned as a standalone family** (kept as a documented cross-check of Y01). Time ~30 min. |
| Y04_GAL_Spombe_HGT | 1. HGT not testable in panel Y: the donor is a Candida (CUG-clade) yeast, represented in panel Y by CANAL, and the only recipient is SCHPO. With donor = CANAL the H-node merge level is lca(CANAL, SCHPO) = Ascomycota, the panel root, which equals the species LCA of every SCHPO-vs-Saccharomycotina pair. The xenolog relation is therefore indistinguishable from orthology in this panel (no xenolog-separation signal); native S. pombe GAL1/GAL7 copies no longer exist to be separated from. Only the native unclustered GAL10 paralog (uge1) would give a contrast, and no source places it. 2. Member identi |

## Built but not headline-eligible

- E06_myosins: mapped 73% < 80%
- E07a_adaptins_large: mapped 72% < 80%
- E07b_adaptins_medium_small: mapped 72% < 80%
- K04a_LysRS_I: no events
- K04b_LysRS_II: no events
- K04c_SOD_FeMn: no events
- K04d_SOD_CuZn: no events
- K04e_PGM_BPGdep: no events
- K04f_PGM_BPGindep: no events
- K04g_SK_Ploop: no events
- K04h_SK_GHMP: no events
- K04i_DHQase_I: no events
- K04j_DHQase_II: no events
- P01_homeobox: mapped 62% < 80%
- V01_globins: mapped 77% < 80%; verification FAIL
- V04_Sox: mapped 46% < 80%
- K01_RPROT units L21, L29, L32, S16: <80 % uniquely mapped (built, non-headline; per shared log)

## V07_Xenopus_LS: candidate units not retained

- 4171 candidates: not_sampled (passes all checks)
- 3119 candidates: not all three genes uniquely mapped
- 819 candidates: names do not form one .L and one .S homeolog
- 36 candidates: L and S proteins identical (excluded, protocol rule 6)
- 5 candidates: homology QC failed (S:nohit)
- 2 candidates: protein or unit id shared with another source row
- 1 candidates: homology QC failed (S:A0A8J1IYJ7 67.2%id 100%cov)
- 1 candidates: homology QC failed (S:A0A8J0PLH8 65.6%id 100%cov)
- 1 candidates: homology QC failed (S:F7A2E2 67.8%id 100%cov)
- 1 candidates: homology QC failed (L:A0A8J0SS71 93.8%id 50%cov)
- 1 candidates: homology QC failed (L:A0A8J0R0Z0 92.9%id 41%cov; S:A0A8J0R0Z0 88.4%id 41%cov)
- 1 candidates: homology QC failed (S:A0A8J0PK63 61.2%id 87%cov)
- 1 candidates: homology QC failed (L:A0A8J1IX76 65.9%id 98%cov; S:A0A8J1IX76 62.4%id 98%cov)
- 1 candidates: homology QC failed (L:A0A8J1J016 63.2%id 100%cov; S:A0A8J1J016 59.6%id 100%cov)
- 1 candidates: homology QC failed (L:A0A6I8QF03 65.8%id 100%cov; S:A0A6I8QF03 62.2%id 100%cov)
- 1 candidates: homology QC failed (L:A0A6I8QVQ0 64.8%id 100%cov)
- 1 candidates: homology QC failed (L:A4QNR3 59.6%id 98%cov; S:A4QNR3 62.4%id 98%cov)
- 1 candidates: homology QC failed (L:Q6DF48 53.3%id 98%cov; S:Q6DF48 57.9%id 100%cov)
- 1 candidates: homology QC failed (L:A0A6I8Q482 67.3%id 100%cov; S:A0A6I8Q482 65.1%id 100%cov)
- 1 candidates: homology QC failed (L:Q6DIX6 62.7%id 96%cov; S:Q6DIX6 63.5%id 96%cov)
- 1 candidates: homology QC failed (L:A0A8J1IQB1 68.7%id 98%cov; S:A0A8J1IQB1 64.1%id 100%cov)
- 1 candidates: homology QC failed (L:F6VBD7 69.9%id 100%cov)
- 1 candidates: homology QC failed (L:A0A803K8J4 64.6%id 99%cov; S:A0A803K8J4 66.1%id 100%cov)
- 1 candidates: homology QC failed (L:B1H3G4 68.0%id 97%cov)
- 1 candidates: homology QC failed (L:A0A803K3Q1 65.8%id 100%cov)
- 1 candidates: homology QC failed (L:A0A8J0PKH7 86.9%id 49%cov)
- 1 candidates: homology QC failed (S:A0A6I8PZ44 56.8%id 87%cov)
- 1 candidates: homology QC failed (S:A0A8J1JIT8 65.9%id 100%cov)
- 1 candidates: homology QC failed (S:A0A9Y6B502 68.9%id 100%cov)
- 1 candidates: homology QC failed (L:A0A8J0QBJ6 57.3%id 75%cov; S:A0A8J0QBJ6 62.4%id 74%cov)
- 1 candidates: homology QC failed (L:A0A8J0PLY1 69.1%id 93%cov)
- 1 candidates: homology QC failed (L:A0A8J0QPA1 35.8%id 80%cov)
- 1 candidates: homology QC failed (S:F6R5I8 65.8%id 100%cov)
- 1 candidates: homology QC failed (S:A0A803JMX5 69.4%id 97%cov)
- 1 candidates: homology QC failed (L:A0A8J0PMK1 95.3%id 19%cov; S:A0A8J0PMK1 95.0%id 19%cov)
- 1 candidates: homology QC failed (L:A0A8J0QKA9 94.5%id 47%cov; S:A0A8J0QKA9 94.5%id 47%cov)
- 1 candidates: homology QC failed (L:A0A8J0QEZ5 93.0%id 18%cov; S:A0A8J0QEZ5 93.7%id 18%cov)
- 1 candidates: homology QC failed (L:F6QGI3 94.1%id 18%cov; S:F6QGI3 96.0%id 18%cov)
- 1 candidates: homology QC failed (L:F7D2U4 95.8%id 24%cov)
- 1 candidates: homology QC failed (L:A0A8J0Q6H2 97.7%id 25%cov; S:A0A8J0Q6H2 97.1%id 25%cov)
- 1 candidates: homology QC failed (L:A0A8J0QA96 93.9%id 24%cov; S:A0A8J0QA96 100.0%id 6%cov)
- 1 candidates: homology QC failed (L:A0A8J0T247 52.6%id 48%cov)
- 1 candidates: homology QC failed (L:Q6P635 63.4%id 100%cov)
- 1 candidates: homology QC failed (L:A0A8J1IZC8 43.8%id 62%cov; S:nohit)
- 1 candidates: homology QC failed (S:Q28H62 93.2%id 39%cov)
- 1 candidates: homology QC failed (L:A0A8J0PKG9 68.7%id 100%cov)
- 1 candidates: homology QC failed (L:A0A8J0SW01 30.9%id 24%cov)
- 1 candidates: homology QC failed (L:A0A8J0QKH5 40.3%id 98%cov; S:A0A8J0QKH5 66.4%id 96%cov)
- 1 candidates: homology QC failed (L:A0A6I8PUM1 53.9%id 93%cov)
- 1 candidates: homology QC failed (S:A0A803KJ62 67.9%id 99%cov)
- 1 candidates: homology QC failed (S:A0A6I8S328 50.4%id 89%cov)
- 1 candidates: homology QC failed (S:A0A6I8RK85 56.2%id 99%cov)
- 1 candidates: homology QC failed (L:A0A8J0Q644 61.6%id 100%cov; S:A0A8J0Q644 60.5%id 99%cov)
- 1 candidates: homology QC failed (S:A0A8J0PKL4 69.9%id 94%cov)
- 1 candidates: homology QC failed (L:A0A8J0SVC5 50.4%id 100%cov; S:A0A8J0SVC5 57.7%id 100%cov)
- 1 candidates: homology QC failed (L:Q28F43 65.5%id 100%cov; S:Q28F43 67.5%id 100%cov)
- 1 candidates: homology QC failed (L:Q5RJV1 96.4%id 50%cov)
- 1 candidates: homology QC failed (S:A0A8J0PLE2 63.9%id 91%cov)
- 1 candidates: homology QC failed (L:A0A8J0PHX1 54.0%id 98%cov; S:A0A8J0PHX1 58.9%id 90%cov)
- 1 candidates: homology QC failed (L:A0A803JL01 69.9%id 94%cov)
- 1 candidates: homology QC failed (L:F6Q943 67.5%id 99%cov)
- 1 candidates: homology QC failed (S:F7B6X6 69.3%id 100%cov)
- 1 candidates: homology QC failed (L:A0A8J0QYD1 67.6%id 100%cov; S:A0A8J0QYD1 66.8%id 96%cov)
- 1 candidates: homology QC failed (L:A0A8J0QU05 68.0%id 97%cov; S:A0A8J0QU05 69.4%id 98%cov)
- 1 candidates: homology QC failed (L:A0A8J1J1D5 65.6%id 61%cov; S:A0A8J1J1D5 69.7%id 72%cov)
- 1 candidates: homology QC failed (L:A0A8J1JAM6 62.3%id 91%cov; S:A0A8J1JAM6 65.3%id 98%cov)
- 1 candidates: homology QC failed (L:A0A8J0R732 69.0%id 100%cov)
- 1 candidates: homology QC failed (S:A0A8J0QB16 68.7%id 97%cov)
- 1 candidates: homology QC failed (S:A0A803JB90 67.5%id 96%cov)
- 1 candidates: homology QC failed (L:Q6DF67 69.6%id 50%cov)
- 1 candidates: homology QC failed (L:A0A8J0PCW7 60.6%id 100%cov)
- 1 candidates: homology QC failed (L:A0A8J0QH50 63.5%id 100%cov; S:A0A8J0QH50 59.0%id 99%cov)
- 1 candidates: homology QC failed (L:A0A8J1J9K7 63.9%id 100%cov; S:A0A8J1J9K7 66.1%id 100%cov)
- 1 candidates: homology QC failed (L:A0A8J0Q9X3 65.7%id 100%cov; S:A0A8J0Q9X3 66.4%id 100%cov)
- 1 candidates: homology QC failed (L:A0A6I8RXC4 68.7%id 99%cov; S:A0A6I8RXC4 68.6%id 99%cov)
- 1 candidates: homology QC failed (L:A0A6I8PZT5 59.9%id 89%cov; S:A0A6I8PZT5 62.5%id 90%cov)
- 1 candidates: homology QC failed (L:A0A6I8QQ56 54.1%id 90%cov)
- 1 candidates: homology QC failed (S:A0A8J1JLH7 54.2%id 72%cov)
- 1 candidates: homology QC failed (L:A0A803K1Y5 59.8%id 98%cov; S:A0A803K1Y5 65.8%id 98%cov)
- 1 candidates: homology QC failed (L:A0A6I8Q579 64.1%id 69%cov; S:A0A6I8Q579 64.3%id 100%cov)
- 1 candidates: homology QC failed (S:A0A8J0QW67 30.3%id 96%cov)
- 1 candidates: homology QC failed (S:A0A6I8PNB0 69.5%id 100%cov)
- 1 candidates: homology QC failed (L:A0A6I8QZZ5 69.6%id 100%cov; S:A0A6I8QZZ5 68.9%id 100%cov)
- 1 candidates: homology QC failed (S:A0A8J0QRK7 64.3%id 99%cov)
- 1 candidates: homology QC failed (L:A0A6I8SH40 66.4%id 100%cov)
- 1 candidates: homology QC failed (S:A0A6I8RHE7 45.5%id 94%cov)
- 1 candidates: homology QC failed (S:A0A8J0SJX7 65.8%id 100%cov)
- 1 candidates: homology QC failed (S:A0A6I8RX68 34.5%id 99%cov)
- 1 candidates: homology QC failed (L:A0A1B8XSL0 85.1%id 22%cov; S:A0A1B8XSL0 81.2%id 22%cov)
- 1 candidates: homology QC failed (L:A0A8J0QP93 61.1%id 100%cov; S:A0A8J0QP93 62.6%id 100%cov)
- 1 candidates: homology QC failed (L:A0A8J1JDH5 51.3%id 100%cov; S:A0A8J1JDH5 66.2%id 100%cov)
- 1 candidates: homology QC failed (S:A0A8J0T5Y1 61.5%id 97%cov)
- 1 candidates: homology QC failed (L:A0A803JHS3 69.1%id 100%cov; S:A0A803JHS3 66.5%id 100%cov)
- 1 candidates: homology QC failed (L:A0A6I8SI97 65.8%id 98%cov; S:A0A6I8SI97 67.9%id 96%cov)
- 1 candidates: homology QC failed (L:A0A803JNP0 69.6%id 100%cov; S:A0A803JNP0 67.9%id 100%cov)
- 1 candidates: homology QC failed (L:A0A8J1J4H7 62.8%id 100%cov; S:A0A8J1J4H7 44.6%id 99%cov)
- 1 candidates: homology QC failed (L:A0A8J0SRT2 46.5%id 79%cov; S:A0A8J0SRT2 39.9%id 99%cov)
- 1 candidates: homology QC failed (S:A0A8J0R854 62.9%id 91%cov)
- 1 candidates: homology QC failed (L:A0A8J1JWC5 64.2%id 100%cov; S:A0A8J1JWC5 65.0%id 100%cov)
- 1 candidates: homology QC failed (L:A0A8J0R499 69.1%id 100%cov)
- 1 candidates: homology QC failed (S:A0A6I8RQM0 51.9%id 95%cov)
- 1 candidates: homology QC failed (L:A0A803J685 66.0%id 38%cov; S:nohit)
- 1 candidates: homology QC failed (L:A0A8J0R4E3 64.5%id 100%cov; S:A0A8J0R4E3 66.8%id 100%cov)
- 1 candidates: homology QC failed (L:A0A8J0SQN3 76.7%id 50%cov)
- 1 candidates: homology QC failed (S:A0A8J1J261 63.1%id 100%cov)
- 1 candidates: homology QC failed (L:A0A8J0QUJ9 58.9%id 100%cov; S:A0A8J0QUJ9 60.0%id 97%cov)
- 1 candidates: homology QC failed (L:A0A8J1ITZ9 69.4%id 99%cov; S:A0A8J1ITZ9 63.1%id 100%cov)
- 1 candidates: homology QC failed (L:B0S4Q5 60.2%id 100%cov; S:B0S4Q5 60.5%id 100%cov)
- 1 candidates: homology QC failed (L:A0A8J0PDS2 40.0%id 95%cov)
- 1 candidates: homology QC failed (L:A0A6I8QND8 91.5%id 38%cov; S:A0A6I8QND8 93.5%id 35%cov)
- 1 candidates: homology QC failed (L:A0A803JEA0 60.2%id 100%cov; S:A0A803JEA0 62.8%id 100%cov)
- 1 candidates: homology QC failed (L:A0A8J0PK02 67.3%id 88%cov)
- 1 candidates: homology QC failed (L:nohit)
- 1 candidates: homology QC failed (S:A0A6I8Q2C1 69.0%id 55%cov)
- 1 candidates: homology QC failed (L:A0A8J0R0U7 96.5%id 43%cov)
- 1 candidates: homology QC failed (L:A0A6I8R3U5 66.8%id 100%cov; S:A0A6I8R3U5 65.4%id 100%cov)
- 1 candidates: homology QC failed (S:A0A8J0QAV6 44.7%id 100%cov)
- 1 candidates: homology QC failed (S:A0A6I8SHL5 60.0%id 100%cov)
- 1 candidates: homology QC failed (L:A0A6I8SLR8 56.4%id 100%cov; S:A0A6I8SLR8 59.3%id 100%cov)
- 1 candidates: homology QC failed (L:A0A8J0Q9R9 82.8%id 47%cov)
- 1 candidates: homology QC failed (L:A0A8J1JGD0 50.1%id 100%cov; S:A0A8J1JGD0 76.3%id 42%cov)
- 1 candidates: homology QC failed (S:A0A6I8RRV5 84.4%id 29%cov)
- 1 candidates: homology QC failed (L:A0A8J0PHL6 96.1%id 34%cov; S:A0A8J0PHL6 88.2%id 34%cov)
- 1 candidates: homology QC failed (L:A0A6I8PNP7 64.4%id 100%cov)
- 1 candidates: homology QC failed (L:A0A8J0Q7M4 91.9%id 36%cov; S:A0A8J0Q7M4 63.7%id 36%cov)
- 1 candidates: homology QC failed (L:A0A8J0Q9V4 66.3%id 87%cov; S:A0A8J0Q9V4 66.1%id 82%cov)
- 1 candidates: homology QC failed (L:A0A8J0PIV9 91.2%id 14%cov; S:A0A8J0PIV9 89.4%id 15%cov)
- 1 candidates: homology QC failed (L:A0A803KIA3 87.5%id 23%cov; S:A0A803KIA3 86.4%id 48%cov)
- 1 candidates: homology QC failed (L:A0A803JF35 64.2%id 100%cov; S:A0A803JF35 53.1%id 100%cov)
- 1 candidates: homology QC failed (S:A0A8J0QQ54 65.0%id 100%cov)
- 1 candidates: homology QC failed (S:A0A6I8QRK3 63.3%id 51%cov)
- 1 candidates: homology QC failed (L:A0A6I8RX42 64.2%id 100%cov)
- 1 candidates: homology QC failed (L:A0A8J0PMD4 95.4%id 50%cov; S:A0A8J0PMD4 96.3%id 50%cov)
- 1 candidates: homology QC failed (L:A0A8J0PIZ1 87.9%id 39%cov; S:A0A8J0PIZ1 90.3%id 39%cov)
- 1 candidates: homology QC failed (L:A0A8J0PD29 68.7%id 52%cov)
- 1 candidates: homology QC failed (S:A0A8J0PKZ3 58.2%id 100%cov)
- 1 candidates: homology QC failed (L:Q0V9A3 66.1%id 98%cov; S:Q0V9A3 65.8%id 98%cov)
- 1 candidates: homology QC failed (L:A0A6I8SF91 68.2%id 100%cov; S:A0A6I8SF91 65.5%id 100%cov)
- 1 candidates: homology QC failed (S:A0A6I8RQD2 66.3%id 97%cov)
- 1 candidates: homology QC failed (L:Q2AB69 68.7%id 100%cov; S:Q2AB69 67.3%id 100%cov)
- 1 candidates: homology QC failed (L:A0A6I8RNY1 62.2%id 100%cov; S:A0A6I8RNY1 61.2%id 100%cov)
- 1 candidates: homology QC failed (L:A4IHD1 67.0%id 94%cov)
- 1 candidates: homology QC failed (L:A0A6I8R7K5 57.2%id 100%cov; S:A0A6I8R7K5 62.4%id 100%cov)
- 1 candidates: homology QC failed (L:F7DPT4 69.6%id 100%cov; S:F7DPT4 66.7%id 100%cov)
- 1 candidates: homology QC failed (L:A0A6I8PW52 69.0%id 100%cov; S:A0A6I8PW52 64.6%id 100%cov)
- 1 candidates: homology QC failed (S:A0A6I8PS95 61.8%id 100%cov)
- 1 candidates: homology QC failed (L:A0A8J0SWU3 61.2%id 85%cov)
- 1 candidates: homology QC failed (L:A0A803K876 59.7%id 100%cov)
- 1 candidates: homology QC failed (L:A0A8J0T1Z7 58.9%id 100%cov; S:A0A8J0T1Z7 56.1%id 100%cov)
- 1 candidates: homology QC failed (L:F7EI49 56.8%id 100%cov; S:F7EI49 57.1%id 99%cov)
- 1 candidates: homology QC failed (S:Q5I030 68.6%id 98%cov)
- 1 candidates: homology QC failed (L:A0A8J0Q9F8 35.2%id 55%cov)
- 1 candidates: homology QC failed (L:A0A1B8XZT7 64.5%id 66%cov; S:A0A1B8XZT7 61.9%id 66%cov)
- 1 candidates: homology QC failed (L:A0A8J0PCM6 91.3%id 41%cov; S:A0A8J0PCM6 88.9%id 29%cov)
- 1 candidates: homology QC failed (L:A0A8J0PM06 89.7%id 45%cov; S:A0A8J0PM06 94.9%id 44%cov)
- 1 candidates: homology QC failed (S:A0A8J1INN5 69.1%id 100%cov)
- 1 candidates: homology QC failed (S:A0A803K5P9 62.9%id 98%cov)
- 1 candidates: homology QC failed (L:A0A6I8SKN1 61.6%id 100%cov; S:A0A6I8SKN1 68.6%id 100%cov)
- 1 candidates: homology QC failed (L:A0A9Y6DVU6 67.9%id 99%cov; S:A0A9Y6DVU6 67.7%id 100%cov)
- 1 candidates: homology QC failed (L:A0A8J0QES8 94.8%id 44%cov; S:A0A8J0QES8 95.4%id 45%cov)

## Y01_YGOB_pillars: candidate units not retained

- 2337 candidates: <4 panel species mapped
- 823 candidates: 
- 10 candidates: no asserted pair
- 3 candidates: mapping 12/18 < 80%
- 2 candidates: mapping 12/16 < 80%
- 2 candidates: mapping 10/13 < 80%
- 2 candidates: mapping 6/8 < 80%
- 2 candidates: mapping 9/12 < 80%
- 2 candidates: mapping 5/8 < 80%
- 1 candidates: mapping 12/19 < 80%
- 1 candidates: mapping 13/19 < 80%
- 1 candidates: mapping 7/10 < 80%
- 1 candidates: mapping 10/16 < 80%
- 1 candidates: mapping 8/12 < 80%
- 1 candidates: mapping 7/11 < 80%
- 1 candidates: mapping 5/9 < 80%
- 1 candidates: mapping 6/9 < 80%
- 1 candidates: mapping 11/19 < 80%
- 1 candidates: mapping 15/19 < 80%
- 1 candidates: mapping 10/18 < 80%
- 1 candidates: mapping 13/17 < 80%
- 1 candidates: mapping 8/15 < 80%
- 1 candidates: excluded: URA1/YKL216W (DHOD) acquired by HGT from bacteria; encoded with its HGT history in families/Y03_URA1_HGT (PI decision)
