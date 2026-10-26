# RESEARCH TREE — Literature-based benchmark for hierarchical orthology inference

Legend: [ ] planned · [~] in progress · [x] done · [-] abandoned/negative · [!] key finding
Start: see START.txt (2026-10-26T13:17:42Z). Budget 24 h → hard stop 2026-10-27T13:17Z; REPORT.md drafted by ~T+21h.

## Status summary (latest milestone first)

### M4 — References frozen (T+10h35, 23:52 UTC)
Frozen: 26 curated families/collections across 5 panels (headline: Y01 YGOB [1,976 units], Y03 URA1-HGT, K02, K04a–j
decoys, K01 r-proteins [control], E01 kinesins, E02 ARF, E03 Rab, E04 tubulins, E05 MCM [control], E07c, V02 opsins,
V03 Hox, V07 X. laevis L/S [500 units], P02 bHLH, P03 GID1/CXE, P04 MIPs; secondary: E06, E07a/b, P01, V01, V04).
Blind verification: 21/22 pass (V01 fails → secondary). Source audit: 0 violations. Human delivered ~45/50 requested
items (globins, Xenopus, Manning kinome...); used where IDs existed. Round-2 review accepted (testable events, baselines
shown with skill, controls track, pooled confirmatory test).

Vertical slice (Y, K full) scored: positive control perfect; negative controls exposed (and we fixed, logged in
results/ERRATA.md) untestable WGD events and a flat-baseline projection bug. First signal: on YGOB, FastOMA beats
OrthoFinder on merge accuracy (0.94 vs 0.88) and membership F1 (0.92 vs 0.83); testable WGD events are rarely recovered
by either (hidden paralogy is hard). Running: E/P/V/Y/K subset runs, then Y robustness (alt tree, fragmentation,
contamination, taxon removal). Next: score all panels, equivalence test, confirmatory test, difficulty analysis, report.

### M3 — Curation wave 1 + Y/K runs (T+10h10, 23:27 UTC)
Second account-limit outage (~20:20→23:10 UTC). Done before it: hobench v1 (88 tests; interval merge levels, two-sided
events, membership F1, OrthoBench measures, paired/two-stage stats; YGOB positive control perfect in 17 s), panels Y/K/E/V/P
(binary named trees from Shen 2018, Coleman 2021, Burki 2020 etc.; levels.tsv with citations; DATA.md), 272 OMAmer hogmaps,
full runs for Y and K (OF default/dendroblast, FastOMA; Y: OF 24 min, FastOMA 84 min), and families: Y01 YGOB (1,976
pillar units, allopolyploid WGD interval), Y03 URA1 HGT, K01 r-proteins (40 units), K02 trGTPases, K04 5 analogous-enzyme
decoy pairs, E01 kinesins, E02 ARF, E03 Rab, E04 tubulins, E05 MCM, E06 myosins (non-headline, 73% mapped), P01 homeobox
(non-headline 62%), P02 bHLH, P03 GID1/CXE, V02 opsins, V03 Hox, V04 Sox. Abandoned with reasons: K03, Y02 (kept as
cross-check, 141/141 agree with Y01), Y04.

Now: scoring Y/K vertical slice; E/P subsets + runs; V05 teleost 3R + V07 X. laevis collections; V01 globins. Scope cut
(budget): salmonid 4R, chemokines, KRAB-ZNF, MADS, SwissTree, PolB deferred. Freeze target moved to T+13.

### M2 — Design v1 + rate-limit recovery (T+5h00, 18:17 UTC)
Review round 1 (7 independent reviewers + synthesis) found real defects: point-valued reference merge levels mis-score
correct WGD predictions (fixed by interval-valued references), trivial predictors scored well on event recovery and F1
(two-sided recovery, F1 never drops families), β-lactamase classes wrongly called non-homologous (replaced by analogous-
enzyme decoys), FastOMA levels silently renamed. DESIGN.md v1 applies all accepted items; CURATION_PROTOCOL.md makes
provenance mechanical (parsers, verified quotes, blind 30-item re-derivation, freeze by hash).

At T+1h05 13 concurrent agents exhausted the account session limit; all died; work resumed at T+4h55 with ≤5 agents.
Running now: hobench v1 metrics, panels/trees/hogmaps, wrapper + Y/K full runs, curation Y01 (YGOB) and K01/K02.
Remaining curation waves (E, V, P, SwissTree) start as slots free. Freeze target T+12.

### M1 — Literature survey + design v0 (T+0h55)
Four parallel surveys produced LITERATURE.md (+ literature/A,B1,B2,C): no existing benchmark scores HOGs per level or
literature events, so per-level/event-aware metrics fill a real gap. ~85 candidate families screened; ~30 have
machine-readable curated data (kinesins, ARF, Rab, tubulins, MCM, myosins, YGOB pillars, teleost 3R ohnologs, salmonid 4R,
opsins, Sox, plant homeobox, NLR, ribosomal proteins, nitrogenase, β-lactamases...). Several curated DBs are offline
(YGOB, KinBase, HomeoDB) → REQUESTED.md lists 50 items for humans. Panels: Y (16 yeasts), V (14–16 vertebrates),
P (13 plants), E (34 eukaryotes + prokaryote outgroups), K (21 prokaryotes), D (Drosophila, stretch).

hobench core (formats, parsers for OrthoXML/OrthoFinder, merge-level metrics, controls) implemented with 53 passing tests;
curation toolkit in progress. DESIGN.md v0 drafted; review round 1 (7 reviewers) running. Next: revise design, launch
curation agents (blind to tool outputs), download proteomes, build subsets, start tool runs on Y/K.

### M0 — Kickoff (T+0h05)
Environment: 12 cores, 503 GB RAM, no scheduler, no GPU, internet OK (UniProt, NCBI, PyPI, bioconda reachable).
Bioconda offers OrthoFinder 3.1.5 and FastOMA 0.5.0 (both emit hierarchical orthologous groups); env being built in conda_envs/ortho.
Git repo initialised. Next: parallel literature survey (benchmarks/metrics; candidate families; genomes/species panels), then DESIGN.md + review rounds.

## Tree
- [x] 0. Setup: START.txt, git, conda env (orthofinder 3.1.5, fastoma 0.5.0, diamond, mafft, iqtree, fasttree, ete3)
- [x] 1. Literature & data survey → LITERATURE.md, requested/REQUESTED.md
  - [x] 1a. Existing benchmarks & hierarchical-orthology metrics
  - [x] 1b. Candidate gene families with literature-resolved histories (phenomenon-tagged)
  - [x] 1c. Genome/proteome sources, species panels, species trees, licenses
- [~] 2. DESIGN.md v1 after review round 1 (reviews/round1/SYNTHESIS.md); round 2 pending
- [ ] 3. Data: species panels + proteomes (DATA.md)
- [x] 4. Curation (frozen refs-frozen-v1; families/VERIFICATION.md, families/AUDIT.md): literature extraction → reference families (assertions with provenance)
- [~] 5. Metrics library + toy/positive/negative controls (hobench core: 53 tests pass)
- [ ] 6. Tool wrapper (OrthoFinder, FastOMA, baselines) + runs
- [ ] 7. Analysis: difficulty, method agreement, perturbation robustness
- [ ] 8. REPORT.md

## Calibration facts (M0)
- [!] Usable RAM ≈ 50 GB (cgroup), 12 cores. OMAmer LUCA db = 9.9 GB file, 14 GB RSS per omamer process.
- OrthoFinder 3.1.5, 4 yeasts (21k proteins): 2.8 min wall (DIAMOND 31 s, MSA/trees 1.8 min). No N0.tsv in 3.1.5 → root = Orthogroups.tsv.
- FastOMA 0.5.0, same input: 11.5 min wall, dominated by fastoma-infer-subhogs and scratch-FS I/O.
- Extrapolation: full 15–20-species vertebrate/plant panels ⇒ several hours per tool-run → need homology-enriched subsets and/or local-NVMe workdirs.

## Decision log
- D0 (T+0h): Work in parallel tracks with subagents; PI integrates. Opus for design/curation/review; Sonnet for bulk retrieval.
- D1 (T+0h50): Reference = per-family reconciled tree with S/D/U/H nodes + assertion tiers (T1 text, T2 classification, T3 tree-derived); unifying quantity = pair merge level. Why: covers orthology, paralogy, event timing and per-level HOGs with one formalism; unresolved parts are simply not asserted.
- D2 (T+0h50): Full proteomes for Y,K; family-neighbourhood subsets for V,P,E (full runs ≈4–8 h each); validate subset≈full on Y,K. Why: 24 h budget, 12 cores.
- D3 (T+0h50): Curators blind to tool outputs; references frozen (git tag) before scoring. Why: prevent metric/label tuning.
- D4 (T+5h): Accept round-1 synthesis recommendations in full (see DESIGN §12). Why: independent convergence of several reviewers on the same defects; toy tests reproduced them.
- D5 (T+5h): Concurrency cap ≤5 agents; resume interrupted agents rather than relaunching. Why: account session limit killed all 13 agents.
- D6 (T+5h50): Protocol v1.1 domain-only mapping (≥250 aa, ≥90% qcov, ≥98% id, RBH, margin). Why: kinesin/myosin sources give only motor domains; same-species identity over a long domain identifies the gene.
- D7 (T+6h): E06 myosins (73%) and P01 homeobox (62%) stay non-headline: 80% rule pre-registered; not relaxed.
- D8 (T+10h): Cut salmonid 4R/chemokines/KRAB/MADS/SwissTree/PolB; ≤4 concurrent agents, Sonnet for evaluation. Why: second account-limit outage; budget is the bottleneck.
- D9 (T+10h30): Freeze now (earlier than T+12) because all launched curation finished; V01 FAIL → secondary (no re-curation; budget).
- D10 (T+10h30): Keep compressed tool outputs in git (≤20 MB/run) despite audit suggestion — brief requires results in git; metrics reproducible without re-running tools.
