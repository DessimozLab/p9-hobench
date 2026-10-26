# RESEARCH TREE — Literature-based benchmark for hierarchical orthology inference

Legend: [ ] planned · [~] in progress · [x] done · [-] abandoned/negative · [!] key finding
Start: see START.txt (2026-10-26T13:17:42Z). Budget 24 h → hard stop 2026-10-27T13:17Z; REPORT.md drafted by ~T+21h.

## Status summary (latest milestone first)

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
- [ ] 4. Curation: literature extraction → reference families (assertions with provenance)
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
