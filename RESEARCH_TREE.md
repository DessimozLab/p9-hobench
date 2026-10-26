# RESEARCH TREE — Literature-based benchmark for hierarchical orthology inference

Legend: [ ] planned · [~] in progress · [x] done · [-] abandoned/negative · [!] key finding
Start: see START.txt (2026-10-26T13:17:42Z). Budget 24 h → hard stop 2026-10-27T13:17Z; REPORT.md drafted by ~T+21h.

## Status summary (latest milestone first)

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
- [~] 2. DESIGN.md draft (v0 done) → review rounds (≤3): round 1 running
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
