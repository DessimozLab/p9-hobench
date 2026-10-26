# RESEARCH TREE — Literature-based benchmark for hierarchical orthology inference

Legend: [ ] planned · [~] in progress · [x] done · [-] abandoned/negative · [!] key finding
Start: see START.txt (2026-10-26T13:17:42Z). Budget 24 h → hard stop 2026-10-27T13:17Z; REPORT.md drafted by ~T+21h.

## Status summary (latest milestone first)

### M0 — Kickoff (T+0h05)
Environment: 12 cores, 503 GB RAM, no scheduler, no GPU, internet OK (UniProt, NCBI, PyPI, bioconda reachable).
Bioconda offers OrthoFinder 3.1.5 and FastOMA 0.5.0 (both emit hierarchical orthologous groups); env being built in conda_envs/ortho.
Git repo initialised. Next: parallel literature survey (benchmarks/metrics; candidate families; genomes/species panels), then DESIGN.md + review rounds.

## Tree
- [~] 0. Setup: START.txt, git, conda env (orthofinder 3.1.5, fastoma 0.5.0, diamond, mafft, iqtree, fasttree, ete3)
- [~] 1. Literature & data survey → LITERATURE.md
  - [~] 1a. Existing benchmarks & hierarchical-orthology metrics
  - [~] 1b. Candidate gene families with literature-resolved histories (phenomenon-tagged)
  - [~] 1c. Genome/proteome sources, species panels, species trees, licenses
- [ ] 2. DESIGN.md draft → review rounds (≤3)
- [ ] 3. Data: species panels + proteomes (DATA.md)
- [ ] 4. Curation: literature extraction → reference families (assertions with provenance)
- [ ] 5. Metrics library + toy/positive/negative controls
- [ ] 6. Tool wrapper (OrthoFinder, FastOMA, baselines) + runs
- [ ] 7. Analysis: difficulty, method agreement, perturbation robustness
- [ ] 8. REPORT.md

## Decision log
- D0 (T+0h): Work in parallel tracks with subagents; PI integrates. Opus for design/curation/review; Sonnet for bulk retrieval.
