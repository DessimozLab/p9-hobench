# Round 2 — PI synthesis and decisions (T+11h30)
Reviewers: 01 PI+statistics (Opus), 02 data/agent audit (Sonnet). Both saw DESIGN v1 and the Y/K vertical-slice results.

| id | item | decision |
|---|---|---|
| R2-1 (01 M1) | untestable events inflate trivial baselines | ACCEPT (modified): event testable iff cross-copy interval bounded above within family range OR ≥1 copy has a reference-together pair; this provably gives 0 to all-singletons and one-big-HOG; report n testable. Stricter "≥2 genes from ≥2 species per copy" kept as sensitivity. |
| R2-2 (01 M2) | show trivial baselines next to every endpoint; skill (s−b)/(1−b) | ACCEPT; S/D-balanced merge accuracy added as post-results secondary (logged in results/ERRATA.md). |
| R2-3 (01 M3) | membership F1 NA without complete cells/decoys; rename one_hog → "one big HOG" | ACCEPT |
| R2-4 (01 M4) | sanity gate: method identical to all-singletons on >50% units | ACCEPT (of_flat projection bug being fixed) |
| R2-5 (01 M5) | controls track (K01 RPROT, K02, E05 MCM) separate from headline means | ACCEPT |
| R2-6 (01 M6, 02 #3) | freeze + ERRATA stating Y/K results were seen; finish blind verification; subset-vs-full equivalence test | ACCEPT; Y/K subset runs queued |
| R2-7 (01 M7) | per-panel descriptive only; confirmatory = one paired permutation test per primary endpoint over pooled items, Holm over 4; relabel min–max columns | ACCEPT |
| R2-8 (02 #1) | write scripts/audit_sources.py | DONE (0 violations) |
| R2-9 (02 #2) | runs/ outputs in git (21 MB) | REJECT (brief: "commit results to git"); keep compressed HOG outputs ≤20 MB/run so metrics are reproducible without re-running tools; raw work dirs never committed. |
| R2-10 (02 #4) | DATA.md row for rice MSU file | DONE (source URL verified by re-download sha256) |
| R2-11 (02 #5) | scoring must reject manifests not `done` | ACCEPT |
| R2-12 (02 #6) | Y03 Gojković quote not T1 | DONE (demoted to supporting T3) |
| R2-13 (02) | COG-based membership (K01, K04) | ACCEPT with disclosure in REPORT limitations |
