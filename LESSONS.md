# LESSONS — process rules for running this project efficiently
(one line per rule, with evidence; pruned at each milestone)

- Prepend conda_envs/ortho/bin to PATH and use absolute paths for nextflow/FastOMA --input/--output_folder. Evidence: M0 env check (java/module not found); first FastOMA attempt failed in 1 s.
- Plan for ~50 GB usable RAM (cgroup), not the 503 GB `free` shows; cap parallel omamer at 2 (14 GB each). Evidence: 4 parallel omamer killed; numpy probe died at 48 GB (M0).
- Put tool work dirs (nextflow -work-dir, OrthoFinder temp) on local NVMe /tmp, not /scratch. Evidence: FastOMA on 4 yeasts, 31 min sys vs 22 min user, 11.5 min wall (M0).
- Never let demo/test scripts write into shared data dirs. Evidence: toolkit demo wiped 109 downloaded proteomes (M1).
- Cap concurrent agents at <=4-5; resume (SendMessage) rather than relaunch; use Sonnet for parsing/engineering; tell agents to grep/head big files. Evidence: 13 agents hit 429 at T+1h05, all died (lost ~3.8 h); second outage ~20:20-23:10 UTC (D5, D8).
- Run negative controls (all-singletons, flat baseline, trivial predictors) on real data before trusting any metric. Evidence: round-1 reviewers + toys found point-level refs and trivial-predictor wins; singletons scored 0.86 YGOB event recovery, exposing untestable events (ERRATA E1, E2).
- Use independent multi-reviewer rounds on design before building at scale. Evidence: round 1 (7 reviewers) found real defects (interval refs, beta-lactamase error, silent FastOMA renaming); fixed in DESIGN v1 (M2).
- Ask humans for papers/data early and check requested/ every milestone (scripts/check_requested.sh). Evidence: REQUESTED.md ~90% fulfilled in 1.5 h but not checked until T+10h; ~45/50 delivered (M1, M4).
- Check "literature" tables for hidden tool provenance (denylist audit early, not at curation end). Evidence: Braasch 2016 ohnologs derive from Ensembl Compara, V05 cancelled (T+11h); audit_sources.py found 0 violations after.
- Make provenance mechanical: parsers, verified quotes, blind re-derivation, freeze by hash/tag; pre-register thresholds and do not relax. Evidence: blind verification 21/22 pass; 80% mapping rule kept (D7, D3, M4).
- Abandon fast with a logged reason when a family lacks panel-species IDs (check >=4 panel species up front). Evidence: K03, E09, E10, P06, Y04 abandoned; E06 and P01 stayed non-headline (M3).
- Cut scope early when budget is the bottleneck, and freeze as soon as launched work finishes. Evidence: D8 cut six items; D9 froze at T+10h30 vs T+12 target (M4).
