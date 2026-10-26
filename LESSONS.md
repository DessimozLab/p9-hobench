# LESSONS — process rules for running this project efficiently
(one line per rule, with evidence; pruned at each milestone)

- Activate tools by prepending conda_envs/ortho/bin to PATH (orthofinder/nextflow fail without it: java/module not found). Evidence: M0 env check.
- Usable RAM is ~50 GB (slurm cgroup), NOT the 503 GB shown by free; omamer+LUCA needs 14 GB/process → cap parallel omamer at 2. Evidence: 4 parallel omamer killed; numpy probe died after 48 GB (M0).
- Put tool work dirs (nextflow -work-dir, OrthoFinder temp) on local NVMe /tmp (679 GB free), not on /scratch: FastOMA on 4 yeasts spent 31 min sys vs 22 min user time (11.5 min wall) from small-file I/O on scratch. Evidence: runs/timing (M0).
- Nextflow/FastOMA needs absolute paths for --input/--output_folder ("Not a valid path value"). Evidence: first FastOMA attempt failed in 1 s.
- Never let demo/test scripts write into shared data dirs (data/proteomes): toolkit demo wiped 109 freshly downloaded proteomes. Evidence: M1, re-download needed.
- Never run >5 subagents concurrently: 13 parallel Opus agents exhausted the account session limit within ~1 h and ALL died mid-task (lost T+1h05→T+4h54). Prefer resuming agents (SendMessage) over relaunching; prefer Sonnet for table parsing/engineering; tell agents to grep/head instead of reading big files. Evidence: rate-limit 429 at ~14:20 UTC, reset 18:10 UTC.
