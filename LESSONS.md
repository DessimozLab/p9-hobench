# LESSONS — process rules for running this project efficiently
(one line per rule, with evidence; pruned at each milestone)

- Activate tools by prepending conda_envs/ortho/bin to PATH (orthofinder/nextflow fail without it: java/module not found). Evidence: M0 env check.
- Usable RAM is ~50 GB (slurm cgroup), NOT the 503 GB shown by free; omamer+LUCA needs 14 GB/process → cap parallel omamer at 2. Evidence: 4 parallel omamer killed; numpy probe died after 48 GB (M0).
- Put tool work dirs (nextflow -work-dir, OrthoFinder temp) on local NVMe /tmp (679 GB free), not on /scratch: FastOMA on 4 yeasts spent 31 min sys vs 22 min user time (11.5 min wall) from small-file I/O on scratch. Evidence: runs/timing (M0).
- Nextflow/FastOMA needs absolute paths for --input/--output_folder ("Not a valid path value"). Evidence: first FastOMA attempt failed in 1 s.
