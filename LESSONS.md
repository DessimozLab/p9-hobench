# LESSONS — process rules for running this project efficiently
(one line per rule, with evidence; pruned at each milestone)

- Activate tools by prepending conda_envs/ortho/bin to PATH (orthofinder/nextflow fail without it: java/module not found). Evidence: M0 env check.
- Usable RAM is ~50 GB (slurm cgroup), NOT the 503 GB shown by free; omamer+LUCA needs 14 GB/process → cap parallel omamer at 2. Evidence: 4 parallel omamer killed; numpy probe died after 48 GB (M0).
