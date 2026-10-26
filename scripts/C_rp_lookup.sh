#!/bin/bash
# usage: C_rp_lookup.sh "regex" ...  : grep species name in UniProt RP 2026_03 table
F=/scratch/nromashc/agentic-ai/workspaces/mywork/literature/raw/C/UniProt_RP_2026_03_species.tsv
for p in "$@"; do echo "## $p"; awk -F'\t' -v p="$p" '$8 ~ p {print $1"\t"$2"\t"$4"\t"$5"\t"$6"\t"$8}' $F | cut -c1-150; done
