#!/bin/bash
# stream-count sequences in UniProt RP canonical FASTA (nothing stored). args: upid taxid
u=$1; t=$2
base=https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes
for sr in Eukaryota Bacteria Archaea; do
  url=$base/$sr/$u/${u}_${t}.fasta.gz
  if curl -sfI -m 20 "$url" >/dev/null; then
    n=$(curl -sf -m 600 "$url" | zcat | grep -c '^>'); echo -e "$u\t$t\t$n"; exit 0
  fi
done
echo -e "$u\t$t\tNA"
