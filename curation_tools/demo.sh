#!/usr/bin/env bash
# End-to-end demonstration of the curation toolkit on real examples. Run from the project root.
# Writes small outputs to curation_tools/demo_output/; downloads go to literature/papers/ and data/proteomes/ (git-ignored parts).
set -euo pipefail
cd "$(dirname "$0")/.."
export PATH="$PWD/conda_envs/ortho/bin:$PATH"
OUT=curation_tools/demo_output
mkdir -p "$OUT"

echo "== 1a. search"
python -m curation_tools.fetch_paper search "whole genome duplication ohnologs gene family retention" -n 8 --oa --suppl --sort cited > $OUT/search.tsv

echo "== 1b. fetch Byrne & Wolfe 2005 (PMID 16169922; in PMC but not OA-subset -> HTML text, supplements blocked)"
python -m curation_tools.fetch_paper fetch 16169922 || echo "(exit $? = some parts not retrievable, see meta.json)"

echo "== 1c. fetch an OA paper with supplementary tables (hexokinase family in Metazoa, GBE 2026, PMC13012417; ~2 min)"
python -m curation_tools.fetch_paper fetch PMC13012417
HK=literature/papers/PMID41793121

echo "== 2. extract identifiers from the YGOB paper (tables + text)"
python -m curation_tools.extract_ids literature/papers/PMID16169922 -o $OUT/ygob_ids.tsv
python -m curation_tools.extract_ids literature/papers/PMID16169922 --unique --min-confidence medium -o $OUT/ygob_ids_unique.tsv
echo "== 2b. identifiers in the hexokinase paper's supplement (xlsx inside a nested zip)"
python -m curation_tools.extract_ids $HK --unique --min-confidence medium -o $OUT/hk_ids_unique.tsv

echo "== 3. map ids to UniProt (a few RefSeq / Ensembl / SGD / AGI / obsolete accession; plus the first ids from YGOB)"
printf 'id\tid_type\ttaxon\nNP_014368.1\trefseq_protein\t\nNP_001234.1\trefseq_protein\t\nXP_001234567.1\trefseq_protein\t\nENSMUSP00000000001\tensembl_protein\t\nENSDARP00000052337\tensembl_protein\t\nFBgn0000008\tflybase_gene\t\nAT1G01010.1\tarabidopsis_agi\t\nAAA34695.1\tgenbank_protein\t\nQ9Y3G6\tuniprot\t\nP12345\tuniprot\t\n' > $OUT/map_ids_input.tsv
cut -f1,2 $OUT/ygob_ids_unique.tsv | grep sgd_systematic | head -8 | awk -F'\t' '{print $1"\t"$2"\t559292"}' >> $OUT/map_ids_input.tsv
cut -f1,2 $OUT/hk_ids_unique.tsv | grep refseq_protein | head -8 | awk -F'\t' '{print $1"\t"$2"\t"}' >> $OUT/map_ids_input.tsv
python -m curation_tools.map_ids $OUT/map_ids_input.tsv -o $OUT/map_ids_output.tsv

echo "== 5. fetch S. cerevisiae reference proteome"
python -m curation_tools.fetch_proteome UP000002311 --out data/proteomes --markdown | tee $OUT/fetch_proteome.stdout
cp data/proteomes/UP000002311_559292/manifest.json $OUT/UP000002311_manifest.json

echo "== 4. map paper sequences onto the proteome (RefSeq/GenBank/SGD/obsolete/human control)"
printf 'id\nNP_014368.1\nNP_009362.1\nYAL001C\nP34111\nQ9Y3G6\nNP_001234.1\nAAA34695.1\n' > $OUT/map_by_seq_ids.tsv
python -m curation_tools.map_by_sequence --ids $OUT/map_by_seq_ids.tsv --target data/proteomes/UP000002311_559292/proteome.fasta \
    --save-query $OUT/map_by_seq_query.fasta -o $OUT/map_by_seq_output.tsv
echo "== 4b. hexokinase paper: first 7 RefSeq proteins of Table S1 (5 human + 2 mouse as cross-species negative control) -> human proteome"
python -m curation_tools.fetch_proteome UP000005640 --out data/proteomes --markdown | tee -a $OUT/fetch_proteome.stdout
cp data/proteomes/UP000005640_9606/manifest.json $OUT/UP000005640_manifest.json
(echo id; cut -f1,2 $OUT/hk_ids_unique.tsv | grep refseq_protein | head -7 | cut -f1) > $OUT/hk_human_ids.tsv
python -m curation_tools.map_by_sequence --ids $OUT/hk_human_ids.tsv --target data/proteomes/UP000005640_9606/proteome.fasta \
    --save-query $OUT/hk_human_query.fasta -o $OUT/hk_human_map_output.tsv
echo "demo finished; see $OUT"
