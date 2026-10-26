# curation_tools

Toolkit for turning papers into reference gene families: find/download a paper and its supplements,
scan them for sequence identifiers, map identifiers to UniProt, cross-reference paper sequences onto a
benchmark proteome, and fetch the proteomes with a reproducibility manifest.
Run everything from the project root with `export PATH=$PWD/conda_envs/ortho/bin:$PATH`
(extra pip deps installed in that env: `pypdf openpyxl python-docx xlrd`; `pdftotext` is not installed, pypdf is used).

```
python -m curation_tools.<module> --help      # fetch_paper | extract_ids | map_ids | map_by_sequence | fetch_proteome
python -m pytest curation_tools/tests -q              # offline tests (fast)
python -m pytest curation_tools/tests -q --network    # + live Europe PMC / PMC / UniProt / NCBI tests (~1 min)
bash curation_tools/demo.sh                           # end-to-end demo (outputs in curation_tools/demo_output/)
```

## 1. fetch_paper
```
python -m curation_tools.fetch_paper search "orthology gene family duplication" -n 10 --oa --suppl --sort cited
python -m curation_tools.fetch_paper resolve 10.1101/gr.3672305          # DOI | PMID | PMCID -> metadata JSON
python -m curation_tools.fetch_paper fetch 16169922 PMC1197285 [--out literature/papers] [--no-supp]
```
`fetch` writes `literature/papers/PMID<id>/` with `meta.json` (ids, title, authors, year, journal, OA flag, license,
`redistributable`, retrieval date, every file with sha256/bytes/source URL, and an `events` log of each retrieval step incl.
failures, plus `not_retrievable`), `fulltext.txt` (headers `#`, figure captions `[FIGURE id]`, tables as
TSV between `[TABLE id]`...`[/TABLE id]`, supplement captions, references), `tables/*.tsv`, and `raw/` (XML/HTML/PDF, supplements; gitignored).
Exit code 2 = something was not retrievable (also printed as `NOT RETRIEVABLE`).

Cascade (all attempts logged): full text: Europe PMC `fullTextXML` -> PMC OA open-data bucket (S3; the modern
replacement of the OA web service) -> legacy `oa.fcgi` (currently returns 404 from NCBI, kept as fallback) -> NCBI efetch
db=pmc (only if it has a body) -> PMC HTML page -> PDF (pdftotext/pypdf). Supplements: Europe PMC `supplementaryFiles` zip
(slow: 1-3 min) -> files linked on the PMC page.
Papers outside the CC-licensed OA subset (`redistributable=false`) get their derived text under `raw/` (not committed).

Facts found while building this (PMID 16169922, Byrne & Wolfe 2005 YGOB): it is in PMC but *not* in the OA subset
(publisher forbids XML), so only the PMC HTML text is obtainable; its supplementary files (`.xls`, `.pdf`) sit behind an
interactive proof-of-work page on PMC that this tool deliberately does not bypass -> reported as NOT RETRIEVABLE; put them
in `requested/` for manual download.

## 2. extract_ids
```
python -m curation_tools.extract_ids literature/papers/PMID16128622 -o ids.tsv [--unique] [--types uniprot,refseq_protein] [--min-confidence medium]
```
Scans txt/tsv/csv/fasta/json/xml/html, xlsx/xls, docx, doc (crude strings; flagged), pdf, and zip/gz (recursive); skips images.
Columns: `id id_type subtype confidence file location column col_header context` (context = +-1 text lines, or the whole table
row; `location` = line / sheet+row / page+line / docx table+row). Types: refseq_protein/nucleotide, ensembl_protein/transcript/gene
(any species prefix), flybase_gene/protein/transcript/annotation (CG), wormbase_gene/protein/seqname, sgd_systematic/sgd_id,
arabidopsis_agi, uniprot (incl. isoforms, 10-char), uniprot_entry_name, genbank_protein, and `other` (+`subtype`: jgi, phytozome, rice,
maize, pombase, plasmodb, toxodb, tritrypdb, dictybase, vectorbase, fungal, ncbi_gene_loc, gi, ...).
Specific patterns claim text first, so e.g. AT1G01010 is never also reported as UniProt-shaped. `confidence`: `low` for
shapes that collide with ordinary words (bare UniProt/GenBank-like accessions, C. elegans cosmid names) unless the row/column
mentions UniProt/GenBank/accession.

## 3. map_ids
```
python -m curation_tools.map_ids ids.tsv -o mapped.tsv [--taxon 559292]       # ids.tsv from extract_ids (id, id_type[, taxon])
python -m curation_tools.map_ids --id NP_014368.1 --id YAL001C --id Q9Y3G6
```
Batched (1000/job) UniProt ID-mapping jobs with polling and 429/5xx backoff. Per-type cascade of from-databases, e.g.
RefSeq_Protein (as-is, then unversioned), EMBL-GenBank-DDBJ_CDS, Ensembl_Protein/Transcript/Ensembl (+Ensembl_Genomes_*), SGD -> Gene_Name+taxon,
Araport, FlyBase, WormBase, UniProtKB_AC-ID (secondary accessions -> primary, noted), and UniParc fallback for deleted accessions
(`uniparc` column, no acc). `Gene_Name` needs a taxon (`--taxon` or per-row `taxon`; defaults for SGD/AGI/CG/worm names).
Output: `original_id id_type uniprot_acc taxid gene_name reviewed length entry_name mapping_method n_hits note uniparc`.
Unmapped ids are kept (`mapping_method=unmapped`); multiple hits all kept, reviewed first. `<out>.meta.json` records the UniProt release.

## 4. map_by_sequence
```
python -m curation_tools.map_by_sequence --query paper.fasta --target data/proteomes/UP000002311_559292/proteome.fasta -o map.tsv
python -m curation_tools.map_by_sequence --ids ids.tsv --target proteome.fasta --save-query fetched.fasta -o map.tsv [--tool mmseqs]
```
Same-species cross-reference only. Exact match by sha256 of the normalised sequence (`exact`; identical targets ->
`exact_ambiguous`), else best DIAMOND (`--very-sensitive`) or MMseqs2 hit by bit score. Columns: `query_id query_len target_id target_len
mapping_identity identity_over_query query_cov target_cov evalue bitscore method status flags n_exact_targets ambiguous_with`.
`status`: exact | exact_ambiguous | unique | ambiguous (another target within 1 pident point and query_cov >= 0.5) | partial
(`flags`: low_query_cov / low_target_cov, i.e. fragment or domain-level hit) | unmapped | sequence_not_fetched. All queries are kept.
`--ids` fetches sequences: UniProt (UniParc for obsolete), NCBI efetch (RefSeq/GenBank), others via map_ids first.

## 5. fetch_proteome
```
python -m curation_tools.fetch_proteome UP000002311 --out data/proteomes --markdown [--source ftp|rest] [--include-isoforms]
```
Writes `<out>/<UP>_<taxid>/{proteome.fasta, header_map.tsv, manifest.json, raw/}`. Headers become `>ACCESSION`;
`header_map.tsv` keeps the original header, entry/protein/gene name, PE, SV, length, sequence sha256, gene2acc symbol.
`manifest.json`: UniProt release + date, source URLs, sha256 of raw and output files, counts, proteome type, assembly, licence;
`--markdown` prints a DATA.md table row.
**Source choice**: default is the Reference_Proteomes FTP canonical file (`<UP>_<taxid>.fasta.gz`) for reference/representative
proteomes because it is exactly one sequence per gene (UP000002311: 6066 = geneCount) and comes with `gene2acc`; the REST
query `(proteome:UP...)` returns one canonical sequence per *entry* (6067 here: +P9WEJ6) and is used as fallback/for non-reference proteomes.
Release is read from FTP `RELEASE.metalink` (REST header `X-UniProt-Release` for the REST path). Only the "current" release is
downloadable, so the sha256 in the manifest is the reproducibility anchor - keep the manifest.

## Demo results (see `curation_tools/demo_output/`)
- `fetch_paper`: PMID 16169922 -> text via PMC HTML (not OA subset), supplements reported NOT RETRIEVABLE (PMC proof-of-work wall);
  PMC13012417 (hexokinases in Metazoa, OA) -> XML full text + `Table S1.xlsx` inside a nested supplementary zip.
- `extract_ids`: 128 RefSeq proteins etc. from that xlsx with sheet/row/column header; 44 SGD systematic names from the YGOB tables.
- `map_ids`: NP_014368.1 -> P02309; YAL001C-type names via Gene_Name+taxon; ENSMUSP00000000001 -> Q9DC51; Q9Y3G6 (deleted) -> UniParc UPI000006D6F7;
  NP_079406.4 stays unmapped in UniProt's RefSeq cross-references but `map_by_sequence` finds its exact human match (Q2TB90).
- `fetch_proteome`: UP000002311 (6066 proteins) and UP000005640 (20652), UniProt release 2026_03 (manifests in demo_output).
- `map_by_sequence`: 5/5 human HK RefSeq proteins exact; the 2 mouse controls map to human paralogs/orthologs at 85-93 % identity
  (`unique`, not `exact`) - a reminder that this tool is a same-species cross-reference: filter on `status=exact` / identity near 100.

## Known limitations
- PMC supplements of non-OA papers (proof-of-work wall), publisher-only supplements (no PMCID), and scanned PDFs (no OCR) cannot be retrieved/read.
- Supplementary images/figures are skipped; trees in images are not parsed. Legacy `.doc` parsing is a crude string scan.
- ID regexes favour recall; low-confidence rows need curator review. Gene symbols without a systematic ID are not extracted.
- UniProt mappings reflect the current UniProt release; Ensembl/RefSeq ids that were retired may be unmapped (then use map_by_sequence).
