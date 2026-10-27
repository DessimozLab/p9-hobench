# Coverage

## Panel x phenomenon (number of catalogue rows tagged; K01 = 1 row, Y01/V07 = 1 row each; tag on any unit counts)

| index | WGD | ALLO | LOSS | TAND/LSE | HGT | DEEP/LECA/LUCA | PROM | FUS | FAST | ANALOG | CONTROL | ANALOG (rows with decoys/distinct_from) | rows |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| E | 0 | 0 | 9 | 6 | 1 | 9 | 4 | 0 | 0 | 0 | 1 | 0 | 9 |
| K | 0 | 0 | 0 | 0 | 0 | 2 | 0 | 0 | 0 | 10 | 2 | 10 | 12 |
| P | 0 | 0 | 2 | 0 | 1 | 3 | 0 | 0 | 0 | 0 | 0 | 0 | 4 |
| V | 4 | 2 | 4 | 3 | 0 | 3 | 0 | 0 | 0 | 0 | 0 | 0 | 5 |
| Y | 1 | 1 | 2 | 0 | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 2 |

Tag vocabulary seen (normalised): ALLO, ANALOG, CONTROL, CONV, DEEP/LECA/LUCA, DOM, HGT, LOSS, NEOF, NESTED, PROM, RETRO, SSD, TAND/LSE, WGD

## Depth coverage (D and H events)

Binning rule (no dated species trees are stored, so bins are an NCBI-rank proxy for age): each event is placed at its youngest-bound level `lvl_lo` (the species-tree node of the LCA mapping/lower interval bound). Level -> NCBI rank via `ncbi_taxid` in panels/<P>/levels.tsv (no-rank/clade nodes inherit the nearest ranked ancestor). **recent** (~<100 Mya): species leaf, or rank species/genus/family (incl. sub-/tribe); **intermediate**: order, class, cohort, subphylum and their super-/infra- ranks; **deep**: phylum, kingdom, domain and above; **LECA/LUCA**: level `LUCA` or `Eukaryota` (panel nodes exactly at LUCA / LECA). Intervals are binned by lower bound, so deep-hi intervals are counted young. Ranks are only a proxy: e.g. Chordata-level (phylum) events are ~550 Mya, family-level may exceed 100 Mya in bacteria.

All D+H events (units counted individually):

| panel | recent | intermediate | deep | LECA/LUCA |
|---|---|---|---|---|
| E | 10 | 0 | 3 | 15 |
| K | 0 | 0 | 0 | 1 |
| P | 0 | 11 | 3 | 0 |
| V | 535 | 49 | 5 | 0 |
| Y | 1591 | 0 | 1 | 0 |
| all | 2136 | 60 | 12 | 16 |

Headline rows only, T1+T2 D events (H events have no tier and are excluded):

| panel | recent | intermediate | deep | LECA/LUCA |
|---|---|---|---|---|
| E | 10 | 0 | 3 | 12 |
| K | 0 | 0 | 0 | 1 |
| P | 0 | 9 | 0 | 0 |
| V | 535 | 29 | 1 | 0 |
| Y | 1591 | 0 | 1 | 0 |
| all | 2136 | 38 | 5 | 13 |

By tier:

| bin | T1 | T2 |
|---|---|---|
| LECA/LUCA | 16 | 0 |
| deep | 12 | 0 |
| intermediate | 47 | 13 |
| recent | 11 | 2125 |

## Comparison with existing benchmarks
Facts for the other benchmarks are taken only from LITERATURE.md and literature/A_benchmarks_metrics.md; "not stated" = not in those files.

| | HOBench (this work) | OrthoBench (Trachana 2011; Emms & Kelly 2020) | SwissTree (Boeckmann 2011) | QfO benchmark service (Altenhoff 2016/2020; Nevers 2022; Altenhoff 2024) |
|---|---|---|---|---|
| Taxonomic scope | 5 panels (E eukaryotes+archaea+bacteria 34 spp; K prokaryote-centred 21; P plants 14; V vertebrates 17; Y yeasts 16) | 12 metazoans; Bilateria only, animal-centred | 19 expert-curated reference gene trees spanning all domains (also described as few "tractable" families) | Reference proteome set of ~78-81 species (48 eukaryotes, 23 bacteria, 7 archaea in 2019-22); test sets vary (VGNC: 8 mammals) |
| Levels scored | Every species-tree node (panel levels.tsv; see depth table) | One (Bilateria); described as a limitation | Pairwise orthology only; no level concept stated | Pairwise; no test scores nesting; OrthoBench group test prototyped, one level |
| Events asserted | Per-family reconciled trees: S, D (exact level or interval), H, with tiers T1/T2/T3 | 70 RefOGs (groups only; no events) | Reconciled/curated reference trees (events implicit in trees; pairwise scored) | None scored: duplication timing, paralogy and nesting are not scored |
| HGT / WGD content | Explicit H nodes (Y03, Y04, K03 etc.) and WGD/ALLO D tags (Y01, V03, V07, ...) | Not stated; challenge tags are duplication/loss, rate, domain architecture, alignment quality | Not stated | YGOB-style WGD "traps" appear only in the separate Salichos 2011 yeast benchmark; no HGT test stated |
| Hierarchical scoring | Per-level HOG P/R/F, B-cubed, three-class relation, stratified by depth (design) | Single-level pair P/R/F per RefOG, (n-1) normalised | Pairwise precision/recall | HOGs flattened to ortholog pairs, nesting discarded |
