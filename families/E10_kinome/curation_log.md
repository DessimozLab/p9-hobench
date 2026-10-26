# E10_kinome curation log
status: abandoned (protocol: <4 panel species with member identifiers)

## Read
- Manning et al. 2002 Science 298:1912 (PMID 12471243): requested/06_* copied to literature/papers/PMID12471243/raw/
  (sha256 as in requested/RECEIVED.tsv: tables_excel e6522eaa…, tables_ascii e8676a83…, text bd87e2b5…, SOM 098af09f…;
  text extracts raw/fulltext.txt, raw/som.txt, raw/tables_ascii/).
- Manning, Plowman, Hunter, Sudarsanam 2002 TIBS 27:514 (PMID 12368087; doi 10.1016/S0968-0004(02)02179-5):
  requested/07_text.pdf -> literature/papers/PMID12368087/raw/ (sha256 1eeb70e0…).

## Why abandoned
- Table S1 (sequences + group/family/subfamily/pseudogene) and S4 (pseudogenes) list HUMAN kinases only.
- Table S6 gives per-subfamily domain COUNTS for yeast, worm, fly, human - no yeast/fly gene identifiers or sequences.
- TIBS paper has no per-gene lists (it points to www.kinase.com, offline; the Wayback KinBase files in
  literature/raw/B1/kinbase/ are again human-only + the same count table).
- => members mappable in panel E: HUMAN only (1 species). Cross-species orthology cannot be asserted; protocol requires
  >=4 panel species.

## What would make it feasible
- KinBase per-species kinome FASTA/classification for S. cerevisiae and D. melanogaster (kinase.com/web/current/kinbase,
  e.g. a Wayback snapshot of the 'Sc'/'Dm' downloads), plus Manning et al. 2008 PNAS Monosiga kinome (PMC, supplement)
  and Goldberg 2006 Dictyostelium kinome. With yeast+fly+MONBE+human the TIBS dating rule (subfamilies 'present in all
  four genomes' -> splits at/above Opisthokonta) could be applied.
