# Curation log: K02_trGTPases (panel K)

Status: **done** (T1 root event + T2 within-group orthology). Curator: K-panel curation agent, 2026-10-26, ~45 min.

## Read
- Baldauf, Palmer & Doolittle 1996 PNAS (PMID 8755547, PMC38819): not in the PMC OA subset; only the PMC HTML page
  (abstract + scanned page images) is retrievable; PDF behind PMC proof-of-work / Europe PMC render failed. Abstract used.
- Iwabe et al. 1989 PNAS (PMID 2531898, PMC298494): same situation; abstract used.
- Creti et al. 1994 PNAS (PMID 8159735): PMC page returned a reCAPTCHA (deleted); abstract stored from the Europe PMC REST
  API (literature/papers/PMID8159735/abstract.txt + epmc_core.json, sha256 in MANIFEST.sha256).
- Gogarten 1989 (PMID 2528146) fetched but concerns V-ATPase subunits, not used.
- NCBI COG 2020 (see K01 logs): COG0050 (EF-Tu), COG5256 (EF-1alpha), COG0480 (EF-G/EF-2).

## Decisions
- Event: Baldauf 1996 abstract "EF-Tu and EF-G are the products of an ancient gene duplication, which appears to predate
  the divergence of all extant organismal lineages" (+ Iwabe 1989) -> root D, T1, no explicit level (defaults to
  [panel root, ABOVE_RANGE]; root_single_copy false) = the two groups never share a HOG within the panel.
- Group membership: COG 2020 membership in the 19 panel genomes covered by COG (BRADU, GEOSL absent). Archaeal
  EF-1alpha (COG5256) is placed in the EF-Tu group and archaeal EF-2 (in COG0480) in the EF-G group per Creti 1994
  ("EF-2(G)", "EF-1 alpha(Tu)" global phylogenies across Archaea/Bacteria/Eukarya). Labels: EF-Tu / EF-G.
- Within each group: single-copy species -> one S polytomy (T2 rule b). Multi-copy species (T2 rule c, unresolved):
  EF-Tu ECOLI tufA/tufB (P0CE47/P0CE48: distinct UniProt entries, near-identical; no source dates the duplication or
  documents conversion -> not tagged CONV); EF-G: PSEAE, MYCTU, SYNY3 (3), DEIRA, THEMA, BACTN, LEPIN (EF-G paralogs
  such as EF-G2/spd; not resolved by any source in hand) -> these genes hang from a group-level U node, so only their
  cross-group pairs (EF-Tu vs EF-G) are scored.
- Identical tuf copies (PSEAE tufA/tufB, DEIRA, AQUAE) are already collapsed to one UniProt entry; the extra COG protein
  is reported 'excluded (identical/collapsed copy)' (protocol rule 6).
- Not included: IF2/eIF5B, LepA/EF4, SelB, RF3, BipA/TypA, aIF2gamma (no stored source asserts their placement); they
  are homologs of the members and may legitimately appear in predicted HOGs (not scored).
- complete_species: none (no source states a genome-wide search for these genes in the panel strains).

## Problems / follow-up
- Full texts of Baldauf 1996 and Iwabe 1989 would allow T3 within-group structure; Atkinson & Baldauf 2011 (MBE, EF-G
  paralog classes) would resolve the EF-G paralogs. Not essential for the headline (root D) assertion.
