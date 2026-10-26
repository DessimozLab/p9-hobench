# E06_myosins curation log
status: abandoned (mapping infeasible under CURATION_PROTOCOL 0.5 coverage rule)

## Read / retrieved
- Sebé-Pedrós et al. 2014 GBE 6:290 (PMID 24443438; journal verified = GBE): DataS1.fasta (725 seqs). Headers are
  Species_index (e.g. Acas_1_10) with NO class labels; classes appear only as clade annotations in the Suppl. Fig. S1 tree
  (PDF vector figure) -> would be T3fig at best.
- Kollmar & Mühlhausen 2017 BMC Evol Biol 17:211 (PMID 28870165): figshare 10.6084/m9.figshare.4565155 (CC BY 4.0)
  downloaded to literature/papers/PMID28870165/raw/figshare_4565155/ (gitignored):
  Myosin_motordomain_alignment.fas md5 2a5ff6c936bf8655bcefd66a83ca7cb8 (7,846 seqs, class-encoded CyMoBase names,
  e.g. HsMyo1A, HsMhc7), Phylogenetic_trees.zip md5 81124aee2761d9ecbc8f61084b15f974.
  sha256: 608858a08b2da1f811dcaf76a7d72cf86597b00a8b58a0f6bca1662db5028d5a / ecfe7224f544739a19d87f69ae89542ebfc69af39eb29c8f8daa4788d16c87b3

## Why abandoned
- Neither source gives accessions for Kollmar names; the only sequences are motor domains (~690 aa) while panel myosins are
  1,000-2,500 aa. Pilot (DIAMOND, same species): HUMAN 38/43, DICDI 13/13, ARATH 17/18 motor domains hit a panel entry at
  >=98% identity over >=90% of the query, but only 0, 1 and 0 reach >=80% target coverage -> fail protocol 0.5
  (>=80% coverage of both sequences). Mapping rate would be <5% -> below the 80% headline criterion.
- CyMoBase (full-length sequences) is interactive only; no bulk export found; not scraped.

## What would make it feasible (for the PI)
- A protocol amendment for domain-defined families: accept >=98% identity over >=90% of a >=300-aa curated domain, RBH,
  >=2-point margin, unique target (same species). Under that rule the pilot maps ~95%, and the Kollmar class names
  (Myo1..Myo79, Mhc = class 2) would give T2 labels; LECA/lineage-specific class dating would need quotes from the
  Kollmar 2017 text (not yet extracted). The same amendment would also recover ~60 motor-domain-only kinesins in E01
  (CHLRE, BATDJ, ENTH1, NAEGR).
- Alternatively request full-length CyMoBase exports from the authors.
