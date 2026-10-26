# Y03_URA1_HGT — curation log (panel Y)
Status: **done** (headline-eligible: 100 % of in-panel source members mapped). Time ~1.5 h.

## Read
- Hall, Brachat & Dietrich 2005 Eukaryot Cell (PMID 15947202; PMC HTML text, not OA subset -> raw/ only; Supplementary
  Table 1 .doc not retrievable). Gives the HGT event (T1), the donor (Lactobacillales), per-species DHOD type statements,
  and GenPept accessions of resequenced DHODs.
- Gojkovic et al. 2004 Mol Genet Genomics (PMID 15014982; closed): abstract only (Europe PMC) -> abstract.txt.
  Full text would add the S. kluyveri family-1A/2 characterisation; not essential (S. kluyveri not in panel).
- Bouwknegt et al. 2021 Fungal Biol Biotechnol (PMID 34656184; OA CC BY): Additional file 1 (331 fungal + 73 bacterial
  Ura9 = Class-II DHOD orthologs, UniProt accessions; HMMER search seeded with L. kluyveri Ura9, ortholog selection with
  proteinortho), Additional file 3 (RAxML-NG tree, 100 bootstraps).
- YGOB v7 Pillars.tab line 4856 (URA1 pillar) and YGOB v1 Synteny_Klactis.tab line 2061.

## Which panel-Y species have which type (sources)
| species | family 1A URA1 (xenolog) | family 2 URA9-type (native) |
|---|---|---|
| YEAST | P28272 YKL216W (YGOB; Hall) | lost (Hall: "endogenous eukaryotic DHOD gene was lost") |
| NAUCC | Q6V3X0 (YGOB; Hall AAQ57200 exact) | none reported (Hall: S. castellii has the type 1a gene) |
| TETPH | G8BZ26 (YGOB) | G8BYI8 (Bouwknegt) |
| VANPO | absent from pillar | A7TRX4 (Bouwknegt) |
| CANGA | absent from pillar | Q6SZS5 (Hall AAR17523 exact; "C. glabrata have the type 2 enzyme") |
| ZYGRC, TORDC | C5DUP9, G8ZXM9 (YGOB) | no source |
| KLULA | Q7Z895 (YGOB; Hall AAQ01773 exact) | present per Hall ("K. lactis ... have both") but no identifier in any source -> not listed |
| LACTH | C5DJY4 (YGOB) | C5DBT0 (Bouwknegt) |
| ASHGO | absent from pillar | Q75CE1 (Bouwknegt; Hall: "family 2 DHOD from A. gossypii") |
| ERECY | absent from pillar | G8JWU1 (Bouwknegt) |
| CANAL | — | Q874I4 (Hall AAP39962 exact; Bouwknegt) |
| YARLI, SCHPO | — | Q6C8K4, P32747 (Bouwknegt) |
| NAUDC, KAZAF | absent from pillar | no source |
Out of panel: K. marxianus AAR17522 (maps to A0ABX6EVC5, KLUMA optional species dropped from panel Y), HENBL (pillar
TBLA0B10140 + Bouwknegt I2H5T2), S. uvarum (Suva_11.4), S. kluyveri, K. waltii, S. paradoxus/mikatae/bayanus.

## Decisions
- Root H (T1), donor OUTSIDE (Lactobacillales not in panel Y) -> URA1 vs URA9-type pairs never in one HOG.
- URA1 subtree: U[ S( S(KLULA,LACTH), S(ZYGRC,TORDC) ), YEAST, NAUCC, TETPH ] — valid on both species_tree.nwk and
  species_tree_alt.nwk. Post-WGD URA1 genes unasserted: YGOB pillar slots are not track-assigned; YGOB v1 scores the
  S. cerevisiae/S. castellii URA1 positions as non-syntenic ('!'); hybrid origin of the WGD lineage.
- URA9 subtree (T3): Bouwknegt tree clades with BS >= 70 restricted to panel members: (Eremothecium pair) BS 94 inside
  (LACTH, Eremothecium, U(VANPO,TETPH)) BS 99 inside (+CANAL) BS 78; CANGA (not in their tree), YARLI, SCHPO at the U root.
- Paper typo: Hall prints "K. lactis (GenPept AAR1773)"; no such accession exists; AAQ01773 (same deposition series,
  NCBI title "dihydroorotate dehydrogenase [Kluyveromyces lactis]") maps exactly to Q7Z895 = the YGOB URA1 gene, so the
  correction is consistent; flagged in mapping_report.tsv.
- Bouwknegt membership used proteinortho (a tool) for ortholog selection; members are only used as family-2 DHOD
  members (corroborated by Hall for CANAL/ASHGO); no orthology is asserted from it beyond T3 clades.

## Problems
- **Gene-ownership overlap with Y01_YGOB_pillars**: unit YGOB_L04856 contains the same 7 URA1 genes (P28272, Q6V3X0,
  G8BZ26, C5DUP9, G8ZXM9, Q7Z895, C5DJY4). validate_panel will report them as members of two families. Recommendation:
  drop unit YGOB_L04856 from Y01 (Y03 is the richer, HGT-aware encoding of the same pillar). Not edited here (not my family).
- Family 2 of KLULA/ZYGRC/TORDC and any DHOD of NAUDC/KAZAF have no source identifiers -> no complete_species.
