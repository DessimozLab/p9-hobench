# V02_opsins curation log (protocol v1.1)
Status: **validated, headline-eligible** (89% of in-panel S1 rows uniquely mapped).

Read: Lagman et al. 2013 BMC Evol Biol (PMID 24180662) full text + Additional file 1 (S1.xlsx, sheet OPN);
Ramirez et al. 2016 GBE Supplemental Table T1 (inspected, not used: names are database protein names, no per-sequence
classification -> no T2 labels; the classification lives only in tree figures).

Decisions
- Scope = vertebrate visual opsins (LWS, SWS1, SWS2, RH2, RH1). OPN3 (outgroup), pinopsin, V/A opsin and Ciona opsins
  (S1 'additional' section) excluded: Lagman's synteny does not place them in the visual-opsin duplications.
- Gar loci mapped by Ensembl locus overlap at the S1 gene start (LepOcu1); the IDs printed for gar in S1 are the
  heterologous query proteins. Retired zebrafish/frog IDs mapped through S1 transcript IDs / the NCBI sequence named in
  the S1 comment, then same-species sequence mapping at 99% (TAND-tagged family).
- Human OPN1MW/OPN1MW2 excluded as identical copies (3 identical accessions in the proteome), so the primate LWS/MWS
  duplication is not encoded.
- Tree: root SWS|LWS local duplication pre-2R (lvl_lo=Vertebrata, lvl_hi=ABOVE_RANGE); SWS -> SWS1/SWS2/RH1/RH2 by 2R
  (interval Gnathostomata..Vertebrata; Lagman places it before the lamprey split); RH1: retrotransposition at
  Neopterygii separates exorh (intron-containing; gar RH1-1, teleost RH1-1) from rho (gar RH1-2, teleost RH1-2);
  zebrafish rho/rhol + medaka RH1-2 = U (paper: 3R origin not certain); teleost local LWS/RH2 duplicates = U (sharing
  across teleosts not stated); gar SWS1-1/-2 = local duplicate on the LEPOC branch.
- complete_species: none (2013 assemblies; chicken genes missing from WASHUC2, frog RH2 possibly an assembly gap).
Problems: medaka S1 gives one Ensembl gene ID for two opsins twice (LWS-2/SWS2-1, RH2-2/RH2-3) -> ambiguous.
Non-visual opsins (Ramirez) and lamprey/shark members not covered.
Time: ~1.5 h.
