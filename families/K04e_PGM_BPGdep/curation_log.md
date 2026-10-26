# K04e_PGM_BPGdep — curation log (K04 isofunctional decoys, panel K)
Status: **done** (decoy/merge-test family; no orthology asserted inside). Time: ~1.5 h for all ten K04 families.

## Read
- Omelchenko et al. 2010 Biol Direct 5:31 (PMID 20433725, OA CC BY): full text (Europe PMC XML) + Additional files 1-3
  (PDF; Table S1 = 1998 list update, Table S2 = final list of confirmed NISE with UniProt example, PDB, SCOP fold, Pfam).
- NCBI COG 2020 (cog-20.def.tab, cog-20.cog.csv restricted to panel K = derived/cog20_panelK.csv).

## Decisions
- Pairs chosen = EC nodes listed in Table S2 as "Two distinct structural folds" whose two forms are each one COG with
  members in >=5 panel-K species: LysRS class I (COG1384) / class II (COG1190); SOD Fe,Mn (COG0605) / Cu,Zn (COG2032);
  PGM BPG-dependent (COG0588) / BPG-independent (COG0696); shikimate kinase P-loop (COG0703) / archaeal GHMP (COG1685);
  3-dehydroquinate dehydratase type I AroD (COG0710) / type II AroQ (COG0757).
- Form <-> COG link: COG definition name + gene symbol (cog-20.def.tab), and the PDB structure cited by both sources
  where present (1BBU LysRS-II, 1IRX LysRS-I, 1ESO Cu/Zn SOD, 1EJJ iPGM). The Table S2 example proteins of E. coli fall
  in the expected COG (cog20_panelK.csv: b2890 lysS -> COG1190, b1656 sodB -> COG0605, b1646 sodC -> COG2032,
  b0755 gpmA -> COG0588, b0388 aroL -> COG0703, b1693 aroD -> COG0710); AROK_METJA = MJ1440 (Table S2) is the METJA
  COG1685 member.
- **Not used:** fructose-bisphosphate aldolase class I vs class II — Table S1 (1998 list update) marks EC 4.1.2.13 as
  "Same fold, same superfamily" (both TIM barrels) and it is absent from the final NISE list (Table S2) -> not a
  non-homologous pair, cannot be a decoy. Carbonic anhydrase alpha (COG3338) has a single panel-K member (HELPY);
  gamma-CA COG0663 mixes acetyltransferases -> skipped. Ni-SOD has no panel-K COG member. Archaeal 2,3-BPG-independent
  PGM (COG3635, AlkP superfamily like iPGM) not placed in either PGM family (its relation to COG0696 is not stated).
- Tree: single U node (no source resolves orthology vs xenology/paralogy; >1 copy in several species, e.g. E. coli
  LysS/LysU, SodA/SodB, AroK/AroL). Decoys = all members of the partner family; distinct_from partner at_level ALL.
- Mapping: COG proteins -> pinned UniProt 2026_03 proteome of the same species (exact sha256 or unique DIAMOND
  >=98 % id / >=80 % cov both / 2-point margin) via literature/papers/COG2020/seqmap_K04.py (output derived/seqmap_K04/,
  separate from the K01/K02 seqmap files). 100 % of non-excluded COG members mapped. COG membership class 3 excluded.
- BRADU, GEOSL: not in COG 2020 -> no members (reported as not_in_panel rows in mapping_report.tsv); complete_species
  empty (COG membership is not a declared genome-wide exhaustive search for each form).

## Problems / caveats
- panels/K not yet available: validated against literature/raw/C/trees/panel_K_ncbi_topology.nwk (internal nodes
  auto-named); no levels used by these families, so nothing depends on level names. Re-run validation once panels/K exists.
- The same five decoy pairs can be re-declared in panel E (eukaryotic members) later.
