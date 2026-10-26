# REQUESTED — full texts / data files we could not retrieve automatically

**How to help:** drop any file into this folder (`requested/`), named with the rank number as prefix,
e.g. `01_Pillars.tab`, `02_Scannell2006.pdf`, `05_Niimura2014_SuppTables.zip`. The project checks this folder at
every milestone (`scripts/check_requested.sh`) and integrates new files; work continues without them.
PDFs and publisher files are git-ignored (not redistributed); only extracted, citable facts (IDs, labels, events, with
page/table provenance) enter the benchmark.

Ranking criterion: expected gain for the benchmark = (family is on the shortlist) × (data cannot be obtained otherwise)
× (phenomenon/taxon coverage it unlocks). Sources: literature/B1_*.md §5, literature/B2_*.md §5.
Status legend: `open` · `received` · `integrated`.

| Rank | Item | Why it matters (family / phenomenon) | Status |
|---|---|---|---|
| 01 | **YGOB `Pillars.tab` (full file, ygob.ucd.ie/browser/ — server refuses our cluster; Wayback copy truncated at 1 MB)** | Genome-wide synteny-curated yeast WGD ortholog/ohnolog truth set; hidden paralogy/pseudo-orthology (WGD, LOSS) | open |
| 02 | Scannell DR et al. 2006 Nature 440:341, doi:10.1038/nature04562 (+SI) | Reciprocal ohnolog loss after yeast WGD → pseudo-orthologs (LOSS/ILS-like) | received |
| 03 | Li F-W et al. 2015 PNAS phytochrome Dryad package doi:10.5061/dryad.5rs50 (Dryad blocks scripts) | Phytochromes: plant WGDs + neochrome HGT hornwort→fern | received |
| 04 | Hoffmann FG et al. 2012 MBE 29:303 `supp_msr207_Supplementary_Table1.pdf` (PMC3245541) and Hoffmann 2010 PNAS `1006756107_sd01.txt` (PMC2922537) | Vertebrate globins: 2R, tandem arrays, differential loss | received |
| 05 | Niimura Y et al. 2014 Genome Res 24:1485, doi:10.1101/gr.169532.113 — supplementary OGG membership tables | Mammalian olfactory receptors: extreme lineage-specific expansion | received |
| 06 | Manning G et al. 2002 Science 298:1912, doi:10.1126/science.1075762 (+ supplementary tables) | Human kinome classification (group/family/subfamily); promiscuous kinase domain | received |
| 07 | Manning G et al. 2002 TIBS 27:514, doi:10.1016/S0968-0004(02)02179-5 | Kinomes yeast→fly→worm→human; subfamily origins | received |
| 08 | Session AM et al. 2016 Nature 538:336, doi:10.1038/nature19840 — Supplementary Tables (L/S homeolog lists) | Xenopus laevis allotetraploidy (homeologs) | received |
| 09 | Shao Z-Q et al. 2016 Plant Physiol 170:2095, doi:10.1104/pp.15.01487 (+SI) | Angiosperm NLR ancestral lineages (promiscuous domains, tandem expansion) | received |
| 10 | Parenicová L et al. 2003 Plant Cell 15:1538, doi:10.1105/tpc.011544 (+SI) | Arabidopsis MADS-box complete list with phylogenetic groups (plant WGDs) | received |
| 11 | Amores A et al. 1998 Science 282:1711, doi:10.1126/science.282.5394.1711 | Zebrafish Hox clusters (teleost 3R) | received |
| 12 | Kuraku S, Meyer A 2009 Int J Dev Biol 53:765, doi:10.1387/ijdb.072533km | Teleost Hox cluster evolution | received |
| 13 | Lin Z et al. 2006 PNAS — SI (RecA/RadA/Rad51/Dmc1 sequences and trees) | Three-domain family; LECA duplication; organellar transfer | received |
| 14 | Elias M et al. 2012 J Cell Sci — supplementary data (Rab GTPases across eukaryotes) | LECA Rab paralogs; small-GTPase domain merging | received |
| 15 | Imbeault M et al. 2017 Nature 543:550, doi:10.1038/nature21683 (+SI) | KRAB-ZNF phylogeny and ages (worst-case merging) | received |
| 16 | Kominek J et al. 2019 Cell — Supplemental Tables 1–4 (PMC7295392) | Fungal HGT (enterobactin operon) | received |
| 17 | Shen X-X et al. 2018 Cell — supplementary xlsx 8–13 (PMC6291210) | Budding-yeast species tree and HGT/gene content | received |
| 18 | Cobbe N, Heck MM 2004 MBE 21:332, doi:10.1093/molbev/msh023 | SMC protein phylogeny (LECA paralogs) | received |
| 19 | Muller J et al. 2005 MBoC 16:5736, doi:10.1091/mbc.e05-06-0508 | Actin-related proteins classification | received |
| 20 | Leung KF, Dacks JB, Field MC 2008 Traffic, doi:10.1111/j.1600-0854.2008.00797.x | ESCRT across eukaryotes | received |
| 21 | Bertrand S et al. 2004 MBE 21:1923, doi:10.1093/molbev/msh200 | Nuclear receptor evolution (2R, nematode expansion) | received |
| 22 | Brindefalk B et al. 2007 MBE, doi:10.1093/molbev/msl202 | Mitochondrial aaRS origins (EGT/HGT) | received |
| 23 | Duchêne A-M et al. 2005 PNAS, doi:10.1073/pnas.0504682102 | Arabidopsis organellar aaRS dual targeting | received |
| 24 | Shiu S-H et al. 2004 Plant Cell 16:1220, doi:10.1105/tpc.020834 (+SI) | RLK subfamily membership (Arabidopsis + rice) | received |
| 25 | Lehti-Shiu MD et al. 2009 Plant Physiol 150:12, doi:10.1104/pp.108.134353 | RLK/Pelle across land plants | received |
| 26 | Becker A, Theissen G 2003 Mol Phylogenet Evol, doi:10.1016/s1055-7903(03)00207-0 | MADS clade definitions | received |
| 27 | Nam J et al. 2003 MBE, doi:10.1093/molbev/msg152 | MADS duplication timing | received |
| 28 | Karlgren A et al. 2011 Plant Physiol, doi:10.1104/pp.111.176206 | PEBP (FT/TFL1/MFT) across land plants | received |
| 29 | Hittinger CT, Carroll SB 2007 Nature 449:677, doi:10.1038/nature06151 | GAL1/GAL3 ohnolog subfunctionalisation | received |
| 30 | Slot JC, Rokas A 2011 Curr Biol, doi:10.1016/j.cub.2010.12.020 | Fungal gene-cluster HGT | received |
| 31 | Marcet-Houben M, Gabaldón T 2010 TIG, doi:10.1016/j.tig.2009.11.007 | Bacteria→fungi HGT catalogue | received |
| 32 | Fitzpatrick DA 2012 FEMS Microbiol Lett, doi:10.1111/j.1574-6968.2011.02465.x | Fungal HGT review/catalogue | received |
| 33 | Raymond J et al. 2004 MBE, doi:10.1093/molbev/msh047 | Nitrogenase natural history (prokaryotic HGT) | received |
| 34 | Hall BG, Barlow M 2004 Drug Resist Updat, doi:10.1016/j.drup.2004.02.003 | Serine β-lactamase phylogeny | received |
| 35 | Brooke NM et al. 1998 Nature 392:920, doi:10.1038/31933 | ParaHox cluster | received |
| 36 | Robinson-Rechavi M et al. 2005 J Mol Evol, doi:10.1007/s00239-004-0175-8 | Nematode HNF4 nuclear-receptor burst | received |
| 37 | Gladyshev EA et al. 2008 Science, doi:10.1126/science.1156407 | Bdelloid rotifer HGT | received |
| 38 | Flot J-F et al. 2013 Nature, doi:10.1038/nature12326 | Bdelloid genome: degenerate tetraploidy + HGT | received |
| 39 | Moran NA, Jarvik T 2010 Science, doi:10.1126/science.1187113 | Fungus→aphid carotenoid HGT | received |
| 40 | Kusserow A et al. 2005 Nature, doi:10.1038/nature03158 | Wnt subfamilies in cnidarians | received |
| 41 | Larroux C et al. 2008 MBE, doi:10.1093/molbev/msn047 | Early metazoan TF classes | received |
| 42 | Haegeman A et al. 2011 MPMI, doi:10.1094/MPMI-03-11-0055 | Nematode HGT (cell-wall enzymes) | received |
| 43 | Leonard G, Richards TA 2012 PNAS — `1210909110_sd01.xls` (PMC3535628) | Gene fusion/fission across fungi | received |
| 44 | Brown CA et al. 2010 Curr Biol — 4 xls (PMC2877759) | Subtelomeric MAL/IMA families | received |
| 45 | Kollmar M, Mühlhausen S 2017 — per-sequence myosin class assignments (CyMoBase export) | Myosin class labels for 7,852 sequences | received |
| 46 | Sackton TB et al. 2007 Nat Genet, doi:10.1038/ng.2007.60 | Drosophila immune gene evolution (fast) | received |
| 47 | Dermauw W et al. 2020 IBMB, doi:10.1016/j.ibmb.2020.103490 | Arthropod P450 clans | received |
| 48 | Archibald JM et al. 1999 Curr Biol, doi:10.1016/s0960-9822(99)80457-6 | Chaperonin/CCT paralogy (archaea, LECA) | received |
| 49 | Figge RM et al. 1999 MBE, doi:10.1093/oxfordjournals.molbev.a026125 | GAPDH endosymbiotic gene transfer | received |
| 50 | Lurin C et al. 2004 Plant Cell, doi:10.1105/tpc.104.022236 | PPR family list (repeat proteins) | received |
