# Part C - Genome/proteome sources, species panels, species trees, compute

Author: agent C (genomes and species). Written 2026-10-26 (project START 2026-10-26T13:17Z). All "verified" statements were checked on 2026-10-26 against the live services named in the text (UniProt REST/FTP, EBI FTP, Ensembl REST/FTP, OMA REST API, NCBI Datasets REST, Europe PMC REST). Anything I could not check is marked **UNVERIFIED**. Raw metadata, scripts and logs are listed in section 8 (`literature/raw/C/`, `scripts/C_*`).

## 0. Executive summary and recommendations

1. **Primary proteome source: UniProt Reference Proteomes (RP) canonical FASTA** (one protein per gene by construction, CC BY 4.0, per-proteome ID-mapping files with RefSeq/Ensembl/GeneID/UniParc/OMA/OrthoDB cross-references, and every species the PI listed exists as a "Reference proteome"). For the 81 Quest-for-Orthologs (QfO) species use the **frozen QfO release 2026_02** (`QfO_release_2026_02.tar.gz`, 2.2 GB); for all other species use the RP files of the current UniProt release (2026_03) and record release + md5 in DATA.md. UniProt does not keep per-proteome archives of older RP releases (only whole-knowledgebase tarballs under `previous_releases/`), so **download once and pin checksums**.
2. **Ensembl/NCBI are secondary** (ID mapping, cross-checks, species missing from UniProt RP). **Phytozome is avoided** (JGI data-release restrictions, no redistribution). **YGOB is not scriptable** (web-only; its host refused connections from the cluster on ports 80 and 443; no documented bulk download) and is only useful as a literature reference for yeast synteny-curated homology. **OMA** is useful as an independent precomputed HOG resource and ID cross-reference (REST API works from the cluster), not as a proteome source (OMA contains all isoforms).
3. **Five panels** (86 unique UniProt reference proteomes in the core tiers; all 111 candidate proteome IDs verified as "Reference proteome" with matching taxid via UniProt REST):

| Panel | Scope | Core species | Core canonical proteins | Optional add-ons |
|---|---|---|---|---|
| V | Vertebrata (+3 invertebrate-deuterostome outgroups) | 14 | 333,315 | +7 (170,716) |
| Y | Saccharomycetaceae / WGD yeasts (+3 outgroups) | 16 | 85,181 | +4 (21,245) |
| P | Viridiplantae | 13 | 389,385 | +5 (187,402) |
| E | Eukaryota-wide (25 eukaryotes + 5 archaea + 4 bacteria) | 34 | 340,209 | +11 (231,643) |
| K | Prokaryotes (14 bacteria + 7 archaea) | 21 | 68,128 | none |

4. **Compute is the binding constraint for V, P, E** (each ~330-390k proteins). A crude upper-bound extrapolation from the team's measured OrthoFinder v3 run (section 5) gives tens of CPU-hours per panel; I recommend the PI run a two-point scaling calibration on V (e.g. 6 and 12 species) before committing, and keep "lite" variants ready (V-lite 256,842 proteins; P-lite 273,684; E-lite 286,588; recipes in section 3.6).
5. **Known data traps** (details in section 6): (a) the UniProt RP README statistics table lists 147,503 canonical entries for human but the FASTA has 20,652 (verified by streaming count; the other 51 non-QfO panel species match the table exactly); (b) UniProt REST `proteinCount` counts all entries (canonical + isoforms + non-gene-collapsed), not canonical genes; (c) assembly drift between releases (zebrafish RP switched from GCF_000002035.6 to GCF_049306965.1 in QfO 2026_02); (d) the yeast WGD and the vertebrate 2R are allopolyploid/hybrid events, so for those gene families the "species tree" is not a faithful description of gene histories.

---

## 1. Sources

### 1.1 Comparison table

| Source | Version / date (verified) | Licence / terms | Canonical one-protein-per-gene | Cross-reference / ID mapping | Notes for this project |
|---|---|---|---|---|---|
| **UniProt Reference Proteomes** (`ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/`) | Release 2026_03, 02-Sep-2026 (REST header `x-uniprot-release`); RP files updated every 8 weeks (RP README); 37,190 species in RP 2026_03 (RP README) | CC BY 4.0 for all copyrightable parts of the databases; free redistribution with the notice (text of `ftp.uniprot.org/.../complete/LICENSE`, verified) | **Yes**: `UPID_taxid.fasta.gz` = canonical sequences, one per gene; `_additional.fasta.gz` = isoforms with header "Isoform of ..."; `.gene2acc` maps gene -> accession(s), col 3 = canonical gene symbol | per-proteome `.idmapping.gz` (3 columns: accession, DB, ID). For S. cerevisiae (UP000002311) it contains e.g. RefSeq 5,950; GeneID 5,884; EnsemblGenome 6,073; OMA 5,423; OrthoDB 5,814; eggNOG 5,477; SGD; STRING; KEGG; UniParc 6,222 (counted by streaming the file). Also `_DNA.fasta` (CDS) | Main source. Proteome IDs UPxxxxxxxxx are stable names for a *proteome record*, but the assembly behind a species' RP can change (see 6.3). |
| **QfO reference proteomes** (`ftp.ebi.ac.uk/pub/databases/reference_proteomes/QfO/`) | Release **2026_02** (10-Jun-2026): UniProt 2026_02, Ensembl 115, Ensembl Genomes 62 (QfO README); yearly; earlier: 2023_03 (79 species), 2024_02, 2025_04 (81 each) under `previous_releases/` | No separate licence file in the QfO directory (verified listing); content is UniProt (CC BY 4.0) plus Ensembl-derived sequences (EMBL-EBI terms: "no additional restrictions", attribution expected) | Yes, same format as RP (canonical `.fasta`, `_additional.fasta`, `gene2acc`, `idmapping`, `_DNA.fasta`); plus `UPKB_PROTEIN_GENE_MAPPING` (columns observed: UPID, UniProt acc, EMBL/Ensembl protein id, genome acc, "-", gene name; content only partly inspected) | same | 81 species; frozen and archived: use for the QfO species. |
| **Ensembl** (+ Ensembl Genomes, now served from the same REST) | Release **116**, EG version 63 (REST `/info/data`, `/info/eg_version`); species per division via REST: Vertebrates 356, Plants 267, Fungi 1,499, Protists 236, Metazoa 386 (Bacteria 31,332) | EMBL-EBI Terms of Use: "no additional restrictions on the use or redistribution ... other than those provided by the original data owners"; attribution expected; third-party rights remain the user's responsibility (verified text) | Not in FASTA. `pep.all.fa.gz` contains **all** transcripts' proteins (header `gene:... transcript:...`). Canonical = transcript carrying the `Ensembl_canonical` tag in the GFF3 (`current_gff3/<species>/*.gff3.gz`) (MANE Select if one exists, else Ensembl algorithm; GENCODE tag description). Verified on Ciona intestinalis KH.116: 16,684 transcripts tagged vs 16,671 gene features | UniProt ID mapping `Ensembl_Protein`, `Ensembl_Transcript`, `Ensembl` (gene) (verified in `configure/idmapping/fields`); IDs carry version suffixes (ENSP00000269305.x) | Good for vertebrates/plants/fungi; archive releases exist (e.g. `jun2026.archive.ensembl.org`). Not available: *Xenopus laevis* (no Ensembl entry for taxid 8355), amphioxus (taxid 7739) not in any Ensembl division (exact-taxid match). Ensembl web pages were unreachable at fetch time (platform-transition notice), FTP/REST fine. |
| **NCBI RefSeq / Datasets** | Datasets API v2, CLI/API version 18.39.0 (verified); RefSeq annotation releases are named in the assembly report (e.g. `GCF_048934315.1-RS_2025_08` for lamprey) | NCBI places no restrictions on use or distribution of the data; submitters may nonetheless claim IP (NCBI policies page, verified) | Not directly. `protein.faa` of a genome package contains all isoforms (XP_/NP_). One per gene: join GFF3 gene -> mRNA -> CDS and keep the longest (NCBI Datasets tutorial "Get one protein per gene", found via search, tutorial URL `ncbi.nlm.nih.gov/datasets/docs/v2/tutorials/ortholog-get-one-isoform/`); the project env also ships `ncbi_primary_transcript` and `primary_transcript` (`conda_envs/ortho/bin/tools/`) | UniProt ID mapping `RefSeq_Protein`, `RefSeq_Nucleotide`, `GeneID`, `EMBL-GenBank-DDBJ_CDS` (verified, see 1.3); versioned accessions (NP_000537.3) | Best fallback for species not in UniProt RP or Ensembl (e.g. gymnosperms: RefSeq has *Pinus tabuliformis* GCF_031772625.1 with 72,033 protein-coding genes and *Cryptomeria japonica* GCF_030272615.1 with 41,227, via Datasets taxon query "Acrogymnospermae"). |
| **YGOB** (ygob.ucd.ie, Byrne & Wolfe) | "version 7" mentioned in secondary sources (UNVERIFIED which species set); cannot be reached from the cluster | No licence statement found; third-party catalogue says "not available for download" (UNVERIFIED) | Curated "pillars" (gene-order based homology) - web browser only | none | `curl` to ygob.ucd.ie: connection refused on 80/443 (verified); the Wayback copy of the home page is a Flash/Ruffle shell (no data). Use POInT-browse (BMC Bioinformatics 2023, `wgd.statgen.ncsu.edu`, orthology for 12 paleopolyploidies; site not tested) and Génolevures/Shen 2018 data (below) instead. Proteomes for YGOB species are in UniProt RP (all 16 core Y species verified). |
| **Phytozome / JGI** | not queried (login-walled); policy text verified from a mirror of the original Phytozome `DataReleasePolicy.txt` | **Restricted (Fort Lauderdale agreement)**: users agree not to publish "Reserved Analyses" (whole-genome gene-family/gene-set/chromosome-scale comparisons) before the genome paper is accepted; "assembly and sequence data should not be redistributed or repackaged without permission from the JGI"; acknowledgement required | n/a (not used) | n/a | **Do not commit or redistribute Phytozome-derived files.** For plants, UniProt RP and Ensembl Plants carry the same genomes (e.g. Chlamydomonas, Selaginella, Physcomitrium, soybean) under open terms. Whether a given Phytozome genome is still embargoed was not checked. Bulk-download access mode (login/Globus) UNVERIFIED. |
| **OMA browser** | **All.May2026**, API v1.11 (`/api/version/`); **2,983 genomes** (`/api/genome/`) | CC BY-SA 2.5 for OMA's own predictions, API and website per the OMA licence page (seen as a web-search snippet; page itself blocked by a Cloudflare challenge from the cluster, so **UNVERIFIED directly**); underlying sequences may carry third-party terms. Note share-alike. | OMA keeps all splice variants (human `nr_entries` 103,200; zebrafish 52,089); each protein record has `is_main_isoform`, `sequence_md5`, `canonicalid` (UniProt-style ID), `oma_hog_id`, chromosome/locus (verified on YEAST) | `/api/xref/?search=...` maps external accessions (UniProt, RefSeq, EMBL...) to OMA entries (verified: P04637 query returns hits); `/api/taxonomy/?type=newick&members=...` returns a taxonomy-based subtree (verified) | FastOMA needs the OMAmer database (`LUCA.h5`, already on disk at `data/raw/omamer/`). Precomputed OMA HOGs (`/api/hog/`) could serve as a third, independent comparator on the QfO/OMA species (not used for truth). |

### 1.2 How to obtain canonical sets (recipes)

* **UniProt RP (preferred).** For each UPID in a panel: `Eukaryota|Bacteria|Archaea/UPID/UPID_TAXID.fasta.gz` (+ `.gene2acc.gz`, `.idmapping.gz`). Headers are UniProt-style (`>sp|P04637|P53_HUMAN ... OX=9606 GN=TP53`). The canonical set was verified by streaming-count of the FASTA for the 51 panel proteomes that are not in the QfO release plus human (nothing stored; `literature/raw/C/fasta_canonical_counts.tsv`). Do not use `rest.uniprot.org/uniprotkb/stream?query=proteome:UPID`: that returns every UniProtKB entry of the proteome (human: 147,520 total; 20,416 reviewed) rather than one protein per gene.
* **QfO**: unpack `QfO_release_2026_02.tar.gz`; `Eukaryota|Bacteria|Archaea/*.fasta` are the canonical sets (51 + 23 + 7 species). `_additional.fasta` (20 GB tarball) is not needed.
* **Ensembl**: download `pep.all.fa.gz` + `gff3.gz`; keep proteins of transcripts with tag `Ensembl_canonical`; or run the OrthoFinder helper `primary_transcript.py` (longest per gene, understands Ensembl headers). Ensembl canonical != UniProt canonical, so literature accessions (typically UniProt or RefSeq) must be mapped through gene, not through protein sequence identity.
* **NCBI**: `datasets download genome accession GCF_... --include protein,gff3`, then longest CDS per gene (`ncbi_primary_transcript`).
* **Isoform caveat for the benchmark.** A curated literature gene may correspond to a non-canonical isoform in the proteome. Resolve literature accession -> UniProt accession (ID mapping) -> gene (`gene2acc`) -> canonical protein of that gene in the pinned proteome, and keep the original accession in the mapping table.

### 1.3 Cross-reference services (tested)

UniProt ID mapping REST pattern (client in `scripts/C_uniprot_idmap.py`): `POST https://rest.uniprot.org/idmapping/run` (fields `from`, `to`, `ids`) -> `GET /idmapping/status/{jobId}` -> `GET /idmapping/uniprotkb/results/stream/{jobId}?format=tsv&fields=...`. Verified calls (2026-10-26):

| from | ids | result |
|---|---|---|
| `RefSeq_Protein` | NP_001230.1, NP_000537.3, NP_033049.2 | CCNH P51946; TP53 -> P04637 **plus two unreviewed duplicates K7PPA8, Q53GA5** (one-to-many); NP_033049.2 -> P12672 |
| `Ensembl_Protein` | ENSP00000269305 | P04637 (an old zebrafish ENSDARP id returned nothing: retired/versioned IDs fail) |
| `EMBL-GenBank-DDBJ_CDS` | AAA32667.1 | P27760 |

Available source/target names include `RefSeq_Protein`, `RefSeq_Nucleotide`, `Ensembl`, `Ensembl_Protein`, `Ensembl_Transcript`, `Ensembl_Genomes*`, `EMBL-GenBank-DDBJ(_CDS)`, `GeneID`, `Gene_Name`, `UniParc`, `UniRef50/90/100`, `HGNC`, `SGD`, `PomBase`, `ZFIN`, `Xenbase` (from `configure/idmapping/fields`). Practical rules: (i) always filter mapping results to the accession set of the pinned proteome (many-to-many; unreviewed duplicates); (ii) strip/record version suffixes; (iii) store the UniParc ID (UPI...) and sequence md5 of each reference sequence - UniParc is a permanent sequence-level identifier and survives accession merges/demerges; OMA exposes `sequence_md5`; (iv) when an accession is obsolete, resolve through UniParc/UniProt secondary accessions rather than by gene name.

---

## 2. Quest for Orthologs (QfO) reference proteome set

* **Latest release: QfO 2026_02** (README: "Based on UniProt Release 2026_02, Ensembl release 115 and Ensembl Genome release 62", dated 10-Jun-2026; directory timestamp 2026-06-26). URL: `https://ftp.ebi.ac.uk/pub/databases/reference_proteomes/QfO/` (files: `README`, `STATS`, `PROTEOME_ASSEMBLY_MAPPING`, `UPKB_PROTEIN_GENE_MAPPING` (249 MB), `QfO_release_2026_02.tar.gz` (2.2 GB), `QfO_release_2026_02_additional.tar.gz` (20 GB), folders `Archaea/ Bacteria/ Eukaryota/`). Earlier releases: `previous_releases/qfo_release-2023_03` (79 species), `-2024_02`, `-2025_04` (81 each).
* **81 species**: 51 eukaryotes (976,393 canonical proteins), 23 bacteria (82,549), 7 archaea (17,317); total 1,076,259 canonical sequences (computed from the README table; `literature/raw/C/QfO_2026_02_species.tsv`). Licence: see 1.1 (UniProt CC BY 4.0 / EMBL-EBI terms).
* Full list (canonical = main FASTA entries; isoforms = `_additional`). Counts are QfO 2026_02 values:

| UniProt proteome | TaxID | Code | Group | Canonical | Isoforms | Species |
|---|---|---|---|---|---|---|
| UP000000554 | 64091 | HALSA | Arc | 2,427 | 0 | Halobacterium salinarum |
| UP000001686 | 374847 | KORCO | Arc | 1,602 | 0 | Korarchaeum cryptofilum |
| UP000000805 | 243232 | METJA | Arc | 1,786 | 0 | Methanocaldococcus jannaschii |
| UP000002487 | 188937 | METAC | Arc | 4,468 | 0 | Methanosarcina acetivorans |
| UP000000792 | 436308 | NITMS | Arc | 1,795 | 0 | Nitrosopumilus maritimus |
| UP000001974 | 273057 | SACS2 | Arc | 2,938 | 0 | Saccharolobus solfataricus |
| UP000000536 | 69014 | THEKO | Arc | 2,301 | 0 | Thermococcus kodakarensis |
| UP000000798 | 224324 | AQUAE | Bac | 1,553 | 0 | Aquifex aeolicus |
| UP000001570 | 224308 | BACSU | Bac | 4,288 | 16 | Bacillus subtilis |
| UP000001414 | 226186 | BACTN | Bac | 4,782 | 0 | Bacteroides thetaiotaomicron |
| UP000002526 | 224911 | BRADU | Bac | 8,253 | 0 | Bradyrhizobium diazoefficiens |
| UP000000431 | 272561 | CHLTR | Bac | 895 | 0 | Chlamydia trachomatis serovar D |
| UP000002008 | 324602 | CHLAA | Bac | 3,850 | 0 | Chloroflexus aurantiacus |
| UP000002524 | 243230 | DEIRA | Bac | 3,085 | 0 | Deinococcus radiodurans |
| UP000007719 | 515635 | DICTD | Bac | 1,743 | 0 | Dictyoglomus turgidum |
| UP000000625 | 83333 | ECOLI | Bac | 4,403 | 23 | Escherichia coli |
| UP000241660 | 190304 | FUSNN | Bac | 1,962 | 0 | Fusobacterium nucleatum subsp. nucleatum |
| UP000000577 | 243231 | GEOSL | Bac | 3,402 | 0 | Geobacter sulfurreducens |
| UP000000557 | 251221 | GLOVI | Bac | 4,406 | 0 | Gloeobacter violaceus |
| UP000000429 | 85962 | HELPY | Bac | 1,554 | 0 | Helicobacter pylori |
| UP000001408 | 189518 | LEPIN | Bac | 3,676 | 0 | Leptospira interrogans serogroup Icterohaemorrhagiae serovar Lai |
| UP000001584 | 83332 | MYCTU | Bac | 3,997 | 7 | Mycobacterium tuberculosis |
| UP000000807 | 243273 | MYCGE | Bac | 483 | 2 | Mycoplasmoides genitalium |
| UP000000425 | 122586 | NEIMB | Bac | 2,003 | 0 | Neisseria meningitidis serogroup B |
| UP000002438 | 208964 | PSEAE | Bac | 5,563 | 2 | Pseudomonas aeruginosa |
| UP000001025 | 243090 | RHOBA | Bac | 7,271 | 0 | Rhodopirellula baltica |
| UP000001973 | 100226 | STRCO | Bac | 8,038 | 2 | Streptomyces coelicolor |
| UP000001425 | 1111708 | SYNY3 | Bac | 3,508 | 2 | Synechocystis sp. |
| UP000000718 | 289376 | THEYD | Bac | 1,982 | 0 | Thermodesulfovibrio yellowstonii |
| UP000008183 | 243274 | THEMA | Bac | 1,852 | 0 | Thermotoga maritima |
| UP000007062 | 7165 | ANOGA | Euk | 13,029 | 1,540 | Anopheles gambiae |
| UP000006548 | 3702 | ARATH | Euk | 27,496 | 14,101 | Arabidopsis thaliana |
| UP000002530 | 330879 | ASPFU | Euk | 9,647 | 1 | Aspergillus fumigatus |
| UP000007241 | 684364 | BATDJ | Euk | 8,610 | 0 | Batrachochytrium dendrobatidis |
| UP000009136 | 9913 | BOVIN | Euk | 22,714 | 29,276 | Bos taurus |
| UP000001554 | 7739 | BRAFL | Euk | 26,616 | 12,034 | Branchiostoma floridae |
| UP000001940 | 6239 | CAEEL | Euk | 19,789 | 8,642 | Caenorhabditis elegans |
| UP000000559 | 237561 | CANAL | Euk | 6,036 | 1 | Candida albicans |
| UP000805418 | 9615 | CANLF | Euk | 20,990 | 22,678 | Canis lupus familiaris |
| UP000006906 | 3055 | CHLRE | Euk | 17,614 | 1,217 | Chlamydomonas reinhardtii |
| UP000008144 | 7719 | CIOIN | Euk | 16,680 | 631 | Ciona intestinalis |
| UP000002149 | 214684 | CRYD1 | Euk | 6,603 | 140 | Cryptococcus deneoformans |
| UP000000437 | 7955 | DANRE | Euk | 26,677 | 42,446 | Danio rerio |
| UP000076858 | 35525 | None | Euk | 26,600 | 0 | Daphnia magna |
| UP000002195 | 44689 | DICDI | Euk | 12,718 | 28 | Dictyostelium discoideum |
| UP000000803 | 7227 | DROME | Euk | 13,817 | 9,669 | Drosophila melanogaster |
| UP000000539 | 9031 | CHICK | Euk | 18,372 | 25,597 | Gallus gallus |
| UP000001548 | 184922 | GIAIC | Euk | 4,900 | 0 | Giardia intestinalis |
| UP000001519 | 9595 | GORGO | Euk | 21,788 | 22,938 | Gorilla gorilla gorilla |
| UP000015101 | 6412 | HELRO | Euk | 23,328 | 0 | Helobdella robusta |
| UP000005640 | 9606 | HUMAN | Euk | 20,652 | 148,985 | Homo sapiens |
| UP001735580 | 6945 | IXOSC | Euk | 25,977 | 5,595 | Ixodes scapularis |
| UP000000542 | 5664 | LEIMA | Euk | 8,038 | 1 | Leishmania major |
| UP000018468 | 7918 | LEPOC | Euk | 18,323 | 4,140 | Lepisosteus oculatus |
| UP000006718 | 9544 | MACMU | Euk | 21,893 | 22,523 | Macaca mulatta |
| UP000002280 | 13616 | MONDO | Euk | 21,225 | 14,996 | Monodelphis domestica |
| UP000001357 | 81824 | MONBE | Euk | 9,156 | 0 | Monosiga brevicollis |
| UP000000589 | 10090 | MOUSE | Euk | 21,853 | 41,471 | Mus musculus |
| UP000000561 | 5270 | MYCMD | Euk | 6,788 | 18 | Mycosarcoma maydis |
| UP000001593 | 45351 | NEMVE | Euk | 24,428 | 17 | Nematostella vectensis |
| UP000001805 | 367110 | NEUCR | Euk | 9,759 | 507 | Neurospora crassa |
| UP000059680 | 39947 | ORYSJ | Euk | 43,678 | 5,545 | Oryza sativa subsp. japonica |
| UP000001038 | 8090 | ORYLA | Euk | 23,617 | 12,521 | Oryzias latipes |
| UP000002277 | 9598 | PANTR | Euk | 23,014 | 25,780 | Pan troglodytes |
| UP000000600 | 5888 | PARTE | Euk | 39,350 | 0 | Paramecium tetraurelia |
| UP000663193 | 321614 | PHANO | Euk | 17,429 | 301 | Phaeosphaeria nodorum |
| UP000006727 | 3218 | PHYPA | Euk | 31,356 | 16,426 | Physcomitrium patens |
| UP000005238 | 164328 | PHYRM | Euk | 15,350 | 0 | Phytophthora ramorum |
| UP000001450 | 36329 | PLAF7 | Euk | 5,361 | 8 | Plasmodium falciparum |
| UP000008783 | 418459 | PUCGT | Euk | 15,688 | 120 | Puccinia graminis f. sp. tritici |
| UP000002494 | 10116 | RAT | Euk | 21,545 | 32,949 | Rattus norvegicus |
| UP000002311 | 559292 | YEAST | Euk | 6,066 | 32 | Saccharomyces cerevisiae |
| UP000002485 | 284812 | SCHPO | Euk | 5,129 | 18 | Schizosaccharomyces pombe |
| UP000001312 | 665079 | SCLS1 | Euk | 14,445 | 0 | Sclerotinia sclerotiorum |
| UP000001449 | 35128 | THAPS | Euk | 11,735 | 18 | Thalassiosira pseudonana |
| UP000007266 | 7070 | TRICA | Euk | 16,568 | 1,937 | Tribolium castaneum |
| UP000001542 | 412133 | TRIV3 | Euk | 50,190 | 0 | Trichomonas vaginalis |
| UP000186698 | 8355 | XENLA | Euk | 36,352 | 25,722 | Xenopus laevis |
| UP000008143 | 8364 | XENTR | Euk | 21,609 | 16,114 | Xenopus tropicalis |
| UP000001300 | 284591 | YARLI | Euk | 6,449 | 5 | Yarrowia lipolytica |
| UP000007305 | 4577 | MAIZE | Euk | 39,346 | 23,916 | Zea mays |

### 2.1 Coverage relative to our needs

| Clade needed | In QfO 2026_02 | Missing in QfO but available as UniProt Reference Proteome (verified) |
|---|---|---|
| **Vertebrate WGDs** | gar (LEPOC), zebrafish (DANRE), medaka (ORYLA), *X. laevis* (XENLA, 36,352 incl. L+S), *X. tropicalis* (XENTR), chicken, 9 mammals, amphioxus (BRAFL), Ciona | lamprey *Petromyzon marinus* UP001735300; hagfish *Eptatretus burgeri* UP000694388; elephant shark *Callorhinchus milii* UP000314986 (listed as "Ghost shark"); coelacanth UP000008672; reedfish UP000694620; Atlantic salmon UP001652741 (also coho, rainbow trout UP000694395, chinook, brown trout, lake trout); common carp UP000694427; stickleback; fugu; anole; platypus; sea urchin |
| **Yeast WGD** | only *S. cerevisiae*, *C. albicans*, *Y. lipolytica*, *S. pombe* (and Neurospora, Aspergillus etc.) - **no** *C. glabrata*, *N. castellii*, *K. lactis*, *Lachancea*, *Z. rouxii*, *E. gossypii* | all of them: CANGA UP000002428, NAUCC UP000001640, KLULA UP000000598, LACTH UP000002036, ZYGRC UP000008536, ASHGO UP000000591 (*Lachancea waltii* has no UniProt proteome) |
| **Plants** | *A. thaliana*, *O. sativa*, *Z. mays*, *P. patens*, *C. reinhardtii* (5) | *Marchantia polymorpha* UP000244005, *Selaginella* UP000001514, *Amborella* UP000017836, *Brassica rapa* (listed as *B. campestris*, taxid 3711) UP000011750, poplar UP000006729, grape UP000009183, soybean UP000008827, *Klebsormidium*, *Chara*, *Ostreococcus*, fern *Ceratopteris*; gymnosperm: only *Taxus chinensis* UP000824469 (44,066) in UniProt RP - no *Picea/Ginkgo/Pinus* (REST search returns only viruses); NCBI RefSeq has two conifers (1.1) |
| **Protists / supergroups** | Amoebozoa 1 (DICDI); SAR 4 (PARTE, PHYRM, PLAF7, THAPS); Excavata 3 (GIAIC, LEIMA, TRIV3); Opisthokonta non-animal/fungal 1 (MONBE); no Archaeplastida algae besides CHLRE, no Haptista, Cryptista, Rhizaria, Apusozoa, Ancyromonads, Hemimastigophora | Amoebozoa: *Acanthamoeba*, *Entamoeba*, *Heterostelium*; SAR: *Tetrahymena*, *Toxoplasma*, *Phaeodactylum*, *Perkinsus*, *Ichthyophthirius*; Haptista: *Emiliania huxleyi*; Cryptista: *Guillardia theta*; Excavata: *Trypanosoma*, *Naegleria*, *Bodo*; Rhodophyta: *Cyanidioschyzon*, *Galdieria*, *Porphyra*, *Chondrus*; Apusozoa: *Thecamonas*; holozoan relatives: *Capsaspora*, *Salpingoeca*, *Sphaeroforma*; non-bilaterians: *Amphimedon*, *Mnemiopsis*, *Trichoplax*, *Hydra*. **No UniProt proteome** at all (REST organism-name search empty): *Cyanophora* (Glaucophyta), *Euglena*, *Mesostigma*, *Chlorokybus*, *Anthoceros*, *Ginkgo*, *Picea*; *Bigelowiella natans* has only a 283-protein (nucleomorph) proteome |
| **Archaea** | 7 (Euryarchaeota x4, Crenarchaeota x1, Thaumarchaeota x1, Korarchaeota x1) | Asgard: *Candidatus Lokiarchaeum ossiferum* UP001208689 (5,112); DPANN: *Nanoarchaeum equitans*, *Micrarchaeum*; many more (RP has hundreds of archaea) |
| **Bacteria** | 23 over about 16 phyla (Proteobacteria alpha/beta/gamma/delta/epsilon, Firmicutes, Actinobacteria, Chlamydiae, Planctomycetes, Cyanobacteria x2, Bacteroidetes, Spirochaetes, Fusobacteria, Deinococcus-Thermus, Thermotogae, Aquificae, Chloroflexi, Dictyoglomi, Nitrospirae, Tenericutes) | essentially unlimited (RP has thousands); CPR/Patescibacteria are not in QfO |

QfO strengths: curated, frozen, widely used by existing benchmarks and by FastOMA/OrthoFinder papers; gaps: yeast WGD relatives, most plant lineages, non-mammalian vertebrates, Archaeplastida/Haptista/Rhizaria/Cryptista, Asgard.

---

## 3. Proposed species panels

Selection rules: every proteome is a UniProt "Reference proteome" with taxid identical to the intended species (all 111 candidate IDs verified by `rest.uniprot.org/proteomes/{UPID}`); canonical protein counts are QfO 2026_02 values for QfO species and RP 2026_03 table values otherwise, the latter validated by streaming the FASTA (51/51 non-QfO panel species match; human taken from the QfO set, see 6.1). Species were chosen for (a) the polyploidy/duplication events the literature families will probe, (b) at least two outgroups per panel and a range of evolutionary rates, (c) preference for compact proteomes and against gene-count inflation (e.g. *Trichomonas*, *Paramecium*, *Tetrahymena*, *Emiliania* and the large polyploid crops are optional), (d) feasible compute. "Set" = QfO (frozen 2026_02) or RP (current UniProt RP 2026_03). "OMA code" and "Ensembl" columns are exact-taxid matches only (strain-level taxids make "-" a weak negative, e.g. *K. lactis*).


#### Panel V: Vertebrata

**Core: 14 proteomes, 333,315 canonical proteins.**

| Code | Species | TaxID | UniProt proteome | Set | Canonical proteins | OMA code (exact taxid) | Ensembl (exact taxid) | Why / note |
|---|---|---|---|---|---|---|---|---|
| BRAFL | Branchiostoma floridae | 7739 | UP000001554 | QfO | 26,616 | BRAFL | - | cephalochordate outgroup (pre-1R/2R) |
| CIOIN | Ciona intestinalis | 7719 | UP000008144 | QfO | 16,680 | CIOIN | Y | tunicate outgroup (pre-2R; fast evolving) |
| PETMA | Petromyzon marinus | 7757 | UP001735300 | RP | 21,315 | PETMA | Y | lamprey (cyclostome; 1R/2R timing) |
| CALMI | Callorhinchus milii | 7868 | UP000314986 | RP | 19,409 | CALMI | Y | elephant/ghost shark (chondrichthyan; slow evolving) |
| LEPOC | Lepisosteus oculatus | 7918 | UP000018468 | QfO | 18,323 | LEPOC | Y | spotted gar (pre-3R actinopterygian) |
| DANRE | Danio rerio | 7955 | UP000000437 | QfO | 26,677 | DANRE | Y | zebrafish (3R teleost) |
| ORYLA | Oryzias latipes | 8090 | UP000001038 | QfO | 23,617 | ORYLA | Y | medaka (3R teleost) |
| SALSA | Salmo salar | 8030 | UP001652741 | RP | 42,272 | SALSA | Y | Atlantic salmon (4R salmonid) |
| LATCH | Latimeria chalumnae | 7897 | UP000008672 | RP | 19,568 | LATCH | Y | coelacanth (sarcopterygian outgroup) |
| XENTR | Xenopus tropicalis | 8364 | UP000008143 | QfO | 21,609 | XENTR | Y | diploid frog |
| XENLA | Xenopus laevis | 8355 | UP000186698 | QfO | 36,352 | XENLA | - | allotetraploid frog (L/S subgenomes) |
| CHICK | Gallus gallus | 9031 | UP000000539 | QfO | 18,372 | CHICK | Y | bird |
| HUMAN | Homo sapiens | 9606 | UP000005640 | QfO | 20,652 | HUMAN | Y | human |
| MOUSE | Mus musculus | 10090 | UP000000589 | QfO | 21,853 | MOUSE | Y | mouse |

**Optional: 7 proteomes, 170,716 canonical proteins.**

| Code | Species | TaxID | UniProt proteome | Set | Canonical proteins | OMA code (exact taxid) | Ensembl (exact taxid) | Why / note |
|---|---|---|---|---|---|---|---|---|
| STRPU | Strongylocentrotus purpuratus | 7668 | UP000007110 | RP | 26,461 | STRPU | Y | echinoderm outgroup |
| EPTBU | Eptatretus burgeri | 7764 | UP000694388 | RP | 16,495 | EPTBU | Y | hagfish (cyclostome) |
| ERPCA | Erpetoichthys calabaricus | 27687 | UP000694620 | RP | 22,211 | - | Y | reedfish (Polypteriformes; basal actinopterygian) |
| CYPCA | Cyprinus carpio | 7962 | UP000694427 | RP | 44,401 | - | Y | common carp (4R cyprinid) |
| GASAC | Gasterosteus aculeatus | 481459 | UP000007635 | RP | 22,263 | - | Y | stickleback (3R teleost) |
| ANOCA | Anolis carolinensis | 28377 | UP000001646 | RP | 21,494 | ANOCA | Y | squamate |
| ORNAN | Ornithorhynchus anatinus | 9258 | UP000002279 | RP | 17,391 | ORNAN | Y | monotreme |


#### Panel Y: Saccharomycetaceae / budding yeasts

**Core: 16 proteomes, 85,181 canonical proteins.**

| Code | Species | TaxID | UniProt proteome | Set | Canonical proteins | OMA code (exact taxid) | Ensembl (exact taxid) | Why / note |
|---|---|---|---|---|---|---|---|---|
| YEAST | Saccharomyces cerevisiae S288c | 559292 | UP000002311 | QfO | 6,066 | YEAST | Y | post-WGD reference |
| CANGA | Candida glabrata (Nakaseomyces glabratus) | 284593 | UP000002428 | RP | 5,198 | CANGA | Y | post-WGD |
| NAUCC | Naumovozyma castellii | 27288 | UP000001640 | RP | 5,565 | - | Y | post-WGD |
| NAUDC | Naumovozyma dairenensis | 1071378 | UP000000689 | RP | 5,523 | NAUDC | Y | post-WGD |
| KAZAF | Kazachstania africana | 1071382 | UP000005220 | RP | 5,354 | KAZAF | Y | post-WGD |
| VANPO | Vanderwaltozyma polyspora | 436907 | UP000000267 | RP | 5,314 | VANPO | Y | post-WGD (early diverging) |
| TETPH | Tetrapisispora phaffii | 1071381 | UP000005666 | RP | 5,235 | TETPH | Y | post-WGD (early diverging) |
| ZYGRC | Zygosaccharomyces rouxii | 559307 | UP000008536 | RP | 4,987 | - | - | pre-WGD (ZT clade) |
| TORDC | Torulaspora delbrueckii | 1076872 | UP000005627 | RP | 4,974 | - | - | pre-WGD (ZT clade) |
| LACTH | Lachancea thermotolerans | 559295 | UP000002036 | RP | 5,093 | LACTC | Y | pre-WGD (Lachancea; unstable node) |
| KLULA | Kluyveromyces lactis | 284590 | UP000000598 | RP | 5,071 | KLULA | - | pre-WGD |
| ASHGO | Eremothecium gossypii | 284811 | UP000000591 | RP | 4,758 | ASHGO | - | pre-WGD (reduced genome) |
| ERECY | Eremothecium cymbalariae | 931890 | UP000006790 | RP | 4,429 | ERECY | Y | pre-WGD |
| CANAL | Candida albicans SC5314 | 237561 | UP000000559 | QfO | 6,036 | CANAL | Y | CTG-clade outgroup |
| YARLI | Yarrowia lipolytica | 284591 | UP000001300 | QfO | 6,449 | YARLI | Y | distant yeast outgroup |
| SCHPO | Schizosaccharomyces pombe | 284812 | UP000002485 | QfO | 5,129 | SCHPO | Y | Taphrinomycotina outgroup |

**Optional: 4 proteomes, 21,245 canonical proteins.**

| Code | Species | TaxID | UniProt proteome | Set | Canonical proteins | OMA code (exact taxid) | Ensembl (exact taxid) | Why / note |
|---|---|---|---|---|---|---|---|---|
| HENBL | Henningerozyma blattae (Tetrapisispora blattae) | 1071380 | UP000002866 | RP | 5,385 | TETBL | Y | post-WGD |
| KAZNA | Huiozyma naganishii (Kazachstania naganishii) | 1071383 | UP000006310 | RP | 5,301 | KAZNA | Y | post-WGD |
| SACUV | Saccharomyces uvarum | 230603 | UP001162085 | RP | 5,489 | - | - | shallow post-WGD (Saccharomyces sensu stricto) |
| KLUMA | Kluyveromyces marxianus | 4911 | UP000422736 | RP | 5,070 | - | Y | pre-WGD, KLE clade |


#### Panel P: Viridiplantae

**Core: 13 proteomes, 389,385 canonical proteins.**

| Code | Species | TaxID | UniProt proteome | Set | Canonical proteins | OMA code (exact taxid) | Ensembl (exact taxid) | Why / note |
|---|---|---|---|---|---|---|---|---|
| CHLRE | Chlamydomonas reinhardtii | 3055 | UP000006906 | QfO | 17,614 | CHLRE | Y | chlorophyte |
| OSTLU | Ostreococcus lucimarinus | 436017 | UP000001568 | RP | 7,403 | OSTLU | Y | prasinophyte (compact genome) |
| KLENI | Klebsormidium nitens | 105231 | UP000054558 | RP | 16,171 | KLENI | - | streptophyte alga |
| MARPO | Marchantia polymorpha | 3197 | UP000244005 | RP | 19,277 | MARPO | Y | liverwort (little paleopolyploidy) |
| PHYPA | Physcomitrium patens | 3218 | UP000006727 | QfO | 31,356 | PHYPA | Y | moss (WGDs) |
| SELMO | Selaginella moellendorffii | 88036 | UP000001514 | RP | 33,131 | SELML | Y | lycophyte |
| AMBTC | Amborella trichopoda | 13333 | UP000017836 | RP | 27,398 | AMBTC | Y | sister to other angiosperms (contested) |
| VITVI | Vitis vinifera | 29760 | UP000009183 | RP | 29,855 | VITVI | Y | gamma only (no lineage WGD) |
| ARATH | Arabidopsis thaliana | 3702 | UP000006548 | QfO | 27,496 | ARATH | Y | brassicaceae (alpha/beta WGD) |
| BRARP | Brassica rapa | 3711 | UP000011750 | RP | 40,808 | - | Y | mesohexaploid |
| SOYBN | Glycine max | 3847 | UP000008827 | RP | 55,852 | SOYBN | Y | legume WGD |
| ORYSJ | Oryza sativa subsp. japonica | 39947 | UP000059680 | QfO | 43,678 | ORYSJ | Y | monocot (rho WGD) |
| MAIZE | Zea mays | 4577 | UP000007305 | QfO | 39,346 | MAIZE | Y | maize tetraploidy |

**Optional: 5 proteomes, 187,402 canonical proteins.**

| Code | Species | TaxID | UniProt proteome | Set | Canonical proteins | OMA code (exact taxid) | Ensembl (exact taxid) | Why / note |
|---|---|---|---|---|---|---|---|---|
| CHABU | Chara braunii | 69332 | UP000265515 | RP | 34,645 | CHABU | Y | charophyte alga |
| CERRI | Ceratopteris richardii | 49495 | UP000825935 | RP | 36,181 | CERRI | - | fern (large genome; WGD) |
| AQUCO | Aquilegia coerulea | 218851 | UP000230069 | RP | 24,741 | - | - | basal eudicot (Ranunculales) |
| POPTR | Populus trichocarpa | 3694 | UP000006729 | RP | 57,605 | POPTR | Y | salicoid WGD |
| BRADI | Brachypodium distachyon | 15368 | UP000008810 | RP | 34,230 | BRADI | Y | grass |


#### Panel E: Eukaryota-wide (+ archaea/bacteria outgroups)

**Core: 34 proteomes, 340,209 canonical proteins.**

| Code | Species | TaxID | UniProt proteome | Set | Canonical proteins | OMA code (exact taxid) | Ensembl (exact taxid) | Why / note |
|---|---|---|---|---|---|---|---|---|
| HUMAN | Homo sapiens | 9606 | UP000005640 | QfO | 20,652 | HUMAN | Y | bilaterian/vertebrate anchor; richest literature |
| DROME | Drosophila melanogaster | 7227 | UP000000803 | QfO | 13,817 | DROME | Y | protostome (ecdysozoa) anchor |
| NEMVE | Nematostella vectensis | 45351 | UP000001593 | QfO | 24,428 | NEMVE | Y | cnidarian (slow-evolving; non-bilaterian) |
| AMPQE | Amphimedon queenslandica | 400682 | UP000007879 | RP | 20,174 | AMPQE | Y | sponge (Porifera) - node vs ctenophore |
| MNELE | Mnemiopsis leidyi | 27923 | UP001729289 | RP | 19,289 | MNELE | Y | ctenophore - contested animal root |
| TRIAD | Trichoplax adhaerens | 10228 | UP000009022 | RP | 11,518 | TRIAD | Y | placozoan (tiny, compact) |
| MONBE | Monosiga brevicollis | 81824 | UP000001357 | QfO | 9,156 | MONBE | - | choanoflagellate (animal sister) |
| CAPO3 | Capsaspora owczarzaki | 595528 | UP000008743 | RP | 8,738 | CAPO3 | Y | filasterean (holozoan outgroup) |
| YEAST | Saccharomyces cerevisiae | 559292 | UP000002311 | QfO | 6,066 | YEAST | Y | fungal anchor; post-WGD |
| NEUCR | Neurospora crassa | 367110 | UP000001805 | QfO | 9,759 | NEUCR | Y | filamentous ascomycete |
| BATDJ | Batrachochytrium dendrobatidis | 684364 | UP000007241 | QfO | 8,610 | BATDJ | Y | chytrid (early-diverging fungus) |
| ENCCU | Encephalitozoon cuniculi | 284813 | UP000000819 | RP | 2,041 | ENCCU | Y | microsporidian (reduced, fast-evolving) |
| DICDI | Dictyostelium discoideum | 44689 | UP000002195 | QfO | 12,718 | DICDI | - | amoebozoa (social amoeba) |
| ENTH1 | Entamoeba histolytica | 294381 | UP000001926 | RP | 7,967 | ENTH1 | Y | amoebozoa (reduced parasite) |
| ARATH | Arabidopsis thaliana | 3702 | UP000006548 | QfO | 27,496 | ARATH | Y | archaeplastida anchor; literature-rich |
| CHLRE | Chlamydomonas reinhardtii | 3055 | UP000006906 | QfO | 17,614 | CHLRE | Y | green alga |
| CYAME | Cyanidioschyzon merolae | 280699 | UP000007014 | RP | 5,018 | CYAM1 | Y | red alga (compact; Rhodophyta) |
| THAPS | Thalassiosira pseudonana | 35128 | UP000001449 | QfO | 11,735 | - | - | stramenopile (diatom) |
| PHYRM | Phytophthora ramorum | 164328 | UP000005238 | QfO | 15,350 | PHYRM | Y | stramenopile (oomycete) |
| PLAF7 | Plasmodium falciparum 3D7 | 36329 | UP000001450 | QfO | 5,361 | PLAF7 | Y | alveolate (apicomplexa; AT-rich, fast) |
| TRYB2 | Trypanosoma brucei | 185431 | UP000008524 | RP | 8,561 | TRYB2 | Y | excavate (kinetoplastid) |
| NAEGR | Naegleria gruberi | 5762 | UP000006671 | RP | 15,634 | - | Y | excavate (heterolobosean; free-living) |
| GIAIC | Giardia intestinalis | 184922 | UP000001548 | QfO | 4,900 | GIAIC | - | excavate (diplomonad; reduced) |
| THETR | Thecamonas trahens | 461836 | UP000054408 | RP | 10,533 | THETB | Y | apusozoan (opimoda; deep branching) |
| SACS2 | Saccharolobus solfataricus | 273057 | UP000001974 | QfO | 2,938 | - | - | archaeon (TACK/Crenarchaeota) |
| THEKO | Thermococcus kodakarensis | 69014 | UP000000536 | QfO | 2,301 | - | - | archaeon (Euryarchaeota) |
| NITMS | Nitrosopumilus maritimus | 436308 | UP000000792 | QfO | 1,795 | - | - | archaeon (Thaumarchaeota) |
| HALSA | Halobacterium salinarum | 64091 | UP000000554 | QfO | 2,427 | - | - | archaeon (halophile Euryarchaeota) |
| ECOLI | Escherichia coli K12 | 83333 | UP000000625 | QfO | 4,403 | - | - | bacterium (gamma-proteobacteria) |
| BACSU | Bacillus subtilis | 224308 | UP000001570 | QfO | 4,288 | - | - | bacterium (Firmicutes) |
| SYNY3 | Synechocystis sp. PCC 6803 | 1111708 | UP000001425 | QfO | 3,508 | - | - | bacterium (Cyanobacteria; plastid ancestry) |
| BRADU | Bradyrhizobium diazoefficiens | 224911 | UP000002526 | QfO | 8,253 | - | - | bacterium (alpha-proteobacteria; mitochondrial ancestry) |
| ICHMU | Ichthyophthirius multifiliis | 5932 | UP000008983 | RP | 8,049 | - | Y | ciliate (compact proteome alternative to Tetrahymena/Paramecium) |
| LOKOS | Candidatus Lokiarchaeum ossiferum | 2951803 | UP001208689 | RP | 5,112 | - | - | Asgard archaeon (closest archaeal relatives of eukaryotes) |

**Optional: 11 proteomes, 231,643 canonical proteins.**

| Code | Species | TaxID | UniProt proteome | Set | Canonical proteins | OMA code (exact taxid) | Ensembl (exact taxid) | Why / note |
|---|---|---|---|---|---|---|---|---|
| CAEEL | Caenorhabditis elegans | 6239 | UP000001940 | QfO | 19,789 | CAEEL | Y | nematode; fast-evolving |
| SCHPO | Schizosaccharomyces pombe | 284812 | UP000002485 | QfO | 5,129 | SCHPO | Y | fission yeast (Taphrinomycotina) |
| CRYD1 | Cryptococcus deneoformans | 214684 | UP000002149 | QfO | 6,603 | CRYNJ | Y | basidiomycete yeast |
| ACACA | Acanthamoeba castellanii | 1257118 | UP000011083 | RP | 14,938 | ACACA | Y | amoebozoa (free-living) |
| PHYPA | Physcomitrium patens | 3218 | UP000006727 | QfO | 31,356 | PHYPA | Y | moss |
| TOXGM | Toxoplasma gondii ME49 | 508771 | UP000001529 | RP | 8,307 | - | Y | alveolate (apicomplexa) |
| TETTS | Tetrahymena thermophila | 312017 | UP000009168 | RP | 26,975 | TETTS | Y | alveolate (ciliate; large proteome) |
| EMIHU | Emiliania huxleyi | 280463 | UP000013827 | RP | 35,603 | EMIHU | Y | haptophyte (large proteome) |
| GUITH | Guillardia theta | 905079 | UP000011087 | RP | 24,715 | GUITC | Y | cryptophyte |
| LEIMA | Leishmania major | 5664 | UP000000542 | QfO | 8,038 | - | - | excavate (kinetoplastid) |
| TRIV3 | Trichomonas vaginalis | 412133 | UP000001542 | QfO | 50,190 | TRIVA | - | excavate (parabasalid; 50k genes) |


#### Panel K: Prokaryotes

**Core: 21 proteomes, 68,128 canonical proteins.**

| Code | Species | TaxID | UniProt proteome | Set | Canonical proteins | OMA code (exact taxid) | Ensembl (exact taxid) | Why / note |
|---|---|---|---|---|---|---|---|---|
| ECOLI | Escherichia coli K12 | 83333 | UP000000625 | QfO | 4,403 | - | - | Gammaproteobacteria; best-annotated bacterium |
| PSEAE | Pseudomonas aeruginosa | 208964 | UP000002438 | QfO | 5,563 | - | - | Gammaproteobacteria (Pseudomonadales); large genome |
| BRADU | Bradyrhizobium diazoefficiens | 224911 | UP000002526 | QfO | 8,253 | - | - | Alphaproteobacteria (mitochondrial lineage); large genome |
| HELPY | Helicobacter pylori | 85962 | UP000000429 | QfO | 1,554 | - | - | Epsilonproteobacteria; small genome |
| GEOSL | Geobacter sulfurreducens | 243231 | UP000000577 | QfO | 3,402 | - | - | Deltaproteobacteria |
| BACSU | Bacillus subtilis | 224308 | UP000001570 | QfO | 4,288 | - | - | Firmicutes (Bacillota) |
| MYCTU | Mycobacterium tuberculosis | 83332 | UP000001584 | QfO | 3,997 | - | - | Actinobacteria |
| CHLTR | Chlamydia trachomatis | 272561 | UP000000431 | QfO | 895 | - | - | Chlamydiae; obligate intracellular, tiny genome |
| SYNY3 | Synechocystis sp. PCC 6803 | 1111708 | UP000001425 | QfO | 3,508 | - | - | Cyanobacteria (plastid lineage) |
| DEIRA | Deinococcus radiodurans | 243230 | UP000002524 | QfO | 3,085 | - | - | Deinococcota |
| THEMA | Thermotoga maritima | 243274 | UP000008183 | QfO | 1,852 | - | - | Thermotogota (extensive HGT with archaea) |
| AQUAE | Aquifex aeolicus | 224324 | UP000000798 | QfO | 1,553 | - | - | Aquificota |
| BACTN | Bacteroides thetaiotaomicron | 226186 | UP000001414 | QfO | 4,782 | - | - | Bacteroidota |
| LEPIN | Leptospira interrogans | 189518 | UP000001408 | QfO | 3,676 | - | - | Spirochaetota |
| SACS2 | Saccharolobus solfataricus | 273057 | UP000001974 | QfO | 2,938 | - | - | Crenarchaeota (Sulfolobales) |
| KORCO | Korarchaeum cryptofilum | 374847 | UP000001686 | QfO | 1,602 | - | - | Korarchaeota |
| NITMS | Nitrosopumilus maritimus | 436308 | UP000000792 | QfO | 1,795 | - | - | Thaumarchaeota (Nitrososphaeria) |
| THEKO | Thermococcus kodakarensis | 69014 | UP000000536 | QfO | 2,301 | - | - | Euryarchaeota (Thermococcales) |
| METJA | Methanocaldococcus jannaschii | 243232 | UP000000805 | QfO | 1,786 | - | - | Euryarchaeota (Methanococci) |
| METAC | Methanosarcina acetivorans | 188937 | UP000002487 | QfO | 4,468 | - | - | Euryarchaeota (Methanosarcinales; bacterial HGT) |
| HALSA | Halobacterium salinarum | 64091 | UP000000554 | QfO | 2,427 | - | - | Euryarchaeota (Halobacteria; HGT-rich) |



### 3.1 Panel V - Vertebrata (events: 1R, 2R, cyclostome-specific expansion, teleost 3R, salmonid 4R, Xenopus allotetraploidy)

* Core (14): outgroups BRAFL, CIOIN (pre-1R); cyclostome PETMA; chondrichthyan CALMI; pre-3R actinopterygian LEPOC (spotted gar; Braasch 2016); teleosts DANRE, ORYLA (3R, ~320 Ma per Lien 2016 abstract) and SALSA (Ss4R, ~80 Ma per Lien 2016; <= 88 Ma per Macqueen & Johnston 2014); LATCH; frogs XENTR vs XENLA (allotetraploid ~17-18 Ma, 56% genes retained as homeologs, Session 2016); CHICK, HUMAN, MOUSE as amniote anchors.
* Optional (7): STRPU (3rd invertebrate outgroup), EPTBU (second cyclostome - tests the lamprey-specific expansion), ERPCA (Polypteriformes, earliest actinopterygian branch), CYPCA (carp, a second 4R), GASAC (second 3R teleost), ANOCA (reptile), ORNAN (monotreme).
* Notes: SALSA UniProt canonical count (42,272) is close to the NCBI RS_2025_05 annotation (43,000 protein-coding genes, assembly Ssal_v3.1, GCF_905237065.1); XENLA RP (36,351) vs NCBI v10.1 annotation 34,476 protein-coding genes; lamprey RP is the new assembly `UKy_Petmar_22M1.pri1.0` (GCF_048934315.1; NCBI RS_2025_08, 22,381 protein-coding genes; whether it is a germline or somatic assembly is **UNVERIFIED** - relevant because the lamprey germline genome has genes eliminated from somatic cells, Smith 2018).

### 3.2 Panel Y - Saccharomycetaceae (events: the yeast WGD with its hybrid origin, KLE/ZT relationships, gene loss in reduced genomes)

* Core (16): post-WGD YEAST, CANGA, NAUCC, NAUDC, KAZAF, VANPO, TETPH (seven post-WGD genera spanning the clade); pre-WGD ZT clade ZYGRC, TORDC; KLE clade LACTH, KLULA, ASHGO, ERECY; outgroups CANAL (CTG clade), YARLI, SCHPO.
* Optional (4): HENBL, KAZNA, SACUV (shallow Saccharomyces comparison for recent-divergence behaviour), KLUMA.
* 85,181 proteins in total: this is cheap; if the PI wants a bigger yeast panel, Y can be extended by 10-20 more Saccharomycetaceae at little cost.
* All UniProt RP yeast proteomes here have 4,429-6,449 canonical proteins; ASHGO/ERECY are reduced-genome species (4,758/4,429).

### 3.3 Panel P - Viridiplantae (events: a seed-plant WGD and an angiosperm WGD (Jiao 2011), Brassica triplication, soybean palaeopolyploidy, maize tetraploidy, moss duplications, grape hexaploidization)

* Core (13): chlorophytes CHLRE, OSTLU (compact prasinophyte genome); streptophyte alga KLENI; liverwort MARPO (low-redundancy genome, "resembling the ancestral land plant", Bowman 2017); moss PHYPA (QfO); lycophyte SELMO; ANGIOSPERM backbone AMBTC (sister to all other angiosperms; one pre-angiosperm WGD and no lineage-specific WGD, Amborella 2013), VITVI (gamma only), ARATH (QfO), BRARP (mesohexaploid; 41,174 genes modelled, Wang 2011), SOYBN, ORYSJ and MAIZE (QfO).
* Optional (5): CHABU (*Chara*, 34,645 canonical; charophyte), CERRI (fern, 36,181), AQUCO (Ranunculales), POPTR (57,605 canonical!), BRADI.
* Gymnosperms: the only UniProt RP is *Taxus chinensis* (not placed in any panel; 44,066 proteins); RefSeq conifers listed in 1.1 are possible but heavy.
* The plant canonical counts are the largest in this survey (SOYBN 55,852; ORYSJ 43,678; BRARP 40,808; MAIZE 39,346; SELMO 33,131). They are what UniProt gene-collapses, i.e. they include low-confidence/TE-related models; I did not inspect whether these are inflated vs high-confidence gene sets (**UNVERIFIED**). If compute forces a cut, drop SOYBN, ORYSJ, KLENI first (P-lite).

### 3.4 Panel E - Eukaryota-wide with prokaryotic outgroups

* Core (34): Opisthokonta - animals HUMAN, DROME, NEMVE (cnidarian), AMPQE (sponge), MNELE (ctenophore), TRIAD (placozoan), unicellular holozoans MONBE, CAPO3; fungi YEAST, NEUCR, BATDJ (chytrid), ENCCU (microsporidian, fast-evolving); Amoebozoa DICDI, ENTH1; Archaeplastida ARATH, CHLRE, CYAME (red alga); SAR/stramenopile-alveolate THAPS, PHYRM, PLAF7, ICHMU (ciliate, compact alternative to *Tetrahymena*); Excavata TRYB2, NAEGR, GIAIC; Apusozoa THETR; archaea SACS2, THEKO, NITMS, HALSA, LOKOS (Asgard); bacteria ECOLI, BACSU, SYNY3 (cyanobacterial/plastid), BRADU (alphaproteobacterial/mitochondrial).
* Optional (11): SCHPO, CRYD1, ACACA, TOXGM, TETTS, TRIV3 (50,190 canonical), EMIHU (35,603), GUITH (24,715), LEIMA, CAEEL, PHYPA (these were left out of the core to hold it at ~340k proteins).
* The panel deliberately includes both ctenophore (MNELE) and sponge (AMPQE) so that the contested animal root is testable (section 4). Haptista (EMIHU) and Cryptista (GUITH) are optional because of proteome size; without them these supergroups are uncovered. Rhizaria cannot be covered (only a nucleomorph proteome for *Bigelowiella*).
* Note on the OMAmer/FastOMA root level: the E panel is the one on which FastOMA's LUCA-level database matters most.

### 3.5 Panel K - Prokaryotes (14 bacteria + 7 archaea, 68,128 proteins, all QfO)

Covers Gamma/Alpha/Epsilon/Delta-proteobacteria, Firmicutes, Actinobacteria, Chlamydiae, Cyanobacteria, Deinococcota, Thermotogota, Aquificota, Bacteroidota, Spirochaetota and 4 archaeal phyla (Euryarchaeota x3, Crenarchaeota, Thaumarchaeota, Korarchaeota). Its purpose is stress-testing HGT-rich, mostly single-copy genomes (e.g. THEMA - bacterial/archaeal exchange; METAC, HALSA - bacterial-origin genes in archaea). Cheap enough to run twice (e.g. with and without LOKOS-type outgroups) and a natural place for the **mitochondrial/plastid-derived eukaryotic gene families** if the literature set includes them (cross-panel E+K).

### 3.6 Lite variants, optional tiers and compute budget

| Variant | Recipe | Proteins |
|---|---|---|
| V-lite | V core minus CIOIN, CHICK, LATCH, MOUSE | 256,842 |
| P-lite | P core minus SOYBN, ORYSJ, KLENI | 273,684 |
| E-lite | E core minus NEMVE, TRIAD, NAEGR, ENCCU | 286,588 |
| Y-ext | Y core + optional 4 | 106,426 |
| V-ext, P-ext, E-ext | core + optional | 504,031 / 576,787 / 571,852 |

All panels with their per-species counts are in `literature/raw/C/panels_verified.tsv` (generated by `scripts/C_verify_panels.py` from `panel_candidates.tsv`; one row per panel/species including REST type, assembly, REST total entries, RP table counts, QfO flag/count, OMA and Ensembl matches).

---

## 4. Reference species trees and controversial nodes

**Generic route (all panels).** (1) NCBI taxonomy topology via `ete3.NCBITaxa.get_topology` (what the FastOMA paper used by default, "the NCBI tree was used by downloading via the ETE3 package"); I generated these for all candidate species: `literature/raw/C/trees/panel_{V,Y,P,E,K}_ncbi_topology.nwk` (script `scripts/C_ncbi_species_trees.py`; leaves named by UniProt code). They contain **polytomies** (counts across all candidates: V 1, Y 2, P 3, E 5, K 4; e.g. Y: a 13-way node over the Saccharomycetaceae; V: HUMAN/MOUSE/ORNAN; E: 10-way over excavates + others). A taxonomy tree is not a defensible species tree for the WGD nodes and yields multifurcations that some methods may not accept (OrthoFinder species-tree input requirements for non-binary trees: **UNVERIFIED**). (2) Preferred: prune a published phylogenomic tree to the panel (section 4.1-4.5) and resolve remaining polytomies explicitly, recording the source tree and the pruning script. (3) TimeTree 5 (Kumar 2022) gives dated species-level trees from a species list (the FastOMA paper used the TimeTree web server for this); useful for branch lengths. (4) OMA exposes a taxonomy-based subtree through `/api/taxonomy/`. (5) OrthoFinder can infer its own species tree (STAG/STRIDE, Emms & Kelly 2018 bioRxiv) - useful as an **independent check**, not as the reference.

### 4.1 Vertebrates (panel V)

* Backbone: Irisarri 2017 (7,189 nuclear genes, 100 jawed-vertebrate taxa); Hughes 2018 and Betancur-R 2017 (ray-finned fishes; Hughes settles the branching order at the base of teleosts); lungfish are the closest living relatives of tetrapods (Meyer 2021) so coelacanth is placed below the lungfish+tetrapod split; tunicates (not cephalochordates) are the closest relatives of vertebrates (Delsuc 2006), so BRAFL is the outermost chordate in V.
* Contested/complex nodes: (a) **number and timing of 1R/2R**: Simakov 2020 - first duplication shared by all vertebrates (autotetraploid, mid/late Cambrian), second only in jawed vertebrates (allotetraploid, mid-late Ordovician); Nakatani 2021 - cyclostomes diverged after 1R but before gnathostome 2R, plus an additional **cyclostome hexaploidization**, 2R allotetraploid with biased loss; Sacerdot 2018 and Dehal & Boore 2005 are earlier treatments. Consequence: lamprey families have a different duplication history than gnathostome ones; at the "Vertebrata" level a HOG may legitimately contain 2-8 gnathostome copies and extra cyclostome copies. (b) **teleost 3R and the pre-3R outgroup**: gar (Braasch 2016) is the benchmarking outgroup; Meyer & Van de Peer 2005 for the "fish-specific genome duplication" concept. (c) **Salmonid Ss4R**: rediploidization lag means duplicate pairs have heterogeneous divergence times (Lien 2016; Macqueen & Johnston 2014). (d) **Xenopus**: L and S subgenomes (Session 2016).
* DOIs: all in section 7.

### 4.2 Yeasts (panel Y)

* Backbone: Shen 2018 (332 species, 12 major clades; trees and data on Figshare doi:10.6084/m9.figshare.5854692 - DOI taken from the paper's Key Resources table; the Figshare landing page itself was not fetched, **UNVERIFIED**) and Shen 2016 (backbone from 25 genomes).
* WGD: Wolfe & Shields 1997; Kellis 2004; Dujon 2004; Gordon 2009 (ancestral genome reconstruction). **Marcet-Houben & Gabaldon 2015**: the WGD was preceded by interspecies hybridization; the post-WGD subgenomes have phylogenetic affinities to the ZT clade (*Zygosaccharomyces*, *Torulaspora*) and the KLE clade (*Kluyveromyces*, *Lachancea*, *Eremothecium*). Wolfe 2015 (PLoS Biol primer): earlier concatenation/supertree analyses put ZT sister to post-WGD with KLE as outgroup; the gene-by-gene analysis finds that ~30% of post-WGD genes instead cluster with KLE or outgroup KLE+ZT. The paper itself states that "the difficulty of correctly determining the position of parent B" arises from gene conversion and long-branch effects and that parent B diverged just before KLE in 33% of trees.
* **Practical consequence for the benchmark:** the position of *Lachancea* inside KLE and where the post-WGD clade attaches are exactly the nodes where a single species tree misrepresents gene histories. WGD-derived ohnolog pairs can legitimately be non-monophyletic with respect to the species tree. Treat ZT/KLE/WGD nodes as **soft** in scoring (report both with and without) and use ZYGRC/TORDC (ZT) and LACTH/KLULA/ASHGO/ERECY (KLE) so that both parental signals are present. YGOB (Byrne & Wolfe 2005, NAR 2006 doi 10.1093/nar/gkj041) and POInT-browse (2023) are the curated synteny-based orthology resources for these genomes.
* Gotcha: *Lachancea waltii* (the classic non-WGD reference in YGOB) is not in UniProt RP; *L. thermotolerans* substitutes.

### 4.3 Plants (panel P)

* Backbone: One Thousand Plant Transcriptomes Initiative 2019 (1,124 transcriptomes; discordance between plastid and nuclear trees at several nodes, polyploidy, rapid speciation); angiosperm-wide phylogenomics Zuntini 2024 ("Phylogenomics and the rise of the angiosperms"); Puttick 2018 resolves Zygnematophyceae as sister to land plants (Wickett 2014 is the earlier phylotranscriptomic study of streptophyte origins; its specific result was not re-checked here).
* Contested nodes: **bryophytes monophyletic vs paraphyletic** - Puttick 2018 (mosses+liverworts "Setaphyta" sister to the rest, or hornworts; bryophyte monophyly supported by gene-tree-based analyses), Harris 2020 (bryophyte monophyly); this determines how MARPO/PHYPA relate to SELMO and angiosperms. **Amborella as the single living species of the sister lineage to all other angiosperms** (Amborella 2013, abstract; I did not confirm how Zuntini 2024 treats this node, **UNVERIFIED**). **Seed-plant and angiosperm ancestral WGDs** (Jiao 2011: gene duplications concentrated around 319 and 192 [Ma per the abstract]) place duplications near the roots of seed-plant and angiosperm HOGs. Lineage WGDs: Physcomitrium (Rensing 2008; WGD events inferred, not checked), Brassica triplication (Wang 2011), soybean palaeopolyploidy (Schmutz 2010), maize (Schnable 2009), poplar (Tuskan 2006), grape hexaploidization (Jaillon 2007); review: Van de Peer 2017.
* Note: the OneKP paper itself flags widespread incomplete lineage sorting and polyploidy at some nodes; benchmark scoring should tolerate polytomies at the Zygnematophyceae/land-plant, bryophyte and early-angiosperm nodes.

### 4.4 Eukaryota-wide (panel E)

* Backbone: Burki 2020 ("The New Tree of Eukaryotes": supergroup labels are largely arbitrary; many deep nodes unresolved); eukaryotic root: Derelle 2015 (root between Unikonta/"Opimoda" and Bikonta with malawimonads and collodictyonids on the unikont side; Excavata not monophyletic). Archaeal ancestry: Williams 2020 (eukaryotes nested within archaea, Asgard/Heimdallarchaeota closest), Eme 2023 (eukaryotes sister to Hodarchaeales within Heimdallarchaeia), Zaremba-Niedzwiedzka 2017 (Asgard).
* **Animal root (ctenophore vs sponge sister)**: ctenophores-first: Whelan 2015 (placement of Ctenophora sister to all other animals, with the caveat of systematic error), Schultz 2023 (chromosome-scale synteny); sponges-first (Porifera-sister): Simion 2017, Redmond & McLysaght 2021; early reference points: Dunn 2008 and Philippe 2009 (sampling of ctenophores/sponges; their individual topologies were not re-checked); Kapli & Telford 2020 argue ctenophore-first signals can come from systematic error. The panel contains one of each, so a "root-agnostic" scoring (treat the Metazoa-internal root as a trifurcation of Ctenophora, Porifera, rest) is advisable. Xenacoelomorpha placement is also contested (Kapli & Telford) but no xenacoelomorph is in the panel.
* Practical: the NCBI tree is polytomous for excavates (LEIMA/TRYB2/NAEGR/GIAIC...), the SAR/Haptista/Cryptista region and archaea. For the E panel, take the supergroup-level topology from Burki 2020/Derelle 2015 and animal/fungal/plant internal topologies from the panel-specific trees above.

### 4.5 Prokaryotes (panel K)

* Backbone: GTDB (Parks 2018 Nat Biotechnol; Parks 2022 NAR) as taxonomic/phylogenomic reference; bacterial root: Coleman 2021 (between Terrabacteria and Gracilicutes; HGT affects 92% of gene families although two-thirds of transmissions are vertical); archaeal/eukaryote tree: Williams 2020, Eme 2023; LUCA: Moody 2024. Contested nodes: root of bacteria; placement of the Candidate Phyla Radiation (Coleman 2021 places CPR within Terrabacteria; CPR is not in the panel); Thermotogota/Aquificota position; TACK/DPANN membership. **HGT is pervasive**, so a bifurcating species tree is only an approximation; scoring should avoid penalizing apparent HGT-induced topology violations or should list the HGT families separately.

---

## 5. Compute estimates (OrthoFinder v3 + FastOMA on 12 cores / 500 GB)

**Published (verified in the papers' text):**

* *FastOMA* (Majidian 2025, Nat Methods): "all 2,086 eukaryotic UniProt reference proteomes in under 24 h, using 300 CPU cores" (<= ~7,200 CPU-h, i.e. at most ~3.5 CPU-h per proteome as an upper bound for that huge, very diverse set); most steps < 10 GB RAM, max 32 GB on the QfO set; the original OMA algorithm managed only 50 genomes in the same time; OrthoFinder and SonicParanoid shown as quadratic. FastOMA uses OMAmer placement to a database (LUCA.h5) and then per-HOG alignment/tree inference with FastTree (default subsample of 20 sequences per HOG at deep levels).
* *OrthoFinder v3* (Emms et al., Nat Methods 2026; bioRxiv 2025): all tools run with 32 threads (16 cores) and up to 200 GB RAM on Ensembl rapid-release proteomes; **<= 64 species use the classic all-vs-all workflow**, beyond that a "core 64 + assign" mode; OrthoFinder v3 ran 1,024 eukaryotic species in 128 h (~8x faster than OrthoFinder v2 above 64 species, ~4-fold lower peak memory than competing methods that completed); on bacterial proteome sets downloaded from Ensembl it completed 2,048 species in 50 h and 4,096 species in 13 d 15 h (504 GB RAM), while FastOMA needed 14 days for the 2,048-species bacterial set and could not finish 4,096. On a 64-proteome set both OrthoFinder v3 and FastOMA gained little beyond 32 threads (runtime plateau), so 12 cores will not be thread-starved by small panels. The text gives no wall-clock numbers for small panels (they are in a figure, Fig. 4a, not checked).
* OrthoFinder 2019 (Emms & Kelly): states that the default DIAMOND workflow gives "equivalent speed and scalability to the fastest score-based heuristic methods"; no runtimes extracted.

**Measured on this machine (by the team, `runs/timing/`; I did not repeat it because the node was already at load average ~14.6 on 12 cores):**

| Run | Input | Wall | CPU (user+sys) | Memory |
|---|---|---|---|---|
| OrthoFinder v3, 12 threads (`of.log`, `of_time.txt`) | 4 yeast proteomes YEAST/CANGA/KLULA/ASHGO = 21,098 proteins (6,067+5,200+5,071+4,760) | 167 s | 729 + 74 = 803 s | not recorded |
| OMAmer 2.1.0 vs LUCA.h5 (`om_time.txt`) | CANGA, 5,200 proteins | 41.5 s | 50 + 4 s | **14.3 GB max RSS** (database load) |
| FastOMA (Nextflow) (`fo.log`) | same 4 yeasts | failed (9 of 19 tasks failed, in `.nextflow.log`: missing `hogmap_in`/`splice` inputs) | - | - |

**Crude scaling model.** All-vs-all DIAMOND cost scales as N^2; calibrating on the single point above gives CPU-s ~ 803 x (N / 21,098)^2. This is an *upper-bound-flavoured* extrapolation: it assumes the 4-yeast run (closely related proteomes, many hits, machine under load, fixed overheads included) is representative, and it ignores superlinear gene-tree costs in giant orthogroups from polyploids. Wall-clock assumes 70% parallel efficiency on 12 cores (8.4 effective cores).

| Panel (core) | N proteins | OrthoFinder v3 CPU-h (model) | Wall on 12 cores (model) | FastOMA upper bound (<= 3.5 CPU-h per proteome) |
|---|---|---|---|---|
| K (21 spp) | 68,128 | 2.3 | ~0.3 h | (bacteria are FastOMA's slow case, see above) |
| Y (16) | 85,181 | 3.6 | ~0.4 h | <= 55 CPU-h (4.6 h wall at 12 cores) |
| V-lite (10) | 256,842 | 33 | ~3.9 h | <= 35 CPU-h |
| V (14) | 333,315 | 56 | ~6.6 h | <= 48 CPU-h (~4 h wall) |
| E (34) | 340,209 | 58 | ~6.9 h | <= 117 CPU-h (~10 h wall; likely much less since 9 genomes are prokaryotic and small) |
| P (13) | 389,385 | 76 | ~9.0 h | <= 45 CPU-h (~3.7 h wall) |
| V + optional (21) | 504,031 | 127 | ~15 h | - |

Interpretation and recommendations: (i) Y and K are trivially cheap and can serve as development/regression panels. (ii) V, E and P are at the edge of the "few hours" budget for OrthoFinder under this pessimistic model; the true cost for diverse panels (E, P) should be below the model because the model is calibrated on very close species, whereas for V (many near-identical vertebrate pairs) it may be close. **Do a quick two-point calibration** (e.g. 5 and 10 vertebrate proteomes, 12 threads, idle node) to fit the exponent before launching full runs; also time OrthoFinder's `-op` (DIAMOND-only) and the tree/MSA stage separately. (iii) FastOMA's measured cost per proteome here is unknown because the run failed; its OMAmer step costs ~40 s per 5k proteins plus ~14 GB RAM per process (so limit concurrent OMAmer jobs to ~20 on 500 GB). (iv) OrthoFinder with > 64 species would switch to the assign mode; all panels stay <= 64 species so this does not apply. (v) The project requirement that each tool finish within a few hours per panel is realistic for Y, K, V-lite and probably P-lite/E-lite; confirm with measurement.

---

## 6. Caveats, anomalies and UNVERIFIED items

### 6.1 Verified anomalies
* **UniProt RP README table, human row**: `#(1)` = 147,503 / `#(2)` = 22,131, whereas the human canonical FASTA (`UP000005640_9606.fasta.gz`, 7.4 MB) has **20,652** sequences and `gene2acc` has 172,773 rows (streamed counts; QfO 2026_02 also says 20,652). Use FASTA counts, not the README table, for validation. The other 51 non-QfO panel proteomes matched the table exactly.
* **REST `proteinCount` is not the canonical count** (human 147,520; mouse 54,855 vs canonical 21,860; zebrafish 69,241 vs 26,567 + 42,555 isoforms).
* Canonical counts in the RP table vs QfO differ slightly because the sources differ by one release (e.g. CAEEL 19,792 RP 2026_03 vs 19,789 QfO 2026_02; DANRE 26,567 vs 26,677; XENTR 21,732 vs 21,609). **Pick one source per species and pin it.**

### 6.2 Redundancy, hybrid and reduced genomes
* *X. laevis* (36,351/36,352) is in the canonical set as L and S homeologs; these are the "same gene" orthologs of one *X. tropicalis* gene, expected 2:1.
* *Saccharomyces pastorianus*/hybrid yeasts are in RP (e.g. UP000501346 with 10,697) - **not** used.
* *Trichomonas vaginalis* (50,190), *Paramecium* (39,350), *Emiliania* (35,603) have gene-rich/repeat-rich annotation and inflate N^2 compute; kept optional.
* Reduced/fast-evolving genomes (ENCCU 2,041, ASHGO 4,758, GIAIC 4,900, PLAF7 5,361, CYAME 5,018) are included on purpose as hard cases for orthology.

### 6.3 Versioning
* UniProt RP changes every 8 weeks; QfO yearly. UniProt keeps only whole-knowledgebase tarballs in `previous_releases/` (verified: `release-2026_01/knowledgebase/` has `knowledgebase2026_01.tar.gz` 151 GB; its `reference_proteomes/` path returns 404 in the listing I probed). The QfO tarballs are the only archived RP-like sets. => Download and checksum now.
* Assembly drift inside QfO: zebrafish moved from GCF_000002035.6 (QfO 2025_04) to **GCF_049306965.1** (2026_02) with isoform count doubling (20,450 -> 42,446) (STATS file). Old literature gene IDs (Ensembl GRCz11 / NCBI GRCz11) may need re-mapping.
* A species can have several RP records (e.g. *Brassica campestris* UP000011750, UP000694005, UP000823674; *Petromyzon marinus* has one Reference proteome and two "Excluded" records, UP000245300, UP001318040; *Salmo salar* old UP000087266 is "Excluded"). Always use the UPID in the panel table.

### 6.4 UNVERIFIED (could not check, or only secondary evidence)
1. OMA licence (CC BY-SA 2.5) - from a search snippet of the OMA licence page; page blocked by Cloudflare from the cluster.
2. YGOB: download availability, version-7 species set, licence; host unreachable.
3. Phytozome bulk-download mechanism; embargo status of any specific genome.
4. Whether the lamprey RP assembly is germline or somatic; whether UniProt/QfO plant canonical counts include low-confidence/TE models.
5. OrthoFinder handling of polytomic user species trees; FastOMA tree requirements beyond "NCBI tree via ete3 or TimeTree tree" (from the paper's Methods).
6. Zuntini 2024 treatment of Amborella; Physcomitrella WGD details; Figshare landing page for Shen 2018.
7. Per-run FastOMA runtimes on this machine (run failed); OrthoFinder time split DIAMOND vs trees.
8. The `UPKB_PROTEIN_GENE_MAPPING` semantics (only the first 600 bytes were read).
9. BioMart canonical attribute names for Ensembl (not tested; GFF3 tag was verified instead).

---

## 7. Verified citations (Europe PMC lookup of each DOI, 2026-10-26)

All rows below returned a record whose title matches the intended paper. Several DOIs I first guessed from memory were wrong and were **corrected before inclusion** (Meyer & Van de Peer 2005 = 10.1002/bies.20293; Hughes 2018 = 10.1073/pnas.1719358115; Byrne & Wolfe 2006 = 10.1093/nar/gkj041; Harris 2020 = 10.1016/j.cub.2020.03.048; Philippe 2009 = 10.1016/j.cub.2009.02.052; Macqueen & Johnston 2014 = 10.1098/rspb.2013.2881; Wolfe & Shields 1997 = 10.1038/42711).

| Key | DOI | Status | Title | Journal | Year |
|---|---|---|---|---|---|
| Simakov2020 | 10.1038/s41559-020-1156-z | VERIFIED | Deeply conserved synteny resolves early events in vertebrate evolution. | Nat Ecol Evol | 2020 |
| Nakatani2021 | 10.1038/s41467-021-24573-z | VERIFIED | Reconstruction of proto-vertebrate, proto-cyclostome and proto-gnathostome genomes provides new insights into  | Nat Commun | 2021 |
| Sacerdot2018 | 10.1186/s13059-018-1559-1 | VERIFIED | Chromosome evolution at the origin of the ancestral vertebrate genome. | Genome Biol | 2018 |
| Braasch2016 | 10.1038/ng.3526 | VERIFIED | The spotted gar genome illuminates vertebrate evolution and facilitates human-teleost comparisons. | Nat Genet | 2016 |
| Lien2016 | 10.1038/nature17164 | VERIFIED | The Atlantic salmon genome provides insights into rediploidization. | Nature | 2016 |
| Session2016 | 10.1038/nature19840 | VERIFIED | Genome evolution in the allotetraploid frog Xenopus laevis. | Nature | 2016 |
| Smith2013lamprey | 10.1038/ng.2568 | VERIFIED | Sequencing of the sea lamprey (Petromyzon marinus) genome provides insights into vertebrate evolution. | Nat Genet | 2013 |
| Venkatesh2014 | 10.1038/nature12826 | VERIFIED | Elephant shark genome provides unique insights into gnathostome evolution. | Nature | 2014 |
| Meyer2005 | 10.1002/bies.20293 | VERIFIED | From 2R to 3R: evidence for a fish-specific genome duplication (FSGD). | Bioessays | 2005 |
| Betancur2017 | 10.1186/s12862-017-0958-3 | VERIFIED | Phylogenetic classification of bony fishes. | BMC Evol Biol | 2017 |
| Hughes2018 | 10.1073/pnas.1719358115 | VERIFIED | Comprehensive phylogeny of ray-finned fishes (Actinopterygii) based on transcriptomic and genomic data. | Proc Natl Acad Sci U S A | 2018 |
| Irisarri2017 | 10.1038/s41559-017-0240-5 | VERIFIED | Phylotranscriptomic consolidation of the jawed vertebrate timetree. | Nat Ecol Evol | 2017 |
| Kumar2022TimeTree5 | 10.1093/molbev/msac174 | VERIFIED | TimeTree 5: An Expanded Resource for Species Divergence Times. | Mol Biol Evol | 2022 |
| Meyer2021lungfish | 10.1038/s41586-021-03198-8 | VERIFIED | Giant lungfish genome elucidates the conquest of land by vertebrates. | Nature | 2021 |
| Delsuc2006 | 10.1038/nature04336 | VERIFIED | Tunicates and not cephalochordates are the closest living relatives of vertebrates. | Nature | 2006 |
| Shen2018 | 10.1016/j.cell.2018.10.023 | VERIFIED | Tempo and Mode of Genome Evolution in the Budding Yeast Subphylum. | Cell | 2018 |
| Shen2016 | 10.1534/g3.116.034744 | VERIFIED | Reconstructing the Backbone of the Saccharomycotina Yeast Phylogeny Using Genome-Scale Data. | G3 (Bethesda) | 2016 |
| WolfeShields1997 | 10.1038/42711 | VERIFIED | Molecular evidence for an ancient duplication of the entire yeast genome. | Nature | 1997 |
| Kellis2004 | 10.1038/nature02424 | VERIFIED | Proof and evolutionary analysis of ancient genome duplication in the yeast Saccharomyces cerevisiae. | Nature | 2004 |
| Dujon2004 | 10.1038/nature02579 | VERIFIED | Genome evolution in yeasts. | Nature | 2004 |
| MarcetHouben2015 | 10.1371/journal.pbio.1002220 | VERIFIED | Beyond the Whole-Genome Duplication: Phylogenetic Evidence for an Ancient Interspecies Hybridization in the Ba | PLoS Biol | 2015 |
| Wolfe2015 | 10.1371/journal.pbio.1002221 | VERIFIED | Origin of the Yeast Whole-Genome Duplication. | PLoS Biol | 2015 |
| Gordon2009 | 10.1371/journal.pgen.1000485 | VERIFIED | Additions, losses, and rearrangements on the evolutionary route from a reconstructed ancestor to the modern Sa | PLoS Genet | 2009 |
| ByrneWolfe2005 | 10.1101/gr.3672305 | VERIFIED | The Yeast Gene Order Browser: combining curated homology and syntenic context reveals gene fate in polyploid s | Genome Res | 2005 |
| ByrneWolfe2006 | 10.1093/nar/gkj041 | VERIFIED | Visualizing syntenic relationships among the hemiascomycetes with the Yeast Gene Order Browser. | Nucleic Acids Res | 2006 |
| Souciet2009 | 10.1101/gr.091546.109 | VERIFIED | Comparative genomics of protoploid Saccharomycetaceae. | Genome Res | 2009 |
| POInTbrowse2023 | 10.1186/s12859-023-05298-w | VERIFIED | POInT<sub>browse</sub>: orthology prediction and synteny exploration for paleopolyploid genomes. | BMC Bioinformatics | 2023 |
| OneKP2019 | 10.1038/s41586-019-1693-2 | VERIFIED | One thousand plant transcriptomes and the phylogenomics of green plants. | Nature | 2019 |
| Amborella2013 | 10.1126/science.1241089 | VERIFIED | The Amborella genome and the evolution of flowering plants. | Science | 2013 |
| Zuntini2024 | 10.1038/s41586-024-07324-0 | VERIFIED | Phylogenomics and the rise of the angiosperms. | Nature | 2024 |
| Puttick2018 | 10.1016/j.cub.2018.01.063 | VERIFIED | The Interrelationships of Land Plants and the Nature of the Ancestral Embryophyte. | Curr Biol | 2018 |
| Harris2020 | 10.1016/j.cub.2020.03.048 | VERIFIED | Phylogenomic Evidence for the Monophyly of Bryophytes and the Reductive Evolution of Stomata. | Curr Biol | 2020 |
| Wickett2014 | 10.1073/pnas.1323926111 | VERIFIED | Phylotranscriptomic analysis of the origin and early diversification of land plants. | Proc Natl Acad Sci U S A | 2014 |
| Bowman2017 | 10.1016/j.cell.2017.09.030 | VERIFIED | Insights into Land Plant Evolution Garnered from the Marchantia polymorpha Genome. | Cell | 2017 |
| Rensing2008 | 10.1126/science.1150646 | VERIFIED | The Physcomitrella genome reveals evolutionary insights into the conquest of land by plants. | Science | 2008 |
| Banks2011 | 10.1126/science.1203810 | VERIFIED | The Selaginella genome identifies genetic changes associated with the evolution of vascular plants. | Science | 2011 |
| Nishiyama2018 | 10.1016/j.cell.2018.06.033 | VERIFIED | The Chara Genome: Secondary Complexity and Implications for Plant Terrestrialization. | Cell | 2018 |
| Merchant2007 | 10.1126/science.1143609 | VERIFIED | The Chlamydomonas genome reveals the evolution of key animal and plant functions. | Science | 2007 |
| Jiao2011 | 10.1038/nature09916 | VERIFIED | Ancestral polyploidy in seed plants and angiosperms. | Nature | 2011 |
| VandePeer2017 | 10.1038/nrg.2017.26 | VERIFIED | The evolutionary significance of polyploidy. | Nat Rev Genet | 2017 |
| Wang2011Brassica | 10.1038/ng.919 | VERIFIED | The genome of the mesopolyploid crop species Brassica rapa. | Nat Genet | 2011 |
| Schmutz2010 | 10.1038/nature08670 | VERIFIED | Genome sequence of the palaeopolyploid soybean. | Nature | 2010 |
| Schnable2009 | 10.1126/science.1178534 | VERIFIED | The B73 maize genome: complexity, diversity, and dynamics. | Science | 2009 |
| Tuskan2006 | 10.1126/science.1128691 | VERIFIED | The genome of black cottonwood, Populus trichocarpa (Torr. & Gray). | Science | 2006 |
| Jaillon2007 | 10.1038/nature06148 | VERIFIED | The grapevine genome sequence suggests ancestral hexaploidization in major angiosperm phyla. | Nature | 2007 |
| Burki2020 | 10.1016/j.tree.2019.08.008 | VERIFIED | The New Tree of Eukaryotes. | Trends Ecol Evol | 2020 |
| Hug2016 | 10.1038/nmicrobiol.2016.48 | VERIFIED | A new view of the tree of life. | Nat Microbiol | 2016 |
| Williams2020 | 10.1038/s41559-019-1040-x | VERIFIED | Phylogenomics provides robust support for a two-domains tree of life. | Nat Ecol Evol | 2020 |
| Eme2023 | 10.1038/s41586-023-06186-2 | VERIFIED | Inference and reconstruction of the heimdallarchaeial ancestry of eukaryotes. | Nature | 2023 |
| Zaremba2017 | 10.1038/nature21031 | VERIFIED | Asgard archaea illuminate the origin of eukaryotic cellular complexity. | Nature | 2017 |
| Dunn2008 | 10.1038/nature06614 | VERIFIED | Broad phylogenomic sampling improves resolution of the animal tree of life. | Nature | 2008 |
| Whelan2015 | 10.1073/pnas.1503453112 | VERIFIED | Error, signal, and the placement of Ctenophora sister to all other animals. | Proc Natl Acad Sci U S A | 2015 |
| Simion2017 | 10.1016/j.cub.2017.02.031 | VERIFIED | A Large and Consistent Phylogenomic Dataset Supports Sponges as the Sister Group to All Other Animals. | Curr Biol | 2017 |
| Schultz2023 | 10.1038/s41586-023-05936-6 | VERIFIED | Ancient gene linkages support ctenophores as sister to other animals. | Nature | 2023 |
| Redmond2021 | 10.1038/s41467-021-22074-7 | VERIFIED | Evidence for sponges as sister to all other animals from partitioned phylogenomics with mixture models and rec | Nat Commun | 2021 |
| Kapli2020 | 10.1126/sciadv.abc5162 | VERIFIED | Topology-dependent asymmetry in systematic errors affects phylogenetic placement of Ctenophora and Xenacoelomo | Sci Adv | 2020 |
| Derelle2015 | 10.1073/pnas.1420657112 | VERIFIED | Bacterial proteins pinpoint a single eukaryotic root. | Proc Natl Acad Sci U S A | 2015 |
| Parks2018 | 10.1038/nbt.4229 | VERIFIED | A standardized bacterial taxonomy based on genome phylogeny substantially revises the tree of life. | Nat Biotechnol | 2018 |
| Parks2022GTDB | 10.1093/nar/gkab776 | VERIFIED | GTDB: an ongoing census of bacterial and archaeal diversity through a phylogenetically consistent, rank normal | Nucleic Acids Res | 2022 |
| Coleman2021 | 10.1126/science.abe0511 | VERIFIED | A rooted phylogeny resolves early bacterial evolution. | Science | 2021 |
| Moody2024 | 10.1038/s41559-024-02461-1 | VERIFIED | The nature of the last universal common ancestor and its impact on the early Earth system. | Nat Ecol Evol | 2024 |
| Altenhoff2016QfO | 10.1038/nmeth.3830 | VERIFIED | Standardized benchmarking in the quest for orthologs. | Nat Methods | 2016 |
| Nevers2022 | 10.1093/nar/gkac330 | VERIFIED | The Quest for Orthologs orthology benchmark service in 2022. | Nucleic Acids Res | 2022 |
| Linard2021 | 10.1093/molbev/msab098 | VERIFIED | Ten Years of Collaborative Progress in the Quest for Orthologs. | Mol Biol Evol | 2021 |
| OMA2024 | 10.1093/nar/gkad1020 | VERIFIED | OMA orthology in 2024: improved prokaryote coverage, ancestral and extant GO enrichment, a revamped synteny vi | Nucleic Acids Res | 2024 |
| FastOMA2025 | 10.1038/s41592-024-02552-8 | VERIFIED | Orthology inference at scale with FastOMA. | Nat Methods | 2025 |
| OrthoFinder2019 | 10.1186/s13059-019-1832-y | VERIFIED | OrthoFinder: phylogenetic orthology inference for comparative genomics. | Genome Biol | 2019 |
| OrthoFinder2015 | 10.1186/s13059-015-0721-2 | VERIFIED | OrthoFinder: solving fundamental biases in whole genome comparisons dramatically improves orthogroup inference | Genome Biol | 2015 |
| STAG2018 | 10.1101/267914 | VERIFIED | STAG: Species Tree Inference from All Genes | bioRxiv | 2018 |
| Macqueen2014 | 10.1098/rspb.2013.2881 | VERIFIED | A well-constrained estimate for the timing of the salmonid whole genome duplication reveals major decoupling f | Proc Biol Sci | 2014 |
| Betancur2013 | 10.1371/currents.tol.53ba26640df0ccaee75bb165c8c26288 | VERIFIED | The tree of life and a new classification of bony fishes. | PLoS Curr | 2013 |
| Dehal2005 | 10.1371/journal.pbio.0030314 | VERIFIED | Two rounds of whole genome duplication in the ancestral vertebrate. | PLoS Biol | 2005 |
| Putnam2008 | 10.1038/nature06967 | VERIFIED | The amphioxus genome and the evolution of the chordate karyotype. | Nature | 2008 |
| Philippe2009 | 10.1016/j.cub.2009.02.052 | VERIFIED | Phylogenomics revives traditional views on deep animal relationships. | Curr Biol | 2009 |
| OrthoFinder3_NatMethods2026 | 10.1038/s41592-026-03126-6 | VERIFIED | OrthoFinder: improved phylogenetic orthology inference with enhanced accuracy and scalability. | Nat Methods | 2026 |
| OrthoFinder3_bioRxiv2025 | 10.1101/2025.07.15.664860 | VERIFIED | OrthoFinder: scalable phylogenetic orthology inference for comparative genomics | bioRxiv | 2025 |
| SartonLoheac2025HOG | 10.1007/s00239-025-10277-1 | VERIFIED | Reconstructing Evolutionary Histories with Hierarchical Orthologous Groups. | J Mol Evol | 2025 |
| Majidian2025QfOreview | 10.1007/s00239-025-10272-6 | VERIFIED | Quest for Orthologs in the era of Data Deluge and AI: Challenges and Innovations in Orthology Prediction and D | J Mol Evol | 2025 |
| Smith2018germline | 10.1038/s41588-017-0036-1 | VERIFIED | The sea lamprey germline genome provides insights into programmed genome rearrangement and vertebrate evolutio | Nat Genet | 2018 |

Additional bibliographic notes: QfO 2022 service paper = Nevers 2022 (row `Nevers2022`); QfO benchmark = Altenhoff 2016; HOG review = Sarton-Loheac 2025; OrthoFinder v3 = Emms 2026 (Nat Methods; bioRxiv 2025 version listed). A hierarchical-orthology paper search in Europe PMC (title words "Quest for Orthologs" + "hierarchical") returned no hit; none of the above is claimed to be a hierarchical-orthology benchmark.

---

## 8. Reproducibility: files created by part C

* Text sources / metadata in `literature/raw/C/`: `QfO_2026_02_{README,STATS,PROTEOME_ASSEMBLY_MAPPING,species.tsv}`, `UniProt_RP_2026_03_README_with_table.txt` + `UniProt_RP_2026_03_species.tsv` (37,190 rows), `oma_genomes_All.May2026.json`, `ensembl_species_Ensembl{Vertebrates,Plants,Fungi,Protists,Metazoa}.json`, `panel_candidates.tsv` (input: panel, code, taxid, UPID, tier, rationale), `panels_verified.tsv` (output of verification), `fasta_canonical_counts.tsv` (streamed counts, nothing stored), `assembly_meta.tsv` (NCBI Datasets, partial: only 44 of ~120 accessions returned), `dois_to_verify.tsv` + `dois_verified.tsv`, `trees/panel_*_ncbi_topology.nwk` (the ete3 `taxa.sqlite`, ~740 MB, was deleted after use; `NCBITaxa(dbfile=...)` rebuilds it from the NCBI taxdump in a few minutes), `md5sums_metadata.txt` (checksums of the QfO README/STATS, the UniProt RP README+table and the OMA genome list as downloaded), `qfo_table.md`, `panel_tables.md`, `doi_table.md`, XML full texts of PMC articles used for statements (FastOMA, OrthoFinder 2019/3, Shen 2018 HTML, Marcet-Houben, Wolfe 2015).
* Scripts (`scripts/`): `C_uniprot_proteomes.py` (proteome records by taxid), `C_rp_lookup.sh` (grep RP table), `C_verify_panels.py`, `C_panel_tables.py`, `C_count_fasta.sh` (stream-count), `C_uniprot_idmap.py`, `C_verify_dois.py`, `C_epmc_title.py`, `C_epmc_abs.py`, `C_ncbi_assembly_meta.py`, `C_ncbi_species_trees.py`, `C_build_md.py` (+ `C_text_template.md`) which regenerate this file.
* To refresh: re-run `C_verify_panels.py` -> `C_panel_tables.py` -> `C_build_md.py` (network required; the extrapolation table in section 5 is computed by hand from the counts and the measured 803 CPU-s at N = 21,098).
