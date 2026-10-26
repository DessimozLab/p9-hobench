# Round 1 review: data scientist (genomics data)

Checked against the files on disk and the UniProt REST API (2026-10-26).

1. **DS-1 (MUST): one pinned release, offline ID mapping.** DESIGN §3 says "RP 2026_03 / QfO 2026_02", but all 111 `data/proteomes/*/manifest.json`, QfO species included, are `current_release` 2026_03.
   - Rewrite §3 to "UniProt RP 2026_03 canonical FASTA for all species".
   - A new UniProt release could appear mid-project, and `map_ids.py` queries the live REST API. Download each proteome's `.idmapping.gz` and map against 2026_03 first, with REST as fallback only.
   - Make `fetch_proteome` fail if the release is not the pinned one.
   - Add a validator: every `members.tsv` gene_id must be in its species' pinned `proteome.fasta`. Store a sequence sha256 and the UniParc UPI for each member.

2. **DS-2 (MUST): one species table per panel, trees validated against it.**
   - `literature/raw/C/trees/panel_Y_ncbi_topology.nwk` has the leaf `4956` instead of ZYGRC (ete3 replaced the strain taxid with the species taxid).
   - Its internal nodes are unnamed, which `speciestree.py` rejects.
   - Codes differ from UniProt mnemonics (LACTH/LACTC, ASHGO/EREGS, NAUCC/NAUCA, SELMO/SELML, LOKOS/9ARCH...).

   Create `panels/<P>/species.tsv` with code, UPID, taxid, mnemonic, OMA code, FASTA sha256 and completeness. Validate that tree leaves = codes = tool FASTA basenames = members' species. Name each level by its NCBI name if the species sets match exactly, else by a cited literature name (PostWGD, ZT, KLE). List them in `levels.tsv`. Give the identical binary, named Newick to both tools and to scoring.

3. **DS-3 (MUST): map to genes, then canonical proteins; stricter thresholds.**
   - Route: ID → UniProt acc → `gene2acc` gene → canonical protein in the pinned proteome. `map_ids` keeps one-to-many hits (TP53 → P04637 plus two TrEMBL duplicates).
   - Retired Ensembl/RefSeq IDs: fetch the archived sequence, then map by sequence.
   - Gene_Name-only mappings need sequence confirmation.
   - §5.3's "≥95% identity" can pick the wrong copy for XENLA L/S, salmonid ohnologs and tandem clusters. Use ≥98% identity (≥99% for XENLA/SALSA/TAND families), ≥80% coverage of query and target, a reciprocal best hit, and a ≥2-point margin to the runner-up (code default 1.0).
   - Record merged/split-model collisions. Never map across species (L. waltii → L. thermotolerans).
   - Write a mapping report per family. Keep a family out of the headline if <80% of members map uniquely.

4. **DS-4 (MUST): define "declared complete"** (used in §10, never defined). A (family, species) cell is complete only if all three hold:
   - (i) the source states a genome-wide search of that species;
   - (ii) ≥90% of its members map uniquely;
   - (iii) a homology screen against the pinned proteome, run from the reference only, finds ≤10% unexplained candidates (proteins scoring at least the family's minimum within-reference bitscore).

   Write `families/<ID>/completeness.tsv`. Count precision penalties only in complete cells. Treat candidates as "unlabelled homologs": no pairs, but counted in contamination.

5. **DS-5 (MUST): fix the subset construction (§7).** A one-hop neighbourhood drops the second-hop proteins with promiscuous domains that cause merging. It under-tests over-merging and inflates precision exactly for the promiscuous-domain families (kinases, KRAB-ZNF, NLR, Rab/ARF).
   - It changes OrthoFinder's per-species-pair length/score normalisation and its MCL graph.
   - FastOMA (per-protein OMAmer placement) is nearly invariant, so the mode effect differs by tool.
   - Change: 2-hop closure with a cap and sizes reported per hop. Search the full panel with fixed `--dbsize` and `--more-sensitive` (LECA families are ~20–30% identity). Freeze the subset with `refs-frozen-v1`. Put a `mode` column on every result row and never pool modes.

6. **DS-6 (MUST): concrete validation.** "Agree within CI" passes almost anything because CIs are wide with few families. Instead:
   - Run per tool, comparing paired per-family Δ in merge-level accuracy and per-level F1 with a ±0.05 equivalence margin.
   - Add a reference-free check: ARI between the full- and subset-mode HOG partitions on reference members, per level.
   - Y and K lack large promiscuous families, so also test on V: 1-hop vs 2-hop and 1k vs 5k background. Results outside the margin make V/P/E results conditional on the mode.

7. **DS-7 (MUST): licensing, git hygiene, DATA.md.**
   - `taxdump.tar.gz` (80 MB) and `literature/raw/B1/{kinbase,ygob}` (scraped, no licence) are untracked and not ignored. Ignore them; commit only derived IDs/labels with provenance.
   - DATA.md does not exist. Give each dataset its URL, release, date, sha256, licence, "redistributed y/n" and regenerating script. Cover UniProt (CC BY 4.0), OMAmer `LUCA.h5` (OMA licence, possibly CC BY-SA, unverified), NCBI taxonomy, and each supplement (OA licence or "facts only").

8. **DS-8 (SHOULD): completeness on the canonical sets.**
   - UniProt BUSCO counts isoforms (pike D=59.3% without 4R, DROME 58.2%), so it is useless for QC. Run OMArk on FastOMA's OMAmer hogmaps; it also reports contamination and fragments.
   - Flag SELML (D=65%, possible haplotype redundancy) and low completeness: EPTBU 57.5%, CIOIN 78.7%, MNELE 80.9%, PETMA 81.1%, ORYSJ 84.3%.
   - GIAIC/ENTH1/TRIV3 look low because they are divergent, not necessarily incomplete.
   - Use these values as covariates in the difficulty analysis.

9. **DS-9 (SHOULD): background size.** 1,000 proteins is 49% of ENCCU but 2.4% of SALSA. Use a fixed fraction (10%, minimum 500, seed from `seeds.py`), with sampled IDs in the manifest.

10. **DS-10 (SHOULD): reconcile §3 with the disk.**
    - Pike (UP000265140) and trout (UP000694395, 119k entries) are not downloaded, so "~370k" is unverified.
    - ZYGRC came via the REST fallback (organism=null, no gene2acc).
    - Generate the panel tables from the manifests.

11. **DS-11 (COULD): D panel.** Good 2014 uses 2007 FlyBase GF IDs, so drift will be severe. Attempt it only if a 20-gene pilot reaches ≥80% unique mapping.
