# CURATION PROTOCOL v1 (binding for all curation agents)

Purpose: turn publications into reference families for HOBench without inventing, leaking or silently dropping anything.
Read first: DESIGN.md, design/core_spec.md, `hobench/reference.py` docstring (on-disk format), `curation_tools/README.md`.

## 0. Hard rules
1. **Blindness.** Never open `runs/`, `results/`, or any OrthoFinder/FastOMA/OMA/Ensembl-Compara/PANTHER/eggNOG/OrthoDB
   output, and never use them to decide membership or events.
2. **No priors as facts.** Every member label, D-node level, phenomenon tag, decoy and `distinct_from` claim must be backed
   by a source locator (file + table row / figure / section) **and** a verbatim quote (≤300 chars) that a script can find in
   the stored source text. If you "know" something but cannot quote it → encode as `U` / leave untagged.
3. **Denylist (not acceptable as orthology truth):** OMA, OrthoFinder, Ensembl Compara/HCOP and Compara-seeded names (VGNC),
   OHNOLOGS (Singh & Isambert; computational), PANTHER, TreeFam, eggNOG, OrthoDB, PLAZA, NCBI "LOC"/"uncharacterized" names,
   unreviewed electronically-projected names, KRAB-ZNF `age_MA` (allowed for integrity metrics only). Allowed: expert
   phylogenies, synteny-curated sets (YGOB, Xenbase L/S, Hox nomenclature with synteny basis), curated classification tables
   from papers, SwissTree. Record `source_method` for every source (`phylogeny`, `synteny`, `expert_classification`,
   `nomenclature`, `HMM_classification`).
4. **Members come from scripts, not typing.** Each family has `families/<ID>/build.py` that reads the stored source files
   (supplement tables, Newick, FASTA) and writes `members.tsv`, `reference_tree.nhx`, `family.yaml`. Hand-entered data is
   allowed only for T1 events (levels, tags) and must carry quotes.
5. **Mapping thresholds** (onto pinned UniProt 2026_03 panel proteomes, `panels/<P>/proteomes/<CODE>.fa`): exact ID mapping
   (UniProt ID mapping / gene name within species with a single hit) or same-species sequence mapping with ≥98% identity
   (≥99% for XENLA, salmonids, TAND-tagged families), ≥80% coverage of both sequences, reciprocal best hit, ≥2-point margin
   over the runner-up. **v1.1 (PI decision T+5h50):** when a source provides only a curated domain sequence (e.g. kinesin/myosin motor domain), `mapping_method` contains `domain`: require query length ≥250 aa, ≥90% query coverage at the identity threshold, reciprocal best hit and ≥2-point margin; target coverage is waived (columns `mapping_qcov`, `mapping_qlen`). Never map across species. Unmapped/ambiguous members are listed in `mapping_report.tsv`, never
   silently dropped. A family enters the headline set only with ≥80% of in-panel members uniquely mapped.
6. **Pseudogenes, identical copies, organelle-encoded genes** are excluded and listed in the mapping report.
7. Keep a `curation_log.md` per family: what was read, decisions, problems, time spent.

## 1. Directory layout per family `families/<ID>/`
```
family.yaml            # metadata (below)
members.tsv            # gene_id (UniProt acc in panel FASTA), species (panel code), original_id_in_paper, label,
                       # mapping_method, mapping_identity, mapping_coverage, source, locator
reference_tree.nhx     # S/D/U/H tree (below)
mapping_report.tsv     # every paper member: status mapped|unmapped|ambiguous|not_in_panel|excluded + reason
assertions.tsv         # one row per T1 claim: assertion_id, type (event|tag|decoy|distinct_from|completeness), target
                       # (node id / gene ids), tier, source_id, locator, quote
build.py               # regenerates the above from literature/papers/<id>/ files + panel FASTAs
curation_log.md
```
`family.yaml` keys: `id, name, panel, collection (optional), phenomena, descriptors, sources: [{id, doi, pmid, title,
year, source_method, files}], root_level, root_single_copy (bool, default false), complete_species: [codes] (only species
where the source states a genome-wide search AND ≥90% of its members map uniquely), decoys, distinct_from: [{family,
at_level, src}], tier_summary, protocol_version: 1, notes`.

## 2. Reference tree encoding (NHX)
Leaves = UniProt accessions. Internal nodes: `[&&NHX:ev=S|D|U|H:id=<node_id>:tier=T1|T2|T3|T3fig:src=<source_id>]`, plus for
duplications `lvl=<level>` (exact) or `lvl_lo=<level>:lvl_hi=<level|ABOVE_RANGE>` (youngest / oldest ancestral lineage in
which the duplication may have happened), `tag=WGD|TAND|LSE|HGT|ALLO|CONV|...`.
- **S polytomy** = all pairs whose LCA is this node are orthologs (allowed only where ≤1 member per species below it AND
  the lineage did not pass through a WGD unless the set is synteny-curated — see rule 3 of §3).
- **D** = the children are paralogs from one duplication (or a WGD event). Several copies → multifurcating D.
- **U** = unresolved: pairs whose LCA is a U node are not scored.
- **H** = HGT: the child subtree flagged `recipient=1` was acquired by transfer; give `donor=<level or OUTSIDE>`.
- Levels must exist in `panels/<P>/levels.tsv`. Express literature events by **panel branch**, not by name; if the
  literature timing is contested (2R vs lamprey, plant WGD placement), use an interval.

## 3. Compilation rules for tiers
- **T1 (text)**: the paper explicitly states the relation/event. Quote required.
- **T2 (classification)**: curated labels (subfamily tables, pillars, L/S names). Rule set: (a) genes with different labels
  whose split the paper dates (quote) → D at that level/interval; (b) same label, ≤1 gene per species → S polytomy, except
  across WGD-descended lineages where a same-label set is not synteny-curated → U (pseudo-orthology risk); (c) same label
  and >1 gene in a species → U unless the paper resolves it.
- **T3 (tree)**: published tree clades with support ≥70% BS / ≥0.95 PP; events by species overlap; label `tier=T3`.
- **T3fig**: read from a figure image; requires two independent extractions that agree; never headline.

## 4. Verification (done by a different agent, blind)
`scripts/verify_provenance.py` (mechanical: quote found in stored source, file hashes match, IDs in panel FASTA,
thresholds) must pass. Then a verifier agent receives only (paper files, locator, gene IDs) for min(all, 30) seeded
assertions + all D-levels, decoys and `distinct_from`, re-derives them, and reports discrepancies.
Kill rule: ≥2/30 discrepancies → re-curate or demote to T3; 1 → resample 30 more.

## 5. Freeze
At T+12 h: `python -m hobench.freeze` writes `families/MANIFEST.sha256`; git tag `refs-frozen-v1`. Post-freeze changes only
via `families/ERRATA.md` with before/after numbers.
