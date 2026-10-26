# Core formalism (PI draft v0, to be folded into DESIGN.md)

## Objects
- **Panel**: a set of species S with a rooted reference species tree T (Newick). Every internal node of T has a unique
  name = taxonomic level (e.g. `Vertebrata`, `Euteleostomi`). Leaves = species codes (UniProt mnemonic style, e.g. `HUMAN`).
- **Gene**: a protein ID (UniProt accession preferred) belonging to exactly one species.
- **Prediction** (from any tool): for each level L of T, a partition of (a subset of) the genes of species below L into
  HOGs. Assumed nested (a HOG at L is contained in exactly one HOG at parent(L), if any). Genes may be absent (singletons).
  Parsers: OrthoXML (FastOMA, OMA), OrthoFinder `Phylogenetic_Hierarchical_Orthogroups/N*.tsv` (node names mapped to our levels).
- **Reference family** F (from literature): a set of member genes M_F (in panel species) plus a rooted, possibly
  multifurcating *reference gene tree* R_F whose internal nodes carry `event ∈ {S, D, U}` (speciation, duplication,
  unresolved/unknown), optional explicit `level` for duplications (the species-tree node in whose ancestral lineage the
  duplication took place; must be ⊇ the LCA mapping), optional tags (WGD, HGT, ILS, ...), and provenance
  (paper DOI, figure/table, quote). Unresolved parts → no assertion.

## Merge level (the unifying quantity)
For genes a, b with species LCA S0 = lca_T(sp(a), sp(b)):
- reference: let v = lca_R(a, b).
  - v is S → `merge(a,b) = S0` (orthologs).
  - v is D at level D → separate at all levels from S0 up to D inclusive; together above D → `merge(a,b) = parent_T(D)`
    (if D is the family's root level and no parent level exists inside the family → `merge = ∞`, i.e. never in one HOG).
  - v is U → pair not asserted (excluded).
- prediction: `merge_pred(a,b)` = lowest level L ⊇ S0 at which a and b are in the same predicted HOG; ∞ if none.
Two genes are predicted orthologs iff merge_pred = S0. They are in the same HOG at level L iff L ⊇ merge.

## Metrics (first set)
1. **Pair merge-level accuracy**: fraction of asserted pairs with merge_pred == merge_ref; signed level error
   (number of species-tree edges; + = predicted too old (under-split), − = too young (over-split), ∞ handled as category).
2. **Per-level HOG pair F1**: at each level L, precision/recall of "same HOG at L" over asserted pairs of genes below L.
   Gives a depth profile.
3. **Orthology (QfO-style) P/R**: ortholog vs paralog classification derived from merge levels.
4. **Event recovery**: for each literature-asserted duplication e (level D, copies C1..Ck): fraction of cross-copy pairs with
   correct merge level; event recovered if ≥ 0.5 (also report continuous); misplacement distance in levels.
5. **Family integrity**: #predicted root-level HOGs needed to cover ≥ 90% of M_F (split); fraction of genes in those HOGs
   belonging to *other* reference families or to declared non-homologous decoys (merge); coverage (members in any HOG).
6. **Gene-level**: per gene, fraction of its asserted pairs correct → per-species "correctly placed genes".
Aggregation: macro-average over families (each family equal weight) and per phenomenon tag; bootstrap over families for CIs.

## Controls
- Positive: prediction = reference-derived HOGs → all metrics perfect.
- Negative: random relabelling, "everything one HOG at root" (all merge at S0), "all singletons" → known values.
- Dose-response: reference with k random leaf swaps / collapsed duplications → monotone degradation.
