# Curation log: Y04_GAL_Spombe_HGT (panel Y)

Status: **abandoned** (no source identifiers for the HGT-relevant members; the HGT signal is not testable in panel Y;
the WGD part duplicates Y01). Time ~25 min.

## Read
- Slot & Rokas 2010 PNAS 107:10136 (PMID 20479238, PMC2890473; verified: PNAS, not OA subset -> PMC HTML text only).
  Key statements: "the GAL cluster of Schizosaccharomyces yeasts was acquired through horizontal gene transfer from a
  Candida yeast"; "All four gene phylogenies robustly supported the placement of Schizosaccharomyces sequences within the
  Candida clade"; "Schizosaccharomyces species also contain an additional unclustered GAL10 paralog".
  Gene identifiers are only in SI Table S1/S3 (0914418107_pnas.200914418SI.pdf), which PMC serves behind the
  proof-of-work wall (NOT RETRIEVABLE by script); main text and Table 1 (Cryptococcus only) contain no IDs.
- Hittinger & Carroll 2007 (GAL1/GAL3) is paywalled; not used.

## Why abandoned
1. HGT not testable in panel Y: the donor is a Candida (CUG-clade) yeast, represented in panel Y by CANAL, and the only
   recipient is SCHPO. With donor = CANAL the H-node merge level is lca(CANAL, SCHPO) = Ascomycota, the panel root, which
   equals the species LCA of every SCHPO-vs-Saccharomycotina pair. The xenolog relation is therefore indistinguishable
   from orthology in this panel (no xenolog-separation signal); native S. pombe GAL1/GAL7 copies no longer exist to be
   separated from. Only the native unclustered GAL10 paralog (uge1) would give a contrast, and no source places it.
2. Member identifiers: none in the retrievable text. SCHPO gal1/gal7/gal10 (Q9HDU2/Q9HDU5/Q9HDU3) and CANAL GAL1/7/10
   could only be added by gene-name lookup without a source table tying them to the paper's phylogenies.
3. The GAL1/GAL3 WGD ohnologs (YBR020W/YDR009W), GAL7 and GAL10 in Saccharomycetaceae are already units of
   Y01_YGOB_pillars (YGOB_Anc_3_219, _217, _218); a separate family would duplicate gene ownership.

## If revisited
Retrieve the PNAS SI PDF manually (Table S1/S3 gene IDs) and add SCHPO/CANAL members to the Y01 GAL units with an
H node (donor=CANAL) only in a panel where the donor-recipient LCA is below the root (e.g. a panel with several
Taphrinomycotina and Pezizomycotina species). Not requested in requested/NEEDED_BY_CURATORS.md because the panel-Y
signal is degenerate even with the IDs.
