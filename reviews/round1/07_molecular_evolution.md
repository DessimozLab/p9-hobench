# Review R1-07 — Molecular evolution / phylogenetics

1. **ME-1 (MUST) Interval-valued duplication levels.** `D` gets `lvl_min`/`lvl_max`, so the reference merge level is an interval [lo, hi]. A prediction inside it is correct, and per-level F1 skips the levels inside it. Needed for ME-2/3/7–9.

2. **ME-2 (MUST) Yeast WGD is allopolyploid (Marcet-Houben & Gabaldón 2015; Wolfe 2015), so "D at the WGD node" is wrong.** Ohnologs coalesce in the parental lineages, partly before the KLE/ZT split. Assert:
   - Ohnolog pairs: separate at the post-WGD node and below; merge ∈ [ZT+WGD, Saccharomycetaceae].
   - Post-WGD vs KLE: orthologs at Saccharomycetaceae (safe).
   - Post-WGD vs ZT: unasserted at the ZT+WGD level, unless an MH2015 topology class places that copy with ZT.
   - YGOB pseudo-orthologs (Scannell 2006): paralogs, merge ≥ ZT+WGD.
   - Converted ohnologs (e.g. RP pairs; Evangelisti & Conant 2010): CONV stratum (ME-10).

3. **ME-3 (MUST) Salmonid LORe.** Robertson 2017 defines LORe relative to the Coregoninae/Thymallinae split. With only salmon, trout and pike in the panel, LORe and AORe both collapse onto the Salmoninae stem. Fixes:
   - Add *Coregonus* or *Thymallus*.
   - Encode per gene from the MOESM2 trees. AORe: `D` on the Salmonidae stem. LORe: an independent `D` per lineage, so cross-lineage pairs are co-orthologs.
   - Exclude genes with ambiguous trees and rediploidisation calls the panel cannot resolve.
   - Do not pool LORe genes into one "Ss4R" event.

4. **ME-4 (MUST) Correct §2.2: HOG semantics do define xenolog levels.**
   - Xenolog vs native recipient homolog: separate at every recipient-lineage level; merge at LCA(donor, recipient), or ∞ if the donor is not in the panel.
   - Xenolog vs donor-lineage gene: same HOG at that LCA; exclude from ortholog/paralog classification.
   - URA1 has no donor in the panel and pairs distant DHOD classes, so it is an easy ∞ case; tag it.
   - Replace GAL1/GAL3 (already in YGOB) with the *S. pombe* GAL cluster HGT (Slot & Rokas 2010; donor lineage in the panel).
   - Add EGT families to E: GapA/B vs GapC, organellar aaRS (donors BRADU/SYNY3). There are none now.

5. **ME-5 (MUST) The β-lactamase entry is factually wrong.** Classes A, C and D are homologous serine PBP-fold enzymes (Hall & Barlow 2004). Only class B (MBL fold) is unrelated. The NCBI catalogue is allele nomenclature for plasmid genes, mostly absent from K. Replace it with non-homologous isofunctional decoys (Omelchenko 2010):
   - LysRS class I vs class II
   - FBP aldolase class I vs class II
   - Cu/Zn vs Fe/Mn SOD

6. **ME-6 (MUST) Declare homology between reference families.** Rab/Arf, kinesin/myosin (Kull 1996), tubulin/FtsZ, MCM vs archaeal MCM, and EF-Tu/EF-G are homologous pairs. The E/K roots are LUCA-level, so merging them at the root is HOG-correct.
   - Count contamination only for merges with declared non-homologous decoys, or merges younger than the documented split.
   - Include all trGTPases (IF2, LepA, BipA, SelB) in the EF family; otherwise it is a trivial ∞ test.

7. **ME-7 (MUST) Restrict the T2 rule "≤1 gene/species ⇒ S".** It fails after reciprocal ohnolog loss: yeast, teleosts (Sémon & Wolfe 2007), salmonids.
   - Teleost a/b suffixes are not cross-species homologous unless synteny-curated, as for Hox (Kuraku & Meyer 2009).
   - Braasch ST23 anchors only to gar, so zebrafish–medaka singletons get merge ∈ [Teleostei, Neopterygii].

8. **ME-8 (SHOULD) Vertebrate 2R.** Model: 1R shared, 2R gnathostome-only, plus a cyclostome triplication (Simakov 2020; Nakatani 2021; Yu 2024 *Nat Ecol Evol*; Marlétaz 2024 *Nature*). Lamprey–gnathostome pairs get merge ∈ [Vertebrata, Olfactores], reported separately. Ignore names (lamprey Hb ≈ Cygb; Hoffmann 2010).

9. **ME-9 (SHOULD) Place plant WGDs by panel branch, not by name.**
   - ζ/ε are contested (Ruprecht 2017), and the panel has no gymnosperm: record one branch event, with no WGD claim.
   - ρ/σ/τ collapse onto the Poaceae stem; α/β onto the Brassicaceae stem.
   - Add *Phaseolus* to split soybean's 59 and 13 Ma WGDs.
   - γ is possibly hybrid (Aköz & Nordborg 2019); it needs an interval if AQUCO is added.
   - Aggregate event recovery per branch.

10. **ME-10 (SHOULD) Homeologs and gene conversion.**
    - *X. laevis*: `D` on the XENLA leaf is correct (progenitors diverged ~34 Ma, after the split from *X. tropicalis* ~48 Ma; Session 2016). XENTR vs L/S are co-orthologs; tag ALLO.
    - Maize: same encoding; needs an interval if sorghum is added (Swigoňová 2004). Do not add wheat.
    - Add a CONV stratum, excluded from the headline: yeast RP ohnologs, *E. coli* tufA/tufB (Lathe & Bork 2001), OPN1LW/MW, teleost RH2. Identical copies (OPN1MW/MW2/MW3) are unmappable.

11. **ME-11 (SHOULD) ILS.** HOGs of single-copy genes are ILS-invariant; the scoreable error is spurious duplication (Hahn 2007).
    - Add a spurious-D rate on confirmed single-copy groups across short internodes (*mel/sim/sec/yak*; Pollard 2006).
    - Give intervals to duplications next to short internodes (hemiplasy; Rasmussen & Kellis 2012).
    - Panel D has no §4 family: add the Good 2014 P450s or drop it.

12. **ME-12 (SHOULD) Poorly suited families.**
    - KRAB-ZNF: `age_MA` is orthology-inferred (Imbeault 2017), so circular; human-centric; repeat arrays. Use for integrity only.
    - NLR: restrict to RNL (ADR1/NRG1), ZAR1 (Adachi 2019) and TNL loss in monocots (Collier 2011).
    - Fusions (Pol B, APOBEC3): declare an anchor domain; leave pairs with the second parent family unasserted, not counted as contamination.

13. **ME-13 (COULD) Coverage.**
    - Make genome-scale sets (YGOB, ST23, Robertson, Xenbase) a separate stratum so they do not swamp the macro-averages.
    - E is over-weighted (6/10 LECA cytoskeleton/GTPase families): swap dynein for EGT.
    - Check K nitrogenase coverage first (likely 3–4 species).
