# ERRATA / freeze log for reference families

## Freeze refs-frozen-v1 (2026-10-27, T+12h)
- Pre-freeze disclosure: the Y and K vertical-slice results (results/Y_full, results/K_full, commit e0f02f4 and the
  metric-defect re-run) were seen by the PI before this freeze, as planned in DESIGN §11 (F1 gate). No reference family
  was changed after those results except Y03 assertion A4 (demoted to supporting T3 on the round-2 audit's request;
  it does not change any scored relation).
- Blind verification (families/VERIFICATION.md): V01_globins FAILED the kill rule (2/37: curator used the 2R interval
  [Gnathostomata, Vertebrata] where the source dates single branches). V01 was already non-headline (77% mapping); it
  is kept as secondary, not re-curated. V02_opsins passed after a resample (1/45).
- Non-headline (mapping <80% or verification): E06_myosins, E07a, E07b, P01_homeobox, V01_globins, V04_Sox.
- Abandoned (logs only): K03, Y02 (kept as cross-check), Y04, E09, V05 (Compara-derived source), V06, V08, V09.

Post-freeze changes must be logged below with before/after metric values.

## 2026-10-27 — V07_Xenopus_LS (collection XLAEVIS)
- Format fix: `mapping_coverage` was percent (e.g. 100.0, 96.4); now a fraction (1.000, 0.964) per hobench/panel.py.
  Identity values and all mapping decisions unchanged.
- Dropped 2 units whose genes are also members of V03_Hox (XL_hoxa4, XL_hoxb7); rule: drop any V07 unit sharing a gene
  with another family. Units 500 -> 498, members 1,500 -> 1,494, T2 D(ALLO) 500 -> 498, T2 S 500 -> 498. No re-sampling.
- Before/after metric values: none — no metric values for panel V had been computed/seen by the curator.
