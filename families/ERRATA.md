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
