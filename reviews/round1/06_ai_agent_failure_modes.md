# Review round 1 — AI-agent failure modes (reviewer 06)

**Evidence of live failure modes already present**
- Committed `literature/raw/C/PMC6291210.xml` is an HTTP-500 JSON body.
- `literature/B1_*.md` states counts "from general knowledge".
- DESIGN §4 lists "β-lactamase classes A/B/C/D (non-homologous ⇒ merge decoys)". Wrong: A, C and D are homologous serine β-lactamases; only B is unrelated.
- `hobench/reference.py::validate_family` only *warns* when a quote is missing. It never checks quotes against stored text or `gene_id` against the panel proteome.
- `hobench/evaluate.py` has `--skip-validation`, uses `reconcile(check=False)`, and truncates `prediction_warnings` at 200. `n_members_absent_from_prediction` is computed but never gated.

## Changes

**AI-1 (MUST) — Verify quotes and locators mechanically.** Add `curation_tools/verify_provenance.py` and call it from `validate_family`, with failures as errors (not warnings). Every members.tsv row and NHX `src` carries `(paper_id, file, sha256, locator, quote)`.
- It re-hashes the file against `meta.json`.
- It matches the quote verbatim after normalisation (Unicode NFKC, ligatures, dashes, whitespace). Anything below exact match fails.
- For T1 assertions, the quote must contain the gene label or ID it supports.
- For table rows, the tool re-reads (file, sheet, row, column) and checks that the cell equals `original_id_in_paper` and `label`.
- This takes about 1–2 h of work.

**AI-2 (MUST) — Extract by script, not by transcription.** T2 members must come from a committed parser (`families/<ID>/build.py`) that reads checksummed supplement files. LLM-typed IDs are allowed only for T1 assertions that pass AI-1. Parsers assert the member counts stated in the paper text.

**AI-3 (MUST) — Check that members exist in the panel.** In the validator, every `gene_id` must be present in the hashed panel FASTA of its `species`. Mapped rows must also respect the DESIGN §5.3 thresholds (≥95% identity, ≥80% coverage, `status ∈ {exact, unique}`). Today `map_by_sequence` defaults to `--min-qcov 0.5` and `map_ids --taxon` is optional.

**AI-4 (MUST) — Ban circular sources explicitly, and audit for them.** Extend DESIGN §4/§10 with a written denylist:
- Ensembl Compara/HCOP-derived data (OHNOLOGS v2 is built on Compara paralogs).
- Gene names in UniProt unreviewed entries, RefSeq `LOC`/"-like" names, and Ensembl projected names.
- OMA, OrthoDB, eggNOG, PANTHER and TreeFam.

HGNC groups and OHNOLOGS may enter only as a separately reported "external-DB" tier, never as headline truth. `scripts/audit_sources.py` greps family sources and build scripts for denylisted hosts.

**AI-5 (MUST) — Ground every event level and tag.** Every D node `lvl=`, every phenomenon tag, and every decoy declaration needs a quote-verified `src`. If none exists, the item becomes `U`/untagged; nothing comes from "known to pre-date L". Fix the β-lactamase decoy now: only class B against A/C/D.

**AI-6 (MUST) — Make verification independent and size it statistically.** This replaces "random 10%, >5%" in §5.5.
- The verifier sees only `(paper, locator, gene ID)`. It re-derives the label, relation or level blind to the curator's answer, and a script compares the two.
- Sample n = min(all, 30) assertions per family (seeded script, `SEEDS[1]`, stratified by tier/type), plus all D-node levels and decoys.
- Kill rule: ≥2 discrepancies out of 30 → re-curate or demote to T3. One discrepancy → resample 30 more.
- Report the per-family Clopper-Pearson upper bound.

**AI-7 (MUST) — Enforce blinding and freezing mechanically.**
- Tool outputs go to a directory curators are told not to read. Afterwards, grep curator logs for `runs/`, `evaluate` and denylisted URLs.
- Before the first real scoring, commit `families/MANIFEST.sha256` and tags `refs-frozen-v1` and `metrics-frozen-v1` (metrics.py, metric list, aggregation, thresholds).
- `evaluate.py` records the git commit and the families/proteome hashes. Analysis scripts refuse mismatches.
- After the freeze, any fix is an erratum entry, made by an agent blind to scores. Both pre- and post-fix numbers are reported.

**AI-8 (MUST) — Remove silent fallbacks in the scoring path.**
- Drop `--skip-validation`, or make it write `UNVALIDATED` into the output.
- Gate `n_members_absent_from_prediction`: above 1% is a hard error (an ID-normaliser failure would look like over-splitting).
- Store full warning counts (not truncated), plus `induced_levels`, `membership_conflicts` and lenient-OrthoXML fallbacks in the results TSV.
- Make fetch/survey scripts reject error bodies and texts under `MIN_TEXT`, and delete PMC6291210.xml.

**AI-9 (MUST) — Add a licence gate before commits.** `fetch_paper` writes non-CC files to `literature/papers/<id>/raw/`, which `.gitignore` does not cover. `literature/raw/B1/{kinbase,hgnc,ohnologs}` are untracked but not ignored. Fix: ignore `literature/papers/*/raw/` and add `scripts/check_licences.py`. It blocks files with `redistributable=false` or no DATA.md entry.

**AI-10 (SHOULD) — Treat figures as a separate tier.** Assertions read from figure images become `T3-fig` and stay out of the headline. They require two independent extractions that agree exactly; disagreement → `U`.

**AI-11 (SHOULD) — Make the tests adversarial.** Tests for the metric code are written or reviewed by a different agent from its author. Add a small mutation check: flip the sign, swap S/D, drop the ∞ handling, and confirm the suite fails each time. Run positive/negative controls on every frozen family as a gate.

**AI-12 (SHOULD) — Keep claims traceable to generated tables.** Every number in REPORT.md comes from tables in `results/` produced by scripts and carries a commit hash. A fresh agent checks every numeric and "X better" claim against them and the paired-bootstrap CIs. Per-tag claims with n < 5 families are labelled descriptive.

**AI-13 (SHOULD) — Prevent context loss and scope creep.** Write `design/CURATION_PROTOCOL.md` (tiers, the AI-1–AI-5 rules, denylist), pass it verbatim to every curator, with its hash recorded in each family's `curation.yaml` (agent/model, source hashes, decisions).
- Set an MVP priority order: Y, K, E, V, then P; D and E3 last.
- Kill rule: at T+10h, freeze whatever has passed AI-6. Never relax thresholds to raise the family count.

**AI-14 (COULD)** — For 3 families, a second curator builds the reference independently; report pair-level agreement as curation noise.
