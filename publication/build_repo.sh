#!/usr/bin/env bash
# Assemble the standalone HOBench MVP repository (publication/hobench) from the project workspace.
# Copies only code, frozen references, panel definitions, compressed tool outputs of the paper's runs and result
# tables. No proteome FASTAs (downloaded by scripts/download_proteomes.py or taken from the release archive),
# no paper PDFs or publisher supplements.
set -euo pipefail
SRC="$(cd "$(dirname "$0")/.." && pwd)"
DST="$SRC/publication/hobench"
mkdir -p "$DST"
rs() { python3 "$SRC/publication/copytree.py" "$@"; }

# code
rs "$SRC/hobench/" "$DST/hobench/"
rs "$SRC/wrapper/" "$DST/wrapper/"
rs --exclude 'demo_output' "$SRC/curation_tools/" "$DST/curation_tools/"
rs "$SRC/tests/" "$DST/tests/"
rs "$SRC/config/" "$DST/config/"
mkdir -p "$DST/scripts/analysis"
for f in evaluate_panel.py evaluate_robustness.py confirmatory.py equivalence.py export_references.py \
         family_catalogue.py audit_sources.py verify_provenance.py make_verification_sheet.py \
         compare_verification.py controls_gate.py mutation_check.py precompute_hogmaps.py run_omark.py; do
  cp "$SRC/scripts/$f" "$DST/scripts/"
done
rs "$SRC/scripts/analysis/" "$DST/scripts/analysis/"
[ -d "$SRC/scripts/tests" ] && rs "$SRC/scripts/tests/" "$DST/scripts/tests/"

# panels: definitions + subset ID lists (no FASTAs)
for P in Y K E V P; do
  d="$DST/panels/$P"; mkdir -p "$d/subset"
  for f in species.tsv levels.tsv species_tree.nwk species_tree_alt.nwk subset_seeds.txt; do
    [ -f "$SRC/panels/$P/$f" ] && cp "$SRC/panels/$P/$f" "$d/"
  done
  for f in manifest.json members.tsv; do
    [ -f "$SRC/panels/$P/subset/$f" ] && cp "$SRC/panels/$P/subset/$f" "$d/subset/"
  done
  mkdir -p "$d/subset/ids"
  for fa in "$SRC/panels/$P/subset/"*.fa; do
    grep '^>' "$fa" | sed 's/^>//; s/ .*//' | sort > "$d/subset/ids/$(basename "$fa" .fa).txt"
  done
done
cp "$SRC/panels/citations.tsv" "$DST/panels/" 2>/dev/null || true

# frozen references (curation build scripts read literature files that are not redistributed)
rs --exclude 'cache' "$SRC/families/" "$DST/families/"
rs "$SRC/release/references/" "$DST/references_export/"

# paper runs: compressed HOG outputs + manifests (main variants only)
for r in Y/full K/full E/subset P/subset V/subset; do
  for v in orthofinder_default orthofinder_dendroblast fastoma_default of_flat_default of_flat_dendroblast; do
    s="$SRC/runs/$r/$v"; [ -d "$s" ] || continue
    rs --exclude raw --exclude 'of_orthologues.tar.gz' --exclude '*.log' --exclude trace.txt "$s/" "$DST/runs/$r/$v/"
  done
done

# results tables (small files only)
mkdir -p "$DST/results"
( cd "$SRC/results" && find . -type f -size -5M ! -path './analysis/cache/*' -print0 ) | \
  while IFS= read -r -d '' f; do mkdir -p "$DST/results/$(dirname "$f")"; cp "$SRC/results/$f" "$DST/results/$f"; done

# documentation and process record
mkdir -p "$DST/docs/process" "$DST/docs/reviews" "$DST/prompts"
cp "$SRC/DESIGN.md" "$SRC/DATA.md" "$SRC/LITERATURE.md" "$DST/docs/"
cp "$SRC/design/CURATION_PROTOCOL.md" "$SRC/design/CURATOR_BRIEF.md" "$SRC/design/core_spec.md" \
   "$SRC/design/ENGINEERING_NOTES.md" "$SRC/design/MUTATION_REPORT.md" "$DST/docs/"
cp "$SRC/REPORT.md" "$SRC/RESEARCH_TREE.md" "$SRC/LESSONS.md" "$SRC/START.txt" "$DST/docs/process/"
cp "$SRC/requested/REQUESTED.md" "$SRC/requested/RECEIVED.tsv" "$DST/docs/process/"
rs "$SRC/reviews/" "$DST/docs/reviews/"
mkdir -p "$DST/docs/literature"; cp "$SRC"/literature/*.md "$DST/docs/literature/"
cp "$SRC/PROMPT.md" "$DST/prompts/01_project_prompt.md"
cp "$SRC/LICENSE" "$DST/LICENSE"
echo "assembled $DST: $(du -sh "$DST" | cut -f1)"

# repository-specific files (README, prompts, environments, user scripts)
cp -r "$SRC/publication/extra/." "$DST/"
cp "$SRC/scripts/evaluate_panel.py" "$DST/scripts/"
echo "final size: $(du -sh --apparent-size "$DST" | cut -f1)"

# manuscript sources, figures and PDFs
mkdir -p "$DST/manuscript/figures"
for f in hobench_manuscript.tex S1_Text.tex S2_Text.tex preamble.tex references.bib plos2025.bst make_figures.py \
         hobench_manuscript.pdf S1_Text.pdf S2_Text.pdf; do
  [ -f "$SRC/publication/manuscript/$f" ] && cp "$SRC/publication/manuscript/$f" "$DST/manuscript/"
done
cp "$SRC"/publication/manuscript/figures/* "$DST/manuscript/figures/"
