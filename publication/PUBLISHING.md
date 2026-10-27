# How to publish (for the human maintainer)

1. **Repository** — `publication/hobench/` is a standalone git repo (one commit, no remote).
   ```bash
   cd publication/hobench
   git remote add origin git@github.com:<you>/hobench.git
   git push -u origin main   # or master, see `git branch`
   ```
   Size ~103 MB (largest files < 20 MB). Regenerate it from the workspace with `bash publication/build_repo.sh`.
2. **Proteome archive** — upload `publication/release_assets/hobench_proteomes_uniprot2026_03.tar.gz` (469 MB,
   sha256 in `release_assets/SHA256SUMS`) as a GitHub release asset or to Zenodo, then put its URL/DOI in
   `README.md` ("Release assets") and in the manuscript ("Data and code availability").
3. **Manuscript** — `publication/manuscript/hobench_manuscript.tex` (+ S1_Text.tex, S2_Text.tex). Fill in authors,
   affiliations, repository URL and archive DOI; decide on authorship/AI-use statements per PLOS policy. Compile with
   `tectonic -X compile hobench_manuscript.tex` (or pdflatex + bibtex). For submission, PLOS wants figures uploaded
   separately (TIFF/EPS) and removed from the .tex; `make_figures.py` regenerates them.
4. **Before publishing**, consider a human spot-check of a few reference families (`families/<ID>/curation_log.md`,
   `assertions.tsv`), since no human reviewed the curated data line by line.
