# HOBench: a literature-based benchmark for hierarchical orthology inference

> **Read this first.** This repository is the output of an experiment that tested whether AI agents can run a
> bioinformatics research project. An AI agent acted as principal investigator and directed AI subagents through
> design, literature curation, software, tool runs, analysis, the accompanying manuscript and this repository. A human
> wrote the prompts (`prompts/`), supplied papers that the agents requested, and restarted the session after usage
> limits. No human checked the curated data or the code line by line. Treat the benchmark as a research prototype: the
> process record (`docs/process/`, `docs/reviews/`) and the error logs (`families/ERRATA.md`, `results/ERRATA.md`)
> are part of the deliverable.

HOBench scores hierarchical orthologous groups (HOGs), or equivalently reconciled gene trees, against evolutionary
histories that publications assert. It provides:

- **5 species panels** (`panels/`): yeasts (Y, 16 species), prokaryotes (K, 21), eukaryotes with prokaryote outgroups
  (E, 34), vertebrates (V, 17) and plants (P, 14). Proteomes are UniProt Reference Proteomes release 2026_03. Each
  panel has a binary species tree with named internal levels (`species_tree.nwk`, `levels.tsv` with citations).
- **Reference families** (`families/`): 32 curated entries from published phylogenies and classifications, with
  member lists, a reference tree with speciation, duplication, HGT and unresolved nodes, a quote for every asserted
  event, mapping reports and blind-verification results. They are frozen (`families/MANIFEST.sha256`). Exports in
  OrthoXML and as asserted gene pairs are in `references_export/`.
- **Metrics** (`hobench/`): per-level HOG B-cubed F1, merge-level accuracy, two-sided duplication-event recovery and
  family-membership F1, each reported next to trivial baselines (all singletons, one big HOG, random relabelling).
- **The paper's runs** (`runs/`): compressed outputs of OrthoFinder 3.1.5 and FastOMA 0.5.0, and the result tables
  (`results/`).

## Quick start: score your method

```bash
conda env create -f environment.yml && conda activate hobench
# 1. proteomes (exact files used in the paper; download the archive from the release page, see below)
python scripts/download_proteomes.py --archive hobench_proteomes_uniprot2026_03.tar.gz
python scripts/build_subsets.py                     # subset proteomes for panels E, V, P (and Y, K)
# 2. run your method on panels/<P>/proteomes/*.fa (full mode) or panels/<P>/subset/*.fa (subset mode),
#    using panels/<P>/species_tree.nwk as the species tree
# 3. score it
python scripts/score_method.py --panel K --mode full --name mytool --format orthoxml --input my_hogs.orthoxml.gz
cat results_user/mytool/K_full/summary.md
```

Use full mode for Y and K and subset mode for E, V and P; these are the modes used in the paper. Accepted formats are
OrthoXML (levels matched to panel levels by species content), an OrthoFinder 3 `Phylogenetic_Hierarchical_Orthogroups`
directory, or a TSV with columns `level, hog_id, gene`. Gene ids must be UniProt accessions or UniProt-style headers.
See `python scripts/score_method.py --help`.

Without the archive, `download_proteomes.py` fetches the proteomes from UniProt. This reproduces the paper only while
UniProt still serves release 2026_03; the script checks every file against its sha256 and refuses mismatches.

## What the numbers mean
For two genes, the *merge level* is the most recent taxonomic level at which they belong to the same HOG. Each
literature reference gives an interval of allowed merge levels for each asserted gene pair. Pairs that the sources do
not resolve are not scored. The four primary endpoints are explained in the manuscript and in `docs/DESIGN.md` §6.
Read every endpoint together with the trivial baselines in `summary.md`: for example, "one big HOG" scores well on
merge-level accuracy in families that consist mostly of speciation nodes.

## Repository layout
| Path | Content |
|---|---|
| `hobench/` | formats, parsers (OrthoXML, OrthoFinder, TSV), metrics, controls, statistics, freeze |
| `scripts/` | `score_method.py`, `evaluate_panel.py`, confirmatory/equivalence/robustness tests, analysis, audits |
| `wrapper/` | runs OrthoFinder/FastOMA with manifests; subset builder; perturbation generator |
| `curation_tools/` | paper retrieval (Europe PMC), identifier extraction and mapping, proteome download |
| `families/` | frozen references; `ERRATA.md`, `VERIFICATION.md`, `AUDIT.md` |
| `panels/` | species tables, trees, levels, subset ID lists |
| `runs/`, `results/` | the paper's tool outputs and result tables |
| `docs/` | design, curation protocol, data sources (`DATA.md`), process record, design reviews (process files are the
live research record and contain intermediate numbers; the manuscript and `results/` supersede them) |
| `prompts/` | the human prompts that started the project and the manuscript |
| `tests/` | unit tests with hand-derived expected values (`pytest -q tests`) |

## Reproducing the paper
`environment-tools.yml` lists the tools. The runs were made with `python -m wrapper.run --tool orthofinder|fastoma
--panel <P> --mode full|subset` (FastOMA needs the OMAmer database `LUCA.h5` from https://omabrowser.org/All/LUCA.h5;
see `docs/DATA.md`). Scoring: `python scripts/evaluate_panel.py --panel <P> --mode <mode>`; analysis:
`bash scripts/analysis/run_all.sh`. Rebuilding references (`families/<ID>/build.py`) requires the source papers and
supplements listed in each `family.yaml`; these are not redistributed.

## Release assets (to be uploaded by the maintainer)
`hobench_proteomes_uniprot2026_03.tar.gz` (469 MB; 90 proteomes; sha256
`d07d917fa72bf54d14d98d02e18bd425a08731b72296a88743c74303cd343561`).

## Licences
Code: MIT (`LICENSE`). Proteomes: UniProt, CC BY 4.0. Reference families contain identifiers, labels and short quotes
from the cited publications (see `DATA_LICENSE.md`).

## Citation
See `CITATION.cff` and the manuscript (`manuscript/`).
