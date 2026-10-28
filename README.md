(Human-written text starts) This is an exercise in assessing the bioinformatics capabilities of agentic AI. It is part of a hackathon run by the Beltrao lab at ETH Zurich, where we used Claude to run bioinformatics projects [LINK TO BE ADDED LATER]. All results in this project were generated entirely by AI with minimal human intervention. This is not intended to be viewed as a traditional publication and bioinformatics tool.

The goal of this specific project was to create a set of species, gene families, and accuracy metrics for benchmarking hierarchical orthology. It was executed by Claude Opus 5.5 with delegation to Opus 5.5/Sonnet 5.5 subagents, and finished within 24 hours on a cluster.

Currently none of the results have been verified by a human. We plan to verify the agent's methodology and data curation, and will update this section accordingly.

The rest of this repository is entirely AI generated. It includes code for data retrieval and curation, and a [AI-written report](https://github.com/DessimozLab/hobench/blob/master/publication/manuscript/hobench_manuscript.pdf) in the style of a journal article, detailing the final method. (Human-written text ends)

# HOBench

HOBench scores hierarchical orthologous group (HOG) inference (OrthoFinder, FastOMA, anything that writes OrthoXML)
against reference gene-family histories asserted in the literature, not against another tool's output.
Five species panels (Y yeasts, K prokaryotes, E eukaryotes, P plants, V vertebrates) fix the proteomes (UniProt Reference
Proteomes 2026_03), a rooted species tree and its named levels. Each reference family is a species-tree-reconciled gene tree whose nodes are
S, D, H (transfer) or U (unresolved, not scored), with duplication levels given as exact levels or intervals and every assertion tied to a source.
Four pre-registered endpoints: per-level HOG B-cubed F1, node-balanced merge-level accuracy, two-sided duplication-event recovery, family-membership F1.
Controls (positive, singletons, one-HOG, relabel) and robustness perturbations are built in; no composite score.

## Layout
| Path | Content |
|---|---|
| `hobench/` | formats, parsers, metrics, controls, freeze, `python -m hobench.evaluate` |
| `families/<ID>/` | frozen references: `family.yaml`, `reference_tree.nhx`, `members.tsv`, `assertions.tsv`, `build.py` |
| `panels/<P>/` | `species.tsv` (pinned sha256), `species_tree.nwk`, `levels.tsv`, `subset/` (family-neighbourhood subsets) |
| `wrapper/` | run tools on a panel with run manifests (`run`), subset builder (`subset`), perturbations |
| `runs/<P>/<mode>/<tool>_<variant>/` | archived tool outputs + `manifest.json` |
| `scripts/` | `evaluate_panel.py` (score all runs of a panel), `export_references.py`, `make_release.py`, audits |
| `results/`, `release/` | result tables; reference export + `MANIFEST.tsv` |
| `tests/` | `python -m pytest tests/` |

## Quickstart: score your own output
Setup: `export PATH=$PWD/conda_envs/ortho/bin:$PATH` (or any Python with `requirements.txt`). Run your tool on the panel
proteomes (`panels/K/proteomes/*.fa`, full mode; `panels/K/subset/` for subset mode) with the panel species tree
(`panels/K/species_tree.nwk`; OrthoFinder needs `-s`, FastOMA `--species_tree`). Then, for panel K:

```bash
# OrthoFinder (3.x) results
tar -xzf runs/K/full/orthofinder_default/of_core.tar.gz -C /tmp/of_K       # or your own results directory
python -m hobench.evaluate --panel-tree panels/K/species_tree.nwk --families families --panel K --panel-dir panels/K \
   --prediction /tmp/of_K/Phylogenetic_Hierarchical_Orthogroups --format orthofinder \
   --of-species-tree runs/K/full/orthofinder_default/SpeciesTree_rooted_node_labels.txt --out /tmp/K_of

# any OrthoXML (FastOMA shown; --tool-species-tree = the tree the tool used; levels are matched by clade species set)
python -m hobench.evaluate --panel-tree panels/K/species_tree.nwk --families families --panel K --panel-dir panels/K \
   --prediction runs/K/full/fastoma_default/FastOMA_HOGs.orthoxml.gz --format orthoxml \
   --tool-species-tree runs/K/full/fastoma_default/species_tree_checked.nwk --out /tmp/K_fastoma
```
Output: `<out>.json` (per family, bootstrap CIs, provenance) and `<out>.tsv`. Add `--no-skill` for speed, `--flat` for flat orthogroup TSVs
(`--format tsv`), `--dup-placement late` for pyHam-style duplication placement. Inconsistent input (unknown levels or species, genes in two groups)
is a hard error, never repaired. Whole-panel scoring with controls and statistics: `python scripts/evaluate_panel.py --panel K --mode full`
(writes `results/K_full/`). Sanity check: the exported positive control scores 1.0 (`release/references/K/reference_hogs.positive_control.orthoxml.gz`).

## Run the tools via the wrapper
```bash
python -m wrapper.run --tool orthofinder --panel K --mode full --variant default --threads 6   # also: fastoma, of_flat
python -m wrapper.run --tool fastoma --panel Y --mode subset --dry-run                          # print the command only
```
The wrapper never reads `families/`; work dirs go to `$HOBENCH_WORK` (default `/tmp/hobench_work`), final outputs and `manifest.json` to `runs/`.

## Rebuild the references
Each family is built offline from its cached literature tables: `python families/<ID>/build.py` (from the project root), then
`python -m hobench.freeze --check` against `families/MANIFEST.sha256`. References are frozen: do not edit `families/` after a run;
changes go through `families/ERRATA.md`. Re-export with `python scripts/export_references.py --validate` and `python scripts/make_release.py`.

## Licences
Code: MIT (`LICENSE`). Sequences: UniProt Reference Proteomes, CC BY 4.0 (FASTAs are distributed only in the release archive, not in git;
checksums in `panels/*/species.tsv`, `release/MANIFEST.tsv`). References: facts taken from the cited papers (identifiers, relations, levels), with short
supporting quotes and locators in `assertions.tsv`; no restricted files are redistributed. Per-file details: `DATA.md`.

## Citation
See `CITATION.cff` (placeholder: HOBench contributors, 2026; add DOI/URL on release).

## More
`DESIGN.md` (goal, definitions, metrics, experiments), `REPORT.md` (results), `DATA.md` (data provenance and licences),
`design/ENGINEERING_NOTES.md` (API and scoring semantics), `release/README.md` (data archive recipe).
