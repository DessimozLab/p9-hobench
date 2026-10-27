# HOBench-core data release

`MANIFEST.tsv` lists every released file (kind, panel, path, bytes, sha256): full and subset input FASTAs of every
panel (including the robustness variants `Y_*`), subset member tables, species trees, `levels.tsv`, `species.tsv`,
perturbation tables and the reference export (`references/<panel>/`).

* `references/<panel>/reference_hogs.orthoxml.gz` - asserted structure only, OrthoXML 0.5 (U subtrees are separate
  top-level groups with `<property name="unresolved" value="true"/>`; interval events carry `lvl_lo`/`lvl_hi`/`merge_lo`/`merge_hi`).
* `references/<panel>/reference_hogs.positive_control.orthoxml.gz` - same, but U subtrees grouped at the level the
  positive control uses; reads back with `hobench.predictions.read_orthoxml` (used for the round-trip check).
* `references/<panel>/asserted_pairs.tsv.gz` - `gene_a gene_b relation(ortholog|paralog|xenolog) merge_lo merge_hi tier family`.

FASTAs are **not committed to git**. Regenerate the export, manifest and file list, then build the archive:

```bash
export PATH=$PWD/conda_envs/ortho/bin:$PATH
python scripts/export_references.py --validate      # release/references/*
python scripts/make_release.py                      # release/MANIFEST.tsv, release/file_list.txt
tar --sort=name --mtime='2026-10-27 00:00Z' --owner=0 --group=0 --numeric-owner \
    -czf HOBench-core-v0.2.1.tar.gz -T release/file_list.txt
sha256sum HOBench-core-v0.2.1.tar.gz > HOBench-core-v0.2.1.tar.gz.sha256
```

Verify an unpacked archive: `cd <unpacked> && tail -n +2 release/MANIFEST.tsv | awk -F'\t' '{print $5"  "$3}' | sha256sum -c -`
(the three trailing entries of `file_list.txt` - MANIFEST.tsv, DATA.md, LICENSE - are not listed in the manifest itself).
Licences: sequences UniProt CC BY 4.0; references are literature-derived facts; code MIT. See `DATA.md`.
