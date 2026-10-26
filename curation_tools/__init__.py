"""Curation toolkit: turn papers into reference gene families.

Modules (all runnable as ``python -m curation_tools.<module> --help``):
  fetch_paper     resolve/search/download papers (+ supplements) via Europe PMC / PMC
  extract_ids     scan supplement files for sequence identifiers
  map_ids         map identifiers to UniProtKB via the UniProt ID-mapping API
  map_by_sequence map paper sequences onto a same-species target proteome
  fetch_proteome  download a UniProt proteome with a reproducibility manifest
"""
__version__ = "0.1.0"
