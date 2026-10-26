"""hobench — literature-based benchmark for hierarchical orthology inference.

Modules
-------
speciestree   rooted species tree with named levels and level arithmetic
reference     reference-family format (family.yaml, members.tsv, reference_tree.nhx), validation
predictions   normalised predicted HOGs + parsers (OrthoXML, OrthoFinder N*.tsv, simple TSV)
metrics       merge-level based metrics (core_spec.md 1–6) and aggregation
controls      positive / negative / dose-response control predictions
evaluate      command-line evaluation (``python -m hobench.evaluate``)
"""
__version__ = "0.1.0"
