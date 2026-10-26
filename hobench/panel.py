"""Pinned panel data (``panels/<P>/species.tsv``) and member checks against it (C4).

``species.tsv`` columns: ``code, upid, taxid, mnemonic, fasta_path, sha256`` (optional
``strict_mapping`` = 1 for species needing the stricter identity threshold). ``fasta_path``
may be absolute, relative to the repository root (the parent of ``panels/``) or relative to
the panel directory — the first existing file is used; none existing is an error.

:func:`validate_members` turns every violation into an ERROR:

* a member ``gene_id`` absent from its species' pinned FASTA (headers are normalised with the
  ``auto`` ID normaliser, i.e. ``sp|ACC|NAME`` → ``ACC``);
* a sequence-based mapping (``mapping_method`` matching :data:`SEQUENCE_METHOD_RE`) with
  identity below 98% (99% for strict species — XENLA, salmonids, or ``strict_mapping`` — and
  for families tagged TAND or with ``strict_mapping: true``) or coverage of either sequence
  below 80%. ``mapping_identity`` is a percentage (0–100); ``mapping_coverage`` a fraction
  (0–1, the minimum over both sequences), or separate ``mapping_qcov``/``mapping_tcov``;
* an unmapped/ambiguous/partial mapping method, or an empty one.

Any other ``mapping_method`` (UniProt ID-mapping source names, ``exact_sha256``,
``UniParc_fallback``, …) is an exact ID/sequence mapping and is exempt from the thresholds.
"""
from __future__ import annotations

import hashlib
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import pandas as pd

from .ids import uniprot
from .reference import ReferenceFamily, ValidationReport

REQUIRED_COLUMNS = ["code", "upid", "taxid", "mnemonic", "fasta_path", "sha256"]
SEQUENCE_METHOD_RE = re.compile(r"^(diamond|mmseqs|blast|blastp|sequence|seq[_-]|ssearch|sw[_-]?align)", re.I)
NOT_MAPPED = {"", "unmapped", "ambiguous", "partial", "exact_ambiguous", "sequence_not_fetched"}
STRICT_SPECIES = {"XENLA", "SALSA", "ONCMY", "ONCKI", "ONCTS", "ONCNE", "SALTR", "SALAL", "HUCHU",
                  "THYTH", "CORCL", "CORLA"}
STRICT_FAMILY_TAGS = {"TAND"}
MIN_IDENTITY, MIN_IDENTITY_STRICT, MIN_COVERAGE = 98.0, 99.0, 0.80


def sha256_file(path: str | Path, chunk: int = 1 << 20) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(chunk), b""):
            h.update(b)
    return h.hexdigest()


@dataclass
class PanelData:
    name: str
    dir: Path
    species: pd.DataFrame
    _ids: dict[str, set[str]] = field(default_factory=dict, repr=False)

    def fasta_path(self, code: str) -> Path:
        row = self.species[self.species["code"] == code]
        if row.empty:
            raise ValueError(f"panel {self.name}: species {code!r} not in species.tsv")
        p = Path(row.iloc[0]["fasta_path"])
        cands = [p] if p.is_absolute() else [self.dir.parents[1] / p, self.dir / p]
        for c in cands:
            if c.exists():
                return c
        raise FileNotFoundError(f"panel {self.name}: FASTA for {code} not found (tried {[str(c) for c in cands]})")

    def ids(self, code: str) -> set[str]:
        if code not in self._ids:
            out = set()
            with open(self.fasta_path(code)) as fh:
                for line in fh:
                    if line.startswith(">"):
                        out.add(uniprot(line[1:]))
            self._ids[code] = out
        return self._ids[code]

    def verify_checksums(self) -> list[str]:
        errs = []
        for _, row in self.species.iterrows():
            try:
                got = sha256_file(self.fasta_path(row["code"]))
            except FileNotFoundError as e:
                errs.append(str(e))
                continue
            if got != str(row["sha256"]).strip().lower():
                errs.append(f"panel {self.name}: sha256 mismatch for {row['code']} ({got} != {row['sha256']})")
        return errs

    def is_strict_species(self, code: str) -> bool:
        if code in STRICT_SPECIES:
            return True
        if "strict_mapping" in self.species.columns:
            row = self.species[self.species["code"] == code]
            return not row.empty and str(row.iloc[0]["strict_mapping"]).strip() in ("1", "true", "True", "yes")
        return False


def load_panel(panel_dir: str | Path) -> PanelData:
    d = Path(panel_dir)
    f = d / "species.tsv"
    if not f.exists():
        raise FileNotFoundError(f"{f} not found")
    df = pd.read_csv(f, sep="\t", dtype=str, keep_default_na=False)
    missing = [c for c in REQUIRED_COLUMNS if c not in df.columns]
    if missing:
        raise ValueError(f"{f}: missing columns {missing}")
    return PanelData(name=d.name, dir=d.resolve(), species=df)


def _float(x) -> Optional[float]:
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def validate_members(fam: ReferenceFamily, panel: PanelData) -> ValidationReport:
    """C4: every member must be in its species' pinned FASTA and meet the mapping thresholds."""
    r = ValidationReport()
    P = f"[{fam.id}]"
    m = fam.members
    strict_family = bool(fam.meta.get("strict_mapping")) or bool(STRICT_FAMILY_TAGS & {t.upper() for t in fam.phenomena})
    codes = set(panel.species["code"])
    for row in m.to_dict("records"):
        g, sp = str(row["gene_id"]), str(row["species"])
        if sp not in codes:
            r.errors.append(f"{P} {g}: species {sp!r} not in panel {panel.name} species.tsv")
            continue
        try:
            ids = panel.ids(sp)
        except FileNotFoundError as e:
            r.errors.append(f"{P} {e}")
            continue
        if g not in ids:
            r.errors.append(f"{P} {g}: not found in the pinned {sp} FASTA of panel {panel.name}")
        method = str(row.get("mapping_method", "") or "").strip()
        if method.lower() in NOT_MAPPED:
            r.errors.append(f"{P} {g}: mapping_method {method!r} is not an accepted mapping")
            continue
        if not SEQUENCE_METHOD_RE.match(method):
            continue  # exact ID / exact sequence mapping
        ident = _float(row.get("mapping_identity"))
        thr = MIN_IDENTITY_STRICT if (strict_family or panel.is_strict_species(sp)) else MIN_IDENTITY
        if ident is None:
            r.errors.append(f"{P} {g}: sequence mapping without mapping_identity")
        elif ident < thr:
            r.errors.append(f"{P} {g}: mapping identity {ident}% < {thr}% (percent scale)")
        covs = []
        if row.get("mapping_coverage") not in (None, ""):
            covs.append(_float(row.get("mapping_coverage")))
        for c in ("mapping_qcov", "mapping_tcov"):
            if row.get(c) not in (None, ""):
                covs.append(_float(row.get(c)))
        if not covs or any(c is None for c in covs):
            r.errors.append(f"{P} {g}: sequence mapping without a numeric mapping_coverage")
        elif any(c > 1.0 for c in covs):
            r.errors.append(f"{P} {g}: mapping_coverage {covs} must be a fraction in [0, 1]")
        elif min(covs) < MIN_COVERAGE:
            r.errors.append(f"{P} {g}: mapping coverage {min(covs):.2f} < {MIN_COVERAGE}")
    return r
