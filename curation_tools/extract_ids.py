"""Scan paper/supplement files for sequence identifiers and report them with context.

CLI:
  python -m curation_tools.extract_ids literature/papers/PMID16128622 [more files/dirs] \
        --out ids.tsv [--unique] [--types uniprot,refseq_protein] [--min-confidence medium]

Supported inputs: text-like (txt, tsv, csv, fasta, newick, json, xml, html, md), XLSX/XLSM (openpyxl),
XLS (xlrd), DOCX (python-docx), legacy DOC (crude string extraction, flagged), PDF (pdftotext, else
pypdf), and .zip / .gz archives (recursed).

Output TSV columns:
  id  id_type  subtype  confidence  file  location  column  col_header  context
where context is the +-1 text lines or, for tables, the whole table row (cells joined by ' | ').
Unreadable files are reported on stderr and in the summary line; nothing is skipped silently.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import io
import re
import struct
import sys
import tempfile
import zipfile
from collections import OrderedDict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator

from . import common
from .common import log

# ----------------------------------------------------------------------------- patterns
# (id_type, subtype, compiled regex, base confidence). Order matters: earlier patterns claim spans first.
_I = re.IGNORECASE
PATTERNS: list[tuple[str, str, re.Pattern, str]] = [
    ("refseq_protein", "", re.compile(r"(?<![A-Za-z0-9_])(?:NP|XP|YP|WP|AP)_\d{6,9}(?:\.\d+)?(?![A-Za-z0-9])"), "high"),
    ("refseq_nucleotide", "", re.compile(r"(?<![A-Za-z0-9_])(?:NM|XM|NR|XR|NC|NG|NW|NT|NZ)_\d{6,12}(?:\.\d+)?(?![A-Za-z0-9])"), "high"),
    ("ensembl_protein", "", re.compile(r"(?<![A-Za-z0-9])ENS[A-Z]{0,8}P\d{11}(?:\.\d+)?(?![A-Za-z0-9])"), "high"),
    ("ensembl_transcript", "", re.compile(r"(?<![A-Za-z0-9])ENS[A-Z]{0,8}T\d{11}(?:\.\d+)?(?![A-Za-z0-9])"), "high"),
    ("ensembl_gene", "", re.compile(r"(?<![A-Za-z0-9])ENS[A-Z]{0,8}G\d{11}(?:\.\d+)?(?![A-Za-z0-9])"), "high"),
    ("flybase_protein", "", re.compile(r"(?<![A-Za-z0-9])FBpp\d{7}(?![A-Za-z0-9])"), "high"),
    ("flybase_transcript", "", re.compile(r"(?<![A-Za-z0-9])FBtr\d{7}(?![A-Za-z0-9])"), "high"),
    ("flybase_gene", "", re.compile(r"(?<![A-Za-z0-9])FBgn\d{7}(?![A-Za-z0-9])"), "high"),
    ("flybase_annotation", "CG", re.compile(r"(?<![A-Za-z0-9])CG\d{3,5}(?:-P[A-Z])?(?![A-Za-z0-9])"), "medium"),
    ("wormbase_gene", "", re.compile(r"(?<![A-Za-z0-9])WBGene\d{8}(?![A-Za-z0-9])"), "high"),
    ("wormbase_protein", "", re.compile(r"(?<![A-Za-z0-9])WP:CE\d{5}(?![A-Za-z0-9])"), "high"),
    ("sgd_systematic", "", re.compile(r"(?<![A-Za-z0-9])Y[A-P][LR]\d{3}[WC](?:-[A-Z])?(?![A-Za-z0-9])"), "high"),
    ("sgd_id", "", re.compile(r"(?<![A-Za-z0-9])S0\d{8}(?![A-Za-z0-9])"), "medium"),
    ("arabidopsis_agi", "", re.compile(r"(?<![A-Za-z0-9])AT[1-5CM]G\d{5}(?:\.\d{1,2})?(?![A-Za-z0-9])", _I), "high"),
    ("wormbase_seqname", "", re.compile(
        r"(?<![A-Za-z0-9.])(?:[A-Z]{1,2}\d{1,2}[A-Z]\d{1,2}[A-Z]?|ZK\d{1,4}|[A-Z]\d{4}|[A-Z]{2}\d{1,3})\.\d{1,2}[a-z]?(?![A-Za-z0-9.])"), "low"),
    # --- "other" families: species/consortium specific identifiers (report as 'other' with subtype)
    ("other", "jgi", re.compile(r"(?<![A-Za-z0-9])jgi\|[A-Za-z0-9_]+\|\d+(?![A-Za-z0-9])", _I), "medium"),
    ("other", "gi", re.compile(r"(?<![A-Za-z0-9])gi\|\d{4,}(?![A-Za-z0-9])"), "medium"),
    ("other", "phytozome", re.compile(r"(?<![A-Za-z0-9])(?:Potri|Glyma|Sobic|Sevir|Pavir|Bradi\d?|Zm\d{5}|Medtr|Solyc|Phvul|Cre\d\d|Mp\d?|Pp\d)[A-Za-z0-9]*[._]\d{1,3}[Gg]?\d{4,9}(?:\.\d+)*(?![A-Za-z0-9])"), "medium"),
    ("other", "maize", re.compile(r"(?<![A-Za-z0-9])(?:Zm\d{5}d\d{6}(?:_[PT]\d+)?|GRMZM\d[GT]\d{6}(?:_[PT]\d+)?)(?![A-Za-z0-9])"), "medium"),
    ("other", "rice", re.compile(r"(?<![A-Za-z0-9])(?:LOC_)?Os\d{2}g\d{5,7}(?:\.\d+)?(?![A-Za-z0-9])"), "medium"),
    ("other", "pombase", re.compile(r"(?<![A-Za-z0-9])SP(?:A|B|C|CC|CP|AC|AP|BC|BP)\d{1,4}[A-Za-z0-9]*\.\d{2}c?(?![A-Za-z0-9])"), "medium"),
    ("other", "plasmodb", re.compile(r"(?<![A-Za-z0-9])(?:PF3D7|PVP01|PBANKA|PCHAS|PY17X|PKNH)_[0-9A-Z]{6,10}(?![A-Za-z0-9])"), "medium"),
    ("other", "toxodb", re.compile(r"(?<![A-Za-z0-9])TGME49_\d{5,7}(?![A-Za-z0-9])"), "medium"),
    ("other", "tritrypdb", re.compile(r"(?<![A-Za-z0-9])(?:Tb\d{2,3}\.\w+\.\d+|LmjF\.\d+\.\d+|TcCLB\.\d+\.\d+)(?![A-Za-z0-9])"), "medium"),
    ("other", "dictybase", re.compile(r"(?<![A-Za-z0-9])DDB_G\d{7}(?![A-Za-z0-9])"), "high"),
    ("other", "vectorbase", re.compile(r"(?<![A-Za-z0-9])(?:AGAP|AAEL|CPIJ|AALF|ASTE|AARA)\d{5,6}(?:-P[A-Z])?(?![A-Za-z0-9])"), "medium"),
    ("other", "fungal", re.compile(r"(?<![A-Za-z0-9])(?:MGG|NCU|CNAG|CNB[A-Z]\d|CAGL0|CAL\d|orf19|ANID|AFUA|FGSG|UM|BC1G|SS1G|VDAG)[_.]?[0-9A-Za-z]{4,10}(?![A-Za-z0-9])"), "low"),
    ("other", "ncbi_gene_loc", re.compile(r"(?<![A-Za-z0-9])LOC\d{5,9}(?![A-Za-z0-9])"), "medium"),
    ("other", "ncbi_geneid", re.compile(r"(?<![A-Za-z0-9])GeneID[:=]\s?\d{1,9}", _I), "medium"),
    ("other", "wgs_protein", re.compile(r"(?<![A-Za-z0-9])[A-Z]{4}\d{2}\d{6,}(?:\.\d+)?(?![A-Za-z0-9])"), "low"),
    ("uniprot_entry_name", "", re.compile(
        r"(?<![A-Za-z0-9_])[A-Z0-9]{1,10}_(?:HUMAN|MOUSE|RAT|YEAST|ARATH|DROME|CAEEL|DANRE|SCHPO|ECOLI|BOVIN|CHICK|XENLA|XENTR|PIG|MACMU|"
        r"DICDI|ORYSJ|MAIZE|PLAF7|CANAL|NEUCR|ASPNG|BACSU|MYCTU|HELPY|SALTY|PSEAE|ANOGA|APIME|BRAFL|NEMVE|STRPU|TRIAD|MONBE|"
        r"CIOIN|CIOSA|LEIMA|TRYB2|TOXGO|SOLLC|SOYBN|VITVI|POPTR|PHYPA|CHLRE|ECHMS)(?![A-Za-z0-9_])"), "medium"),
    # generic accession formats last (lowest specificity)
    ("uniprot", "", re.compile(
        r"(?<![A-Za-z0-9_.\-])(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})(?:-\d{1,3})?(?![A-Za-z0-9_])"), "low"),
    ("genbank_protein", "", re.compile(r"(?<![A-Za-z0-9_.\-])[A-Z]{3}\d{5}(?:\d{2})?\.\d{1,2}(?![A-Za-z0-9])"), "medium"),
    ("genbank_protein", "no_version", re.compile(r"(?<![A-Za-z0-9_.\-])[A-Z]{3}\d{5}(?:\d{2})?(?![A-Za-z0-9.])"), "low"),
]

# 3-letter prefixes that look like GenBank protein accessions but are something else
NOT_GENBANK = {"PMC", "SRR", "SRX", "SRS", "SRP", "ERR", "ERX", "ERS", "DRR", "DRX", "DRS", "PRJ", "GSE", "GSM",
               "GPL", "GCA", "GCF", "ISO", "CAS", "EFO", "EC0", "RFC", "MIM", "ICD", "NCT", "ISRCTN", "WGS",
               "MPI", "HGNC", "PDB", "DOI", "ISS", "ISB", "MAP", "SDS", "MSG", "TAB", "FIG", "SEQ", "ORF", "NGS",
               "CDS", "UTR", "BAC", "YAC", "PAC", "RNA", "DNA", "EST", "SNP", "GWA", "TIF", "PNG", "JPG", "PDF"}
CONF_ORDER = {"low": 0, "medium": 1, "high": 2}
KEYWORDS = re.compile(r"uniprot|swiss-?prot|trembl|accession|genbank|ncbi|refseq|protein id|ensembl|entry", _I)


@dataclass
class Hit:
    id: str
    id_type: str
    subtype: str
    confidence: str
    span: tuple[int, int]


def find_ids(text: str, context_hint: str = "") -> list[Hit]:
    """Find identifiers in one string. Specific patterns claim character spans first; generic ones
    (UniProt / GenBank accession shapes) only match what is left."""
    claimed: list[tuple[int, int]] = []
    hits: list[Hit] = []

    def free(a: int, b: int) -> bool:
        return all(b <= x or a >= y for x, y in claimed)

    for id_type, subtype, rx, conf in PATTERNS:
        for m in rx.finditer(text):
            a, b = m.span()
            if not free(a, b):
                continue
            s = m.group(0)
            if id_type == "genbank_protein" and s[:3] in NOT_GENBANK:
                continue
            c = conf
            if id_type in ("uniprot", "genbank_protein") and KEYWORDS.search(context_hint or text):
                c = "medium" if c == "low" else c
            if id_type == "uniprot":
                if re.fullmatch(r"[A-Z]\d[A-Z0-9]{3}\d(?:-\d+)?", s) and not re.search(r"[A-Z]", s[1:]):
                    pass  # e.g. P12345: plausible
                if re.fullmatch(r"[A-Z]\d{5}", s) and c == "low" and s[0] not in "OPQ":
                    continue  # letter+5 digits not starting with O/P/Q is almost always something else
            if id_type == "arabidopsis_agi":
                s = s.upper()
            if id_type == "other" and subtype == "ncbi_geneid":
                s = re.sub(r"\s", "", s)
            claimed.append((a, b))
            hits.append(Hit(s, id_type, subtype, c, (a, b)))
    hits.sort(key=lambda h: h.span)
    return hits


# ----------------------------------------------------------------------------- readers
@dataclass
class Record:
    """One unit of scanning: a text line or a table row."""
    location: str
    text: str                      # string to regex-scan (cells joined)
    cells: list[str] | None = None  # for table rows
    header: list[str] | None = None
    prev: str = ""
    next: str = ""


def _clip(s: str, a: int | None = None, b: int | None = None, width: int = 300) -> str:
    s = s.replace("\t", " ").replace("\r", " ").replace("\n", " ")
    if len(s) <= width:
        return s
    if a is None:
        return s[:width] + "..."
    lo = max(0, a - width // 2)
    return ("..." if lo else "") + s[lo:lo + width] + "..."


TEXT_EXT = {".txt", ".tsv", ".csv", ".tab", ".fa", ".fasta", ".faa", ".pep", ".aa", ".nwk", ".tre", ".newick",
            ".json", ".xml", ".html", ".htm", ".md", ".fas", ".phy", ".list", ".ids", ".gff", ".gff3", ".gtf",
            ".nex", ".out", ".log", ".dat", ".text", ""}
FASTA_EXT = {".fa", ".fasta", ".faa", ".pep", ".aa", ".fas"}


def _text_records(lines: list[str], prefix: str = "line", fasta: bool = False) -> Iterator[Record]:
    for i, ln in enumerate(lines):
        if fasta and not ln.startswith(">"):
            continue
        yield Record(f"{prefix} {i + 1}", ln,
                     prev=lines[i - 1] if i > 0 else "", next=lines[i + 1] if i + 1 < len(lines) else "")


def _table_records(rows: Iterator[list[str]], prefix: str) -> Iterator[Record]:
    header = None
    for i, row in enumerate(rows):
        row = ["" if c is None else str(c) for c in row]
        if not any(x.strip() for x in row):
            continue
        if header is None:
            header = row
        yield Record(f"{prefix} row {i + 1}", " \t ".join(row), cells=row, header=header)


def _delimited(path: Path, delim: str) -> Iterator[Record]:
    csv.field_size_limit(1 << 27)
    with open(path, newline="", errors="replace") as f:
        yield from _table_records(csv.reader(f, delimiter=delim), "")


def _xlsx(path: Path) -> Iterator[Record]:
    import openpyxl
    wb = openpyxl.load_workbook(path, read_only=True, data_only=True)
    for ws in wb.worksheets:
        yield from _table_records((list(r) for r in ws.iter_rows(values_only=True)), f"sheet '{ws.title}'")


def _xls(path: Path) -> Iterator[Record]:
    import xlrd
    wb = xlrd.open_workbook(path)
    for ws in wb.sheets():
        yield from _table_records(([c for c in ws.row_values(r)] for r in range(ws.nrows)), f"sheet '{ws.name}'")


def _docx(path: Path) -> Iterator[Record]:
    import docx
    d = docx.Document(str(path))
    paras = [p.text for p in d.paragraphs]
    yield from _text_records(paras, "para")
    for ti, t in enumerate(d.tables, 1):
        yield from _table_records(([c.text for c in r.cells] for r in t.rows), f"table {ti}")


def _doc_legacy(path: Path) -> Iterator[Record]:
    """Very crude: pull printable ASCII and UTF-16LE runs out of the OLE binary. Flagged in output."""
    data = path.read_bytes()
    runs = re.findall(rb"[\x20-\x7e\t]{4,}", data)
    u16 = re.findall(rb"(?:[\x20-\x7e]\x00){4,}", data)
    lines = [r.decode("ascii", "replace") for r in runs] + [r.decode("utf-16le", "replace") for r in u16]
    log.warning("%s: legacy .doc parsed with crude string extraction; convert to .docx/.txt for reliable results", path)
    yield from _text_records(lines, "doc-string")


def _pdf(path: Path) -> Iterator[Record]:
    text, method = common.pdf_to_text(path)
    for pi, page in enumerate(text.split("\f"), 1):
        lines = page.split("\n")
        yield from _text_records(lines, f"page {pi} line")


SKIP_EXT = {".jpg", ".jpeg", ".gif", ".png", ".tif", ".tiff", ".bmp", ".svg", ".eps", ".ps", ".mp4", ".mov",
            ".avi", ".pptx", ".ppt", ".dmnd", ".h5", ".pdb", ".cif", ".ico", ".mp3"}


def iter_records(path: Path) -> Iterator[Record]:
    ext = path.suffix.lower()
    if ext in (".xlsx", ".xlsm"):
        yield from _xlsx(path)
    elif ext == ".xls":
        yield from _xls(path)
    elif ext == ".docx":
        yield from _docx(path)
    elif ext == ".doc":
        yield from _doc_legacy(path)
    elif ext == ".pdf":
        yield from _pdf(path)
    elif ext in (".tsv", ".tab"):
        yield from _delimited(path, "\t")
    elif ext == ".csv":
        yield from _delimited(path, ",")
    else:
        if ext not in TEXT_EXT:   # unknown extension (.nhx, .treefile, .aln, ...): scan if it looks like text
            with open(path, "rb") as f:
                head = f.read(4096)
            if b"\x00" in head:
                raise ValueError(f"unsupported binary file type {ext!r}")
        lines = path.read_text(errors="replace").splitlines()
        yield from _text_records(lines, fasta=ext in FASTA_EXT)


def _walk(inputs: list[Path]) -> Iterator[tuple[Path, str]]:
    """Yield (real_path, display_name), expanding dirs, .zip and .gz."""
    tmp_dirs: list[tempfile.TemporaryDirectory] = []
    for p in inputs:
        if p.is_dir():
            for q in sorted(p.rglob("*")):
                if q.is_file() and q.name != "meta.json" and not q.name.endswith(".dmnd") \
                        and not _redundant(q):
                    yield from _walk_file(q, str(q), tmp_dirs)
        elif p.is_file():
            yield from _walk_file(p, str(p), tmp_dirs)
        else:
            log.error("input not found: %s", p)
            yield p, str(p)
    for t in tmp_dirs:
        t.cleanup()


def _redundant(q: Path) -> bool:
    """Inside a fetch_paper folder, skip copies of content that is scanned in another form."""
    d = q.parent
    if q.name == "supplementary.zip" and (d / "supp").is_dir():
        return True                                   # members already extracted to raw/supp/
    if q.name in ("fulltext.xml", "article.html", "fulltext.pdf") and \
            ((d / "fulltext.txt").exists() or (d.parent / "fulltext.txt").exists()):
        return True
    return False


def _walk_file(p: Path, disp: str, tmp_dirs: list):
    ext = p.suffix.lower()
    if ext == ".zip":
        td = tempfile.TemporaryDirectory()
        tmp_dirs.append(td)
        try:
            with zipfile.ZipFile(p) as z:
                for i in z.infolist():
                    if i.is_dir():
                        continue
                    out = Path(td.name) / f"{len(list(Path(td.name).iterdir()))}_{common.safe_name(i.filename)}"
                    out.write_bytes(z.read(i))
                    yield from _walk_file(out, f"{disp}!{i.filename}", tmp_dirs)
        except zipfile.BadZipFile as e:
            log.error("bad zip %s: %s", p, e)
            yield p, disp
    elif ext == ".gz":
        td = tempfile.TemporaryDirectory()
        tmp_dirs.append(td)
        out = Path(td.name) / p.stem
        with gzip.open(p, "rb") as f:
            out.write_bytes(f.read())
        yield from _walk_file(out, f"{disp}!{p.stem}", tmp_dirs)
    else:
        yield p, disp


# ----------------------------------------------------------------------------- scan
COLS = ["id", "id_type", "subtype", "confidence", "file", "location", "column", "col_header", "context"]


def _col_letter(i: int) -> str:
    s = ""
    i += 1
    while i:
        i, r = divmod(i - 1, 26)
        s = chr(65 + r) + s
    return s


def scan_file(path: Path, disp: str, min_conf: str = "low", types: set[str] | None = None) -> Iterator[dict]:
    for rec in iter_records(path):
        if rec.cells is not None:
            # scan per cell so the column is known; context = whole row
            row_ctx = _clip(" | ".join(rec.cells), width=500)
            for ci, cell in enumerate(rec.cells):
                if not cell:
                    continue
                for h in find_ids(cell, " ".join(rec.header or [])):
                    yield _row(h, disp, rec.location, _col_letter(ci), (rec.header[ci] if rec.header and ci < len(rec.header) else ""),
                               row_ctx, min_conf, types)
        else:
            if not rec.text:
                continue
            for h in find_ids(rec.text):
                a, b = h.span
                ctx = " // ".join(x for x in (_clip(rec.prev), _clip(rec.text, a, b), _clip(rec.next)) if x)
                yield _row(h, disp, rec.location, "", "", ctx, min_conf, types)


def _row(h: Hit, disp, location, col, col_header, ctx, min_conf, types):
    if CONF_ORDER[h.confidence] < CONF_ORDER[min_conf] or (types and h.id_type not in types):
        return None
    return dict(id=h.id, id_type=h.id_type, subtype=h.subtype, confidence=h.confidence, file=disp,
                location=location.strip(), column=col, col_header=_clip(col_header, width=80), context=ctx)


SKIPPED: list[str] = []


def scan(inputs: list[str | Path], min_conf: str = "low", types: set[str] | None = None) -> tuple[list[dict], list[str]]:
    rows, errors = [], []
    for real, disp in _walk([Path(x) for x in inputs]):
        if not real.exists():
            errors.append(f"{disp}: not found")
            continue
        if real.suffix.lower() in SKIP_EXT:
            log.debug("skipping non-text file %s", disp)
            SKIPPED.append(disp)
            continue
        try:
            n0 = len(rows)
            for r in scan_file(real, disp, min_conf, types):
                if r:
                    rows.append(r)
            log.info("%s: %d identifier hits", disp, len(rows) - n0)
        except Exception as e:  # noqa: BLE001  - report and continue; never silent
            msg = f"{disp}: could not be scanned ({e.__class__.__name__}: {e})"
            log.error(msg)
            errors.append(msg)
    return rows, errors


def unique(rows: list[dict]) -> list[dict]:
    seen: "OrderedDict[tuple, dict]" = OrderedDict()
    for r in rows:
        k = (r["id"], r["id_type"])
        if k in seen:
            seen[k]["n_occurrences"] += 1
        else:
            seen[k] = dict(r, n_occurrences=1)
    return list(seen.values())


def write_tsv(rows: list[dict], out) -> None:
    cols = COLS + (["n_occurrences"] if rows and "n_occurrences" in rows[0] else [])
    w = csv.writer(out, delimiter="\t", lineterminator="\n")
    w.writerow(cols)
    for r in rows:
        w.writerow([str(r.get(c, "")).replace("\t", " ").replace("\n", " ") for c in cols])


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("inputs", nargs="+", help="files or directories")
    ap.add_argument("--out", "-o", help="output TSV (default: stdout)")
    ap.add_argument("--unique", action="store_true", help="one row per distinct id (first context, with count)")
    ap.add_argument("--types", help="comma-separated id_types to keep")
    ap.add_argument("--min-confidence", choices=list(CONF_ORDER), default="low")
    ap.add_argument("-v", "--verbose", action="store_true")
    args = ap.parse_args(argv)
    common.setup_logging(args.verbose)
    types = set(args.types.split(",")) if args.types else None
    SKIPPED.clear()
    rows, errors = scan(args.inputs, args.min_confidence, types)
    if args.unique:
        rows = unique(rows)
    out = open(args.out, "w") if args.out else sys.stdout
    write_tsv(rows, out)
    if args.out:
        out.close()
    from collections import Counter
    c = Counter(r["id_type"] for r in rows)
    log.info("%d rows; by type: %s; %d file(s) with errors; %d non-text file(s) skipped (images etc.)",
             len(rows), dict(c), len(errors), len(SKIPPED))
    return 3 if errors else 0


if __name__ == "__main__":
    sys.exit(main())
