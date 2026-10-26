"""Shared helpers: polite HTTP with retries, hashing, FASTA I/O, logging."""
from __future__ import annotations

import hashlib
import logging
import os
import re
import sys
import threading
import time
from collections import defaultdict
from pathlib import Path
from typing import Iterable, Iterator

import requests

log = logging.getLogger("curation_tools")

USER_AGENT = "curation_tools/0.1 (orthology-benchmark curation; contact: %s)" % os.environ.get(
    "CURATION_EMAIL", "unset"
)

# Minimum seconds between requests to the same host (be polite).
HOST_MIN_INTERVAL = {
    "eutils.ncbi.nlm.nih.gov": 0.4,   # NCBI: 3 req/s without API key
    "www.ncbi.nlm.nih.gov": 0.4,
    "pmc.ncbi.nlm.nih.gov": 0.5,
    "rest.uniprot.org": 0.1,
    "www.ebi.ac.uk": 0.1,
}
_last_hit: dict[str, float] = defaultdict(float)
_lock = threading.Lock()
_session: requests.Session | None = None


class HttpError(RuntimeError):
    """Raised when a request ultimately fails (after retries)."""

    def __init__(self, url: str, status: int | None, msg: str = ""):
        super().__init__(f"{status} {url} {msg}".strip())
        self.url, self.status, self.msg = url, status, msg


def setup_logging(verbose: bool = False) -> None:
    logging.basicConfig(
        level=logging.DEBUG if verbose else logging.INFO,
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
        stream=sys.stderr,
    )


def session() -> requests.Session:
    global _session
    if _session is None:
        _session = requests.Session()
        _session.headers.update({"User-Agent": USER_AGENT})
    return _session


def _throttle(url: str) -> None:
    host = re.sub(r"^https?://([^/]+).*$", r"\1", url)
    gap = HOST_MIN_INTERVAL.get(host, 0.0)
    if gap:
        with _lock:
            wait = _last_hit[host] + gap - time.time()
            if wait > 0:
                time.sleep(wait)
            _last_hit[host] = time.time()


def http(method: str, url: str, *, params=None, data=None, headers=None, stream=False,
         timeout: float = 120, retries: int = 5, ok_status: Iterable[int] = (200,),
         backoff: float = 2.0) -> requests.Response:
    """HTTP request with throttling, retry on 429/5xx/network errors, Retry-After support.

    Non-retryable statuses not in ``ok_status`` raise HttpError immediately.
    """
    ok = set(ok_status)
    last: Exception | None = None
    for attempt in range(retries + 1):
        _throttle(url)
        try:
            r = session().request(method, url, params=params, data=data, headers=headers,
                                  stream=stream, timeout=timeout)
        except (requests.ConnectionError, requests.Timeout, requests.exceptions.ChunkedEncodingError) as e:
            last = e
            wait = backoff ** attempt
            log.warning("network error on %s (%s); retry %d/%d in %.0fs", url, e.__class__.__name__,
                        attempt + 1, retries, wait)
            time.sleep(wait)
            continue
        if r.status_code in ok:
            return r
        if r.status_code in (429, 500, 502, 503, 504):
            ra = r.headers.get("Retry-After")
            try:
                wait = float(ra) if ra else backoff ** attempt
            except ValueError:
                wait = backoff ** attempt
            wait = min(wait, 120)
            last = HttpError(url, r.status_code)
            # a persistent 500 (e.g. Europe PMC for publisher-restricted XML) is not worth 5 retries
            if r.status_code == 500 and attempt >= 1:
                r.close()
                raise HttpError(url, 500, "server error (persistent)")
            log.warning("HTTP %s on %s; retry %d/%d in %.0fs", r.status_code, url, attempt + 1, retries, wait)
            r.close()
            time.sleep(wait)
            continue
        status = r.status_code
        body = " ".join(re.sub(r"<[^>]+>", " ", r.text[:600]).split())[:120] if not stream else ""
        r.close()
        raise HttpError(url, status, body)
    raise HttpError(url, getattr(last, "status", None), f"gave up after {retries} retries: {last}")


def get(url: str, **kw) -> requests.Response:
    return http("GET", url, **kw)


def post(url: str, **kw) -> requests.Response:
    return http("POST", url, **kw)


def sha256_bytes(b: bytes) -> str:
    return hashlib.sha256(b).hexdigest()


def sha256_file(path: str | Path, chunk: int = 1 << 20) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while True:
            b = f.read(chunk)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def md5_str(s: str) -> str:
    return hashlib.md5(s.encode()).hexdigest()


def read_fasta(path: str | Path) -> Iterator[tuple[str, str]]:
    """Yield (header_without_>, sequence) from a FASTA file (plain or .gz)."""
    import gzip

    op = gzip.open if str(path).endswith(".gz") else open
    hdr, seq = None, []
    with op(path, "rt") as f:
        for line in f:
            line = line.rstrip("\r\n")
            if line.startswith(">"):
                if hdr is not None:
                    yield hdr, "".join(seq)
                hdr, seq = line[1:], []
            elif line.strip():
                seq.append(line.strip())
        if hdr is not None:
            yield hdr, "".join(seq)


def write_fasta(records: Iterable[tuple[str, str]], path: str | Path, width: int = 60) -> int:
    n = 0
    with open(path, "w") as f:
        for h, s in records:
            f.write(f">{h}\n")
            for i in range(0, len(s), width):
                f.write(s[i:i + width] + "\n")
            n += 1
    return n


def norm_seq(s: str) -> str:
    """Normalise a protein sequence for exact-match hashing (upper-case, strip gaps/stop/space)."""
    return re.sub(r"[\s\-*.]", "", s).upper()


def seq_sha256(s: str) -> str:
    return hashlib.sha256(norm_seq(s).encode()).hexdigest()


def safe_name(s: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]+", "_", s).strip("_")[:120]


def now_iso() -> str:
    return time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())


def have_exe(name: str) -> str | None:
    """Locate an executable on PATH or next to the running python (conda env bin)."""
    import shutil

    p = shutil.which(name)
    if p:
        return p
    cand = Path(sys.executable).parent / name
    return str(cand) if cand.exists() else None


def pdf_to_text(path: str | Path) -> tuple[str, str]:
    """Extract text from a PDF. Returns (text, method). Uses pdftotext if installed, else pypdf.

    Raises RuntimeError if neither backend is available.
    """
    import subprocess

    exe = have_exe("pdftotext")
    if exe:
        r = subprocess.run([exe, "-layout", str(path), "-"], capture_output=True, text=True, timeout=300)
        if r.returncode == 0:
            return r.stdout, "pdftotext"
        log.warning("pdftotext failed on %s: %s", path, r.stderr[:200])
    try:
        from pypdf import PdfReader
    except ImportError as e:  # pragma: no cover
        raise RuntimeError("no PDF backend: install poppler-utils (pdftotext) or `pip install pypdf`") from e
    rd = PdfReader(str(path))
    return "\n".join((pg.extract_text() or "") + "\n\f" for pg in rd.pages), "pypdf"
