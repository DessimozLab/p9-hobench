"""Convert JATS XML (Europe PMC / PMC OA) and PMC HTML article pages to plain text.

Output keeps: title, authors, abstract, section headers ('#', '##', ...), paragraphs,
figure captions ([FIGURE id] ...), tables as TSV blocks ([TABLE id] ... [/TABLE id]),
supplementary-material captions, footnotes and the reference list.
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field

from lxml import etree, html as lhtml

XLINK = "{http://www.w3.org/1999/xlink}href"


@dataclass
class Converted:
    text: str = ""
    tables: dict = field(default_factory=dict)       # id -> TSV string
    n_figures: int = 0
    n_tables: int = 0
    supp_hrefs: list = field(default_factory=list)   # hrefs mentioned in the XML
    warnings: list = field(default_factory=list)


def _ln(el) -> str:
    t = el.tag
    if not isinstance(t, str):
        return ""
    return t.rsplit("}", 1)[-1] if "}" in t else t


def _ws(s: str) -> str:
    return re.sub(r"\s+", " ", s).strip()


# ----------------------------------------------------------------------------- tables
def _grid_from_rows(rows) -> str:
    """rows: list of list of (text, colspan, rowspan). Expand spans into a rectangular TSV."""
    grid: list[list] = []
    pending: dict[tuple[int, int], str] = {}
    for ri, row in enumerate(rows):
        out: list[str] = []
        ci = 0
        cells = list(row)
        while cells or (ri, ci) in pending:
            if (ri, ci) in pending:
                out.append(pending.pop((ri, ci)))
                ci += 1
                continue
            txt, cs, rs = cells.pop(0)
            for k in range(cs):
                out.append(txt if k == 0 else "")
                for r2 in range(1, rs):
                    pending[(ri + r2, ci + k)] = txt if k == 0 else ""
            ci += cs
        grid.append(out)
    return "\n".join("\t".join(c.replace("\t", " ").replace("\n", " ") for c in r) for r in grid)


def _span(el, name: str) -> int:
    try:
        return max(1, int(el.get(name, "1")))
    except ValueError:
        return 1


def table_to_tsv(table_el, inline) -> str:
    rows = []
    for tr in table_el.iter():
        if _ln(tr) != "tr":
            continue
        row = []
        for c in tr:
            if _ln(c) in ("td", "th"):
                row.append((_ws(inline(c)), _span(c, "colspan"), _span(c, "rowspan")))
        rows.append(row)
    return _grid_from_rows(rows)


# ----------------------------------------------------------------------------- JATS
BLOCK_CHILDREN = {"fig", "fig-group", "table-wrap", "table-wrap-group", "list", "disp-formula",
                  "boxed-text", "disp-quote", "supplementary-material", "sec", "def-list", "preformat"}


class _Jats:
    def __init__(self, root):
        self.root = root
        self.out: list[str] = []
        self.c = Converted()

    # inline text -------------------------------------------------------
    def inline(self, el, deferred: list | None = None) -> str:
        parts: list[str] = []

        def rec(e):
            if e.text:
                parts.append(e.text)
            for ch in e:
                tag = _ln(ch)
                if tag in BLOCK_CHILDREN and deferred is not None:
                    deferred.append(ch)
                elif tag == "sup":
                    t = _ws("".join(ch.itertext()))
                    parts.append(("^" + t) if t and len(t) < 12 else (" " + t))
                elif tag == "sub":
                    t = _ws("".join(ch.itertext()))
                    parts.append(("_" + t) if t else "")
                elif tag == "break":
                    parts.append(" ")
                elif tag == "ext-link":
                    t = _ws("".join(ch.itertext()))
                    href = ch.get(XLINK, "")
                    parts.append(t if (not href or href in t) else f"{t} <{href}>")
                elif tag in ("object-id", "label") and deferred is not None and _ln(e) in ("fig", "table-wrap"):
                    pass
                elif tag in ("graphic", "inline-graphic", "media"):
                    pass
                else:
                    rec(ch)
                if ch.tail:
                    parts.append(ch.tail)

        rec(el)
        return _ws("".join(parts))

    def emit(self, s: str = ""):
        self.out.append(s)

    # blocks ------------------------------------------------------------
    def block(self, el, level: int):
        tag = _ln(el)
        if tag in ("sec", "app", "boxed-text", "ack", "notes", "glossary", "fn-group", "bio", "abstract",
                   "trans-abstract", "app-group", "disp-quote", "fig-group", "table-wrap-group", "def-list"):
            title = el.find("title")
            if title is None:
                lab = el.find("label")
            else:
                lab = None
            head = None
            if title is not None:
                lbl = el.find("label")
                head = (self.inline(lbl) + " " if lbl is not None else "") + self.inline(title)
            elif tag in ("ack", "fn-group", "notes", "boxed-text", "glossary"):
                head = {"ack": "Acknowledgements", "fn-group": "Footnotes", "notes": "Notes",
                        "boxed-text": "Box", "glossary": "Glossary"}[tag]
            if head and tag not in ("fig-group", "table-wrap-group", "disp-quote", "def-list"):
                self.emit()
                self.emit("#" * min(level, 6) + " " + head)
            for ch in el:
                if _ln(ch) in ("title", "label"):
                    continue
                self.block(ch, level + 1 if tag in ("sec", "app", "boxed-text", "notes") else level)
        elif tag == "p":
            deferred: list = []
            t = self.inline(el, deferred)
            if t:
                self.emit(t)
            for d in deferred:
                self.block(d, level)
        elif tag == "title":
            t = self.inline(el)
            if t:
                self.emit("#" * min(level, 6) + " " + t)
        elif tag == "fig":
            self.figure(el)
        elif tag == "table-wrap":
            self.table_wrap(el)
        elif tag in ("list",):
            for i, li in enumerate(el.findall("list-item"), 1):
                deferred = []
                txt = " ".join(self.inline(p, deferred) for p in li.findall("p")) or self.inline(li)
                self.emit(f"- {txt}")
        elif tag == "supplementary-material":
            self.supp(el)
        elif tag in ("disp-formula", "preformat", "fn", "def-item"):
            t = self.inline(el)
            if t:
                self.emit(t)
        elif tag == "ref-list":
            self.emit()
            title = el.find("title")
            self.emit("#" * min(level, 6) + " " + (self.inline(title) if title is not None else "References"))
            for ref in el.iter("ref"):
                rid = ref.get("id", "")
                self.emit(f"[{rid}] " + self.citation(ref))
        elif tag in ("label", "caption", "object-id", "kwd-group", "permissions", "alternatives"):
            return
        else:
            for ch in el:
                self.block(ch, level)

    def citation(self, ref) -> str:
        mc = ref.find("mixed-citation")
        if mc is not None:
            return self.inline(mc)
        ec = ref.find("element-citation")
        if ec is None:
            ec = ref.find("citation")
        if ec is None:
            return self.inline(ref)
        names, rest = [], []
        for ch in ec:
            if _ln(ch) == "person-group":
                for n in ch.iter("name"):
                    names.append(_ws(f"{n.findtext('surname') or ''} {n.findtext('given-names') or ''}"))
                if ch.find(".//collab") is not None:
                    names.append(_ws("".join(ch.find(".//collab").itertext())))
            else:
                t = self.inline(ch)
                if t:
                    rest.append(t)
        return ", ".join(names) + (". " if names else "") + " ".join(rest)

    def caption_text(self, el) -> str:
        cap = el.find("caption")
        if cap is None:
            return ""
        paras = []
        for ch in cap:
            t = self.inline(ch)
            if t:
                paras.append(t)
        return " ".join(paras) if paras else self.inline(cap)

    def figure(self, el):
        self.c.n_figures += 1
        fid = el.get("id", f"fig{self.c.n_figures}")
        lab = el.find("label")
        lab = self.inline(lab) if lab is not None else ""
        self.emit()
        self.emit(f"[FIGURE {fid}] {lab} {self.caption_text(el)}".strip())

    def table_wrap(self, el):
        self.c.n_tables += 1
        tid = el.get("id", f"table{self.c.n_tables}")
        lab = el.find("label")
        lab = self.inline(lab) if lab is not None else ""
        self.emit()
        self.emit(f"[TABLE {tid}] {lab} {self.caption_text(el)}".strip())
        tbl = next((t for t in el.iter() if _ln(t) == "table"), None)
        if tbl is None:
            g = next((t for t in el.iter() if _ln(t) == "graphic"), None)
            self.emit(f"[table is an image: {g.get(XLINK) if g is not None else '?'}; contents not extractable]")
            self.c.warnings.append(f"table {tid} is an image")
        else:
            tsv = table_to_tsv(tbl, self.inline)
            self.c.tables[tid] = tsv
            self.emit(tsv)
        foot = el.find("table-wrap-foot")
        if foot is not None:
            for p in foot.iter():
                if _ln(p) in ("p",):
                    t = self.inline(p)
                    if t:
                        self.emit("[table note] " + t)
        self.emit(f"[/TABLE {tid}]")

    def supp(self, el):
        href = el.get(XLINK, "")
        if href:
            self.c.supp_hrefs.append(href)
        lab = el.find("label")
        lab = self.inline(lab) if lab is not None else ""
        cap = self.caption_text(el)
        self.emit(f"[SUPPLEMENT {el.get('id', '')}] {lab} {cap} (file: {href})".strip())

    # whole article -----------------------------------------------------
    def run(self) -> Converted:
        art = self.root if _ln(self.root) == "article" else next(self.root.iter("article"), self.root)
        front = art.find("front")
        body = art.find("body")
        back = art.find("back")
        if front is not None:
            self.front(front)
        for sub in art.findall("sub-article"):
            pass
        if body is not None:
            self.emit()
            self.emit("# Main text")
            for ch in body:
                self.block(ch, 2)
        else:
            self.c.warnings.append("no <body> in XML (publisher may not allow full text)")
        if back is not None:
            for ch in back:
                self.block(ch, 2)
        for sa in art.findall("sub-article"):
            self.emit()
            self.emit("# Sub-article")
            b = sa.find("body")
            if b is not None:
                for ch in b:
                    self.block(ch, 2)
        txt = "\n".join(self.out)
        self.c.text = re.sub(r"\n{3,}", "\n\n", txt).strip() + "\n"
        return self.c

    def front(self, front):
        am = front.find("article-meta")
        jm = front.find("journal-meta")
        if am is None:
            return
        t = am.find(".//article-title")
        if t is not None:
            self.emit("# " + self.inline(t))
        auth = []
        for c in am.iter("contrib"):
            if c.get("contrib-type", "author") != "author":
                continue
            sn, gn = c.findtext(".//surname"), c.findtext(".//given-names")
            if sn:
                auth.append(f"{gn or ''} {sn}".strip())
        if auth:
            self.emit("Authors: " + "; ".join(auth))
        if jm is not None:
            jt = jm.findtext(".//journal-title")
            yr = am.findtext(".//pub-date/year")
            self.emit(f"Journal: {jt or ''}  Year: {yr or ''}")
        for aid in am.findall("article-id"):
            if aid.get("pub-id-type") not in ("pmid", "pmcid", "doi"):
                continue
            self.emit(f"{aid.get('pub-id-type', 'id').upper()}: {_ws(aid.text or '')}")
        for ab in am.findall("abstract"):
            self.emit()
            ttl = ab.find("title")
            atype = ab.get("abstract-type")
            self.emit("## " + (self.inline(ttl) if ttl is not None else
                               "Abstract" + (f" ({atype})" if atype else "")))
            for ch in ab:
                if _ln(ch) == "title":
                    continue
                if _ln(ch) == "sec":
                    t = ch.find("title")
                    body = " ".join(self.inline(p) for p in ch.findall("p"))
                    self.emit(((self.inline(t) + ": ") if t is not None else "") + body)
                else:
                    self.block(ch, 3)
        kws = [self.inline(k) for k in am.iter("kwd")]
        if kws:
            self.emit("Keywords: " + "; ".join(kws))


def jats_to_text(xml: bytes | str) -> Converted:
    if isinstance(xml, str):
        xml = xml.encode("utf-8")
    parser = etree.XMLParser(resolve_entities=False, load_dtd=False, no_network=True, recover=True,
                             huge_tree=True)
    root = etree.fromstring(xml, parser)
    if root is None:
        raise ValueError("unparseable XML")
    if _ln(root) == "pmc-articleset":
        root = root.find("article")
    return _Jats(root).run()


# ----------------------------------------------------------------------------- PMC HTML
class _Html:
    SKIP_CLASSES = ("usa-link--external", "img-box", "tileshop")

    def __init__(self, root):
        self.root = root
        self.out: list[str] = []
        self.c = Converted()

    def inline(self, el) -> str:
        parts: list[str] = []

        def rec(e):
            if e.text:
                parts.append(e.text)
            for ch in e:
                tag = _ln(ch)
                cls = ch.get("class", "") or ""
                if tag in ("script", "style", "img", "button", "svg"):
                    pass
                elif tag == "sup":
                    t = _ws("".join(ch.itertext()))
                    parts.append(("^" + t) if t and len(t) < 12 else (" " + t))
                elif tag == "sub":
                    t = _ws("".join(ch.itertext()))
                    parts.append("_" + t if t else "")
                elif tag == "br":
                    parts.append(" ")
                elif tag == "a" and "Open in a new tab" in "".join(ch.itertext()):
                    pass
                else:
                    rec(ch)
                if ch.tail:
                    parts.append(ch.tail)

        rec(el)
        return _ws("".join(parts))

    def emit(self, s=""):
        self.out.append(s)

    def block(self, el, level: int):
        tag = _ln(el)
        cls = el.get("class", "") or ""
        if tag in ("script", "style", "nav", "button", "hr", "img"):
            return
        if re.fullmatch(r"h[1-6]", tag):
            t = self.inline(el)
            if t:
                self.emit()
                self.emit("#" * min(level, 6) + " " + t)
        elif tag == "p" or (tag == "div" and "p" in cls.split() and "text-right" not in cls):
            if "img-box" in cls:
                return
            t = self.inline(el)
            if t:
                self.emit(t)
        elif tag == "figure":
            self.c.n_figures += 1
            fid = el.get("id", f"fig{self.c.n_figures}")
            lab = el.find(".//h2") if el.find(".//h2") is not None else el.find(".//h3")
            cap = el.find(".//figcaption")
            self.emit()
            self.emit(f"[FIGURE {fid}] " + (self.inline(lab) if lab is not None else "") + " " +
                      (self.inline(cap) if cap is not None else ""))
        elif tag == "section" and "tw" in cls.split() or (tag == "div" and "tw" in cls.split()):
            self.c.n_tables += 1
            tid = el.get("id", f"table{self.c.n_tables}")
            lab = next((e for e in el if re.fullmatch(r"h[1-6]", _ln(e))), None)
            cap = el.find(".//div[@class='caption p']")
            self.emit()
            self.emit(f"[TABLE {tid}] " + (self.inline(lab) if lab is not None else "") + " " +
                      (self.inline(cap) if cap is not None else ""))
            tbl = next(el.iter("table"), None)
            if tbl is None:
                self.emit("[table is an image; contents not extractable]")
                self.c.warnings.append(f"table {tid} is an image")
            else:
                tsv = table_to_tsv(tbl, self.inline)
                self.c.tables[tid] = tsv
                self.emit(tsv)
            for fn in el.xpath(".//*[contains(@class,'fn')]"):
                t = self.inline(fn)
                if t:
                    self.emit("[table note] " + t)
            self.emit(f"[/TABLE {tid}]")
        elif tag == "table":
            self.c.n_tables += 1
            tid = f"table{self.c.n_tables}"
            tsv = table_to_tsv(el, self.inline)
            self.c.tables[tid] = tsv
            self.emit(f"[TABLE {tid}]")
            self.emit(tsv)
            self.emit(f"[/TABLE {tid}]")
        elif tag in ("ul", "ol"):
            for li in el.findall("li"):
                self.emit("- " + self.inline(li))
        elif tag == "section" and ("ref-list" in cls.split()):
            self.emit()
            self.emit("## References")
            for li in el.iter("li"):
                cite = li.find("cite")
                self.emit(f"[{li.get('id', '')}] " + self.inline(cite if cite is not None else li))
        elif tag == "section" and "supplementary-materials" in cls:
            self.emit()
            self.emit("## Supplementary Material")
            for a in el.iter("a"):
                h = a.get("href")
                if h:
                    self.c.supp_hrefs.append(h)
                    self.emit(f"[SUPPLEMENT] {self.inline(a)} (file: {h})")
        else:
            deeper = level + 1 if tag == "section" else level
            for ch in el:
                self.block(ch, deeper)

    def run(self) -> Converted:
        root = self.root
        h1 = root.find(".//h1")
        if h1 is not None:
            self.emit("# " + self.inline(h1))
        cit = root.xpath("//section[@aria-label='Article citation and metadata']")
        meta = root.xpath("//meta[@name='citation_author']/@content")
        if meta:
            self.emit("Authors: " + "; ".join(meta))
        for n, label in (("citation_journal_title", "Journal"), ("citation_publication_date", "Date"),
                         ("citation_doi", "DOI"), ("citation_pmid", "PMID")):
            v = root.xpath(f"//meta[@name='{n}']/@content")
            if v:
                self.emit(f"{label}: {v[0]}")
        for ab in root.xpath("//section[contains(@class,'abstract')]")[:1]:
            self.block(ab, 1)
        body = root.xpath("//section[contains(@class,'body')]")
        if not body:
            self.c.warnings.append("no article body in HTML page")
        for b in body:
            self.emit()
            self.emit("# Main text")
            for ch in b:
                if "abstract" in (ch.get("class") or ""):
                    continue
                self.block(ch, 1)
        txt = "\n".join(self.out)
        self.c.text = re.sub(r"\n{3,}", "\n\n", txt).strip() + "\n"
        return self.c


def pmc_html_to_text(page: bytes | str) -> Converted:
    root = lhtml.fromstring(page)
    art = root.find(".//main")
    return _Html(root if art is None else root).run()


def tsv_blocks_from_text(text: str) -> dict:
    """Utility: re-extract {table_id: tsv} from converted text."""
    out = {}
    for m in re.finditer(r"\[TABLE (\S+)\][^\n]*\n(.*?)\n\[/TABLE \1\]", text, re.S):
        out[m.group(1)] = m.group(2)
    return out
