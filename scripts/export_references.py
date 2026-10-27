#!/usr/bin/env python
"""Export the frozen references of a panel (read-only on families/).

    python scripts/export_references.py --panel Y [--out release/references]

Writes release/references/<panel>/
  reference_hogs.orthoxml.gz   OrthoXML 0.5, asserted structure only.  One top-level orthologGroup per family
                               (or per resolved subtree).  S node = orthologGroup at the node's level; D/H node =
                               orthologGroup at the representative merge level wrapping a paralogGroup, with
                               properties lvl_lo/lvl_hi/merge_lo/merge_hi (names; ABOVE_RANGE/NEVER for unbounded)
                               and event=HGT for H.  Children of a U node are exported as separate top-level
                               groups carrying <property name="unresolved" value="true"/> (relations between them
                               are not asserted).  Pairs never merged (merge = NEVER) are likewise separate groups.
  asserted_pairs.tsv.gz        gene_a gene_b relation(ortholog|paralog|xenolog) merge_lo merge_hi tier family
"""
from __future__ import annotations
import argparse, gzip, sys
from pathlib import Path
import numpy as np
from lxml import etree

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from hobench.reference import load_families, reconcile, INF, NOT_ASSERTED, ABOVE   # noqa: E402
from hobench.metrics import reference_matrices, representative_matrix            # noqa: E402
from hobench.speciestree import SpeciesTree                                      # noqa: E402

NS = "http://orthoXML.org/2011/"


def lvl_name(st, i):
    return "NEVER" if i == INF else st.names[i]


def prop(parent, k, v):
    etree.SubElement(parent, "property", name=k, value=str(v))


def fam_groups(rf, st, rep, gid, grouped=False):
    """Return list of top-level orthologGroup elements for one reconciled family."""
    ev_of = {id(e.node): e for e in rf.events}
    leaf_idx = {g: i for i, g in enumerate(rf.genes)}
    detached: list = []
    sp_of_gid = {gid[g]: int(rf.species_idx[i]) for g, i in leaf_idx.items()}

    def leaves(nd):
        return [l.name for l in nd.preorder() if l.is_leaf]

    def wrap(el, level_i, unresolved=False):
        """Make sure el is an orthologGroup (geneRef -> single-gene group at species level)."""
        if el.tag == "geneRef":
            level_i = sp_of_gid[int(el.get("id"))]
            g = etree.Element("orthologGroup", taxonId=str(level_i))
            g.append(el)
            el = g
        if unresolved:
            el.insert(0, etree.Element("property", name="unresolved", value="true"))
        return el

    def emit(nd):
        if nd.is_leaf:
            return etree.Element("geneRef", id=str(gid[nd.name]))
        e = ev_of[id(nd)]
        kids = [emit(c) for c in nd.children]
        if e.ev == "S":
            kids = [k for k in kids if k is not None]
            if not kids:
                return None
            el = etree.Element("orthologGroup", taxonId=str(e.mapping))
            el.extend(kids)
            return el
        if e.ev == "U" and grouped:
            # positive-control variant: U subtree grouped at the oldest level implied inside it (as controls.py does)
            kk = [(c, k) for c, k in zip(nd.children, kids) if k is not None]
            idxs = [np.array([leaf_idx[x] for x in leaves(c)]) for c, _ in kk]
            lv = [int(v) for a in range(len(idxs)) for b in range(a + 1, len(idxs))
                  for v in np.unique(rep[np.ix_(idxs[a], idxs[b])])]
            if len(kk) >= 2 and lv and INF not in lv and min(lv) >= 0:
                M = min(lv, key=lambda v: st.depth[v])
                el = etree.Element("orthologGroup", taxonId=str(M))
                prop(el, "unresolved", "true")
                el.extend([k for _, k in kk])
                return el
            if len(kk) == 1:
                return kk[0][1]
        if e.ev == "U":
            for c, k in zip(nd.children, kids):
                if k is not None:
                    sp = rf.species_idx[leaf_idx[leaves(c)[0]]]
                    detached.append(wrap(k, int(sp), unresolved=True))
            return None
        # D / H
        i0, j0 = leaf_idx[leaves(nd.children[0])[0]], leaf_idx[leaves(nd.children[1])[0]]
        M = int(rep[i0, j0])
        if M == INF or M < 0:
            for c, k in zip(nd.children, kids):
                if k is not None:
                    sp = rf.species_idx[leaf_idx[leaves(c)[0]]]
                    detached.append(wrap(k, int(sp)))
            return None
        kids = [k for k in kids if k is not None]
        el = etree.Element("orthologGroup", taxonId=str(M))
        prop(el, "lvl_lo", lvl_name(st, e.lvl_lo) if e.lvl_lo != ABOVE else "ABOVE_RANGE")
        prop(el, "lvl_hi", "ABOVE_RANGE" if e.lvl_hi == ABOVE else st.names[e.lvl_hi])
        prop(el, "merge_lo", lvl_name(st, e.merge_lo))
        prop(el, "merge_hi", lvl_name(st, e.merge_hi))
        if e.ev == "H":
            prop(el, "event", "HGT")
            prop(el, "donor", e.donor)
        pg = etree.SubElement(el, "paralogGroup")
        pg.extend(kids)
        return el

    root = emit(rf.family.tree)
    out = [wrap(root, int(rf.species_idx[0]))] if root is not None else []
    out += detached
    for k, el in enumerate(out):
        el.set("id", f"{rf.id}|{k}")
        prop(el, "family", rf.id)
    return out


def export_panel(P: str, famdir: Path, outdir: Path):
    _export(P, famdir, outdir, True, write_pairs=False)   # positive-control variant (U subtrees grouped)
    return _export(P, famdir, outdir, False)


def _export(P, famdir, outdir, grouped, write_pairs=True):
    st = SpeciesTree.from_newick(str(ROOT / "panels" / P / "species_tree.nwk"))
    fams = load_families(famdir, panel=P)
    rfs = [reconcile(f, st, check=False) for f in fams]
    outdir.mkdir(parents=True, exist_ok=True)
    root = etree.Element("orthoXML", nsmap={None: NS}, version="0.5", origin="HOBench-reference", originVersion="1")
    # genes per species
    gid, per_sp, n = {}, {}, 0
    for f in fams:
        for g, s in f.gene_species.items():
            if g not in gid:
                n += 1; gid[g] = n; per_sp.setdefault(s, []).append(g)
    for s in sorted(per_sp):
        se = etree.SubElement(root, "species", name=s, NCBITaxId="0")
        db = etree.SubElement(se, "database", name="HOBench", version="1")
        ge = etree.SubElement(db, "genes")
        for g in per_sp[s]:
            etree.SubElement(ge, "gene", id=str(gid[g]), protId=g)
    tx = etree.SubElement(root, "taxonomy")

    def tax(i, parent):
        t = etree.SubElement(parent, "taxon", id=str(i), name=st.names[i])
        for c in np.flatnonzero(st.parent_idx == i):
            tax(int(c), t)
    tax(int(np.flatnonzero(st.parent_idx < 0)[0]), tx)
    groups = etree.SubElement(root, "groups")
    rows = []
    for rf in rfs:
        rep = representative_matrix(rf)
        for el in fam_groups(rf, st, rep, gid, grouped):
            groups.append(el)
        if not write_pairs:
            continue
        R = reference_matrices(rf)
        iu, ju = np.triu_indices(len(rf.genes), 1)
        keep = R.lo[iu, ju] != NOT_ASSERTED
        iu, ju = iu[keep], ju[keep]
        for i, j in zip(iu, ju):
            e = rf.events[R.node[i, j]]
            rel = "xenolog" if (R.xeno[i, j] or e.ev == "H") else ("ortholog" if e.ev == "S" else "paralog")
            rows.append((rf.genes[i], rf.genes[j], rel, lvl_name(st, int(R.lo[i, j])), lvl_name(st, int(R.hi[i, j])),
                         e.tier or "", rf.id))
    with gzip.open(outdir / ("reference_hogs.positive_control.orthoxml.gz" if grouped else "reference_hogs.orthoxml.gz"), "wb") as fh:
        fh.write(etree.tostring(root, xml_declaration=True, encoding="UTF-8", pretty_print=True))
    if not write_pairs:
        return len(fams), 0
    with gzip.open(outdir / "asserted_pairs.tsv.gz", "wt") as fh:
        fh.write("gene_a\tgene_b\trelation\tmerge_lo\tmerge_hi\ttier\tfamily\n")
        for r in rows:
            fh.write("\t".join(r) + "\n")
    return len(fams), len(rows)


def validate(P: str, famdir: Path, outdir: Path, n: int = 3) -> int:
    """Read the export back with hobench.predictions.read_orthoxml and score n families against their references."""
    from hobench.predictions import read_orthoxml
    from hobench.metrics import evaluate_family, family_scalars, panel_context, PRIMARY
    st = SpeciesTree.from_newick(str(ROOT / "panels" / P / "species_tree.nwk"))
    fams = load_families(famdir, panel=P)
    keep = set().union(*[set(f.member_ids) for f in fams])
    cs = set().union(*[set(f.complete_species) for f in fams])
    pred = read_orthoxml(outdir / P / "reference_hogs.positive_control.orthoxml.gz", st, normalize=None, keep_genes=keep, count_species=cs)
    ctx = panel_context(fams, st)
    ok = 0   # number of exactly-1.0 families
    for f in fams[:n]:
        r = evaluate_family(f, st, pred, ctx, per_gene=False)
        sc = family_scalars(r)
        vals = {k: sc.get(k) for k in PRIMARY}
        good = sum(v is not None for v in vals.values()) > 0 and all(v is None or abs(v - 1.0) < 1e-9 for v in vals.values())
        ok += good
        print(f"  validate {P} {f.id}: " + " ".join(f"{k}={v:.3f}" if isinstance(v, float) else f"{k}={v}" for k, v in vals.items()),
              "OK" if good else "differs (interval-event representative placement / U grouping)")
    return ok


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--panel", nargs="+", default=["E", "K", "P", "V", "Y"])
    ap.add_argument("--families", default=str(ROOT / "families"))
    ap.add_argument("--out", default=str(ROOT / "release" / "references"))
    ap.add_argument("--validate", action="store_true", help="round-trip check: export scores 1.0 on 3 families per panel")
    ap.add_argument("--no-export", action="store_true")
    a = ap.parse_args()
    n_ok = 0
    for P in a.panel:
        if not a.no_export:
            nf, npairs = export_panel(P, Path(a.families), Path(a.out) / P)
            print(f"{P}: {nf} families, {npairs} asserted pairs", flush=True)
        if a.validate:
            n_ok += validate(P, Path(a.families), Path(a.out))
    if a.validate:
        print(f"round-trip exact (all primary endpoints 1.0): {n_ok} families of {3 * len(a.panel)} checked")
    return 0 if (not a.validate or n_ok >= 3) else 1


if __name__ == "__main__":
    sys.exit(main())
