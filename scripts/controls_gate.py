#!/usr/bin/env python
"""Gate: run hobench.controls on real frozen-candidate families.
positive = perfect on every pair metric; all-singletons / one-big-HOG = 0 recovery of testable events."""
import sys
from pathlib import Path
ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from hobench import controls
from hobench.metrics import evaluate_family, family_scalars, panel_context
from hobench.reference import load_families
from hobench.speciestree import SpeciesTree

FAMS = {"Y03_URA1_HGT": "Y", "E02_ARF": "E", "V03_Hox": "V"}
PERFECT = ("merge_accuracy_node", "merge_accuracy_pairs", "level_bcubed_f1", "membership_f1", "orthobench_f")


def main():
    bad = []
    for fid, panel in FAMS.items():
        st = SpeciesTree.from_newick(str(ROOT / "panels" / panel / "species_tree.nwk"))
        fams = [f for f in load_families(ROOT / "families", panel=panel) if f.id == fid or f.id.startswith(fid + "/")]
        assert fams, fid
        ctx = panel_context(fams, st)
        for name, pred in (("positive", controls.reference_prediction(fams, st)),
                           ("singletons", controls.all_singletons(fams)),
                           ("one_big", controls.one_big_hog(fams, st))):
            for f in fams:
                r = evaluate_family(f, st, pred, ctx)
                sc = family_scalars(r)
                ev = [e for e in r["events"] if e["testable"]]
                rec = [e["recovered_0.5"] for e in ev]
                line = f"{f.id:28s} {name:10s} events={len(r['events'])} testable={len(ev)} recovered={sum(rec)}"
                if name == "positive":
                    nbad = [k for k in PERFECT if sc[k] != 1.0] + ([] if all(rec) else ["events"])
                    line += f"  inconsistent_pairs={pred.meta['inconsistent_pairs']} " + ("PASS" if not nbad and not pred.meta["inconsistent_pairs"] else f"FAIL {nbad}")
                    if nbad or pred.meta["inconsistent_pairs"]: bad.append(line)
                else:
                    # one-big-HOG merges same-species pairs at parent(species): a duplication AT a species
                    # leaf whose reference merge is parent(species) is then correct by construction
                    # (single-gene copies, vacuous cohesion). Reported separately; any other recovery fails.
                    leaf = {e["event"] for e in ev if st.is_species(e["level"])}
                    other = [e["event"] for e in ev if e["recovered_0.5"] and e["event"] not in leaf]
                    nleaf = sum(1 for e in ev if e["recovered_0.5"] and e["event"] in leaf)
                    ok = not other and (name == "one_big" or nleaf == 0)
                    line += f"  (of which species-leaf dups: {nleaf})  " + ("PASS" if ok else f"FAIL {other[:3]}")
                    if not ok: bad.append(line)
                print(line, flush=True)
    print("GATE", "FAILED" if bad else "PASSED")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
