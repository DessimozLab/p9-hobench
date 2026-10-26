#!/usr/bin/env python
"""Compare blind verifier answers with curator claims (CURATION_PROTOCOL s4 kill rule).
Reads families/<dir>/verification/{sheet,key}[_r2].tsv, answers[_r2].tsv, optional review.tsv (manual verdict overrides:
item_id<TAB>verdict(AGREE|DISAGREE)<TAB>note). Writes verification/compare.tsv and prints a per-family table.
Auto rules: member labels agree when normalised strings contain one another or share a COG id / subfamily token;
events agree when the verifier's interval ([lo, hi] or lvl=X) equals the claimed (lvl | lvl_lo,lvl_hi), or tag-only equal;
decoy/distinct_from agree when verifier answered 'non-homologous'. Everything else -> REVIEW (needs manual review.tsv)."""
import csv, glob, os, re, sys
from scipy.stats import beta
csv.field_size_limit(10**9)
def rd(p):
    return list(csv.DictReader(open(p), delimiter='\t')) if os.path.exists(p) else []
def norm(s): return re.sub(r'[^a-z0-9]', '', s.lower())
def parse_claim(c):
    d = dict(kv.split('=', 1) for kv in c.split(';') if '=' in kv)
    if 'lvl' in d: return (d['lvl'], d['lvl'])
    if 'lvl_lo' in d: return (d['lvl_lo'], d.get('lvl_hi', ''))
    return None
def parse_ans(a):
    m = re.search(r'\[(\w+), ?(\w+)\]', a)
    if m: return (m.group(1), m.group(2))
    m = re.search(r'lvl=(\w+)', a)
    if m: return (m.group(1), m.group(1))
    return None
def member_agree(cl, an):
    cn, an_ = norm(cl), norm(an)
    if not an_ or 'notfound' in an_: return False
    for cog in re.findall(r'COG\d+', an):
        if cog in cl: return True
    if cn and (an_ in cn or cn in an_): return True
    return False
def judge(kind, claimed, ans):
    if kind == 'member': return 'AGREE' if member_agree(claimed, ans) else 'REVIEW'
    if kind in ('decoy', 'distinct_from'): return 'AGREE' if ans.lower().startswith('non-homologous') else 'REVIEW'
    if kind == 'event':
        pc, pa = parse_claim(claimed), parse_ans(ans)
        if pc and pa and pc == pa: return 'AGREE'
        return 'REVIEW'
    return 'REVIEW'
def cp(k, n, a=0.05):  # two-sided 95% Clopper-Pearson upper bound
    return 1.0 if k >= n else float(beta.ppf(1 - a / 2, k + 1, n - k))
rows = []
for vd in sorted(glob.glob('families/*/verification') + glob.glob('families/_verif_*')):
    fam = os.path.basename(vd) if os.path.basename(vd).startswith('_verif_') else os.path.basename(os.path.dirname(vd))
    for rnd in ('', '_r2'):
        sh, ky, an = rd(f'{vd}/sheet{rnd}.tsv'), {r['item_id']: r for r in rd(f'{vd}/key{rnd}.tsv')}, {r['item_id']: r for r in rd(f'{vd}/answers{rnd}.tsv')}
        rv = {r['item_id']: r for r in rd(f'{vd}/review{rnd}.tsv')}
        out = []
        for s in sh:
            i = s['item_id']
            if i not in an: continue
            k = ky[i]; v = judge(k['kind'], k['claimed'], an[i]['answer'])
            if i in rv: v = rv[i]['verdict']
            out.append((i, s['family'], k['kind'], v, k['claimed'], an[i]['answer'], rv.get(i, {}).get('note', '')))
        if out:
            with open(f'{vd}/compare{rnd}.tsv', 'w') as f:
                f.write('item_id\tfamily\tkind\tverdict\tclaimed\tanswer\tnote\n')
                for o in out: f.write('\t'.join(o).replace('\n', ' ') + '\n')
            rows.append((fam, rnd or 'r1', out))
print('family\tround\tn\tdisagree\treview_pending\tCP95_upper\tverdict')
for fam, rnd, out in rows:
    n = len(out); dis = sum(o[3] == 'DISAGREE' for o in out); pend = sum(o[3] == 'REVIEW' for o in out)
    ub = cp(dis, n)
    verdict = 'PENDING' if pend else ('PASS' if dis == 0 else ('RESAMPLE' if dis == 1 else 'FAIL'))
    print(f'{fam}\t{rnd}\t{n}\t{dis}\t{pend}\t{ub:.3f}\t{verdict}')
