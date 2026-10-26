#!/usr/bin/env python
"""Blind verification sheet generator (CURATION_PROTOCOL s4).
Usage: make_verification_sheet.py [--extra]   (--extra: draw 30 further members, resample round 2)
For each family / collection: seeded sample of min(all,30) member-label rows (round-robin over labels),
plus ALL D/H events (from reference_tree.nhx), decoys and distinct_from rows. Outputs verification/sheet.tsv
(blinded: no label/level/tag) and key.tsv (claimed values; do not open until answers filled)."""
import csv, os, re, sys, random, glob
F = 'families'; SEED = 20261026
EXTRA = '--extra' in sys.argv
ROUND = 2 if EXTRA else 1
csv.field_size_limit(10**9)

def rd(p):
    if not os.path.exists(p): return []
    with open(p) as f: return list(csv.DictReader(f, delimiter='\t'))

def parse_nhx(s):
    """return list of nodes: dict(attrs, leaves)"""
    s = s.strip().rstrip(';'); nodes = []; stack = [[]]; i = 0; n = len(s)
    while i < n:
        c = s[i]
        if c == '(':
            stack.append([]); i += 1
        elif c in ',':
            i += 1
        elif c == ')':
            kids = stack.pop(); i += 1
            m = re.match(r'([^\[\(\),:;]*)(\[&&NHX:([^\]]*)\])?(:[0-9.eE+-]+)?', s[i:])
            name, _, attrs, _ = m.groups(); i += m.end()
            if attrs is None:
                m2 = re.match(r'\[&&NHX:([^\]]*)\]', s[i:])
                if m2: attrs = m2.group(1); i += m2.end()
            leaves = [l for k in kids for l in k]
            d = dict(kv.split('=', 1) for kv in (attrs or '').split(':') if '=' in kv)
            d['_name'] = name
            nodes.append((d, leaves)); stack[-1].append(leaves)
        else:
            m = re.match(r'([^\[\(\),:;]+)(\[&&NHX:[^\]]*\])?(:[0-9.eE+-]+)?', s[i:])
            stack[-1].append([m.group(1)]); i += m.end()
    return nodes

def mask_loc(l):
    parts = [p for p in re.split(r';\s*|\s*\|\s*', l)]
    keep = [p for p in parts if '=' not in p]
    return '; '.join(keep[:3]) if keep else re.split(r'[:,]', l)[0]

def family_items(dirp, fam, nmem, rng, unit=None, exclude=()):
    mem = rd(f'{dirp}/members.tsv'); ass = rd(f'{dirp}/assertions.tsv')
    sp = {m['gene_id']: m for m in mem}
    items = []
    # members: round robin over labels
    by = {}
    for m in mem: by.setdefault(m['label'], []).append(m)
    for v in by.values(): rng.shuffle(v)
    labs = list(by); rng.shuffle(labs); picked = []
    while len(picked) < min(nmem + len(exclude), len(mem)) and any(by.values()):
        for l in labs:
            if by[l]: picked.append(by[l].pop())
            if len(picked) >= nmem + len(exclude): break
    picked = [m for m in picked if m['gene_id'] not in exclude][:nmem]
    for m in picked:
        items.append(dict(kind='member', genes=m['gene_id'], species=m['species'], orig=m['original_id_in_paper'],
                          source=m['source'], locator=mask_loc(m['locator']), node='', claimed=m['label'], unit=unit or ''))
    # events from tree
    tp = f'{dirp}/reference_tree.nhx'
    if os.path.exists(tp):
        for d, leaves in parse_nhx(open(tp).read().strip().splitlines()[0]):
            if d.get('ev') not in ('D', 'H'): continue
            nid = d.get('id', '')
            locs, srcs = [], set()
            for a in ass:
                if a.get('type') == 'event' and nid in (a.get('assertion_id'), a.get('target'), a.get('what')):
                    locs.append(a.get('locator', '')); srcs.add(a.get('source_id', ''))
            cl = ';'.join(f'{k}={v}' for k, v in d.items() if k in ('ev', 'lvl', 'lvl_lo', 'lvl_hi', 'tag', 'donor', 'recipient') )
            items.append(dict(kind='event', genes=','.join(leaves[:8]) + (f',...(+{len(leaves)-8})' if len(leaves) > 8 else ''),
                              species=','.join(sorted({sp[g]['species'] for g in leaves if g in sp}))[:60], orig=','.join(sp[g]['original_id_in_paper'] for g in leaves[:6] if g in sp),
                              source=','.join(sorted(srcs)) or d.get('src', ''), locator=' || '.join(locs[:3]) or '(none)', node=nid, claimed=cl, unit=unit or ''))
    for a in ass:
        t = a.get('type')
        if t == 'decoy':
            items.append(dict(kind='decoy', genes=a.get('genes') or a.get('target', ''), species='', orig='', source=a.get('source_id', ''),
                              locator=a.get('locator', ''), node=a['assertion_id'], claimed='nonhomologous', unit=unit or ''))
        elif t == 'distinct_from':
            items.append(dict(kind='distinct_from', genes=a.get('target', '').split(' at_level')[0], species='', orig='', source=a.get('source_id', ''),
                              locator=a.get('locator', ''), node=a['assertion_id'], claimed=a.get('target', ''), unit=unit or ''))
    return items

def write(outdir, fam, items):
    os.makedirs(outdir, exist_ok=True)
    sp = f'{outdir}/sheet.tsv'; kp = f'{outdir}/key.tsv'
    rnd = '' if ROUND == 1 else '_r2'
    sp = f'{outdir}/sheet{rnd}.tsv'; kp = f'{outdir}/key{rnd}.tsv'
    cols = ['item_id', 'family', 'unit', 'kind', 'genes', 'species', 'original_id_in_paper', 'source', 'locator', 'node_id']
    with open(sp, 'w') as s, open(kp, 'w') as k:
        s.write('\t'.join(cols) + '\n'); k.write('item_id\tkind\tclaimed\n')
        for i, it in enumerate(items, 1):
            iid = f'{fam}.{ROUND}.{i:03d}'
            s.write('\t'.join([iid, fam, it['unit'], it['kind'], it['genes'], it['species'], it['orig'], it['source'], it['locator'], it['node']]).replace('\n', ' ') + '\n')
            k.write(f"{iid}\t{it['kind']}\t{it['claimed']}\n")
    print(fam, len(items), 'items', sum(i['kind'] == 'member' for i in items), 'members')

def prev_members(outdir, fam):
    p = f'{outdir}/sheet.tsv'
    return {r['genes'] for r in rd(p) if r['kind'] == 'member'} if EXTRA else set()

def rng_for(name): return random.Random(f'{SEED}:{name}:{ROUND}')

fams = sorted(d for d in os.listdir(F) if os.path.isdir(f'{F}/{d}') and not d.startswith('_') and os.path.exists(f'{F}/{d}/reference_tree.nhx') or d.startswith('Y01'))
singles = [d for d in fams if not d.startswith(('K01_RPROT', 'K04', 'Y01'))]
for d in singles:
    out = f'{F}/{d}/verification'; r = rng_for(d)
    write(out, d, family_items(f'{F}/{d}', d, 30, r, exclude=prev_members(out, d)))
# collections: sample 30 units, 1 member per unit, all specials in sampled units
def collection(name, unitdirs, out, nunits=30):
    r = rng_for(name); ud = sorted(u for u in unitdirs if len(rd(f'{u}/members.tsv'))>0)
    if EXTRA:
        old = {x['unit'] for x in rd(f'{out}/sheet.tsv')}
        ud = [u for u in ud if os.path.basename(u) not in old]
    r.shuffle(ud); items = []
    for u in ud[:nunits]:
        items += family_items(u, os.path.basename(u), 1, r, unit=os.path.basename(u))
    write(out, name, items)
collection('Y01_YGOB_pillars', glob.glob(f'{F}/Y01_YGOB_pillars/units/*'), f'{F}/Y01_YGOB_pillars/verification')
collection('K01_RPROT', glob.glob(f'{F}/K01_RPROT_*'), f'{F}/_verif_K01')
# K04a-j: one pooled sample
r = rng_for('K04'); k4 = sorted(glob.glob(f'{F}/K04?_*')); pool = []
allm = []
for u in k4:
    mem = rd(f'{u}/members.tsv')
    for m in mem: m['_u'] = os.path.basename(u)
    allm += mem
r.shuffle(allm)
old = {x['genes'] for x in rd(f'{F}/_verif_K04/sheet.tsv')} if EXTRA else set()
sel = [m for m in allm if m['gene_id'] not in old][:30]
items = []
for u in k4:
    it = family_items(u, os.path.basename(u), 0, r, unit=os.path.basename(u))
    items += [i for i in it if not EXTRA or True]
for m in sel:
    items.append(dict(kind='member', genes=m['gene_id'], species=m['species'], orig=m['original_id_in_paper'], source=m['source'],
                      locator=mask_loc(m['locator']), node='', claimed=m['label'], unit=m['_u']))
if EXTRA: items = [i for i in items if i['kind'] == 'member']
write(f'{F}/_verif_K04', 'K04', items)
