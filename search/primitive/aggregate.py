"""Combine the per-method outputs of the primitive-group SRG search into one status per group.

Usage: python aggregate.py  (run in search/primitive/)
Inputs: prim_ranks.csv, open_srg.csv, dfs_results.jsonl.gz (DFS), cpsat_results.jsonl (fixed-parameter CP-SAT),
milp_results.jsonl (commutative eigenvalue MILP), hybrid_results.jsonl (non-commutative MILP).
Groups of degree p^2 are settled by the multiplier argument in primitive_note.md (all are affine, checked in GAP).
A group counts as settled for an open parameter set if some method gives a complete answer for it (DFS: always;
CP-SAT/MILP: status INFEASIBLE, FOUND or 'k not a sum of unit valencies').
Output: summary.csv (v, i, rank, method, settled, hits).
"""
import csv, json, gzip
from collections import defaultdict

SQ = {121, 169, 289, 361, 529, 841, 961}
opens = defaultdict(set)
for r in csv.DictReader(open('open_srg.csv')):
    v, k, l, m = int(r['v']), int(r['k']), int(r['l']), int(r['m'])
    opens[v] |= {(v, k, l, m), (v, v - k - 1, v - 2 - 2 * k + m, v - 2 * k + l)}
groups = [(int(r['v']), int(r['i']), int(r['rank'])) for r in csv.DictReader(open('prim_ranks.csv'))]
settled = defaultdict(set); method = defaultdict(set); hits = defaultdict(list)
for l in gzip.open('dfs_results.jsonl.gz', 'rt'):
    r = json.loads(l); key = (r['v'], r['i'])
    settled[key] |= opens[r['v']]; method[key].add('dfs')
    for g in r['graphs']:
        if 'srg' in g and tuple(g['srg']) in opens[r['v']]: hits[key].append(g['srg'])
for fn, name in [('cpsat_results.jsonl', 'cpsat'), ('milp_results.jsonl', 'milp'), ('hybrid_results.jsonl', 'hybrid')]:
    try:
        lines = open(fn).read().splitlines()
    except FileNotFoundError:
        continue
    for l in lines:
        r = json.loads(l)
        ents = [r['partial']] if 'partial' in r else r.get('results', [])
        key = (r['v'], r['i'])
        for e in ents:
            if e['status'] in ('INFEASIBLE', 'FOUND', 'OPTIMAL', 'k not a sum of unit valencies'):
                settled[key].add(tuple(e['param'])); method[key].add(name)
            if e['status'] in ('FOUND', 'OPTIMAL') and tuple(e['param']) in opens[r['v']]:
                hits[key].append(e['param'])
rows = []
for v, i, rk in groups:
    key = (v, i)
    if v in SQ:
        settled[key] |= opens[v]; method[key].add('multiplier')
    ok = opens[v] <= settled[key]
    rows.append({'v': v, 'i': i, 'rank': rk, 'method': '+'.join(sorted(method[key])), 'settled': ok,
                 'hits': json.dumps(hits[key])})
with open('summary.csv', 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
print('groups', len(rows), 'settled', sum(r['settled'] for r in rows),
      'unsettled', [(r['v'], r['i']) for r in rows if not r['settled']][:30],
      'hits', [(r['v'], r['i'], r['hits']) for r in rows if r['hits'] != '[]'])
