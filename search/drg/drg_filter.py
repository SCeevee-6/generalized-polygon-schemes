"""Vertex-primitive distance-regular graphs: which arrays from Brouwer's tables can a primitive group carry?

Step 1 (this file): cheap filter. A distance-regular graph of diameter d whose automorphism group contains the
primitive group G is a union of orbitals of G, and so is each distance graph. So G needs rank >= d + 1 and every
valency k_1, ..., k_d has to be a sum of nontrivial suborbit lengths of G; in fact the nontrivial suborbits must
split into d groups with sums k_1, ..., k_d (tested when G has at most 61 suborbits). A vertex-primitive graph is
neither bipartite nor antipodal, so only the primitive arrays matter.

Inputs: primitive_arrays.json (from fetch_arrays.py), suborbits.txt.gz (from export_suborbits.g; one line
"G v i |G| [suborbit lengths]" per primitive group other than Sym(v), Alt(v)).
Output: survivors.csv (v, i, order, rank, array, hamming) for the pairs that pass. The last column is H(d,q) when
the array is that of a Hamming graph with q != 4; those graphs are determined by their array (Egawa; BCN 9.2.5),
so step 2 skips them.

Usage: python drg_filter.py
"""
import os, re, csv, json, gzip
from collections import defaultdict
from fractions import Fraction as Fr

HERE = os.path.dirname(os.path.abspath(__file__))


def valencies(b, c):
    k = [Fr(1)]
    for i in range(len(b)):
        k.append(k[-1] * b[i] / c[i])
    return [int(x) for x in k]


def subset_sum(lengths, target):
    bits = 1
    for L in lengths:
        bits |= bits << L
    return (bits >> target) & 1 == 1


def splits(lengths, targets):
    """can the multiset `lengths` be split into parts with sums `targets` (in order)? First each target has to be
    a subset sum. Then, for at most 60 suborbits, a DP over reachable tuples of the first len(targets)-1 partial
    sums (the last part gets the rest). Returns True, False, or None if the joint test was skipped."""
    assert sum(lengths) == sum(targets)
    if not all(subset_sum(lengths, k) for k in targets):
        return False
    if len(lengths) > 60:
        return None
    d = len(targets)
    states = {(0,) * (d - 1)}
    for L in sorted(lengths, reverse=True):
        new = set()
        for st in states:
            new.add(st)                       # put L into the last part
            for j in range(d - 1):
                if st[j] + L <= targets[j]:
                    new.add(st[:j] + (st[j] + L,) + st[j + 1:])
        states = new
        if len(states) > 300000:
            return None
    return tuple(targets[:-1]) in states


def hamming(b, c):
    d = len(b)
    q1 = b[-1]
    if c == list(range(1, d + 1)) and b == [(d - i) * q1 for i in range(d)] and q1 + 1 != 4:
        return 'H(%d,%d)' % (d, q1 + 1)
    return ''


def main():
    arrays = defaultdict(list)
    for a in json.load(open(os.path.join(HERE, 'primitive_arrays.json'))):
        arrays[a['v']].append(a)
    rows = []; stats = defaultdict(int)
    for l in gzip.open(os.path.join(HERE, 'suborbits.txt.gz'), 'rt'):
        if not l.startswith('G '):
            continue
        m = re.match(r'G (\d+) (\d+) (\d+) \[(.*)\]', l.strip())
        v, i, order = int(m.group(1)), int(m.group(2)), int(m.group(3))
        sub = sorted(int(x) for x in m.group(4).replace(' ', '').split(','))
        assert sub[0] == 1 and sum(sub) == v
        stats['groups'] += 1
        for a in arrays[v]:
            stats['pairs'] += 1
            k = valencies(a['b'], a['c'])
            d = a['d']
            if len(sub) < d + 1:
                stats['rank too small'] += 1; continue
            ok = splits(sub[1:], k[1:])
            if ok is False:
                stats['valencies do not split'] += 1; continue
            stats['pass' if ok else 'pass (split undecided)'] += 1
            rows.append({'v': v, 'i': i, 'order': order, 'rank': len(sub),
                         'array': '{' + ','.join(map(str, a['b'])) + ';' + ','.join(map(str, a['c'])) + '}',
                         'hamming': hamming(a['b'], a['c'])})
    with open(os.path.join(HERE, 'survivors.csv'), 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=['v', 'i', 'order', 'rank', 'array', 'hamming']); w.writeheader(); w.writerows(rows)
    print(dict(stats))
    print('survivor pairs', len(rows), 'groups', len({(r['v'], r['i']) for r in rows}),
          'arrays', len({(r['v'], r['array']) for r in rows}))
    rest = [r for r in rows if not r['hamming']]
    print('without Hamming arrays (q != 4): pairs', len(rest), 'groups', len({(r['v'], r['i']) for r in rest}),
          'arrays', len({(r['v'], r['array']) for r in rest}))


if __name__ == '__main__':
    main()
