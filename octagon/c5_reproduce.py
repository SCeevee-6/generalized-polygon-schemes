"""Reproduce the C5-stable count 2,167,056 -> 251 -> 6 -> 1 with an independent method.
The rotation acts on edges of the base (no spanning-tree gauge); covers are built from the invariant subspaces."""
import sys, os, json, numpy as np
from collections import Counter
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import c5enum, recog
from base import *
cl, cF, cN, RNT, Fix, Nn, eF, sF, sN = c5enum.main()
tot, surv = c5enum.count(cF, cN, RNT)
print('labelled covers', 1 + 31*69905, 'surviving girth 8:', tot // 15 + (1 if (cF != 0).all() else 0))
def span4(f):
    fs = [f]
    for _ in range(3): fs.append(c5enum.compose(fs[-1], RNT))
    s = [0]
    for b in fs: s += [x ^ b for x in s]
    return fs, min(s[1:])
covers = {}
for f, g in surv:
    fs, key = span4(f); covers.setdefault((key, g), (fs, g))
def voltages(fs, g):
    vol = [0]*40
    for j, k in enumerate(CHORDS):
        x = 1 << j; xf = c5enum.apply(eF, x); xn = xf ^ x
        bits = [bin(g & sF(xf)).count('1') & 1] + [bin(fj & sN(xn)).count('1') & 1 for fj in fs]
        vol[k] = sum(b << i for i, b in enumerate(bits))
    return vol
hist = Counter(); kept = []
for key, (fs, g) in sorted(covers.items()):
    vol = voltages(fs, g); F = fibre_dist(cover(vol)); mc = mate_counts(vol, F)
    hist[(min(mc), max(mc))] += 1
    if not (max(mc) <= 1 and min(mc) == 0): kept.append(vol)
print('mate-count histogram', dict(hist))
for vol in kept:
    E = [(32*VI[a]+t, 32*VI[b]+(t ^ vol[k])) for k, (a, b) in enumerate(EDGES) for t in range(32)]
    r = recog.recognise(E)
    print('cover', r['ok'], r.get('why'), Counter(r['row']) if 'row' in r else '')
json.dump([[int(v) for v in vol] for vol in kept], open('c5_kept_covers.json', 'w'))
