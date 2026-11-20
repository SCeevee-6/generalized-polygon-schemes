"""Primitive intersection arrays from Brouwer's machine-readable tables (https://aeb.win.tue.nl/drg/drgtables.html).

Reads IA3P, IA4P, IA5P, IA6P, IA7P and IA8 ... IA13, keeps the arrays that are neither bipartite nor antipodal,
and writes primitive_arrays.json (b, c, v, d). These are the feasible arrays of BCN (1989) for at most 4096
vertices; the files carry no information on which arrays are known to exist or to be impossible.

Usage: python fetch_arrays.py
"""
import os, re, json
from fractions import Fraction as Fr
import requests

HERE = os.path.dirname(os.path.abspath(__file__))
FILES = ['IA3P', 'IA4P', 'IA5P', 'IA6P', 'IA7P', 'IA8', 'IA9', 'IA10', 'IA11', 'IA12', 'IA13']


def antipodal(b, c):
    d = len(b)
    return all(b[i] == c[d - 1 - i] for i in range(d) if i != d // 2)


def main():
    out = []; seen = 0
    for f in FILES:
        r = requests.get('https://aeb.win.tue.nl/drg/' + f, timeout=60)
        r.raise_for_status()
        for m in re.finditer(r'\(([\d,\s]+);([\d,\s]+)\)', r.text):
            b = [int(x) for x in m.group(1).split(',')]
            c = [int(x) for x in m.group(2).split(',')]
            seen += 1
            assert len(b) == len(c)
            k = [Fr(1)]
            for i in range(len(b)):
                k.append(k[-1] * b[i] / c[i])
            v = sum(k)
            assert v.denominator == 1
            d = len(b)
            a = [b[0] - b[i] - (c[i - 1] if i > 0 else 0) for i in range(d)] + [b[0] - c[-1]]
            if all(x == 0 for x in a) or antipodal(b, c):
                continue
            out.append({'file': f, 'b': b, 'c': c, 'v': int(v), 'd': d})
    json.dump(out, open(os.path.join(HERE, 'primitive_arrays.json'), 'w'), indent=0)
    print('arrays read', seen, 'primitive', len(out), 'distinct v', len({a['v'] for a in out}))


if __name__ == '__main__':
    main()
