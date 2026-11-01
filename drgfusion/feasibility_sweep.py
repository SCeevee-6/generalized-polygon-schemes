"""Krein conditions, absolute bound and the Delsarte-Goethals-Seidel s-distance bound
for every intersection array in Brouwer's machine-readable DRG tables.

Usage: python feasibility_sweep.py CACHE_DIR OUT_CSV
CACHE_DIR holds the table files (IA3P.txt, ..., IA4a_ext.txt, IA6b_ext.txt, ...).

Floating point (numpy) with tolerances; every flagged array is re-checked exactly
with sympy by exact_check.py.
"""
import sys, os, re, csv
from math import comb
import numpy as np

TOL = 1e-7


def parse(cache):
    arrs = {}
    for fn in sorted(os.listdir(cache)):
        if not fn.startswith('IA'):
            continue
        tab = fn.split('.')[0]
        for l in open(os.path.join(cache, fn)):
            m = re.match(r'^\s*\w*\(([^;]+);([^)]+)\)', l)
            if m:
                b = tuple(int(x) for x in m.group(1).split(','))
                c = tuple(int(x) for x in m.group(2).split(','))
                arrs.setdefault((b, c), []).append(tab)
    return arrs


def spectrum(b, c):
    d = len(b); k = b[0]
    bb = list(b) + [0]; cc = [0] + list(c)
    a = [k - bb[i] - cc[i] for i in range(d + 1)]
    kk = [1.0]
    for i in range(d):
        kk.append(kk[-1] * bb[i] / cc[i + 1])
    v = sum(kk)
    T = np.zeros((d + 1, d + 1))
    for i in range(d + 1):
        T[i, i] = a[i]
        if i < d:
            T[i, i + 1] = T[i + 1, i] = np.sqrt(bb[i] * cc[i + 1])
    th = np.sort(np.linalg.eigvalsh(T))[::-1]
    U = []
    for t in th:
        u = [1.0, t / k]
        for i in range(1, d):
            u.append(((t - a[i]) * u[i] - cc[i] * u[i - 1]) / bb[i])
        U.append(u)
    U = np.array(U); kk = np.array(kk)
    m = v / (U ** 2 @ kk)
    return v, kk, th, U, m


def analyse(b, c):
    v, kk, th, U, m = spectrum(b, c)
    d = len(b)
    res = {'v': round(v), 'krein': [], 'absolute': [], 'sdist': []}
    W = U * kk  # W[i,l] = k_l u_i(l)
    for i in range(1, d + 1):
        for j in range(i, d + 1):
            q = m[i] * m[j] / v * (W[i] * U[j]) @ U.T  # q[h] = q^h_ij
            scale = m[i] * m[j]
            neg = [h for h in range(d + 1) if q[h] < -TOL * scale]
            if neg:
                res['krein'].append((i, j, neg, [float(q[h]) for h in neg]))
            supp = [h for h in range(d + 1) if abs(q[h]) > TOL * scale]
            s = sum(m[h] for h in supp)
            bound = m[i] * (m[i] + 1) / 2 if i == j else m[i] * m[j]
            if s > bound + 1e-6:
                res['absolute'].append((i, j, supp, round(s), bound))
    for i in range(1, d + 1):
        vals = []
        for l in range(1, d + 1):
            if not any(abs(U[i][l] - x) < 1e-9 for x in vals):
                vals.append(U[i][l])
        if any(abs(x - 1) < 1e-9 for x in vals):
            continue  # representation not injective
        s = len(vals); mi = round(m[i])
        bnd = comb(mi + s - 1, s) + comb(mi + s - 2, s - 1)
        if v > bnd + 1e-6:
            res['sdist'].append((i, float(th[i]), mi, s, bnd))
    res['mult_integral'] = all(abs(x - round(x)) < 1e-6 for x in m)
    return res


def main(cache, out):
    arrs = parse(cache)
    rows = []
    for (b, c), tabs in sorted(arrs.items(), key=lambda x: (len(x[0][0]), x[0])):
        r = analyse(b, c)
        if r['krein'] or r['absolute'] or r['sdist'] or not r['mult_integral']:
            rows.append({'array': '{%s;%s}' % (','.join(map(str, b)), ','.join(map(str, c))),
                         'd': len(b), 'v': r['v'], 'tables': ' '.join(sorted(set(tabs))),
                         'mult_integral': r['mult_integral'],
                         'krein_fail': r['krein'], 'absolute_fail': r['absolute'], 'sdist_fail': r['sdist']})
    with open(out, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()) if rows else ['array'])
        w.writeheader(); w.writerows(rows)
    print(len(arrs), 'arrays;', len(rows), 'flagged')


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
