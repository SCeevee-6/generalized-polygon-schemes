"""Exact SRG test for primitive affine groups, replacing the floating-point eigenvalue MILP.

G = N:G_0 primitive affine of degree v = p^n, N elementary abelian and regular. Every G-invariant graph is
Cay(N, D) with D = -D a union of G_0-orbits on N minus 0. Its eigenvalue at the character chi_y(x) = zeta^(y.x) is
sum_c N_c(y) zeta^c with N_c(y) = #{x in D : y.x = c}. This is a rational integer iff N_1 = ... = N_{p-1}, and then
it equals N_0 - N_1. So Cay(N, D) is an SRG(v,k,lambda,mu) with integral restricted eigenvalues r > s iff |D| = k
and for every y != 0:  N_1(y) = N_c(y) for all c, and N_0(y) - N_1(y) is r or s. All of this is linear in the 0/1
indicators of the symmetric units with integer coefficients, so CP-SAT decides it exactly (INFEASIBLE is a proof).

Input: affine/aff_V_I.txt from export_affine.g (one line "PT point suborbit [exponents]" per element of N).
Usage: python affine_exact.py V I [V I ...]      (appends to affine_exact_results.jsonl)
"""
from ortools.sat.python import cp_model
import os, sys, csv, json, time, re, math
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))


def load(v, i):
    lab, co = {}, {}
    for l in open(os.path.join(HERE, 'affine', f'aff_{v}_{i}.txt')):
        if l.startswith('PT '):
            m = re.match(r'PT (\d+) (\d+) \[(.*)\]', l.strip())
            x = int(m.group(1)); lab[x] = int(m.group(2))
            co[x] = tuple(int(a) for a in m.group(3).replace(' ', '').split(','))
    assert len(lab) == v
    n = len(co[1]); p = round(v ** (1.0 / n))
    assert p ** n == v and all(c == 0 for c in co[1])
    return p, n, lab, co


def units(p, n, lab, co):
    """symmetric units: suborbit D_a merged with -D_a; returns list of arrays of vectors."""
    pts = sorted(lab)
    X = np.array([co[x] for x in pts], dtype=np.int64)
    L = np.array([lab[x] for x in pts])
    idx = {tuple(r): j for j, r in enumerate(X)}
    neg = np.array([idx[tuple((-r) % p)] for r in X])
    par = {}
    for a in set(L.tolist()):
        par[a] = a
    def find(a):
        while par[a] != a:
            a = par[a]
        return a
    for j in range(len(pts)):
        a, b = find(int(L[j])), find(int(L[neg[j]]))
        if a != b:
            par[max(a, b)] = min(a, b)
    zero = int(L[idx[(0,) * n]])
    groups = {}
    for j in range(len(pts)):
        if int(L[j]) != zero:
            groups.setdefault(find(int(L[j])), []).append(j)
    return X, [np.array(g) for _, g in sorted(groups.items())]


def char_rows(p, n, X, U):
    """for every y != 0: E[y,u] = N_0 - N_1 and, for p > 3, the differences N_1 - N_c (c = 2..(p-1)/2).
    Rows are deduplicated."""
    Y = np.array(np.unravel_index(np.arange(1, p ** n), (p,) * n)).T.astype(np.int64)
    D = (Y @ X.T) % p                                   # D[y, x]
    half = (p - 1) // 2
    cnt = np.zeros((len(Y), len(U), p), dtype=np.int64)
    for a, g in enumerate(U):
        sub = D[:, g]
        for c in range(p):
            cnt[:, a, c] = (sub == c).sum(axis=1)
    if p > 2:
        for c in range(1, half + 1):
            assert np.array_equal(cnt[:, :, c], cnt[:, :, p - c])
    N1 = cnt[:, :, 1] if p > 2 else cnt[:, :, 1]
    E = cnt[:, :, 0] - N1
    F = [cnt[:, :, 1] - cnt[:, :, c] for c in range(2, half + 1)]
    full = np.concatenate([E] + F, axis=1)
    _, keep = np.unique(full, axis=0, return_index=True)
    keep = np.sort(keep)
    return E[keep], [f[keep] for f in F]


def srg_eigs(v, k, l, mu):
    D = (l - mu) ** 2 + 4 * (k - mu)
    q = math.isqrt(D)
    if q * q != D or (l - mu + q) % 2:
        return None
    return (l - mu + q) // 2, (l - mu - q) // 2


def decide(E, F, uv, v, k, lam, mu, tl):
    rs = srg_eigs(v, k, lam, mu)
    if rs is None:
        return 'NONINTEGRAL', None
    r, s = rs
    M = cp_model.CpModel()
    x = [M.NewBoolVar(f'x{a}') for a in range(len(uv))]
    M.Add(sum(int(uv[a]) * x[a] for a in range(len(uv))) == k)
    for j in range(E.shape[0]):
        z = M.NewBoolVar(f'z{j}')
        M.Add(sum(int(E[j, a]) * x[a] for a in range(len(uv))) == s + (r - s) * z)
        for f in F:
            if f[j].any():
                M.Add(sum(int(f[j, a]) * x[a] for a in range(len(uv))) == 0)
    S = cp_model.CpSolver()
    S.parameters.max_time_in_seconds = tl
    S.parameters.num_search_workers = 8
    st = S.Solve(M)
    name = {cp_model.INFEASIBLE: 'INFEASIBLE', cp_model.OPTIMAL: 'FOUND', cp_model.FEASIBLE: 'FOUND'}.get(st, 'LIMIT')
    sel = [a for a in range(len(uv)) if S.Value(x[a])] if name == 'FOUND' else None
    return name, sel


def main(args, tl=1800.0):
    params, ctrl = {}, {}
    for rr in csv.DictReader(open(os.path.join(HERE, 'open_srg.csv'))):
        v, k, l, mm = int(rr['v']), int(rr['k']), int(rr['l']), int(rr['m'])
        params.setdefault(v, set()).update({(k, l, mm), (v - k - 1, v - 2 - 2 * k + mm, v - 2 * k + l)})
    cf = os.path.join(HERE, 'affine', 'controls.csv')
    if os.path.exists(cf):
        for rr in csv.DictReader(open(cf)):
            ctrl.setdefault((int(rr['v']), int(rr['i'])), []).append((int(rr['k']), int(rr['l']), int(rr['m'])))
    out = os.path.join(HERE, 'affine_exact_results.jsonl')
    with open(out, 'a') as fo:
        for v, i in zip(args[0::2], args[1::2]):
            v, i = int(v), int(i)
            t0 = time.time()
            p, n, lab, co = load(v, i)
            X, U = units(p, n, lab, co)
            uv = [len(g) for g in U]
            E, F = char_rows(p, n, X, U)
            res = []
            for kind, plist in (('open', sorted(params.get(v, []))), ('control', ctrl.get((v, i), []))):
                for (k, lam, mu) in plist:
                    reach = {0}
                    for w in sorted(uv):
                        reach |= {a + w for a in reach if a + w <= k}
                    if k not in reach:
                        res.append({'param': [v, k, lam, mu], 'kind': kind, 'status': 'k not a sum of unit sizes'}); continue
                    st, sel = decide(E, F, uv, v, k, lam, mu, tl)
                    e = {'param': [v, k, lam, mu], 'kind': kind, 'status': st}
                    if sel is not None:
                        e['units'] = sel
                    res.append(e)
            rec = {'v': v, 'i': i, 'p': p, 'n': n, 'suborbits': len(set(lab.values())), 'units': len(U),
                   'distinct_character_rows': int(E.shape[0]), 'results': res, 'seconds': round(time.time() - t0, 1)}
            fo.write(json.dumps(rec) + '\n'); fo.flush()
            from collections import Counter
            print(v, i, 'units', len(U), 'rows', E.shape[0],
                  dict(Counter((r_['kind'], r_['status']) for r_ in res)), rec['seconds'], flush=True)


if __name__ == '__main__':
    main(sys.argv[1:])
