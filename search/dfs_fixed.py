"""Exact depth-first search for a union of orbitals S = S^T with prescribed local parameters.

SRG mode (v, k, lambda, mu): sum of valencies over S is k, c_j(S) = sum_{a,b in S} p^j_ab equals lambda for every
class j in S and mu for every class j != 0 outside S.
DRG mode (k, a_1, c_2): the same, except that c_j(S) is allowed to be 0 or c_2 outside S (vertices at distance 2
have c_2 common neighbours, vertices further apart have none). Candidates are checked afterwards.

c_j is monotone in S and so is the valency, so every partial assignment gives intervals for the final values and a
branch is cut as soon as one target falls outside its interval. Integer arithmetic only (int64); works for
non-commutative configurations as well. Used for the groups where the generic search in dfs_numba.py is too slow.

Usage: python dfs_fixed.py GENS_DIR V I PARAMS_CSV OUT_JSONL      (SRG mode, all parameter sets with that v)
"""
import sys, os, json, time, csv
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from numba import njit


@njit(cache=True)
def dfs_fixed(Rsym, Qd, ucl, uval, order, K, LAM, MU, drg, maxsol, maxnodes):
    U = Rsym.shape[0]; r = Rsym.shape[2]
    cS = np.zeros((U + 1, r), np.int64)
    cW = np.zeros((U + 1, r), np.int64)
    kS = np.zeros(U + 1, np.int64)
    kW = np.zeros(U + 1, np.int64)
    st = np.zeros((U + 1, U), np.uint8)    # 1 in, 2 out, 0 undecided
    state = np.zeros(U + 1, np.int64)
    for a in range(U):
        kW[0] += uval[a]
        for k in range(r):
            cW[0, k] += Qd[a, k]
            for b in range(a):
                cW[0, k] += Rsym[a, b, k]
    sols = np.zeros((maxsol, U), np.uint8)
    nsol = 0; nodes = 0; t = 0
    while t >= 0:
        if t == U:
            if nsol < maxsol:
                for a in range(U):
                    sols[nsol, a] = st[t, a]
            nsol += 1
            t -= 1
            state[t] += 1
            continue
        if state[t] >= 2:
            t -= 1
            if t >= 0:
                state[t] += 1
            continue
        nodes += 1
        if nodes > maxnodes:
            return sols[:min(nsol, maxsol)], nodes, -1
        x = order[t]
        for a in range(U):
            st[t + 1, a] = st[t, a]
        if state[t] == 0:
            st[t + 1, x] = 1
            kS[t + 1] = kS[t] + uval[x]; kW[t + 1] = kW[t]
            for k in range(r):
                d = Qd[x, k]
                for a in range(U):
                    if st[t, a] == 1:
                        d += Rsym[a, x, k]
                cS[t + 1, k] = cS[t, k] + d
                cW[t + 1, k] = cW[t, k]
        else:
            st[t + 1, x] = 2
            kS[t + 1] = kS[t]; kW[t + 1] = kW[t] - uval[x]
            for k in range(r):
                d = Qd[x, k]
                for a in range(U):
                    if st[t, a] != 2 and a != x:
                        d += Rsym[a, x, k]
                cW[t + 1, k] = cW[t, k] - d
                cS[t + 1, k] = cS[t, k]
        ok = kS[t + 1] <= K and kW[t + 1] >= K
        if ok:
            for a in range(U):
                sa = st[t + 1, a]
                for c in range(2):
                    k = ucl[a, c]
                    if k < 0:
                        continue
                    lo = cS[t + 1, k]; hi = cW[t + 1, k]
                    if sa == 1:
                        if lo > LAM or hi < LAM:
                            ok = False
                    elif sa == 2:
                        if drg:
                            if lo > MU or (lo > 0 and hi < MU):
                                ok = False
                        elif lo > MU or hi < MU:
                            ok = False
                    else:
                        # undecided: must be able to reach lambda (if taken) or mu / 0 (if left out)
                        can_in = lo <= LAM and hi >= LAM
                        can_out = (lo <= MU and hi >= MU) or (drg and lo == 0)
                        if not (can_in or can_out):
                            ok = False
                if not ok:
                    break
        if not ok:
            state[t] += 1
            continue
        t += 1
        state[t] = 0
    return sols[:min(nsol, maxsol)], nodes, nsol


def prepare(p):
    """p[i,j,k] = p^k_ij. Returns valencies, units, and the unit-level sums used by the search."""
    import srg_drg_search as sds
    val, pair, units = sds.structure(p)
    U = len(units); r = p.shape[0]
    Q = np.zeros((U, U, r), dtype=np.int64)
    for a, ua in enumerate(units):
        for b, ub in enumerate(units):
            Q[a, b] = p[np.ix_(ua, ub)].sum(axis=(0, 1))
    Rsym = Q + Q.transpose(1, 0, 2)
    Qd = np.array([Q[a, a] for a in range(U)])
    ucl = -np.ones((U, 2), dtype=np.int64)
    for a, u in enumerate(units):
        for c, k in enumerate(u):
            ucl[a, c] = k
    uval = np.array([int(val[u].sum()) for u in units], dtype=np.int64)
    return val, units, Rsym, Qd, ucl, uval


def search(p, K, LAM, MU, drg=False, maxsol=1000, maxnodes=10 ** 11, prep=None):
    val, units, Rsym, Qd, ucl, uval = prep if prep is not None else prepare(p)
    order = np.array(sorted(range(len(units)), key=lambda a: -uval[a]), dtype=np.int64)
    sols, nodes, nsol = dfs_fixed(Rsym, Qd, ucl, uval, order, K, LAM, MU, drg, maxsol, maxnodes)
    S = [sorted(int(i) for a in range(len(units)) if s[a] == 1 for i in units[a]) for s in sols]
    return S, int(nodes), int(nsol)


def main(gdir, v, i, pcsv, out):
    import prim_tools as pt
    v, i = int(v), int(i)
    params = set()
    for rr in csv.DictReader(open(pcsv)):
        if int(rr['v']) == v:
            k, l, m = int(rr['k']), int(rr['l']), int(rr['m'])
            params |= {(k, l, m), (v - k - 1, v - 2 - 2 * k + m, v - 2 * k + l)}
    t0 = time.time()
    p = pt.scheme(gdir, v, i)
    prep = prepare(p)
    res = []
    for (k, l, m) in sorted(params):
        t1 = time.time()
        S, nodes, nsol = search(p, k, l, m, prep=prep)
        res.append({'param': [v, k, l, m], 'status': 'LIMIT' if nsol < 0 else ('FOUND' if nsol else 'NONE'),
                    'solutions': S[:5], 'nodes': nodes, 'seconds': round(time.time() - t1, 1)})
        print(v, i, (k, l, m), res[-1]['status'], nodes, res[-1]['seconds'], flush=True)
    with open(out, 'a') as fo:
        fo.write(json.dumps({'v': v, 'i': i, 'rank': int(p.shape[0]), 'units': len(prep[1]),
                             'commutative': bool(np.array_equal(p, p.transpose(1, 0, 2))), 'results': res,
                             'seconds': round(time.time() - t0, 1)}) + '\n')


if __name__ == '__main__':
    main(*sys.argv[1:6])
