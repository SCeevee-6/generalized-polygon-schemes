"""Exhaustive depth-first search for unions of orbitals S (S = S^T) that satisfy the necessary conditions
for a connected distance-regular graph of diameter >= 2:
  (a) c_k(S) = sum_{i,j in S} p^k_ij is the same number a_1 for every class k in S;
  (b) for every class k outside S and k != 0, c_k(S) is 0 or one common number c_2 >= 1, and c_2 occurs.
c_k is monotone in S, so with S' the units chosen so far and W = S' plus the undecided units, the final
value lies in [c_k(S'), c_k(W)]; a branch is cut when the intervals for a_1 or for c_2 become empty.
Leaves are checked exactly with srg_drg_search.distance_regular. Integer arithmetic (int64).

Usage: python dfs_numba.py <intersection_numbers.csv> <out.json> [order] [prefix]
order: 'desc' (default, by unit valency), 'desc2' (same, ties reversed) or 'asc'.
prefix: optional string of 0/1 fixing the first decisions (for splitting into shards).
"""
import sys, os, json, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import srg_drg_search as sds
from numba import njit


@njit(cache=True)
def dfs(Q, Rsym, Qd, ucl, order, prefix, maxsol):
    U = Q.shape[0]; r = Q.shape[2]
    cS = np.zeros((U + 1, r), np.int64)
    cW = np.zeros((U + 1, r), np.int64)
    inS = np.zeros((U + 1, U), np.uint8)   # 1 included, 2 excluded, 0 undecided
    state = np.zeros(U + 1, np.int64)
    for a in range(U):
        for b in range(U):
            for k in range(r):
                cW[0, k] += Q[a, b, k]
    sols = np.zeros((maxsol, U), np.uint8)
    nsol = 0
    nodes = 0
    t = 0
    while t >= 0:
        if t == U or state[t] == 2:
            # leaf or both branches done
            if t == U:
                # final check
                S_any = False; S_all = True
                for a in range(U):
                    if inS[t, a] == 1:
                        S_any = True
                    else:
                        S_all = False
                if S_any and not S_all:
                    a1 = -1; ok = True
                    for a in range(U):
                        if inS[t, a] == 1:
                            for c in range(2):
                                k = ucl[a, c]
                                if k >= 0:
                                    if a1 < 0:
                                        a1 = cS[t, k]
                                    elif cS[t, k] != a1:
                                        ok = False
                    c2 = 0
                    if ok:
                        for a in range(U):
                            if inS[t, a] == 2:
                                for c in range(2):
                                    k = ucl[a, c]
                                    if k >= 0 and cS[t, k] > 0:
                                        if c2 == 0:
                                            c2 = cS[t, k]
                                        elif cS[t, k] != c2:
                                            ok = False
                    if ok and c2 > 0 and nsol < maxsol:
                        for a in range(U):
                            sols[nsol, a] = inS[t, a]
                        nsol += 1
            t -= 1
            if t >= 0:
                state[t] += 1
            continue
        x = order[t]
        if t < prefix.shape[0] and prefix[t] >= 0:
            # forced decision: branch 0 = include (prefix 1), branch 1 = exclude (prefix 0)
            allowed = 0 if prefix[t] == 1 else 1
            if state[t] < allowed:
                state[t] = allowed
            elif state[t] > allowed:
                state[t] = 2
                continue
        if state[t] >= 2:
            continue
        nodes += 1
        # build child t+1
        for a in range(U):
            inS[t + 1, a] = inS[t, a]
        if state[t] == 0:
            inS[t + 1, x] = 1
            for k in range(r):
                s = Qd[x, k]
                for a in range(U):
                    if inS[t, a] == 1:
                        s += Rsym[a, x, k]
                cS[t + 1, k] = cS[t, k] + s
                cW[t + 1, k] = cW[t, k]
        else:
            inS[t + 1, x] = 2
            for k in range(r):
                s = 0
                for a in range(U):
                    if inS[t, a] != 2:
                        s += Rsym[a, x, k]
                cW[t + 1, k] = cW[t, k] - s + Qd[x, k]
                cS[t + 1, k] = cS[t, k]
        # bounds
        lo1 = 0; hi1 = 1 << 60; lo2 = 1; hi2 = 1 << 60
        for a in range(U):
            st = inS[t + 1, a]
            if st == 0:
                continue
            for c in range(2):
                k = ucl[a, c]
                if k < 0:
                    continue
                if st == 1:
                    if cS[t + 1, k] > lo1: lo1 = cS[t + 1, k]
                    if cW[t + 1, k] < hi1: hi1 = cW[t + 1, k]
                elif cS[t + 1, k] > 0:
                    if cS[t + 1, k] > lo2: lo2 = cS[t + 1, k]
                    if cW[t + 1, k] < hi2: hi2 = cW[t + 1, k]
        if lo1 > hi1 or lo2 > hi2 or cW[t + 1, 0] == 0:
            state[t] += 1
            continue
        t += 1
        state[t] = 0
    return sols[:nsol], nodes


def prepare(p):
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
    return val, units, Q, Rsym, Qd, ucl


def run(inp, out, order='desc', prefix=''):
    t0 = time.time()
    p = sds.load_p(inp)
    val, units, Q, Rsym, Qd, ucl = prepare(p)
    uv = [int(val[u].sum()) for u in units]
    idx = sorted(range(len(units)), key=lambda a: (uv[a] if order == 'asc' else -uv[a], -a if order == 'desc2' else a))
    pre = np.array([int(ch) for ch in prefix] if prefix else [], dtype=np.int64)
    sols, nodes = dfs(Q, Rsym, Qd, ucl, np.array(idx, dtype=np.int64), pre, 100000)
    n = int(val.sum())
    closed = sds.closed_subsets(p)
    res = []
    for s in sols:
        S = sorted(i for a in range(len(units)) if s[a] == 1 for i in units[a])
        comp = tuple(sorted(set(range(p.shape[0])) - set(S)))
        e = {'S': S, 'k': int(val[S].sum()), 'complete_multipartite': comp in closed}
        d = sds.distance_regular(p, S)
        e['drg'] = d if isinstance(d, (str, type(None))) else {kk: d[kk] for kk in ('diameter', 'b', 'c')}
        res.append(e)
    info = {'input': inp, 'n': n, 'rank': int(p.shape[0]), 'units': len(units), 'order': order, 'prefix': prefix,
            'nodes': int(nodes), 'candidates': res, 'seconds': round(time.time() - t0, 1)}
    with open(out, 'w') as f:
        json.dump(info, f, indent=1)
    drg = [(e['k'], e['drg']) for e in res if isinstance(e['drg'], dict) and not e['complete_multipartite']]
    print(f"n={n} rank={p.shape[0]} units={len(units)} nodes={nodes} candidates={len(res)} "
          f"multipartite={sum(e['complete_multipartite'] for e in res)} other_drg={drg} t={info['seconds']}s")


if __name__ == '__main__':
    run(sys.argv[1], sys.argv[2], sys.argv[3] if len(sys.argv) > 3 else 'desc', sys.argv[4] if len(sys.argv) > 4 else '')
