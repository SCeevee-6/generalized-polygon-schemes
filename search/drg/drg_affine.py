"""Vertex-primitive distance-regular graphs, affine groups of large rank (step 2b).

For G = N:G_0 of degree p^n with p = 2 or 3 every invariant graph is Cay(N, D) and its eigenvalues are rational
integers, sum_{x in D} chi_y(x) = N_0(y) - N_1(y) (see ../primitive/affine_exact.py). A graph with intersection
array {b; c} has the spectrum determined by the array, so: if some eigenvalue of the array is not an integer there
is no such graph; otherwise D must have size b_0, every character sum must be one of the eigenvalues, and each
eigenvalue must occur with its multiplicity. Solutions of this integer model are enumerated with CP-SAT (no-good
cuts), each one is built as a Cayley graph and tested for the array directly; the enumeration stops at the first
graph with the right array.

Usage: python drg_affine.py V I [V I ...]      (appends to drg_affine_results.jsonl)
"""
from ortools.sat.python import cp_model
import os, sys, csv, json, time
import numpy as np
import sympy as sp
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'primitive'))
import affine_exact as ae


def spectrum(b, c):
    d = len(b)
    a = [b[0] - b[i] - (c[i - 1] if i > 0 else 0) for i in range(d)] + [b[0] - c[-1]]
    X = sp.Symbol('X')
    u = [sp.Integer(1), X / b[0]]
    for j in range(1, d):
        u.append(sp.expand(((X - a[j]) * u[j] - c[j - 1] * u[j - 1]) / b[j]))
    mp = sp.Poly(sp.expand((X - a[d]) * u[d] - c[d - 1] * u[d - 1]), X)
    k = [sp.Integer(1)]
    for j in range(d):
        k.append(k[-1] * sp.Rational(b[j], c[j]))
    rts = sp.roots(mp, X)
    if sum(rts.values()) != d + 1 or not all(r.is_Integer for r in rts):
        return None
    out = {}
    for r in rts:
        m = sum(k) / sum(k[j] * u[j].subs(X, r) ** 2 for j in range(d + 1))
        if not m.is_Integer:
            return None
        out[int(r)] = int(m)
    return out


def weighted_rows(p, n, X, U):
    Y = np.array(np.unravel_index(np.arange(1, p ** n), (p,) * n)).T.astype(np.int64)
    D = (Y @ X.T) % p
    E = np.zeros((len(Y), len(U)), dtype=np.int64)
    for a, g in enumerate(U):
        sub = D[:, g]
        E[:, a] = (sub == 0).sum(axis=1) - (sub == 1).sum(axis=1)
    rows, cnt = np.unique(E, axis=0, return_counts=True)
    return rows, cnt


def is_drg(p, n, X, D, b, c):
    """distance partition from 0 in Cay(N, D) and the intersection numbers along it (the graph is vertex-transitive)."""
    v = len(X)
    code = lambda V: (V * (p ** np.arange(n))).sum(axis=1)
    idx = np.zeros(p ** n, dtype=np.int64); idx[code(X)] = np.arange(v)
    nb = np.stack([idx[code((X + d_) % p)] for d_ in D], axis=1)
    dist = -np.ones(v, dtype=np.int64); zero = idx[0]; dist[zero] = 0
    frontier = [zero]; h = 0
    while frontier:
        nxt = np.unique(nb[frontier].ravel()); nxt = nxt[dist[nxt] < 0]
        if len(nxt) == 0:
            break
        h += 1; dist[nxt] = h; frontier = list(nxt)
    d = len(b)
    if dist.max() != d or (dist < 0).any():
        return False
    dn = dist[nb]
    for j in range(d + 1):
        w = np.where(dist == j)[0]
        cj = (dn[w] == j - 1).sum(axis=1); bj = (dn[w] == j + 1).sum(axis=1)
        if len(set(cj.tolist())) > 1 or len(set(bj.tolist())) > 1:
            return False
        if j > 0 and cj[0] != c[j - 1]:
            return False
        if j < d and bj[0] != b[j]:
            return False
    return True


def main(args, tl=1500.0):
    todo = {}
    for r in csv.DictReader(open(os.path.join(HERE, 'survivors.csv'))):
        if not r['hamming']:
            todo.setdefault((int(r['v']), int(r['i'])), []).append(r['array'])
    with open(os.path.join(HERE, 'drg_affine_results.jsonl'), 'a') as fo:
        for v, i in zip(args[0::2], args[1::2]):
            v, i = int(v), int(i)
            p, n, lab, co = ae.load(v, i)
            assert p in (2, 3)
            X, U = ae.units(p, n, lab, co)
            uv = [len(g) for g in U]
            rows, cnt = weighted_rows(p, n, X, U)
            res = []
            for arr in todo.get((v, i), []):
                t0 = time.time()
                b, c = [[int(x) for x in part.split(',')] for part in arr.strip('{}').split(';')]
                spec = spectrum(b, c)
                if spec is None:
                    res.append({'array': arr, 'status': 'spectrum not integral'}); continue
                th = sorted(t_ for t_ in spec if t_ != b[0])
                M = cp_model.CpModel()
                x = [M.NewBoolVar(f'x{a}') for a in range(len(U))]
                M.Add(sum(uv[a] * x[a] for a in range(len(U))) == b[0])
                z = {}
                for j in range(rows.shape[0]):
                    zs = [M.NewBoolVar(f'z{j}_{q}') for q in range(len(th))]
                    M.AddExactlyOne(zs)
                    M.Add(sum(int(rows[j, a]) * x[a] for a in range(len(U))) == sum(th[q] * zs[q] for q in range(len(th))))
                    z[j] = zs
                for q, t_ in enumerate(th):
                    M.Add(sum(int(cnt[j]) * z[j][q] for j in range(rows.shape[0])) == spec[t_])
                S = cp_model.CpSolver(); S.parameters.num_search_workers = 8
                n_cand = 0; hit = None; status = 'NONE'
                while True:
                    S.parameters.max_time_in_seconds = max(1.0, tl - (time.time() - t0))
                    st = S.Solve(M)
                    if st == cp_model.INFEASIBLE:
                        break
                    if st not in (cp_model.OPTIMAL, cp_model.FEASIBLE):
                        status = 'LIMIT'; break
                    sel = [a for a in range(len(U)) if S.Value(x[a])]
                    n_cand += 1
                    D = np.concatenate([X[U[a]] for a in sel])
                    if is_drg(p, n, X, D, b, c):
                        hit = sel; status = 'REALISED'; break
                    M.AddBoolOr([x[a].Not() for a in sel] + [x[a] for a in range(len(U)) if a not in sel])
                res.append({'array': arr, 'status': status, 'cospectral_candidates_tested': n_cand, 'hit': hit,
                            'seconds': round(time.time() - t0, 1)})
                print(v, i, arr, status, n_cand, res[-1]['seconds'], flush=True)
            fo.write(json.dumps({'v': v, 'i': i, 'units': len(U), 'results': res}) + '\n'); fo.flush()


if __name__ == '__main__':
    main(sys.argv[1:])
