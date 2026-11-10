"""SRG(v,k,lambda,mu) as a union of orbitals, for non-commutative orbital configurations.

Exact part (same as srg_fixed_cpsat.py): binary x_a for the transposition-closed units, y_ab = x_a x_b
(linearised), and for every unit u the walk count e_u = sum_{i,j in S} p^h_ij (h in u) equals lambda if u is in S
and mu if not. Extra cuts: every one-dimensional representation phi of the orbital algebra (a constituent of
multiplicity one in the permutation character) gives an eigenvalue phi(A_S) = sum_a x_a phi(A_a) of the graph,
which must be r or s. These cuts are necessary conditions, so the MILP is still exact as a search; they are
computed in floating point (labelled numerical), and every solution is re-checked exactly.

Usage: python srg_hybrid_milp.py GENS_DIR JOBS_FILE PARAMS_CSV OUT_JSONL [time_limit] [nocomp]
('nocomp': test only the listed parameter sets, not also their complements)
"""
import sys, os, json, time, csv
import numpy as np
from scipy.optimize import milp, LinearConstraint, Bounds
from scipy.sparse import lil_matrix
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import prim_tools as pt
import srg_drg_search as sds


def one_dim_reps(p, units):
    """Values phi(A_a) on units for the non-trivial 1-dim representations of the algebra."""
    r = p.shape[0]
    B = [p[l].T.astype(float) for l in range(r)]      # B_l[h,j] = p^h_lj : A_l A_j = sum_h p^h_lj A_h
    rng = np.random.default_rng(7)
    C = sum(c * b for c, b in zip(rng.normal(size=r), B))
    w, V = np.linalg.eig(C.T)                          # left eigenvectors of C
    reps = []
    for idx in range(r):
        if np.sum(np.abs(w - w[idx]) < 1e-6) != 1:
            continue
        vv = V[:, idx]
        vals = []
        ok = True
        for l in range(r):
            u = B[l].T @ vv
            piv = np.argmax(np.abs(vv)); lam = u[piv] / vv[piv]
            if not np.allclose(u, lam * vv, atol=1e-6 * max(1, abs(lam))):
                ok = False; break
            vals.append(lam)
        if ok:
            vals = np.array(vals)
            if abs(vals.imag).max() < 1e-6:
                reps.append(vals.real)
    val = p[:, :, 0].sum(axis=1)
    out = []
    for f in reps:
        if np.allclose(f, val):
            continue  # trivial representation
        out.append(np.array([f[u].sum() for u in units]))
    # drop duplicates
    uniq = []
    for f in out:
        if not any(np.allclose(f, g) for g in uniq):
            uniq.append(f)
    return uniq


def solve(p, units, val, phis, k, lam, mu, tl):
    U = len(units); r = p.shape[0]
    uv = [int(val[u].sum()) for u in units]
    pairs = [(a, b) for a in range(U) for b in range(a + 1, U)]
    nx, ny, nz = U, len(pairs), len(phis)
    N = nx + ny + nz
    D = (lam - mu) ** 2 + 4 * (k - mu); rr = ((lam - mu) + D ** 0.5) / 2; ss = ((lam - mu) - D ** 0.5) / 2
    Qs = {}
    for a in range(U):
        for b in range(U):
            Qs[a, b] = p[np.ix_(units[a], units[b])].sum(axis=(0, 1))
    rows = []; lo = []; hi = []
    def add(coefs, l_, h_):
        rows.append(coefs); lo.append(l_); hi.append(h_)
    add({a: uv[a] for a in range(U)}, k, k)
    for t, (a, b) in enumerate(pairs):
        add({nx + t: 1, a: -1}, -np.inf, 0); add({nx + t: 1, b: -1}, -np.inf, 0)
        add({nx + t: 1, a: -1, b: -1}, -1, np.inf)
    M = k * k + 1
    for u in range(U):
        h = units[u][0]
        e = {}
        for a in range(U):
            if Qs[a, a][h]: e[a] = e.get(a, 0) + int(Qs[a, a][h])
        for t, (a, b) in enumerate(pairs):
            q = int(Qs[a, b][h] + Qs[b, a][h])
            if q: e[nx + t] = q
        # e - lam <= M(1-x_u), e - lam >= -M(1-x_u); e - mu <= M x_u, e - mu >= -M x_u
        c = dict(e); c[u] = c.get(u, 0) + M; add(c, -np.inf, lam + M)
        c = dict(e); c[u] = c.get(u, 0) - M; add(c, lam - M, np.inf)
        c = dict(e); c[u] = c.get(u, 0) - M; add(c, -np.inf, mu)
        c = dict(e); c[u] = c.get(u, 0) + M; add(c, mu, np.inf)
    for j, f in enumerate(phis):
        c = {a: f[a] for a in range(U) if abs(f[a]) > 1e-12}; c[nx + ny + j] = -(rr - ss)
        add(c, ss - 1e-6, ss + 1e-6)
    A = lil_matrix((len(rows), N))
    for i, c in enumerate(rows):
        for j, vv in c.items():
            A[i, j] = vv
    res = milp(c=np.zeros(N), constraints=LinearConstraint(A.tocsr(), lo, hi), integrality=np.ones(N),
               bounds=Bounds(0, 1), options={'time_limit': tl})
    return res


def main(gdir, jobs, pcsv, out, tl=3600.0, nocomp=False):
    params = {}
    for row in csv.DictReader(open(pcsv)):
        v, k, l, mm = int(row['v']), int(row['k']), int(row['l']), int(row['m'])
        params.setdefault(v, set()).add((k, l, mm))
        if not nocomp:
            params[v].add((v - k - 1, v - 2 - 2 * k + mm, v - 2 * k + l))
    with open(out, 'a') as fo:
        for line in open(jobs):
            if not line.strip(): continue
            v, i = map(int, line.split())
            t0 = time.time()
            p = pt.scheme(gdir, v, i)
            val, pair, units = sds.structure(p)
            phis = one_dim_reps(p, units)
            uv = sorted(int(val[u].sum()) for u in units)
            res_all = []
            for (k, lam, mu) in sorted(params.get(v, [])):
                reach = {0}
                for w in uv:
                    reach |= {a + w for a in reach if a + w <= k}
                if k not in reach:
                    res_all.append({'param': [v, k, lam, mu], 'status': 'k not a sum of unit valencies'}); continue
                t1 = time.time()
                res = solve(p, units, val, phis, k, lam, mu, tl)
                e = {'param': [v, k, lam, mu], 'status': {0: 'FOUND', 1: 'LIMIT', 2: 'INFEASIBLE'}.get(res.status, str(res.status)),
                     'seconds': round(time.time() - t1, 1)}
                if res.status == 0:
                    S = sorted(c for a in range(len(units)) if res.x[a] > 0.5 for c in units[a])
                    e['S'] = S; e['check'] = str(sds.distance_regular(p, S))[:200]
                res_all.append(e)
                fo.write(json.dumps({'v': v, 'i': i, 'partial': e}) + '\n'); fo.flush()
            fo.write(json.dumps({'v': v, 'i': i, 'rank': int(p.shape[0]), 'units': len(units), 'one_dim_reps': len(phis),
                                 'results': res_all, 'seconds': round(time.time() - t0, 1)}) + '\n'); fo.flush()


if __name__ == '__main__':
    main(*sys.argv[1:5], *(float(a) for a in sys.argv[5:6]), *([True] if sys.argv[6:7] == ['nocomp'] else []))
