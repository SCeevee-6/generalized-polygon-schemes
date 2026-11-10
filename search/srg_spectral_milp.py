"""SRG(v,k,lambda,mu) as a union of orbitals of a group whose orbital configuration is commutative.

For a commutative scheme with eigenmatrix P, a transposition-closed union S is an SRG(v,k,lambda,mu) iff
sum_{a in S} val_a = k and every non-principal eigenvalue theta_j(S) = sum_{c in S} P[j,c] is r or s, the two
restricted eigenvalues of the parameter set (then A_S^2 - (r+s)A_S + rs I is a multiple of J). This is linear in
binary x_a (units) and y_j (theta_j(S) = r): a MILP solved with HiGHS (scipy.optimize.milp). The eigenmatrix and
the MILP are floating point, so an infeasible answer is numerical; a solution is re-checked exactly with
srg_drg_search.distance_regular.

Usage: python srg_spectral_milp.py GENS_DIR JOBS_FILE PARAMS_CSV OUT_JSONL [time_limit]
"""
import sys, os, json, time, csv
import numpy as np
from scipy.optimize import milp, LinearConstraint, Bounds
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import prim_tools as pt
import srg_drg_search as sds
import coset_schemes as cs


def eig_units(p, units):
    n = int(p[:, :, 0].sum(axis=1).sum())
    P, Q, m = cs.eigenmatrix(np.ascontiguousarray(p.transpose(2, 0, 1)), n)  # expects pk[k,i,j] = p^k_ij
    th = np.array([[P[j, u].sum() for u in units] for j in range(P.shape[0])])
    assert np.abs(th.imag).max() < 1e-6
    return th.real, np.real(m)


def srg_eigs(v, k, l, mu):
    D = (l - mu) ** 2 + 4 * (k - mu)
    r = ((l - mu) + D ** 0.5) / 2; s = ((l - mu) - D ** 0.5) / 2
    f = ((v - 1) * (-s) - k) / (r - s)  # multiplicity of r
    return r, s, f


def solve(th, m, uv, v, k, l, mu, tl):
    U = th.shape[1]; R = th.shape[0] - 1
    r, s, f = srg_eigs(v, k, l, mu)
    nv = U + R
    A = []; lo = []; hi = []
    A.append(np.r_[uv, np.zeros(R)]); lo.append(k); hi.append(k)
    for j in range(1, R + 1):
        row = np.zeros(nv); row[:U] = th[j]; row[U + j - 1] = -(r - s)
        A.append(row); lo.append(s - 1e-6); hi.append(s + 1e-6)
    row = np.r_[np.zeros(U), m[1:]]; A.append(row); lo.append(f - 1e-3); hi.append(f + 1e-3)
    res = milp(c=np.zeros(nv), constraints=LinearConstraint(np.array(A), lo, hi), integrality=np.ones(nv),
               bounds=Bounds(0, 1), options={'time_limit': tl, 'presolve': True})
    return res


def main(gdir, jobs, pcsv, out, tl=600.0):
    params = {}
    for rr in csv.DictReader(open(pcsv)):
        v, k, l, mm = int(rr['v']), int(rr['k']), int(rr['l']), int(rr['m'])
        params.setdefault(v, set()).update({(k, l, mm), (v - k - 1, v - 2 - 2 * k + mm, v - 2 * k + l)})
    done = set()
    if os.path.exists(out):
        for l in open(out): r = json.loads(l); done.add((r['v'], r['i']))
    with open(out, 'a') as fo:
        for l in open(jobs):
            if not l.strip(): continue
            v, i = map(int, l.split())
            if (v, i) in done: continue
            t0 = time.time()
            p = pt.scheme(gdir, v, i)
            val, pair, units = sds.structure(p)
            if not np.array_equal(p, p.transpose(1, 0, 2)):
                fo.write(json.dumps({'v': v, 'i': i, 'error': 'not commutative'}) + '\n'); fo.flush(); continue
            uv = np.array([int(val[u].sum()) for u in units], dtype=float)
            th, m = eig_units(p, units)
            out_r = []
            for (k, lam, mu) in sorted(params.get(v, [])):
                reach = {0}
                for w in sorted(uv.astype(int)):
                    reach |= {a + w for a in reach if a + w <= k}
                if k not in reach:
                    out_r.append({'param': [v, k, lam, mu], 'status': 'k not a sum of unit valencies'}); continue
                res = solve(th, m, uv, v, k, lam, mu, tl)
                e = {'param': [v, k, lam, mu], 'status': {0: 'FOUND', 1: 'LIMIT', 2: 'INFEASIBLE'}.get(res.status, str(res.status))}
                if res.status == 0:
                    sel = [a for a in range(len(units)) if res.x[a] > 0.5]
                    S = sorted(c for a in sel for c in units[a])
                    e['S'] = S; e['check'] = str(sds.distance_regular(p, S))[:200]
                out_r.append(e)
            fo.write(json.dumps({'v': v, 'i': i, 'rank': int(p.shape[0]), 'units': len(units), 'results': out_r,
                                 'seconds': round(time.time() - t0, 1)}) + '\n'); fo.flush()


if __name__ == '__main__':
    main(*sys.argv[1:5], *(float(a) for a in sys.argv[5:6]))
