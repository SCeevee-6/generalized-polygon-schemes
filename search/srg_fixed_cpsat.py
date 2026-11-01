"""Is there an SRG(v,k,lambda,mu) that is a union of orbitals of a given transitive group?

Usage: python srg_fixed_cpsat.py GENS_DIR JOBS_FILE PARAMS_CSV OUT_JSONL [time_limit_per_param]
JOBS_FILE lines "v i"; PARAMS_CSV has columns v,k,l,m (the parameter sets to test; both a set and its
complement are tried). For a union S of transposition-closed units with A_S of valency k:
  S gives an SRG(v,k,lambda,mu)  <=>  sum_{i,j in S} p^h_ij = lambda for h in S and = mu for h not in S, h != 0.
With k, lambda, mu fixed these are linear constraints in x_a and y_ab = x_a x_b, solved by CP-SAT (exact,
integers). A solution is re-checked with srg_drg_search.distance_regular.
"""
import sys, os, json, time, csv, inspect
from ortools.sat.python import cp_model  # import before numpy-heavy modules (Windows DLL issue)
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, '..', 'tits1755'))
import srg_drg_search as sds
import coset_schemes as cs

# orbital label matrix with int16 labels (coset_schemes uses int8, enough only for rank <= 127)
_src = inspect.getsource(cs.orbitals).replace('dtype=np.int8', 'dtype=np.int16')
_ns = dict(cs.__dict__); exec(_src, _ns); orbitals16 = _ns['orbitals']


def scheme(gdir, v, i):
    gens = [np.array([int(x) - 1 for x in l.split(',') if x], dtype=np.int64)
            for l in open(os.path.join(gdir, f'{v}_{i}.txt')).read().split('\n') if l.strip()]
    L, sizes, reps, _ = orbitals16(gens, v)
    pk = cs.intersection_numbers(L, reps)
    return np.ascontiguousarray(pk.transpose(1, 2, 0))


def solve(p, units, val, k, lam, mu, tl):
    U = len(units)
    uv = [int(val[u].sum()) for u in units]
    m = cp_model.CpModel()
    x = [m.NewBoolVar(f'x{a}') for a in range(U)]
    m.Add(sum(uv[a] * x[a] for a in range(U)) == k)
    y = {}
    for a in range(U):
        for b in range(a + 1, U):
            t = m.NewBoolVar(f'y{a}_{b}')
            m.AddImplication(t, x[a]); m.AddImplication(t, x[b]); m.AddBoolOr([x[a].Not(), x[b].Not(), t])
            y[a, b] = t
    for u in range(U):
        h = units[u][0]
        terms = []
        for a in range(U):
            q = int(p[np.ix_(units[a], units[a])].sum(axis=(0, 1))[h])
            if q: terms.append(q * x[a])
        for (a, b), t in y.items():
            q = int(p[np.ix_(units[a], units[b])].sum(axis=(0, 1))[h] + p[np.ix_(units[b], units[a])].sum(axis=(0, 1))[h])
            if q: terms.append(q * t)
        e = sum(terms) if terms else 0
        m.Add(e == lam).OnlyEnforceIf(x[u])
        m.Add(e == mu).OnlyEnforceIf(x[u].Not())
    s = cp_model.CpSolver(); s.parameters.num_workers = 1; s.parameters.max_time_in_seconds = tl
    st = s.Solve(m)
    name = s.StatusName(st)
    sol = [a for a in range(U) if s.Value(x[a])] if st in (cp_model.OPTIMAL, cp_model.FEASIBLE) else None
    return name, sol, round(s.WallTime(), 2)


def main(gdir, jobs, pcsv, out, tl=600.0):
    params = {}
    for r in csv.DictReader(open(pcsv)):
        v, k, l, mm = int(r['v']), int(r['k']), int(r['l']), int(r['m'])
        params.setdefault(v, set()).update({(k, l, mm), (v - k - 1, v - 2 - 2 * k + mm, v - 2 * k + l)})
    done = set()
    if os.path.exists(out):
        for l in open(out):
            r = json.loads(l); done.add((r['v'], r['i']))
    with open(out, 'a') as f:
        for l in open(jobs):
            if not l.strip(): continue
            v, i = map(int, l.split())
            if (v, i) in done or v not in params: continue
            t0 = time.time()
            p = scheme(gdir, v, i)
            val, pair, units = sds.structure(p)
            uvals = sorted(int(val[u].sum()) for u in units)
            res = []
            for (k, lam, mu) in sorted(params[v]):
                # subset-sum pre-check on k
                reach = {0}
                for w in uvals:
                    reach |= {a + w for a in reach if a + w <= k}
                if k not in reach:
                    res.append({'param': [v, k, lam, mu], 'status': 'k not a sum of unit valencies'}); continue
                name, sol, wt = solve(p, units, val, k, lam, mu, tl)
                e = {'param': [v, k, lam, mu], 'status': name, 'seconds': wt}
                if sol is not None:
                    S = sorted(c for a in sol for c in units[a])
                    e['S'] = S; e['check'] = str(sds.distance_regular(p, S))[:200]
                res.append(e)
            f.write(json.dumps({'v': v, 'i': i, 'rank': int(p.shape[0]), 'units': len(units), 'results': res,
                                'seconds': round(time.time() - t0, 1)}) + '\n'); f.flush()


if __name__ == '__main__':
    main(*sys.argv[1:5], *(float(a) for a in sys.argv[5:6]))
