"""SRG/DRG search over primitive groups (GAP primitive groups library), one group per line of a job list.

Usage: python prim_search.py GENS_DIR JOBS_FILE OUT_JSONL
JOBS_FILE lines: "v i" ; GENS_DIR/v_i.txt holds the generators (one permutation per line, images of 1..v).
For each group: orbitals and exact p^k_ij (tits1755/coset_schemes.py), then the depth-first search of
dfs_numba.py for unions of orbitals with constant a_1 and walk counts off S in {0, c_2}; every candidate is
checked exactly with srg_drg_search.distance_regular. Writes one JSON line per group.
"""
import sys, os, json, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, '..', 'tits1755'))
import srg_drg_search as sds
import dfs_numba as dn
import coset_schemes as cs


def one(gdir, v, i):
    t0 = time.time()
    gens = [np.array([int(x) - 1 for x in l.split(',') if x], dtype=np.int64)
            for l in open(os.path.join(gdir, f'{v}_{i}.txt')).read().split('\n') if l.strip()]
    L, sizes, reps, _ = cs.orbitals(gens, v)
    pk = cs.intersection_numbers(L, reps)          # pk[k,i,j] = p^k_ij
    p = np.ascontiguousarray(pk.transpose(1, 2, 0))  # p[i,j,k]
    val, units, Q, Rsym, Qd, ucl = dn.prepare(p)
    uv = [int(val[u].sum()) for u in units]
    idx = np.array(sorted(range(len(units)), key=lambda a: (-uv[a], a)), dtype=np.int64)
    sols, nodes = dn.dfs(Q, Rsym, Qd, ucl, idx, np.zeros(0, dtype=np.int64), 1000000)
    found = []
    for s in sols:
        S = sorted(c for a in range(len(units)) if s[a] == 1 for c in units[a])
        d = sds.distance_regular(p, S)
        if isinstance(d, dict):
            e = {'S': S, 'k': int(val[S].sum()), 'diameter': d['diameter'], 'b': d['b'], 'c': d['c']}
            if d['diameter'] == 2:
                k, b1, c2 = d['b'][0], d['b'][1], d['c'][1]
                e['srg'] = [v, k, k - b1 - 1, c2]
            found.append(e)
    return {'v': v, 'i': i, 'rank': int(p.shape[0]), 'units': len(units), 'nodes': int(nodes),
            'candidates': int(len(sols)), 'graphs': found, 'seconds': round(time.time() - t0, 2)}


if __name__ == '__main__':
    gdir, jobs, out = sys.argv[1:4]
    done = set()
    if os.path.exists(out):
        for l in open(out):
            r = json.loads(l); done.add((r['v'], r['i']))
    with open(out, 'a') as f:
        for l in open(jobs):
            if not l.strip():
                continue
            v, i = map(int, l.split())
            if (v, i) in done:
                continue
            try:
                r = one(gdir, v, i)
            except Exception as e:
                r = {'v': v, 'i': i, 'error': repr(e)}
            f.write(json.dumps(r) + '\n'); f.flush()
