"""Vertex-primitive distance-regular graphs, step 2: search the pairs that survive drg_filter.py.

For a primitive group G (generators in GENS_DIR/V_I.txt, one comma-separated image list per line) the orbital
configuration is computed exactly, and for every surviving array {b_0,...;c_1,...} with that v (Hamming arrays
with q != 4 excluded) all unions S = S^T of orbitals with valency b_0, a_1 = b_0 - b_1 - 1 common neighbours on S
and 0 or c_2 off S are enumerated (search/dfs_fixed.py, DRG mode). Each candidate is then tested for
distance-regularity and its array compared.

Usage: python drg_search.py GENS_DIR OUT_JSONL [max_nodes] [skip v_i,v_i,...]
"""
import os, sys, csv, json, time
from collections import defaultdict
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..'))
import prim_tools as pt
import srg_drg_search as sds
import dfs_fixed as df


def main(gdir, out, maxnodes=3 * 10 ** 8, skip=''):
    skip = {tuple(int(x) for x in s.split('_')) for s in skip.split(',') if s}
    todo = defaultdict(list)
    for r in csv.DictReader(open(os.path.join(HERE, 'survivors.csv'))):
        if not r['hamming']:
            todo[(int(r['v']), int(r['i']))].append(r['array'])
    done = set()
    if os.path.exists(out):
        for l in open(out):
            r = json.loads(l); done.add((r['v'], r['i']))
    with open(out, 'a') as fo:
        for (v, i), arrs in sorted(todo.items()):
            if (v, i) in done or (v, i) in skip or not os.path.exists(os.path.join(gdir, f'{v}_{i}.txt')):
                continue
            t0 = time.time()
            p = pt.scheme(gdir, v, i)
            prep = df.prepare(p)
            res = []
            for arr in arrs:
                b, c = [[int(x) for x in part.split(',')] for part in arr.strip('{}').split(';')]
                S, nodes, nsol = df.search(p, b[0], b[0] - b[1] - 1, c[1], drg=True, maxsol=10000,
                                           maxnodes=maxnodes, prep=prep)
                hits = []
                for s_ in S:
                    d = sds.distance_regular(p, s_)
                    if isinstance(d, dict) and list(d['b']) == b and list(d['c']) == c:
                        hits.append(s_)
                res.append({'array': arr, 'status': 'LIMIT' if nsol < 0 else 'DONE', 'candidates': int(max(nsol, 0)),
                            'hits': hits[:20], 'n_hits': len(hits), 'nodes': nodes})
            rec = {'v': v, 'i': i, 'rank': int(p.shape[0]), 'units': len(prep[1]), 'results': res,
                   'seconds': round(time.time() - t0, 1)}
            fo.write(json.dumps(rec) + '\n'); fo.flush()
            print(v, i, rec['rank'], [(e['array'], e['status'], e['n_hits']) for e in res], rec['seconds'], flush=True)


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2], *([int(sys.argv[3])] if len(sys.argv) > 3 else []), *sys.argv[4:5])
