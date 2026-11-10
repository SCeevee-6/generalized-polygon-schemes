"""Cross-check of the CP-SAT model (search_cpsat_cuts) against the exhaustive Gray-code search.

For each scheme, the Gray-code candidates (constant a_1) are filtered by the extra conditions of the CP-SAT
model (walk counts off S in {0, c_2} with one value c_2 > 0, and S in no maximal closed subset). The two
lists must coincide. Usage: python validate_cuts.py <csv> [<csv> ...]"""
import os, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ortools.sat.python import cp_model  # noqa: F401  (load before numpy-heavy imports on Windows)
import srg_drg_search as sd

for fp in sys.argv[1:]:
    p = sd.load_p(fp)
    val, pair, units, cand, _ = sd.search(p)
    mcs = sd.maximal_closed_subsets(p)
    g = []
    for us in cand:
        S = [i for a in us for i in units[a]]
        c = sd.walk(p, S, S)
        nz = {int(c[k]) for k in range(1, p.shape[0]) if k not in S and c[k] > 0}
        if len(nz) == 1 and not any(set(S) <= set(C) for C in mcs):
            g.append(tuple(us))
    t = time.time()
    b = sd.search_cpsat_cuts(p, workers=4, time_limit=3000)
    cp = sorted(map(tuple, b[3]))
    print(os.path.basename(fp), 'gray', sorted(g), 'cpsat', cp, b[4], 'agree' if sorted(g) == cp else 'DIFFER',
          f'{time.time() - t:.0f}s', flush=True)
