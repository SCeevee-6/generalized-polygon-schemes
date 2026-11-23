"""Exhaustive search over connected regular elementary-abelian 32-sheeted trace covers (no lifted symmetry).
python elem_run.py FP [first last]  runs the level-8 prefixes [first, last) of case FP (0: same-edge mate orbit at base
edge 4, 1: cross orbit between edges 4 and 17) and writes elem_fp{FP}_{first}_{last}.json with per-prefix node counts,
leaf statistics (0 rejected by the H-sibling test, 1/2 survivors, 3 disconnected) and the surviving voltage sequences."""
import sys, os, json, time, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import elem_setup as S
fp = int(sys.argv[1]); A = S.setup_multi2([fp]); G = S.pair_stabiliser(fp)
code, n, pre, cnt, st = S.run_fast2(A, None, 1, use_sym=True, group=G, stop=8, s2=False, maxleaves=10**6)
lo, hi = (int(sys.argv[2]), int(sys.argv[3])) if len(sys.argv) > 3 else (0, n)
res = {}
for j in range(lo, min(hi, n)):
    t = time.time(); c, m, lv, cnt, st = S.run_fast2(A, list(pre[j, :A['ncol']]), 8, use_sym=True, group=G, maxleaves=20000)
    res[j] = dict(code=int(c), n=int(m), surv=lv[:, :A['ncol']].tolist(), counts=cnt.tolist(), stat=st.tolist(), sec=time.time()-t)
json.dump(dict(fp=fp, prefixes=int(n), results=res), open('elem_fp%d_%d_%d.json' % (fp, lo, hi), 'w'))
