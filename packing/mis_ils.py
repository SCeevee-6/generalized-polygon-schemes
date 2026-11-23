import numpy as np, numba as nb
@nb.njit(cache=True)
def ils(adjp, adji, n, iters, seed, init):
    np.random.seed(seed)
    insol = np.zeros(n, np.bool_); tight = np.zeros(n, np.int64)
    for v in init:
        if tight[v] == 0 and not insol[v]:
            insol[v] = True
            for k in range(adjp[v], adjp[v+1]): tight[adji[k]] += 1
    best = insol.copy(); bsize = insol.sum(); size = bsize
    for it in range(iters):
        # greedy fill of free vertices
        perm = np.random.permutation(n)
        for v in perm:
            if not insol[v] and tight[v] == 0:
                insol[v] = True; size += 1
                for k in range(adjp[v], adjp[v+1]): tight[adji[k]] += 1
        # (1,2)-swaps
        improved = True
        while improved:
            improved = False
            for x in np.random.permutation(n):
                if not insol[x]: continue
                # candidates: 1-tight neighbours of x
                cnt = 0
                cand = np.empty(adjp[x+1]-adjp[x], np.int64)
                for k in range(adjp[x], adjp[x+1]):
                    w = adji[k]
                    if tight[w] == 1: cand[cnt] = w; cnt += 1
                if cnt < 2: continue
                found = False
                for i in range(cnt):
                    for j in range(i+1, cnt):
                        u = cand[i]; w = cand[j]; adjuw = False
                        for k in range(adjp[u], adjp[u+1]):
                            if adji[k] == w: adjuw = True; break
                        if not adjuw:
                            insol[x] = False
                            for k in range(adjp[x], adjp[x+1]): tight[adji[k]] -= 1
                            for z in (u, w):
                                insol[z] = True
                                for k in range(adjp[z], adjp[z+1]): tight[adji[k]] += 1
                            size += 1; found = True; break
                    if found: break
                if found:
                    improved = True
                    for v in np.random.permutation(n):
                        if not insol[v] and tight[v] == 0:
                            insol[v] = True; size += 1
                            for k in range(adjp[v], adjp[v+1]): tight[adji[k]] += 1
        if size > bsize:
            bsize = size; best = insol.copy()
        elif size < bsize - 1:
            insol = best.copy(); size = bsize
            tight[:] = 0
            for v in range(n):
                if insol[v]:
                    for k in range(adjp[v], adjp[v+1]): tight[adji[k]] += 1
        # perturb: force in 1-2 random non-solution vertices
        for _ in range(1 + (np.random.random() < 0.3)):
            v = np.random.randint(n)
            if insol[v]: continue
            for k in range(adjp[v], adjp[v+1]):
                w = adji[k]
                if insol[w]:
                    insol[w] = False; size -= 1
                    for kk in range(adjp[w], adjp[w+1]): tight[adji[kk]] -= 1
            insol[v] = True; size += 1
            for k in range(adjp[v], adjp[v+1]): tight[adji[k]] += 1
    return best, bsize
