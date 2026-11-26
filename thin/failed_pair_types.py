"""Enumerate the middle-tip metrics and check the stabilizer and shell counts."""
from collections import Counter
from itertools import combinations
import numpy as np

from grid_path_closure import LineConstraints, extend_paths
from theta_profiles import GROUP, LINES, MIDDLE, seed_distances


def completions(constraints, masks, lines):
    if not constraints.propagate(masks, lines):
        return
    pairs = np.argwhere((masks & (masks - 1)) != 0)
    if not len(pairs):
        yield np.log2(masks).astype(np.int8)
        return
    i, j = min(pairs, key=lambda pair: int(masks[tuple(pair)]).bit_count())
    for distance in range(5):
        bit = 1 << distance
        if masks[i, j] & bit:
            trial = masks.copy()
            trial[i, j] = trial[j, i] = bit
            yield from completions(constraints, trial, lines)


def on_slots(g, definitions):
    perm = list(g)
    for kind, a, b, k in definitions:
        x, y = perm[a], perm[b]
        if kind == "path" and x > y:
            x, y, k = y, x, 3-k
        matches = [21+i for i,(s,u,v,h) in enumerate(definitions)
                   if (s,h) == (kind,k) and
                   ((u,v) == (x,y) or (kind == "third" and (u,v) == (y,x)))]
        target, = matches
        perm.append(target)
    assert sorted(perm) == list(range(len(perm)))
    return perm


def shell_bounds(d, lines):
    # Distinct named observers grouped by distance from two named anchors.
    # Intersection numbers are entries of distance polynomials on the
    # quotient by the distance partition around one anchor.
    slots = [next(i for i in range(len(d)) if d[x,i] == 0) for x in range(len(d))]
    pts = sorted(set(slots))
    unique = sorted({tuple(sorted(slots[x] for x in line)) for line in lines})
    dl = np.array([[0 if a==b else 1+d[np.ix_(a,b)].min() for b in unique] for a in unique])
    for m,s,t in [(d[np.ix_(pts,pts)],2,4),(dl,4,2)]:
        b = [s*(t+1),s*t,s*t,s*t,0]
        c = [0,1,1,1,t+1]
        a = [b[0]-x-y for x,y in zip(b,c)]
        q = np.diag(a)+np.diag(b[:-1],1)+np.diag(c[1:],-1)
        polys = [np.eye(5,dtype=int),q]
        for h in range(1,4):
            v=(q-a[h]*np.eye(5,dtype=int))@polys[h]-b[h-1]*polys[h-1]
            assert not (v%c[h+1]).any()
            polys.append(v//c[h+1])
        for i in range(len(m)):
            for j in range(i+1):
                for (h,k),n in Counter(zip(m[i],m[j])).items():
                    assert n <= polys[h][m[i,j],k]


def main():
    masks, lines, definitions = extend_paths((1<<seed_distances()).astype(np.uint8),
                                             LINES, list(combinations(MIDDLE,2)))
    perms = [on_slots(g,definitions) for g in GROUP]
    seen = set()
    total = 0
    for d in completions(LineConstraints(),masks,lines):
        total += 1
        if d.tobytes() in seen:
            continue
        labels = tuple(sorted(map(int,(d[19,24],d[12,29],d[5,34]))))
        orbit = set()
        for p in perms:
            image = d[np.ix_(p,p)]
            orbit.add(image.tobytes())
            if labels != (1,1,1) and np.array_equal(d,image):
                arms = [p[j]//7 for j in (0,7,14)]
                odd = sum(a>b for a,b in combinations(arms,2))%2
                assert (p[0]%7 == 3) == odd
        shell_bounds(d, lines)
        seen.update(orbit)
        print(labels,len(orbit))
    assert len(seen) == total == 33


if __name__ == "__main__":
    main()
