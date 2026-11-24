"""Count point triples at mutual distance three by the shape of their geodesics.

Checks in the Ree-Tits model the three counts derived in grid_tree_stars.md.
"""
from collections import Counter
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "octagon"))
from geom import octagon  # noqa: E402


def main():
    _, A, D, lines = octagon()
    on = {}
    for i, line in enumerate(lines):
        for p in line:
            on.setdefault(int(p), set()).add(i)

    def start(x, y):
        (u,) = [w for w in np.flatnonzero(A[x]) if D[w, y] == D[x, y] - 1]
        (h,) = on[x] & on[int(u)]
        return h, int(u)

    def shape(x, y, z):
        (hy, uy), (hz, uz) = start(x, y), start(x, z)
        return "c" if uy == uz else "b" if hy == hz else "a"

    def through(u, v):
        (h,) = on[int(u)] & on[int(v)]
        return h

    def middles(c1, c2):
        # The three middle lines of opposite lines, with their attachment
        # points on c1, c2 and their own two used points.
        out = []
        for p in lines[c1]:
            (q,) = [q for q in lines[c2] if D[p, q] == 3]
            (u,) = [w for w in np.flatnonzero(A[p]) if D[w, q] == 2]
            (v,) = [w for w in np.flatnonzero(A[q]) if D[w, p] == 2]
            out.append((through(u, v), int(p), int(q), {int(u), int(v)}))
        return out

    def preimages(t1, t2, t3):
        # Opposite pairs U whose unused middle-line points are t1, t2, t3.
        found = 0
        for c1 in on[t1]:
            for c2 in on[t2]:
                if min(D[p, q] for p in lines[c1] for q in lines[c2]) != 3:
                    continue
                m = middles(c1, c2)
                if [r[2] for r in m if r[1] == t1] != [t2]:
                    continue
                a, b = [r[0] for r in m if r[1] != t1]
                (c3,) = [r for r in middles(a, b) if r[0] not in (c1, c2)]
                (unused,) = set(map(int, lines[c3[0]])) - c3[3]
                found += unused == t3
        return found

    x = 0
    pre = Counter()
    shell = np.flatnonzero(D[x] == 3)
    p33 = int(((D[shell[0]] == 3) & (D[x] == 3)).sum())
    names = {"ccc": "centered", "aab": "16-cycle", "aaa": "18-cycle"}
    count = Counter()
    for y in shell:
        for z in np.flatnonzero((D[x] == 3) & (D[y] == 3)):
            if z > y:
                key = "".join(sorted(shape(*t) for t in ((x, y, z), (y, x, z), (z, x, y))))
                count[names[key]] += 1
                if key == "ccc":
                    pre[preimages(x, int(y), int(z))] += 1
    total = 1755 * 640 * p33 // 6
    found = {k: 1755 * v // 3 for k, v in count.items()}  # the group is point-transitive
    predicted = {"centered": 2925 * 512, "16-cycle": 1755 * 5 * 2048}
    predicted["18-cycle"] = total - sum(predicted.values())
    print("p33", p33, "triples", total)
    for k in ("centered", "16-cycle", "18-cycle"):
        print(k, found[k], "predicted", predicted[k])
    assert found == predicted
    print("opposite pairs over centered triples through point 0", dict(sorted(pre.items())))


if __name__ == "__main__":
    main()
