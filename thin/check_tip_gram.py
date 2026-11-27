"""Check one complete tip-incidence row and five matchings in the Ree–Tits model."""
import numpy as np
from check_opposite_seed import ROOT, octagon


def main():
    pd, lines = octagon(np.load(ROOT / "data/generators.npy"))
    x = 0
    to_line = pd[:, lines].min(axis=2)
    near = np.minimum.reduce([to_line[lines[:, i]] for i in range(3)]) == 1
    counts = np.zeros(len(pd), dtype=int)
    for c in np.flatnonzero((lines == x).any(axis=1)):
        y, z = lines[c][lines[c] != x]
        rows = np.flatnonzero(near[c] & (to_line[y] == 1))
        cols = np.flatnonzero(near[c] & (to_line[z] == 1))
        assert len(rows) == len(cols) == 32
        uses = np.zeros(len(lines), dtype=int)
        for a in rows:
            for b in cols:
                mids = np.flatnonzero(near[a] & near[b])
                assert len(mids) == 3 and c in mids
                for j in mids[mids != c]:
                    tips = lines[j]
                    t, = tips[(to_line[tips, a] != 1) & (to_line[tips, b] != 1)]
                    counts[t] += 1
                    uses[j] += 1
        cross = ((to_line[x] == 3) & (to_line[y] == 3) & (to_line[z] == 3))
        assert cross.sum() == 2048 and np.array_equal(uses, cross)
    assert np.array_equal(counts, 16*(pd[x] == 3))
    print("Five matchings of 1024 edges and all off-diagonal entries of one tip Gram row verified")


if __name__ == "__main__":
    main()
