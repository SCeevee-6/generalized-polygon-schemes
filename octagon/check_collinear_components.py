"""Check the collinear-root component counts on all lines in the stored model."""
from pathlib import Path
import sys
import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "flags"))
from o2_build import octagon


def main():
    d, lines = octagon(np.load(ROOT / "data/generators.npy"))
    graph = csr_matrix(d == 1)
    components = np.full(d.shape, -1, dtype=np.int8)
    for x, row in enumerate(d):
        shell = np.flatnonzero(row == 4)
        count, ids = connected_components(graph[shell][:, shell], directed=False)
        assert count == 2 and np.array_equal(np.bincount(ids), [512, 512])
        components[x, shell] = ids
    # Every product and partial sum here is an exactly representable integer
    # in binary64 (entries of F are 0,+1,-1 and a dot product has 1755 terms).
    f = np.where(components < 0, 0, 1-2*components).astype(np.float64)
    gram = f @ f.T
    signed = f*f.T
    assert np.array_equal(gram, 1024*np.eye(len(d)) + 64*signed)
    assert np.array_equal(signed @ signed, 1024*np.eye(len(d)) + 48*signed)
    print("Full model contrast Gram matrix verified; signed opposition has roots 64 and -16")
    point_to_line = d[:, lines].min(axis=2)
    for roots in lines:
        for x, y in ((roots[0], roots[1]), (roots[1], roots[2]), (roots[2], roots[0])):
            a, b = components[[x, y]]
            shared = (a >= 0) & (b >= 0)
            assert np.array_equal(np.bincount(2*a[shared]+b[shared]), [128]*4)
        opposite = lines[(point_to_line[roots] == 3).all(axis=0)]
        assert len(opposite) == 2048
        code = np.zeros(len(opposite), dtype=int)
        for x in roots:
            labels = np.sort(components[x, opposite], axis=1)
            assert (labels[:, 0] == -1).all() and (labels[:, 1] == labels[:, 2]).all()
            code = 2*code + labels[:, 2]
        assert np.array_equal(np.bincount(code), [256]*8)
    print("2925 lines: four pair cells of 128 at each pair, eight triple cells of 256 per line")


if __name__ == "__main__":
    main()
