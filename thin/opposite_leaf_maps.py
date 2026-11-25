"""Check the 32-line bijections at a fixed flag in the Ree–Tits model."""
from itertools import combinations
import numpy as np

from check_opposite_seed import ROOT, octagon


def main():
    distance, lines = octagon(np.load(ROOT / "data/generators.npy"))
    through = {tuple(sorted(edge)): i for i, line in enumerate(lines)
               for edge in combinations(line, 2)}

    def leaves(point, parent):
        # Three incidence edges from point, with first line different from parent.
        near = distance[point, lines]
        ids = np.flatnonzero(near.min(axis=1) == 1)
        points = lines[ids, near[ids].argmin(axis=1)]
        first = np.array([through[tuple(sorted((point, p)))] for p in points])
        keep = first != parent
        assert keep.sum() == 32
        return ids[keep], first[keep], points[keep]

    def match(source, target):
        src, first, second = source
        dst, other_first, other_second = target
        square = distance[lines[src][:, None, :, None], lines[dst][None, :, None, :]]
        good = square.min(axis=(2, 3)) == 1
        assert np.all(good.sum(axis=0) == 1) and np.all(good.sum(axis=1) == 1)
        image = good.argmax(axis=1)
        for i, j in combinations(range(32), 2):
            k = int(first[i] == first[j]) + int(second[i] == second[j])
            h, t = image[i], image[j]
            l = int(other_first[h] == other_first[t]) + int(other_second[h] == other_second[t])
            assert k+l <= 2
        return image

    parent = 0
    u, a, other_a = map(int, lines[parent])
    source, target = leaves(a, parent), leaves(other_a, parent)
    total = np.zeros((32, 32), dtype=int)
    blocks = 0
    for middle in np.flatnonzero((lines == u).any(axis=1)):
        if middle == parent:
            continue
        for v in lines[middle]:
            if v == u:
                continue
            columns = []
            for c_line in np.flatnonzero((lines == v).any(axis=1)):
                if c_line == middle:
                    continue
                images = []
                for c in lines[c_line]:
                    if c != v:
                        dst = leaves(c, c_line)
                        f, g = match(source, dst), match(target, dst)
                        images.append(np.argsort(g)[f])
                assert np.array_equal(*images)  # the Ree–Tits repeat
                columns += images
            values = np.array(columns).T
            assert values.shape == (32, 8)
            for d, row in enumerate(values):
                _, mult = np.unique(row, return_counts=True)
                assert list(mult) == [2, 2, 2, 2]
                np.add.at(total[d], row, 1)
            blocks += 1
    assert blocks == 8 and np.all(total == 2)
    print("source leaves 32, blocks 8, bijections 128")
    print("each block: four distinct outputs per source, each twice")
    print("each source-target pair counted twice across all eight blocks")


if __name__ == "__main__":
    main()
