"""Try two squares at the same opposite roots against the three-line test."""

from itertools import combinations, permutations

import numpy as np

from grid_path_closure import LineConstraints, extend_paths
from opposite_link_surfaces import failed_metrics
from theta_profiles import LINES, MIDDLE, seed_distances


def check_overlaps():
    constraints = LineConstraints()
    distances = failed_metrics(constraints)[0]
    requests = list(combinations(MIDDLE, 2))
    _, lines, definitions = extend_paths((1 << seed_distances()).astype(np.uint8), LINES, requests)
    assert len(distances) == 36
    centers = [next(line for line in LINES if tip in line) for tip in MIDDLE]
    # Use the first two arms for the second vertex of the cubic edge.
    # Its third middle line is the middle of the first requested tip path.
    assert definitions[:2] == [("path", *requests[0], k) for k in (1, 2)]
    x = np.array(centers[2])
    y = np.array(next(line for line in lines if 21 in line and 22 in line))
    projection = distances[np.ix_(x, y)]
    assert sorted(projection.flatten()) == [3]*3 + [4]*6
    assert (projection == 3).sum(axis=0).tolist() == [1, 1, 1]
    y = y[projection.argmin(axis=1)]
    assert np.array_equal(distances[np.ix_(x, y)], 4 - np.eye(3, dtype=int))
    corners = [lines.index(line) for line in [*LINES[:2], *centers[:2]]]

    def line_distance(a, b):
        return 1 + int(distances[np.ix_(a, b)].min())

    left = [lines[i] for i in corners[:2]]
    right = [lines[i] for i in corners[2:]]
    assert all((line_distance(x, a), line_distance(y, a)) == (2, 4) for a in left)
    assert all((line_distance(x, b), line_distance(y, b)) == (4, 2) for b in right)
    assert all(line_distance(a, b) == 2 for a in left for b in right)

    metric = (1 << distances).astype(np.uint8)
    distinct = sum(1 << k for k, square in enumerate(constraints.squares)
                   if np.count_nonzero(square == 0) != 3)
    for perm in permutations(range(3)):
        maps = [np.arange(36)]
        # Second square, then the failed-pair metric at the common roots.
        for old, named in [(np.r_[x, y], np.r_[x[list(perm)], y[list(perm)]]),
                           (np.r_[LINES[0], LINES[1]], np.r_[x, y])]:
            slots = np.full(36, -1)
            slots[old] = named
            start = 1 + int(maps[-1].max())
            slots[slots < 0] = np.arange(start, start + 30)
            maps.append(slots)
        size = 1 + int(maps[-1].max())
        masks = np.full((size, size), 31, dtype=np.uint8)
        joined = []
        for slots in maps:
            masks[np.ix_(slots, slots)] &= metric
            joined.extend(tuple(slots[list(line)]) for line in lines)
        domains = np.full((len(joined),)*2, (1 << 39) - 1, dtype=np.uint64)
        for a in corners:
            for b in corners:
                domains[a, len(lines)+b] = domains[len(lines)+b, a] = distinct
        assert constraints.propagate(masks, joined, line_domains=domains)
    print("Two squares and a failed root metric pass three-line propagation in all six gluings")


if __name__ == "__main__":
    check_overlaps()
