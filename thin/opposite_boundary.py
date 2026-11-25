"""Verify the local counts used in the global opposite-pair matrix identities."""

from collections import Counter
from fractions import Fraction as Q
from itertools import combinations, product

import numpy as np

from failed_pair_types import completions
from grid_path_closure import LineConstraints, extend_paths
from theta_profiles import LINES, MIDDLE, seed_distances


def local_counts():
    requests = list(combinations(MIDDLE, 2))
    initial, lines, _ = extend_paths((1 << seed_distances()).astype(np.uint8), LINES, requests)
    index_by_line = {tuple(sorted(line)): i for i, line in enumerate(lines)}
    # Each unused tip singles out its own original middle line.
    centers = []
    for tip in MIDDLE:
        line, = [line for line in LINES if tip in line]
        centers.append(line)
    types = Counter()
    for distances in completions(LineConstraints(), initial, lines):
        new_centers = []
        labels = []
        outer_distances = []
        for i, tip in enumerate(MIDDLE):
            pair = tuple(x for x in MIDDLE if x != tip)
            slot = 21 + 5*requests.index(pair)
            midline = (slot, slot + 1, slot + 3)
            assert midline in index_by_line
            new_centers.append(midline)
            labels.append(int(distances[tip, slot + 3]))
            outer_distances.append(1 + int(distances[np.ix_(centers[i], midline)].min()))
        if labels == [1, 1, 1]:
            assert outer_distances == [2, 2, 2]
            types["grid"] += 1
            continue
        internal = sorted(1 + int(distances[np.ix_(a, b)].min()) for a, b in combinations(new_centers, 2))
        assert sorted(labels) in ([3, 3, 3], [3, 4, 4])
        kind = sorted(labels) == [3, 4, 4]
        assert internal == ([2, 3, 3] if kind else [3, 3, 3])
        assert outer_distances == [4, 4, 4]
        types["16-cycle" if kind else "18-cycle"] += 1
    assert types == {"grid": 1, "16-cycle": 24, "18-cycle": 8}
    print("All three omitted line distances checked in each labelled matrix", dict(types))


def moments(eigenvalue, t):
    cosine = [Q(1), Q(eigenvalue, 12)]
    for _ in range(3):
        cosine.append(((eigenvalue - 3)*cosine[-1] - cosine[-2])/8)
    a2, a3, a4 = cosine[2:]
    g = 1024*(3 + 6*a3) + 2048*t*(a2 - a3)
    return [2048*(1 + a4),
            12288*a2,
            12288*(1 + 2*a4),
            12288*(8*a2 + a4),
            36*2048*(1 + a4) + 24*6144*a4 + 4*g]


def determinant(matrix):
    copy = [list(row) for row in matrix]
    result = Q(1)
    for i in range(len(copy)):
        pivot = copy[i][i]
        if not pivot:
            return Q(0)
        result *= pivot
        for j in range(i+1, len(copy)):
            ratio = copy[j][i]/pivot
            for k in range(i+1, len(copy)):
                copy[j][k] -= ratio*copy[i][k]
    return result


def verify_positivity():
    # Both matrices are affine in the fraction of 16-cycle pairs.
    for e, t in product((7, 3, -1, -3), (0, 1)):
        p = moments(e, t)
        for matrix in ([[p[i+j] for j in range(3)] for i in range(3)],
                       [[9*p[i+j] - p[i+j+2] for j in range(2)] for i in range(2)]):
            for n in range(1, len(matrix) + 1):
                assert determinant([row[:n] for row in matrix[:n]]) > 0
    print("Both moment tests hold for every 16-cycle fraction in [0, 1]")


if __name__ == "__main__":
    local_counts()
    verify_positivity()
