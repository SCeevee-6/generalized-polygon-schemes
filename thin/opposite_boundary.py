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
    centers = [lines[4 + 3*i] for i in range(3)]
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
            continue
        internal = sorted(1 + int(distances[np.ix_(a, b)].min()) for a, b in combinations(new_centers, 2))
        kind = int(sorted(labels) == [3, 4, 4])
        assert internal == ([2, 3, 3] if kind else [3, 3, 3])
        for label, outer in zip(labels, outer_distances):
            assert outer in ((4,) if label == 3 else (3, 4))
        types[kind, outer_distances.count(3)] += 1
    assert types == {(0, 0): 8, (1, 0): 6, (1, 1): 12, (1, 2): 6}
    print("Labelled failed matrices by shape and number of marked edges", dict(sorted(types.items())))


def moments(eigenvalue, t, q):
    cosine = [Q(1), Q(eigenvalue, 12)]
    for _ in range(3):
        cosine.append(((eigenvalue - 3)*cosine[-1] - cosine[-2])/8)
    a2, a3, a4 = cosine[2:]
    g = 1024*(3 + 6*a3) + 2048*t*(a2 - a3)
    return [2048*(1 + a4),
            12288*a2,
            12288*(1 + 2*a4),
            12288*(8*a2 + q*a3 + (1 - q)*a4),
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
    # The required matrices are affine on the triangle 0 <= q <= 2t/3.
    for e, (t, q) in product((7, 3, -1, -3), ((0, Q(0)), (1, Q(0)), (1, Q(2, 3)))):
        p = moments(e, t, q)
        for matrix in ([[p[i+j] for j in range(3)] for i in range(3)],
                       [[9*p[i+j] - p[i+j+2] for j in range(2)] for i in range(2)]):
            for n in range(1, len(matrix) + 1):
                assert determinant([row[:n] for row in matrix[:n]]) > 0
    print("Both moment tests hold throughout the parameter triangle")


if __name__ == "__main__":
    local_counts()
    verify_positivity()
