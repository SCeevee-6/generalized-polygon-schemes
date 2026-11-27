"""Verify the nine lines in J^3 D^T and the sixth walk distribution."""
from collections import Counter
from itertools import combinations
from fractions import Fraction

import numpy as np

from nb_m5_cross import REQUESTS, prepare
from nb_moment_sdp import known_distributions, walk_expansion
from opposite_link_surfaces import failed_metrics
from grid_path_closure import LineConstraints, extend_paths
from theta_profiles import LINES, seed_distances


def sixth_distribution(p, u, v):
    """Divide absolute totals by the number of ambient lines, as in Q_0..Q_5."""
    return [1024*x for x in (24, 0, 32*p+8*u, 96-32*p+8*v, 264-8*(u+v))]


def main():
    constraints = LineConstraints()
    _, lines, _ = extend_paths((1 << seed_distances()).astype(np.uint8), LINES, REQUESTS)
    # In REQUESTS order the omitted arms are 2, 1, 0.
    middle = [LINES[3 + 3*i] for i in (2, 1, 0)]
    distance = 1 + constraints.squares.min(axis=(1, 2))
    distance[(constraints.squares == 0).sum(axis=(1, 2)) == 3] = 0
    shape_totals = Counter()
    side_totals = Counter()
    for metric in failed_metrics(constraints):
        masks, _, extra = prepare(constraints, lines, metric)

        def values(first, second):
            agree = ((1 << constraints.squares) & masks[np.ix_(first, second)]) != 0
            return set(distance[agree.all(axis=(1, 2))])

        close_sides = []
        for i in range(3):
            for side, root in enumerate("AB"):
                y = extra[i, root]
                assert values(y, LINES[side]) == {2}
                assert values(y, LINES[1-side]) == {4}
                central_values = values(y, middle[i])
                assert central_values in ({2}, {3})
                if central_values == {2}:
                    close_sides.append(root)
                for j in range(3):
                    if i == j:
                        continue
                    assert values(y, middle[j]) == {4}
                    assert values(y, extra[j, root]) == {4}
            assert values(extra[i, "A"], extra[i, "B"]) in ({2}, {3}, {2, 3})
        possible_twos = 0
        for i, j in combinations(range(3), 2):
            for x, y in ((i, j), (j, i)):
                allowed = values(extra[x, "A"], extra[y, "B"])
                assert allowed in ({3, 4}, {2, 3, 4})
                possible_twos += 2 in allowed
        shape = tuple(sorted(metric[[5, 12, 19], [34, 29, 24]]))
        shape_totals[shape, possible_twos] += 1
        side_totals[shape, len(close_sides), len(set(close_sides))] += 1
    assert shape_totals == {((3, 3, 3), 0): 8, ((3, 4, 4), 2): 12, ((3, 4, 4), 3): 12}

    assert side_totals == {((3, 3, 3), 0, 0): 8, ((3, 4, 4), 2, 1): 12, ((3, 4, 4), 2, 2): 12}

    coeff = walk_expansion(6)[6]
    assert coeff == [87, 0, 47, 0, 11, 0, 1]
    for p, u, v in ((0, 0, 0), (1, 0, 0), (0, 6, 0), (0, 0, 6)):
        # Ordered line pairs at the middle vertex, weights 14 on each C and 2 on each Y.
        cc = [3, 0, 0, 0, 6]
        cy = [0, 0, 2*p, 6-2*p, 12]  # one orientation
        yy = [6, 0, 2*p+2*u, 6-2*p+2*v, 12+2*(6-u-v)]
        gram = [1024*(196*a+56*b+4*c) for a, b, c in zip(cc, cy, yy)]
        low = known_distributions(p)
        result = [gram[i]-sum(coeff[k]*low[k][i] for k in (0, 2, 4)) for i in range(5)]
        assert result == sixth_distribution(p, u, v) and sum(result) == 384*1024
    print("32 failed metrics; six-step returns per line 24576, distance-one entries 0")


def positive_definite(matrix):
    # Exact symmetric elimination. Every pivot is a positive rational.
    for i, row in enumerate(matrix):
        pivot = row[i]
        assert pivot > 0
        for j in range(i+1, len(matrix)):
            for k in range(i+1, len(matrix)):
                matrix[j][k] -= matrix[j][i]*row[k]/pivot


def moment_witness():
    from nb_moment_sdp import EIGENVALUES, q

    p, t, s, u, v = 0, Fraction(1, 2), 2, 0, 3
    rows = known_distributions(p) + [[0, 0, 98304, 12288*s, 98304-12288*s],
                                     sixth_distribution(p, u, v),
                                     [0, 0, 250244, 458591, 77597],
                                     [55296, 0, 0, 12288, 1505280]]
    assert 0 <= s <= 4-2*p and 0 <= u <= 3*p and v >= 0 and u+v <= 6
    for k in range(1, 9):
        assert min(rows[k]) >= 0 and sum(rows[k]) == 12288*2**(k-1)
    assert rows[7][0] <= 16384*p
    assert rows[7][2] >= 36864+12288*(1-t) and sum(rows[7][3:]) >= 12288*t
    assert rows[8][0] >= 1024*(60-12*t) and rows[8][4] >= 96*1024
    assert sum(rows[8][1:4]) >= 12288*t and sum(rows[8][2:4]) >= 6144*t
    coeff = walk_expansion(8)
    for e in EIGENVALUES:
        values = [q(row, e) for row in rows]
        m = [sum(a*b for a, b in zip(c, values)) for c in coeff]
        if e == 12:
            assert m == [m[0]*3**j for j in range(9)]
            continue
        positive_definite([[m[i+j] for j in range(5)] for i in range(5)])
        for sign in (-1, 1):
            positive_definite([[3*m[i+j]+sign*m[i+j+1] for j in range(4)] for i in range(4)])
        positive_definite([[9*m[i+j]-m[i+j+2] for j in range(4)] for i in range(4)])
    print("Rational moment witness through order eight at p=0, t=1/2")


if __name__ == "__main__":
    main()
    moment_witness()
