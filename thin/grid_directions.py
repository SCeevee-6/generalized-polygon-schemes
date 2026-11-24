"""Enumerate local directions in a nine-tip grid with no centered transversal."""

from collections import Counter
from fractions import Fraction
from itertools import combinations, combinations_with_replacement, permutations, product
import json
from math import lcm
from pathlib import Path

import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import shortest_path

from theta_profiles import distance_polynomials
from grid_geodesics import check_geodesics

CELLS = range(9)
OFF = [[y for y in CELLS if x // 3 != y // 3 and x % 3 != y % 3] for x in CELLS]
SYMMETRIES = [
    (
        tuple(3 * c[j] + r[i] if swap else 3 * r[i] + c[j] for i in range(3) for j in range(3)),
        swap,
    )
    for r in permutations(range(3))
    for c in permutations(range(3))
    for swap in (False, True)
]


def matching_type(y, z):
    if y // 3 == z // 3:
        return 0
    if y % 3 == z % 3:
        return 1
    return 2


def cases():
    seen = set()
    for directions in product(range(3), repeat=9):
        if directions in seen:
            continue
        orbit = set()
        for g, swap in SYMMETRIES:
            word = tuple(directions[i] for i in g)
            if swap:
                word = tuple((1, 0, 2)[i] for i in word)
            orbit.add(word)
        assert not seen.intersection(orbit)
        seen.update(orbit)
        yield directions, len(orbit)
    assert len(seen) == 3**9


def skeleton(directions):
    """Return the 81 tip-distance profiles and the 51 forced lines, or a collision."""
    lines = [tuple(9 + 2 * (3 * i + j) for j in range(3)) for i in range(3)]
    lines += [tuple(10 + 2 * (3 * i + j) for i in range(3)) for j in range(3)]
    profiles = np.full((81, 9), 4, dtype=np.int8)
    profiles[:9] = 3
    np.fill_diagonal(profiles[:9], 0)
    partner = {}
    neighbour = {}
    for x in CELLS:
        row, col = 9 + 2 * x, 10 + 2 * x
        lines.append((x, row, col))
        profiles[[row, col], x] = 1
        for y in CELLS:
            if x == y:
                continue
            if x // 3 == y // 3:
                profiles[[row, col], y] = (2, 3)
            if x % 3 == y % 3:
                profiles[[row, col], y] = (3, 2)
        for j, y in enumerate(OFF[x]):
            (z,) = (k for k in OFF[x] if k != y and matching_type(y, k) == directions[x])
            partner[x, y] = z
            neighbour[x, y] = 27 + 4 * x + j
            profiles[neighbour[x, y], [x, y, z]] = (1, 2, 3)
        for y in OFF[x]:
            z = partner[x, y]
            if y < z:
                lines.append((x, neighbour[x, y], neighbour[x, z]))
    mid = 63
    for x in CELLS:
        for y in OFF[x]:
            if x >= y:
                continue
            v, w = partner[x, y], partner[y, x]
            if v == w:
                return None
            profiles[mid] = 3
            profiles[mid, [x, y, v, w]] = (2, 2, 4, 4)
            lines.append((neighbour[x, y], neighbour[y, x], mid))
            mid += 1
    assert mid == 81 and len(lines) == 51
    return profiles, lines


def third_domains():
    # An entry is the bit mask of possible distances to the third point
    # of a line, given the masks at the first two points.
    table = np.zeros((32, 32), dtype=np.uint8)
    for d in range(4):
        for i, j, k in set(permutations((d, d + 1, d + 1))):
            for a in range(1, 32):
                for b in range(1, 32):
                    if (a & (1 << i)) and (b & (1 << j)):
                        table[a, b] |= 1 << k
    return table


def distance_domains(profiles, lines, table):
    adj = np.zeros((81, 81), dtype=int)
    for line in lines:
        for i, j in combinations(line, 2):
            adj[i, j] = adj[j, i] = 1
    path = shortest_path(csr_matrix(adj), unweighted=True)
    low = np.abs(profiles[:, None, :] - profiles[None, :, :]).max(axis=2)
    high = (profiles[:, None, :] + profiles[None, :, :]).min(axis=2)
    masks = np.zeros((81, 81), dtype=np.uint8)
    for d in range(5):
        allowed = (low <= d) & (d <= high) & ((path == d) | (path > 4))
        if d == 0:
            allowed &= np.eye(81, dtype=bool)
        masks |= allowed.astype(np.uint8) << d
    for i in CELLS:
        known = (1 << profiles[:, i]).astype(np.uint8)
        masks[i, :] &= known
        masks[:, i] &= known
    while True:
        old = masks.copy()
        for i, j, k in lines:
            masks[:, i] &= table[masks[:, j], masks[:, k]]
            masks[:, j] &= table[masks[:, i], masks[:, k]]
            masks[:, k] &= table[masks[:, i], masks[:, j]]
        masks &= masks.T
        if not masks.all():
            return None
        if np.array_equal(masks, old):
            return masks


def star_patterns(degree):
    patterns = []
    for d in range(5):
        allowed = []
        if d < 4:
            allowed.append((d + 1, d + 1))
        if d > 0:
            allowed.extend(((d - 1, d), (d, d - 1)))
        for slots in product(allowed, repeat=degree):
            if d < 4 and sum(min(pair) < d for pair in slots) > 1:
                continue
            patterns.append((d,) + sum(slots, ()))
    return np.left_shift(1, np.asarray(patterns)).astype(np.uint8)


def star_domains(masks, lines):
    # Fix an observer in each row and project the joint distance pattern
    # at each point and all neighbours on its named lines.
    stars = []
    for p in range(81):
        slots = [tuple(i for i in line if i != p) for line in lines if p in line]
        stars.append(((p,) + sum(slots, ()), star_patterns(len(slots))))
    masks = masks.copy()
    while True:
        old = masks.copy()
        for vertices, patterns in stars:
            allowed = ((masks[:, vertices][:, None, :] & patterns) != 0).all(axis=2)
            masks[:, vertices] &= np.bitwise_or.reduce(allowed[:, :, None] * patterns, axis=1)
        masks &= masks.T
        if not masks.all():
            return None
        if np.array_equal(masks, old):
            return masks


def psd_rank(matrix):
    # Fraction-free symmetric elimination. A zero diagonal in a PSD
    # matrix forces its entire row to vanish.
    a = matrix.astype(object)
    previous, rank = 1, 0
    while a.size:
        diag = np.diag(a)
        assert (diag >= 0).all() and not a[diag == 0].any()
        a = a[np.ix_(diag > 0, diag > 0)]
        if not a.size:
            break
        pivot = a[0, 0]
        numerator = pivot * a[1:, 1:] - np.outer(a[1:, 0], a[0, 1:])
        assert not (numerator % previous).any()
        a = numerator // previous
        previous, rank = pivot, rank + 1
    return rank


def verify_metric(masks, profiles, lines, digits):
    # Digits specify only ambiguous upper-triangular entries, in row order.
    digits = iter(digits)
    dist = np.zeros_like(masks, dtype=int)
    for i, j in combinations(range(81), 2):
        mask = int(masks[i, j])
        d = mask.bit_length() - 1 if mask.bit_count() == 1 else int(next(digits))
        assert 1 <= d <= 4 and mask & (1 << d)
        dist[i, j] = dist[j, i] = d
    assert next(digits, None) is None
    assert np.array_equal(dist[:, :9], profiles)
    for line in lines:
        values = np.sort(dist[:, line], axis=1)
        assert (values[:, 1] == values[:, 0] + 1).all() and (values[:, 1] == values[:, 2]).all()
    for i in range(81):
        inward = sum(
            (dist[:, line].min(1) < dist[:, i]).astype(int) for line in lines if i in line
        )
        assert ((dist[:, i] == 4) | (inward <= 1)).all()
    capacity = np.array(distance_polynomials()).transpose(1, 0, 2)
    for i, j in combinations_with_replacement(range(81), 2):
        shell = np.bincount(dist[i] * 5 + dist[j], minlength=25).reshape(5, 5)
        assert (shell <= capacity[dist[i, j]]).all()
    for eigenvalue, multiplicity in ((-5, 78), (-3, 675), (1, 650), (5, 351)):
        seq = [Fraction(1), Fraction(eigenvalue, 10)]
        for _ in range(3):
            seq.append(((eigenvalue - 1) * seq[-1] - seq[-2]) / 8)
        scale = lcm(*(x.denominator for x in seq))
        cosine = np.asarray([int(scale * x) for x in seq])
        rank = psd_rank(cosine[dist])
        assert rank <= multiplicity
        if eigenvalue == -5:
            assert rank == 30

    check_geodesics(dist, lines)


def exclude_middle_line_case(reduced, lines):
    # For (r,a,b,c) below, c is the third point on the middle line from r to b.
    # Thus a neighbour a of r must satisfy d(a,b) = d(a,c)+1.
    forbidden = np.array([[0, 1, 3, 2], [1, 0, 3, 3], [3, 3, 0, 2], [2, 3, 2, 0]])
    assert reduced[28, 37] == (1 << 3 | 1 << 4)
    for d, four in ((3, [35, 40, 49, 78]), (4, [36, 44, 50, 76])):
        branch = reduced.copy()
        branch[28, 37] = branch[37, 28] = 1 << d
        after = star_domains(branch, lines)
        assert after is not None and np.array_equal(after[np.ix_(four, four)], 1 << forbidden)


def main():
    totals = Counter()
    classes = Counter()
    table = third_domains()
    witnesses = json.loads(Path(__file__).with_name("grid_metrics.json").read_text())
    used = set()
    for directions, size in cases():
        word = "".join(map(str, directions))
        data = skeleton(directions)
        if data is None:
            result = "coincident partners"
        else:
            profiles, lines = data
            masks = distance_domains(profiles, lines, table)
            if masks is None:
                result = "inconsistent line distances"
            elif word in witnesses:
                verify_metric(masks, profiles, lines, witnesses[word])
                used.add(word)
                result = "verified partial metrics"
            else:
                reduced = star_domains(masks, lines)
                if reduced is None:
                    result = "inconsistent inward lines"
                elif word == "001001220":
                    exclude_middle_line_case(reduced, lines)
                    result = "inconsistent middle line"
                else:
                    # A two-branch refutation for the remaining negative case.
                    assert word == "001010211" and reduced[27, 41] == (1 << 3 | 1 << 4)
                    for d in (3, 4):
                        branch = reduced.copy()
                        branch[27, 41] = branch[41, 27] = 1 << d
                        assert star_domains(branch, lines) is None
                    result = "inconsistent after distance split"
        totals[result] += size
        classes[result] += 1
    assert used == set(witnesses)
    print("directions", sum(totals.values()), "symmetry classes", sum(classes.values()))
    for key in (
        "coincident partners",
        "inconsistent line distances",
        "inconsistent inward lines",
        "inconsistent after distance split",
        "inconsistent middle line",
        "verified partial metrics",
    ):
        print(key, totals[key], "classes", classes[key])


if __name__ == "__main__":
    main()
