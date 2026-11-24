"""Divisibility of profile-orbit totals from integral line counts."""

from functools import reduce
from math import gcd, lcm

import numpy as np

from theta_incidence_stars import reduced
from theta_line_margins import system


def divisors(A, b, sizes, n):
    """For every mixed row, derive a necessary divisor of its point total."""
    scale = 12 // sizes
    out = np.ones(n, dtype=np.int64)
    for i, target in enumerate(b):
        row = A.getrow(i)
        points = [(int(j), int(v * scale[j])) for j, v in zip(row.indices, row.data) if j < n]
        lines = [int(v * scale[j]) for j, v in zip(row.indices, row.data) if j >= n]
        if not lines:
            continue
        g = reduce(gcd, lines, 0)
        if not points:
            assert (12 * int(target)) % g == 0
        else:
            # Degree and inward-count equations have one point coefficient
            # and homogeneous right-hand side.
            assert len(points) == 1 and target == 0
            j, a = points[0]
            out[j] = lcm(int(out[j]), g // gcd(a, g))
    return out


def parity_constraints(A, b, sizes, n, v):
    """Eliminate integral line orbit totals over F_2.

    The remaining equations are on U_i = T_i/v_i for the point orbit
    totals T_i. Normalize by the gcd of line coefficients before reducing.
    """
    m = len(sizes) - n
    scale = 12 // sizes
    basis = {}
    for i, target in enumerate(b):
        row = A.getrow(i)
        point = [(int(j), int(c * scale[j] * v[j])) for j, c in zip(row.indices, row.data) if j < n]
        line = [(int(j - n), int(c * scale[j])) for j, c in zip(row.indices, row.data) if j >= n]
        g = reduce(gcd, (c for j, c in line or point), 0)
        if not g:
            assert target == 0
            continue
        rhs = 12 * int(target)
        assert rhs % g == 0 and all(c % g == 0 for j, c in point)
        bits = (rhs // g % 2) << (m + n)
        for j, c in line:
            bits |= (c // g % 2) << j
        for j, c in point:
            bits |= (c // g % 2) << (m + j)
        # Increasing pivot index eliminates line variables first.
        while bits:
            pivot = (bits & -bits).bit_length() - 1
            if pivot in basis:
                bits ^= basis[pivot]
            else:
                assert pivot < m + n, "inconsistent parity system"
                basis[pivot] = bits
                break
    rows = []
    rhs = []
    for pivot, bits in sorted(basis.items()):
        if pivot >= m:
            rows.append([(bits >> (m + i)) & 1 for i in range(n)])
            rhs.append((bits >> (m + n)) & 1)
    return np.asarray(rows, dtype=np.int64), np.asarray(rhs, dtype=np.int64)


if __name__ == '__main__':
    data = reduced()
    A, b, sizes = system(data)
    v = divisors(A, b, sizes, int(data[3].max()) + 1)
    nontrivial = [(i, int(x)) for i, x in enumerate(v) if x > 1]
    assert nontrivial == [(0, 3), (1, 2), (3, 2), (6, 2), (9, 2), (10, 2),
                          (21, 2), (22, 2), (23, 2)]
    print("nontrivial profile orbit divisors", nontrivial)
    C, d = parity_constraints(A, b, sizes, len(v), v)
    point_only = A[:, len(v):].getnnz(axis=1) == 0
    C0, d0 = parity_constraints(A[point_only], b[point_only], sizes, len(v), v)
    assert (len(d), len(d0)) == (52, 43)
    print("independent parity equations", len(d), "point-only", len(d0),
          "additional", len(d) - len(d0), flush=True)
