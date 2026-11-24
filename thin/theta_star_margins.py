"""Star-mask partition constraints on the grid-free opposite-pair seed."""

from functools import lru_cache
import json
from pathlib import Path

import numpy as np
from scipy.sparse import coo_matrix, hstack, vstack

from theta_incidence_stars import reduced
from theta_line_margins import system
from theta_divisibility import divisors, parity_constraints


def partitions(active, masks):
    @lru_cache(None)
    def cover(remaining, slots):
        if not remaining:
            return ((0,) * slots,) if slots == 0 or 0 in masks else ()
        if not slots:
            return ()
        bit = remaining & -remaining
        return tuple((m,) + tail for m in masks if m & bit and m & remaining == m
                     for tail in cover(remaining ^ m, slots - 1))
    return cover(int(active), 5)


def star_system(with_metadata=False):
    data = reduced()
    D, P, L, pl, ll, pp, lp = data
    A, b, sizes = system(data)
    n = int(pl.max()) + 1
    v = divisors(A, b, sizes, n)
    parity, parity_rhs = parity_constraints(A, b, sizes, n, v)
    starts = len(sizes)
    m = int(ll.max()) + 1
    row, col, val = [], [], []
    nr, nc = 0, starts
    bits = 1 << np.arange(pp.shape[1], dtype=np.int64)
    total = 0
    for rep in range(n):
        i = int(np.flatnonzero(pl == rep)[0])
        p = pp[i]
        active = (p > 0) & (p < 8)
        involved = (L == i).sum(axis=1)
        incident = np.flatnonzero(involved)
        line_masks = ((lp[incident] < p) & active) @ bits
        masks = sorted(set(map(int, line_masks)))
        patterns = partitions(int(active @ bits), masks)
        assert patterns
        total += len(patterns)
        # Each point belongs to one star pattern, and every line place with
        # a specified inward mask appears in a slot of such a pattern.
        row.append(nr)
        col.append(rep)
        val.append(-1)
        for j, slots in enumerate(patterns):
            row.append(nr)
            col.append(nc + j)
            val.append(1)
        nr += 1
        for mask in masks:
            coeff = np.zeros(m, dtype=int)
            inds = incident[line_masks == mask]
            np.add.at(coeff, ll[inds], involved[inds])
            for j in np.flatnonzero(coeff):
                row.append(nr)
                col.append(n + int(j))
                val.append(int(coeff[j]))
            for j, slots in enumerate(patterns):
                k = slots.count(mask)
                if k:
                    row.append(nr)
                    col.append(nc + j)
                    val.append(-k)
            nr += 1
        nc += len(patterns)
    print("star patterns", total, "new rows", nr, flush=True)
    C = coo_matrix((val, (row, col)), shape=(nr, nc), dtype=np.int64).tocsr()
    A = vstack([hstack([A, coo_matrix((A.shape[0], nc - starts), dtype=np.int64)]), C]).tocsr()
    b = np.r_[b, np.zeros(nr, dtype=np.int64)]
    print("full matrix", A.shape, A.nnz, flush=True)
    if with_metadata:
        return A, b, sizes, v, parity, parity_rhs
    return A, b


def main():
    A, b, sizes, v, parity, parity_rhs = star_system(with_metadata=True)
    data = json.loads(Path(__file__).with_name("theta_star_margins.json").read_text())
    denominator = data["denominator"]
    assert type(denominator) is int and denominator > 0
    weights = [0] * A.shape[1]
    for key, value in data["weights"].items():
        i = int(key)
        assert 0 <= i < len(weights) and weights[i] == 0
        assert type(value) is int and value > 0
        weights[i] = value
    for i, target in enumerate(b):
        start, end = A.indptr[i:i + 2]
        total = sum(int(v) * weights[j] for j, v in zip(A.indices[start:end], A.data[start:end]))
        assert total == int(target) * denominator, i
    # These totals, not necessarily the common multiplicities of individual
    # profiles, have to be integers in any averaged actual completion.
    n = len(v)
    scaled = [weights[i] * int(sizes[i]) for i in range(n)]
    assert all(total % denominator == 0 for total in scaled)
    totals = np.asarray([total // denominator for total in scaled], dtype=np.int64)
    assert totals.sum() == 1755 and not (totals % v).any()
    assert sum(weights[i] * int(sizes[i]) for i in range(n, len(sizes))) == 2925 * denominator
    print("exact rational star witness", sum(w > 0 for w in weights), "occupied columns")
    print("integer point orbit totals, all individual divisibility conditions satisfied")
    failures = np.count_nonzero((parity @ (totals // v) - parity_rhs) % 2)
    assert failures == 0
    print("failed parity equations", failures, "of", len(parity_rhs))


if __name__ == '__main__':
    main()
