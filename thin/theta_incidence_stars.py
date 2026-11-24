"""Complete incidence directions relative to the 32 named vertices."""

from functools import lru_cache
from itertools import combinations

import numpy as np

from theta_line_margins import geometry


def line_rule(pp, lines, lp):
    """Distances from every named vertex to the three points of a line."""
    values = pp[lines]
    inward = (values == lp[:, None, :] - 1).sum(axis=1)
    outward = (values == lp[:, None, :] + 1).sum(axis=1)
    return np.all(np.where(lp == 0, outward == 3,
                           np.where(lp == 8, inward == 3, (inward == 1) & (outward == 2))), axis=1)


def compatible(pp, lines, lp, dist):
    """The initial one or two edges of short paths must fit the pair distances."""
    good = np.ones(len(lines), dtype=bool)
    values = pp[lines]
    # First edges of geodesics from the line itself, where all three
    # possible neighbours are known from its type.
    nearest = values.argmin(axis=1)
    active = (lp > 0) & (lp < 8)
    for y, z in combinations(range(pp.shape[1]), 2):
        condition = active[:,y] & active[:,z]
        length = lp[:,y] + lp[:,z]
        d = dist[y, z]
        same = nearest[:,y] == nearest[:,z]
        good &= ~(condition & same & (d > length-2))
        invalid = np.where(length <= 8, d != length, d < 16-length)
        good &= ~(condition & ~same & invalid)
    for i in range(3):
        h = values[:, i]
        active = (h > 0) & (h < 8)
        inward = active & (lp == h - 1)
        others = [j for j in range(3) if j != i]
        direction = (values[:, others[0]] == h - 2).astype(int) + 2 * (values[:, others[1]] == h - 2)
        assert np.all(direction[inward & (h > 1)] > 0)
        assert np.all(direction[inward & (h > 1)] < 3)
        for y, z in combinations(range(pp.shape[1]), 2):
            a, b = inward[:, y], inward[:, z]
            condition = active[:, y] & active[:, z]
            length = h[:, y] + h[:, z]
            d = dist[y, z]
            # If both inward and one target is the line itself, the common
            # path has length one: the other target is exactly length-2 away.
            terminal = a & b & ((h[:, y] == 1) | (h[:, z] == 1))
            good &= ~(condition & terminal & (d != length - 2))
            same_point = a & b & ~terminal & (direction[:, y] == direction[:, z])
            good &= ~(condition & same_point & (d > length - 4))
            same_line = a & b & ~terminal & (direction[:, y] != direction[:, z])
            separate = a ^ b
            for category, steps in ((same_line, length - 2), (separate, length)):
                invalid = np.where(steps <= 8, d != steps, d < 16 - steps)
                good &= ~(condition & category & invalid)
    return good


def supported(pp, lines, lp):
    bits = 1 << np.arange(pp.shape[1], dtype=np.int64)
    keep = np.ones(len(pp), dtype=bool)
    for i, p in enumerate(pp):
        active = (p > 0) & (p < 8)
        incident = np.flatnonzero((lines == i).any(axis=1))
        masks = set(((lp[incident] < p) & active) @ bits)

        @lru_cache(None)
        def cover(remaining, slots):
            if not remaining:
                return slots == 0 or 0 in masks
            if not slots:
                return False
            bit = remaining & -remaining
            return any(cover(remaining ^ m, slots - 1) for m in masks if m & bit and m & remaining == m)

        keep[i] = cover(int(active @ bits), 5)
    return keep


def reduced():
    D, P, L, pl, ll, pp, lp = geometry()
    original_point_sizes = np.bincount(pl)
    original_line_sizes = np.bincount(ll)
    keep = line_rule(pp,L,lp)
    print("incidence line rule", int(keep.sum()), "of", len(L), flush=True)
    L, ll, lp = L[keep], ll[keep], lp[keep]
    keep = compatible(pp, L, lp, D)
    print("incidence compatible lines", int(keep.sum()), "of", len(L), flush=True)
    L, ll, lp = L[keep], ll[keep], lp[keep]
    while True:
        keep = supported(pp, L, lp)
        if keep.all():
            break
        print("unsupported point profiles", int((~keep).sum()), flush=True)
        lkeep = keep[L].all(axis=1)
        L, ll, lp = L[lkeep], ll[lkeep], lp[lkeep]
        L = (keep.cumsum() - 1)[L]
        P, pl, pp = P[keep], pl[keep], pp[keep]
    # The restrictions must delete entire orbits, not impose asymmetry.
    for labels, sizes in ((pl, original_point_sizes), (ll, original_line_sizes)):
        counts = np.bincount(labels, minlength=len(sizes))
        assert ((counts == 0) | (counts == sizes)).all()
    _, pl = np.unique(pl, return_inverse=True)
    _, ll = np.unique(ll, return_inverse=True)
    print("remaining points", len(P), len(set(pl)), "lines", len(L), len(set(ll)), flush=True)
    return D, P, L, pl, ll, pp, lp
