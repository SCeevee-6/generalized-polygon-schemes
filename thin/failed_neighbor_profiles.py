"""Check the four overlapping seeds at a failed pair and its three neighbors."""

from collections import defaultdict
from itertools import combinations, combinations_with_replacement

import numpy as np
from scipy.optimize import Bounds, LinearConstraint, milp
from scipy.sparse import coo_matrix

from failed_pair_types import on_slots
from grid_path_closure import LineConstraints, extend_paths
from opposite_link_surfaces import failed_metrics
from theta_incident_lines import unsupported
from theta_line_margins import quotient
from theta_profiles import GROUP, LINES, MIDDLE, expand, profiles, seed_distances


def neighboring_seeds(metric, lines, seed):
    # Equal slots in the failed-pair metric refer to one ambient point.
    rep = metric.argmin(axis=1)
    named = [tuple(rep[list(line)]) for line in lines]
    third = {frozenset((a, b)): c for line in named
             for a, b, c in (line, line[1:] + line[:1], line[2:] + line[:2])}
    seeds = [tuple(rep[:21])]
    for i, j in combinations(range(3), 2):
        slots = []
        for a in named[3+3*i]:
            b, = [x for x in named[3+3*j] if metric[a, x] == 3]
            path = [a]
            for h in (1, 2):
                internal = set(rep[(metric[a] == h) & (metric[b] == 3-h)])
                point, = internal
                path.append(point)
            path.append(b)
            slots.extend(path)
            slots.extend(third[frozenset(edge)] for edge in zip(path, path[1:]))
        assert np.array_equal(metric[np.ix_(slots, slots)], seed)
        seeds.append(tuple(slots))
    return seeds, named


def joined_profiles(metric, lines, seed, candidates):
    seeds, named = neighboring_seeds(metric, lines, seed)
    points = candidates
    ids = list(seeds[0])
    for slots in seeds[1:]:
        common = [i for i, x in enumerate(slots) if x in ids]
        extra = [i for i, x in enumerate(slots) if x not in ids]
        tails = defaultdict(list)
        for p in candidates:
            tails[p[common].tobytes()].append(p[extra])
        overlap = [ids.index(slots[i]) for i in common]
        points = np.array([np.concatenate((p, tail)) for p in points
                           for tail in tails[p[overlap].tobytes()]])
        ids.extend(slots[i] for i in extra)
    assert len(ids) == len(set(ids)) == len(set(metric.argmin(axis=1)))
    named = np.array(sorted({tuple(sorted(ids.index(i) for i in line)) for line in named}))
    d = metric[np.ix_(ids, ids)]
    mixed = 1 + 2*d[:, named].min(axis=2)
    dual = 1 + mixed[named].min(axis=1)
    np.fill_diagonal(dual, 0)
    full = np.block([[2*d, mixed], [mixed.T, dual]])
    observers = np.column_stack((2*points, 1+2*points[:, named].min(axis=2)))
    keep = np.ones(len(points), dtype=bool)
    for i, j in combinations(range(len(full)), 2):
        a, b = observers[:, i], observers[:, j]
        keep &= (abs(a-b) <= full[i, j]) & (a+b >= full[i, j])
    for i, j, k in combinations(range(len(full)), 3):
        a = observers[:, i] + full[j, k]
        b = observers[:, j] + full[i, k]
        c = observers[:, k] + full[i, j]
        # In the short four-point test the two largest sums must be equal.
        keep &= (((a+b >= 16) & (a+c >= 16) & (b+c >= 16))
                 | ((a == b) & (a >= c)) | ((a == c) & (a >= b)) | ((b == c) & (b >= a)))
    return points[keep], observers[keep], full, named


def binary_margins(observers, full, n):
    # Store the right-hand side in the lowest bit. An inconsistent row reduces to 1.
    pivots = {}

    def add(row, rhs):
        row = (row << 1) | int(rhs % 2)
        while row > 1:
            pivot = row.bit_length()
            if pivot not in pivots:
                pivots[pivot] = row
                return
            row ^= pivots[pivot]
        assert row == 0

    shells = [[int.from_bytes(np.packbits(col == h, bitorder="little").tobytes(), "little")
               for h in range(9)] for col in observers.T]
    radial = [quotient(0), quotient(1)]
    for i, j in combinations_with_replacement(range(len(full)), 2):
        kind = int(i >= n)
        for h in range(kind, 9, 2):
            for k in range(int(j >= n), 9, 2):
                add(shells[i][h] & shells[j][k], radial[kind][k, full[i, j], h])
    pair_rank = len(pivots)
    # The scaled least-eigenspace cosine is odd precisely at point distance four.
    for i, j, k in combinations_with_replacement(range(n), 3):
        add(shells[i][8] & shells[j][8] & shells[k][8], 0)
    return pair_rank, len(pivots)


def joint_stars(points, observers, full, named):
    n = points.shape[1]
    d = full[:n, :n] // 2
    lookup = {p.tobytes(): i for i, p in enumerate(points)}
    near = np.flatnonzero((points <= 1).any(axis=1))
    variable = dict(zip(near, range(len(near))))
    lower = [int((points[i] == 0).any()) for i in near]
    upper = [1 if fixed or (points[i] == 1).sum() >= 2 else 10 for i, fixed in zip(near, lower)]
    rows, cols, coeff, rhs = [], [], [], []

    def equation(terms, target):
        for j, value in terms.items():
            rows.append(len(rhs))
            cols.append(j)
            coeff.append(value)
        rhs.append(target)

    def new_variable(limit):
        lower.append(0)
        upper.append(limit)
        return len(lower)-1

    for x in range(n):
        neighbors = np.flatnonzero(points[:, x] == 1)
        places, inward = defaultdict(dict), defaultdict(dict)
        choices = {}
        line_variables = {}
        h = full[x]
        active = np.flatnonzero((h > 0) & (h < 8))
        for i in neighbors:
            if (abs(points[i]-d[x]) > 1).any():
                continue
            z = np.where(points[i] == d[x], points[i]-1, np.maximum(points[i], d[x]))
            j = lookup.get(z.tobytes())
            if j is None or j < i:
                continue
            values = np.array([h, observers[i], observers[j]])
            dist = 1 + values.min(axis=0)
            dist[n:][(values[:, n:] == 1).all(axis=0)] = 0
            toward = (values == dist-1).sum(axis=0)
            away = (values == dist+1).sum(axis=0)
            if not np.where(dist == 0, away == 3, np.where(dist == 8, toward == 3,
                            (toward == 1) & (away == 2))).all():
                continue
            v = new_variable(5)
            line_variables[(i, j)] = v
            choices[v] = places[i][v] = 1
            places[j][v] = 1 + int(i == j)
            for k in active[dist[active] < h[active]]:
                inward[k][v] = 1
        # Every already named line is an actual distinct line through x.
        # Identical profile triples still consume separate line slots.
        for line in named:
            if x not in line:
                continue
            others = list(line)
            others.remove(x)
            pair = tuple(sorted(lookup[d[y].tobytes()] for y in others))
            if pair not in line_variables:
                return None
            lower[line_variables[pair]] += 1
        equation(choices, 5)
        for i in neighbors:
            places[i][variable[i]] = -1
            equation(places[i], 0)
        for k in active:
            equation(inward[k], 1)

    for x, y in combinations(range(n), 2):
        if d[x, y] != 3:
            continue
        left = np.flatnonzero((points[:, x] == 1) & (points[:, y] == 2))
        right = np.flatnonzero((points[:, x] == 2) & (points[:, y] == 1))
        matches = defaultdict(dict)
        for i in left:
            for j in right:
                if (abs(points[i]-points[j]) > 1).any():
                    continue
                third = np.where(points[i] == points[j], points[i]-1, np.maximum(points[i], points[j]))
                if third.tobytes() in lookup:
                    v = new_variable(1)
                    matches[i][v] = matches[j][v] = 1
        for i in (*left, *right):
            matches[i][variable[i]] = -1
            equation(matches[i], 0)

    matrix = coo_matrix((coeff, (rows, cols)), shape=(len(rhs), len(lower))).tocsr()
    result = milp(np.zeros(len(lower)), integrality=1, bounds=Bounds(lower, upper),
                  constraints=LinearConstraint(matrix, rhs, rhs))
    assert result.x is not None, result.message
    counts = np.rint(result.x).astype(np.int64)
    assert (counts >= lower).all() and (counts <= upper).all()
    assert np.array_equal(matrix @ counts, rhs)
    return matrix.shape


def main():
    seed = seed_distances()
    candidates, _ = expand(profiles(seed))
    candidates = np.delete(candidates, unsupported(candidates), axis=0)
    _, lines, defs = extend_paths((1 << seed).astype(np.uint8), LINES, list(combinations(MIDDLE, 2)))
    perms = [on_slots(g, defs) for g in GROUP]
    metrics = failed_metrics(LineConstraints())
    seen = set()
    for metric in metrics:
        if metric.tobytes() in seen:
            continue
        orbit = {metric[np.ix_(p, p)].tobytes() for p in perms}
        seen.update(orbit)
        points, observers, full, named = joined_profiles(metric, lines, seed, candidates)
        ranks = binary_margins(observers, full, points.shape[1])
        shape = joint_stars(points, observers, full, named)
        print("orbit", len(orbit), "points", points.shape[1], "profiles", len(points),
              "binary ranks", ranks, "integral star system", shape, flush=True)
    assert seen == {m.tobytes() for m in metrics}


if __name__ == "__main__":
    main()
