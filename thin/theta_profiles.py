"""Check point-distance marginals on the 21 points forced by two opposite lines.

python thin/theta_profiles.py              # verify the shipped integer certificate
python thin/theta_profiles.py --solve      # find another with scipy.optimize.milp
"""

import argparse
import json
from itertools import combinations, combinations_with_replacement, permutations, product
from pathlib import Path

import numpy as np
from scipy.optimize import Bounds, LinearConstraint, milp
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import shortest_path

LINES = [(0, 7, 14), (3, 10, 17)]
for a in (0, 7, 14):
    LINES.extend((a + i, a + i + 1, a + 4 + i) for i in range(3))
ENDS = [4, 6, 11, 13, 18, 20]
MIDDLE = [5, 12, 19]
COS = np.array([16, -8, 4, -2, 1], dtype=np.int64)
GROUP = [
    tuple(7 * leg + q for leg in perm for q in order)
    for perm in permutations(range(3))
    for order in (range(7), (3, 2, 1, 0, 6, 5, 4))
]


def seed_distances():
    # Core point paths of at most four steps in the three-arm diagram are geodesics.
    adj = np.zeros((12, 12), dtype=int)
    for a in range(3):
        for j in range(3):
            u = 4 * a + j
            adj[u, u + 1] = adj[u + 1, u] = 1
        for b in range(a):
            for j in (0, 3):
                adj[4 * a + j, 4 * b + j] = adj[4 * b + j, 4 * a + j] = 1
    core = shortest_path(csr_matrix(adj), unweighted=True).astype(int)
    assert core.max() == 4
    coeff = np.zeros((21, 12), dtype=int)
    for i in range(3):
        for j in range(4):
            coeff[7 * i + j, 4 * i + j] = 1
        for j in range(3):
            coeff[7 * i + 4 + j] = -coeff[7 * i + j] - coeff[7 * i + j + 1]
    assert np.array_equal(coeff[ENDS].sum(0) @ COS[core], coeff[MIDDLE].sum(0) @ COS[core])
    gram = coeff @ COS[core] @ coeff.T
    decode = {value: d for d, value in enumerate(COS)}
    return np.array([[decode[v] for v in row] for row in gram])


def third(a, b):
    return a - 1 if a == b else max(a, b)


def profiles(dist):
    end_types = sorted({t for v in range(4) for t in permutations((v, v + 1, v + 1))})
    arms = {}
    for a, b in product(range(5), repeat=2):
        arms[a, b] = []
        for l, r in product(range(5), repeat=2):
            chain = (a, l, r, b)
            if any(abs(x - y) > 1 or (x == y == 0) for x, y in zip(chain, chain[1:])):
                continue
            arms[a, b].append(chain + tuple(third(x, y) for x, y in zip(chain, chain[1:])))
    representatives = set()
    for left, right in product(end_types, repeat=2):
        for pieces in product(*(arms[a, b] for a, b in zip(left, right))):
            profile = sum(pieces, ())
            if all(profile[t] == 3 for t in ENDS):
                continue  # Hypothesis under test: no seventh point.
            representatives.add(min(tuple(profile[j] for j in g) for g in GROUP))
    P = np.array(sorted(representatives), dtype=np.int8)
    ok = np.ones(len(P), dtype=bool)
    for i, j in combinations(range(21), 2):
        ok &= (np.abs(P[:, i] - P[:, j]) <= dist[i, j]) & (P[:, i] + P[:, j] >= dist[i, j])
    # On a point at distances 1 to 3 there is just one inward line toward the observer.
    for i in range(21):
        inward = sum((P[:, line].min(1) < P[:, i]).astype(int) for line in LINES if i in line)
        ok &= (P[:, i] == 4) | (inward <= 1)
    return P[ok]


def subset_orbits(k):
    seen = set()
    for tup in combinations(range(21), k):
        if tup not in seen:
            orbit = {tuple(sorted(g[i] for i in tup)) for g in GROUP}
            seen.update(orbit)
            yield sorted(orbit)


def distance_polynomials():
    q = np.diag([0, 1, 1, 1, 5]) + np.diag([10, 8, 8, 8], 1) + np.diag([1, 1, 1, 5], -1)
    polys = [np.eye(5, dtype=int), q]
    for i in (1, 2, 3):
        m = (q - np.eye(5, dtype=int)) @ polys[-1] - (10 if i == 1 else 8) * polys[-2]
        c = 5 if i == 3 else 1
        assert np.all(m % c == 0)
        polys.append(m // c)
    return polys


def expand(reps):
    named, labels = [], []
    for i, rep in enumerate(reps):
        orbit = sorted({tuple(rep[list(g)]) for g in GROUP})
        named.extend(orbit)
        labels.extend([i] * len(orbit))
    return np.asarray(named, dtype=np.int8), np.asarray(labels)


def solve(D, P, labels):
    # One unknown per orbit is the common multiplicity of EACH named profile.
    # This symmetry is only an ansatz for constructing a feasible witness.
    n = int(labels.max()) + 1
    rows, rhs = [np.bincount(labels)], [1755]
    polys = distance_polynomials()
    seen = set()
    for i, j in permutations(range(21), 2):
        if (i, j) in seen:
            continue
        seen.update((g[i], g[j]) for g in GROUP)
        for h, k in product(range(5), repeat=2):
            mask = (P[:, i] == h) & (P[:, j] == k)
            rows.append(np.bincount(labels[mask], minlength=n))
            rhs.append(polys[h][D[i, j], k])
    scaled = COS[P]
    for orbit in subset_orbits(3):
        i, j, k = orbit[0]
        row = np.zeros(n, dtype=np.int64)
        np.add.at(row, labels, scaled[:, i] * scaled[:, j] * scaled[:, k])
        rows.append(row)
        rhs.append(0)
    A = csr_matrix(rows)
    result = milp(
        np.zeros(n),
        integrality=np.ones(n),
        bounds=Bounds(0, 1755),
        constraints=LinearConstraint(A, rhs, rhs),
        options={"time_limit": 180},
    )
    if result.x is None:
        raise RuntimeError(result.message)
    return np.rint(result.x).astype(np.int64)


def verify(D, P, weights):
    # Check every named pair and triple without orbit reduction.
    polys = distance_polynomials()
    assert (weights >= 0).all() and weights.sum() == 1755
    assert not np.any(np.all(P[:, ENDS] == 3, axis=1))
    for i, j in product(range(21), repeat=2):
        counts = np.zeros((5, 5), dtype=np.int64)
        np.add.at(counts, (P[:, i], P[:, j]), weights)
        assert np.array_equal(counts, [poly[D[i, j]] for poly in polys])
    scaled = COS[P]
    for i, j in combinations_with_replacement(range(21), 2):
        assert not np.any((weights * scaled[:, i] * scaled[:, j]) @ scaled)
    for row in D:
        assert weights[(P == row).all(axis=1)].sum() == 1


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--solve", action="store_true")
    parser.add_argument(
        "--certificate", type=Path, default=Path(__file__).with_name("theta_no_extension.json")
    )
    args = parser.parse_args()
    D = seed_distances()
    reps = profiles(D)
    P, labels = expand(reps)
    if args.solve:
        x = solve(D, P, labels)
    else:
        index = {tuple(row): j for j, row in enumerate(reps.tolist())}
        x = np.zeros(len(reps), dtype=np.int64)
        for row, count in json.loads(args.certificate.read_text()):
            assert type(count) is int and 0 < count <= 1755
            j = index[tuple(row)]
            assert x[j] == 0
            x[j] = count
    weights = x[labels]
    verify(D, P, weights)
    print("admissible named profiles", len(P), "orbits", len(reps))
    print("occupied named profiles", np.count_nonzero(weights), "points", weights.sum())
    print("ordered pair tables 441, cubic entries 1771, exact residual 0")
    print("points at distance three from all six tips", 0)


if __name__ == "__main__":
    main()
