"""Check forced first lines of geodesics in partial 24-cell assignments."""
from itertools import combinations
import numpy as np

from grid_path_closure import LineConstraints
from partial_frame_metrics import Classes, Reference, metric_seed


def determined(masks):
    d = np.full(masks.shape, -1, dtype=np.int8)
    for h in range(5):
        d[masks == 1 << h] = h
    return d


def forced_directions(d, p):
    """Return equality classes of first lines for known targets at distances 1..3."""
    h = d[p]
    targets = np.flatnonzero((h > 0) & (h < 4))
    a = h[targets, None]
    b = a.T
    k = d[np.ix_(targets, targets)]
    classes = Classes(len(targets))
    equal = (k >= 0) & (a+b+k < 8) & (a+b > k)
    for x, y in zip(*np.nonzero(np.tril(equal, -1))):
        classes.merge(x, y)
    names = [classes.rep(i) for i in range(len(targets))]
    graph = {i: set() for i in names}
    distinct = (k == a+b)
    # At perimeter eight, sharing a first edge still forces a tripod.
    borderline = (k >= 2) & (k <= 3) & (a+b+k == 8) & (a+b > k)
    for i, j in zip(*np.nonzero(np.tril(borderline, -1))):
        x, y = targets[i], targets[j]
        interior = (d[x] > 0) & (d[y] > 0) & (d[x]+d[y] == k[i, j]) & (h >= 0)
        prediction = h[x]+h[y]-k[i, j] + np.abs(h[x]-h[y]+k[i, j]-2*d[x, interior])
        if np.any(2*h[interior] != prediction):
            distinct[i, j] = True
    for x, y in zip(*np.nonzero(np.tril(distinct, -1))):
        i, j = names[x], names[y]
        if i == j:
            return None
        graph[i].add(j)
        graph[j].add(i)

    # Four pairwise distinct forbidden first lines leave just the fifth.
    def four_clique(vertices):
        return any(graph[x] & graph[y] & graph[z] & vertices
                   for x in vertices for y in graph[x] & vertices if x < y
                   for z in graph[x] & graph[y] & vertices if y < z)

    for _ in range(len(graph)):
        match = next(((i, j) for i, j in combinations(sorted(graph), 2)
                      if j not in graph[i] and len(graph[i] & graph[j]) >= 4
                      and four_clique(graph[i] & graph[j])), None)
        if match is None:
            break
        i, j = match
        graph[i].update(graph.pop(j))
        for v in graph[i]:
            graph[v].discard(j)
            graph[v].add(i)
        names = [i if label == j else label for label in names]
    return dict(zip(map(int, targets), names))


def rectangle(d, witness):
    p, r, q, s = witness
    assert d[p, r] == d[q, s] == 4
    for x, y, z in ((p, q, s), (r, q, s), (q, p, r), (s, p, r)):
        assert d[x, y] == d[x, z] == 3
        directions = forced_directions(d, x)
        if directions is None:
            return  # Contradictory first-line equalities already exclude this seed.
        assert directions[y] == directions[z]


def find_rectangle(d, stop):
    # The verifier can check just four roots once their slots have been saved.
    same = {}
    for p in range(stop):
        directions = forced_directions(d, p)
        if directions is None:
            return [p]
        pairs = set()
        by_label = {}
        for q, label in directions.items():
            if q < stop and d[p, q] == 3:
                by_label.setdefault(label, []).append(q)
        for row in by_label.values():
            pairs.update((q, s) for q, s in combinations(row, 2) if d[q, s] == 4)
        same[p] = pairs
    for p, pairs in same.items():
        for q, s in pairs:
            for u, v in same[q] & same[s]:
                if p in (u, v) and (q, s) in same[u+v-p]:
                    return [p, u+v-p, q, s]
    return None


def model_control():
    from check_opposite_seed import ROOT, octagon

    pd, lines = octagon(np.load(ROOT / "data/generators.npy"))
    # Include all ten neighbors of each tested root, so every true first
    # line has a target. Other named targets and interior points are sampled.
    subset = set(range(0, len(pd), 13))
    roots = [0, 117, 754]
    for p in roots:
        subset.update(map(int, np.flatnonzero(pd[p] == 1)))
    subset = sorted(subset)
    d = pd[np.ix_(subset, subset)]
    for p in roots:
        through = lines[(lines == p).any(axis=1)]
        assert len(through) == 5
        directions = forced_directions(d, subset.index(p))
        assert directions is not None
        actual = {}
        for pos, label in directions.items():
            q = subset[pos]
            which, = np.flatnonzero(pd[q, through].min(axis=1) == pd[p, q] - 1)
            assert actual.setdefault(label, which) == which
        assert len(set(actual.values())) == 5
    print("First-line equalities agree with the Ree–Tits model at three sampled roots")


def twelve_controls(constraints):
    from indexed_twelve_metrics import WITNESSES, check_case
    from opposite_link_surfaces import failed_metrics

    metrics = (1 << failed_metrics(constraints)).astype(np.uint8)
    for word, choices in WITNESSES.items():
        masks, _ = check_case(word, choices, metrics, constraints)
        assert find_rectangle(determined(masks), len(masks)) is None
    print("All thirteen twelve-cycle witnesses pass the first-line rectangle test")


def witnesses():
    from pathlib import Path
    rows = Path(__file__).with_name("frame_direction_rectangles.txt").read_text().splitlines()
    parsed = [list(map(int, row.split())) for row in rows if row and not row.startswith("#")]
    assert all(len(row) == 6 for row in parsed) and len(parsed) == 96
    return {tuple(row[:2]): row[2:] for row in parsed}


def verify(ref, constraints, omitted, witness):
    masks, lines = metric_seed(ref, omitted)
    assert constraints.propagate(masks, lines)
    rectangle(determined(masks), witness)


def main():
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pair", nargs=2, type=int, metavar=("X", "Y"))
    parser.add_argument("--model-control", action="store_true")
    parser.add_argument("--twelve", action="store_true", help="test existing twelve-cycle witnesses")
    args = parser.parse_args()
    if args.model_control:
        model_control()
        return
    if args.twelve:
        twelve_controls(LineConstraints())
        return
    ref = Reference()
    certs = witnesses()
    adjacent = {p for p in combinations(range(24), 2)
                if len(ref.incident[p[0]] & ref.incident[p[1]]) == 3}
    assert set(certs) == adjacent
    constraints = LineConstraints()
    if args.pair:
        omitted = tuple(sorted(args.pair))
        verify(ref, constraints, omitted, certs[omitted])
        print("Verified rectangle for", omitted)
        return
    for x in range(24):
        assert all(ref.known({x}))
    for omitted in combinations(range(24), 2):
        if omitted in adjacent:
            verify(ref, constraints, omitted, certs[omitted])
        else:
            assert all(ref.known(omitted))
    print("All 24 single omissions and 276 double omissions excluded;",
          "96 adjacent pairs checked by first-line rectangles")
    pairs = 2925 * 2048 // 2
    assert 6*pairs//24 == 748800 and 748800//128 == 5850


if __name__ == "__main__":
    main()
