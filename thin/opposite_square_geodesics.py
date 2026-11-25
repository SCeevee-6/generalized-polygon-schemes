"""Replay the six label-4 square tests at the roots of each failed metric."""

from itertools import combinations, permutations

import numpy as np

from failed_pair_types import on_slots
from grid_path_closure import LineConstraints, extend_paths
from opposite_link_surfaces import failed_metrics
from theta_profiles import GROUP, LINES, MIDDLE, seed_distances


def path_lines(constraints, masks, lines, requests):
    masks, lines, _ = extend_paths(masks, lines, requests)
    assert constraints.propagate(masks, lines)
    paths = []
    for a, b in requests:
        path = [a]
        for k in (1, 2):
            candidates = np.flatnonzero((masks[a] == 1 << k) & (masks[b] == 1 << (3-k)))
            path.append(int(candidates[0]))
        path.append(b)
        segments = []
        for x, y in zip(path, path[1:]):
            line = next(line for line in lines if x in line and y in line)
            spare, = set(line) - {x, y}
            segments.append((line, spare))
        paths.append(segments)
    return masks, lines, paths


def replay(constraints, metric, lines):
    colors = list(permutations(range(3), 2))
    first = [(7*i+4, 7*j+6) for i, j in colors]
    masks, lines, paths = path_lines(constraints, (1 << metric).astype(np.uint8), lines, first)
    second = []
    for (i, j), path in zip(colors, paths):
        k = 3-i-j
        second.extend([(path[0][1], 7*k+3), (path[2][1], 7*k)])
    assert all(masks[a, b] == 8 for a, b in second)
    masks, lines, continuations = path_lines(constraints, masks, lines, second)

    n = len(masks)
    enlarged = np.full((n+6, n+6), 31, dtype=np.uint8)
    enlarged[:n, :n] = masks
    np.fill_diagonal(enlarged, 1)
    squares = []
    for pos, path in enumerate(paths):
        (e, p), (b, q) = [arm[1] for arm in continuations[2*pos:2*pos+2]]
        connector = (p, q, n+pos)
        lines.append(connector)
        for u, v in combinations(connector, 2):
            enlarged[u, v] &= 2
            enlarged[v, u] &= 2
        squares.append((path[0][0], b, path[2][0], e))
    # Each square has four different auxiliary vertices, disjoint from the
    # corners of every other square. Exclude coincident line patterns.
    corners = [lines.index(line) for square in squares for line in square]
    assert len(set(corners)) == 24
    distinct = sum(1 << i for i, s in enumerate(constraints.squares) if (s == 0).sum() != 3)
    domains = np.full((len(lines),)*2, (1 << 39)-1, dtype=np.uint64)
    for i, j in combinations(corners, 2):
        domains[i, j] = domains[j, i] = distinct
    assert constraints.propagate(enlarged, lines, line_domains=domains)

    for a, b, d, e in squares:
        requests = []
        for left, right, centers in [(a, b, (d, e, LINES[0])), (d, e, (a, b, LINES[1]))]:
            tips = []
            for line in centers:
                used = {np.flatnonzero((enlarged[np.ix_(line, end)] == 2).any(axis=1)).item()
                        for end in (left, right)}
                spare, = set(range(3)) - used
                tips.append(line[spare])
            requests.extend(combinations(tips, 2))
        for x, y in requests:
            enlarged[x, y] &= 8
            enlarged[y, x] &= 8
        enlarged, lines, _ = extend_paths(enlarged, lines, requests)
        assert constraints.propagate(enlarged, lines, line_domains=domains)
    return len(lines)


def main():
    constraints = LineConstraints()
    metrics = failed_metrics(constraints)
    _, lines, defs = extend_paths((1 << seed_distances()).astype(np.uint8), LINES,
                                  list(combinations(MIDDLE, 2)))
    group = [on_slots(g, defs) for g in GROUP]
    seen = set()
    sizes = []
    for metric in metrics:
        if metric.tobytes() in seen:
            continue
        orbit = {metric[np.ix_(p, p)].tobytes() for p in group}
        seen.update(orbit)
        sizes.append(len(orbit))
        print("orbit", len(orbit), "line slots after six square closures", replay(constraints, metric, lines))
    assert sorted(sizes) == [2, 6, 6, 6, 12] and seen == {d.tobytes() for d in metrics}


if __name__ == "__main__":
    main()
