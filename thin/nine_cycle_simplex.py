"""Check the centered triples and first lines of the nine-cycle's nine tips."""
from collections import Counter
from itertools import combinations, product

import numpy as np

from grid_path_closure import LineConstraints, extend_paths
from opposite_cycles import traces
from opposite_nine_cycles import attachment_choices, named_lines


def simplex_seed(constraints, distances, attachments):
    masks, lines = named_lines(distances, attachments)
    assert constraints.propagate(masks, lines)
    assert ((masks & (masks-1)) == 0).all()
    # Identify only known coincidences, before naming new paths.
    classes, index = np.unique(masks.argmin(axis=1), return_inverse=True)
    lines = sorted({tuple(sorted(index[list(line)])) for line in lines})
    masks = masks[np.ix_(classes, classes)]
    degree = Counter(p for line in lines for p in line)
    tips = sorted(p for p, d in degree.items() if d == 1)
    assert len(classes) == 36 and len(lines) == 21 and len(tips) == 9
    assert all(masks[x, y] == 8 for x, y in combinations(tips, 2))
    masks, lines, _ = extend_paths(masks, lines, list(combinations(tips, 2)))
    assert constraints.propagate(masks, lines)
    first = {}
    for x, y in product(tips, repeat=2):
        if x != y:
            slots = np.flatnonzero((masks[x] == 2) & (masks[y] == 4))
            assert len(slots) and (masks[slots[0], slots] == 1).all()
            first[x, y] = int(slots[0])
    return masks, lines, tips, first


def possible_centers(constraints, masks, lines, tips, first):
    candidates = []
    tested = 0
    for x, y, z in combinations(tips, 3):
        # In a centered triple, from each tip both paths go through
        # the same projection point on the center line.
        equal = [(first[a, b], first[a, c]) for a, b, c in
                 ((x, y, z), (y, z, x), (z, x, y))]
        if any(not masks[p, q] & 1 for p, q in equal):
            continue
        tested += 1
        trial = masks.copy()
        for p, q in equal:
            trial[p, q] = trial[q, p] = 1
        if constraints.propagate(trial, lines):
            candidates.append((x, y, z))
    assert tested == 30 and len(candidates) == 3
    assert sorted(p for triple in candidates for p in triple) == tips
    return candidates


def unused_lines(masks, lines, tips):
    start = len(masks)
    big = np.full((start+18, start+18), 31, dtype=np.uint8)
    big[:start, :start] = masks
    np.fill_diagonal(big, 1)
    unused = []
    for i, x in enumerate(tips):
        a, b = start+2*i, start+2*i+1
        unused.append((x, a, b))
        big[a, b] = big[b, a] = 2
        for t in tips:
            for p in (a, b):
                big[p, t] = big[t, p] = 2 if t == x else 16
    for a, b in combinations(unused, 2):
        big[np.ix_(a, b)] &= 24
        big[np.ix_(b, a)] &= 24
    return big, lines+unused


def symmetric_witness(constraints, masks, lines, tips, first):
    # The two tips reached through the original line at each tip make a
    # nine-cycle. Find it from distances to the original 36 points.
    adjacent = {x: [y for y in tips if y != x and first[x, y] < 36] for x in tips}
    assert all(len(row) == 2 for row in adjacent.values())
    cycle = [tips[0], min(adjacent[tips[0]])]
    while len(cycle) < 9:
        y, = set(adjacent[cycle[-1]]) - {cycle[-2]}
        cycle.append(y)
    assert len(set(cycle)) == 9 and cycle[0] in adjacent[cycle[-1]]
    # The first line pairs offsets -1,+1. This candidate pairs the six
    # other targets on three more lines. It is a compatibility test.
    trial, named = unused_lines(masks, lines, tips)
    for i, x in enumerate(cycle):
        for a, b in ((1, 8), (2, 3), (4, 5), (6, 7)):
            y, z = cycle[(i+a) % 9], cycle[(i+b) % 9]
            for target, partner in ((y, z), (z, y)):
                p = first[x, target]
                for t in tips:
                    d = 1 if t == x else 2 if t == target else 3 if t == partner else 4
                    trial[p, t] &= 1 << d
                    trial[t, p] &= 1 << d
    assert constraints.propagate(trial, named)


def main():
    constraints = LineConstraints()
    count = 0
    for _, distances in traces(9):
        for attachments in product(*(attachment_choices(i, distances) for i in range(9))):
            masks, lines, tips, first = simplex_seed(constraints, distances, attachments)
            possible_centers(constraints, masks, lines, tips, first)
            symmetric_witness(constraints, masks, lines, tips, first)
            count += 1
    assert count == 2
    print("Both nine-cycle traces: 84 triples, three disjoint candidates,"
          " all excluded by even degree; joint first-line witness passes")


if __name__ == "__main__":
    main()
