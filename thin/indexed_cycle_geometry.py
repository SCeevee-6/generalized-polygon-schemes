"""Check which eight-cycle attachment types pass the named-line constraints."""

from collections import Counter
from itertools import combinations, product

import numpy as np

from grid_path_closure import LineConstraints
from indexed_cycle_blocks import build_relaxation


def line_distances(a, b, matches_at_leaf, blocks):
    """Necessary line distances in the two rooted families."""
    depth = {"parent": 1, "leaf": 2}
    if a == b:
        return (0,)
    if ("link",) in (a, b):
        other = b if a == ("link",) else a
        return (1 if other[0] == "root" else depth[other[0]] + 1,)
    if "root" in (a[0], b[0]):
        if a[0] == b[0]:
            return (2,)
        root, other = (a, b) if a[0] == "root" else (b, a)
        return (depth[other[0]] + 2*(root[1] != other[1]),)
    if a[1] == b[1]:
        prefix = 0
        for x, y in zip(a[2:], b[2:]):
            if x != y:
                break
            prefix += 1
        return (depth[a[0]] + depth[b[0]] - prefix,)
    if a[1] == 1:
        a, b = b, a
    if a[0] == b[0] == "parent":
        return (4,)
    if a[0] == b[0] == "leaf":
        row, col = a[1:], b[1:]
        fixed = blocks.get((row[:4], col[:4]))
        if fixed is not None:
            return (4 - int(row == fixed[0]) - int(col == fixed[1]),)
        return (3, 4) if row in matches_at_leaf or col in matches_at_leaf else (2, 3, 4)
    leaf, parent = (a[1:], b[1:]) if a[0] == "leaf" else (b[1:], a[1:])
    if leaf not in matches_at_leaf:
        return (3, 4)
    return (3 if any(other[:3] == parent for other in matches_at_leaf[leaf]) else 4,)


def replay(groups, edges, constraints):
    # A leaf is named by side, branch, parent line, ancestor point and slot.
    # Every leaf occurring in edges has both of its matches supplied.
    leaves = [group + (child, slot) for group in groups
              for child, slot in product(range(2), range(4))]
    named = [("root", 0), ("root", 1), ("link",)]
    named += [("parent", *group) for group in groups]
    named += [("leaf", *leaf) for leaf in leaves]
    blocks = {}
    matches = {}
    for a, b in edges:
        block = a[:4], b[:4]
        if block in blocks:
            return "block", None
        blocks[block] = a, b
        matches.setdefault(a, set()).add(b)
        matches.setdefault(b, set()).add(a)
    assert all({other[1] for other in partners} == {0, 1} and len(partners) == 2
               for partners in matches.values())

    point_ids = {}

    def index(label):
        return point_ids.setdefault(label, len(point_ids))

    def points(line):
        kind, *label = line
        if kind == "root":
            return [("branch", label[0], branch) for branch in range(3)]
        if kind == "link":
            return [("branch", 0, 2), ("branch", 1, 2), ("link-tip",)]
        if kind == "parent":
            return [("branch", *label[:2]), ("point", *label, 0), ("point", *label, 1)]
        return [("point", *label[:4]), ("outer", *label, 0), ("outer", *label, 1)]

    lines = [tuple(index(point) for point in points(line)) for line in named]
    for a, b in edges:
        lines.append(tuple(index(point) for point in
                           (("outer", *a, b[1]), ("outer", *b, a[1]), ("extra", a, b))))
    masks = np.full((len(point_ids),)*2, 31, dtype=np.uint8)
    np.fill_diagonal(masks, 1)
    for line in lines:
        for a, b in combinations(line, 2):
            masks[a, b] = masks[b, a] = 2

    squares = constraints.squares
    distance = 1 + squares.min(axis=(1, 2))
    distance[(squares == 0).sum(axis=(1, 2)) == 3] = 0
    types = [int(constraints.bits[distance == d].sum()) for d in range(5)]
    domains = np.empty((len(named), len(named)), dtype=np.uint64)
    for i, a in enumerate(named):
        for j, b in enumerate(named[:i+1]):
            allowed = line_distances(a, b, matches, blocks)
            domains[i, j] = domains[j, i] = sum(types[d] for d in allowed)
    if not constraints.propagate(masks, lines, line_domains=domains):
        return "paths", None
    return "open", (masks, point_ids, lines, named)


def cycle(relations):
    """Name an eight-cycle with alternating branches on each side."""
    sides = []
    for side in range(2):
        first, second = [], []
        for branch in range(2):
            relation = relations[2*side+branch]
            first.append((side, branch, 0, 0, 0))
            second.append((side, branch, int(relation == 0), int(relation == 1), int(relation == 2)))
        sides.append(first + second)
    rows, cols = sides
    edges = sorted((rows[i], cols[j % 4]) for i in range(4) for j in (i, i - 1))
    groups = sorted({leaf[:3] for leaf in rows + cols})
    return groups, edges


def threshold_witness(constraints):
    matching, _ = build_relaxation()
    slots = list(product(range(2), range(4), range(2), range(4)))
    groups = list(product(range(2), range(2), range(4)))
    edges = [((0, *slots[i]), (1, *slots[j])) for i, j in np.argwhere(matching)]
    result, data = replay(groups, edges, constraints)
    assert result == "open"
    masks, ids, lines, _ = data
    center = ids["point", 0, 0, 0, 0]
    witness = ids["extra", (0, 1, 0, 0, 0), (1, 1, 0, 1, 0)]
    star = [line for line in lines if center in line]
    neighbors = [p for line in star for p in line if p != center]
    assert len(star) == 5 and len(set(neighbors)) == 10
    assert not np.any(masks[np.ix_(neighbors, neighbors)] & 1 & ~np.eye(10, dtype=bool))
    assert masks[center, witness] == 8
    assert np.all(masks[witness, neighbors] == 16)
    print("Rank-48 matching table excluded at a full five-line star")


def main():
    constraints = LineConstraints()
    counts = Counter()
    survivors = set()
    for relation in product(range(3), repeat=4):
        result, _ = replay(*cycle(relation), constraints)
        counts[result] += 1
        if result == "open":
            survivors.add(relation)
    expected = {(0, 0, 0, 0), (1, 1, 1, 1)}
    expected.update(row for row in product(range(2), repeat=4)
                    if sum(row) == 1 or (sum(row[:2]) == sum(row[2:]) == 1))
    assert counts == {"block": 25, "paths": 46, "open": 10} and survivors == expected
    print("Eight-cycle attachment cases", dict(counts))
    threshold_witness(constraints)


if __name__ == "__main__":
    main()
