"""An abstract rank-48 matrix satisfying the local block and count rules."""

from collections import Counter
from itertools import combinations, product

import numpy as np


def build_relaxation():
    # A slot records its branch, intermediate line, ancestor point, and
    # position among the four lines at that point.
    labels = list(product(range(2), range(4), range(2), range(4)))
    indices = {slot: i for i, slot in enumerate(labels)}
    location = np.asarray(labels)
    adjacency = np.zeros((64, 64), dtype=bool)
    for row, (branch, group, child, slot) in enumerate(labels):
        for other_branch in range(2):
            other_group = (group + slot) % 4
            flip = branch & other_branch
            col = indices[other_branch, other_group, child ^ flip, slot]
            adjacency[row, col] = True

    distance = np.full((64, 64), 4, dtype=int)
    for a, b in product(range(16), repeat=2):
        rows = np.arange(4*a, 4*a + 4)
        cols = np.arange(4*b, 4*b + 4)
        edges = np.argwhere(adjacency[np.ix_(rows, cols)])
        if len(edges):
            (r, c), = edges
            distance[rows[r], cols] = 3
            distance[rows, cols[c]] = 3
            distance[rows[r], cols[c]] = 2
        else:
            group_x = location[rows[0], :2]
            group_y = location[cols[0], :2]
            # The marked row and column already have a distance-three entry
            # through the other ancestor of their respective groups.
            in_group_y = np.all(location[:, :2] == group_y, axis=1)
            in_group_x = np.all(location[:, :2] == group_x, axis=1)
            remaining_rows = rows[~adjacency[rows][:, in_group_y].any(axis=1)]
            remaining_cols = cols[~adjacency[:, cols][in_group_x].any(axis=0)]
            assert len(remaining_rows) == len(remaining_cols) == 3
            distance[remaining_rows, remaining_cols] = 3
    return adjacency, distance


def verify(adjacency, distance):
    assert np.array_equal(adjacency, distance == 2)
    for matrix in (distance, distance.T):
        for row in matrix:
            assert Counter(row) == {2: 2, 3: 18, 4: 44}
            for branch in row.reshape(2, 32):
                assert Counter(branch) == {2: 1, 3: 9, 4: 22}
                for group in branch.reshape(4, 8):
                    expected = {2: 1, 3: 3, 4: 4} if 2 in group else {3: 2, 4: 6}
                    assert Counter(group) == expected
    for a, b in product(range(16), repeat=2):
        block = distance[4*a:4*a+4, 4*b:4*b+4]
        if (block == 2).any():
            (r, c), = np.argwhere(block == 2)
            assert (block[r, :] <= 3).all() and (block[:, c] <= 3).all()
            assert Counter(block.flat) == {2: 1, 3: 6, 4: 9}
        else:
            assert Counter(block.flat) == {3: 3, 4: 13}
            assert ((block == 3).sum(axis=0) <= 1).all()
            assert ((block == 3).sum(axis=1) <= 1).all()

    edges = [np.flatnonzero(row) + 64 for row in adjacency]
    edges += [np.flatnonzero(col) for col in adjacency.T]
    assert all(len(neighbors) == 2 for neighbors in edges)
    unseen = set(range(128))
    lengths = Counter()
    while unseen:
        pending = [unseen.pop()]
        size = 0
        while pending:
            v = pending.pop()
            size += 1
            neighbors = unseen.intersection(edges[v])
            unseen.difference_update(neighbors)
            pending.extend(neighbors)
        lengths[size] += 1
    assert lengths == {8: 16}
    return lengths

def expanded_metric(adjacency, cross):
    """Adjoin the two roots, their common line, and the two ancestor trees."""
    nodes = [("root", 0), ("root", 1), ("link",)]
    for side in range(2):
        nodes.extend(("parent", side, i, j) for i, j in product(range(2), range(4)))
        nodes.extend(("leaf", side, *coordinates)
                     for coordinates in product(range(2), range(4), range(2), range(4)))
    lines = len(nodes)
    nodes.extend(("point", *coordinates)
                 for coordinates in product(range(2), range(2), range(4), range(2)))
    depths = {"parent": 2, "point": 3, "leaf": 4}

    def leaf_index(node):
        _, _, branch, group, child, slot = node
        return 32*branch + 8*group + 4*child + slot

    def block(node):
        _, _, branch, group, *children = node
        stride = 8 if not children else 4
        start = 32*branch + 8*group + (4*children[0] if children else 0)
        return slice(start, start + stride)

    def distance(a, b):
        if a == b:
            return 0
        if a[0] == "link" or b[0] == "link":
            other = b if a[0] == "link" else a
            return 2 if other[0] == "root" else 2 + depths[other[0]]
        if "root" in (a[0], b[0]):
            if a[0] == b[0]:
                return 4
            root, other = (a, b) if a[0] == "root" else (b, a)
            return depths[other[0]] + (0 if root[1] == other[1] else 4)
        if a[1] == b[1]:
            common = 0
            for x, y in zip(a[2:], b[2:]):
                if x != y:
                    break
                common += 1
            return depths[a[0]] + depths[b[0]] - 2*common
        if a[1] == 1:
            a, b = b, a
        if a[0] == b[0] == "leaf":
            return 2 * cross[leaf_index(a), leaf_index(b)]
        if "leaf" not in (a[0], b[0]):
            if a[0] == b[0] == "point":
                return 6 if adjacency[block(a), block(b)].any() else 8
            return 8 if a[0] == b[0] else 7
        row = leaf_index(a) if a[0] == "leaf" else block(a)
        col = leaf_index(b) if b[0] == "leaf" else block(b)
        ancestor = b if a[0] == "leaf" else a
        if ancestor[0] == "point":
            return 5 if adjacency[row, col].any() else 7
        return 6 if adjacency[row, col].any() else 8

    return np.asarray([[distance(a, b) for b in nodes] for a in nodes], dtype=int), lines


def verify_metric(distance, line_count):
    """Check exact necessary conditions for the mixed incidence distance table."""
    from grid_directions import psd_rank
    from theta_line_margins import quotient

    order = len(distance)
    assert np.array_equal(distance, distance.T)
    for row in distance:
        assert (distance <= row[:, None] + row).all()
    for i in range(order - 3):
        j, k, l = np.array(list(combinations(range(i + 1, order), 3))).T
        x = distance[i, j] + distance[k, l]
        y = distance[i, k] + distance[j, l]
        z = distance[i, l] + distance[j, k]
        largest = np.maximum(x, np.maximum(y, z))
        smallest = np.minimum(x, np.minimum(y, z))
        total = x + y + z
        assert ((total - largest >= 16) | (total - smallest == 2*largest)).all()
    tables = [quotient(1), quotient(0)]
    for i, j in combinations(range(order), 2):
        expected = tables[int(i >= line_count)][:, distance[i, j], :]
        for section in (slice(None, line_count), slice(line_count, None)):
            counts = np.bincount(9*distance[j, section] + distance[i, section],
                                 minlength=81).reshape(9, 9)
            assert (counts <= expected).all()
    ll = distance[:line_count, :line_count] // 2
    pp = distance[line_count:, line_count:] // 2
    pl = (distance[line_count:, :line_count] - 1) // 2
    for i, j in combinations(range(line_count), 2):
        if ll[i, j] != 4:
            continue
        middle = np.flatnonzero((ll[i] == 2) & (ll[j] == 2))
        assert len(middle) <= 3 and all(ll[x, y] == 4 for x, y in combinations(middle, 2))
        if len(middle) == 3:
            completes = (ll[i] == 4) & (ll[j] == 4) & (ll[middle] == 2).all(axis=0)
            assert not completes.any()
    point_kernels = [
        [64, 32, 8, 0, -1],
        [320, 32, -40, -4, 5],
        [320, -96, 8, 8, -5],
        [16, -8, 4, -2, 1],
    ]
    line_kernels = [
        [192, 112, 32, 2, -3],
        [64, 16, -8, -2, 1],
        [192, -16, -16, 10, -3],
        [256, -64, 16, -4, 1],
    ]
    ranks = []
    for point, line in zip(point_kernels, line_kernels):
        p, l = np.array(point), np.array(line)
        cross = l[0] * (p[pl] + 2*p[pl + 1])
        gram = np.block([[3*(p[0] + 2*p[1])*l[ll], cross.T],
                         [cross, l[0]*p[pp]]])
        ranks.append(psd_rank(gram))
    assert ranks == [128, 131, 131, 20]
    assert psd_rank(np.array(line_kernels[-1])[ll]) == 111
    print("Mixed projector ranks", ranks, "line kernel rank 111")


if __name__ == "__main__":
    adjacency, cross = build_relaxation()
    lengths = verify(adjacency, cross)
    distance, line_count = expanded_metric(adjacency, cross)
    assert distance.shape == (179, 179) and line_count == 147
    verify_metric(distance, line_count)
    print("Cycle lengths", dict(lengths))
    print("Rank over any field", 64 - sum(lengths.values()))
