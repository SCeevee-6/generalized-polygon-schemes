"""Check the attachment choices for an eight-face surface component."""

from itertools import combinations, product

import numpy as np

from grid_path_closure import LineConstraints


def replay(relations, constraints):
    triples = list(product(range(2), repeat=3))
    root = ("root",)
    middle = [("middle", i, s) for i, s in product(range(3), range(2))]
    partners = [("partner", *t) for t in triples]
    point_ids = {}

    def index(label):
        return point_ids.setdefault(label, len(point_ids))

    named = {line: tuple(index((line, k)) for k in range(3))
             for line in (root, *partners)}
    connectors = []
    close = {}
    for i, shared in enumerate(relations):
        parents = [index(("parent", i, s)) for s in (0, 1)]
        for s in (0, 1):
            line = ("middle", i, s)
            named[line] = (parents[s], index((line, 1)), index((line, 2)))
            close[root, line] = i, 0
        if shared:
            connectors.append((named[root][i], *parents))
        else:
            for s in (0, 1):
                connectors.append((named[root][i], parents[s], index(("tip", i, s))))

    for t in triples:
        b = ("partner", *t)
        for i, s in enumerate(t):
            c = ("middle", i, s)
            # Opposite neighbors of a cube edge force a checkerboard
            # on the two available attachment points of each middle line.
            j = 1 + (sum(t) - s) % 2
            close[c, b] = j, i
            connectors.append((named[c][j], named[b][i], index(("tip", i, *t))))

    lines = list(named.values()) + connectors
    masks = np.full((len(point_ids),)*2, 31, dtype=np.uint8)
    np.fill_diagonal(masks, 1)
    for line in lines:
        for a, b in combinations(line, 2):
            masks[a, b] = masks[b, a] = 2

    for (a, b), (i, j) in close.items():
        square = 1 << (1 + (np.arange(3)[:, None] != i).astype(int)
                       + (np.arange(3)[None, :] != j))
        masks[np.ix_(named[a], named[b])] &= square.astype(np.uint8)
        masks[np.ix_(named[b], named[a])] &= square.T.astype(np.uint8)
    opposite = [(root, b) for b in partners]
    opposite += [(a, b) for a, b in combinations(partners, 2)
                 if sum(x != y for x, y in zip(a[1:], b[1:])) == 1]
    opposite += [(a, b) for a, b in combinations(middle, 2) if a[1] != b[1]]
    for a, b in opposite:
        masks[np.ix_(named[a], named[b])] &= 24
        masks[np.ix_(named[b], named[a])] &= 24
    return constraints.propagate(masks, lines)


def main():
    constraints = LineConstraints()
    survivors = [choice for choice in product(range(2), repeat=3)
                 if replay(choice, constraints)]
    assert survivors == [(0, 0, 0), (0, 0, 1), (0, 1, 0), (1, 0, 0)]
    print("Eight-face attachment choices", survivors)


if __name__ == "__main__":
    main()
