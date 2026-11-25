"""Check intermediate lines and face diagonals of the smallest rooted surfaces."""

from itertools import combinations, product

import numpy as np

from failed_pair_types import completions
from grid_path_closure import LineConstraints, extend_paths
from theta_profiles import LINES, MIDDLE, seed_distances


def suspension(size):
    faces = [(pole, edge // 2, ((edge + 1) // 2) % size)
             for pole, edge in product(range(2), range(2 * size))]
    parity = [(face[0] + pos) % 2 for pos, face in enumerate(faces)]
    return faces, parity


def joined_octahedra():
    faces = [face for side in (1, 2) for face in product((0, side), repeat=3) if any(face)]
    parity = [(sum(value != 0 for value in face) + (2 in face)) % 2 for face in faces]
    return faces, parity


def surface_seed(shape, groups):
    """Build the face partners and their paths from one root line.

    A color's group labels specify its intermediate lines. These small
    surfaces have distinct ancestor points within each color.
    """
    faces, parity = shape
    root = ("root",)
    middle = [("middle", i, j) for i, group in enumerate(groups)
              for j in range(len(group))]
    partners = [("partner", *face) for face in faces]
    point_ids = {}

    def index(label):
        return point_ids.setdefault(label, len(point_ids))

    named = {line: tuple(index((line, k)) for k in range(3))
             for line in (root, *partners)}
    connectors = []
    close = {}
    for i, group in enumerate(groups):
        parents = [index(("parent", i, j)) for j in range(len(group))]
        for j, parent in enumerate(parents):
            line = ("middle", i, j)
            named[line] = (parent, index((line, 1)), index((line, 2)))
            close[root, line] = i, 0
        for label in sorted(set(group)):
            points = [parent for j, parent in enumerate(parents) if group[j] == label]
            if len(points) == 1:
                points.append(index(("tip", i, label)))
            assert len(points) == 2
            connectors.append((named[root][i], *points))

    for pos, t in enumerate(faces):
        b = ("partner", *t)
        for i, s in enumerate(t):
            c = ("middle", i, s)
            # At one middle line its link alternates between the two
            # outward points. This parity alternates across every face edge.
            j = 1 + parity[pos]
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
    neighbors = {b: set() for b in partners}
    for a, b in combinations(partners, 2):
        if sum(x != y for x, y in zip(a[1:], b[1:])) == 1:
            opposite.append((a, b))
            neighbors[a].add(b)
            neighbors[b].add(a)
    parity_at = dict(zip(partners, parity))
    assert all(len(adjacent) == 3 and all(parity_at[b] != parity_at[v] for v in adjacent)
               for b, adjacent in neighbors.items())
    opposite += [(a, b) for a, b in combinations(middle, 2) if a[1] != b[1]]
    for a, b in opposite:
        masks[np.ix_(named[a], named[b])] &= 24
        masks[np.ix_(named[b], named[a])] &= 24
    diagonals = {tuple(sorted((a, b))) for adjacent in neighbors.values()
                 for a, b in combinations(adjacent, 2)}
    # Each face supplies all 36 slots for its opposite pair. The three
    # geodesics between its unused tips go through its adjacent partners.
    third = {}
    for line in lines:
        for i in range(3):
            edge = tuple(sorted((line[i-1], line[i-2])))
            assert edge not in third
            third[edge] = line[i]

    def spare(x, y):
        return third[tuple(sorted((x, y)))]

    templates = []
    for b in partners:
        slots = []
        for i, label in enumerate(b[1:]):
            parent, far, tip = named["middle", i, label]
            if parity_at[b]:
                far, tip = tip, far
            a_i, b_i = named[root][i], named[b][i]
            slots.extend([a_i, parent, far, b_i, spare(a_i, parent), tip, spare(far, b_i)])
        for j, k in combinations(range(3), 2):
            other, = [v for v in neighbors[b] if v[j+1] == b[j+1] and v[k+1] == b[k+1]]
            x, y = named[other][j], named[other][k]
            slots.extend([x, y, spare(slots[MIDDLE[j]], x), spare(x, y),
                          spare(y, slots[MIDDLE[k]])])
        templates.append(np.array(slots))
    return masks, lines, [(named[a], named[b]) for a, b in sorted(diagonals)], templates


def diagonal_domains(constraints, masks, lines, diagonals):
    """Try every ordered line square of distance two or three at each diagonal."""
    squares = (1 << constraints.squares).astype(np.uint8)
    squares = squares[np.isin(squares.min(axis=(1, 2)), (2, 4))]
    cells = [np.ix_(a, b) for a, b in diagonals]
    for cell in cells:
        masks[cell] &= 30  # the two lines are disjoint
    masks &= masks.T
    if not constraints.propagate(masks, lines):
        return None
    while True:
        previous = masks.copy()
        answer = []
        for cell in cells:
            possible = squares[((squares & masks[cell]) != 0).all(axis=(1, 2))]
            if len(possible) > 1:
                supported = []
                for square in possible:
                    trial = masks.copy()
                    trial[cell] &= square
                    trial &= trial.T
                    if constraints.propagate(trial, lines):
                        supported.append(square)
                possible = supported
            if not len(possible):
                return None
            masks[cell] &= np.bitwise_or.reduce(possible, axis=0)
            masks &= masks.T
            answer.append(frozenset(int(sq.min()).bit_length() for sq in possible))
        if np.array_equal(previous, masks):
            return answer
        if not constraints.propagate(masks, lines):
            return None


def parent_groups(size):
    distinct = tuple(range(size))
    return [distinct] + [tuple(a if k == b else k for k in distinct)
                         for a, b in combinations(distinct, 2)]


def failed_metrics(constraints):
    template, lines, definitions = extend_paths(
        (1 << seed_distances()).astype(np.uint8), LINES, list(combinations(MIDDLE, 2)))
    for i, (j, k) in enumerate(combinations(MIDDLE, 2)):
        x, y = 21 + 5*i, 22 + 5*i
        assert definitions[5*i:5*i+5] == [
            ("path", j, k, 1), ("path", j, k, 2),
            ("third", j, x, 1), ("third", x, y, 1), ("third", y, k, 1)]
    metrics = sorted((d for d in completions(constraints, template, lines)
                      if d[5, 34] != 1), key=lambda d: d.tobytes())
    assert len(metrics) == 32
    return np.array(metrics)


def first_face_counts(constraints, metrics, shape, groups):
    masks, lines, _, templates = surface_seed(shape, groups)
    cell = np.ix_(templates[0], templates[0])
    counts = [0, 0]
    for distances in metrics:
        trial = masks.copy()
        np.bitwise_and.at(trial, cell, (1 << distances).astype(np.uint8))
        if constraints.propagate(trial, lines):
            has_four = (distances[[5, 12, 19], [34, 29, 24]] == 4).any()
            counts[int(has_four)] += 1
    return counts


def check_assignment(constraints, metrics, shape, groups, indices):
    """Replay one simultaneous assignment to all faces of a surviving type."""
    masks, lines, diagonals, templates = surface_seed(shape, groups)
    assert len(indices) == len(templates)
    for slots, index in zip(templates, indices):
        np.bitwise_and.at(masks, np.ix_(slots, slots),
                          (1 << metrics[index]).astype(np.uint8))
    assert diagonal_domains(constraints, masks, lines, diagonals) is not None


def main():
    constraints = LineConstraints()
    metrics = failed_metrics(constraints)
    for size in (2, 3):
        shape = suspension(size)
        successes = {}
        choices = list(product(parent_groups(2), parent_groups(size), parent_groups(size)))
        for groups in choices:
            masks, lines, diagonals, _ = surface_seed(shape, groups)
            result = diagonal_domains(constraints, masks, lines, diagonals)
            if result is not None:
                successes[groups] = result
        distinct = tuple(range(size))
        separate = ((0, 1), distinct, distinct)
        if size == 2:
            assert set(successes) == {separate}
            assert successes[separate] == [{3}] * 12
            assert first_face_counts(constraints, metrics, shape, separate) == [8, 0]
        else:
            shared = ((0, 0), distinct, distinct)
            assert set(successes) == {separate, shared}
            assert successes[separate] == [{2, 3}] * 24
            assert successes[shared].count({3}) == successes[shared].count({2, 3}) == 12
            assert first_face_counts(constraints, metrics, shape, shared) == [0, 0]
            assert first_face_counts(constraints, metrics, shape, separate) == [8, 6]
        check_assignment(constraints, metrics, shape, separate, [0] * (4 * size))
        print(4 * size, "faces:", len(choices), "labelled attachment choices, only",
              separate, "passes both tests")

    shape = joined_octahedra()
    remaining = set()
    for groups in product(parent_groups(3), repeat=3):
        masks, lines, diagonals, _ = surface_seed(shape, groups)
        if diagonal_domains(constraints, masks, lines, diagonals) is not None:
            remaining.add(groups)
    assert remaining == set(product(((0, 1, 2), (0, 1, 1)), repeat=3))
    # Permuting the colors leaves only the number of shared pairs to choose.
    witnesses = [
        [0] * 14,
        [5, 0, 0, 0, 0, 0, 0] * 2,
        [5, 10, 1, 0, 12, 12, 0] * 2,
        [5, 10, 1, 28, 13, 12, 2] * 2,
    ]
    for count, indices in enumerate(witnesses):
        groups = ((0, 1, 2),) * (3-count) + ((0, 1, 1),) * count
        for face, index in zip(shape[0], indices):
            touches_shared_edge = sum(label != 0 for label in face) == 1 and any(face[3-count:])
            has_four = (metrics[index, [5, 12, 19], [34, 29, 24]] == 4).any()
            assert bool(has_four) == touches_shared_edge
        check_assignment(constraints, metrics, shape, groups, indices)
    print("14 faces: 64 labelled attachment choices, 8 pass the diagonal test")


if __name__ == "__main__":
    main()
