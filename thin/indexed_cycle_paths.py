"""Replay attachment restrictions on open paths in an indexed cycle."""

from itertools import product

from grid_path_closure import LineConstraints
from indexed_cycle_geometry import replay


def window(relations):
    """Name a path with four periodic side and branch classes."""
    classes = []
    for side, branch in ((0, 0), (1, 0), (0, 1), (1, 1)):
        position = len(classes)
        leaves = [(side, branch, 0, 0, 0)]
        if position < len(relations):
            relation = relations[position]
            leaves.append((side, branch, int(relation == 0),
                           int(relation == 1), int(relation == 2)))
        classes.append(leaves)
    path = [classes[i % 4][i // 4] for i in range(4 + len(relations))]
    edges = sorted((a, b) if a[0] == 0 else (b, a)
                   for a, b in zip(path, path[1:]))
    return sorted({leaf[:3] for leaf in path}), edges


def main():
    constraints = LineConstraints()
    for size in (1, 3):
        survivors = set()
        for relations in product(range(3), repeat=size):
            result, _ = replay(*window(relations), constraints)
            if result == "open":
                survivors.add(relations)
        expected = set(product(range(2), repeat=size)) - {(1, 0, 1)}
        assert survivors == expected
        print("Path length", size + 3, "surviving comparison words", sorted(survivors))

    # For each of four colors on a twelve-cycle, choose the location
    # of the one possible shared parent or choose no shared parent.
    total = 0
    for choices in product((None, 0, 1, 2), repeat=4):
        word = [int(choices[color] == visit)
                for visit in range(3) for color in range(4)]
        read = word + word[:2]
        if all(read[i:i+3] != [1, 0, 1] for i in range(12)):
            total += 1
    assert total == 136
    print("Indexed twelve-cycle attachment choices", total)


if __name__ == "__main__":
    main()
