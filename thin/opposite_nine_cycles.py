"""Replay the forced distance data of every nine-cycle trace."""

from collections import Counter
from itertools import combinations, product

import numpy as np

from grid_path_closure import LineConstraints, extend_paths
from opposite_cycles import line_pair, traces


def attachment_choices(vertex, distances):
    """Enumerate local projection choices modulo permutations of the line."""
    neighbors = sorted(b if a == vertex else a for (a, b), d in distances.items()
                       if d == 2 and vertex in (a, b))

    def assign(index, chosen):
        if index == len(neighbors):
            yield chosen
            return
        v = neighbors[index]
        forbidden = {color for u, color in chosen.items()
                     if distances.get(line_pair(u, v)) == 4}
        for color in range(min(3, 1 + len(set(chosen.values())))):
            if color not in forbidden:
                yield from assign(index + 1, {**chosen, v: color})

    yield from assign(0, {})


def named_lines(distances, attachments):
    count = len(attachments)
    lines = [tuple(range(3*i, 3*i + 3)) for i in range(count)]
    size = 3*count + sum(d == 2 for d in distances.values())
    masks = np.full((size, size), 31, dtype=np.uint8)
    np.fill_diagonal(masks, 1)
    slot = 3*count
    for (a, b), d in distances.items():
        if d == 4:
            square = np.full((3, 3), 24, dtype=np.uint8)
        else:
            x, y = attachments[a][b], attachments[b][a]
            square = 1 << (1 + (np.arange(3)[:, None] != x).astype(int)
                           + (np.arange(3)[None, :] != y))
            lines.append((3*a + x, 3*b + y, slot))
            slot += 1
        masks[np.ix_(lines[a], lines[b])] = square
        masks[np.ix_(lines[b], lines[a])] = square.T
    return masks, lines


def unused_tips(path, attachments):
    data = []
    for i, (a, b) in enumerate(path):
        neighbors = (path[i - 1], path[(i + 1) % len(path)])
        tips = {}
        for c in set(neighbors[0]) | set(neighbors[1]):
            unused, = {0, 1, 2} - {attachments[c][a], attachments[c][b]}
            tips[c] = 3*c + unused
        data.append((tips, neighbors))
    return data


def replay(path, attachments, masks, lines, constraints):
    assert constraints.propagate(masks, lines)
    assert np.all((masks & (masks - 1)) == 0)
    representatives = [np.flatnonzero(row == 1)[0] for row in masks]
    distinct_lines = {tuple(sorted(representatives[p] for p in line)) for line in lines}
    degrees = Counter(point for line in distinct_lines for point in line)
    assert len(set(representatives)) == 36 and len(distinct_lines) == 21
    assert Counter(degrees.values()) == {1: 9, 2: 27}
    dangling_points = [point for point, degree in degrees.items() if degree == 1]
    assert all(masks[a, b] == 8 for a, b in combinations(dangling_points, 2))

    tips_by_vertex = unused_tips(path, attachments)
    requests = sorted({line_pair(a, b) for tips, _ in tips_by_vertex
                       for a, b in combinations(tips.values(), 2)})
    enlarged, extended_lines, definitions = extend_paths(masks, lines, requests)
    assert constraints.propagate(enlarged, extended_lines)
    # Extend_paths writes the connecting line on the middle edge as a
    # named line, either newly inserted or present already.
    new_points = {definition: len(masks) + i for i, definition in enumerate(definitions)}

    def internal_point(a, b, k):
        known = np.flatnonzero((masks[a] == (1 << k)) & (masks[b] == (1 << (3-k))))
        return int(known[0]) if len(known) else new_points["path", a, b, k]

    for tips, neighbors in tips_by_vertex:
        for middle_pair in combinations(sorted(tips), 2):
            a, b = line_pair(*(tips[i] for i in middle_pair))
            x, y = [internal_point(a, b, k) for k in (1, 2)]
            opposite_tip, = set(tips) - set(middle_pair)
            third_points = [next(iter(set(line) - {x, y})) for line in extended_lines
                            if x in line and y in line]
            expected = 4 if middle_pair in neighbors else 3
            assert third_points and all(enlarged[tips[opposite_tip], z] == 1 << expected
                                        for z in third_points)


def main():
    constraints = LineConstraints()
    count = 0
    for path, distances in traces(9):
        labels = 1 + max(max(pair) for pair in path)
        assert labels == 9
        for attachments in product(*(attachment_choices(i, distances) for i in range(labels))):
            masks, lines = named_lines(distances, attachments)
            replay(path, attachments, masks, lines, constraints)
            count += 1
    assert count == 2
    print("Nine-cycle traces checked", count)
    print("Each forces 36 points, 21 lines and nine cycle edges of label 4")


if __name__ == "__main__":
    main()
