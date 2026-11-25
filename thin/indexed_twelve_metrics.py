"""Replay joint failed-pair metric assignments on the indexed twelve-cycle."""

from collections import Counter
from itertools import product

import numpy as np

from grid_path_closure import LineConstraints
from indexed_cycle_geometry import replay
from opposite_link_surfaces import failed_metrics


# A comparison word and the metric index at each of its twelve vertices.
WITNESSES = {
    "000000000000": [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    "000000000001": [10, 4, 21, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    "000000000011": [6, 5, 21, 0, 0, 0, 0, 0, 0, 0, 0, 10],
    "000000000111": [7, 5, 29, 0, 0, 0, 0, 0, 0, 0, 26, 6],
    "000000001001": [10, 4, 21, 0, 0, 0, 0, 0, 0, 10, 4, 21],
    "000000001111": [7, 5, 29, 0, 0, 0, 0, 0, 0, 26, 6, 7],
    "000000100001": [10, 4, 21, 0, 0, 0, 0, 10, 4, 21, 0, 0],
    "000001000001": [10, 4, 21, 0, 0, 0, 10, 4, 21, 0, 0, 0],
    "000001000011": [6, 5, 21, 0, 0, 0, 10, 4, 21, 0, 0, 10],
    "000001001001": [10, 4, 21, 0, 0, 0, 10, 4, 21, 10, 4, 21],
    "000010000111": [7, 5, 29, 0, 0, 10, 4, 21, 0, 0, 26, 6],
    "000011000011": [6, 5, 21, 0, 0, 10, 6, 5, 21, 0, 0, 10],
    "001001001001": [10, 4, 21, 10, 4, 21, 10, 4, 21, 10, 4, 21],
}


def words():
    for choices in product(range(4), repeat=4):
        word = [0]*12
        for color, pos in enumerate(choices):
            if pos < 3:
                word[4*pos+color] = 1
        s = "".join(map(str, word))
        if "101" not in s+s[:2]:
            yield min(v[i:] + v[:i] for v in (s, s[::-1]) for i in range(12))


def placements(word):
    leaves = []
    for t in range(12):
        color, pos = t % 4, t // 4
        share = next((k for k in range(3) if word[4*k+color] == "1"), None)
        parent, ancestor = pos, 0
        if share is not None and pos == (share+1) % 3:
            parent, ancestor = share, 1
        leaves.append((color % 2, color // 2, parent, ancestor, 0))
    edges = []
    for i, a in enumerate(leaves):
        b = leaves[(i+1) % 12]
        edges.append((b, a) if i % 2 else (a, b))
    return sorted({p[:3] for p in leaves}), edges, leaves


def check_case(word, choices, metrics, constraints):
    groups, edges, leaves = placements(word)
    result, data = replay(groups, edges, constraints)
    assert result == "open"
    masks, _, lines, names = data
    named = dict(zip(names, lines))
    third = {(a, b): c for line in lines for a in line for b in line if a != b
             for c in line if c != a and c != b}

    def attachments(a, b):
        x, y = np.nonzero(masks[np.ix_(a, b)] == 2)
        assert len(x) == len(y) == 1
        return a[x[0]], b[y[0]]

    def leaf(i):
        return named[("leaf", *leaves[i % 12])]

    n = len(masks)
    expanded = np.full((n+60, n+60), 31, dtype=np.uint8)
    expanded[:n, :n] = masks
    np.fill_diagonal(expanded, 1)
    lines = list(lines)
    for i, choice in enumerate(choices):
        side = i % 2
        a, b = named["root", 1-side], leaf(i)
        middles = [named["root", side], leaf(i-1), leaf(i+1)]
        slots = []
        for line in middles:
            p, x = attachments(a, line)
            y, q = attachments(line, b)
            slots.extend((p, x, y, q, third[p, x], third[x, y], third[y, q]))
        # For the first two tip pairs the third arm goes through a cycle leaf.
        for delta, k in [(-1, 1), (1, 2)]:
            far = leaf(i+2*delta)
            x, u = attachments(middles[0], far)
            v, y = attachments(far, middles[k])
            assert (x, y) == (slots[5], slots[7*k+5])
            slots.extend((u, v, third[x, u], third[u, v], third[v, y]))
        # Leave the five extra slots free to coincide with other named points.
        x, y, left, mid, right = range(n+5*i, n+5*i+5)
        lines.extend(((slots[12], x, left), (x, y, mid), (y, slots[19], right)))
        slots.extend((x, y, left, mid, right))
        assert len(slots) == 36
        np.bitwise_and.at(expanded, np.ix_(slots, slots), metrics[choice])
    assert constraints.propagate(expanded, lines)


def main():
    classes = Counter(words())
    assert set(classes) == set(WITNESSES) and sum(classes.values()) == 136
    constraints = LineConstraints()
    metrics = (1 << failed_metrics(constraints)).astype(np.uint8)
    for word, choices in sorted(WITNESSES.items()):
        assert len(choices) == 12
        check_case(word, choices, metrics, constraints)
    print(len(classes), "dihedral classes cover 136 labelled twelve-cycle patterns")
    print("All joint metric witnesses pass three-line propagation")


if __name__ == "__main__":
    main()
