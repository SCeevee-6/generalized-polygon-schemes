"""Replay the line-label and attachment checks for lengths 11 through 13."""

from collections import Counter
from itertools import combinations, product

import numpy as np

from grid_path_closure import LineConstraints, extend_paths
from opposite_cycles import line_pair, traces
from opposite_nine_cycles import attachment_choices, named_lines, unused_tips


def trace_key(trace):
    def rename(path):
        labels = {}
        for pair in path:
            for line in pair:
                labels.setdefault(line, len(labels))
        return tuple(line_pair(*(labels[line] for line in pair)) for pair in path)

    variants = []
    for direction in (trace, trace[::-1]):
        for start in range(len(trace)):
            path = direction[start:] + direction[:start]
            # The first two vertices are disjoint. After those, each
            # vertex introduces at most one new line label.
            for first, second in product((False, True), repeat=2):
                prefix = (path[0][::-1] if first else path[0],
                          path[1][::-1] if second else path[1])
                variants.append(rename(prefix + path[2:]))
    return min(variants)


def classified_traces(length):
    representatives = {}
    for path, distances in traces(length):
        representatives.setdefault(trace_key(path), (path, distances))
    return representatives.values()


def edge_labels(path, attachments, masks, lines, constraints):
    tips_at_vertex = unused_tips(path, attachments)
    requests = sorted({line_pair(a, b) for tips, _ in tips_at_vertex
                       for a, b in combinations(tips.values(), 2)})
    # The three unused points of the middle lines of each opposite pair
    # have mutual point distance three.
    for a, b in requests:
        assert masks[a, b] & 8
        masks[a, b] = masks[b, a] = 8
    expanded, extended_lines, definitions = extend_paths(masks, lines, requests)
    assert constraints.propagate(expanded, extended_lines)
    new_slots = {meaning: len(masks) + i for i, meaning in enumerate(definitions)}

    def internal_point(a, b, k):
        existing = np.flatnonzero((masks[a] == 1 << k) & (masks[b] == 1 << (3-k)))
        return int(existing[0]) if len(existing) else new_slots["path", a, b, k]

    labels = []
    for tips, neighbors in tips_at_vertex:
        at_vertex = {}
        for pair in combinations(sorted(tips), 2):
            a, b = line_pair(*(tips[i] for i in pair))
            x, y = [internal_point(a, b, k) for k in (1, 2)]
            omitted, = set(tips) - set(pair)
            third = [z for line in extended_lines if x in line and y in line
                     for z in line if z not in (x, y)]
            values = {int(expanded[tips[omitted], z]) for z in third}
            value, = values
            assert value in (8, 16)
            at_vertex[pair] = value.bit_length() - 1
        assert sorted(at_vertex.values()) == [3, 4, 4]
        labels.append(tuple(at_vertex[pair] for pair in neighbors))
    assert all(labels[i][1] == labels[(i + 1) % len(path)][0] for i in range(len(path)))
    return tuple(after for before, after in labels)


def main():
    constraints = LineConstraints()
    expected = {11: (8, 0, 1), 12: (6, 1, 5), 13: (71, 0, 1)}
    for length, result in expected.items():
        counts = Counter()
        for path, distances in classified_traces(length):
            labels = 1 + max(max(pair) for pair in path)
            choices = [list(attachment_choices(i, distances)) for i in range(labels)]
            if not all(choices):
                continue
            counts["traces"] += 1
            assert all(len(options) == 1 for options in choices)
            attachments = [options[0] for options in choices]
            masks, lines = named_lines(distances, attachments)
            if not constraints.propagate(masks, lines):
                continue
            indexed = length % 2 == 0 and all(
                set.intersection(*(set(pair) for pair in path[side::2])) for side in range(2)
            )
            counts["indexed" if indexed else "other"] += 1
            if length % 2:
                labels = edge_labels(path, attachments, masks, lines, constraints)
                assert labels.count(3) == {11: 2, 13: 4}[length]
                print("Odd-cycle edge labels", length, labels, flush=True)
        counts = tuple(counts[k] for k in ("traces", "indexed", "other"))
        assert counts == result
        print(length, *counts, flush=True)


if __name__ == "__main__":
    main()
