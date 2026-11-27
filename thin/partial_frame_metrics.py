"""Necessary point distances in a partially consistent 24-cell frame assignment."""
from collections import defaultdict
from itertools import combinations

import numpy as np

from all_residue_marks import assignment_tree, filled_edges, setup
from indexed_eight_obstruction import coordinates, ports


def pair(a, b):
    return tuple(sorted((a, b)))


class Reference:
    def __init__(self):
        _, self.incident, self.moves, self.roots, self.frames, faces, self.edges = setup()
        self.tree, self.internal = assignment_tree(self.moves, self.roots, self.frames, self.edges)
        self.number = {edge: i for i, edge in enumerate(self.edges)}
        self.unconditional = [[self.number[e] for e in face] for conditional, _, face in faces
                              if not conditional]
        self.face_of = [[] for _ in self.edges]
        for i, face in enumerate(self.unconditional):
            for e in face:
                self.face_of[e].append(i)
        vertices, normals, _, dual = coordinates()
        self.port = ports(normals, vertices, self.incident, dual)

    def known(self, omitted):
        return filled_edges(self.tree, self.internal, set(range(24)) - set(omitted),
                            self.unconditional, self.face_of, len(self.edges))


class Classes:
    def __init__(self, size):
        self.parent = list(range(size))

    def rep(self, i):
        while i != self.parent[i]:
            self.parent[i] = self.parent[self.parent[i]]
            i = self.parent[i]
        return i

    def merge(self, i, j):
        i, j = self.rep(i), self.rep(j)
        if i == j:
            return False
        self.parent[max(i, j)] = min(i, j)
        return True


def line_relations(ref, omitted):
    known = ref.known(omitted)
    classes = Classes(5 * len(ref.frames))
    # On move 1 the new partner and the last middle line can change.
    preserved = [(1, 0, 2, 3, 4), (0, None, 2, 3, None), (0, 1, 2, 4, 3), (0, 1, 3, 2, 4)]
    for i, row in enumerate(ref.moves):
        for color, j in enumerate(row):
            if i < j and known[ref.number[i, j]]:
                for k, target in enumerate(preserved[color]):
                    if target is not None:
                        classes.merge(5*i+k, 5*j+target)

    while True:
        labels = [tuple(classes.rep(5*i+k) for k in range(5)) for i in range(len(ref.frames))]
        middles = defaultdict(set)
        for i, (a, b, c, d, e) in enumerate(labels):
            middles[pair(a, b)].add(frozenset((c, d, e)))
            j = ref.moves[i][1]
            if known[ref.number[pair(i, j)]]:
                middles[pair(c, d)].add(frozenset((a, b, labels[j][1])))
        changed = False
        for sets in middles.values():
            for first, second in combinations(sets, 2):
                if len(first & second) == 2:
                    x, = first - second
                    y, = second - first
                    changed |= classes.merge(x, y)
        if not changed:
            break

    near = defaultdict(set)
    opposite = set(middles)
    for ends, triples in middles.items():
        for triple in triples:
            assert len(triple) == 3
            opposite.update(pair(x, y) for x, y in combinations(triple, 2))
            for x in ends:
                near[x].update(triple)
                for y in triple:
                    near[y].add(x)

    # Use the two proven diagonal distances only at guaranteed small roots.
    distance_three = set()
    for x in sorted(set(ref.roots) - set(omitted)):
        partners, facets = defaultdict(set), defaultdict(set)
        for i, (r, y, c, d, e) in enumerate(ref.frames):
            if r == x:
                partners[y].add(labels[i][1])
                for k, f in enumerate((c, d, e), 2):
                    facets[f].add(labels[i][k])
        assert all(len(v) == 1 for v in (*partners.values(), *facets.values()))
        for y, z in combinations(partners, 2):
            if len(ref.incident[x] & ref.incident[y] & ref.incident[z]) == 1:
                distance_three.add(pair(min(partners[y]), min(partners[z])))
        for c, d in combinations(facets, 2):
            if ref.port[x, c] == ref.port[x, d]:
                distance_three.add(pair(min(facets[c]), min(facets[d])))
    return near, opposite, distance_three


def attachment_colorings(neighbors, opposite):
    """Enumerate colorings up to renaming the three points of the central line."""
    edges = {v: {w for w in neighbors if pair(v, w) in opposite} for v in neighbors}
    order = sorted(neighbors, key=lambda v: -len(edges[v]))
    assigned, rows = {}, []

    def visit(pos):
        if pos == len(order):
            rows.append([assigned[v] for v in neighbors])
            return
        v = order[pos]
        used = {assigned[w] for w in edges[v] if w in assigned}
        for color in range(min(3, 2 + max(assigned.values(), default=-1))):
            if color not in used:
                assigned[v] = color
                visit(pos+1)
                del assigned[v]

    visit(0)
    assert rows
    return np.asarray(rows, dtype=np.uint8)


def metric_seed(ref, omitted):
    near, opposite, distance_three = line_relations(ref, omitted)
    count = 0
    central, attachment, extra = {}, {}, []
    for a in sorted(near):
        neighbors = sorted(near[a])
        colors = attachment_colorings(neighbors, opposite)
        # Ports used in every coloring by the same point get a common slot.
        bank = {}
        for j, b in enumerate(neighbors):
            key = colors[:, j].tobytes()
            if key not in bank:
                bank[key] = count
                count += 1
            attachment[a, b] = bank[key]
        # Fix up to three pairwise distinct ports and fill any unused slots.
        basis = next(group for size in (3, 2, 1) for group in combinations(neighbors, size)
                     if all(pair(v, w) in opposite for v, w in combinations(group, 2)))
        line = [attachment[a, b] for b in basis]
        while len(line) < 3:
            line.append(count)
            count += 1
        central[a] = tuple(line)
        extra.extend((p, line) for p in bank.values() if p not in line)

    lines = list(central.values())
    for a, b in sorted(attachment):
        if a < b:
            lines.append((attachment[a, b], attachment[b, a], count))
            count += 1
    masks = np.full((count, count), 31, dtype=np.uint8)
    np.fill_diagonal(masks, 1)
    for line in lines:
        for p, q in combinations(line, 2):
            masks[p, q] = masks[q, p] = 2
    for port, line in extra:
        for p in line:
            masks[port, p] &= 3
            masks[p, port] &= 3
    for relations, allowed in ((opposite, 24), (distance_three, 28)):
        for a, b in relations:
            masks[np.ix_(central[a], central[b])] &= allowed
            masks[np.ix_(central[b], central[a])] &= allowed
    return masks, lines
