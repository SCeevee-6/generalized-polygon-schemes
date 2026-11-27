"""Verify the 24-residue count in all_residue_marks.md."""

from collections import Counter
from itertools import combinations

from indexed_eight_obstruction import coordinates, frame_graph, two_color_faces


# Move-0 joins added after the internal residue trees. Each entry is the
# frame (x, y, c, d, e); the other end of the join is (y, x, c, d, e).
JOINS = (
    (0, 6, 1, 9, 2),
    (1, 2, 5, 0, 1),
    (1, 5, 6, 0, 5),
    (1, 6, 1, 2, 10),
    (1, 10, 10, 5, 1),
    (3, 8, 4, 2, 12),
    (3, 16, 8, 12, 6),
    (4, 17, 13, 8, 7),
    (6, 8, 16, 9, 2),
    (6, 18, 15, 9, 16),
    (7, 9, 3, 9, 17),
    (7, 12, 17, 3, 11),
    (8, 13, 18, 4, 12),
    (10, 15, 11, 19, 5),
    (10, 20, 19, 11, 15),
    (11, 16, 20, 6, 12),
    (12, 17, 7, 21, 13),
    (12, 20, 11, 21, 17),
    (14, 16, 14, 20, 6),
    (18, 19, 23, 16, 15),
    (19, 23, 20, 23, 19),
    (21, 22, 23, 22, 18),
    (22, 23, 22, 23, 21),
)


def setup():
    vertices, _, incident, dual = coordinates()
    moves, roots = frame_graph(incident, dual)
    frames = [(x, y, c, d, e) for x in range(24) for y in range(24)
              if x != y and len(incident[x] & incident[y]) == 3
              for c, d, e in _orders(sorted(incident[x] & incident[y]))]
    assert [f[0] for f in frames] == roots
    faces = two_color_faces(moves)
    edges = sorted({tuple(sorted((i, j))) for i, row in enumerate(moves) for j in row})
    return vertices, incident, moves, roots, frames, faces, edges


def _orders(t):
    a, b, c = t
    return [(a, b, c), (a, c, b), (b, a, c), (b, c, a), (c, a, b), (c, b, a)]


def assignment_tree(moves, roots, frames, edges):
    index = {f: i for i, f in enumerate(frames)}
    eid = {e: k for k, e in enumerate(edges)}
    internal = {x: sorted(eid[tuple(sorted((i, j)))] for i in range(len(moves))
                          if roots[i] == x for j in moves[i][1:])
                for x in range(24)}
    joins = []
    for x, y, c, d, e in JOINS:
        i, j = index[x, y, c, d, e], index[y, x, c, d, e]
        assert moves[i][0] == j
        joins.append(eid[tuple(sorted((i, j)))])
    parent = list(range(len(moves)))

    def rep(i):
        while i != parent[i]:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    tree = []
    for k in sorted(e for x in range(24) for e in internal[x]) + joins:
        a, b = edges[k]
        ra, rb = rep(a), rep(b)
        if ra != rb:
            parent[ra] = rb
            tree.append(k)
    assert len(tree) == len(moves) - 1
    assert set(joins) <= set(tree)
    return tree, internal


def filled_edges(tree, internal, small, unconditional, face_of, n):
    known = bytearray(n)
    for k in tree:
        known[k] = 1
    for x in small:
        for k in internal[x]:
            known[k] = 1
    missing = [sum(1 for e in f if not known[e]) for f in unconditional]
    stack = [k for k, m in enumerate(missing) if m == 1]
    while stack:
        k = stack.pop()
        if missing[k] != 1:
            continue
        e = next(e for e in unconditional[k] if not known[e])
        known[e] = 1
        for k2 in face_of[e]:
            missing[k2] -= 1
            if missing[k2] == 1:
                stack.append(k2)
    return known


def fills(tree, internal, small, unconditional, face_of, n):
    return all(filled_edges(tree, internal, small, unconditional, face_of, n))


def reduce(vectors, basis):
    basis = dict(basis)
    for v in vectors:
        while v:
            h = v.bit_length() - 1
            if h not in basis:
                basis[h] = v
                break
            v ^= basis[h]
    return basis


def main():
    vertices, incident, moves, roots, frames, faces, edges = setup()
    eid = {e: k for k, e in enumerate(edges)}
    ip = lambda a, b: sum(p*q for p, q in zip(vertices[a], vertices[b]))
    unconditional = [[eid[e] for e in f] for cond, _, f in faces if not cond]
    face_of = [[] for _ in edges]
    for k, f in enumerate(unconditional):
        for e in f:
            face_of[e].append(k)
    tree, internal = assignment_tree(moves, roots, frames, edges)

    # Every marked set must meet the complement of each filling omitted set.
    everyone = set(range(24))
    singles = [x for x in range(24)
               if fills(tree, internal, everyone - {x}, unconditional, face_of, len(edges))]
    assert singles == list(range(24))
    pairs = Counter()
    for x, y in combinations(range(24), 2):
        ok = fills(tree, internal, everyone - {x, y}, unconditional, face_of, len(edges))
        pairs[ip(x, y), ok] += 1
    assert all(ok for (prod, ok) in pairs if prod != 1)
    print("one assignment tree; all 24 single omissions fill; pairs", dict(pairs))

    pairs_total = 2925*2048//2
    incidences = -(-6*pairs_total // 36)
    incidences += -incidences % 4
    indices = -(-incidences // 128)
    assert (incidences, indices) == (499200, 3900)
    print("24a >= 2 and a <= 3t give", incidences, "long-cycle incidences in at least",
          indices, "indices")

    # Limits of the twenty-residue test, over GF(2).
    vec = lambda face: sum(1 << eid[e] for e in face)
    base = reduce([vec(f) for cond, _, f in faces if not cond], {})
    octagons = [(roots[f0], vec(f)) for cond, f0, f in faces if cond]
    full = len(reduce([v for _, v in octagons], base))
    assert (len(base), full) == (1056, 1153) == (1056, len(edges) - len(moves) + 1)
    rank = lambda omit: len(reduce([v for r, v in octagons if r not in omit], base))
    good4 = [om for om in combinations(range(24), 4) if rank(set(om)) == full]
    good4set = set(good4)
    good5 = [om for om in combinations(range(24), 5)
             if all(s in good4set for s in combinations(om, 4))]
    assert len(good4) == 72 and not good5
    patterns = Counter(tuple(sorted(ip(a, b) for a, b in combinations(om, 2))) for om in good4)
    assert set(patterns) == {(-1, -1, -1, -1, 0, 0)}
    # The two half-cubes of each vertex figure are the other two 16-cells.
    cell = {}
    for x in range(24):
        cell.setdefault(min(y for y in range(24) if ip(x, y) in (0, 2, -2)), []).append(x)
    cell_of = {x: k for k, xs in enumerate(cell.values()) for x in xs}
    neighbors = lambda x: [y for y in range(24) if ip(x, y) == 1]
    for x in range(24):
        for u, v in combinations(neighbors(x), 2):
            diagonal = len(incident[x] & incident[u] & incident[v]) == 1
            assert diagonal == (cell_of[u] == cell_of[v])
    print("GF(2): 72 omitted 4-sets keep full rank, no omitted 5-set does;",
          "half-cubes are 16-cells")


if __name__ == "__main__":
    main()
