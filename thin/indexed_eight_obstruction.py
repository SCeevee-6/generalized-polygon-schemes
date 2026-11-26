"""Verify the 24-cell obstruction in indexed_eight_obstruction.md."""

from collections import deque
from itertools import combinations, permutations, product


def coordinates():
    vertices = [v for v in product((-1, 0, 1), repeat=4) if sum(x*x for x in v) == 2]
    axial = [tuple(s if k == i else 0 for k in range(4))
             for s in (-2, 2) for i in range(4)]
    normals = sorted(axial + list(product((-1, 1), repeat=4)))
    incident = [{j for j, w in enumerate(normals) if sum(a*b for a, b in zip(v, w)) == 2}
                for v in vertices]
    dual = [{i for i, row in enumerate(incident) if j in row} for j in range(24)]
    return vertices, normals, incident, dual


def distances(graph, start):
    depth = {start: 0}
    queue = deque([start])
    while queue:
        i = queue.popleft()
        for j in graph[i]:
            if j not in depth:
                depth[j] = depth[i] + 1
                queue.append(j)
    return depth


def frame_graph(incident, dual):
    frames = [(x, y, c, d, e) for x, y in permutations(range(24), 2)
              if len(incident[x] & incident[y]) == 3
              for c, d, e in permutations(sorted(incident[x] & incident[y]))]
    index = {f: i for i, f in enumerate(frames)}
    assert len(frames) == len(index) == 1152
    moves = []
    for x, y, c, d, e in frames:
        z, = (dual[c] & dual[d]) - {x, y}
        f, = (incident[x] & incident[z]) - {c, d}
        moves.append(tuple(index[row] for row in (
            (y, x, c, d, e), (x, z, c, d, f), (x, y, c, e, d), (x, y, d, c, e))))
    assert all(moves[j][k] == i for i, row in enumerate(moves) for k, j in enumerate(row))
    for slot, excluded in [(0, 0), (2, 3)]:
        graph = [tuple(j for k, j in enumerate(row) if k != excluded) for row in moves]
        for label in range(24):
            residue = {i for i, f in enumerate(frames) if f[slot] == label}
            assert set(distances(graph, min(residue))) == residue
    return moves, [frame[0] for frame in frames]


def two_color_faces(moves):
    faces = []
    for a, b in combinations(range(4), 2):
        order = {(0, 1): 3, (1, 2): 4, (2, 3): 3}.get((a, b), 2)
        seen = set()
        for f in range(len(moves)):
            if f in seen:
                continue
            i, face = f, set()
            for step in range(2*order):
                assert i not in seen
                seen.add(i)
                j = moves[i][(a, b)[step % 2]]
                face.add(tuple(sorted((i, j))))
                i = j
            assert i == f and len(face) == 2*order
            faces.append((order == 4, f, face))

    return faces


def braid_checks(moves, faces):
    """Fill flag edges from a fixed tree, taking unconditional relations first."""
    depth = distances(moves, 0)
    assert len(depth) == len(moves)
    known = set()
    for i in range(1, len(moves)):
        j = next(j for j in moves[i] if depth[j] < depth[i])
        known.add(tuple(sorted((i, j))))
    pending = sorted(faces, key=lambda item: item[:2])
    tests = []
    all_edges = {tuple(sorted((i, j))) for i, row in enumerate(moves) for j in row}
    while known != all_edges:
        for k, (conditional, f, face) in enumerate(pending):
            missing = face - known
            if len(missing) == 1:
                known.update(missing)
                if conditional:
                    tests.append(f)
                del pending[k]
                break
        else:
            raise AssertionError("some flag edges have no filling relation")
    assert len(tests) == len(set(tests)) == 97
    return len(tests)


# Joins between tested residues. Their order only chooses a spanning tree.
ROOT_JOINS = (
    (1, 2), (1, 14), (2, 4), (3, 8), (3, 16), (6, 7), (6, 8),
    (8, 9), (10, 14), (10, 20), (11, 14), (12, 17), (12, 20),
    (13, 17), (13, 21), (14, 15), (14, 16), (17, 23), (18, 21),
)


def rooted_piece_checks(moves, roots, faces, incident):
    """Certify a spanning tree and fill edges using 20 rooted octahedra."""
    omitted = {0, 5, 19, 22}
    chosen = set(roots) - omitted
    known = {tuple(sorted((i, j))) for i, row in enumerate(moves) if roots[i] in chosen
             for j in row[1:]}
    edges = {tuple(sorted((i, j))) for i, row in enumerate(moves) for j in row}
    parent = list(range(len(moves)))

    def representative(i):
        while i != parent[i]:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    # Construct the assignment tree inside the chosen residues first.
    tree = []
    joins_by_root = {}
    for i, row in enumerate(moves):
        j = row[0]
        if i < j:
            key = tuple(sorted((roots[i], roots[j])))
            joins_by_root.setdefault(key, (i, j))
    prescribed = [joins_by_root[e] for e in ROOT_JOINS]
    for a, b in sorted(known) + prescribed + sorted(edges - known):
        x, y = representative(a), representative(b)
        if x != y:
            parent[x] = y
            tree.append((a, b))
    assert len(tree) == len(moves) - 1 and len(chosen) == 20
    # After contracting the twenty residue trees, these joins form a tree
    # on the tested roots. Each join exchanges the endpoints of one pair.
    joins = [(i, j) for i, j in tree if roots[i] != roots[j]
             and roots[i] in chosen and roots[j] in chosen]
    assert len(joins) == 19 and all(moves[i][0] == j for i, j in joins)
    assert set(joins) == set(prescribed)
    known.update(tree)
    unconditional = [face for conditional, _, face in faces if not conditional]
    while known != edges:
        size = len(known)
        for face in unconditional:
            missing = face - known
            if len(missing) == 1:
                known.update(missing)
        assert len(known) > size
    # At a small root the partners form a cube. Two partner faces at
    # cube distance two retain exactly one of their three middle labels.
    diagonals = 0
    for x in sorted(chosen):
        neighbors = {b if a == x else a for a, b in ROOT_JOINS if x in (a, b)}
        graph = {u: {v for v in neighbors if len(incident[x] & incident[u] & incident[v]) == 1}
                 for u in neighbors}
        # Take a spanning tree of each component of the diagonal graph.
        while neighbors:
            component = set(distances(graph, min(neighbors)))
            diagonals += len(component) - 1
            neighbors -= component
    assert diagonals == 17
    return len(chosen), len(joins), diagonals


def check_square_links():
    for high in product((0, 1), repeat=4):
        if sum(high) < 2:
            continue  # excluded in a large piece by the all-low triangle lemma
        unpaired = sum(high[i] and not high[(i+2) % 4] for i in range(4))
        face_counts = [high[i] + high[(i+1) % 4] for i in range(4)]
        assert unpaired <= 3*face_counts.count(2) - face_counts.count(0)


def ports(points, centers, incident, dual):
    """Number the three attachment points of each central line."""
    port = {}
    for i, neighbors in enumerate(incident):
        classes = {}
        for j in sorted(neighbors):
            paired = tuple(2*a-b for a, b in zip(centers[i], points[j]))
            key = min(points[j], paired)
            port[i, j] = 3*i + classes.setdefault(key, len(classes))
        values = [port[i, j] for j in neighbors]
        assert sorted(values) == [3*i, 3*i, 3*i+1, 3*i+1, 3*i+2, 3*i+2]
        for j, k in combinations(neighbors, 2):
            assert (len(dual[j] & dual[k]) == 3) == (port[i, j] != port[i, k])
        triples = {frozenset(neighbors & row) for h, row in enumerate(incident)
                   if h != i and len(neighbors & row) == 3}
        assert len(triples) == 8
        assert all(len({port[i, j] for j in triple}) == 3 for triple in triples)
    return port


def check_point_graph(vertices, normals, incident, dual):
    left = ports(normals, vertices, incident, dual)
    doubled = [tuple(2*a for a in v) for v in vertices]
    right = ports(doubled, normals, dual, incident)
    graph = [set() for _ in range(144)]
    for i in range(48):
        for a, b in combinations(range(3*i, 3*i+3), 2):
            graph[a].add(b)
            graph[b].add(a)
    for i, row in enumerate(incident):
        for j in row:
            a, b = left[i, j], 72 + right[j, i]
            graph[a].add(b)
            graph[b].add(a)
    assert all(len(row) == 4 for row in graph)
    groups = [(0, 15, 69, 54), (104, 107, 113, 110)]
    dist = {i: distances(graph, i) for group in groups for i in group}
    for group, targets in (groups, groups[::-1]):
        for i, p in enumerate(group):
            nxt = group[(i+1) % 4]
            assert dist[p][nxt] == 2 and dist[p][group[(i+2) % 4]] == 4
            midpoint, = graph[p] & graph[nxt]
            for q in targets:
                assert dist[q][midpoint] == 4
                assert all(dist[q][x] == 4 for x in graph[p])


def main():
    v, w, incident, dual = coordinates()
    moves, roots = frame_graph(incident, dual)
    faces = two_color_faces(moves)
    tests = braid_checks(moves, faces)
    rooted, joins, diagonals = rooted_piece_checks(moves, roots, faces, incident)
    check_square_links()
    check_point_graph(v, w, incident, dual)
    pairs = 2925*2048//2
    systems = (6*pairs + 128*tests-1)//(128*tests)
    assert systems == 1448
    print("1152 reference frames,", tests, "eight-cycle tests")
    print("Two four-point cycles verified on 144 named points")
    denominator = 3*rooted
    stronger = (6*pairs + 128*denominator-1)//(128*denominator)
    assert denominator == 60 and stronger == 2340
    print(rooted, "rooted-piece tests give at least", stronger, "systems with a long cycle")
    # (rooted - diagonals)*a + 2*diagonals*t >= 1 + (joins - diagonals)*b,
    # where a <= 3t. Count in units of four edge incidences.
    denominator = 3*(rooted - diagonals) + 2*diagonals
    assert rooted > diagonals and (denominator, joins - diagonals) == (43, 2)
    units = (6*pairs + 4*denominator - 1)//(4*denominator)
    indices = (4*units + 127)//128
    assert (4*units, indices) == (417936, 3266)
    print(diagonals, "cube diagonals force at least", 4*units,
          "long-cycle edge incidences in at least", indices, "indices")


if __name__ == "__main__":
    main()
