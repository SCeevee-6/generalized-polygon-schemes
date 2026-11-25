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
    return moves


def braid_checks(moves):
    """Fill flag edges from a fixed tree, taking unconditional relations first."""
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
    tests = braid_checks(frame_graph(incident, dual))
    check_point_graph(v, w, incident, dual)
    pairs = 2925*2048//2
    systems = (6*pairs + 128*tests-1)//(128*tests)
    assert systems == 1448
    print("1152 reference frames,", tests, "eight-cycle tests")
    print("Two four-point cycles verified on 144 named points")
    print("At least", systems, "indexed systems contain a cycle of length at least 12")


if __name__ == "__main__":
    main()
