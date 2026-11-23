"""Independent final two lifts, starting from four.txt, without a spectral filter.

Only the base-vertex and voltage-bit numbering are shared with the published
records. Python integers implement GF(2) elimination. Candidate lifts are
checked by counting visits in the actual lifted graphs.
"""
from collections import Counter, defaultdict, deque
from itertools import combinations
from pathlib import Path


def record_numbering():
    words = [i for i in range(32) if i.bit_count() % 2 == 0]
    edges = [(i, j) for i, j in combinations(range(16), 2)
             if (words[i] ^ words[j]).bit_count() == 4]
    # The records use a breadth-first tree with increasing neighbour indices.
    seen, queue, tree = {0}, [0], set()
    for u in queue:
        for e in edges:
            if u in e:
                v = sum(e) - u
                if v not in seen:
                    seen.add(v)
                    queue.append(v)
                    tree.add(e)
    chords = [edge for edge in edges if edge not in tree]
    return words, edges, queue, [0 if e in tree else 1 << chords.index(e) for e in edges]


def read_quotient(cyclic, low, high):
    _, edges, _, masks = record_numbering()
    lifted = []
    for e, mask in zip(edges, masks):
        voltage = int(bool(low & mask)) + 2 * int(bool(high & mask))
        for s in range(4):
            dest = (s + voltage) % 4 if cyclic else s ^ voltage
            lifted.append((4 * e[0] + s, 4 * e[1] + dest))
    return lifted


def adjacency(edges):
    graph = [[] for _ in range(2 * len(edges) // 5)]
    for idx, (u, v) in enumerate(edges):
        graph[u].append((v, idx))
        graph[v].append((u, idx))
    return graph


def chord_bits(graph, edges, breadth_first=False):
    # Use a depth-first spanning tree by default.
    tree, visited = set(), set()
    pending = deque([(0, None)])
    take = pending.popleft if breadth_first else pending.pop
    while pending:
        vertex, incoming = take()
        if vertex not in visited:
            visited.add(vertex)
            if incoming is not None:
                tree.add(incoming)
            pending.extend(graph[vertex])
    assert len(visited) == len(graph)
    free = [i for i in range(len(edges)) if i not in tree]
    bits = [0] * len(edges)
    for j, i in enumerate(free):
        bits[i] = 1 << j
    return bits, len(free)


def visit_groups(edges, graph, bits):
    for idx, (u, v) in enumerate(edges):
        layer = [(u, idx, 0), (v, idx, bits[idx])]
        by_vertex = defaultdict(list)
        for depth in range(4):
            for node, _, sign in layer:
                by_vertex[node].append(sign)
            if depth < 3:
                layer = [(w, e, sign ^ bits[e]) for node, prev, sign in layer
                         for w, e in graph[node] if e != prev]
        yield from by_vertex.values()


def solutions(rows, degree):
    # Rows encode (coefficient mask << 1) | right-hand side.
    rows = list(set(rows))
    pivot_rows = {}
    for j in range(1, degree + 1):
        bit = 1 << j
        index = next((i for i, row in enumerate(rows) if row & bit), None)
        if index is None:
            continue
        pivot = rows.pop(index)
        rows = [row ^ pivot if row & bit else row for row in rows]
        pivot_rows[j] = pivot
    if 1 in rows:
        return
    free = [j for j in range(1, degree + 1) if j not in pivot_rows]
    for assignment in range(1 << len(free)):
        vector = 1 + sum(((assignment >> i) & 1) << j for i, j in enumerate(free))
        for j, row in reversed(pivot_rows.items()):
            if (row & vector).bit_count() % 2:
                vector |= 1 << j
        yield vector >> 1


def admissible_lifts(edges):
    graph = adjacency(edges)
    capacity = 256 // len(graph)
    bits, degree = chord_bits(graph, edges)
    equations = set()
    for group in visit_groups(edges, graph, bits):
        if len(group) > capacity:
            return
        if len(group) == capacity:
            mask = 0
            for bit in group:
                mask ^= bit
            equations.add((mask << 1) | ((capacity // 2) % 2))
    for solution in solutions(equations, degree):
        lifted = []
        for (u, v), mask in zip(edges, bits):
            flip = (solution & mask).bit_count() % 2
            lifted.extend(((2*u, 2*v + flip), (2*u+1, 2*v + 1-flip)))
        actual = visit_groups(lifted, adjacency(lifted), [0] * len(lifted))
        if all(len(group) <= capacity // 2 for group in actual):
            yield lifted


def bipartite(edges):
    graph = adjacency(edges)
    colours = {}
    for root in range(len(graph)):
        if root in colours:
            continue
        colours[root] = 0
        queue = [root]
        for v in queue:
            for w, _ in graph[v]:
                if w in colours:
                    if colours[w] == colours[v]:
                        return False
                else:
                    colours[w] = 1 - colours[v]
                    queue.append(w)
    return True


if __name__ == "__main__":
    count = Counter()
    for record in Path(__file__).with_name("four.txt").read_text().splitlines():
        cyclic, low, high = map(int, record.split())
        count[cyclic, 4] += 1
        for eight in admissible_lifts(read_quotient(cyclic, low, high)):
            count[cyclic, 8] += 1
            for sixteen in admissible_lifts(eight):
                count[cyclic, 16] += 1
                count[cyclic, "bipartite"] += bipartite(sixteen)
    print("quotient four eight sixteen bipartite")
    for cyclic, name in enumerate(("elementary", "cyclic")):
        print(name, *(count[cyclic, k] for k in (4, 8, 16, "bipartite")))
