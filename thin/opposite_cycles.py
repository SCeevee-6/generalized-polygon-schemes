"""Enumerate line-label traces of induced short cycles of opposite pairs."""

from itertools import combinations, product


def line_pair(a, b):
    return tuple(sorted((a, b)))


def constrain(path, closed):
    """Return required distances, or None if a trace is inconsistent."""
    distances = {}
    size = len(path)
    if len(set(path)) != size:
        return None

    def require(pair, distance):
        return pair[0] != pair[1] and distances.setdefault(pair, distance) == distance

    for pair in path:
        if not require(pair, 4):
            return None
    for i, j in combinations(range(size), 2):
        separation = j - i
        if closed:
            separation = min(separation, size - separation)
        if separation == 1:
            if any(not require(line_pair(a, b), 2) for a, b in product(path[i], path[j])):
                return None
        elif separation == 2:
            delta = sorted(set(path[i]) ^ set(path[j]))
            if len(delta) != 2 or not require(tuple(delta), 4):
                return None
    for i, j in combinations(range(size), 2):
        # The ends of an open path might be adjacent when the cycle closes.
        if j == i + 1 or (i, j) == (0, size - 1):
            continue
        if all(distances.get(line_pair(a, b)) == 2 for a, b in product(path[i], path[j])):
            return None
    return distances


def traces(size):
    def extend(path):
        closed = len(path) == size
        distances = constrain(path, closed)
        if distances is None:
            return
        if closed:
            yield path, distances
            return
        fresh = 1 + max(max(pair) for pair in path)
        for inherited, other in product(path[-2], range(fresh + 1)):
            yield from extend(path + (line_pair(inherited, other),))

    yield from extend(((0, 1), (2, 3), (0, 4)))


def attachment_possible(path, distances):
    """Opposite lines at distance two from a reference line project differently."""
    for a in sorted(set().union(*path)):
        nearby = [b for b in sorted(set().union(*path)) if distances.get(line_pair(a, b)) == 2]
        graph = {b: {c for c in nearby if distances.get(line_pair(b, c)) == 4} for b in nearby}
        order = sorted(nearby, key=lambda b: -len(graph[b]))

        def color(todo, assigned):
            if not todo:
                return True
            vertex = todo[0]
            unavailable = {assigned[b] for b in graph[vertex] if b in assigned}
            return any(color(todo[1:], {**assigned, vertex: value})
                       for value in range(3) if value not in unavailable)

        if not color(order, {}):
            return False
    return True


def main():
    expected = {5: (0, 0), 6: (2, 0), 7: (0, 0), 8: (6, 2), 9: (2, 2), 10: (38, 0)}
    for size, counts in expected.items():
        candidates = list(traces(size))
        survivors = [path for path, d in candidates if attachment_possible(path, d)]
        assert (len(candidates), len(survivors)) == counts
        if size == 8:
            assert all(len(set.intersection(*(set(pair) for pair in path[i::2]))) == 1
                       for path in survivors for i in range(2))
        print(size, len(candidates), len(survivors))


if __name__ == "__main__":
    main()
