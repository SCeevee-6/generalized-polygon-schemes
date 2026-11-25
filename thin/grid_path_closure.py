"""Check bad-grid direction words by completing short paths and comparing full lines.

Run without arguments for the whole enumeration. Pass a single word for one
replay, or --control for tests against the stored Ree–Tits model.
"""
import argparse
from collections import Counter
from itertools import combinations, permutations, product
import json
from pathlib import Path

import numpy as np

from grid_directions import cases, skeleton


def pair_patterns():
    """The point-distance squares of two ordered three-point lines."""
    squares = []
    for perm in permutations(range(3)):
        unequal = np.array([[j != perm[i] for j in range(3)] for i in range(3)], dtype=int)
        squares.extend((unequal, 3 + unequal))  # coincident or opposite
    for i, j, h in product(range(3), repeat=3):
        squares.append(h + (np.arange(3)[:, None] != i).astype(int) + (np.arange(3) != j))
    assert len({tuple(s.flat) for s in squares}) == 39
    return np.array(squares)


def composition_table(squares):
    """Bit c is set at (a,b) if AB=a, BC=b, AC=c passes the small-tree tests."""
    triples = np.array(list(product(range(39), repeat=3)))
    # Each line contributes its three points and its line vertex, in that order.
    distance = np.zeros((12, 12, len(triples)), dtype=np.int8)
    for start in (0, 4, 8):
        for i, j in combinations(range(4), 2):
            distance[start+i, start+j] = distance[start+j, start+i] = 1 if j == 3 else 2
    for t, (a, b) in enumerate(((0, 4), (4, 8), (0, 8))):
        square = squares[triples[:, t]]
        block = np.empty((4, 4, len(triples)), dtype=np.int8)
        block[:3, :3] = 2 * square.transpose(1, 2, 0)
        block[3, :3] = 1 + 2 * square.min(axis=1).T
        block[:3, 3] = 1 + 2 * square.min(axis=2).T
        equal_lines = (square == 0).sum(axis=(1, 2)) == 3
        block[3, 3] = np.where(equal_lines, 0, 2 + 2 * square.min(axis=(1, 2)))
        distance[a:a+4, b:b+4] = block
        distance[b:b+4, a:a+4] = block.transpose(1, 0, 2)

    possible = np.ones(len(triples), dtype=bool)
    for i, j, k in combinations(range(12), 3):
        x, y, z = distance[i, j], distance[i, k], distance[j, k]
        possible &= (x+y >= z) & (x+z >= y) & (y+z >= x)
    for i, j, k, l in combinations(range(12), 4):
        sums = np.sort([distance[i, j] + distance[k, l],
                        distance[i, k] + distance[j, l],
                        distance[i, l] + distance[j, k]], axis=0)
        possible &= (sums[0] + sums[1] >= 16) | (sums[1] == sums[2])

    table = np.zeros((39, 39), dtype=np.uint64)
    a, b, c = triples[possible].T
    np.bitwise_or.at(table, (a, b), np.uint64(1) << c.astype(np.uint64))
    return table


class LineConstraints:
    def __init__(self):
        self.squares = pair_patterns()
        self.compose = composition_table(self.squares)
        self.bits = np.uint64(1) << np.arange(39, dtype=np.uint64)
        self.support = np.empty((9, 32), dtype=np.uint64)
        for k, values in enumerate(self.squares.reshape(39, 9).T):
            for mask in range(32):
                self.support[k, mask] = self.bits[((1 << values) & mask) != 0].sum()

    def propagate(self, masks, lines, line_domains=None):
        """Contract point distances and line-pair domains to stability or failure."""
        lines = np.asarray(lines)
        pairs = np.full((len(lines), len(lines)), (1 << 39) - 1, dtype=np.uint64)
        if line_domains is not None:
            n = len(line_domains)
            assert line_domains.shape == (n, n) and n <= len(lines)
            pairs[:n, :n] &= line_domains
        cells = [np.ix_(lines[:, i], lines[:, j]) for i, j in product(range(3), repeat=2)]
        while True:
            old_pairs, old_masks = pairs.copy(), masks.copy()
            for k, cell in enumerate(cells):
                pairs &= self.support[k, masks[cell]]
            if not pairs.all():
                return False
            for middle in range(len(lines)):
                values, index = np.unique(np.r_[pairs[:, middle], pairs[middle, :]],
                                          return_inverse=True)
                members = (values[:, None] & self.bits) != 0
                # Compose all left and right domains that occur at this middle line.
                right = np.array([np.bitwise_or.reduce(self.compose[:, s], axis=1)
                                  for s in members]).T
                joined = np.array([np.bitwise_or.reduce(right[s], axis=0) for s in members])
                pairs &= joined[np.ix_(index[:len(lines)], index[len(lines):])]
                if not pairs.all():
                    return False
            for k, cell in enumerate(cells):
                allowed = sum(((pairs & self.support[k, 1 << d]) != 0).astype(np.uint8) << d
                              for d in range(5))
                np.bitwise_and.at(masks, cell, allowed)
            masks &= masks.T
            if not masks.all():
                return False
            if np.array_equal(old_pairs, pairs) and np.array_equal(old_masks, masks):
                return True


def extend_paths(masks, lines, requests):
    """Insert full lines on shortest paths forced before this batch."""
    n = len(masks)
    definitions = []
    segments = []
    result = list(lines)
    edges = {tuple(sorted(edge)) for line in lines for edge in combinations(line, 2)}
    for a, b in requests:
        assert 0 <= a < n and 0 <= b < n and masks[a, b] in (4, 8)
        d = int(masks[a, b]).bit_length() - 1
        path = [a]
        for k in range(1, d):
            known = np.flatnonzero((masks[a] == 1 << k) & (masks[b] == 1 << (d-k)))
            if len(known):
                path.append(int(known[0]))
            else:
                path.append(n + len(definitions))
                definitions.append(("path", a, b, k))
        path.append(b)
        for i, j in zip(path, path[1:]):
            edge = tuple(sorted((i, j)))
            if edge not in edges:
                edges.add(edge)
                result.append((i, j, n + len(definitions)))
                definitions.append(("third", i, j, 1))
        for i, j in combinations(range(len(path)), 2):
            segments.append((path[i], path[j], j-i))

    enlarged = np.full((n + len(definitions),) * 2, 31, dtype=np.uint8)
    enlarged[:n, :n] = masks
    np.fill_diagonal(enlarged, 1)
    for i, j, d in segments:
        enlarged[i, j] &= 1 << d
        enlarged[j, i] &= 1 << d
    return enlarged, result, definitions


def replay(directions, constraints, batches):
    seed = skeleton(directions)
    if seed is None:
        assert not batches
        return "coincident partners", 0, 0
    profiles, lines = seed
    assert len(np.unique(profiles, axis=0)) == len(profiles) == 81
    masks = np.full((81, 81), 30, dtype=np.uint8)  # the 81 seed points are distinct
    np.fill_diagonal(masks, 1)
    masks[:, :9] = 1 << profiles
    masks[:9, :] = masks[:, :9].T
    if not constraints.propagate(masks, lines):
        assert not batches
        return "seed contradiction", 81, 51
    for stage, requests in enumerate(batches, 1):
        masks, lines, _ = extend_paths(masks, lines, requests)
        if not constraints.propagate(masks, lines):
            assert stage == len(batches)
            return "path contradiction", len(masks), len(lines)
    raise AssertionError(f"no contradiction for {directions}")


def model_control(constraints):
    from check_opposite_seed import ROOT, octagon

    pd, lines = octagon(np.load(ROOT / "data/generators.npy"))
    # Fix AB in every line-distance class, permute both, then try every C.
    code = (5 ** np.arange(9)).reshape(3, 3)
    lookup = np.full(5**9, -1, dtype=int)
    lookup[(constraints.squares * code).sum(axis=(1, 2))] = np.arange(39)

    def against_all(line):
        distances = pd[line[:, None, None], lines[None, :, :]].transpose(1, 0, 2)
        ids = lookup[(distances * code).sum(axis=(1, 2))]
        assert (ids >= 0).all()
        return ids

    a = lines[0]
    heights = pd[a[:, None, None], lines].min(axis=(0, 2)) + 1
    heights[0] = 0
    seen = set()
    trials = 0
    for h in range(5):
        b = lines[np.flatnonzero(heights == h)[0]]
        for perm in permutations(range(3)):
            first = a[list(perm)]
            ac = against_all(first)
            for other in permutations(range(3)):
                second = b[list(other)]
                bc = against_all(second)
                ab = lookup[(pd[np.ix_(first, second)] * code).sum()]
                assert ab >= 0 and ((constraints.compose[ab, bc] & constraints.bits[ac]) != 0).all()
                seen.add(int(ab))
                trials += len(lines)
    assert len(seen) == 39

    # Two geodesics from p share their inner point u, omitted from the seed.
    p = 0
    u = int(np.flatnonzero(pd[p] == 1)[0])
    dest = [int(l[0] if l[0] != u else l[1]) for l in lines if u in l and p not in l][:2]
    selected = [next(l for l in lines if v in l and u not in l) for v in [p, *dest]]
    selected.append(selected[0][::-1])
    actual = sorted(set(np.asarray(selected).flat))
    old_lines = [tuple(actual.index(v) for v in l) for l in selected]
    masks = (1 << pd[np.ix_(actual, actual)]).astype(np.uint8)
    assert constraints.propagate(masks, old_lines)
    requests = [(actual.index(p), actual.index(v)) for v in dest]
    masks, expanded, definitions = extend_paths(masks, old_lines, requests)
    for kind, i, j, k in definitions:
        x, y = actual[i], actual[j]
        other_d = 1 if kind == "third" else int(pd[x, y]) - k
        choices = np.flatnonzero((pd[x] == k) & (pd[y] == other_d))
        assert len(choices) == 1
        actual.append(int(choices[0]))
    assert actual.count(u) == 2
    assert constraints.propagate(masks, expanded)
    assert (masks & (1 << pd[np.ix_(actual, actual)])).all()
    print("model line triples", trials, "ordered pair patterns", len(seen),
          "shared geodesic point preserved")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("word", nargs="?")
    parser.add_argument("--control", action="store_true")
    args = parser.parse_args()
    constraints = LineConstraints()
    if args.control:
        model_control(constraints)
        return
    recipes = json.loads(Path(__file__).with_name("grid_path_certificates.json").read_text())
    if args.word:
        assert len(args.word) == 9 and set(args.word) <= set("012")
        print(args.word, *replay(tuple(map(int, args.word)), constraints, recipes.get(args.word, [])))
        return
    orbits, totals, used = Counter(), Counter(), set()
    for directions, size in cases():
        word = "".join(map(str, directions))
        reason, _, _ = replay(directions, constraints, recipes.get(word, []))
        if reason == "path contradiction":
            used.add(word)
        orbits[reason] += 1
        totals[reason] += size
    assert used == set(recipes)
    assert sum(orbits.values()) == 374 and sum(totals.values()) == 3**9
    for reason in orbits:
        print(reason, "classes", orbits[reason], "labelings", totals[reason])


if __name__ == "__main__":
    main()
