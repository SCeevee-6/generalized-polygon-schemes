"""Couple full stars by distances between their neighbours."""

from itertools import combinations, product
import json
from pathlib import Path

from grid_directions import SYMMETRIES, third_domains
from grid_star_extension import Neighbourhoods, OBSTRUCTIONS, line_pattern
from grid_middle_extension import ROOTS, add_middle
from grid_rup import verify


ROOT_PAIRS = {
    "000002002": (32, 37),
    "001110221": (33, 36),
    "002012102": (39, 52),
}


class Coupled(Neighbourhoods):
    def __init__(self, word, table):
        super().__init__(word, table)
        self.stars = {}
        self.distance = dict(self.named)
        self.n = 81

    def add_star(self, root):
        assert root not in self.stars
        known = [p for line in self.lines if root in line for p in line if p != root]
        rows = self.add_root(root)
        add_middle(self, root, rows)
        vertices = known + list(range(self.n, self.n + 10 - len(known)))
        self.n += 10 - len(known)
        for a, row in zip(vertices[len(known):], rows[len(known):]):
            for y, v in enumerate(row):
                self.distance[a, y] = self.distance[y, a] = v
            self.distance[a, a] = self.choose([0])
        for i, j in combinations(range(10), 2):
            a, b = vertices[i], vertices[j]
            d = 1 if i // 2 == j // 2 else 2
            if (a, b) in self.distance:
                self.formula.append([self.distance[a, b][d]] if d in self.distance[a, b] else [])
            else:
                self.distance[a, b] = self.distance[b, a] = self.choose([d])
        for s, others in self.stars.items():
            self.couple(root, vertices, s, others)
        self.stars[root] = vertices

    def couple(self, root, vertices, s, others):
        for a in vertices:
            for b in others:
                if (a, b) in self.distance:
                    continue
                # Both are fresh; they may represent the same ambient point.
                values = set(range(5))
                for y in range(81):
                    values &= {d for d in range(5) if any(abs(p - q) <= d <= p + q
                                for p, q in product(self.distance[a, y], self.distance[b, y]))}
                v = self.choose(sorted(values))
                self.distance[a, b] = self.distance[b, a] = v
                for y in range(81):
                    self.require([self.distance[a, y], self.distance[b, y], v],
                                 lambda p, q, d: abs(p - q) <= d <= p + q)
                self.require([self.named[root, s], self.distance[a, s], self.distance[b, root], v],
                             geodesic_neighbours)
        for verts, r, obs in ((vertices, root, others), (others, s, vertices)):
            for a in obs:
                if a < 81:
                    continue  # already imposed for every named observer
                h = self.distance[a, r]
                for i in range(0, 10, 2):
                    self.require([h, self.distance[a, verts[i]], self.distance[a, verts[i + 1]]],
                                 line_pattern)
                for d, selected in h.items():
                    if 0 < d < 4:
                        inward = [self.distance[a, b][d - 1] for b in verts
                                  if d - 1 in self.distance[a, b]]
                        self.formula.append([-selected, *inward])
                        for x, y in combinations(inward, 2):
                            self.formula.append([-selected, -x, -y])


def geodesic_neighbours(d, p, q, h):
    """Distances of neighbours of two roots, when the roots are near."""
    # p=d(a,s), q=d(b,r), with a neighbour of r and b neighbour of s.
    if d == 1:
        # The geodesic r,s is a single line. Relative to the other root,
        # the inward and lateral neighbours are precisely the two named
        # points of this line; other neighbours are on disjoint branches.
        if p == 0 or q == 0:
            return h == 1
        if p == q == 1:
            return h == 0
        return h == p + q - 1
    if d in (2, 3):
        length = p + q - d
        if length <= 4:
            return h == length
        return h >= 8 - length
    return True


def main():
    saved = json.loads(Path(__file__).with_name("grid_metrics.json").read_text())
    proofs = json.loads(Path(__file__).with_name("grid_coupled_refutations.json").read_text())
    excluded = {word for word, (_, condition) in OBSTRUCTIONS.items() if not condition} | ROOTS.keys()
    assert proofs.keys() == ROOT_PAIRS.keys() <= saved.keys() - excluded
    table = third_domains()
    for word, roots in ROOT_PAIRS.items():
        model = Coupled(word, table)
        for root in roots:
            model.add_star(root)
        steps = verify(model.formula, proofs[word])
        print(word, "excluded at roots", *roots, "RUP additions", steps, flush=True)
    remaining = set(saved) - excluded - ROOT_PAIRS.keys()
    labelled = set()
    for word in remaining:
        for perm, swap in SYMMETRIES:
            labelled.add(tuple((1, 0, 2)[int(word[i])] if swap else int(word[i]) for i in perm))
    print("remaining classes", len(remaining), "labelings", len(labelled))


if __name__ == "__main__":
    main()
