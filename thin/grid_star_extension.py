"""Refute grid direction classes and partial distance assignments by five-line stars.

Run python thin/grid_star_extension.py with python-sat installed.
"""

from itertools import combinations, product
import json
from pathlib import Path

import numpy as np
from grid_rup import CNF, IDPool

try:
    from pysat.solvers import Solver
except ModuleNotFoundError:
    Solver = None

from grid_directions import SYMMETRIES, distance_domains, skeleton, star_domains, third_domains


# Each condition is a conjunction of (point, point, distance) triples.
# Empty conditions exclude the entire direction class.
OBSTRUCTIONS = {
    "000000111": ([27], [(27, 33, 4), (27, 37, 3), (27, 41, 3), (33, 56, 3)]),
    "000011011": ([39], []),
    "000012102": ([77, 36, 37], []),
    "001002111": ([37], [(27, 33, 4), (27, 37, 4), (29, 37, 4), (37, 60, 3)]),
    "001010100": ([29], []),
    "001010102": ([31], []),
    "001010122": ([31], []),
    "001010200": ([29], []),
    "001010202": ([72, 30], []),
    "001020120": ([31], []),
    "001021120": ([31], []),
}


def metric(word, digits, table):
    profiles, lines = skeleton(tuple(map(int, word)))
    masks = distance_domains(profiles, lines, table)
    digits = iter(digits)
    dist = np.zeros((81, 81), dtype=int)
    for i, j in combinations(range(81), 2):
        mask = int(masks[i, j])
        d = mask.bit_length() - 1 if mask.bit_count() == 1 else int(next(digits))
        assert mask & (1 << d)
        dist[i, j] = dist[j, i] = d
    assert next(digits, None) is None
    return dist, lines


def line_pattern(a, b, c):
    low, mid, high = sorted((a, b, c))
    return mid == high == low + 1


class Neighbourhoods:
    """One symmetric partial metric shared by the observed five-line stars."""

    def __init__(self, word, table):
        profiles, self.lines = skeleton(tuple(map(int, word)))
        self.pool = IDPool()
        self.formula = CNF()
        self.named = {}
        masks = star_domains(distance_domains(profiles, self.lines, table), self.lines)
        for i in range(81):
            for j in range(i, 81):
                values = [d for d in range(5) if masks[i, j] & (1 << d)]
                self.named[i, j] = self.named[j, i] = self.choose(values)
        for i in range(81):
            for line in self.lines:
                self.require([self.named[i, j] for j in line], line_pattern)

    def choose(self, values):
        """Represent an integer with one Boolean per allowed value."""
        options = {d: self.pool.id() for d in values}
        self.formula.append(list(options.values()))
        for a, b in combinations(options.values(), 2):
            self.formula.append([-a, -b])
        return options

    def require(self, variables, predicate):
        for assignment in product(*(v.items() for v in variables)):
            if not predicate(*(d for d, lit in assignment)):
                self.formula.append([-lit for d, lit in assignment])

    def add_root(self, root):
        known = [p for line in self.lines if root in line for p in line if p != root]
        assert len(set(known)) == len(known) <= 10
        rows = [[self.named[p, y] for y in range(81)] for p in known]
        for _ in range(10 - len(known)):
            row = []
            for y in range(81):
                h = self.named[root, y]
                row.append(self.choose([d for d in range(5) if (y != root or d == 1)
                                        and any(abs(d - k) <= 1 for k in h)]))
            for line in self.lines:
                self.require([row[y] for y in line], line_pattern)
            for i, j in combinations(range(81), 2):
                self.require([row[i], row[j], self.named[i, j]],
                             lambda a, b, d: abs(a - b) <= d <= a + b)
            rows.append(row)
        for y in range(81):
            h = self.named[root, y]
            for u in range(0, 10, 2):
                self.require([h, rows[u][y], rows[u + 1][y]], line_pattern)
            for d, selected in h.items():
                if 0 < d < 4:
                    inward = [row[y][d - 1] for row in rows if d - 1 in row[y]]
                    self.formula.append([-selected, *inward])
                    for a, b in combinations(inward, 2):
                        self.formula.append([-selected, -a, -b])

        return rows


def check_refutation(formula, proof):
    """Replay clause additions by reverse unit propagation in MiniSat.

    Keeping deleted clauses can only strengthen unit propagation.
    """
    with Solver(name="m22", bootstrap_with=formula) as checker:
        for line in proof:
            if line.startswith("d "):
                continue
            clause = list(map(int, line.split()))
            assert clause.pop() == 0
            consistent, _ = checker.propagate(assumptions=[-v for v in clause])
            assert not consistent, "clause addition failed reverse unit propagation"
            if not clause:
                return
            checker.add_clause(clause)
    raise AssertionError("missing empty clause")


def possible(formula):
    if Solver is None:
        raise RuntimeError("Install python-sat to run the original star refutations.")
    with Solver(name="g4", bootstrap_with=formula, with_proof=True) as solver:
        if solver.solve():
            model = set(solver.get_model())
            assert all(any(lit in model for lit in clause) for clause in formula)
            return True
        check_refutation(formula, solver.get_proof())
        return False


def main():
    saved = json.loads(Path(__file__).with_name("grid_metrics.json").read_text())
    table = third_domains()
    excluded = set()
    for word, (roots, condition) in OBSTRUCTIONS.items():
        model = Neighbourhoods(word, table)
        for root in roots:
            model.add_root(root)
        if condition:
            # Reading a saved assignment is only a cross-check that the cut hits it.
            dist, _ = metric(word, saved[word], table)
            for i, j, d in condition:
                assert dist[i, j] == d
                model.formula.append([model.named[i, j][d]])
        else:
            excluded.add(word)
        assert not possible(model.formula)
        print(word, "excluded" if not condition else "forbidden", condition, flush=True)
    labelled = set()
    for word in excluded:
        for perm, swap in SYMMETRIES:
            labelled.add(tuple((1, 0, 2)[int(word[i])] if swap else int(word[i]) for i in perm))
    print("excluded classes", len(excluded), "labelings", len(labelled))
    print("remaining classes", len(saved) - len(excluded), "labelings", 788 - len(labelled))


if __name__ == "__main__":
    main()

