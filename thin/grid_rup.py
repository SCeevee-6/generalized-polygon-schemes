"""Reverse unit propagation, with no solver or learned assumptions."""

from collections import defaultdict, deque


class CNF(list):
    pass


class IDPool:
    def __init__(self):
        self.top = 0

    def id(self):
        self.top += 1
        return self.top


class Formula:
    def __init__(self, clauses):
        self.clauses = []
        self.occurrences = defaultdict(list)
        self.units = []
        self.empty = False
        for clause in clauses:
            self.add(clause)

    def add(self, clause):
        c = tuple(dict.fromkeys(clause))
        i = len(self.clauses)
        self.clauses.append(c)
        for lit in c:
            self.occurrences[lit].append(i)
        if len(c) == 1:
            self.units.append(c[0])
        if not c:
            self.empty = True

    def conflict(self, assumptions=()):
        if self.empty:
            return True
        # Each count is the number of literals not yet assigned false.
        remaining = [len(c) for c in self.clauses]
        assigned = set()
        queue = deque(self.units)
        queue.extend(assumptions)
        while queue:
            lit = queue.popleft()
            if lit in assigned:
                continue
            if -lit in assigned:
                return True
            assigned.add(lit)
            for i in self.occurrences[-lit]:
                remaining[i] -= 1
                if remaining[i] == 0:
                    return True
                if remaining[i] == 1:
                    (unit,) = (v for v in self.clauses[i] if -v not in assigned)
                    queue.append(unit)
        return False


def verify(formula, proof):
    """Replay every clause addition and finish with a unit contradiction."""
    bound = max(abs(lit) for c in formula for lit in c)
    checker = Formula(formula)
    for i, c in enumerate(proof, 1):
        assert all(type(lit) is int and 0 < abs(lit) <= bound for lit in c)
        assert checker.conflict(-lit for lit in c), (i, c, "not RUP")
        checker.add(c)
    assert checker.conflict(), "missing contradiction"
    return len(proof)

