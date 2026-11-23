"""Independent double- and four-sheet enumeration. Requires numpy and python-sat.

The exact polynomial check uses principal submatrices (Berkowitz recurrence).
The four-sheet search assigns the upper voltage bits directly with SAT.
Compare only at the end with the published lists of records.
"""
from collections import Counter, defaultdict
from itertools import combinations, permutations
from math import comb
from pathlib import Path

import numpy as np
from pysat.solvers import Solver

from check_terminal import record_numbering

WORDS, EDGES, ORDER, BITS = record_numbering()
CHORDS = [j for j, bit in enumerate(BITS) if bit]
NEIGHBOURS = [[] for _ in WORDS]
for k, (u, v) in enumerate(EDGES):
    NEIGHBOURS[u].append((v, k, 1))
    NEIGHBOURS[v].append((u, k, -1))


def orbit_table():
    # Evaluate pullback on the 25 fundamental cycles under every base isometry.
    paths = [0] * 16
    for v in ORDER[1:]:
        u, edge, _ = next((u, k, s) for u, k, s in NEIGHBOURS[v]
                          if BITS[k] == 0 and ORDER.index(u) < ORDER.index(v))
        paths[v] = paths[u] ^ (1 << edge)
    cycles = [paths[u] ^ paths[v] ^ (1 << k) for k in CHORDS for u, v in [EDGES[k]]]
    maps = []
    for p in permutations(range(5)):
        perm = [sum(((v >> i) & 1) << p[i] for i in range(5)) for v in WORDS]
        for tr in WORDS:
            g = [WORDS.index(v ^ tr) for v in perm]
            images = [BITS[EDGES.index(tuple(sorted((g[u], g[v]))))] for u, v in EDGES]
            masks = [sum(images[k] for k in range(40) if cycle >> k & 1) for cycle in cycles]
            maps.append([sum(1 << i for i, mask in enumerate(masks) if mask >> j & 1)
                         for j in range(25)])
    cols = np.array(maps, dtype=np.int32).T
    table = np.zeros((4, 256, len(maps)), dtype=np.int32)
    for j, col in enumerate(cols):
        size = 1 << (j % 8)
        table[j//8, size:2*size] = table[j//8, :size] ^ col
    return table


def orbit(table, word):
    return (table[0, word & 255] ^ table[1, (word >> 8) & 255]
            ^ table[2, (word >> 16) & 255] ^ table[3, word >> 24])


def charpolys(matrices):
    # det(tI-M) = (t-a)det(tI-B) - row adj(tI-B) col for the trailing
    # principal block B. The series coefficients of row (tI-B)^(-1) col
    # are row B^k col. All matrices here have operator norm at most five.
    n = matrices.shape[1]
    co = np.ones((len(matrices), 1), dtype=np.int64)
    for i in reversed(range(n)):
        size = n - i
        block = matrices[:, i+1:, i+1:]
        row = matrices[:, i, i+1:]
        power = row.copy()
        new = np.zeros((len(matrices), size+1), dtype=np.int64)
        new[:, :size] = co
        new[:, 1:] -= matrices[:, i, i, None] * co
        for j in range(2, size+1):
            scalar = (power * row).sum(axis=1)
            new[:, j:] -= scalar[:, None] * co[:, :size+1-j]
            if j < size:
                power = (power[:, None, :] @ block)[:, 0, :]
        co = new
    return co


def roots_above(co, a, b=1):
    # For a real-rooted p, sign variations of p(a/b + x) count positive
    # roots exactly. Multiply by b^degree; the first nonzero index gives
    # the multiplicity at the endpoint. Python integers for all shifts.
    asc = list(reversed(co))
    shifted = [sum(asc[k] * comb(k, j) * a**(k-j) * b**(len(co)-1-k)
                   for k in range(j, len(co))) for j in range(len(co))]
    nz = [v for v in shifted if v]
    return sum((a > 0) != (b > 0) for a, b in zip(nz, nz[1:])), next(i for i, v in enumerate(shifted) if v)


def spectral_ok(co):
    lower, at_lower = roots_above(co, -4)
    after_gap, at_end = roots_above(co, -1)
    return (lower + at_lower == 16 and roots_above(co, -2)[0] == after_gap + at_end
            and roots_above(co, 69, 20)[0] == 0)


def tree_visits(low, cyclic):
    # Walk directly on the oriented base. With upper voltage b, the sheet
    # reached is constant + 2*parity(mask & b) in the cyclic case. The xor
    # case has the same variable upper bit and no carry.
    for edge, (u, v) in enumerate(EDGES):
        levels = [(u, edge, 0, 0), (v, edge, int(bool(low & BITS[edge])), BITS[edge])]
        groups = defaultdict(list)
        for depth in range(4):
            for node, _, const, mask in levels:
                groups[node, const % 2].append((mask, const // 2))
            if depth < 3:
                levels = [(w, e, (const + sign * bool(low & BITS[e])) % 4 if cyclic
                           else const ^ bool(low & BITS[e]), mask ^ BITS[e])
                          for node, prev, const, mask in levels
                          for w, e, sign in NEIGHBOURS[node] if e != prev]
        yield from groups.values()


def first_stage(table):
    seen = np.zeros(1 << 25, dtype=bool)
    reps = []
    for word in range(len(seen)):
        if not seen[word]:
            seen[orbit(table, word)] = True
            reps.append(word)
    assert seen.all()
    passed = []
    spectrum = {}
    spectral_records = []
    for start in range(0, len(reps), 256):
        batch = np.array(reps[start:start+256])
        matrix = np.zeros((len(batch), 16, 16), dtype=np.int64)
        for (u, v), mask in zip(EDGES, BITS):
            matrix[:, u, v] = matrix[:, v, u] = 1 - 2*((batch & mask) != 0)
        for word, coeff in zip(batch, charpolys(matrix).tolist()):
            co = tuple(coeff)
            if co not in spectrum:
                spectrum[co] = spectral_ok(co)
            if spectrum[co]:
                spectral_records.append(int(word))
                if all(len(group) <= 8 for group in tree_visits(int(word), False)):
                    passed.append(int(word))
    print("double orbits, spectral, capacity", len(reps), len(spectral_records), len(passed))
    return passed


def upper_signings(groups, degree, maximum=4):
    # A group of visits on one sheet of the double quotient must split into
    # at most maximum visits of either upper bit. Forbid larger monochromatic subsets.
    variables = {1 << i: i+1 for i in range(degree)}
    variables[0] = degree + 1
    clauses = {(-(degree + 1),)}

    def parity(mask):
        if mask not in variables:
            bit = mask & -mask
            u, v = parity(mask ^ bit), parity(bit)
            z = len(variables) + 1
            variables[mask] = z
            clauses.update(((u, v, -z), (-u, -v, -z), (u, -v, z), (-u, v, z)))
        return variables[mask]

    for group in groups:
        literals = [(-1 if const else 1) * parity(mask) for mask, const in group]
        for subset in combinations(literals, maximum + 1):
            subset = set(subset)
            if any(-x in subset for x in subset):
                continue
            clauses.add(tuple(sorted(subset)))
            clauses.add(tuple(sorted(-x for x in subset)))
    with Solver(name="cadical195", bootstrap_with=clauses) as solver:
        while solver.solve():
            assignment = solver.get_model()[:degree]
            yield sum(1 << i for i, literal in enumerate(assignment) if literal > 0)
            solver.add_clause([-v for v in assignment])


if __name__ == "__main__":
    directory = Path(__file__).parent
    table = orbit_table()
    first = first_stage(table)
    classes = {int(b): a for a in first for b in orbit(table, a)}
    print("double signings passing capacity", len(classes))
    rows = []
    for a in first:
        for cyclic in (0, 1):
            for b in upper_signings(tree_visits(a, cyclic), 25):
                if not cyclic:
                    c = a ^ b
                    if b >= c or min(classes.get(b, -1), classes.get(c, -1)) < a:
                        continue
                rows.append((cyclic, a, b))
    print("four-sheet records", sorted(Counter(c for c, _, _ in rows).items()))
    assert first == list(map(int, (directory / "first.txt").read_text().split()))
    saved = [tuple(map(int, line.split())) for line in (directory / "four.txt").read_text().splitlines()]
    assert sorted(rows) == sorted(saved)
    print("both record lists agree")
