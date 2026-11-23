"""Allow an arbitrary second double cover of every admissible double quotient."""
import argparse
from collections import Counter
from pathlib import Path

import numpy as np

from check_quotients import BITS, EDGES, upper_signings
from check_terminal import adjacency, admissible_lifts, bipartite, chord_bits, visit_groups


def double(edges, masks, signing):
    return [(2*u+s, 2*v+(s ^ int(bool(signing & mask))))
            for (u, v), mask in zip(edges, masks) for s in (0, 1)]


def forbidden_eigenvalue(edges, bits, signing):
    """Find an integer residual certificate inside a forbidden spectral interval."""
    n = 2 * len(edges) // 5
    matrix = np.zeros((n, n), dtype=np.int64)
    for (u, v), bit in zip(edges, bits):
        matrix[u, v] = matrix[v, u] = 1 - 2 * bool(signing & bit)
    eigenvalues, vectors = np.linalg.eigh(matrix)
    scale = 1000
    for value, vector in zip(eigenvalues, vectors.T):
        t = int(round(value * scale))
        # Open intervals below -4, between -2 and -1, and between 3.45 and 5.
        margin = max(-4000-t, min(t+2000, -1000-t), min(t-3450, 5000-t))
        if margin <= 0:
            continue
        w = np.rint(1000000 * vector).astype(np.int64)
        error = scale * (matrix @ w) - t * w
        # An eigenvalue lies within ||error|| / (scale ||w||) of t/scale.
        # Python integers in both square sums. A missed witness keeps the lift.
        if sum(int(e)**2 for e in error) < margin**2 * sum(int(z)**2 for z in w):
            return True
    return False


def sat_lifts(edges):
    graph = adjacency(edges)
    bits, degree = chord_bits(graph, edges, breadth_first=True)
    groups = ([(mask, 0) for mask in group] for group in visit_groups(edges, graph, bits))
    for signing in upper_signings(groups, degree, 128 // len(graph)):
        yield double(edges, bits, signing)


def census(low, sat_terminal=False):
    base = double(EDGES, BITS, low)
    graph = adjacency(base)
    bits, degree = chord_bits(graph, base, breadth_first=sat_terminal)
    assert degree == 49
    groups = ([(m, 0) for m in group] for group in visit_groups(base, graph, bits))
    counts = Counter()
    lifts = sat_lifts if sat_terminal else admissible_lifts
    for high in upper_signings(groups, degree):
        four = double(base, bits, high)
        counts['four'] += 1
        if forbidden_eigenvalue(base, bits, high):
            continue
        counts['window'] += 1
        for eight in lifts(four):
            counts['eight'] += 1
            for sixteen in lifts(eight):
                counts['sixteen'] += 1
                if bipartite(sixteen):
                    counts['bipartite'] += 1
                else:
                    raise AssertionError(f'nonbipartite terminal graph above {low}, {high}')
    return counts


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sat-terminal", action="store_true", help="use breadth-first tree gauges and SAT at all splits")
    args = parser.parse_args()
    total = Counter()
    fields = ('four', 'window', 'eight', 'sixteen', 'bipartite')
    print('first four window eight sixteen bipartite', flush=True)
    for low in map(int, Path(__file__).with_name('first.txt').read_text().split()):
        counts = census(low, args.sat_terminal)
        print(low, *(counts[k] for k in fields), flush=True)
        total.update(counts)
    print('total', *(total[k] for k in fields))
