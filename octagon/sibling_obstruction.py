"""Test whether every outer line of a bare graph can belong to an inward sibling pair.

Input: the saved 40 voltages for each labelled trace cover. The exclusion uses
only the bare edge list after building the graph from the voltages.
Run from any directory: python octagon/sibling_obstruction.py
"""
import argparse
import json
from itertools import combinations
from pathlib import Path

import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import shortest_path

from base import EDGES, VI


def outer_geometry(voltage):
    lines = np.array([(32 * VI[a] + s, 32 * VI[b] + (s ^ voltage[k]))
                      for k, (a, b) in enumerate(EDGES) for s in range(32)])
    u, v = lines.T
    graph = csr_matrix((np.ones(len(lines)), (u, v)), shape=(512, 512))
    distance = shortest_path(graph + graph.T, unweighted=True)
    if not np.isfinite(distance).all():
        raise ValueError("expected a connected cover")
    flat = lines.ravel()
    delta = distance.astype(np.int16)[np.ix_(flat, flat)].reshape(1280, 2, 1280, 2).min((1, 3))
    # Half the incidence distance between two line vertices in the subdivided graph.
    line_distance = delta + 1
    np.fill_diagonal(line_distance, 0)
    return lines, line_distance


def sibling_blocks(distance):
    """A row (a,b,c,d) pairs lines at two boundary points on one inward line."""
    blocks = []
    left, right = np.nonzero(np.triu(distance >= 7, k=1))
    for i, (a, b) in enumerate(zip(left, right)):
        far = (distance[a] >= 6) & (distance[b] >= 6)
        partners = np.flatnonzero(far[left[i+1:]] & far[right[i+1:]]) + i + 1
        blocks.extend((a, b, left[j], right[j]) for j in partners)
    return np.array(blocks, dtype=int).reshape(-1, 4)


def missing_line(anchor, blocks, distance):
    """Add the anchor's inward tree. Find an outer line with no possible block."""
    possible = ~np.isin(blocks, anchor).any(axis=1)
    for i, j in combinations(range(4), 2):
        # The new tree has half-distance one within a mate pair, two across.
        via = 1 if i // 2 == j // 2 else 2
        for a, b in combinations(range(4), 2):
            required = 7 if a // 2 == b // 2 else 6
            left, right = blocks[:, a], blocks[:, b]
            for u, v in ((anchor[i], anchor[j]), (anchor[j], anchor[i])):
                length = distance[left, u] + via + distance[v, right]
                possible &= length >= required
    covered = np.zeros(len(distance), dtype=bool)
    covered[anchor] = True
    covered[blocks[possible].ravel()] = True
    missing = np.flatnonzero(~covered)
    return int(missing[0]) if len(missing) else None


def verify_missing(lines, anchor, blocks, line):
    """Check a negative witness by BFS on the literal augmented incidence graph."""
    alternatives = blocks[(blocks == line).any(axis=1) & ~np.isin(blocks, anchor).any(axis=1)]
    if len(alternatives) == 0:
        return
    n = int(lines.max()) + 1
    boundary = n + len(lines)
    sources = np.unique(alternatives)
    u = np.r_[lines.ravel(), n + anchor, boundary, boundary+1]
    v = np.r_[np.repeat(n + np.arange(len(lines)), 2),
              boundary + np.array([0, 0, 1, 1]), boundary+2, boundary+2]
    graph = csr_matrix((np.ones(len(u)), (u, v)), shape=(boundary+3, boundary+3))
    dist = shortest_path(graph + graph.T, indices=n+sources, unweighted=True)[:, n:n+len(lines)]
    index = np.searchsorted(sources, alternatives)
    legal = np.ones(len(alternatives), dtype=bool)
    for a, b in combinations(range(4), 2):
        threshold = 14 if a // 2 == b // 2 else 12
        legal &= dist[index[:, a], alternatives[:, b]] >= threshold
    assert not legal.any(), (anchor, line)


def test(voltage, verify=False):
    lines, distance = outer_geometry(voltage)
    blocks = sibling_blocks(distance)
    # These 40 lines suffice on the supplied inputs. Failure at any single
    # line excludes a graph, irrespective of the choice or order of anchors.
    for line in range(0, len(distance), 32):
        anchors = blocks[(blocks == line).any(axis=1)]
        missing = [missing_line(q, blocks, distance) for q in anchors]
        if None not in missing:
            if verify:
                for q, e in zip(anchors, missing):
                    verify_missing(lines, q, blocks, e)
            return len(blocks), line, missing
    return len(blocks), None, []


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--verify", action="store_true", help="recheck each obstruction by incidence BFS")
    args = parser.parse_args()
    data = json.loads(Path(__file__).with_name("stage4_fp0.json").read_text())
    print("row bare_blocks impossible_anchor first_uncovered_line_per_choice")
    survivors = []
    for row, record in enumerate(data):
        count, line, missing = test(record["vol"], args.verify)
        if line is None:
            survivors.append(row)
        print(row, count, "-" if line is None else line, *missing)
    print("surviving rows", survivors)
