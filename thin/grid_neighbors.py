"""Find additional nine-tip grids among the 51 prescribed lines."""

from collections import Counter
from itertools import combinations
import json
from pathlib import Path

import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import shortest_path

from grid_directions import OFF, SYMMETRIES, matching_type, skeleton


def named_grids(lines):
    on_line = [set(l) for l in lines]
    links = {}
    neighbours = [set() for _ in lines]
    for i, j in combinations(range(len(lines)), 2):
        if on_line[i] & on_line[j]:
            continue
        join = [l for l in on_line if l & on_line[i] and l & on_line[j]]
        if len(join) == 1:
            link = join[0]
            (a,) = link & on_line[i]
            (b,) = link & on_line[j]
            (t,) = link - {a, b}
            links[i, j], links[j, i] = (a, b, t), (b, a, t)
            neighbours[i].add(j)
            neighbours[j].add(i)
    for rows in combinations(range(len(lines)), 3):
        candidates = set.intersection(*(neighbours[i] for i in rows))
        for cols in combinations(sorted(candidates), 3):
            if rows >= cols or set(rows).intersection(cols):
                continue
            nine = [links[r, c] for r in rows for c in cols]
            if len(set().union(*nine)) == 27:
                yield rows, cols, [t for a, b, t in nine]


def restrictions(tips, lines, dist):
    fixed = {}
    for x, tip in enumerate(tips):
        for line in lines:
            if tip not in line:
                continue
            a, b = set(line) - {tip}
            for y, z in combinations(OFF[x], 2):
                if (
                    dist[a, tips[y]] == dist[b, tips[z]] == 2
                    or dist[b, tips[y]] == dist[a, tips[z]] == 2
                ):
                    required = matching_type(y, z)
                    assert fixed.get(x, required) == required
                    fixed[x] = required
    return fixed


def main():
    reps = json.loads(Path(__file__).with_name("grid_metrics.json").read_text())
    labelled = {}
    for word in reps:
        for perm, swap in SYMMETRIES:
            key = tuple((1, 0, 2)[int(word[i])] if swap else int(word[i]) for i in perm)
            labelled[key] = word
    assert len(labelled) == 788
    counts = {}
    allowed = set()
    for word in reps:
        _, lines = skeleton(tuple(map(int, word)))
        adjacency = np.zeros((81, 81), dtype=int)
        for line in lines:
            for i, j in combinations(line, 2):
                adjacency[i, j] = adjacency[j, i] = 1
        dist = shortest_path(csr_matrix(adjacency), unweighted=True)
        grids = list(named_grids(lines))
        assert any(tips == list(range(9)) for _, _, tips in grids)
        for _, _, tips in grids:
            if tips == list(range(9)):
                continue
            fixed = restrictions(tips, lines, dist)
            works = Counter(
                w for key, w in labelled.items() if all(key[i] == v for i, v in fixed.items())
            )
            assert len(fixed) == 4 and sum(works.values()) == 22
            allowed.add(tuple(sorted(works.items())))
        if len(grids) > 1:
            counts[word] = len(grids) - 1
    (multiplicities,) = allowed
    assert set(dict(multiplicities)) == set(counts)
    print("direction class, additional grids, allowed labelings of each additional grid")
    for word, mult in multiplicities:
        print(word, counts[word], mult)


if __name__ == "__main__":
    main()
