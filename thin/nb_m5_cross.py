"""Cross distances behind the exact fifth moment (nb_moment_sdp.md)."""

from collections import Counter
from itertools import combinations

import numpy as np

import failed_pair_types as fpt
import grid_path_closure as gpc
import theta_profiles as tp

REQUESTS = list(combinations(tp.MIDDLE, 2))


def zline(idx):
    s = 21 + 5*idx
    return (s, s + 1, s + 3)


def prepare(constraints, base_lines, d):
    """Add the geodesics a_k -> m and b_k -> m whose middle lines are Y_A, Y_B."""
    masks = (1 << d.astype(int)).astype(np.uint8)
    requests, meta = [], []
    for idx, (ti, tj) in enumerate(REQUESTS):
        k = 3 - tp.MIDDLE.index(ti) - tp.MIDDLE.index(tj)
        m = zline(idx)[2]
        for side, point in (("A", 7*k), ("B", 7*k + 3)):
            assert masks[point, m] == 8
            requests.append((point, m))
            meta.append((idx, side))
    big, lines, defs = gpc.extend_paths(masks, base_lines, requests)
    n, ys = len(masks), {}
    for (a, b), key in zip(requests, meta):
        path = []
        for k in (1, 2):
            known = np.flatnonzero((masks[a] == 1 << k) & (masks[b] == 1 << (3 - k)))
            path.append(int(known[0]) if len(known) else
                        n + [q for q, x in enumerate(defs) if x == ("path", a, b, k)][0])
        ys[key], = [line for line in lines if path[0] in line and path[1] in line]
    assert constraints.propagate(big, lines)
    return big, lines, ys


def line_values(constraints, big, lines, first, second):
    """Line distances 3 and 4 that survive propagation (1 + min point distance)."""
    found = set()
    trial = big.copy()
    for a in first:
        for b in second:
            trial[a, b] = trial[b, a] = trial[a, b] & np.uint8(8 | 16)
    if constraints.propagate(trial, lines):
        found.add(4)
    for a in first:
        for b in second:
            if big[a, b] & 4:
                trial = big.copy()
                trial[a, b] = trial[b, a] = 4
                if constraints.propagate(trial, lines):
                    found.add(3)
                    return found
    return found


def main():
    constraints = gpc.LineConstraints()
    masks, base_lines, _ = gpc.extend_paths((1 << tp.seed_distances()).astype(np.uint8),
                                            tp.LINES, REQUESTS)
    pattern = Counter()
    for d in fpt.completions(constraints, masks.copy(), base_lines):
        labels = tuple(int(x) for x in (d[19, 24], d[12, 29], d[5, 34]))
        if labels == (1, 1, 1):
            continue
        big, lines, ys = prepare(constraints, base_lines, d)
        for idx in range(3):
            for side in "AB":
                # Y_A and Y_B are middle lines of pairs containing Z_idx.
                y = ys[idx, side]
                assert 1 + min(int(big[a, b]).bit_length() - 1
                               for a in zline(idx) for b in y) == 2
            values = tuple(tuple(sorted(line_values(constraints, big, lines, zline(o), ys[idx, s])))
                           for o in range(3) if o != idx for s in "AB")
            pattern[sorted(labels) == [3, 4, 4], values.count((4,))] += 1
    assert pattern == {(False, 0): 24, (True, 0): 24, (True, 3): 48}
    print("open cross distances per directed edge:", dict(pattern))
    # Exact fifth moment: Q_5 = 36 Q_1 + 24 (RS + RS^T) + 4 S^T J S - 9 Q_3 - 29 Q_1.
    for s in (0, 1, 4):
        q5 = [0, 0, 98304, 12288*s, 98304 - 12288*s]
        q1, q3 = [0, 0, 12288, 0, 0], [0, 0, 36864, 0, 12288]
        stjs = [0, 0, 12288, 3072*s, 3072*(5 - s)]
        rs = [0, 0, 12288, 0, 6144]
        built = [36*a + 24*b + 4*c - 9*x - 29*y for a, b, c, x, y in zip(q1, rs, stjs, q3, q1)]
        assert built == q5 and sum(q5) == 12*2**4*1024
    print("Q_5 = [0, 0, 98304, 12288 s, 98304 - 12288 s] per line")


if __name__ == "__main__":
    main()
