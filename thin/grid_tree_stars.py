"""Tree patterns at the neighbours of a completed star.

For a root r and named b, c with d(r,b), d(r,c) <= 3 and perimeter at most
seven, the three geodesics form a tree. The tree fixes which pairs
(d(a,b)-d(r,b), d(a,c)-d(r,c)) are possible for a neighbour a of r.

With no arguments, replay the saved refutations by reverse unit propagation.
To test one instance: python grid_tree_stars.py WORD ROOT [ROOT2] [--coupled]
[--limit SEC] [--out FILE] (needs python-sat).
"""

from itertools import product
import json
from pathlib import Path
import sys
import time

from grid_directions import third_domains

# Allowed (delta_b, delta_c) for a neighbour of the root.
SAME_POINT = {(-1, -1), (0, 0), (1, 1)}
TWO_LINES = {(-1, 1), (0, 1), (1, -1), (1, 0), (1, 1)}
ONE_LINE = {(-1, 0), (0, -1), (1, 1)}


def tree_case(i, j, k):
    """Case at the root for d(r,b)=i, d(r,c)=j, d(b,c)=k, or None if no tree."""
    if not (1 <= i <= 3 and 1 <= j <= 3 and 1 <= k) or i + j + k > 7:
        return None
    if not (abs(i - j) <= k <= i + j):
        return None
    if (i + j + k) % 2 == 0:
        arm = (i + j - k) // 2
        return SAME_POINT if arm >= 1 else TWO_LINES
    arm = (i + j - k - 1) // 2
    return SAME_POINT if arm >= 1 else ONE_LINE


def add_tree(model, root, rows, named=None):
    """Impose the tree patterns at every neighbour row of root.

    rows[u][y] is the distance variable from the u-th neighbour to named y.
    Returns the number of clauses added.
    """
    named = model.named if named is None else named
    before = len(model.formula)
    for b in range(81):
        if b == root:
            continue
        for c in range(b + 1, 81):
            if c == root:
                continue
            for (i, li), (j, lj), (k, lk) in product(
                named[root, b].items(), named[root, c].items(), named[b, c].items()
            ):
                allowed = tree_case(i, j, k)
                if allowed is None:
                    continue
                antecedent = [-li, -lj, -lk]
                for row in rows:
                    for (x, u), (y, v) in product(row[b].items(), row[c].items()):
                        if (x - i, y - j) not in allowed:
                            model.formula.append(antecedent + [-u, -v])
    return len(model.formula) - before


def add_tree_general(model, root, verts, targets):
    """Tree patterns at the star of root against arbitrary target vertices."""
    dist = model.distance
    before = len(model.formula)
    for x, b in enumerate(targets):
        for c in targets[x + 1:]:
            for (i, li), (j, lj), (k, lk) in product(
                dist[root, b].items(), dist[root, c].items(), dist[b, c].items()
            ):
                allowed = tree_case(i, j, k)
                if allowed is None:
                    continue
                antecedent = [-li, -lj, -lk]
                for a in verts:
                    if a in (b, c):
                        continue
                    for (p, u), (q, v) in product(dist[a, b].items(), dist[a, c].items()):
                        if (p - i, q - j) not in allowed:
                            model.formula.append(antecedent + [-u, -v])
    return len(model.formula) - before


def build(word, roots, coupled=False, cross=True):
    from grid_star_extension import Neighbourhoods
    from grid_middle_extension import add_middle

    table = third_domains()
    if coupled:
        from grid_coupled_stars import Coupled

        model = Coupled(word, table)
        for r in roots:
            model.add_star(r)
        for r, verts in model.stars.items():
            targets = [y for y in range(81) if y != r]
            if cross:
                targets += [b for s, other in model.stars.items() if s != r
                            for b in other if b >= 81]
            add_tree_general(model, r, verts, targets)
    else:
        model = Neighbourhoods(word, table)
        for r in roots:
            rows = model.add_root(r)
            add_middle(model, r, rows)
            add_tree(model, r, rows)
    return model


def solve(formula, limit):
    from pysat.solvers import Solver

    with Solver(name="cd19", bootstrap_with=formula, with_proof=False) as s:
        if limit:
            import threading

            timer = threading.Timer(limit, s.interrupt)
            timer.start()
            res = s.solve_limited(expect_interrupt=True)
            timer.cancel()
        else:
            res = s.solve()
        return res


ROOT_PAIRS = {"000000222": (29, 73)}
REMAINING_BEFORE = """000000000 000000022 000000111 000000222 000011022
001001222 001002111 001012102 002002002""".split()


def replay():
    from grid_directions import SYMMETRIES
    from grid_rup import verify

    proofs = json.loads(Path(__file__).with_name("grid_tree_refutations.json").read_text())
    assert proofs.keys() == ROOT_PAIRS.keys()
    for word, roots in ROOT_PAIRS.items():
        model = build(word, roots, coupled=True)
        steps = verify(model.formula, proofs[word])
        print(word, "excluded at roots", *roots, "RUP additions", steps, flush=True)
    remaining = set(REMAINING_BEFORE) - ROOT_PAIRS.keys()
    labelled = set()
    for word in remaining:
        for perm, swap in SYMMETRIES:
            labelled.add(tuple((1, 0, 2)[int(word[i])] if swap else int(word[i]) for i in perm))
    print("remaining classes", len(remaining), "labelings", len(labelled))
    print(*sorted(remaining))


def main(argv):
    if not argv:
        return replay()
    out = None
    limit = 0
    coupled = False
    args = []
    it = iter(argv)
    for a in it:
        if a == "--out":
            out = next(it)
        elif a == "--limit":
            limit = float(next(it))
        elif a == "--coupled":
            coupled = True
        else:
            args.append(a)
    word, roots = args[0], list(map(int, args[1:]))
    t = time.time()
    model = build(word, roots, coupled)
    nclauses = len(model.formula)
    tb = time.time() - t
    res = solve(model.formula, limit)
    rec = dict(word=word, roots=roots, coupled=coupled, clauses=nclauses,
               build=round(tb, 1), solve=round(time.time() - t - tb, 1),
               result={True: "sat", False: "unsat", None: "unknown"}[res])
    print(json.dumps(rec), flush=True)
    if out:
        with open(out, "a") as f:
            f.write(json.dumps(rec) + "\n")


if __name__ == "__main__":
    main(sys.argv[1:])
