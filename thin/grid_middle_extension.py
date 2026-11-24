"""Extend the middle-line identity to the unnamed neighbours of a root."""

from itertools import product
import json
from pathlib import Path

from grid_directions import SYMMETRIES, third_domains
from grid_star_extension import Neighbourhoods, OBSTRUCTIONS
from grid_rup import verify


ROOTS = {
    "000000002": 36,
    "000001220": 59,
    "001010222": 29,
    "002002111": 52,
}


def add_middle(model, root, rows):
    # d(root,b)=3, d(root,c)=d(b,c)=2 implies that every neighbour
    # a of root satisfies d(a,b)=d(a,c)+1, including unnamed a.
    named = model.named
    for b in range(81):
        if 3 not in named[root, b]:
            continue
        for c in range(81):
            if 2 not in named[root, c] or 2 not in named[b, c]:
                continue
            antecedent = [-named[root, b][3], -named[root, c][2], -named[b, c][2]]
            for row in rows:
                for (x, u), (y, v) in product(row[b].items(), row[c].items()):
                    if x != y + 1:
                        model.formula.append(antecedent + [-u, -v])


def main():
    saved = json.loads(Path(__file__).with_name("grid_metrics.json").read_text())
    proofs = json.loads(Path(__file__).with_name("grid_middle_refutations.json").read_text())
    old = {word for word, (_, condition) in OBSTRUCTIONS.items() if not condition}
    assert ROOTS.keys() <= saved.keys() - old
    assert proofs.keys() == ROOTS.keys()
    table = third_domains()
    for word, root in ROOTS.items():
        model = Neighbourhoods(word, table)
        rows = model.add_root(root)
        add_middle(model, root, rows)
        steps = verify(model.formula, proofs[word])
        assert steps == len(proofs[word])
        print(word, "excluded at root", root, "RUP additions", steps, flush=True)
    labelled = set()
    remaining = set(saved) - old - ROOTS.keys()
    for word in remaining:
        for perm, swap in SYMMETRIES:
            labelled.add(tuple((1, 0, 2)[int(word[i])] if swap else int(word[i]) for i in perm))
    print("remaining classes", len(remaining), "labelings", len(labelled))


if __name__ == "__main__":
    main()
