"""Check the marked symmetry on the entire forbidden-walk set, before search."""
from collections import Counter
from functools import reduce
from operator import xor

import elem_setup as S


def pullback(mask, images):
    return reduce(xor, (col for j, col in enumerate(images) if mask >> j & 1), 0)


for fp in (0, 1):
    setup = S.setup_multi2([fp])
    group = S.pair_stabiliser(fp)
    supports, _ = S.sym_arrays(setup, group)
    order = list(setup['order'])
    marked = set(S.PAIR_CONS[fp])
    row_iter = iter(supports)
    shifts = Counter()
    old_valid = 0
    for g in group:
        h = S.ginv(g)
        images = [S.chord_coords(sum(1 << S.eimg(h, e) for e in range(40) if cycle >> e & 1))
                  for cycle in S.FUND]
        image = {pullback(walk, images) for walk in marked}
        shift = S.mate_shift(fp, h)
        assert {v ^ shift for v in image} == marked
        # Check the actual search columns, in the search's own coordinate order.
        if g != tuple(range(32)):
            # The nonidentity rows appear in group order.
            row = next(row_iter)
            for pos, col in enumerate(order):
                mask = (shift << 1) | 1 if col == 0 else images[col-1] << 1
                value = sum(1 << j for j, o in enumerate(order) if mask >> o & 1)
                assert int(row[pos]) == value
        shifts[shift] += 1
        old_valid += image == marked
    print(f'case {fp}: stabilizers {len(group)}, old walk-set symmetries {old_valid}, '
          f'corrected {sum(shifts.values())}, nonzero shifts {len(group)-shifts[0]}')
