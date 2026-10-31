"""Orbital scheme (integer p^k_ij) of a transitive action given by two generator images.
Usage: python build_scheme.py <action.txt> <out_prefix>
Writes <out_prefix>_intersection_numbers.csv and <out_prefix>.json (rank, valencies, pairing, commutativity).
Uses read_action, orbitals and intersection_numbers from tits1755/coset_schemes.py."""
import os, sys, json
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'tits1755'))
from coset_schemes import read_action, orbitals, intersection_numbers

def main(inp, outp):
    (gord, mord, n), gens = read_action(inp)
    assert all(len(g) == n for g in gens) and gord // mord == n
    L, val, reps, nsch = orbitals(gens, n)
    r = len(val)
    p = intersection_numbers(L, reps)
    pair = [int(L[y, 0]) for y in reps]
    for i in range(r):
        for j in range(r):
            assert p[0, i, j] == (val[i] if j == pair[i] else 0)
        for k in range(r):
            assert p[k, i, :].sum() == val[i]
    comm = bool(all(np.array_equal(p[:, i, j], p[:, j, i]) for i in range(r) for j in range(r)))
    with open(outp + '_intersection_numbers.csv', 'w') as f:
        f.write('k,i,j,p_k_ij\n')
        for k in range(r):
            for i in range(r):
                for j in range(r):
                    if p[k, i, j]:
                        f.write(f'{k},{i},{j},{p[k, i, j]}\n')
    info = dict(input=inp, n=n, group_order=gord, subgroup_order=mord, rank=r, valencies=[int(v) for v in val],
                pairing=pair, self_paired=sum(1 for i in range(r) if pair[i] == i), commutative=comm)
    json.dump(info, open(outp + '.json', 'w'), indent=1)
    print(f"n={n} rank={r} self_paired={info['self_paired']} commutative={comm}")

if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
