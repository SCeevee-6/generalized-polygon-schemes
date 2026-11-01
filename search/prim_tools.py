"""Orbital scheme of a permutation group given by generator files (no OR-Tools import)."""
import os, sys, inspect
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, '..', 'tits1755'))
import coset_schemes as cs

# coset_schemes.orbitals stores labels as int8 (rank <= 127); this copy uses int16
_src = inspect.getsource(cs.orbitals).replace('dtype=np.int8', 'dtype=np.int16')
_ns = dict(cs.__dict__); exec(_src, _ns); orbitals16 = _ns['orbitals']


def read_gens(gdir, v, i):
    return [np.array([int(x) - 1 for x in l.split(',') if x], dtype=np.int64)
            for l in open(os.path.join(gdir, f'{v}_{i}.txt')).read().split('\n') if l.strip()]


def scheme(gdir, v, i):
    """p[i,j,k] = p^k_ij of the orbital configuration."""
    L, sizes, reps, _ = orbitals16(read_gens(gdir, v, i), v)
    pk = cs.intersection_numbers(L, reps)
    return np.ascontiguousarray(pk.transpose(1, 2, 0))
