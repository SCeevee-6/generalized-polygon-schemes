# Radius-one packings in the dual GO(4,2)

`dual_packing_133.txt` lists 133 lines of O(2) (zero-based, lines ordered by their sorted point triples in the
labelling of generators.npy, the same convention as the earlier 129-centre list). Their pairwise distances in the
dual are at least 3 and their 133 closed balls cover 1729 = 133·13 of the 2925 dual points. It was found by an
iterated local search with (1,2)-swaps (`mis_ils.py`, best of 12 runs; the other runs gave 127–130). A HiGHS MILP
on the full model stopped at 112 after 25 minutes.

Upper bound: the Delsarte LP on the distance distribution gives 225, exactly the perfect-code value, so it adds
nothing to the bound 223 that comes from the nonexistence of a perfect packing and the defect-one argument. The
current interval is 133 ≤ max ≤ 223.
