# The Ree-Tits family GO(q, q²): what carries over from q = 2 to q = 8

The octagon O(8) has 1 210 323 465 points, so it cannot be built and searched the way O(2) was. What follows are
statements for the whole family that are computed exactly as functions of s = q, or at (s, t) = (8, 64), and
checked against the q = 2 data.

## Flags (`hecke_family.py`)

For q = 2^(2e+1) the group 2F4(q) has a BN-pair with Weyl group D16 = I2(8). Its orbitals on flags are therefore
exactly the 16 Weyl distances, and the flag algebra is the Iwahori-Hecke algebra H(I2(8); t, s) with t = s². This
holds for every group between 2F4(q) and the full automorphism group, since all of them preserve Weyl distance
(known; standard BN-pair theory).

The script computes the seven irreducible representations of H(I2(8); s², s), checking the quadratic and braid
relations symbolically, and the generic degree D_χ = P_W χ(1) / Σ_w χ(T_w)χ(T_w⁻¹)/ind(w) as a function of s:

| H-irrep | dim (= multiplicity in C[flags]) | degree of the G-character |
|---|---:|---|
| index | 1 | 1 |
| T_σ ↦ t, T_τ ↦ −1 | 1 | s(s²−s+1)(s⁴−s²+1) |
| T_σ ↦ −1, T_τ ↦ s | 1 | s⁵(s²−s+1)(s⁴−s²+1) |
| sign | 1 | s¹² |
| 2-dim, j = 1 | 2 | [s²(s+1)²(s²−s+1)(s⁴+s³−2s²+s+1) − r s(s−1)²(s+1)³(s²−s+1)]/4 |
| 2-dim, j = 2 | 2 | s²(s²+1)²(s⁴−s²+1)/2 |
| 2-dim, j = 3 | 2 | [s²(s+1)²(s²−s+1)(s⁴+s³−2s²+s+1) + r s(s−1)²(s+1)³(s²−s+1)]/4 |

Here r = √(2s³), which is an integer exactly for s = 2^(2e+1).

- s = 2: 1, 78, 1248, 4096 and 351, 650, 675 (each twice). This is the 2F4(2) flag decomposition, and it matches
  `flags/out` and the GAP 4.15.1 run in `gap_runs/2026-11-22`, where T has constituents
  1 + 78 + 2·351 + 624 + 624' + 2·650 + 2·675 + 2048 + 2048'.
- s = 8: 1, 1839048, 7532740608, 68719476736 and 201301200, 545261600, 461921616 (each twice). The sum
  with multiplicities is 78671025225, the number of flags.
- s = 32: in `out/hecke_family_output.txt`.

So the flag scheme of O(q) has rank 16 and is H(D16; q², q) for every q. The rank-18 refinement found at q = 2
comes from T = 2F4(2)', which has index 2 and splits the two linear characters of degree s⁵(…) = 1248 and
s¹² = 4096 into 624 + 624' and 2048 + 2048'. For q > 2 the group 2F4(q) is simple, so this mechanism is gone.
Whether some other flag-transitive subgroup of Aut(O(8)) refines the scheme I haven't checked. I haven't compared
the degree formulas with the unipotent character tables in the literature (Carter's book lists them); the only
check here is q = 2.

## Triples (`triple_family.py`)

For every type (U, V, W) the script sets up the exact LP (marginals + Krein equation q⁴₄₄ = 0, which holds at both
(2,4) and (8,64)) and adds the constraints proved for every GO(s, t) in `triple/triple_relations_note.md` (G: zeros
from the geodesic-tree test with the given s, t, and the counts [313] = t, [133] = t). It then finds the implicit
zeros of the LP with nonnegativity (scaled HiGHS LP; the smallest positive scaled maximum is at least 0.12, so
the zero/non-zero split is clear) and computes the exact free dimension that is left.

Free dimension left after LP + G + implicit zeros:

| type | s = 2 | s = 8 | Ree-Tits O(2): distinct vectors |
|---|---:|---:|---:|
| (1,1,1) | 0 | 0 | 1 |
| (1,1,2) | 0 | 0 | 1 |
| (1,2,2) | 1 | 1 | 1 |
| (1,2,3) | 0 | 0 | 1 |
| (1,3,3) | 3 | 3 | 1 |
| (1,3,4) | 2 | 2 | 1 |
| (1,4,4) | 5 | 5 | 1 |
| (2,2,2) | 4 | 4 | 1 |
| (2,2,3) | 4 | 5 | 1 |
| (2,2,4) | 2 | 2 | 1 |
| (2,3,3) | 10 | 10 | 4 |
| (2,3,4) | 10 | 11 | 2 |
| (2,4,4) | 9 | 9 | 4 |
| (3,3,3) | 18 | 19 | 6 |
| (3,3,4) | 14 | 14 | 14 |
| (3,4,4) | 12 | 12 | 18 |
| (4,4,4) | 10 | 10 | 44 |

- Types (1,1,1), (1,1,2), (1,2,3) are determined at both s = 2 and s = 8 by general constraints. So in every
  GO(2,4), and in every GO(8,64), the triple numbers of these types take a single value. At s = 2 that value is
  the Ree-Tits one.
- For the other seven types where O(2) has one vector, (1,2,2), (1,3,3), (1,3,4), (1,4,4), (2,2,2), (2,2,3),
  (2,2,4), general constraints leave 1 to 5 free dimensions at both s. The rigidity of O(2) on these types is not
  explained by the constraints proved so far, at either s.
- The three differences between s = 2 and s = 8 ((2,2,3), (2,3,4), (3,3,3)) come from the degree check of the tree
  test. At s = 2 a line has only 3 points, which rules out a few configurations that exist for s = 8.
- Nothing here says which vectors O(8) actually has. That would need the orbits of a point stabiliser of 2F4(8)
  on pairs of points, which is out of reach with the explicit permutation approach used for
  q = 2.

## Triples: exact tree counts (`tree_counts.py`, `triple_tree_family.py`)

A stronger general constraint is to count, exactly and for every GO(s, t), the points x that hang off the base
configuration as a tree.

- **Base configuration.** Take the union B of the geodesics between the non-opposite pairs of u, v, w. With at most
  one opposite pair, B is a tripod whose legs are fixed by the distances. The one exception is an opposite pair,
  where the shared leg can vary; each choice is a separate base shape.
- **Lemma.** Let i, j, k ≤ 3 and suppose 2(i_a + i_b) + diam(B) < 16 for every pair a, b. Then every x with these
  distances hangs off B. Every way of hanging a path off B in the free biregular tree that ends at the right tree
  distances is realised by exactly one x. A circuit in the union would be shorter than 16, a shortcut would close
  a circuit shorter than 12, and geodesics are unique.
- **Consequence.** [ijk] equals an explicit polynomial in s and t that depends only on the base shape. At s = 2 these
  counts agree with every Ree-Tits orbit vector whose base is a tripod.
- **Not covered.** Types with two opposite pairs have no tree base and are skipped. For all-non-opposite types,
  the apartment-shaped and longer-circuit bases are not covered either.

Free dimension after LP + G + tree counts + implicit zeros, per tripod shape:

| type | s = 2 | s = 8 |
|---|---:|---:|
| (1,1,1), (1,1,2), (1,2,3) | 0 | 0 |
| (1,2,2) | 1 | 1 |
| (1,3,3) | 3 | 3 |
| (1,3,4) | 2 | 2 |
| (2,2,2) | 3 | 3 |
| (2,2,3) | 4 | 5 |
| (2,2,4) | 2 | 2 |
| (2,3,3) tripod | 7 | 7 |
| (2,3,4), both shapes | 10 | 11 |
| (3,3,3) tripod | 6 | 6 |
| (3,3,4), all three shapes | 14 | 14 |

- Where the free dimension is 0, the unique solution at s = 2 equals the Ree-Tits vector.
- Compared with LP + G alone, the tree counts lower the free dimension for (2,2,2) (4 → 3), the (2,3,3) tripod
  (10 → 7) and the (3,3,3) tripod (18 → 6). They change nothing for the types with an opposite pair.
- So the remaining rigidity of O(2) on the small types, for example (1,2,2) with 1 free dimension and (2,2,2)
  with 3, comes from points that close a circuit with the base: apartments, or longer circuits. These counts are
  not given by the tree lemma, and an argument for them is still missing at every q.
- For the opposite-pair types the tree counts do not tell the base shapes apart in the LP. Every Ree-Tits vector
  of type (2,3,4) or (3,3,4) satisfies the counts of every shape of its type, so "orbits of this shape" in the
  output means orbits consistent with that shape's counts.

Item 7 therefore ends here. The flag side is settled for the whole family. On the triple side, the general
constraints fix (1,1,1), (1,1,2) and (1,2,3) for every q. The rest needs circuit counting, which I haven't done.

## Files

- `hecke_family.py` → `out/hecke_family_output.txt`
- `triple_family.py` → `out/triple_family_output.txt`, `out/triple_family.csv`
- `tree_counts.py` (library), `triple_tree_family.py` → `out/triple_tree_family_output.txt`, `out/triple_tree_family.csv`
