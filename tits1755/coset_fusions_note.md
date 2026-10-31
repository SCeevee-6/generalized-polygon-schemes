# Fusions and algebra structure of the orbital schemes of T = 2F4(2)' on cosets of its maximal subgroups

Labels used below: **known** (in the literature, with reference), **new here** (not found in the sources I
checked; the check is a search, not a full review), **check** (a computational confirmation of something
that follows from known data).

## Method

1. `gap/coset_actions.g` (GAP 4.16.1) builds the action of T on the cosets of max i (ATLAS SLP
   TF42G1-max{i}W1 on the standard generators of TF42G1-p1755B0) with `FactorCosetAction`, and prints the
   two generator images. The outputs are in `cosets/max{i}.txt`, for i = 1, 2, 4, 5, 6, 7, 8. Max 3 is the
   1755-point action, which is treated elsewhere in the repository.
2. `coset_schemes.py` works from those two permutations only. It builds a Schreier tree from point 0, takes
   the orbits of the stabiliser from all Schreier generators (so the orbitals are exact, not sampled),
   fills the full orbital matrix, and counts the intersection numbers p^k_ij as integers. It checks
   p^0_ij = delta(j, i*) k_i and sum_j p^k_ij = k_i, tests commutativity, and computes the centre of the
   orbital algebra exactly (sympy rank of the commutator system).
3. For the commutative cases it tests every partition of the non-trivial classes for closure, with exact
   integers. The criterion is the one in `fusion/verify_fusions.py`: for all blocks A, B, C the sum of
   p^h_ij over i in A, j in B is the same for every h in C. It also records whether the partition is
   closed under transposition. The eigenmatrix, multiplicities and Krein parameters are floating point
   (the multiplicities round to the integers given below with error below 1e-6). P- and Q-polynomial
   orderings are searched for every choice of first class or first idempotent.

Run time: about 2 minutes for all seven actions. `python tits1755/coset_schemes.py tits1755/cosets/max5.txt tits1755/out_cosets/max5`
writes `out_cosets/max5.json` and `out_cosets/max5_intersection_numbers.csv`.

## Checks against the GAP results

For all seven actions the rank and the number of self-paired orbitals equal the values from the
character-table computation in `coherent_configurations_note.md` (double cosets and Frobenius-Schur sums).
For the commutative cases the multiplicities equal the degrees of the constituents of 1_M^T. **check**

| max | M | n | rank | valencies | multiplicities |
|---|---|---|---|---|---|
| 1, 2 | L3(3):2 | 1600 | 4 | 1, 312, 351, 936 | 1, 351, 624, 624' |
| 4 | L2(25) | 2304 | 7 | 1, 78, 300, 300, 325, 325, 975 | 1, 27, 27', 351, 624, 624', 650 |
| 5 | 2^2.[2^8].S3 | 2925 | 6 | 1, 12, 96, 768, 1024, 1024 | 1, 351, 624, 624', 650, 675 |

## Commutative cases: closed fusions

**Degree 1600 (both classes of L3(3):2).** 3 proper partitions tested. One closes: {1,3}|{2}. It gives a
strongly regular graph with parameters (1600, 351, 94, 72) (complement (1600, 1248, 968, 992)), eigenvalue
multiplicities 351 and 1248 = 624 + 624'. **known**: Crnković, Rukavina and Švob, J. Algebraic Combin. 51
(2020), Theorem 8 and Table 12, list this graph with full automorphism group T; it also appears among the
exceptional rank-four strongly regular graphs in arXiv:2510.25410. The parameters are those of a Latin square
graph LS_9(40); I did not check whether the graph has 40-cliques.

**Degree 2304 (L2(25)).** 201 proper partitions tested (all set partitions of the six non-trivial classes,
excluding the trivial two). Three close, all closed under transposition:

| blocks | rank | symmetric | multiplicities | origin |
|---|---|---|---|---|
| {1} {2} {3} {4,5} {6} | 6 | no | 1, 27, 27', 351, 650, 1248 | orbitals of T.2 on cosets of L2(25).2_3 |
| {1} {2,3} {4} {5} {6} | 6 | yes | 1, 54, 351, 624, 624', 650 | symmetrisation |
| {1} {2,3} {4,5} {6} | 5 | yes | 1, 54, 351, 650, 1248 | symmetrisation of the T.2 scheme |

The T.2 identification: the library permutation character of 2F4(2) on L2(25).2_3 has rank 6, with 27 and 27'
still separate and 624 + 624' fused to an irreducible of degree 1248 (`gap/extra_checks_output.txt`). The T.2
orbital scheme is a rank-6 fusion of the T scheme, and the symmetrisation is the only other rank-6 fusion, but
it has an eigenspace of dimension 54; so the T.2 scheme is the first row. No fusion is P- or Q-polynomial and
there is no rank-3 fusion, which agrees with Theorem 8 of Crnković–Rukavina–Švob (only one strongly regular
graph with a primitive T-action up to 12000 vertices). The fusion list itself is **new here** in the sense that
I found it in no source, but it is a routine computation.

**Degree 2925 (2^2.[2^8].S3, the line stabiliser of GO(2,4)).** 50 proper partitions tested. One closes:
{1}|{2}|{3}|{4,5}. The fused scheme has valencies 1, 12, 96, 768, 2048, is P-polynomial in the natural order,
and is the distance scheme of the dual octagon GO(4,2) (lines of GO(2,4); array {12,8,8,8;1,1,1,3}). It is
the orbital scheme of T.2 = 2F4(2) on 2925 points (library rank 5, constituents 1, 351, 650, 675, 1248). It is
not Q-polynomial, as De Bruyn–Vanhove (Combinatorica 35 (2015) 181–208) predicts for thick regular near
octagons. **known** as a graph (the collinearity graph of the dual Ree–Tits octagon); the statement that T itself splits the distance-4 class into two orbitals of size 1024
is a **check**.

**Q-polynomiality.** None of the commutative T schemes on 1600, 2304 or 2925 points is P- or
Q-polynomial. Among the fusions only the rank-3 graph (trivially) and the GO(4,2) scheme (P only) are.

## Non-commutative cases: structure of the orbital algebra

For a transitive group the orbitals already form a coherent configuration (one fibre), so the coherent
closure is the orbital configuration itself; there is nothing further to close. The orbital algebra is the
centraliser (Hecke) algebra End_T(1_M^T), which is semisimple over C and isomorphic to the direct sum of
matrix algebras M_m(C) over the constituents of 1_M^T, with m the multiplicity.

| max | M | n | rank | self-paired | centre dim (exact) | Wedderburn type forced by rank and centre | from the character |
|---|---|---|---|---|---|---|---|
| 6, 7 | A6.2^2 | 12480 | 24 | 16 | 13 | 11 x C + M_2(C) + M_3(C) | 351 three times, 675 twice |
| 8 | 5^2:4A4 | 14976 | 25 | 17 | 13 | 9 x C + 4 x M_2(C) | 351, 675, 1300, 1300' twice each |

The centre dimension was computed from the integer structure constants of the orbital algebra, so this is
independent of the character table. With rank r and centre dimension z, the block sizes m_1..m_z satisfy
sum m_i^2 = r; for (24, 13) and (25, 13) the only solutions are the ones above, and they agree with the
multiplicities from the library character table. **check**

So these are non-commutative association schemes in the usual sense (homogeneous coherent configurations, not
commutative), not a different kind of algebra. All constituents with multiplicity above 1 have indicator +1, so
over R the algebra is still a sum of full matrix algebras over R. Valencies:

- max 6/7: 1, 20, 45, 90, 120, 144, 180, 360 (x5), 720 (x10), 1440 (x2).
- max 8: 1, 75, 100, 150 (x2), 300, 400 (x4), 600 (x9), 1200 (x6).

Fusions of the non-commutative schemes were not searched (there are too many partitions for the brute-force
test: 23 and 24 non-trivial classes).

## Corrections to coherent_configurations_note.md

- The pair 624/624' is **not** a complex-conjugate pair. Both characters are real (indicator +1), and
  complex conjugation fixes each of them (`gap/extra_checks_output.txt`). The pairs that are complex conjugate
  are 26/26', 27/27' and two of the three 351s (library X.10, X.11).
- 27/27' is a complex-conjugate pair (checked), which explains the one pair of non-self-paired orbitals on
  2304 points.

## Novelty of the commutative cases (Breuer–Lux)

The Breuer–Lux database shipped with ctbllib 1.3.11 (`tst/multfree.dat`, variable MULTFREEINFO) has 38 groups:
the 26 sporadic groups and 12 automorphism groups. The Tits group and 2F4(2) are not among them, and the
Breuer–Müller endomorphism-ring file `tst/mferctbl.gap` does not mention 2F4 either. So the multiplicity-free
actions of T are not covered by those tables. They are not new as such: Crnković–Rukavina–Švob (2020, Section 3.4) tabulate the ranks of T on cosets of its subgroups of
index up to 7000 (the rows I could see are 1600, 1755, 2925 and 3510; I did not read the full table), and the 1755-point action is a Gelfand pair in arXiv:2510.25790.
