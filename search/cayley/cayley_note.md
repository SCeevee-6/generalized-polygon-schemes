# The 18 open arrays as Cayley graphs on abelian groups

## Question

`search/drg/drg_note.md` leaves 18 primitive arrays of diameter 3 (v < 4096) with no known graph, no
nonexistence proof that I could find, and no vertex-primitive realisation. Can any of them be a Cayley graph on an
abelian group?

## What is known

- If a distance-regular graph is a Cayley graph on an abelian group, its distance scheme is a translation scheme.
  The dual scheme lives on the character group and its intersection numbers are the Krein parameters, so these must
  be nonnegative integers. This is standard (BCN 2.10, and it is used in the recent paper on Fourier analysis on
  distance-regular Cayley graphs over abelian groups, Electron. J. Combin. 33(1) (2026) #P1.14).
- An integral Cayley graph on an abelian group has a connection set that is a union of atoms
  [g] = {x : <x> = <g>} (Bridges and Mena).
- Shi, Krotov and Sole (arXiv:2206.00031, Discrete Math. 2022) show that a distance-regular Cayley graph of diameter
  at least 3 on an elementary abelian p-group is the coset graph of a linear completely regular code, and they rule
  that out for {52,42,16;1,6,28} (q = 5) and {56,42,20;1,6,28} (q = 3), among others.
- Schur-ring results (Schur, Wielandt, Bercov; see Wielandt's book, Sect. 25, and Ponomarenko-Ryabov,
  arXiv:2410.22998): over many abelian groups every primitive S-ring is trivial, and then every primitive
  distance-regular Cayley graph on the group is complete. Bercov's theorem is quoted in the result below.
- Distance-regular Cayley graphs on Z_p^s + Z_p, on cyclic, dihedral and dicyclic groups are classified in other
  papers. I found no treatment of these 18 arrays on general abelian groups. Keyword searches only.

## Method

1. `krein_filter.py`: exact Krein parameters of the 18 arrays (`krein_filter.csv`).
2. `abelian_cayley.py`: for each surviving array and each abelian group G of order v, a 0/1 model over the atoms of
   G. Each atom goes to one class S_1, S_2, S_3, each nontrivial character atom to one eigenspace E_1, E_2, E_3, the
   character sums of the S_i must be the eigenvalue matrix entries P_ji, and the dual sums over the E_j must be the
   entries Q_ij. Solved with CP-SAT; every solution is checked directly for distance-regularity. `abelian_cayley_mip.py`
   is the same model in SCIP.
3. `quotient_lift.py`: for N <= G write w_i(g + N) = |S_i cap (g + N)|. Characters of G/N are characters of G that
   are trivial on N, so
   - sum over g in G/N of chi(g) w_i(g) = P_ji for chi in E_j, and
   - |N| sum over chi in E_j cap N^perp of chi(g) = sum_i Q_ij w_i(g) + Q_0j [g in N],
   with 0 <= w_i <= |N|, sum_i w_i(g) = |N| (minus 1 on N itself), and w_i constant on atoms of G/N. This model only
   depends on G/N and |N|. If it has no solution, no abelian group with a quotient G/N of that size carries the
   graph. `quotient_mip.py` is the SCIP version.
4. `quotient_table.py` runs the quotient checks below with both solvers (`quotient_table.csv`). The runs of the
   atom model on all 36 (group, array) pairs other than Z_3^6 are in `atom_runs.csv` (`run_all.py`, 600 s per pair).

Checks of the code: the atom model finds the Clebsch graph on Z_2^4, the Brouwer-Haemers graph on Z_3^4, the folded
7-cube on Z_2^6 and H(3,4) on Z_4^3, and the quotient model finds the Brouwer-Haemers graph over Z_3^2 (fibre 9) and
H(3,4) over Z_2^3 (fibre 8), with the lifts checked as graphs. The Aut(G)-orbit types used for symmetry breaking in
`abelian_cayley.py --symbreak` were compared with a brute-force orbit computation on six small groups.

## Result

None of the 18 arrays is the array of a Cayley graph on an abelian group.

- 13 fail the Krein test: {18,14,5;1,2,14}, {44,36,5;1,9,40}, {27,16,9;1,3,12}, {48,35,9;1,7,40}, {42,32,9;1,6,24},
  {78,50,9;1,15,60}, {39,32,16;1,4,24}, {32,28,9;1,2,28}, {80,63,11;1,9,70}, {63,50,8;1,5,56}, {21,20,10;1,1,12},
  {40,33,19;1,3,22}, {99,80,12;1,10,88}.
- The other five are {52,42,16;1,6,28} (v = 625), {56,42,20;1,6,28} (729), {49,40,22;1,5,28} (750),
  {54,42,11;1,3,44} and {63,48,24;1,6,28} (1000). They are primitive, and their distance classes are unions of
  atoms, so the distance module is a rational primitive S-ring over G (Miklavic-Potocnik, JCTB 97 (2007),
  Prop. 3.6). Each abelian group of these orders is excluded by at least one of:
  1. Bercov (Canad. J. Math. 21 (1969), Theorem 2.2, which restates his 1965 result). If G = A x B x C with A cyclic
     of order p^a != 4, B of exponent < p^a and |C| prime to p, then a nontrivial rational primitive S-ring over G
     forces G = H_1 x ... x H_d with |H_i| = q > 2 for all i, and the union of the H_i minus the identity, which has
     d(q-1) elements, is a basis set. None of the class sizes of the five arrays has the form d(q-1) with q^d = v.
     This covers every group of these orders except Z25^2, Z5^4, Z27^2, Z9^3, Z9^2 x Z3^2, Z3^6, Z4 x Z2 x Z5^3 and
     Z2^3 x Z5^3.
  2. A short counting argument (below), for groups whose Sylow p-subgroup is Z_{p^e} x Z_p^r with e >= 2. It
     reproves most of the Bercov cases and also covers Z4 x Z2 x Z5^3.
  3. The atom bound (below): if the Sylow p-subgroup has exponent p^e with e >= 2, then
     p^(e-1)(p-1) <= max(a_1, c_2). This covers every group with a cyclic factor of order 25, 125, 625, 27 or
     more, in particular Z25^2 (20 > 9) and Z27^2 (18 > 13), which before needed the solver.
  4. Shi-Krotov-Sole for Z5^4 and Z3^6 (coset graphs of linear codes; their result is computational).
  5. Computation: the quotient Z9 x Z9 x Z3 (fibre 3) for Z9^3 and Z9^2 x Z3^2, and the quotient Z5^3 (fibre 8)
     for Z2^3 x Z5^3 = Z10^3. Both solvers give infeasible in each case. The quotient Z5^3 also gives a second
     proof for Z5^4. The atom model had also excluded Z25^2 and Z27^2.

  `atom_bound_output.txt` applies the atom bound to every abelian group of the five orders.

  `coset_lemma_output.txt` lists every group with the argument that excludes it.

In particular Z9 x Z3^4, the group that was still open for {56,42,20;1,6,28}, is excluded by Bercov's theorem (A = Z9,
B = Z3^4; d(q-1) would be 52, 24 or 12, while the class sizes are 56, 392 and 280) and by the counting argument
(24 > 13). So this case was never really open; I only missed the Schur-ring literature.

## A short argument

Let Cay(G, S) be distance-regular on an abelian group G with integral eigenvalues k = theta_0 > theta_1 > ..., and
let the Sylow p-subgroup of G be Z_{p^e} x Z_p^r with e >= 2. Let A be the subgroup of index p containing every
element whose p-part has order at most p, let K be the subgroup of order p inside p^{e-1}G, and let lambda be a
character with kernel A.

- S is a union of atoms. If s is not in A, its order n is divisible by p^2, the numbers 1 + j n/p are units mod n, and
  (n/p)s generates K, so s + K lies in the atom of s. Hence for h in K, h != 0, every s in S \ A has s + h in S, and
  |S \ A| <= #{s in S : s + h in S} = p^I_{11}, where h lies in S_I.
- lambda takes every nontrivial p-th root of unity equally often on an atom outside A, so
  lambda(S) = |S cap A| - |S \ A|/(p-1). This is an eigenvalue theta_J, and so |S \ A| = (p-1)(k - theta_J)/p.

So (p-1)(k - theta_1)/p <= max(a_1, c_2). For {56,42,20;1,6,28} and p = 3 this says 24 <= 13. The bound is
attained by the 3-cube on Z4 x Z2, the 4-cube, the Clebsch graph and the folded 7-cube on Z4 x Z2^r
(`coset_lemma.py` checks the identity and the inequality on graphs found by `abelian_cayley.py`, and the two
group facts by brute force).

## The atom bound

This is simpler than the argument above and works when the largest cyclic factor occurs more than once.
Let s in S have order n with p^2 | n, and put h = (n/p)s. For every unit u mod n, u + n/p is again a unit,
because every prime dividing n divides n/p. So the atom [s] is mapped to itself by translation by h. Since
[s] is inside S, every element of [s] is a common neighbour of 0 and h, and h != 0, so

  phi(n) <= |S cap (S - h)| <= max(a_1, c_2).

If the Sylow p-subgroup of G has exponent p^e with e >= 2, the elements whose p-part has order less than p^e form
a proper subgroup, so a connected S contains such an element and phi(p^e) <= max(a_1, c_2). The eigenvalues of
the five arrays are integers (`atom_bound.py` checks this), so S is a union of atoms.

| array | group | forced order | phi | max(a_1, c_2) |
|---|---|---:|---:|---:|
| {52,42,16;1,6,28} | Z25^2 | 25 | 20 | 9 |
| {56,42,20;1,6,28} | Z27^2 | 27 | 18 | 13 |

The general form (an integral distance-regular Cayley graph on an abelian p-group of exponent p^e, e >= 2, has
p^(e-1)(p-1) <= max(a_1, c_2)) is elementary and very likely known in some form; I haven't found it stated.
For {56,42,20;1,6,28} on Z9^3 and Z9^2 x Z3^2 it gives only 6 <= 13.

A short proof for Z9^3 doesn't come out of this either. On Z9^3 the order-9 atoms in S correspond to points of
PG(2,3) (the atom of s goes to <3s>) and the index-3 subgroups to lines. The bound above says at most two
atoms over a point whose order-3 subgroup lies in S, at most one if it is at distance 2, none at distance 3.
Each line must miss 4, 6 or 7 of the atoms (the eigenvalues 20, 2, -7), and MacWilliams on the subgroup 3G ties
the class counts on 3G to the line counts. `occupancy_z93.py` finds 468 configurations that satisfy all of this
(all with nine order-9 atoms), so the finer quotient model is really needed there.

## Nonabelian groups: normal Cayley graphs

Here S is a union of conjugacy classes. Then every distance class is too, and each eigenspace E_k is a sum of the
central idempotents e_phi over a set of irreducible characters phi. If you compare the coefficient of phi in
E_i o E_j, you get

  q^k_ij phi(1) = sum over chi in E_i, psi in E_j of chi(1) psi(1) <chi psi, phi>,

an integer, for every phi in E_k. So every q^k_ij is rational, and phi(1) is a multiple of D_k, the lcm of the
denominators of the q^k_ij. Since phi(1) divides v and the phi(1)^2 over E_k sum to m_k, `normal_krein_filter.py`
excludes all 13 arrays that fail the abelian Krein test, on every group of order v (`normal_krein_filter.csv`):

- {18,14,5;1,2,14} and {40,33,19;1,3,22} have irrational Krein parameters;
- for ten more arrays, D_1 or D_3 has a prime factor not dividing v (for example D_1 = 32 for v = 280);
- {78,50,9;1,15,60}: D_1 = 7, so the irreducibles in E_1 have degree 7 or 14, but m_1 = 26 < 49.

If phi is linear the condition says q^k_ij is an integer, which is the abelian Krein test again. I didn't find this
filter in the literature, but I only did one keyword search. Integrality of normal Cayley graphs (S power-closed) is
Godsil-Spiga.

That leaves the five arrays that pass the abelian test. For them, `export_normal.g` exports the rational character
data of the nonabelian groups of orders 625, 729, 750 and 1000. It skips groups with a normal cyclic Sylow subgroup,
because every nontrivial central S-ring over such a group is imprimitive (Chen-Muzychuk-Ponomarenko,
arXiv:1505.01352, Thm 1.3). `normal_cayley.py` is then the atom model with rational classes and Galois orbits of
irreducibles in place of atoms. It finds the Clebsch graph, the Brouwer-Haemers graph, the 3-cube, H(3,4) and the
4-cube on abelian groups, K_{3,3} on S_3 and K_{4x2} on Q_8. On Q_8 it finds no 3-cube, which is correct.

Result of the run over all groups of these orders (GAP 4.16.1 SmallGroups library and character tables;
`out_normal/tables_*.txt`, `out_normal/results.jsonl`, summary in `out_normal/normal_summary.csv`):

| v | array | groups | abelian | skipped (normal cyclic Sylow) | modelled | candidates |
|---:|---|---:|---:|---:|---:|---:|
| 625 | {52,42,16;1,6,28} | 15 | 5 | 0 | 10 | 0 |
| 729 | {56,42,20;1,6,28} | 504 | 11 | 0 | 493 | 0 |
| 750 | {49,40,22;1,5,28} | 39 | 3 | 29 | 7 | 0 |
| 1000 | {54,42,11;1,3,44} | 199 | 9 | 13 | 177 | 0 |
| 1000 | {63,48,24;1,6,28} | 199 | 9 | 13 | 177 | 0 |

Every model is infeasible (CP-SAT, 8 workers; the longest took 187 s, no time-outs), so there is no candidate
to check in GAP. Together with the Krein-degree filter and the abelian result: none of the 18 arrays is the
array of a normal Cayley graph on any group of order v. This says nothing about Cayley graphs whose
connection set is not a union of conjugacy classes, where a new graph is still possible.

## What rests on what

- The Krein lemma and the atom lemma above.
- CP-SAT and SCIP infeasibility answers. They agree everywhere I ran both, but there is no proof log.
- For the quotient argument, only the identities in step 3, which are short to check by hand.
- Bercov's Theorem 2.2 as stated in the 1969 paper; I did not read the 1965 proof.
- For the normal case: GAP's character tables of the SmallGroups of orders 625, 729, 750, 1000, CP-SAT
  infeasibility (no proof log), and Chen-Muzychuk-Ponomarenko Thm 1.3 for the 42 skipped groups.

## Summary

None of the 18 open arrays from the primitive-group census is realised by a Cayley graph on an abelian group. Most
of this follows from known results: the Krein test, Bercov's theorem on primitive S-rings and Shi-Krotov-Sole. What
is mine is the computation for Z9^3, Z9^2 x Z3^2 and Z10^3, the atom bound, which settles Z25^2 and Z27^2 by
hand, and the short counting argument, which is an elementary special case of Bercov's theorem. No normal Cayley graph on any group of order v has one of these arrays either.
Not a paper. A graph with one of these arrays would need a group that
is neither vertex-primitive nor contains a regular abelian subgroup, and if it is a Cayley graph its
connection set is not a union of conjugacy classes.
