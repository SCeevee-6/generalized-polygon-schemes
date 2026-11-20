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

- 13 of the 18 arrays fail the Krein test, so none of them is a Cayley graph on an abelian group:
  {18,14,5;1,2,14}, {44,36,5;1,9,40}, {27,16,9;1,3,12}, {48,35,9;1,7,40}, {42,32,9;1,6,24}, {78,50,9;1,15,60},
  {39,32,16;1,4,24}, {32,28,9;1,2,28}, {80,63,11;1,9,70}, {63,50,8;1,5,56}, {21,20,10;1,1,12}, {40,33,19;1,3,22},
  {99,80,12;1,10,88}.
- Four of the other five are ruled out by one quotient of order 125. Every abelian group of order 625 has a quotient
  of order 125, and every abelian group of order 750 or 1000 has its Sylow 5-subgroup as a quotient. For each array
  the quotient model is infeasible for all three groups Z_5^3, Z_25 x Z_5 and Z_125 (fibre v/125), with both CP-SAT
  and SCIP:
  - {52,42,16;1,6,28}, v = 625 (fibre 5);
  - {49,40,22;1,5,28}, v = 750 (fibre 6);
  - {54,42,11;1,3,44} and {63,48,24;1,6,28}, v = 1000 (fibre 8).

  So none of these four is a Cayley graph on an abelian group. For {52,42,16;1,6,28} on Z_5^4 this gives a second
  proof of the Shi-Krotov-Sole result. The direct atom model agrees on the 29 (group, array) pairs where it finishes.
- {56,42,20;1,6,28} on 729 vertices: the atom model rules out Z_729, Z_243 x Z_3, Z_81 x Z_9, Z_81 x Z_3^2, Z_27^2,
  Z_27 x Z_9 x Z_3 and Z_27 x Z_3^3 (both solvers). Z_3^6 is the Shi-Krotov-Sole case. That leaves Z_9^3,
  Z_9^2 x Z_3^2 and Z_9 x Z_3^4, where the atom model doesn't finish. The quotient model for Z_9 x Z_9 x Z_3 with
  fibre 3 is infeasible (both solvers), and that group is a quotient of Z_9^3 and of Z_9^2 x Z_3^2, so those two are
  out as well. Z_9 x Z_3^4 is still open: its quotients of order 243 are Z_9 x Z_3^3 and Z_3^5, and for the first
  neither solver finished in 90 minutes. All quotients of order at most 81 that I tried have solutions.

## What rests on what

- The Krein lemma and the atom lemma above.
- CP-SAT and SCIP infeasibility answers. They agree everywhere I ran both, but there is no proof log.
- For the quotient argument, only the identities in step 3, which are short to check by hand.

## Summary

Of the 18 arrays, 17 are not Cayley graphs on any abelian group (13 by the Krein test, four by a quotient of order
125). For the last one, {56,42,20;1,6,28}, only Z_9 x Z_3^4 is left. Together with the primitive-group census, a
graph with one of the 17 arrays would need an automorphism group that is not vertex-primitive and has no regular
abelian subgroup. No new graph.
