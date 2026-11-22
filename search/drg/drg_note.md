# Distance-regular graphs with a primitive automorphism group (v < 4096)

## Question

Brouwer's machine-readable versions of the BCN tables list the feasible intersection arrays of distance-regular
graphs on at most 4096 vertices. Which of the primitive arrays (not bipartite, not antipodal) with v < 4096 and
diameter at least 3 are realised by a graph whose automorphism group is primitive on vertices?

## What is known

- Distance-transitive graphs with a primitive group are nearly classified (see the introduction of
  'Imprimitive distance-transitive graphs with primitive core of diameter at least 3', Michigan Math. J. 58 (2009)). That covers rank d + 1 only; a vertex-primitive graph can have a
  group of much larger rank.
- GRAPE's `VertexTransitiveDRGs` finds the distance-regular graphs for one given group. Crnkovic, Rukavina and Svob search the actions of single groups.
- I did not find a run over the whole primitive groups library against the BCN arrays. Two keyword searches only.

## Method

1. `fetch_arrays.py`: 260 arrays in IA3P, IA4P, ..., IA13; 250 are neither bipartite nor antipodal, on 157
   values of v. Four of them have v = 4096 and are out of range of the library.
2. `export_suborbits.g` (GAP 4.16.1, PrimGrp): suborbit lengths of all 4244 primitive groups on those 156 degrees
   other than Sym(v) and Alt(v) (`suborbits.txt.gz`).
3. `drg_filter.py`: an invariant graph is a union of orbitals and so is each distance graph, so the group needs
   rank >= d + 1 and its nontrivial suborbits must split into parts of sizes k_1, ..., k_d. Of the 12173
   (group, array) pairs 3102 fail on rank and 6882 on the split; 2189 pass. Arrays of H(d,q) with q != 4 are
   then dropped, since those graphs are determined by their array. That leaves 806 pairs, 663 groups and 117
   arrays (`survivors.csv`, `drg_filter_output.txt`).
4. `drg_search.py` (uses `search/dfs_fixed.py` in DRG mode): for 661 of the 663 groups, all symmetric unions S of
   orbitals with valency b_0, a_1 common neighbours on S and 0 or c_2 off S, each candidate then checked exactly
   for distance-regularity and its array. Integer arithmetic, no node limit reached, 74 s in total.
5. `drg_affine.py` for the two affine groups of rank 94 and 63 on 1024 points: a graph invariant under
   2^10:G_0 is a Cayley graph with integer eigenvalues, and the spectrum is fixed by the array, so every
   character sum must be one of the eigenvalues with the right multiplicities. CP-SAT on that model; a solution
   is built as a Cayley graph and tested directly. One case, {231,160,6;1,48,210} on 1024/1, hit the time limit;
   that array is realised by other groups, so it doesn't change anything below.

`drg_summary.py` puts the results together (`drg_summary.csv`, one line per array).

## Result

Of the 117 arrays, 86 are realised by a graph with a primitive group (47 of diameter 3, 25 of diameter 4, 9 of
diameter 5, 4 of diameter 6, 1 of diameter 7). All 86 are known graphs: 76 are in the DistanceRegular.org
catalogue, and the other ten are J(10,4), the thin octagon GO(4,1), the thin hexagons GH(7,1), GH(8,1),
GH(9,1), the hexagon GH(2,8), H(5,4), the thin dodecagon of order (3,1), GO(2,4) and GO(4,2). So no new graph.

For the other 31 arrays no primitive group of degree v carries a graph with that array. So if a graph with
one of these arrays exists, its automorphism group is not primitive on vertices: it is either intransitive, or
transitive with a system of imprimitivity. Both are possible for distance-regular graphs in general; a DRG
does not have to be vertex-transitive (the twisted Grassmann graphs of van Dam and Koolen have the parameters
of Grassmann graphs and are not vertex-transitive, Fujisaki-Koolen-Tagami, Innov. Incidence Geom. 3 (2006)
81-87). Status of the 31 (`array_status.csv` has the sources):

- does not exist (11): {13,10,7;1,2,7}, {22,16,5;1,2,20}, {36,28,4;1,2,24}, {55,36,11;1,4,45},
  {104,66,8;1,12,88}, {119,96,18;1,16,102} (all in the nonexistence lists of Vidali's sage-drg),
  {25,24,3;1,3,20} (Zhurtov-Gerieva 2025), {42,30,12;1,6,28} (Belousov-Makhnev 2018) and {20,16,5;1,1,16}
  (Makhnev-Golubyatnikov, Diskret. Mat. 34:1 (2022) 76-87, English version Discrete Math. Appl. 33:5 (2023)
  273-281; the paper treats the family {mn-1,(m-1)(n+1),n-m+1;1,1,(m-1)(n+1)} and also rules out
  {39,36,4;1,1,36} and {90,84,7;1,1,84}, which are not in this census. A search for the exact array misses it
  because the array only appears inside the paper), {44,36,5;1,9,40} (Chen-Makhnev-Klimin, Proc. Inst. Math.
  NAS Belarus 32:1 (2024) 57-63; abstract read) and {49,40,22;1,5,28} (Gavrilyuk-Vidali, arXiv:2611.02986,
  Thm 2.6(a), posted 2 Nov 2026);
- does not exist if Koval's October 2026 preprint holds (1): {39,30,4;1,5,36};
- exists (3), each with an imprimitive automorphism group, as the search says it must be:
  - {33,30,15;1,2,15}, the Shi-Krotov-Sole graph (Des. Codes Cryptogr. 87 (2019)); its group has order 55296
    according to DistanceRegular.org, and no primitive group of degree 1024 has that order;
  - {26,24,19;1,3,8}, the Brouwer graph Br(3) on F_3^3 x F_3^3, (a,b) ~ (c,d) iff a != c and b - d = a x c;
  - {31,30,17;1,2,15}, the coset graph of the binary cyclic code of length 31 with zeros alpha and alpha^3
    (zeros alpha and alpha^5 give the same array).
  The last two are built in `known_graphs.py`, which checks the array exactly at every vertex, and
  `aut_check.g` (GRAPE) gives |Aut| = 8188128 with blocks of size 27 for Br(3) and |Aut| = 158720 with
  blocks of size 32 for the coset graphs (`known_graphs_output.txt`);
- no status found (16): {18,14,5;1,2,14}, {27,16,9;1,3,12}, {48,35,9;1,7,40},
  {42,32,9;1,6,24}, {78,50,9;1,15,60}, {39,32,16;1,4,24}, {32,28,9;1,2,28}, {52,42,16;1,6,28},
  {56,42,20;1,6,28}, {80,63,11;1,9,70}, {63,50,8;1,5,56}, {21,20,10;1,1,12},
  {40,33,19;1,3,22}, {54,42,11;1,3,44}, {63,48,24;1,6,28}, {99,80,12;1,10,88}.

For those 16 I know of no graph and no nonexistence proof. Some have partial results: automorphisms of a
graph with {48,35,9;1,7,40} were studied by Makhnev, Bitkina and Gutnova (Vladikavkaz Mat. Zh. 22:2 (2020)),
and arXiv:2206.00031 shows that several of these arrays can't come from a coset graph of a linear completely
regular code. The new restriction here is that none of the 16 can have a vertex-primitive group. The three
existing graphs above, which the search correctly misses, are the only test of that statement I have.

Status check of 22 Nov 2026. The 18 arrays of the first version were 'no status found' after keyword searches;
two of them turned out to be settled: {44,36,5;1,9,40} by Chen, Makhnev and Klimin (2024; the earlier keyword
searches missed it) and {49,40,22;1,5,28} by Gavrilyuk and Vidali
(2026). For the other 16 a Crossref and OpenAlex search on each array string found only the automorphism
papers on {48,35,9;1,7,40} mentioned above and a paper that lists {18,14,5;1,2,14} without a result
(Efimov-Makhnev, Ural Math. J. 2020). The four antipodal diameter-4 arrays that Cesarz and Williford rule
out (arXiv:2610.05420) are outside this census, which only has primitive arrays, so they don't overlap.

## What rests on what

- The library of primitive groups of degree < 4096 (Coutts, Quick, Roney-Dougal) and GAP's suborbit lengths.
- The DFS and the exact distance-regularity check; `search/dfs_fixed.py` was checked on the T-actions
  (`search/out/dfs_fixed_validation.txt`).
- For the two groups on 1024 points, CP-SAT's infeasibility answers, with no independent certificate.
- The statuses come from the sources named, not from a review. Two of the three existing graphs were
  found because a pasted list pointed to them; the catalogue I checked first doesn't have them.

## Summary

The search over the primitive groups library gives no new distance-regular graph, as expected, and shows that
16 arrays with no known graph and no known nonexistence proof can't have a vertex-primitive realisation. That is
a table for a note, not a paper by itself.
