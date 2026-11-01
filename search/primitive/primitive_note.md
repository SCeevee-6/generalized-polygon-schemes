# SRGs on open parameter sets from primitive groups (v <= 1300)

Labels: **known**, **new** (not found in the sources I checked), **incomplete** (computation not finished).

## Question

Brouwer's SRG table (fetched 1 Nov 2026) has 1367 open parameter sets with v <= 1300 (counting a set and its
complement once), on 430 distinct values of v (`open_srg.csv`). Is any of them realised by a graph with a
primitive group of automorphisms?

## What is already known

- Rank 3 groups give rank 3 graphs, all classified and all in the table, so only rank >= 4 matters.
- Brouwer–Van Maldeghem, *Strongly regular graphs* (2022), Table 11.9 lists SRGs whose automorphism group is
  primitive of rank at most 10 (as cited in arXiv:2606.10652, which redoes degree 120 with GRAPE's
  `VertexTransitiveDRGs`). I have not seen the table itself.
- Crnković, Švob and coauthors searched transitive actions of many specific groups (U(3,3), U(3,4), U4(3), L(4,3),
  O±(6,2), M11, L(2,q) for q <= 23, L(3,q) for q <= 5, U(5,2), the Mathieu groups, G2(4), G2(5), O(7,3), the Tits
  group, A5 on 81 points; OpenAlex list of their papers, 1 Nov 2026). Those are group-by-group, not degree-by-degree.
- I found no systematic run over the whole primitive groups library restricted to open degrees.

## Method

1. GAP 4.16.1, primitive groups library: for each of the 430 degrees, every primitive group except Sym(v) and
   Alt(v) with rank >= 4 (`export_primitive.g`). That gives 2164 groups (`prim_ranks.csv`); most degrees have only
   Sym and Alt, and most groups are affine (degrees 289, 343, 361, 529, 625, 729, 841, 961, 1024) or product
   actions (1296, 1000).
2. For each group: orbitals and exact intersection numbers from the generators (`tits1755/coset_schemes.py`), then
   the exact depth-first search of `search/dfs_numba.py` over unions of transposition-closed units
   (`search/prim_search.py`). This finds every SRG and every DRG that is a union of orbitals.
3. Groups where the search is too slow (rank above about 30): for each open parameter set of that degree, a CP-SAT
   model with k, lambda, mu fixed, where the SRG conditions are linear in the unit variables and their products
   (`search/srg_fixed_cpsat.py`). INFEASIBLE is a proof that no union has those parameters.

Controls: on degree 120 the search gives SRG(120,56,28,24) from 8 groups, plus (120,51,18,24), (120,42,8,18) and
T(16), in line with arXiv:2606.10652. The fixed-parameter model finds srg(121,40,15,12) and srg(1024,495,238,240)
on the groups where the DFS found them.

## Results so far

- DFS complete on 2038 of the 2164 groups (every group of rank <= 30 and some above; `dfs_results.jsonl.gz`).
  Every SRG found has parameters that Brouwer's table lists as existing (309 parameter sets), and none lies on an
  open set. **new** (search range), empty answer.
- The other 126 groups (rank 31 to 281) are running with the fixed-parameter model. After the first 36 groups:
  178 (group, parameter set) pairs fail at once because k is not a sum of unit valencies, 138 are INFEASIBLE, and
  8 are UNKNOWN after 15 minutes (degrees 861, 903, 946 and 1024). No solution. **incomplete**.

## Assessment

If the rest also comes out empty, the statement is: no open parameter set with v <= 1300 is realised by a graph
admitting a primitive group of rank >= 4 from the GAP library. That is a negative result. It is useful for the
table and might be a short section of a note, but it is not a discovery.
