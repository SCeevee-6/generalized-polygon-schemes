# Strongly regular and distance-regular graphs with a transitive action of T = 2F4(2)'

Labels: **known** (with reference), **new** (not found in the sources I checked; the check is a search, not a
full review), **pending** (computation not finished).

## Starting point

Crnković, Rukavina and Švob (J. Algebraic Combin. 51 (2020), arXiv:1809.10197, Section 3.4, Theorems 8 and 9)
determined all SRGs and DRGs of diameter at least 3 on which T acts primitively with at most 12000 vertices,
and imprimitively with at most 7000 vertices. Their answer: SRG(1600,351,94,72) (primitive), SRG(3510,693,180,126)
(imprimitive, Aut = Fi22.2), and the octagons GO(2,4) on 1755 and GO(4,2) on 2925 vertices. **known.**
The three papers that cite theirs (OpenAlex, 31 Oct 2026) do not extend the range for T.

## Actions

From the table of marks of 2F4(2)' in TomLib (434 classes of subgroups; 10 classes of index at most 7000,
which matches the count in CRS 2020). Index in (7000, 20000]:

| TomLib class | index | subgroup order | rank | self-paired | units | how searched |
|---|---|---|---|---|---|---|
| 424 | 7020 | 2560 | 16 | 10 | 12 | Gray code, all 2^12 unions |
| 423 | 8775 | 2048 | 18 | 12 | 14 | Gray code |
| 422 | 11700 | 1536 | 30 | 16 | 22 | Gray code |
| 420, 421 (maximal A6.2^2) | 12480 | 1440 | 24 | 16 | 19 | Gray code (ATLAS max6, max7) |
| 419 (maximal 5^2:4A4) | 14976 | 1200 | 25 | 17 | 20 | Gray code (ATLAS max8) |
| 418 | 17550 | 1024 | 33 | 19 | 25 | Gray code |
| 416 | 17550 | 1024 | 45 | 21 | 32 | CP-SAT |
| 417 | 17550 | 1024 | 54 | 22 | 37 | depth-first search (`dfs_numba.py`), complete |

"Units" are the transposition-closed sets of non-trivial orbitals (a self-paired orbital, or a pair {R, R^T}).
Ranks from the computed orbitals agree with the double-coset counts in GAP. All of these actions are
non-commutative. The coset actions are in `tits1755/cosets_tom/` (GAP: `tits1755/gap/tom_coset_actions.g`,
generators of `UnderlyingGroup(TableOfMarks("2F4(2)'"))`); the maximal ones come from the ATLAS SLPs as before.

## Method

`build_scheme.py` computes the orbitals and integer p^k_ij (from `tits1755/coset_schemes.py`).
`srg_drg_search.py`:

- `gray`: runs through all unions S of units with a Gray code, keeps those with a_1 constant
  (sum_{i,j in S} p^k_ij the same for all k in S), and for those builds the distance partition from the
  structure constants and checks distance-regularity layer by layer. Exact integers.
- `cuts`: the same question as a CP-SAT model (OR-Tools), solved repeatedly with no-good cuts until infeasible.
  Its constraints are necessary conditions for a connected DRG (constant a_1 on S; walk counts off S in {0, c_2};
  S not inside a proper closed subset, i.e. the graph is connected). Every solution is then checked with the exact
  distance-regularity test. The model was checked against the Gray-code search on four schemes (7020, 8775,
  11700, 14976): same solution sets (`out/validate_cuts_output.txt`).
- `dfs_numba.py` (used for class 417, where CP-SAT ran 3 h without an answer): the same necessary conditions
  (constant a_1 on S; walk counts off S in {0, c_2}, c_2 occurring), decided unit by unit in a depth-first
  search. The walk counts c_k are monotone in S, so after each decision the possible values of a_1 and c_2 are
  intervals [c_k(chosen), c_k(chosen + undecided)], and a branch is cut when one becomes empty. Every surviving
  union goes through the exact distance-regularity test. Checked against the earlier searches on classes 424, 422,
  418 (Gray code) and 416 (CP-SAT): the same SRGs and no others (`out/dfs4*.json`).

Controls: on the commutative actions the search finds exactly the CRS graphs: SRG(1600,351,94,72) and its
complement on both 1600-point actions, GO(2,4) on 1755, GO(4,2) on 2925, nothing on 2304 (`out/prim_max*.json`).

## Results

- Degrees 12480 (both classes) and 14976: no SRG and no DRG of diameter >= 3. Together with CRS 2020 this
  covers every primitive action of T, since 14976 is the largest index of a maximal subgroup. **new** (the search
  range; the answer is empty).
- Imprimitive degrees 7020, 8775, 11700, 17550 (classes 416, 418): the only SRGs are complete multipartite
  graphs (mu = k, for instance (7020,7016,7012,7016), parts of size 4), which come from the block systems, and
  there is no DRG of diameter >= 3. **new** (range), empty answer.
- Class 417 (17550 points, rank 54, 37 units): 1.1e8 search nodes, about 6.5 minutes. 14 unions pass the necessary
  conditions: the four complete multipartite graphs coming from the block systems (blocks of size 2, 5, 6, 10,
  i.e. over the 8775 flags, 3510 points, 2925 lines and 1755 points of GO(2,4)) and ten disconnected unions
  inside a block system. So no other SRG and no DRG of diameter >= 3. A second run with a different unit order
  (`out/dfs417_desc2.json`, a different search tree) gives the same 14 unions. **new** (range), empty answer.

So, up to 17550 vertices, the only non-trivial SRGs and DRGs of diameter at
least 3 with a transitive action of T are the four found by Crnković, Rukavina and Švob. Every action of
2F4(2) in which T is transitive is one of these T-actions, so this also covers 2F4(2) on up to 17550 points.

## Assessment

This finishes the primitive case and pushes the imprimitive range from 7000 to 17550, which is now complete.
It produced no new graph. As a stand-alone paper it is thin; it could be a section of a note together with
something else.
