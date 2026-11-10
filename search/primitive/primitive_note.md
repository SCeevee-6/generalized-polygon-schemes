# SRGs on open parameter sets from primitive groups (v <= 1300)

## Question

Brouwer's SRG table (fetched 1 Nov 2026) has 1367 open parameter sets with v <= 1300 (counting a set and its
complement once), on 430 distinct values of v (`open_srg.csv`). Is any of them realised by a graph with a
primitive group of automorphisms?

## What is already known

- Rank 3 groups give rank 3 graphs, all classified and all in the table, so only rank >= 4 matters.
- Brouwer–Van Maldeghem, *Strongly regular graphs* (2022), Table 11.9 lists SRGs whose automorphism group is
  primitive of rank at most 10 (as cited in arXiv:2606.10652, which redoes degree 120 with GRAPE's
  `VertexTransitiveDRGs`). The table itself was not seen.
- Crnković, Švob and coauthors searched transitive actions of many specific groups (U(3,3), U(3,4), U4(3), L(4,3),
  O±(6,2), M11, L(2,q) for q <= 23, L(3,q) for q <= 5, U(5,2), the Mathieu groups, G2(4), G2(5), O(7,3), the Tits
  group, A5 on 81 points; OpenAlex list of their papers, 1 Nov 2026). Those are group-by-group, not degree-by-degree.
- A systematic run over the whole primitive groups library restricted to open degrees was not found.

## Method

1. GAP 4.16.1, primitive groups library (complete below degree 4096): for each of the 430 degrees, every primitive
   group except Sym(v) and Alt(v) with rank >= 4 (`export_primitive.g`). That gives 2164 groups (`prim_ranks.csv`).
   Most degrees have only Sym and Alt; most of the groups are affine (degrees 121, ..., 961, 343, 625, 729, 1024)
   or product actions (1000, 1296).
2. Exact DFS (`search/prim_search.py`, using `search/dfs_numba.py`): orbitals and exact p^k_ij from the generators,
   then all unions of transposition-closed units with constant a_1 and walk counts off S in {0, c_2}, each checked
   exactly for distance-regularity. This finds every SRG and every DRG that is a union of orbitals. Used for all
   2038 groups of rank <= 30 and five more (861/1, 903/1, 946/1, 1128/1, 1176/1).
3. Degree p^2 (p = 11, 13, 17, 19, 23, 29, 31): no computation needed. GAP confirms that all 460 groups of these
   degrees in the list have an elementary abelian regular socle N = Z_p^2, so every invariant graph is a Cayley graph
   Cay(N, D). For an SRG on p^2 vertices both restricted eigenvalues are integers (a conference graph has
   eigenvalues (-1 +- p)/2 here), so chi(D) is a rational integer for every character chi of N. For t in F_p^* the
   Galois automorphism zeta_p -> zeta_p^t gives chi(tD) = chi(D) for all chi, hence tD = D. So D is a union of
   punctured lines through 0 of AG(2,p), and the graph has Latin square parameters (p^2, m(p-1), p-2+(m-1)(m-2),
   m(m-1)), which exist for every m. None of them is open. (Standard multiplier argument; see Ma's survey of
   partial difference sets.) As a check, all 116114 SRGs the DFS found on these degrees are of this type.
4. Groups of rank above 30 that are not of degree p^2: first the fixed-parameter CP-SAT model
   (`search/srg_fixed_cpsat.py`; exact, INFEASIBLE is a proof), which settled 14 groups. For the other 14 two
   floating-point MILPs were used at first (`search/srg_spectral_milp.py` for commutative configurations,
   `search/srg_hybrid_milp.py` for 1015/1 and 1225/13). Their answers are kept in `milp_results.jsonl` and
   `hybrid_results.jsonl`, but the result no longer depends on them; see 5.
5. Exact redo of those 14 groups.
   - The twelve affine ones (625/109, 625/530, 729/296, 729/330, 729/331, 729/337, 729/346, 729/347, 1024/1,
     1024/2, 1024/3, 1024/6): `affine_exact.py`. The group is N:G_0 with N elementary abelian of order p^n, an
     invariant graph is Cay(N, D) with D a union of G_0-orbits, and its eigenvalue at the character chi_y is
     sum_c N_c(y) zeta^c with N_c(y) = #{x in D : y.x = c}. That is a rational integer iff N_1 = ... = N_{p-1}, and
     then it is N_0 - N_1. So "SRG with eigenvalues r, s" is a system of linear equations with integer coefficients
     in the 0/1 indicators of the units, which CP-SAT decides exactly. The coordinates come from GAP
     (`export_affine.g`, files in `affine/`). All open parameter sets: infeasible, or k is not a sum of unit sizes.
     As a control the model finds srg(1024,495,238,240) on 1024/1, and the graph it returns was rebuilt as a
     Cayley graph and checked directly (k = 495, lambda = 238, mu = 240).
   - 1015/1 (socle PSL(2,29), rank 52, 43 units) and 1225/13 (PSL(2,49), rank 39, 32 units), both non-commutative:
     `search/dfs_fixed.py`, a depth-first search with k, lambda, mu fixed. The numbers sum_{a,b in S} p^j_ab and
     the valency are monotone in S, so each partial choice gives intervals and a branch is cut when a target
     leaves its interval. Integer arithmetic only. No solution for any open parameter set (1225/13: 46 sets,
     94 minutes; 1015/1: 6 sets, 97 minutes, the longest one 1.9e9 nodes). The search was checked on the T-action
     on 1600 points, where it finds srg(1600,351,94,72) and its complement, and on TomLib class 424
     (`search/out/dfs_fixed_validation.txt`).
6. A parameter set whose k is not a sum of unit valencies is ruled out for that group at once.

`aggregate.py` combines everything into `summary.csv` (one line per group: method, settled, hits).

Controls: on degree 120 the DFS gives SRG(120,56,28,24) from 8 groups, plus (120,51,18,24), (120,42,8,18) and
T(16), in line with arXiv:2606.10652. The MILPs find srg(121,40,15,12), srg(1024,495,238,240) and T(48) on
1128/1 where the DFS finds them.

## Result

All 2164 groups are settled (`summary.csv`). Every SRG the searches produce has parameters that Brouwer's table
already lists as existing (309 parameter sets from the DFS), and no open parameter set with v <= 1300 is
realised by an SRG admitting a primitive group of automorphisms. Rank 3 groups and 2-transitive groups are left
out because their graphs are known (rank 3) or trivial.

Status: not found in the sources checked (BVM22 Table 11.9 covers primitive groups
of rank <= 10 only, as cited by arXiv:2606.10652). Every group is now settled without floating point
(`summary.csv`, column `exact`): DFS 2043, multiplier argument for degree p^2 93 more, CP-SAT with walk counts 14,
CP-SAT with character sums 12, fixed-parameter DFS 2. The CP-SAT answers rely on the solver's infeasibility proofs;
no independent certificate was produced.

## Summary

None of the open SRG parameter sets in Brouwer's table up to 1300 vertices can be filled by a graph with a primitive automorphism group. That is a negative result, worth a line in the table or a section in a note rather than a paper. Finding an actual graph would need imprimitive groups or groups that are not transitive at all. The transitive groups library stops at degree 47 (48 with the extra package), below the smallest open v in the table (88), so the imprimitive case cannot be done from a library.
