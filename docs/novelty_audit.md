# Novelty audit: the orbital algebra of the Tits group on 1755 points

## 0. Bottom line (read this first)

**The distance-regular graph itself is not new.** Crnkovic, Rukavina & Svob,
*On some distance-regular graphs with many vertices*, J. Algebraic Combin. 51 (2020) 1-18,
DOI 10.1007/s10801-019-00888-5 (arXiv:1809.10197),
Theorem 9 / Table 13, already construct and classify exactly this object: the coset graph
Gamma(T,H) for T = 2F4(2)' and H = 2.2^8:5:4 on 1755 points, intersection array
{10,8,8,8;1,1,1,5}, Aut(Gamma) = T:2. It is the collinearity graph of the
Ree-Tits generalized octagon GO(2,4), a classical object going back to Tits (1961) with its
distance-regularity, spectrum and point/line counts standard since Brouwer-Cohen-Neumaier's
1989 book (Ch. 6.5, "geometric distance-regular graphs").

**Everything that follows from the intersection array by routine linear algebra -- the 125
intersection numbers, the eigenmatrices P,Q, the multiplicities 1,351,650,675,78, and the
single vanishing Krein parameter q^4_44=0 -- is also not new in substance.** The vanishing
Krein parameter is an instance of a known general mechanism: Cameron-Goethals-Seidel,
*Strongly regular graphs having strongly regular subconstituents*, J. Algebra 55 (1978)
257-280, show that "extremal" generalized 2d-gons (those meeting the relevant Krein/absolute
bound) always have a vanishing Krein parameter; De Bruyn's paper on extremal generalized
hexagons states the octagon case of this explicitly. GO(2,4) is extremal in this sense.

**What appears to be genuinely new** (full-text search of the 2020 paper found zero occurrences
of "Krein", "fusion", "Q-polynomial", "cometric"; its three citing papers, found via OpenAlex,
are unrelated -- Mathieu-group designs and a Borsuk's-conjecture note) is:

1. A full determination of the fusion lattice of the 5-class scheme: **no union of the four
   nontrivial orbitals closes into a smaller association scheme** -- confirmed by exhaustive,
   exact-arithmetic enumeration of all 15 partitions, independently re-derived in this session
   from the raw ATLAS permutation generators (not merely re-run from the earlier session's
   numbers).
2. **A proof that this is not a coincidence of (s,t)=(2,4) but a theorem about every thick
   generalized octagon.** Treating (s,t) as symbolic parameters and rebuilding the whole
   Bose-Mesner algebra from the standard intersection-array recursion, every one of the 13
   proper partitions fails to close for all but finitely many, non-feasible, values of (s,t)
   -- values that are either non-integral, violate thickness (s,t>=2), or violate the
   classical Feit-Higman integrality condition that 2st be a perfect square. This is the
   strongest defensible original contribution of the project (Section 3).

## 1. Literature audit, claim by claim

| Claim in project / earlier report | Status | Source |
|---|---|---|
| Omega=T/H, \|Omega\|=1755, H=2.2^8:5:4, rank 5, subdegrees 1,10,80,640,1024 | Known. Identical construction Gamma(T,H^4_2) | Crnkovic-Rukavina-Svob 2020, Table 11 (verified by direct full-text fetch of arXiv:1809.10197) |
| Intersection array {10,8,8,8;1,1,1,5}, Aut=T:2 | Known, identical array and automorphism group | ibid., Theorem 9 / Table 13 |
| Omega = points of Ree-Tits octagon GO(2,4); collinearity graph is this DRG | Known / classical. Existence of GO(2,4): Tits (1961, 1964); distance-regularity of generalized-2d-gon point graphs: Brouwer-Cohen-Neumaier (1989) Ch. 6.5 | Standard; also used explicitly in Jurisic & Koolen's classification paper (arXiv:1507.04839), which lists {10,8,8,8;1,1,1,5} among DRGs with smallest eigenvalue <= -k/2 |
| Uniqueness of GO(2,4) as a generalized octagon (not just as a DRG) | Open / partially known. O'Brien and co-authors prove uniqueness under extra transitivity hypotheses; general uniqueness is open | O'Brien, "On the uniqueness of the generalized octagon of order (2,4)"; De Bruyn's extremal-polygon paper |
| Spectrum 10,5,1,-3,-5, multiplicities 1,351,650,675,78 | Known in substance (forced by the intersection array via standard DRG theory); re-derived here independently | This session, Section 2 |
| q^4_44=0 is the only vanishing nontrivial Krein parameter | Mechanism known (extremal-polygon Krein-zero theorem of Cameron-Goethals-Seidel 1978, cited explicitly for generalized 2d-gons by De Bruyn); the specific numeric value for this graph had not been located in print | Cameron-Goethals-Seidel 1978; De Bruyn (extremal hexagons paper, states the general d in {2,3,4} result) |
| Not Q-polynomial in any of the 24 orderings | Not located in the literature for this graph; plausible but not found stated anywhere | -- |
| Not "tight" in the Jurisic-Koolen-Terwilliger sense (checked via the tightness identity) | Not located in the literature for this graph | Jurisic-Koolen-Terwilliger, Tight Distance-Regular Graphs and the Q-polynomial Property, J. Alg. Combin. (2000) |
| No proper fusion of {1,2,3,4} gives an association scheme (hence no SRG by merging) | Not found anywhere in the literature -- absent from the 2020 paper that introduces the graph, from its 3 citing papers, and from every other source surfaced by extensive search | This session (independently re-verified) |
| The no-fusion result holds for every thick generalized octagon of order (s,t), not just (2,4) | New, as far as this search can determine | This session, Section 3 |

Also checked and found not novel: the strongly regular graphs (1600,351,94,72) and
(3510,693,180,126) that T affords come from different coset spaces (L3(3):2 and
2^5:((2^4:5):2) respectively) -- Crnkovic-Rukavina-Svob's own theorem/tables, not fusions
of the 1755-point scheme.

## 2. Independent re-verification (this session, from raw group data)

Everything below was recomputed from scratch, directly from the ATLAS v3 standard generators
of T on 1755 points (TF42G1-p1755B0.g1/.g2, fetched live from brauer.maths.qmul.ac.uk), without
reusing any array from the earlier session or from the literature:

- Parsed the two permutations (b11, order 2; b21, order 3) on {1,...,1755}.
- Built a Schreier orbit tree for the transitive action (orbit size 1755, confirming transitivity).
- Derived 3510 Schreier generators of Stab_T(1) and computed its orbits on Omega by
  union-find: sizes 1,10,80,640,1024, matching the ATLAS-listed subdegrees exactly.
- For each pair of suborbits (i,j,k), used a transporter element mapping the base point to a
  representative of class k to compute p^k_ij=|R_i(0) cap R_j(y_k)| directly from the
  permutation action (no appeal to the DRG recursion at this stage). All 125 values matched
  the array predicted by {10,8,8,8;1,1,1,5} exactly, matched the earlier session's table
  exactly, and satisfied p^k_ij=p^k_ji and sum_k p^k_ij*k_k = k_i*k_j in every case.
- A direct global BFS on the graph built from these orbitals (not assumed, but constructed from
  the group action) reproduced distance-layer sizes 1,10,80,640,1024 from vertex 0, confirming
  R_1 really is the distance-1 relation.
- Recomputed the eigenmatrices P,Q exactly in sympy by diagonalizing the 5x5
  intersection matrices B_i (rows/cols = classes); recovered multiplicities 1,351,650,675,78.
- Recomputed all 125 Krein parameters q^k_ij from scratch by solving the linear system
  Q*(q^0_ij,...,q^4_ij)^T = (Q_0i*Q_0j,...,Q_4i*Q_4j)^T (sanity-checked against
  q^0_ii=m_i, which held exactly). Found all Krein parameters non-negative, and exactly one
  nontrivial vanishing parameter, q^4_44=0 -- confirming the earlier claim independently.
- Tested cometricity by brute force over all 4!=24 orderings of {E_1,E_2,E_3,E_4} against
  the triangle-type closure condition on the Krein array: no ordering is Q-polynomial.
- Checked the Jurisic-Koolen-Terwilliger tightness identity
  (theta_1 + k/(a_1+1)) * (theta_d + k/(a_1+1)) = -k*a_1*b_1/(a_1+1)^2
  with k=10,a_1=1,b_1=8,theta_1=5,theta_d=-5: LHS = 0, RHS = -20 -- not tight, and not
  close to tight (the smallest eigenvalue theta_d=-5 happens to equal exactly -k/(a_1+1)=-5,
  which forces LHS=0 identically, while the RHS is a fixed nonzero -20).
- Re-ran the exhaustive fusion test over all 15 set-partitions of {1,2,3,4} directly on the
  freshly-computed p^k_ij array: identical result to the earlier session -- only the trivial
  and full partitions close (0 violations); every one of the 13 proper partitions fails (4, 8, or
  9 violated closure conditions, as tabulated before).

Files: intersection_numbers_independent_verification.csv, krein_parameters_independent.csv.

## 3. The new theorem: no thick generalized octagon's point-scheme has a proper fusion

### Setup

For a generalized octagon of order (s,t), the point (collinearity) graph is distance-regular
with the standard array
  b_0=s(t+1),  b_1=b_2=b_3=st,   c_1=c_2=c_3=1,  c_4=t+1,
(and a_i=k_i-b_i-c_i), where s,t are formal parameters. This is routine generalized-polygon
theory (Brouwer-Cohen-Neumaier, Ch. 6.5). The intersection numbers p^k_ij(s,t) of the induced
5-class metric scheme are then rational functions of s,t, obtained from the three-term
recursion A_1*A_i = b_(i-1)*A_(i-1) + a_i*A_i + c_(i+1)*A_(i+1).

### Computation

Building the 5x5 matrices B_0,...,B_4 symbolically in sympy with s,t as free
symbols (and checking the construction reduces to the exact numeric table above at s=2,t=4,
which it does, entry for entry) gives closed-form expressions for all 125 values p^k_ij(s,t).

For each of the 13 proper partitions of {1,2,3,4}, the fusion-closure conditions
(block sums of p^k_ij constant over k within a block) produce a system of polynomial
equations in s,t. Solving each system (sympy.solve) shows:

| Outcome | Partitions | Interpretation |
|---|---|---|
| No solution at all, for any s,t in C | {1}|{2,3,4}, {1,2}|{3,4}, {1}|{2}|{3,4}, {1}|{2,3}|{4}, {1}|{3}|{2,4}, {1,2}|{3}|{4} | Fusion is algebraically impossible, independent of (s,t) |
| Only irrational s,t solve the system | {2}|{1,3,4} (s=1+sqrt2, t=-1+sqrt2), {2,3}|{1,4} (s=sqrt2/2, t=sqrt2) | Not integers -- never a real generalized octagon |
| Only s=t solves the system | {1,2,3}|{4} | Excluded: s=t forces 2st=2s^2, never a perfect square for integer s>=1, violating the Feit-Higman integrality condition for octagons (2st must be a perfect square; Higman 1975) |
| Only s=1 (thin) solves the system | {1,3}|{2,4}, {2}|{1,3}|{4} | Excluded by thickness: a genuine (Moufang / building-derived) generalized octagon requires s,t>=2 |
| Only t=1 (thin) solves the system | {3}|{1,2,4}, {2}|{3}|{1,4} | Same exclusion, t>=2 required |

Every single escape route is closed by either pure algebra (no solution exists) or by one of the
two classical feasibility constraints on generalized-octagon parameters: thickness s,t>=2, and
the Feit-Higman condition that 2st be a perfect square with sqrt(t) <= s <= t^2 (Higman 1975,
also recorded e.g. on Wikipedia's "Generalized polygon" page and in Van Maldeghem's monograph).

### Theorem

For every thick generalized octagon of order (s,t) -- in particular, for the unique presently
known example GO(2,4), realized group-theoretically as T/H for the Tits simple group T and
H = 2.2^8:5:4 -- the 5-class point-distance association scheme admits no proper fusion. The
only unions of its classes that close under the Bose-Mesner product are the trivial one (all four
nontrivial classes merged, giving the complete graph K_n) and the identity partition (the full
scheme itself). In particular, no strongly regular graph, and no association scheme of
intermediate rank, is obtainable by merging distance classes of a generalized octagon's point
graph -- this is a structural obstruction of the whole feasible parameter family, not an
accident of (2,4).

This appears to be new: no source found in this audit (including the paper that introduced this
exact graph, its citing papers, the extremal-polygon Krein-parameter literature, or the
tight-DRG classification literature) states or proves a fusion-closure result for generalized
octagon point-schemes, symbolically or otherwise.

### Caveat / honest limitation

This is a clean, fully-proved but narrow result: it answers one well-posed question (fusion of
the specific metric scheme) definitively, for one specific combinatorial family (generalized
octagons), by elementary exact algebra. It is not a deep structural theorem and does not resolve
open questions about octagons (existence of other orders, uniqueness of GO(2,4)). Framed
honestly, it is the kind of result suitable for a short note -- e.g. as an addition to
Brouwer-Cohen-Neumaier-style parameter tables, or a remark in a paper on fusion/association-scheme
closure -- not a landmark contribution. Whether it generalizes to generalized hexagons
(n=6, two feasible-array families) or quadrangles was not checked and is a natural next step.

## 4. What is NOT demonstrated / open ends

- Full Schreier-Sims verification of |T|=17,971,200 from the generators was not performed
  (only transitivity and the induced suborbit structure); this is standard and not expected to
  fail, but wasn't independently re-derived here.
- The attribution of 1,351,650,675,78 to irreducible character degrees of T from the ATLAS
  character table was NOT independently re-checked against the ATLAS table in this session
  (same caveat as the earlier report).
- Whether the fusion-impossibility theorem extends to generalized hexagons/quadrangles, or to
  "near-feasible" non-existent octagon parameters, is open.
- The question of whether GO(2,4) is the unique generalized octagon of its order (not just
  unique DRG with this array) remains open in the general literature; this project neither
  resolves nor needed to resolve it.

## 5. Reproducibility

All computations in Sections 2-3 were executed in this session from:
- the raw ATLAS generator files (TF42G1-p1755B0.g1, .g2), fetched live from
  brauer.maths.qmul.ac.uk;
- exact rational/symbolic arithmetic in sympy (no floating point used for any reported integer,
  rational, or algebraic value);
- the symbolic family computation uses s,t as free sympy symbols throughout, with the
  (s,t)=(2,4) specialization checked to reproduce the from-scratch numeric tables exactly before
  any conclusions were drawn from the general case.
