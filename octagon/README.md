# Opposite components and GO(2,4) uniqueness

[Current status: sibling obstruction, symmetry correction and open cases](status.md).

Fix a point x of a generalized octagon of order (2,4), let B be the 640 points at distance 3 and O the 1024 points
opposite x, and let the outer graph be collinearity on O. Each component is a cover of the Clebsch graph (even
five-bit words, adjacency at Hamming distance 4) through the trace words, and the sheet numbers add up to 64.

## 1. One-component recognition, checked again

`verify_recognition.py` rebuilds O(2) from the ATLAS generators of T on 1755 points (TF42G1-p1755B0.g1/.g2;
set ATLAS_DIR) and runs `recog.py`, a separate implementation of the recognition steps that only sees the uncoloured
edge list of a component. For both 512-vertex components:

* every edge has exactly one other edge at outer distance δ ≥ 6, the forced pairing agrees with the true boundary
  points, and all six inward merges (640 → 320 → 80 → 40 → 10 → 5 → 1, thresholds 14, 12, 10, 8, 6, 4 on incidence
  distance through the partial half H) are forced and agree with the true ancestors;
* the reconstructed half contains a full suboctagon of order (2,1) (found by geodesic closure from random 4-point
  seeds containing x), so De Bruyn's theorem (Discrete Math. 338 (2015) 2125–2142: the Ree–Tits octagon is the only
  GO(2,4) with a full suboctagon of order (2,1)) applies.

The two components are isomorphic (BLISS); K has girth 8, diameter 6 and |Aut K| = 10240.

On the step from outer distance to incidence distance: the mate test only uses δ, which is an outer-graph distance;
the inward merges use distances in the partial incidence graph H (points of K, its 1280 lines, the 640 boundary
points). H is a subgraph of the incidence graph of any completion, so H-distances are upper bounds for true ones,
and "same parent ⇒ H-distance ≥ 16 − L" is what the girth argument gives. The direction of every test is the one
needed. A path in H between two boundary points never uses an inner element, so it closes up with the tree path to a
circuit.

## 2. The C5-invariant class, reproduced

`c5_reproduce.py` does the count with the cycle space inside F_2^40 and the rotation acting on edges (no spanning-tree
gauge): 2,167,056 invariant covers, 251 with girth ≥ 8, mate histogram {(0,0): 217, (0,1): 28, (1,1): 5, (5,5): 1}.
Of the five forced-pair covers, four have no admissible sibling at any of the 640 boundary points and the fifth is K.
The 5-mate cover has 320 eligible pairs and no sibling-compatible pair of pairs in each of the ten categories.

## 3. New necessary conditions

**Boundary count (any GO(2,4), any component C).** Put f = 1_C − (|C|/1024) 1_O. Every point of O has its five
neighbours in O inside its own component and its other neighbours outside O, so f^T A f = 5 f^T f, and f is
orthogonal to the constants. Since 5 is the largest eigenvalue of A off the constants, Af = 5f. At a boundary point b
this reads Σ f over its 8 outer neighbours = 0, and the two outer points of each of its 4 outward lines lie in one
component, so b lies on exactly |C|/256 lines of C. Hence |C| ∈ {256, 512, 768, 1024}, sheet numbers are 16, 32, 48 or
64, and in a 512-point component every boundary point carries exactly two of its edges, mate graph triangles or not.
In O(2) the count is 2 + 2 at all 640 boundary points.

**Least eigenvalue (any GO(2,4)).** With e_y the normalised vectors of the 78-dimensional −5 eigenspace
(Gram matrix (−1/2)^d, frame operator (1755/78) I), the Gram matrix of the 1024 vectors at distance 4 has spectrum
4, 12^45, 15^32 (checked in O(2); the general argument goes through S_0 − S_1/2 + S_2/4 − S_3/8 + S_4/16 = 0). So a
vector supported on O keeps at most 15/22.5 = 2/3 of its norm in the −5 eigenspace, the rest of the spectrum is ≥ −3,
and every outer component has least eigenvalue ≥ −13/3. In particular no component is bipartite. In O(2) the least
eigenvalue of each component is −3.9624.

For any component with m sheets the same counting gives, in every category, Σ 1/(min(deg e, 3) + 1) ≤ 64 over its
4m lifted edges e, where deg is the degree in the mate graph. For m > 16 this forces mates in every category, and for
a single 64-sheeted component (256 lifted edges per category) every lifted edge has at least three mates and every
boundary point carries four edges of the component.

Let K be a component with 32 sheets (512 vertices, 1280 edges). Each of the 10 categories (coordinate, bit) has 128
lifted edges and 64 boundary points, and each boundary point lies on 4 outer lines.

* Edges sharing a boundary point are pairwise at δ ≥ 6 and in one category, so the K-edges of a category split into at
  most 64 cliques of the mate graph, each of size ≤ 4. Hence every category has at least 64 mate edges, and if
  mc(ε) is the number of mates of a lift of the base edge ε, then Σ_ε 1/(min(mc(ε),3)+1) ≤ 2 over the four base edges
  of the category.
* If the mate graph of a category is triangle-free, every one of its 64 boundary points carries exactly two K-edges,
  and their pairs form a perfect matching of the mate graph.
* Two K-edges whose boundary points have lowest common ancestor b, N, p2, M in the inward tree satisfy δ ≥ 6, 5, 4, 3
  respectively (circuit of length 2δ + 4 + L ≥ 16). In particular sibling boundary points need δ ≥ 5 across.
* For a regular cover with abelian deck group the mate graph is invariant under the deck group, so mates come in
  orbits (ε, ε′, t). An orbit with ε′ ≠ ε and t = 0 never occurs (some walk of length ≤ 5 already closes it up).

## 4. Elementary abelian 32-sheeted covers without the five-cycle

Here K is assumed to be a connected regular cover with deck group F_2^5, and nothing is assumed about lifted
automorphisms. The cover is a voltage map on the 25 chords of a BFS spanning tree, taken up to GL(5,2).

Reduction. Every category has a mate orbit, and Aut(Clebsch) (order 1920, coordinate permutations and even
translations) is transitive on edges, so some mate orbit can be moved onto base edge 4. Its stabiliser (order 48) has
one orbit on the other three edges of that category, so there are two cases: a same-edge orbit (4, 4, t) and a cross
orbit (4, 17, t). The translate t is nonzero and becomes the first basis vector.

Search (`elem_setup.py`, `elem_fast.py`): columns are t and the 25 chord voltages in a greedy order; GL(5,2) is broken
by first-occurrence echelon form, the stabiliser of the forced orbit (48 and 32 elements) by a lex-leader test that
compares canonical forms of the images as soon as their columns are determined. A node is cut when a simple base cycle
of length 4–7 gets voltage 0, when the forced orbit stops being at δ ≥ 6, or when the clique bound of §3 fails in some
category (the bound only uses constraints whose walks are fully assigned). Leaves are tested twice more:

* stage 2 (inside the search): in every triangle-free category take the forced pairs (mate components of size 2),
  build the partial H from them, and require each forced boundary point to have a possible sibling at H-distance ≥ 14
  (other forced points, or a still undetermined mate pair, using 1 + the smaller distance to its two lines);
* stage 3 (`stage3.py`): for the remaining leaves, a category in which no two disjoint mate pairs are at δ ≥ 5 across
  is impossible (siblings carry two mate pairs each, by the boundary count).

Controls: K is found (as the lex-minimal image of its voltage map) and passes all three stages; in a sample of 2000
leaves 1997 failed stage 2; the lex-leader pruning was checked against unpruned subtrees.

Historical run at 2d1ca19. The cross row used a fixed-translate symmetry which fails for 28 of its 32 group elements. The saved leaf tests reproduce, but that row does not establish completeness.

| case | prefixes at depth 8 | nodes | leaves rejected at stage 2 | disconnected | left for stage 3 | rejected at stage 3 | survivors |
|---|---:|---:|---:|---:|---:|---:|---:|
| cross orbit (4,17,t) | 9914 | 1,720,076,943 | 4776 | 18 | 152 | 152 | 0 |
| same-edge orbit (4,4,t) | 5925 | 10,779,630,523 | 5,646,997 | 8 | 912 | 866 | 46 |

In the same-edge case one of the 46 survivors passed stage 2 outright and is isomorphic to K (BLISS). The other 45 are
labelled covers falling into 7 isomorphism classes of uncoloured graphs (girth 8, diameter 7 or 8, automorphism
groups of orders 128, 512, 2048, 2048, 3072, 4096, 6144). None of them is K. They pass stage 3 in every category, and
`stage4.py` (unit propagation of the boundary pairing, then the H-sibling test) retains them. Their voltages
are in `stage4_fp0.json`. The subsequent uncoloured sibling test in `sibling_obstruction.py` excludes all 45.
Only row 6 remains and is isomorphic to both known components. The corrected cross search has 15,086
depth-eight prefixes. The archived cross results do not exhaust the corrected search. See [status.md](status.md).


Stage 3 is now run on every category, not only triangle-free ones. This is valid because, by the boundary count,
every boundary point carries exactly two K-edges.

Without the stages, girth alone leaves of the order of 10^11–10^12 covers (Knuth estimates), so the mate and sibling
conditions do the work.

## 5. What is still open

* Component sizes are 256, 512, 768 or 1024 (§3), so the patterns are 1024, 768 + 256, 512 + 512, 512 + 256 + 256 and
  4 × 256. `regular256/` excludes every 256-point component that is a regular cover (all deck groups of order 16),
  through the spectral window [−4, −2] ∪ [−1, 1 + √6] for nonprincipal eigenvalues of such a component and a tower of
  double covers. A second end-to-end implementation now checks all stages, including the 168 order-four records
  and the terminal bipartiteness tests. See [regular256/README.md](regular256/README.md). The [double-tower extension](regular256/PROOF.md) also excludes every 256-point cover with 2-group monodromy.
  The other 256-point components, 768- and 1024-point components, and 512-point components
  that are not regular elementary abelian covers are all still open.
* The Gram representation gives the eigenvalue bound of §3 but I have no case where it rejects a 512-point cover that
  the distance tests keep. For two mates at one boundary point all four endpoint distances are 2, which already makes
  the two boundary vectors equal.
* A proof that the axioms force a 512-point component, or a thin suboctagon directly.
