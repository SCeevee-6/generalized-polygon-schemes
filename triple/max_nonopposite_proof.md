# Maximum sets of pairwise non-opposite points in a generalized octagon of order (2,4)

Status: **new** as far as I can tell (I haven't found it in the EKR-for-buildings literature I've read so far, which deals with polar spaces, projective spaces and flags). Hand proof, not refereed. Every lemma below has been checked numerically on the Ree-Tits octagon (`max_set_checks.py`), and the end result for the Ree-Tits octagon has been checked by two independent clique searches (`line_cliques.py` with CP-SAT, `clique_certificate.py` with igraph/Cliquer and a separate branch and bound). The proof itself uses only the octagon axioms and the parameters (s,t) = (2,4), so it applies to **any** GO(2,4), not just the Ree-Tits one.

**Theorem.** Let S be a finite generalized octagon of order (2,4) and let C be a set of points, no two of them opposite. Then |C| ≤ 27, and |C| = 27 if and only if C is the set of points at distance at most 1 from some line L (the three points of L and the 24 points collinear with one of them). L is unique.

## Setup

Distances are in the collinearity graph, so they run from 0 to 4 and "opposite" means distance 4. The graph is distance-regular with intersection array {10,8,8,8; 1,1,1,5}, v = 1755, k_i = 1, 10, 80, 640, 1024, a_1 = a_2 = a_3 = 1. A_i is the distance-i matrix. The eigenvalues of A_1 are 10, 5, 1, −3, −5, and A_4 = P_4(A_1) takes the values

| θ | 10 | 5 | 1 | −3 | −5 |
|---|---|---|---|---|---|
| P_4(θ) | 1024 | −16 | 16 | −16 | 64 |

(exact, printed by `max_set_checks.py`). Two facts about generalized octagons are used throughout:

* (G1) A path x_0 ∼ x_1 ∼ … ∼ x_r with r ≤ 4 in which consecutive steps use different lines is a geodesic, so d(x_0, x_r) = r. (It is a non-backtracking path of length 2r ≤ 8 in the incidence graph, which has girth 16.)
* (G2) If d(x,z) ≤ 3 the geodesic from x to z is unique. For a line M on x, say z is *behind M* (from x) if that geodesic starts along M.

For x ∈ C write n_i(x) = |C ∩ Γ_i(x)| and a(x) = n_1(x). For y ∉ C write m_i(y) = |C ∩ Γ_i(y)|. For x ∈ C and a line M on x write d_M(x), c_M(x) for the number of points of C at distance 2, resp. 3, from x behind M, and for a point y ∈ M \ {x} write d_y(x), c_y(x) for those whose geodesic from x passes through y.

## Step 1: the ratio bound and what equality gives

C is a coclique of the opposition graph Γ_4, which is 1024-regular with least eigenvalue −16, so |C| ≤ 1755·16/(1024+16) = 27. If |C| = 27, equality in the ratio bound puts 1_C − (27/1755)·1 in the −16-eigenspace of A_4, which by the table is the sum of the A_1-eigenspaces for θ = 5 and θ = −3. So

* A_4 1_C = 16(1 − 1_C), and
* (A_1 − 5I)(A_1 + 3I) 1_C is a multiple of 1. Using A_1² = 10I + A_1 + A_2 this operator is A_2 − A_1 − 5I, and applying it to the constant part gives (27/1755)(80 − 10 − 5) = 1. So (A_2 − A_1 − 5I) 1_C = 1.

Reading these off at a point (together with n_0 + … + n_4 = 27):

* for x ∈ C: (n_1, n_2, n_3, n_4) = (a, a + 6, 20 − 2a, 0);
* for y ∉ C: (m_1, m_2, m_3, m_4) = (m, m + 1, 10 − 2m, 16), so in particular m_1(y) ≤ 5.

## Step 2: local lemmas

**Lemma A.** Let x ∈ C, let M be a line on x that is not contained in C, and let y ∈ M \ C. Then c_M(x) = 4 − 2a(x).

*Proof.* Take z ∈ C. If d(x,z) ≤ 2, or d(x,z) = 3 and z is behind M, then d(y,z) ≤ 3. If d(x,z) = 3 and z is not behind M, the path y, x, … , z is non-backtracking of length 4, so d(y,z) = 4 by (G1). Since C ∩ Γ_4(x) = ∅, this gives m_4(y) = n_3(x) − c_M(x), i.e. 16 = 20 − 2a(x) − c_M(x). ∎

**Lemma B.** For x ∈ C, a(x) ∈ {0, 2, 10}. If a(x) = 2, exactly one line on x lies in C and the other four lines on x meet C only in x. If a(x) = 10, all five lines on x lie in C.

*Proof.* Let f be the number of lines on x inside C and p the number that meet C in exactly two points, so a = 2f + p. If some line on x is not inside C, Lemma A gives c_M = 4 − 2a ≥ 0 for it, so a ≤ 2. Hence a ≥ 3 forces f = 5, a = 10. If a = 1, all five lines are not inside C and Σ_M c_M = 5·2 = 10 ≠ 18 = n_3. If a = 2 and f = 0, then Σ_M c_M = 0 ≠ 16 = n_3. ∎

So if no point of C has a = 10, every collinear pair in C lies on a line contained in C, and every line through a point y ∉ C carries at most one point of C.

**Lemma C.** Let x ∈ C. If two points z, z' of C are at distance 2 from x and behind different lines on x, then z and z' are opposite. Hence, when a(x) ≤ 2 and the points of C collinear with x have a ≠ 10, all points of C ∩ Γ_2(x) are behind one line on x.

*Proof.* z, y, x, y', z' is non-backtracking of length 4 (z is not on xy because d(x,z) = 2), so d(z,z') = 4 by (G1). For the second sentence: a point of C at distance 2 from x behind a line inside C would be collinear with a point p of that line off the line, giving p two lines meeting C, so a(p) = 10 by Lemma B. ∎

**Lemma D.** Let x ∈ C with a(x) = 0, let M be a line on x and y, y' its other points (both outside C). Then m_1(y) = 1 + d_y(x) and c_y(x) = 2 + d_y(x) − d_{y'}(x).

*Proof.* The points of C collinear with y are x and the d_y(x) points at distance 2 from x through y. Points of C at distance 3 from y: those at distance 2 from x not behind M, and those at distance 3 from x through y'. So m_3(y) = (6 − d_M(x)) + c_{y'}(x). Put this into m_3(y) = 10 − 2m_1(y) and swap y and y'. ∎

## Step 3: the case analysis

**Case I: some x ∈ C has a(x) = 10.** Then n_2(x) = 16 and n_3(x) = 0. All ten neighbours of x are in C. A point of C at distance 2 from x has a unique middle point y, which is in C and is collinear with it off the line xy, so a(y) ≥ 3 and a(y) = 10. Each neighbour y with a(y) = 10 contributes its 8 neighbours off xy, and these sets are disjoint, so exactly two neighbours y_1, y_2 have a = 10. If they were on different lines on x, a neighbour w of y_2 off xy_2 would satisfy d(y_1, w) = 3 by (G1), against n_3(y_1) = 0. So y_1, y_2 and x form a line L, every point of L has all its neighbours in C, and C contains the 27 points at distance ≤ 1 from L. Since |C| = 27, C is that set.

**Case II: no point has a = 10, some x ∈ C has a(x) = 2.** Let L be the line on x inside C (Lemma B). By Lemma C the 8 points of C ∩ Γ_2(x) are behind one line M ≠ L; call its other points y, y' (outside C). Lemma A gives c_M(x) = 0. Counting as in Lemma D:

* m_1(y) = 1 + d_y(x),
* m_2(y) = 2 + d_{y'}(x) (the two other points of L, the points through y', and c_y(x) = 0),

so m_2 = m_1 + 1 forces d_y(x) = d_{y'}(x) = 4 and m_1(y) = 5. Take z ∈ C collinear with y, z ≠ x, and let N = zy with third point w.

* a(z) = 0: if z were on a line L_z inside C, the other points of L_z would be at distance 3 from x behind M (path x, y, z, q, (G1)), against c_M(x) = 0.
* From z, the point x is behind N, so by Lemma C all 6 points of C ∩ Γ_2(z) are behind N. Through y there are the 4 other points of C collinear with y (one on each line on y other than N), so d_y(z) = 4 and d_w(z) = 2.
* Lemma D at z gives c_y(z) = 2 + 4 − 2 = 4. But the two other points of L (path z, y, x, p) and the 4 points of C collinear with y' off M (path z, y, y', z'') are all at distance 3 from z through y by (G1), so c_y(z) ≥ 6.

Contradiction, so Case II does not occur.

**Case III: a(x) = 0 for every x ∈ C.** Take x ∈ C. By Lemma C its 6 points of C at distance 2 are behind one line M_0(x), so d_y(x) + d_{y'}(x) = 6 for the two other points y, y' of M_0(x). Pick y with d_y(x) ≥ 3 and put k = m_1(y) = 1 + d_y(x) ≥ 4; Step 1 gives k ≤ 5. The k points x_1, …, x_k of C collinear with y lie on distinct lines N_i through y. From x_i there are k − 1 ≥ 1 points of C at distance 2 through y, so N_i = M_0(x_i), and the third point y_i of N_i has d_{y_i}(x_i) = 6 − (k − 1) = 7 − k. Those 7 − k points are collinear with y_i off N_i, so they are at distance 2 from y by (G1), and they are different for different i. So m_2(y) ≥ k(7 − k) ≥ 10, but m_2(y) = k + 1 ≤ 6. Contradiction.

Only Case I is possible, which proves the theorem. ∎

## Remarks

* Uniqueness of L: the points of C with a = 10 are exactly the three points of L.
* For the Ree-Tits octagon this gives exactly 2925 maximum sets, one per line, matching the two computer searches.
* The ratio bound (s+1)(s³+1) and the line-neighbourhood example exist for every order (s, s²), and the proof above is written with the numbers for s = 2. I have not checked whether each step goes through for general s; that is the obvious next thing to try.

## Files

* `clique_certificate.py`, `out_orbits/clique_certificate_output.txt`: second and third independent searches (igraph 1.0.0 largest_cliques, and our own colouring branch and bound) for the maximum sets through point 0; both give exactly the 45 line neighbourhoods. Transitivity of T on points is checked from the generators. Hoffman equality checked on all 2925 lines.
* `max_set_checks.py`, `out_orbits/max_set_checks_output.txt`: exact eigenvalues, the identities of Step 1 and Lemma A on all 2925 line neighbourhoods (280800 (C, x, M) triples).
