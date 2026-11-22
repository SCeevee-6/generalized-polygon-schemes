# Maximum sets of pairwise non-opposite points in finite generalized octagons

Status: new as far as I can tell. I haven't found it in the EKR-for-buildings literature I've read so far (polar spaces, projective spaces, flags; De Beule-Mattheus-Metsch, JCTA 2022), but a proper literature check is still to be done. Hand proof, not refereed.

The first version of this file only treated order (s, s²). The same argument works for every order once t is kept separate from s; the only new ingredient is a short integrality argument that rules out |s − t| = 1.

Checks:
- `max_set_all_orders_checks.py` (sympy) verifies every identity and inequality below as an exact statement in s and t, and the |s − t| = 1 non-integrality for s ≤ 10⁴ as well as symbolically.
- `max_set_general_checks.py` is the older check for t = s².
- `max_set_checks.py` checks the lemmas numerically on the Ree-Tits octagon O(2), order (2,4), on all 2925 line neighbourhoods.
- `escape_bk.py` enumerates all inclusion-maximal non-opposite sets of size ≥ 15 through a point of O(2), and of size ≥ 25 through a point of its dual, order (4,2). The maximum ones are exactly the line neighbourhoods in both cases, so both sign cases t > s and t < s of the proof are confirmed on a real octagon. See `escape_note.md`.

**Theorem.** Let S be a finite generalized octagon of order (s, t) with s, t ≥ 2, and let C be a set of points, no two of them opposite. Then |C| ≤ (s+1)(st+1). Equality holds if and only if C is the set of points at distance at most 1 from some line L. L is unique.

Applied to the dual: a family of pairwise non-opposite lines has at most (t+1)(st+1) members, with equality only for the lines through a point p together with the lines meeting one of them.

The numbers in brackets are the values for O(2), (s,t) = (2,4).

## Setup

Distances are in the collinearity graph, so they run from 0 to 4 and "opposite" means distance 4. The graph is distance-regular with intersection array {s(t+1), st, st, st; 1, 1, 1, t+1}, so

- k = s(t+1) [10] and k_4 = s⁴t³ [1024];
- a_1 = a_2 = a_3 = s − 1;
- v = (s+1)(st+1)(s²t²+1) [1755].

A_i is the distance-i matrix. Put r = √(2st). The eigenvalues of A_1 are k, θ_1 = s−1+r, θ_2 = s−1, θ_3 = s−1−r, θ_4 = −(t+1) [10, 5, 1, −3, −5], and A_4 = P_4(A_1) takes the values

| θ | k | θ_1 | θ_2 | θ_3 | θ_4 |
|---|---|---|---|---|---|
| P_4(θ) | s⁴t³ | −s²t | s²t | −s²t | t³ |

Two facts about generalized octagons are used throughout:

* (G1) Let x_0 ∼ x_1 ∼ … ∼ x_r be a path with r ≤ 4 in which consecutive steps use different lines. Then it is a geodesic, so d(x_0, x_r) = r. (It is a non-backtracking path of length 2r ≤ 8 in the incidence graph, which has girth 16.)
* (G2) If d(x,z) ≤ 3, the geodesic from x to z is unique. For a line M on x, z is *behind M* (from x) if that geodesic starts along M.

Notation:
- For x ∈ C write n_i(x) = |C ∩ Γ_i(x)| and a(x) = n_1(x).
- For y ∉ C write m_i(y) = |C ∩ Γ_i(y)|.
- For x ∈ C and a line M on x, d_M(x) and c_M(x) are the numbers of points of C at distance 2, resp. 3, from x behind M.
- For y ∈ M \ {x}, d_y(x) and c_y(x) count the points among those whose geodesic from x passes through y.

Put K = (s+1)(st+1) [27].

## Step 1: the ratio bound and what equality gives

C is a coclique of the opposition graph Γ_4. That graph is s⁴t³-regular with least eigenvalue −s²t (the other values s²t and t³ are positive), so |C| ≤ v·s²t/(s⁴t³ + s²t) = v/(1 + s²t²) = K. The K points at distance ≤ 1 from a line have pairwise distance ≤ 3, so the bound is attained. Now suppose |C| = K.

Equality in the ratio bound puts 1_C − (K/v)·1 in the −s²t-eigenspace of A_4, which is the sum of the A_1-eigenspaces for θ_1 and θ_3. So:

* A_4 1_C = s²t(1 − 1_C).
* (A_1 − θ_1)(A_1 − θ_3) 1_C is a multiple of 1. Using A_1² = kI + (s−1)A_1 + A_2, this operator is A_2 − (s−1)A_1 + (s² − s + 1 − st)I [A_2 − A_1 − 5I]. Its value on 1 is (k − θ_1)(k − θ_3) = s²t² + 1, so it maps 1_C to (K/v)(s²t² + 1)·1 = 1.

Reading these off at a point, together with n_0 + … + n_4 = K:

* for x ∈ C: n_2 = (s−1)a + s(t − s + 1) [a + 6], n_3 = s²t + s² − s·a [20 − 2a], n_4 = 0;
* for y ∉ C: m_2 = (s−1)m_1 + 1 [m_1 + 1], m_3 = st + s − s·m_1 [10 − 2m_1], m_4 = s²t [16]. In particular m_1(y) ≤ t + 1.

## Step 2: local lemmas

**Lemma A.** Let x ∈ C, let M be a line on x that is not contained in C, and let y ∈ M \ C. Then c_M(x) = s(s − a(x)) [4 − 2a].

*Proof.* Take z ∈ C.
- If d(x,z) ≤ 2, or d(x,z) = 3 and z is behind M, then d(y,z) ≤ 3.
- If d(x,z) = 3 and z is not behind M, the path y, x, … , z is non-backtracking of length 4, so d(y,z) = 4 by (G1).

Since C ∩ Γ_4(x) = ∅, m_4(y) = n_3(x) − c_M(x), i.e. s²t = s²t + s² − s·a(x) − c_M(x). ∎

**Lemma B.** For x ∈ C, a(x) ∈ {0, s, k}.
- If a(x) = s, exactly one line on x lies in C, and the other lines on x meet C only in x.
- If a(x) = k, all lines on x lie in C.

*Proof.* Let f be the number of lines on x inside C. Every other line M on x has r_M = |M ∩ C| − 1 ∈ [0, s−1], and a = s·f + Σ r_M.
- If some line on x is not inside C, Lemma A gives s(s − a) ≥ 0, so a ≤ s. Hence a > s forces every line on x into C, and a = k.
- If f ≥ 1 and some line is not inside C, then a ≥ s and a ≤ s, so a = s, f = 1 and every r_M = 0.
- If f = 0, Lemma A applies to all t+1 lines, and Σ_M c_M = n_3 reads (t+1)·s(s − a) = s²t + s² − s·a, i.e. sta = 0. ∎

So if no point of C has a = k, every collinear pair in C lies on a line contained in C. It follows that every line through a point y ∉ C carries at most one point of C.

**Lemma C.** Let x ∈ C. If two points z, z' of C are at distance 2 from x and behind different lines on x, then z and z' are opposite. Hence, if a(x) ≤ s and no point of C collinear with x has a = k, all points of C ∩ Γ_2(x) are behind one line on x.

*Proof.* z, y, x, y', z' is non-backtracking of length 4 (z is not on xy because d(x,z) = 2), so d(z,z') = 4 by (G1).

For the second sentence, suppose a point z of C at distance 2 from x lies behind a line inside C, through a point p of that line. Then z is collinear with p off the line, so two lines on p meet C in at least two points. Then a(p) > s, so a(p) = k by Lemma B. ∎

**Lemma D.** Let z ∈ C with a(z) = 0, let N be a line on z, and let w ∈ N \ {z} (all points of N \ {z} are outside C). Then m_1(w) = 1 + d_w(z) and c_w(z) = s + s·d_w(z) − d_N(z).

*Proof.* The points of C collinear with w are z and the d_w(z) points at distance 2 from z through w. The points of C at distance 3 from w are of two kinds:
- those at distance 2 from z not behind N;
- those at distance 3 from z through a point of N other than w.

So m_3(w) = (n_2(z) − d_N(z)) + (c_N(z) − c_w(z)). By Step 1 n_2(z) = s(t − s + 1), and by Lemma A c_N(z) = s². Put these into m_3(w) = st + s − s·m_1(w). ∎

**Lemma E.** There is no finite generalized octagon of order (s, t) with s, t ≥ 2 and |s − t| ≤ 1.

*Proof.* s = t is excluded by Feit-Higman (2st must be a square). For θ_4 = −(t+1) the standard sequence is u_i = (−1/s)^i, and the multiplicity v/Σ k_i u_i² simplifies to

  s⁴(st+1)(s²t²+1) / ((s+t)(s²+t²)).

If |s − t| = 1 then d = s² + t² = 2st + 1 is odd, prime to s (gcd(s,t) = 1) and prime to st + 1 (2(st+1) = d + 1). Integrality forces d | s²t² + 1, hence d | 4(s²t² + 1) = (d − 1)² + 4, hence d | 5. But d ≥ 13. ∎

(Feit-Higman's square condition alone does not exclude |s − t| = 1: (s,t) = (8,9) has 2st = 144.)

## Step 3: the case analysis

By Lemma E, |s − t| ≥ 2: either t ≥ s + 2 or t ≤ s − 2.

**Case I: some x ∈ C has a(x) = k.** Then n_2(x) = s²t and n_3(x) = 0, and all k neighbours of x are in C.
- A point of C at distance 2 from x has a unique middle point y. That point is in C and is collinear with it off the line xy, so a(y) > s and therefore a(y) = k.
- Each neighbour y with a(y) = k contributes its k − s = st neighbours off xy, and these sets are disjoint. So exactly s²t/st = s neighbours of x have a = k.
- If two of them, y_1 and y_2, were on different lines on x, a neighbour w of y_2 off xy_2 would satisfy d(y_1, w) = 3 by (G1). That contradicts n_3(y_1) = 0.

So these s neighbours and x form a line L, and every point of L has all its neighbours in C. Then C contains the K points at distance ≤ 1 from L, so C is that set.

**Case II: no point has a = k, and some x ∈ C has a(x) = s.** Let L be the line on x inside C (Lemma B).
- n_2(x) = st. By Lemma C these points are behind one line M ≠ L. Its other s points are outside C, and Lemma A gives c_M(x) = 0.
- For y ∈ M \ {x}, m_1(y) = 1 + d_y(x), and m_2(y) = s + (st − d_y(x)): the s points of L \ {x}, the points through the other points of M, and c_y(x) = 0. So m_2 = (s−1)m_1 + 1 forces d_y(x) = t, and m_1(y) = t + 1.

Now take z ∈ C collinear with y, z ≠ x, and let N = zy.
- a(z) = 0. Otherwise z lies on a line L_z inside C, and the other points of L_z would be at distance 3 from x behind M (path x, y, z, q, (G1)), against c_M(x) = 0.
- x is at distance 2 from z, so n_2(z) = s(t − s + 1) ≥ 1. This already excludes t ≤ s − 2.
- From z, the point x is behind N. So by Lemma C all s(t − s + 1) points of C ∩ Γ_2(z) are behind N, i.e. d_N(z) = s(t − s + 1). Through y there are the t other points of C collinear with y, one on each line on y other than N, so d_y(z) = t.
- Lemma D at z gives c_y(z) = s + st − s(t − s + 1) = s².
- But the following points are all at distance 3 from z through y by (G1): the s points of L \ {x} (path z, y, x, p), and for each of the s − 1 points y' ∈ M \ {x, y}, the t points of C collinear with y' off M (path z, y, y', z''). So c_y(z) ≥ s + (s−1)t, which exceeds s² by (s−1)(t−s) > 0 when t ≥ s + 2.

So Case II does not occur.

**Case III: a(x) = 0 for every x ∈ C.** Then n_2 = s(t − s + 1) ≥ 0, which excludes t ≤ s − 2. So t ≥ s + 2.
- Take x ∈ C. By Lemma C its s(t − s + 1) points of C at distance 2 are behind one line M_0(x). So the s points y ∈ M_0(x) \ {x} have Σ d_y(x) = s(t − s + 1).
- Pick y with d_y(x) ≥ t − s + 1, and put h = m_1(y) = 1 + d_y(x) ≥ t − s + 2 ≥ 4. Step 1 gives h ≤ t + 1.
- The h points x_1, …, x_h of C collinear with y lie on distinct lines N_i through y. From x_i there are h − 1 ≥ 1 points of C at distance 2 through y, so N_i = M_0(x_i). Therefore the other s − 1 points of N_i carry, between them, s(t − s + 1) − (h − 1) points of C at distance 2 from x_i.
- Those points are collinear with a point of N_i \ {x_i, y} off N_i, so they are at distance 2 from y by (G1). They are different for different i, by uniqueness of geodesics.

So m_2(y) ≥ h(s(t − s + 1) + 1 − h), while m_2(y) = (s−1)h + 1. The difference is h(st − s² + 2 − h) − 1, and since h ≤ t + 1 it is at least h(s−1)(t−s−1) − 1 ≥ 4 − 1 > 0. This is a contradiction.

Only Case I is possible, which proves the theorem. ∎

## Remarks

* Uniqueness of L: the points of C with a = k are exactly the s + 1 points of L.
* For O(2) this gives exactly 2925 maximum sets, one per line, and for its dual exactly 1755, one per point of O(2). Both match the computer searches.
* What is used about the geometry: (G1), (G2) and the parameters. What makes Step 1 work is that the least eigenvalue −s²t of A_4 occurs at θ_1 and θ_3 and nowhere else, for every order.
* Only the known orders (q, q²) and (q², q) are realised, but the proof doesn't use that.

## Files

* `max_set_all_orders_checks.py`, `out_orbits/max_set_all_orders_checks_output.txt`: exact symbolic checks in s and t.
* `max_set_general_checks.py`, `out_orbits/max_set_general_checks_output.txt`: the older checks for t = s².
* `max_set_checks.py`, `out_orbits/max_set_checks_output.txt`: numerical checks of Step 1 and Lemma A on O(2) (all 2925 line neighbourhoods, 280800 (C, x, M) triples).
* `clique_certificate.py`, `out_orbits/clique_certificate_output.txt`: independent clique searches for O(2) (igraph 1.0.0 `largest_cliques`, and our own colouring branch and bound). Both find exactly the 45 line neighbourhoods through point 0. T is transitive on points, so there are 2925 in all.
* `escape_bk.py`, `out_escape/`: all maximal non-opposite sets of size ≥ 15 through a point of O(2) and of size ≥ 25 through a point of the dual (`escape_note.md`).
