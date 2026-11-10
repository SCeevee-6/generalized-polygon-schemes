# Fusions of regular near octagons

A regular near octagon with parameters (s, t2, t3, t) has a distance-regular point graph with c_i = t_i + 1 and
b_i = s(t - t_i) (t_1 = 0, t_4 = t). This note lists all integer parameter sets with s >= 2 and
0 <= t2 <= t3 < t for which the distance scheme has a proper fusion, up to three surfaces that I could not finish.
It replaces the box scan in `near_polygon_note.md`.

## What is known

- Diameter 3: BCN Proposition 4.2.17.
- Diameter 4, Gamma_4 strongly regular: Fiol (arXiv:1407.1392) gives the condition b_3 = a_4 + 1, b_1 = b_3 c_3 for any
  distance-regular graph of diameter 4 (as quoted in the survey of van Dam, Koolen and Tanaka, Section 13.7).
- Diameter 4, Gamma_2 or Gamma_1 u Gamma_2 strongly regular: spectral conditions in arXiv:1608.00373, with H(4,3) as
  an example.
- Dual polar graphs of type B_d/C_d: merging distances 1 and 2 gives the Ustimenko graphs
  (Ivanov, Muzychuk, Ustimenko; BCN 9.4). Fusions of Hamming schemes: Muzychuk.

So every example below is known, and the two-class cases could also be read off from the general criteria. What I
did not find is the list for all 13 partitions in terms of (s, t2, t3, t). Keyword searches only.

## Result

For integers s >= 2, 0 <= t2 <= t3 < t, a partition of {1,2,3,4} closes exactly in these cases:

| parameters (s, t2, t3, t) | partitions | graphs |
|---|---|---|
| A: (s, m, s(m+1), (m+1)(s^2+1) - 1), m >= 0 | {1,2}\|{3,4} | m = s: dual polar graphs B_4(s), C_4(s) |
| B: (s, 0, 0, s) | {1,2,3}\|{4} | generalized octagon with s = t; none exists (Feit-Higman) |
| C: (2, i, 2i, 2i+1), i >= 0 | {1,2,4}\|{3} and {1,4}\|{2}\|{3} | i = 0: thin octagon of order (2,1); i = 1: H(4,3) |
| (2, 1, 2, 3) | also {1,4}\|{2,3} and {1,3,4}\|{2} | H(4,3) |
| D: (3, 1, 2, 3) | {1,3}\|{2,4} | H(4,4) |

with one gap: for each of {1,4}|{2,3}, {1,3}|{2,4} and {1,2,4}|{3} the solution set has one more component, a
surface of high degree, and for those three surfaces I only have a search (below), not a proof. Everything else
is proved, including the six partitions into three blocks and the four two-block partitions {1}|{2,3,4},
{1,2}|{3,4}, {1,3,4}|{2}, {1,2,3}|{4}.

Which of these parameter sets can occur:

- Family C: the eigenvalues are 4(i+1), 2i+3, i+1, -1, -2i-2 and the multiplicity of 2i+3 is 36(i+1)/(i+2)^2, an
  integer only for i = 0, 1, 4. For i = 4 that multiplicity is 5 with valency 20, which Godsil's bound
  k <= (m-1)(m+2)/2 excludes. So only the thin octagon (2,1) and H(4,3) remain.
- Family A: for s <= 12 and m <= 200 the valencies and multiplicities are all integers only for m = s
  (`near_octagon_families_output.txt`, exact). That is a computation, not a proof for all s and m.
- B does not occur. D is H(4,4) or a Doob graph.

So for thick regular near octagons the fusions that actually occur, as far as this goes, are the ones already in
the literature: {1,2}|{3,4} for B_4(q)/C_4(q) and the Hamming fusions for H(4,3) and H(4,4).

## How it is proved

`near_octagon_decompose.py` writes the closure conditions of each partition as polynomials in (s, t2, t3, t) and
splits their common zeros by factoring and eliminating: a factor that has one sign on the region s >= 1,
0 <= t2 <= t3, t >= t3 + 1 is dropped, a reducible equation gives one branch per factor, and an equation that is
linear in a variable is solved for it. Every branch ends as a rational parametrisation, sometimes with a residual
equation, and is substituted back into the original conditions as a check (`near_octagon_decompose_output.txt`,
about 4 minutes). The integer points on the branches are then found by hand:

- {1,2}|{3,4}: t2 = (t - s^2)/(s^2+1), t3 = s(t+1)/(s^2+1). Integrality of t2 gives t = s^2 + m(s^2+1), family A.
- {1,3,4}|{2}: s t (t2+1) = (t3+1)(t2+t3+1) and t + 1 = (s-1)(t2+t3+1). Put n = t2+t3+1. Then n divides s(t2+1),
  say s(t2+1) = jn, and t3 + 1 = jt; since t3 < t, j = 1, and then (s-1)^2 (t2+1) = 2, so (s, t2) = (2, 1).
- {1,2,3}|{4}: s = t and t3 = 0 (or s = 1), in line with Fiol's condition.
- {1,4}|{2,3}: two branches have t3 = 1/(s^3 - 2s^2 + 2s - 2), never an integer for s >= 2. On the third,
  t3 = (D(t-s+1) + t + 1)/(sD) with D = (s-1)t - 1, so D divides t + 1; that leaves (s,t) in
  {(2,2), (2,3), (3,1), (3,2), (4,1)}, and only (2,3) gives an admissible point, H(4,3). The fourth branch is
  the residual surface.
- {1,3}|{2,4}: apart from s = 1 and branches with t2 = t or t = 0, the branches lie on the curve
  t^2 + t(t3+1) = t3(t3+1)^2 with s = (t3+1)^2/t. Its integer points are t3 = n(n+1), t = n(n^2+n+1), and then
  s = (n^2+n+1)/n is an integer only for n = 1: H(4,4). One branch has t3^3 = t3 + 1, no integer root. The last
  branch is the residual surface.
- {1,2,4}|{3}: one branch is t = (s-1)(t3+1), t2 = ((s-1)^2 (t3+1) - 1)/s. Then s divides t3, t3 = sj, and
  t2 <= t3 forces s = 2: family C. The other branch is the residual surface.
- {1,4}|{2}|{3}: one branch is t2 = (t+1-s)/s, t3 = (s-1)t - 1; t3 < t forces s = 2, family C. The other branches
  have t = (2s-1)/(cubic in s), which is an integer only for s = 2 and gives (2,1,2,3) and (2,0,0,1) again.
- {1,3}|{2}|{4}: only s = 1. The other five partitions never close.

Checks (`near_octagon_check.py`, about 12 minutes):

- the families satisfy the closure conditions identically;
- all 44640 integer points with 2 <= s <= 10, 0 <= t2 <= t3 < t <= 30 were tested directly with exact fractions
  (no integrality filter): 53 closing (partition, point) pairs, all of them the ones in the table.

## The three residual surfaces

Each is given by t2 = N/D and one equation R(s, t, t3) = 0 of degree 3 or 4 in t (`near_octagon_decompose.json`).
They do contain real points of the region, so a sign argument cannot work, and I did not find an arithmetic
argument. What was done: for all 2 <= s <= 200 and 0 <= t3 <= 2000 the polynomial R(s, ., t3) was factored exactly
and every integer root t > t3 tested.

- {1,4}|{2,3}: 43 integer roots, all with t3 = s(s+1), t = s(s^2+s+1) and t2 = s^4 + s^3 - 1 > t3, so outside
  the region.
- {1,3}|{2,4}: one integer root, H(4,4), where D = N = 0.
- {1,2,4}|{3}: no integer root.

So there is no further admissible point with s <= 200 and t3 <= 2000. Unproved beyond that.

## Summary

The classification of fusions of regular near octagons by parameters is complete for ten of the thirteen
partitions and reduced to one surface each for the other three. Nothing new exists inside the searched range:
the only feasible parameter sets with a fusion are those of B_4(q)/C_4(q), H(4,3) and H(4,4). This would make a
section in the fusion note once the three surfaces are settled; on its own it is too close to the known
diameter-4 criteria to be a paper.
