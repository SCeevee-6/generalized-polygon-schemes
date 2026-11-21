# Fusions of regular near octagons

A regular near octagon with parameters (s, t2, t3, t) has a distance-regular point graph with c_i = t_i + 1 and
b_i = s(t - t_i) (t_1 = 0, t_4 = t). This note lists all integer parameter sets with s >= 2 and
0 <= t2 <= t3 < t for which the distance scheme has a proper fusion. It replaces the box scan in
`near_polygon_note.md`. The three surfaces left open in the first version are settled in
`near_octagon_residual.py` (last section).

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

and in no other case. This is proved for all 13 partitions. For {1,4}|{2,3}, {1,3}|{2,4} and {1,2,4}|{3}
the proof needs one extra step, for a residual surface of high degree; that step is the last section.

In particular, for the distance-3 graph: Gamma_3 is strongly regular if and only if (s, t2, t3, t) =
(2, i, 2i, 2i+1), and only i = 0 (the octagon of order (2,1), where Gamma_3 has the parameters of T(10)) and
i = 1 (H(4,3)) pass the multiplicity test below.

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
- {1,2,4}|{3}: one branch is t = (s-1)(t3+1), t2 = ((s-1)^2 (t3+1) - 1)/s. Here t2 <= t3
  means (s^2 - 3s + 1)(t3+1) <= 1 - s < 0, so s = 2, t3 = 2t2, t = t3 + 1: family C. The other branch is the
  residual surface.
- {1,4}|{2}|{3}: one branch is t2 = (t+1-s)/s, t3 = (s-1)t - 1; t3 < t forces s = 2, family C. The other branches
  have t = (2s-1)/(cubic in s), which is an integer only for s = 2 and gives (2,1,2,3) and (2,0,0,1) again.
- {1,3}|{2}|{4}: only s = 1. The other five partitions never close.

Checks (`near_octagon_check.py`, about 12 minutes):

- the families satisfy the closure conditions identically;
- all 44640 integer points with 2 <= s <= 10, 0 <= t2 <= t3 < t <= 30 were tested directly with exact fractions
  (no integrality filter): 53 closing (partition, point) pairs, all of them the ones in the table.

## The three residual surfaces

Each is given by t2 = N/D and one equation R(s, t3, t) = 0 of degree 3 or 4 in t (`near_octagon_decompose.json`).
The first version only had a search: for 2 <= s <= 200 and 0 <= t3 <= 2000 every integer root t > t3 of R was
tested, and none is admissible (for {1,4}|{2,3} the 43 integer roots have t3 = s(s+1), t = s(s^2+s+1),
t2 = s^4 + s^3 - 1 > t3, and for {1,3}|{2,4} the only one is H(4,4), where N = D = 0). The surfaces do contain real
admissible points, for s = 3 and 4 only, so a single sign argument can't work. `near_octagon_residual.py` now
proves that none of the three has an admissible integer point (output in `near_octagon_residual_output.txt`,
about 5 minutes).

### {1,2,4}|{3} by hand

Gamma_3 is strongly regular iff p^1_33 = p^2_33 = p^4_33. Directly from the intersection array, the first
equation is s(t - t2)(t2 + 1) E = 0 with E linear in t2, and on E = 0 the second becomes

  (t3 + 1)(t2 + 1) u R = 0,   u = (s-1)(t3+1) - t,

with R cubic in t and leading coefficient s(s^2 - s + 1) (the factor t2 + 1 appears through the numerator of
t2 + 1 on that branch). u = 0 is family C. For R = 0 put n = t3 + 1 and M = (s-1)n + 1, so t + 1 = M - u. Then

  s n t2 = n((s-1)^2 n - 1) - s u (t+1),

so t2 <= t3 is the same as u(t+1) >= pL := n((s^2-3s+1)n + s - 1)/s. For s >= 3 this is positive, so u > 0,
u >= pL/M and t <= M - 1 - pL/M < T := M - pL/M. Now

- R(t3) = -n^2 ((n-1)(s^3 - 3s^2 + 2s - 1) - s(s-1)), negative for s >= 4, n >= 2 and for s = 3, n >= 3;
- 8 R(t3/2) has positive coefficients after s = 3 + z, n = 2 + w, so R(t3/2) > 0;
- R(T) = -Nc/(s^2 M^3) with Nc of degree 16; Nc has positive coefficients after s = 5 + z, n = 5 + w, and in one
  variable for s = 4, n >= 7 and s = 3, n >= 122.

A cubic with positive leading coefficient that is positive at t3/2 and negative at t3 and at T has no root in
[t3, T]: two roots there would force a third beyond T, and then it would be negative on all of (-inf, t3].
What is left is finite or one-dimensional. For s = 2, R(2, t3+1+x, t3) has positive coefficients. For t3 <= 3,
t2 runs over 0..t3 and the resultant in t leaves no admissible point. For s = 3, t3 <= 120 and s = 4, t3 <= 5 the
integer roots of R were checked exactly. None is admissible.

### All three surfaces: a region argument

For fixed (s, t3) count the real roots t of R with t > t3 and 0 <= N/D <= t3. This number can only change where
one of lc_t(R), disc_t(R), R(s,t3,t3), Res_t(R,N), Res_t(R, N - t3 D), Res_t(R,D) vanishes. The script factors
these and checks that every factor has constant sign on a region by substituting a parametrisation of the region
and checking that all coefficients have one sign. Then one sample point gives the count everywhere in the
region. For fixed s (or fixed t3) the same thing is done in one variable: beyond the largest real root of the
boundary factors one sample decides, and the integer points below it are checked exactly.

| partition | regions certified (count 0) | fixed lines checked exactly |
|---|---|---|
| {1,2,4}\|{3} | s >= 4, t3 >= 2 | s = 2, 3; t3 = 0, 1 |
| {1,4}\|{2,3} | s >= 4, t3 >= 1 | s = 2, 3; t3 = 0 |
| {1,3}\|{2,4} | s >= 6 split into t3 >= s^2-s-2, s <= t3 <= s^2-s-3, 1 <= t3 <= s | s = 2..5; t3 = 0 |

For {1,3}|{2,4} two things needed care. The factor s(s-1)t3 - s - (t3+1)^2 is where R, N and D have a common
root; its large root in t3 lies strictly between s^2-s-3 and s^2-s-2 for s >= 4, so it never vanishes inside
the regions. And a degree-20 factor of the discriminant has leading form s^12 t3^6 (s - t3)^2 for t3 of order s,
so plain coefficient checks fail on any region containing the line t3 = s; splitting there fixes it. The only
integer point found on a line is (s, t3, t) = (3, 2, 3) with N = D = 0, which is H(4,4) on the branch
s = (t3+1)^2/t treated above.

This second argument also re-proves {1,2,4}|{3}. Both arguments rest on exact sympy arithmetic (resultants,
factorisations, real root isolation) and on the factor-and-split case analysis of `near_octagon_decompose.py`
for the other branches.

## Summary

The fusions of regular near octagons with s >= 2 are now classified by parameters for all thirteen partitions.
The only parameter sets with a fusion that pass the multiplicity test are those of B_4(q)/C_4(q) (for family A
only checked for s <= 12, m <= 200), H(4,3), H(4,4) and the octagon of order (2,1); none of these is new. The
statement that matters most is the distance-3 one: Gamma_3 is strongly regular only for (2, i, 2i, 2i+1). I
haven't found this, or the two partitions {1,4}|{2,3} and {1,3}|{2,4}, treated for regular near octagons of
arbitrary order (keyword searches only); Fiol's criterion covers Gamma_4, and BCN 4.2.17 covers diameter 3.
It is a classification result with an empty answer, so it would fit as a section of the fusion note rather than as
a paper of its own.
