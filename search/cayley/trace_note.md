# Character-trace test for Cayley graphs with array {52,42,16;1,6,28}

`cayley_trace.py` is a second implementation of the trace model. It takes the rational class sizes, the Galois orbit
sizes and the rationalized character values from the GAP export and solves for

* w[K][i] = |K ∩ S_i| for every nonidentity rational class K and distance i = 1,2,3,
* y[r][j] = multiplicity of the nonprincipal eigenvalue θ_j in one character of the Galois orbit r,

subject to

    Σ_K T_r(K) w[K][i] = n_r Σ_j P_ji y[r][j],     Σ_r n_r d_r y[r][j] = m_j,
    Σ_i w[K][i] = |K|,   Σ_K w[K][i] = k_i,   Σ_j y[r][j] = d_r.

The derivation: in the regular representation S_i acts as p_i(A_1), A_1 restricted to a nontrivial irreducible block
is Hermitian with eigenvalues among 17, 2, −8 (the eigenvalue 52 only occurs on the trivial block since the graph is
connected), and conjugate characters have the same multiplicities because the eigenvalues are rational. Nothing is
assumed about S beyond inverse closure, so S may split ordinary classes. The S3 6-cycle, whose generators split a
class, comes out feasible.

With the HiGHS MILP inside SciPy all ten nonabelian groups of order 625 are infeasible (each in under 2 s). The LP
relaxations are feasible and the equality system has integer solutions once the bounds are dropped, so the
contradiction really uses the box constraints. The same model over the exported tables at orders 729, 750 and 1000 is still running. The groups that export_normal.g skipped (normal cyclic Sylow subgroup) need `export_skipped.g` run in GAP first.

## Groups 6, 9, 10 and 14

These four tables have one Galois orbit of 20 characters of degree 5 whose orbit sum is 100 at the identity, −25 on
the four nontrivial elements of a central subgroup Z of order 5 and 0 elsewhere. Every other nontrivial character is
linear or lies in orbits of 4 characters of degree 5.

Write (a, b, c) for the multiplicities of 17, 2, −8 in one character of the 20-orbit, so a + b + c = 5. In the
multiplicity equations Σ_r n_r d_r y[r][j] = m_j the 20-orbit contributes 100 times its multiplicity and everything
else contributes multiples of 1 (linear characters) or 20 (the 4-orbits of degree 5; only present for 9 and 10).
For θ = 17 the total is 52 < 100, so a = 0. For θ = −8 the total is 208, so c ≤ 2, and for θ = 2 the total is 364,
so b ≤ 3. With b + c = 5 this forces b = 3, c = 2. (For groups 9 and 10 there is also a 4-orbit of degree 5, with multiplicities y = (y17, y2, y−8), and 24 linear
characters. The θ = 17 equation reads c17 + 20 y17 + 100a = 52 with c17 ≤ 24, so a = 0, y17 = 2, c17 = 12. Then
y−8 ≤ 3, and the θ = −8 equation c−8 + 20 y−8 + 100c = 208 with c−8 ≤ 12 leaves only c = 2, y−8 = 0, c−8 = 8.)

The trace equation for this orbit at i = 1 reads −25 |S_1 ∩ Z#| = 20 (2·3 − 8·2) = −200, so |S_1 ∩ Z#| = 8.
But |Z#| = 4. So none of these four groups carries a Cayley graph with this array.

## Groups 7 and 8

The two tables coincide. G' has order 25 and consists of the identity, a class K1 of 4 central elements and a class
K2 of 20. There are 24 linear characters (six orbits of 4) and six orbits of 4 characters of degree 5. The sum of all
24 degree-5 characters is 100 at the identity, −5 on G'# and 0 elsewhere (read off the table).

The multiplicity equations determine the totals Y_θ of the multiplicities over the six degree-5 orbits: with c_θ
linear characters (a multiple of 4, at most 24), c_17 + 20 Y_17 = 52 gives Y_17 = 2, c_−8 + 20 Y_−8 = 208 gives
Y_−8 = 10, and then Y_2 = 18. Summing the trace equations of the six orbits at distance i,

    −5 |S_i ∩ G'#| = 4 Σ_j P_ji Y_j.

For i = 3 the right side is 4(−32·2 + 8·18 − 7·10) = 40, so |S_3 ∩ G'#| = −8. (At i = 1 and 2 it gives 8 and 24,
which add up correctly to 24 with the impossible −8.)

## Groups 3 and 4

The tables coincide. The centre Z has order 25: the identity and six classes K1, …, K6 of 4 elements. There are five
orbits of 4 degree-5 characters; the orbit sum of each is 20 at e, 20 on one class K_a, −5 on the other five central
classes and 0 off Z, and a runs over K1, K3, K4, K5, K6. There are five orbits of 20 linear characters whose orbit
sums are 20 on e and K2, −5 on the other central classes, 20 on one of the five classes of size 20, −5 on the other
four and 0 on the classes of size 100, and six orbits of 4 linear characters. Only distance 3 is needed. Put w_K =
|K ∩ S_3|, z = |S_3 ∩ Z#| and t = |S_3 ∩ (classes of size 20)|.

* A degree-5 orbit with multiplicities (y17, y2, y−8) gives 25 w_a − 5z = 4(−32 y17 + 8 y2 − 7 y−8), i.e.
  5 w_a − z = 32 − 32 y17 − 12 y−8. Summing over the five orbits: 32 Y17 + 12 Y−8 = 160 + 5 w_K2, so w_K2 ∈ {0, 4}.
* A linear 20-orbit with eigenvalue θ gives 5 w_K2 + 5 w_c − z − t = 4 P_θ3; summing over the five orbits:
  25 w_K2 − 5z = 4 Σ P_θ3 over those orbits.
* The multiplicity of −8: if u of the 4-orbits and v of the 20-orbits are at −8 then 4u + 20v + 20 Y−8 = 208, so
  u = 2 and v + Y−8 = 10.

If w_K2 = 0 then 8 Y17 + 3 Y−8 = 40 with Y−8 ≥ 5 forces Y−8 = 8 and v = 2; the 20-orbit sum then allows only z = 24,
which needs w_K2 = 4. If w_K2 = 4 then 8 Y17 + 3 Y−8 = 45 forces Y−8 = 7, Y17 = 3, v = 3, and the 20-orbit sum again
forces z = 24, so every w_a = 4. The first bullet now reads 8 y17 + 3 y−8 = 9 for every degree-5 orbit, so y17 = 0
in each, contradicting Y17 = 3.

## Groups 12 and 13

These two tables coincide too and have 37 rational classes. I have not finished a hand argument; their exclusion rests
on the solver (HiGHS MILP; the earlier CP-SAT and Z3 runs agree).

## Groups of order 750 and 1000 with a normal cyclic Sylow subgroup

The normal-Cayley search skipped the 29 nonabelian groups of order 750 and the 13 of order 1000 that have a normal
cyclic Sylow subgroup, because Chen–Muzychuk–Ponomarenko already rule out a primitive distance-regular Cayley graph on
them. As an independent check I ran the trace model on them anyway. `export_skipped.g` produced their rational tables
(`tables_skipped.txt`, GAP 4.12.1, all 42 groups). Each table was tested against {49,40,22;1,5,28} at order 750 and
against both {54,42,11;1,3,44} and {63,48,24;1,6,28} at order 1000, first with HiGHS (600 s) and then with CP-SAT
(900 s); results in `trace_skipped.json`.

| order | array | tables | both solvers infeasible | undecided within the limits |
|---|---|---:|---:|---:|
| 750 | {49,40,22;1,5,28} | 29 | 19 | 10 |
| 1000 | {54,42,11;1,3,44} | 13 | 11 | 2 |
| 1000 | {63,48,24;1,6,28} | 13 | 13 | 0 |

No run found a feasible point. The undecided ones (750: 7, 13, 31, 32, 33, 34, 35, 36, 37, 38; 1000: 22, 23 for the
first array) were stopped at 1700 s with neither solver finished. They are still excluded by the
Chen–Muzychuk–Ponomarenko theorem, so this only means the computational cross-check is incomplete for them.
