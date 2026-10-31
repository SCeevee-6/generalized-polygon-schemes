# Fusions of regular near hexagons and near octagons (angle 3)

Labels: **known**, **not found** (not in the sources I checked), **unverified**.

## Literature

- Diameter 3 is settled in general. For any DRG of diameter 3, BCN Proposition 4.2.17 (also Brouwer–Van Maldeghem,
  *Strongly regular graphs*, Prop. 1.3.12) says exactly when Gamma_2 or Gamma_3 is strongly regular. A regular near
  hexagon is a special case, so the near-hexagon statement below is a direct corollary. **known**
- Dual polar graphs: the B_d(q)/C_d(q) dual polar graph and the half dual polar graph D_{d+1,d+1}(q) live on the same
  vertex set (BCN §9.4), which explains the strongly regular fusion found below for d = 4. **known** (identification
  from BCN, not re-derived here: **unverified** by me in detail).
- Hamming graphs: fusions of H(n,q) are classical. **known**

## Computation

`near_polygon_fusions.py` computes the p^k_ij of a regular near 2d-gon with parameters (s, t_2, ..., t_{d-1}, t)
(c_i = t_i + 1, b_i = s(t - t_i)) as rational functions, from the tridiagonal intersection matrix. At (2,0,0,4) it gives
the GO(2,4) valencies, and at three random integer points all 125 numbers agree with the independent Fraction code in
`drgfusion/drg_table_fusions.py`. For every partition it prints a Groebner basis (grevlex, QQ) of the closure ideal,
and it scans all integer points with s <= 8, t <= 24, t_i nondecreasing and integral valencies. Output:
`near_polygon_fusions_output.txt` (the d = 4 bases are long and not simplified).

## Results

Near hexagons (d = 3), from the factored bases:

- {1,2}|{3}: the basis is s(t - t_2)(s t_2 + s - t), so for s >= 1 and t > t_2 it closes iff t = s(t_2 + 1). On t = s(t_2+1) with integer multiplicities in the box:
  the generalized hexagon arrays (s,0,s), s = 2..8, the dual polar graphs C_3(2), C_3(3), C_3(4) (t_2 = q,
  t = q^2 + q), and (5,1,10), whose feasibility I did not check further. **known** (BCN 4.2.17).
- {1,3}|{2} closes iff s = 1 or s t (t_2 + 1) = (t + 1)(t + t_2 + 1). For t_2 = 0 this is the hexagon curve
  s = (t+1)^2/t. The only thick integer point in the box, (3,1,2), has non-integral valencies.
- {1}|{2,3} never closes.

Near octagons (d = 4), integer scan only (s <= 8, t <= 24), thick points (s, t >= 2):

| partition | thick parameter points with integral valencies | multiplicities |
|---|---|---|
| {1,2,3}\|{4} | (s,0,0,s), s = 2..8 (generalized octagons with s = t) | non-integral (Feit–Higman) |
| {1,3}\|{2,4} | (3,1,2,3) = H(4,4) | integral |
| {1,2}\|{3,4} | (2,2,6,14) = C_4(2) dual polar; (2,1,4,9) | integral; non-integral |
| {1,4}\|{2,3}, {1,3,4}\|{2} | (2,1,2,3) = H(4,3) | integral |
| {1,2,4}\|{3}, {1,4}\|{2}\|{3} | (2,i,2i,2i+1), i = 1..11 | integral only for i = 1 (H(4,3)) and i = 4; i = 4 fails Godsil's bound (multiplicity 5 < valency 20 needs k <= 14) |

A direct check of the rank-4 dual polar families shows {1,2}|{3,4} closing for B_4(q)/C_4(q) (e = 1) at every q tested
(2,3,4,5,7,8,9), only bipartite-type fusions for D_4(q) (e = 0), and none for e = 1/2, 3/2, 2.

## Assessment

Everything thick that closes in the box is a generalized octagon with s = t (excluded by Feit–Higman), a Hamming
graph, or a B_4/C_4 dual polar graph. All of these are known. Turning the box scan into a theorem for all
(s, t_2, t_3, t) would need a proper decomposition of the d = 4 ideals, which the grevlex bases do not give directly.
I stopped here: the expected result is "only the known families", and the diameter-3 case is already in BCN.
