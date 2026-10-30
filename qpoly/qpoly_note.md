# Q-polynomial orderings of point schemes of generalized 2d-gons, d = 2, 3, 4

Point graph of a generalized 2d-gon of order (s,t): intersection array b = (s(t+1), st, ..., st), c = (1, ..., 1, t+1).
Real parameters s, t >= 1. Eigenvalue indices follow the order below.

| d | eigenvalues θ_0..θ_d | parameter |
|---|---|---|
| 2 | s(t+1), s-1, -t-1 | — |
| 3 | s(t+1), s-1+m, s-1-m, -t-1 | m = √(st) |
| 4 | s(t+1), s-1+m, s-1, s-1-m, -t-1 | m = √(2st) |

## Result (exact; `qpoly_check.py`)

An ordering E_0, E_{σ1}, ..., E_{σd} is Q-polynomial iff two conditions hold:
- q^j_{σ1,i} = 0 whenever |pos(i) − pos(j)| ≥ 2;
- q^{next}_{σ1,i} ≠ 0 for consecutive idempotents.

All d! orderings were tested.

- **d = 2 (GQ):** (E_1 = θ_1) is always Q-polynomial, since q^2_{11} = t(s+1)(t+1)(st+1)/(s+t)^2 > 0. (E_1 = θ_2) is Q-polynomial iff s ≠ 1.
- **d = 3 (hexagon):** Q-polynomial iff **s = 1** or **t = s³**.
  - s = 1: orderings (0,1,2,3), and (0,2,1,3) when t ≠ 1.
  - t = s³ (s ≥ 2): exactly one ordering, (0,3,1,2), i.e. E_1 is the θ = −t−1 idempotent.
  - The dual extremal case s = t³ is *not* Q-polynomial.
- **d = 4 (octagon):** Q-polynomial **only for the ordinary octagon s = t = 1**. This covers all real s,t ≥ 1, including the thin orders (1,t) and (s,1). In particular GO(2,4) has none of its 24 orderings Q-polynomial.

### Proof method
1. Write each required-zero Krein numerator in QQ[s,m] (t = m²/s, resp. m²/(2s)).
2. Factor each numerator into irreducibles. Drop factors that are monomials or have single-sign coefficients, since these have no zero with s,m > 0.
3. Compute a lex Gröbner basis over QQ on every branch. This is exact, with no `sympy.solve`.
4. Check the nonzero conditions modulo each surviving curve.

One curve, ordering (0,2,3,1), d = 3, needs a separate argument:

C = m⁵ − m⁴ + m³s − m²s² + ms³ − s³.

Substitute s = m²/t, then m = 1+u, t = 1+v. This gives

t³C/m⁴ = u³ − u²v + u² + uv³ + 4uv² + 3uv + u + v² + v.

By AM-GM, u³ + 4uv² ≥ 4u²v, so the expression is > 0 for (u,v) ≠ 0. Hence C has no zero with s,t ≥ 1 except s = t = 1, and there the ordering fails anyway.

### Independent check
A float Krein-graph test on all Feit–Higman-admissible integer (s,t) ≤ 40 returns exactly these Q-polynomial cases:
- (1, t) hexagons;
- (2, 8) and (3, 27);
- no octagon at all.

## Literature status — not new
- **Thick octagons:** follows from De Bruyn & Vanhove, *On Q-polynomial regular near 2d-gons*, Combinatorica (2015), doi:10.1007/s00493-014-3039-x. For d ≥ 4, the only thick Q-polynomial regular near polygons are Hamming near polygons and dual polar spaces.
- **Thin octagons (1,t):** also follow from Lewis's girth bound for Q-polynomial DRGs (girth ≤ 6).
- **Hexagons:**
  - Extremal (t = s³) hexagons are known to be Q-polynomial; see De Bruyn, EJC 29(4) (2022) P4.46.
  - The (1,t) case is in Miklavič, arXiv:2501.12820.
  - The full "iff" statement is very likely classical (BCN Ch. 6.5/8), but this has not been confirmed from a primary source.
- **Remark 10 of the 1755-point note:** the octagon Q-polynomial question is closed by the citation above. It is not an open problem.
