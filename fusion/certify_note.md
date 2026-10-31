# Exact re-derivation of the fusion result for generalized hexagons and octagons

`certify_fusions.py` replaces the old `sympy.solve` derivation. It never uses eigenvalues or square roots.

1. The intersection numbers p^k_ij are computed as rational functions of s and t straight from the
   intersection array {s(t+1), st, ..., st; 1, ..., 1, t+1}: A_i is written as a polynomial v_i(A_1) by the
   three-term recurrence, products v_i v_j are reduced modulo the minimal polynomial of A_1, and the result
   is expanded back in the basis v_0, ..., v_d. At (s,t) = (2,4) the valencies are 1, 10, 80, 640, 1024.
2. For each partition of {1, ..., d} the closure conditions (block sums of p^h_ij constant on each block)
   are cleared of denominators. This gives a finite set of polynomials in Q[s,t].
3. A lex Groebner basis over QQ is computed with `sympy.groebner` (exact rational arithmetic). For
   zero-dimensional ideals the last basis element is univariate in t; its real roots are isolated exactly
   (`real_roots`) and s is recovered from the remaining basis elements, exactly when t is rational and to
   60 digits otherwise. For positive-dimensional ideals the factored basis is read by hand below.
4. As an independent check that does not use the Groebner bases, every partition is evaluated at all integer
   pairs 1 <= s, t <= 30.

Output: `certify_fusions_output.txt` (about 15 s).

## Result

For real s, t >= 1 a proper fusion of the point-distance scheme closes only in these cases.

Hexagons (d = 3):

- {1,2}|{3}: basis s t (s - t), so exactly when s = t.
- {1,3}|{2}: basis (s - 1)(st - t^2 - 2t - 1), so s = 1 or s = (t+1)^2/t. The second curve has an integer
  point only for t = 1 (s = 4), so it contains no thick hexagon.
- {1}|{2,3}: basis {1}, never.

Octagons (d = 4):

- {1,2,3}|{4}: basis s t^2 (s - t), so s = t, which the Feit–Higman condition (2st a square) rules out for
  integers s, t >= 1.
- {1,3}|{2,4} and {1,3}|{2}|{4}: s = 1 (the other branches have t = -1).
- {1,2,4}|{3} and {1,4}|{2}|{3}: the only real point with s, t >= 1 is (2,1).
- {1,3,4}|{2}, {1,4}|{2,3}, {1,2}|{3,4}: the real solutions all have s < 1 or t < 1.
- {1}|{2,3,4}, {1}|{2}|{3,4}, {1}|{2,3}|{4}, {1}|{2,4}|{3}, {1,2}|{3}|{4}: none (basis {1}, or s^2 t).

The integer scan agrees: for thick parameters (s, t >= 2) the only closing pairs are s = t, for the
partition {1,2}|{3} (hexagons) and {1,2,3}|{4} (octagons).

So for any distance-regular graph with the intersection array of a generalized hexagon or octagon of order
(s,t), s, t >= 2, the distance scheme has a proper fusion exactly when s = t (hexagons) and never (octagons,
using Feit–Higman). The thin cases s = 1 and (2,1) match `thin/thin_note.md`.

## What this changes

- The theorem in `fusion_note.tex` (Section 3 and appendix) now has an exact proof. The tex file has not
  been edited yet and still describes the old derivation.
- The real exceptional set differs from the table in `docs/novelty_audit.md`: there is no "t = 1 family" for
  octagons, and the hexagon partition {1,3}|{2} has the extra real curve s = (t+1)^2/t, which has no thick
  integer points.

## Status

- Hexagon part, s = t: **known**. For the split Cayley hexagon H(q) the points are those of Q(6,q) and
  distance at most 2 is collinearity in Q(6,q).
- Octagon part (no proper fusion for any thick octagon, in fact for any graph with that array): not found in
  the Crossref/OpenAlex searches in `docs/literature_check_2026-10-31.md`. It is a short consequence of the
  intersection numbers, so it may well be stated somewhere (BCN 1989, Van Maldeghem 1998) in another form. I
  have not checked those books page by page. Treat it as **not known to me**, not as new.
