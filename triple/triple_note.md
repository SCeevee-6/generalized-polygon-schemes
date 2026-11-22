# Triple intersection numbers for GO(2,4) and the vanishing Krein parameter

## Literature

- Uniqueness of GO(2,4) is open. Partial results: De Bruyn, Discrete Math. 338 (2015) 2125-2142, doi 10.1016/j.disc.2015.04.022 (the Tits octagon is the only GO(2,4) containing a sub-GO(2,1)); Cohen, O'Brien, Shpectorov, J. Algebra 421 (2015) 369-393, doi 10.1016/j.jalgebra.2014.08.034 (unique under a condition on a point stabiliser); and uniqueness among those with a vertex-transitive automorphism group (cited in the same paper).
- Bamberg and Lansdown, Algebr. Comb. 6 (2023), doi 10.5802/alco.246, Thm 8.1, use the same zero (written q^2_22 in their ordering) to show that a generalized octagon of order (s,s^2) has no nontrivial m-ovoid.
- Petit and Van Maldeghem, arXiv:2607.18560 (July 2026), characterise the natural embedding of the Ree-Tits octagons in metasymplectic spaces of type F4. For a pseudo-isometric embedding, distances 1, 2, 3, 4 in the octagon are collinear, symplectic, special and opposite in the ambient space.
- q^4_44 = 0 for extremal generalized octagons, and the triple-intersection-number method of Coolsaet and Jurisic for vanishing Krein parameters, are both in the literature (De Bruyn, EJC 29(4) P4.46 (2022), item (4) of the introduction; van Dam-Koolen-Tanaka survey 6.3.1). Vidali's `sage-drg` automates the method.

## Computation

`go24_triple_lp.py` sets up, for each of the 17 triple types (d(u,v), d(u,w), d(v,w)), the linear system for the numbers [i j k] (marginals equal to intersection numbers, plus the Krein equation sum Q_{i4} Q_{j4} Q_{k4} [i j k] = 0), and computes the range of every [i j k] by LP, with and without the Krein equation. It then counts the same numbers for one triple of each type in the actual octagon (ATLAS 1755-point action).

## Result

- The Krein equation cuts the solution space by exactly one dimension for every type, leaving 2 to 19 free parameters.
- It fixes values that are otherwise free only for the type (1,1,2) (a path v - u - w with v, w not collinear): [3 3 3] = 0, [3 3 4] = [3 4 3] = 64, [3 4 4] = 384, and [4 3 3] = [4 3 4] = [4 4 3] = [4 4 4] = 256. These are what the near-polygon property of the octagon gives (a point at distance 3 from u reaches u through one line only), and the octagon itself has these values.
- For all other types the Krein equation only narrows the LP ranges; the octagon's values lie inside them.

## Which Krein parameter (`krein_convention.py`)

Bamberg and Lansdown write the zero as q^2_22. Their E_2 is the eigenspace for which a line (inner distribution (1, s, 0, 0, 0)) is a {2}-design, and (aQ)_j = m_j(1 + s theta_j/k) vanishes only for theta_j = -(t+1). Here the eigenvalues are ordered 10, 5, 1, -3, -5, so their E_2 is my E_4 (eigenvalue -5, multiplicity 78), and their q^2_22 is my q^4_44. For general (s,t) and theta = -(t+1), sum_l k_l u_l^3 = (s-1)(s+1)(s^2-t)(s^4+t^2)/s^8, the factors they give. At (2,4) it is the only vanishing Krein parameter with i, j, h nonzero. So both papers use the same zero and nothing new comes from it here.

## Triple orbits (`triple_orbits.py`)

The LP above uses one triple per distance type. Here every orbit of T, and of 2F4(2), on ordered triples is computed (orbits of the point stabiliser, order 10240 or 20480 checked, on ordered pairs), with the full vector [i j k] for each orbit.

T has 415 orbits on ordered triples of points and 2F4(2) has 258. On the 17 types with U <= V <= W (LP free dimension counted with the Krein equation):

| type | unknowns | LP free dim | T-orbits | 2F4(2)-orbits | distinct vectors | affine span |
|---|---:|---:|---:|---:|---:|---:|
| (1, 1, 1) | 25 | 2 | 1 | 1 | 1 | 0 |
| (1, 1, 2) | 29 | 2 | 1 | 1 | 1 | 0 |
| (1, 2, 2) | 37 | 6 | 1 | 1 | 1 | 0 |
| (1, 2, 3) | 37 | 5 | 1 | 1 | 1 | 0 |
| (1, 3, 3) | 40 | 7 | 1 | 1 | 1 | 0 |
| (1, 3, 4) | 36 | 5 | 1 | 1 | 1 | 0 |
| (1, 4, 4) | 34 | 5 | 1 | 1 | 1 | 0 |
| (2, 2, 2) | 49 | 14 | 3 | 3 | 1 | 0 |
| (2, 2, 3) | 50 | 14 | 1 | 1 | 1 | 0 |
| (2, 2, 4) | 44 | 10 | 1 | 1 | 1 | 0 |
| (2, 3, 3) | 54 | 17 | 6 | 5 | 4 | 2 |
| (2, 3, 4) | 48 | 13 | 3 | 2 | 2 | 1 |
| (2, 4, 4) | 44 | 11 | 8 | 5 | 4 | 2 |
| (3, 3, 3) | 57 | 19 | 18 | 9 | 6 | 5 |
| (3, 3, 4) | 51 | 15 | 28 | 17 | 14 | 8 |
| (3, 4, 4) | 46 | 12 | 40 | 20 | 18 | 9 |
| (4, 4, 4) | 42 | 10 | 86 | 49 | 44 | 10 |

- Every orbit's vector satisfies the LP equations exactly and lies in the LP ranges, as it must.
- The ten types up to (2,2,4) are realised by a single vector in the octagon, although the LP leaves 2 to 14 free dimensions. For (2,2,2) there are three T-orbits (also three 2F4(2)-orbits) with the same vector. So for these types the octagon is far more rigid than the LP sees, and what fixes the vectors is the octagon geometry (no short circuits), not the scheme.
- For (2,3,3) up to (3,4,4) the octagon has several vectors, but they span an affine space of smaller dimension than the LP solution set. So there are linear relations among the [i j k] that hold for every triple of that type in the Ree-Tits octagon and are not consequences of the intersection numbers and the Krein equation: 15 for (2,3,3), 12 for (2,3,4), 9 for (2,4,4), 14 for (3,3,3), 7 for (3,3,4), 3 for (3,4,4). Their exact integer form, and which of them hold in every GO(2,4) (14 of the 60, so far), are in `triple_relations_note.md`.
- For (4,4,4) the 44 distinct vectors span the full 10-dimensional LP solution set, so no further linear relation holds for all opposite triples.
- Distance types and group orbits differ a lot from (2,3,3) on, e.g. 86 T-orbits of type (4,4,4).

These are facts about the known octagon. A hypothetical other GO(2,4) need not have a big group, so the orbit refinement suggests constraints but proves nothing about uniqueness.

On the F4 embedding: 2F4(2) acts on the ambient metasymplectic space and preserves its relations, so on the Ree-Tits octagon any partition of triples by ambient relations is coarser than (or equal to) the 2F4(2)-orbit partition. The embedding can't split these orbits further; it could only matter for a hypothetical octagon, and Petit-Van Maldeghem do not show that every GO(2,4) embeds.

## Terwilliger algebra

T(x), generated by A and the dual idempotents E*_0, ..., E*_4 at a point x, lies in the centraliser algebra of the stabiliser of x. In orbit coordinates (`terwilliger_exact` in `triple_orbits.py`): T(x) has dimension 174, exact over Q (basis picked modulo a prime, then independence and closure under all six generators checked with python-flint). The centraliser algebra of 2F4(2)_x has dimension 258 and that of T_x 415. So T(x) is a proper subalgebra even of the 2F4(2)_x centraliser. Its Wedderburn decomposition (15 simple components, 4 non-thin irreducible modules) and the comparison with the centraliser are in `terwilliger_note.md`.

## Lines from the 'not opposite' graph (`line_cliques.py`, `ratio_bound_family.py`)

For a line L, the 27 points at distance <= 1 from L are pairwise at distance <= 3, so they are a coclique of the opposition graph Gamma_4.

- The eigenvalues of A_4 are 1024, -16, 16, -16, 64, and the ratio bound for cocliques of Gamma_4 is 1755 * 16 / (1024 + 16) = 27.
- CP-SAT: the largest clique of the distance-<=3 graph through a point has size 27, and with a cut excluding each of the 45 line neighbourhoods through that point the optimum drops to 26 (both OPTIMAL, a few seconds, no proof log). By transitivity the 27-cliques are exactly the 2925 line neighbourhoods, and L is recovered from its neighbourhood C as the points of C within distance 2 of all of C.
- For order (s, s^2) in general the A_4 eigenvalues are s^10, -s^4, s^4, -s^4, s^6, and the ratio bound is (s+1)^2(s^2-s+1) = (s+1)(st+1), the size of a line neighbourhood. So line neighbourhoods are maximum sets of pairwise non-opposite points in every generalized octagon of order (s,s^2); the bound uses only the parameters. That they are the only maximum sets is now proved for every finite GO of order (s,s^2), s >= 2 (`max_nonopposite_proof.md`).


Update 2026-11-22. The CP-SAT result is now backed by two more searches in `clique_certificate.py` (igraph 1.0.0 `largest_cliques`, and a separate colouring branch and bound written here): both return exactly the 45 line neighbourhoods through point 0, and T is checked to be transitive on points, so there are 2925. The Hoffman identity A_4 1_C = 16(1 - 1_C) holds on all 2925 line neighbourhoods, not just four.

There is now also a hand proof (`max_nonopposite_proof.md`) that in **any** finite generalized octagon of order (s,s^2), s >= 2, every set of (s+1)(s^3+1) pairwise non-opposite points is a line neighbourhood. It uses only the ratio-bound equality and the octagon axioms. The identities are checked symbolically in s (`max_set_general_checks.py`) and each lemma numerically on the Ree-Tits octagon O(2) (`max_set_checks.py`).
I could not find the source of the statement that the distance-1,2,3 graph of GO(2,4) has exactly 2925 cliques of size 27; it is checked here directly. Partial ovoids (pairwise opposite points) of GO(2,4) also have at most 27 points (Coolsaet-Van Maldeghem, J. Algebraic Combin. 12 (2000) 107-113), which is a different question.

## Not done

Quadruple counts, semidefinite (moment) strengthenings of the LP, and higher-dimensional coherent closure of the collinearity graph were suggested as next steps. None of them is done here.

## Summary

At the level of triple intersection numbers the vanishing Krein parameter reproduces the octagon geometry for one triple type and gives no further restriction, so it doesn't move the uniqueness question. The orbit computation shows where the actual octagon sits inside the LP: unique vectors for the ten small types, and extra linear relations for six of the remaining seven. Those relations, if they can be proved from the axioms, are the obvious thing to feed back into the LP.
