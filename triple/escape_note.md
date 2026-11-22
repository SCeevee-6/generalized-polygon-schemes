# Non-opposite sets that lie in no line neighbourhood

`max_nonopposite_proof.md` shows that the largest sets of pairwise non-opposite points in a GO(s,t) are the line neighbourhoods N[L] (points at distance ≤ 1 from L). The natural next question is a stability one: how large can a non-opposite set be if it is contained in no N[L]? Call such a set *escaping*.

The earlier clique searches (`line_cliques.py`, `clique_certificate.py`) don't answer this. They only exclude choosing all of some N[L], so N[L] minus a point is still a feasible 26-set there. Here the condition is containment.

Everything below is for the Ree-Tits octagon O(2), order (2,4), and its dual, order (4,2). It is computational and says nothing yet about other octagons of these orders or about O(q) for q ≥ 8.

## Results

**O(2), points.** `escape_bk.py --geom points --min 15` enumerates every inclusion-maximal non-opposite set of size ≥ 15 that contains point 0: 12800 of size 15, 45 of size 27, and none of any other size (1262329 search nodes, about 25 s). The 45 are the line neighbourhoods through 0. T is transitive on points, so:

- every non-opposite set of at least 16 points of O(2) lies in a line neighbourhood;
- there are exactly 12800·1755/15 = 1497600 escaping 15-sets, all of them maximal.

Each escaping 15-set has the same shape. It contains exactly one full line L, and its other 12 points are at distance 2 from L, pairwise non-collinear. Counting how many of the 12 project onto each point of L gives two profiles:

| projection counts on L | sets per line L | in O(2) | T-orbits | 2F4(2)-orbits | largest overlap with a line neighbourhood |
|---|---:|---:|---|---|---:|
| (4,4,4) | 128 | 374400 | 187200 + 187200 | 374400 | 7 |
| (0,4,8) | 384 | 1123200 | 1123200 | 1123200 | 11 |

The orbits are orbits of the stabiliser of L (Schreier generators), times the 2925 lines; the 2F4(2) column is the same run with `--gens atlas` (`out_escape/escape_bk_points_min15_T2.json`). So the outer automorphism swaps the two (4,4,4) orbits of T. In the (0,4,8) sets, the 8 points behind one point p of L together with L lie in N[M] for a line M ≠ L on p, which is where the 11 comes from.

**Dual of O(2), lines.** `escape_bk.py --geom lines --min 25` does the same for the 2925 lines, with "opposite" meaning distance 4 in the line graph and the neighbourhoods N*[p] = lines meeting the closed neighbourhood of a point p (45 lines each). Through line 0 there are 128 maximal families of size 25, 27 of size 45 (the N*[p] containing line 0), and nothing else of size ≥ 25. So

- every family of at least 26 pairwise non-opposite lines of O(2) lies in some N*[p];
- the escaping families of size 25 number 128·2925/25 = 14976 and form one T-orbit (`escape_dual_orbit.py`).

These 25 lines are pairwise at distance 3 in the line graph, so no two meet and no two meet a common line; they cover 75 distinct points, and no N*[p] contains more than 5 of them. The stabiliser of a family has order |T|/14976 = 1200 and is transitive on its 25 lines. A subgroup of order 1200 of 2F4(2)' can only lie in the maximal subgroup 5²:4A4 (the other maximal subgroups have orders not divisible by 1200), so the stabiliser is 5²:4A4. In 50 random elements of order 5 (seeded), each fixes no line, fixes exactly one family, and permutes its lines in five 5-cycles. So the Sylow 5-subgroup of the stabiliser acts regularly on the family.

The two geometries behave differently: the extremal escaping sets in O(2) contain a line, while those in the dual contain no two collinear points.

## A local model for the O(2) sets

`escape_local.py` looks only at sets that contain a full line L and whose other points are at distance 2 from L. For p ∈ L and a line M_p ≠ L on p, let X(p, M_p) be the s²t = 16 points x at distance 2 from L with geodesic x, q, p and q ∈ M_p.

- Two points with the same projection p are non-opposite iff their middle points lie on one line through p (otherwise the path between them is non-backtracking of length 4). So each projection class lies in a single X(p, M_p), and the points inside one X are pairwise non-opposite.
- Points with different projections are non-opposite iff they are at distance 3. Call this compatible.
- Such a set escapes every line neighbourhood iff at least two projection classes are nonempty.

So the escaping sets of this shape are L together with a set that is a clique for compatibility across three classes of 16. The script checks, for all 4³ = 64 choices of (M_p):

- each x ∈ X(p, M_p) is compatible with exactly one point of N' \ {q'} for every q' ∈ M_p' \ {p'} and every line N' ≠ M_p' on q'. This one is a general fact: q' is opposite x, and every line on a point opposite x has a unique point at distance 3 from x. So x is compatible with st = 8 points of X(p', M_p');
- compatibility between two X's has twin classes of size 4 (points with the same compatible set), and the graph between twin classes is an 8-cycle, in all 384 ordered cases;
- the largest compatible selection with at least two nonempty classes has 12 points for every one of the 64 choices. The 12-point selections have profiles (4,4,4) (128 of them) and (0,4,8) (1536, since the line over the empty class is free).

For two classes the 12 follows from the 8-cycle: one twin class is compatible with its two neighbouring classes (4 + 8), two classes at distance 2 share one neighbour (8 + 4), and nothing larger is complete. The three-class maximum is only checked by computer.

What's not explained is why every escaping set of size 15 must have this shape in the first place.

## Open

- A proof, for every GO(2,4), that an escaping set has at most 15 points, or at least that a 15-point one contains a line.
- The analogue for O(8), whose 1.2·10⁹ points rule out this kind of search. The local model above only needs the geometry near one line, so it may be the right thing to compute there.
- A formula in (s,t). The values 15 at (2,4) and 25 at (4,2) are the only data so far, and the extremal sets have different shapes in the two cases.

## Files

- `escape_common.py`: builds O(2) or its dual from generators (`--gens atlas` fetches c, d of 2F4(2); `gap:DIR` reads TF42G1-p1755B0.g1/.g2 for T).
- `escape_bk.py`: Bron-Kerbosch with pivoting and a colouring bound on bitsets (numba). The search gets no information about lines; the sets are compared with the line neighbourhoods afterwards.
- `escape_local.py`, `escape_dual_orbit.py`: the local model and the dual orbit checks.
- `out_escape/`: logs, JSON summaries and the full lists of sets (`N` = in a line neighbourhood, `E` = escaping). The committed runs use the T generators.

A separate C++ Bron-Kerbosch, run on a graph built from the edge orbit of {0,72} under c, d rather than from the point stabiliser, gave the same counts 12800 and 45 and the same profiles.
