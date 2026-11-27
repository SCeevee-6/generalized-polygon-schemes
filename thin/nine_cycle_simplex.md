# The nine distinguished points of a cubic nine-cycle

Suppose the opposite-pair graph of a generalized octagon of order $(2,4)$
has a nine-cycle. Its seed has 36 points on 21 lines. The nine points
incident with just one of these lines have mutual point distance three.
Write $S$ for this nine-set. **No triple of $S$ has a center line.** For
each $x\in S$ there is a unique line through $x$ unused by the eight
geodesics to $S\setminus\{x\}$. These nine unused lines are pairwise
opposite. [Nine-cycle trace and seed](opposite_cycles.md#nine-cycles-and-ten-cycles)

## Excluding centers

For $x,y\in S$ distinct, write $f(x,y)$ for the first point after
$x$ on the unique point geodesic to $y$. If $x,y,z$ are centered,
$f(x,y)=f(x,z)$ is the projection of $x$ on their center line. There
are two further such equalities at $y,z$.

Start with the two nine-cycle traces and their unique attachment choices.
Identify slots forced equal in the complete seed metric. Complete the
unique geodesics of all 36 unordered tip pairs, including the third point
of each segment line. All fresh names may coincide when the distances
allow it. For each of the $\binom93=84$ triples, impose the three necessary
equalities and propagate the full-line distance patterns.

| Result per trace | Triples |
|---|---:|
| A necessary equality is already impossible | 54 |
| Joint equalities give a contradiction in propagation | 27 |
| Survive the equality test | 3 |

In **each** trace the three surviving triples partition $S$. Every tip
therefore lies in at most one centered triple.

In any nine-set of mutually distance-three points, every tip lies in an
*even* number of centered triples. For a line $L=\{x,p,q\}$ through a
tip, write $a,b$ for the numbers of other tips reached through $p,q$
on that line. The nine least-eigenspace vectors sum to zero. Pairing
with $e_p$ gives $-8+4a-2b+(8-a-b)=0$, hence $a=b$. Each is at most
two since at most two mutually distance-three points can be at distance
two from $p$. Two targets reached through $p$ form a centered triple
with $x$. Indeed $p$ is the third point on the middle line of the
geodesic between those two targets, and $x$ reaches that line through
$p$ at incidence distance three. Conversely a center for the two
targets and $x$ gives the common first point. Thus $L$ contributes
either zero or two centered triples at $x$. Every centered triple at
$x$ is counted once. The degree bound of one now excludes all centers.
[Cosine identity and middle-line lemma](opposite_seed.md#one-centered-transversal-forces-thirty-lines)

Each of the eight targets from $x$ consequently uses a different first
point. They use four lines in pairs, leaving a fifth $L_x$. From $L_x$
to $L_y$ the route through $x,y$ has incidence length eight and does
not backtrack, proving opposition.

## Counts in the whole octagon

Let $n_j$ count the elements of $S$ at point distance $j$ from an ambient
point. The point types and their multiplicities are forced.

| Point | $(n_0,n_1,n_2,n_3,n_4)$ | Number |
|---|---|---:|
| Tip | $(1,0,0,8,0)$ | 9 |
| Neighbor on a used line | $(0,1,1,1,6)$ | 72 |
| Neighbor on the unused line | $(0,1,0,0,8)$ | 18 |
| At distance two from two tips and adjacent to none | $(0,0,2,5,2)$ | 36 |
| At distance two from one tip and adjacent to none | $(0,0,1,4,4)$ | 576 |
| At distance at least three from every tip | $(0,0,0,3,6)$ | 1,044 |

The cosine identity gives $16n_0-8n_1+4n_2-2n_3+n_4=0$. Each pair
of tips has a unique common distance-two point and $n_2\le2$, giving
36 different such points. None is a neighbor of a tip by the used-line
counts. The total distance-two count is $9\cdot80$, which gives 576
in the fifth row after subtracting $72+2\cdot36$. Subtraction from
1,755 gives the last row.

## Remaining compatibility

The original line at each tip pairs two targets. Those pairs give the
two cyclic neighbors in an order on the nine tips. The verifier also
checks an assignment at all nine tips at once. In that cyclic order
pair the target offsets $(1,8),(2,3),(4,5),(6,7)$ at every tip.
The first neighbor toward each target has distance two to that target,
three to its partner, and four to all other targets. Impose the resulting
distances and insert the nine full unused lines. The joint assignment
passes three-line propagation in both traces. Existence of a nine-cycle
in an octagon remains open.

Run `python thin/nine_cycle_simplex.py` from the repository root.
[Verifier](nine_cycle_simplex.py)
