# The non-grid components of the opposite-pair graph

An edge of the cubic graph on opposite-line pairs has an intrinsic label
in $\{1,3,4\}$. Grid components have label 1 throughout. In every other
component the edges labelled 4 form vertex-disjoint cycles. The statement
holds in an arbitrary generalized octagon of order $(2,4)$.

Use the three arms and middle lines $C_1,C_2,C_3$ of an opposite pair
$U=\{A,B\}$ from [opposite_seed.md](opposite_seed.md). Let $t_i$ be the
unused point on $C_i$. For the edge from $U$ to $V=\{C_j,C_k\}$ let $i$
be the omitted index. On the middle line of the unique geodesic from
$t_j$ to $t_k$, let $m$ be the third point. Label the edge by $d(t_i,m)$
in the point graph.

This definition is symmetric. The three middle lines for $V$ are $A,B$
and the middle line of $t_j,t_k$: the paths through $A,B$ use the two
other points on each of $C_j,C_k$, so the third path goes through $t_j,t_k$.
The unused points on $A,B$ are $a_i,b_i$, whose geodesic has middle
line $C_i$ and third point $t_i$. Reading the edge from $V$ therefore
interchanges $t_i$ and $m$.

## The three possible label multisets

Take the three incidence geodesics between the mutually distance-three
points $t_1,t_2,t_3$. Each has length six. Two from a common tip share
an initial segment and cannot meet again after separating, by girth 16.
If their union is a tree, the branch vertex is a line at distance three
from each tip. This is precisely the grid case. The tip omitted from a
pair is adjacent to the third point on its middle line, giving labels
$(1,1,1)$.

Otherwise the three paths contain a simple cycle. Removing each shared
initial segment from the two paths that use it leaves that cycle. Its length
is 18 minus twice the total length removed at the three tips. Girth leaves
two possibilities.

- **Length 18.** For each omitted tip, the two used points on the
  opposite pair's middle line are at incidence distance eight along the
  cycle. Those paths are geodesics. The unique closest point of the line
  is therefore its third point, at incidence distance six. All labels are 3.
- **Length 16.** Exactly one tip, say $t_1$, has a shared initial edge.
  It attaches to the cycle at a line $H$. For the pair $t_2,t_3$, the
  two used points on its middle line are seven steps from $H$ around the
  cycle. Prepending $t_1 H$ gives geodesics of length eight. The third
  point again has point distance three from $t_1$. For the omitted
  $t_2$, the path on the cycle from $t_2$ to $H$ has length five.
  From $H$ it takes one more edge to the nearer used point on the middle
  line of $t_1,t_3$ and two to that line. This seven-edge geodesic
  projects to a used point, so the unused point is at incidence distance
  eight from $t_2$. Exchanging $t_2,t_3$ gives labels $(3,4,4)$.

Thus each vertex has zero or two incident 4-edges. Finiteness and
symmetry of the label give the cycle assertion. Components entirely of
3-edges are still possible under this necessary condition.

## A five-line test on the 21-point seed

For a point profile $p$ on the 21 seed coordinates, the profiles on a
possible line through it can be enumerated exactly. Once the first two
profiles are chosen the third coordinate is forced everywhere: equal
values $a,a$ require $a-1$, and consecutive values require their maximum.
Retain only triples from the support being tested.

Associate to each incident triple the subset of coordinates at distances
one to three for which that line goes inward from $p$. Five such subsets
must partition all coordinates of $p$ at distances one to three. Empty
subsets may repeat. At a coordinate of value four every line goes inward
automatically. Allowing repeated line types only relaxes the test. The
verifier tries every cover by choosing the mask covering the first
uncovered coordinate and recursing.

Of the 756 candidate profiles under the no-seventh-point assumption,
three fail even using the full support. They are the three relabellings
of the profile with seven entries on each arm given below.

```text
3 2 2 3 3 1 3 | 4 3 3 4 4 2 4 | 4 3 3 4 4 2 4
```

Restricting to the 393 occupied profiles of `theta_no_extension.json`
rejects 117 named profiles in 14 orbits, with total multiplicity 249.
Thus that particular solution of the pair and cubic marginal equations
cannot be realized. Other multiplicities can change the available
lines. Neither this calculation nor the edge labels force a grid somewhere.

Run `python thin/theta_incident_lines.py` to verify the counts. [Verifier](theta_incident_lines.py)
