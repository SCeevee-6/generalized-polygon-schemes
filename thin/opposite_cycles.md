# Short cycles of the opposite-pair graph

In any generalized octagon of order $(2,4)$, each component of the
opposite-pair graph is either a grid component $K_{3,3}$ or has girth
at least eight. Every eight-cycle in a non-grid component alternates
between pairs containing a common line $A$ and pairs containing a
common line $C$, where $A$ and $C$ have line distance two.

## Short cycle test

The vertices of the graph are unordered pairs of opposite lines. Two
vertices are adjacent exactly when all four cross-distances are two in
the line graph. The three neighbors of a pair are the three pairs
from its three middle lines. A four-cycle gives a grid and the six pairs
of rows and columns exhaust one component. These facts are proved in
[the opposite-pair reduction](opposite_seed.md).

The three middle lines of any vertex are disjoint from the vertex
itself. Two distinct neighbors share one line and have their remaining
lines opposite. A triangle is therefore impossible: two neighbors of
a vertex share a line and so cannot be adjacent.

An induced cycle is specified by a sequence of pairs of line labels.
Each consecutive pair of vertices gives four distances of two.
Vertices two steps apart share exactly one label, and their unshared
lines are at distance four. These rules also apply across the join
where the cycle closes. Any two prescribed distances for the same
pair of lines must agree. If all four distances across a nonconsecutive
pair of vertices have been prescribed as two, there is a chord.

There is a further test on the line labels. Fix a label $A$. Any
opposite lines $B,B'$ at distance two from $A$ must use different
attachment points on $A$. Otherwise their paths to the attachment
point on $A$ give a walk of incidence length six between $B$ and $B'$,
contrary to opposition. Thus the graph of prescribed oppositions among
the labels at distance two from $A$ is three-colorable.

[The enumeration](opposite_cycles.py) fixes the initial three vertices
as $\{0,1\},\{2,3\},\{0,4\}$, which is possible by relabeling.
For each new vertex it chooses one of the two labels at distance two
back in the sequence, and either a previously used label or the next
unused label. This covers all possibilities because the two vertices
two steps apart share a line. The finite enumeration gives:

| Length | Distance traces | Traces passing the attachment test |
|---:|---:|---:|
| 5 | 0 | 0 |
| 6 | 2 | 0 |
| 7 | 0 | 0 |
| 8 | 6 | 2 |
| 9 | 2 | 2 |
| 10 | 38 | 0 |

The two length-six traces have a common label in each parity class.
At one of these labels, the opposite shared label and the three
unshared labels of its own parity class induce a $K_4$ in the
attachment test. Each surviving length-eight trace likewise has one
common label in each parity class. These are necessary label traces.
The enumeration does not assert that they are realized.

If a non-grid component had any cycle of length below eight, a shortest
one would be induced, so the test proves the girth bound. A chord in
an eight-cycle would also give a shorter cycle. The test applies to
every eight-cycle in a non-grid component.

## Cycles indexed by two lines

For any two lines $A,C$ at line distance two, let $J_{A,C}$ have the
following two vertex sets.

- Opposite pairs containing $A$ and having $C$ as a middle line.
- Opposite pairs containing $C$ and having $A$ as a middle line.

Keep the opposite-pair graph edges between these two sets. The sets
are disjoint since $A,C$ are at distance two. For an opposite pair
containing $A$ on the first side, exactly two of its three neighbors
contain $C$. Both are on the other side, since $A$ is one of their
middle lines. The reverse statement holds on the other side. Thus
$J_{A,C}$ is a disjoint union of cycles.

Each side has 64 vertices. To count the first side, start on the unique
four-edge incidence geodesic from $A$ to $C$. Continue for four edges
with $2\cdot4\cdot2\cdot4=64$ choices. The resulting paths are
nonbacktracking of length eight, so the other endpoint is opposite $A$.
Distinct continuations give distinct endpoints: two geodesics between
opposite lines cannot share the segment from $A$ to $C$, by uniqueness
of paths shorter than eight. Every pair on the first side arises in
this way.

Each cycle has length divisible by four. At a vertex on the side
containing $C$, its two neighbors on the side containing $A$ have form
$\{A,B\},\{A,B'\}$. The lines $A,B,B'$ are the three middle lines
of that vertex. Their attachment points on $C$ are distinct. Thus
as one walks round the cycle, the attachment points of successive
$B$ labels alternate between the two points of $C$ other than its
attachment toward $A$. The number of vertices on each side is even.

Every edge of the full opposite-pair graph belongs to exactly four
such cycle systems, since one chooses one line from each endpoint
pair. Every eight-cycle outside the grid components belongs to exactly
one system. Existence follows from the trace enumeration and uniqueness
from the intersection of the four pairs in each parity class. In
particular, one way to force a grid is to exclude
a simultaneous decomposition of all 128-vertex systems into cycles
of lengths $8,12,\ldots,128$ with the shared edge incidences just
described.

As a further count, let $N$ be the number of vertices of the
opposite-pair graph and $g$ its number of grid components. Every
non-grid vertex lies on 87 closed walks of length six, the number in
the cubic tree. Each grid vertex lies on 243. Hence
$\operatorname{tr}(J^6)=87N+936g$. This is a second trace criterion for
a grid, alongside the fourth trace in the original reduction.

## Nine-cycles and ten-cycles

A chord in a nine- or ten-cycle would split it into two cycles
whose lengths sum to at most twelve. In a non-grid component both
lengths would have to be at least eight, so these cycles are induced.
The enumerated absence of an induced ten-cycle therefore excludes
every ten-cycle.

The two nine-cycle traces name nine distinct lines and determine all
36 distances between them. Eighteen distances are two and eighteen
are four. Enumerate attachment point choices on each line up to a
permutation of its three points. Each trace has one assignment.
Name the intermediate line for each pair at distance two and apply
the short-path propagation described in
[the grid proof](grid_path_closure.md).

After identifying slots at distance zero, the result has 36 distinct
points and 21 distinct lines. Nine points are incident with one named
line, and the other 27 with two. The nine points of degree one are
mutually at point distance three.

The two neighbors of a vertex on the cycle name all three middle
lines, so their three unused tips are known. Complete the geodesics
between these tips. In [the edge labeling](opposite_edge_labels.md)
the two cycle edges at each vertex are forced to have label 4, and
the third edge has label 3. Thus every nine-cycle is an entire
component of the subgraph on edges labeled 4. Distinct nine-cycles
are vertex-disjoint. The computation checks implications of the
trace and leaves its existence open.

Run [the replay](opposite_nine_cycles.py) from the repository root
using `python thin/opposite_nine_cycles.py`.

## A local rank criterion

Fix $A,C$ at line distance two and take the two 64-vertex sides of
$J_{A,C}$. Its $64\times64$ biadjacency matrix $M$ has entry one for
the pairs $\{A,B\},\{C,D\}$ exactly when $B,D$ are at line distance two.
On a component of length $2m$, permuting the rows and columns gives
the block $I_m+P_m$, with $P_m$ the permutation matrix of an $m$-cycle.
Here $m$ is even. Over any field, the kernel is one-dimensional,
since its coordinates alternate in sign around the cycle. In
characteristic two the alternating signs both equal one.

The rank of $M$ is $64$ minus the number of cycles of $J_{A,C}$.
If there is no grid in this system, every cycle has length at least
eight and there are at most sixteen cycles. Therefore a rank below
48 for even one line pair $A,C$ forces a grid. Rank 32 means all
32 cycles have length four. This reduces another grid test to a
matrix on just 128 local lines. A universal rank bound below 48
remains to be proved.

The block counts and an exact rank-48 partial distance example are in
[Distance blocks in the local rank test](indexed_cycle_blocks.md).

## Cycles of lengths eleven through thirteen

In a non-grid component a chord of a cycle of length at most
thirteen would split it into two cycles whose lengths sum to
at most fifteen. This contradicts the girth bound. Hence the
induced trace enumeration also applies to every cycle of these
lengths.

The [long-cycle replay](opposite_long_cycles.py) identifies traces
under cyclic shifts, reversal and line relabelling. The first two
vertices use four distinct labels. Once their labels are fixed,
each subsequent vertex introduces at most one new label, which
fixes its renaming. For each trace passing the attachment test,
the replay enumerates projection assignments up to independent
permutations of the three points on each line. Each line has one
such assignment. The replay names the full intermediate line of
every distance-two connection and applies the three-line
distance propagation.

| Length | Traces after the attachment test | Surviving indexed traces | Other surviving traces | Force a grid at the first pair |
|---:|---:|---:|---:|---:|
| 11 | 8 | 0 | 1 | 0 |
| 12 | 6 | 1 | 3 | 2 |
| 13 | 71 | 0 | 1 | 0 |

Two of the six length-twelve traces that pass distance propagation
force a grid through the first opposite pair. Its whole component
is then $K_{3,3}$, which excludes a twelve-cycle through that pair.
The 128-vertex indexed systems also obey
[additional restrictions on open paths](indexed_cycle_paths.md).

The eleven-cycle and thirteen-cycle survivors also force the
edge labels. Complete the short geodesics between the three unused
middle-line tips at every vertex, as in the nine-cycle replay.
Every vertex of either survivor has incident labels $(3,4,4)$.
Up to a change of starting point and direction, the cycle labels
are as follows.

| Length | Edge labels in order |
|---:|---|
| 11 | $3,4,4,4,3,4,4,4,4,4,4$ |
| 13 | $3,4,4,4,3,4,3,4,4,4,3,4,4$ |

These are necessary traces. A surviving trace is a set of
distance domains and gives no ambient completion. In particular,
excluding eight-cycles would still leave the systems indexed by
two lines with possible cycle lengths $12,16,\ldots,128$.

The indexed systems through a fixed line are links of a single triangle
complex. [Compatibility and the Euler restriction](opposite_link_surfaces.md)
