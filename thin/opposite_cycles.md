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
