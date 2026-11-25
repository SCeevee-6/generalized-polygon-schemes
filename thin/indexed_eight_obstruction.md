# Every opposite-pair component needs a longer indexed cycle

**Theorem.** In a grid-free generalized octagon of order $(2,4)$, every
component of the cubic opposite-pair graph contains an indexed cycle of
length at least twelve. In each component at least $1/60$ of its indexed
edge incidences belong to such cycles. Globally at least **2,340 distinct
unordered indices** $\{X,C\}$ contain a cycle of length at least twelve.

## A finite system of frame relations

Write $(X,Y;C,D,E)$ for an ordering of an opposite line pair and its
three middle lines. There are twelve frames over each opposite pair.
Use the following four involutions.

| Move | Result |
|---|---|
| 0 | $(Y,X;C,D,E)$ |
| 1 | $(X,Z;C,D,F)$ |
| 2 | $(X,Y;C,E,D)$ |
| 3 | $(X,Y;D,C,E)$ |

Here $X,Y,Z$ are the middle lines of $\{C,D\}$, and $C,D,F$ the
middle lines of $\{X,Z\}$. Move 1 passes through $\{C,D\}$, so each
move stays in the component of $\{X,Y\}$. Nonconsecutive moves commute.
Moves 0 and 1 act as transpositions on $X,Y,Z$ with $C,D$ fixed.
Moves 2 and 3 act as transpositions on $C,D,E$. Alternating either
of these adjacent pairs for six steps returns to the starting frame.

For moves 1 and 2, fix $X,C$. A frame specifies the edge from $\{X,Y\}$
to $\{C,D\}$ of $J_{X,C}$. The moves change alternating ends of that
edge. Eight alternating steps return to the frame exactly when the
edge is on an eight-cycle, since the possible grid-free lengths are
$8,12,16,\ldots$. [Indexed cycles](opposite_cycles.md)

Take as reference vertices the 24 vectors $v$ with two entries $\pm1$
and two zeros. The 24 facet normals $w$ are the sixteen sign vectors and
the eight vectors $\pm2e_i$. Incidence means $v\cdot w=2$. A reference
frame $(x,y;c,d,e)$ is a pair of distinct vertices with exactly three
common facets, in the order $c,d,e$. The four moves above act on these
1,152 frames. In particular $c,d$ have a unique third common vertex $z$
besides $x,y$.

Fix a starting reference frame and a spanning tree of its flag graph.
From any actual starting frame, the tree words assign an actual frame to
every reference frame. An edge is consistent when applying its move
takes one assigned frame to the other. Tree edges are consistent. Around
a two-color polygon with just one edge still to check, the alternating
relation at a fixed flag of the polygon forces consistency of the last
edge. Follow the two known paths from that flag to the ends of the missing
edge. For colors 1 and 2 this requires the assigned base flag to lie on
an indexed eight-cycle. All other polygons use unconditional relations.

The verifier builds the tree by distance from frame 0 and fills one
missing edge at a time, taking unconditional polygons first. It fills
all edges using exactly **97** polygons of colors 1 and 2. Consequently,
if all 97 fixed tree words take an actual starting frame to eight-cycles,
the assignment respects all four moves.

## What the assigned frames would force

The moves fixing the first vertex $x$ connect all frames with that first
entry. Likewise the moves fixing the first facet $c$ connect all frames
with that third entry. Both are finite connectivity checks. Each vertex
and each facet therefore receives a well-defined ambient line. Swaps
identify the other entries with the corresponding labels. Incident
vertex–facet labels are at line distance two. Vertices sharing three
facets are opposite, and facets sharing three vertices are opposite.
In the latter case the three vertices label their three middle lines by
move 1. We may interchange vertices and facets below.

At a central line the six incident labels form three pairs. Labels from
different pairs are opposite, so they attach to distinct points of the
central line. This complete tripartite graph has only one partition
into three independent sets. Its classes give the three attachment
points, called ports below. At $v$ paired facets satisfy $w+w'=2v$.
At $w$ paired vertices satisfy $v+v'=w$.

To see that the two connector lines at each port are distinct, fix a
central label. Its eight neighboring face triples use all choices of one
label from each port. Two triples differing in one entry have opposite
partner lines, since those partners are two middle lines of the retained
pair. The changed labels must differ: otherwise two different
opposite-pair vertices would have the same middle triple, forcing a
grid. Thus all six labels and all eight triples are distinct. Every
edge of these faces already has both incident faces. They form an
eight-face rooted piece. By the eight-face attachment calculation, the
two intermediate lines in each color are distinct. This applies at
central labels on either side. [Eight-face result](opposite_link_surfaces.md)

Form a graph on the 144 reference ports. Join the three ports of each
central line in a triangle. For each incident vertex–facet pair join
its two facing ports. Each port has two triangle neighbors and two
connector neighbors. A shortest path of length at most four lifts to
a nonbacktracking ambient incidence walk of length at most eight. Two
successive steps in one triangle would shorten the path, and the two
connector lines at any port are distinct from each other and from the
central line. By girth sixteen the walk is geodesic. Thus every such
graph distance is its actual point distance, even if distant reference
labels have coincident images.

## Two configurations of four points

In the verifier's numbering take the cyclic orders
\[
 P=(0,15,69,54),\qquad Q=(104,107,113,110).
\]
Breadth-first search verifies the following facts.

- Consecutive entries of either group have distance two. Diametrically
  opposite entries have distance four.
- The midpoint of a consecutive pair has distance four from every
  point of the other group.
- All four graph neighbors of each $p$ have distance four from every
  $q$ in the other group.

The two other points of the central line through $p$ thus have distance
four from $q$. Projection to that line gives $d(p,q)=3$. On each connector
through $p$ one named neighbor has distance four from $q$. The three
distances on a line are a permutation of $(h,h+1,h+1)$, so neither
connector goes inward. Of the five lines through $p$, only two remain
for the four cross geodesics.

From $p$, consecutive targets $q,r$ must use different first lines.
Their three incidence geodesics have lengths $6,6,4$. If the two from
$p$ shared a first edge, their union with the third path would have
fewer than sixteen edges and be a tree. Its median would be the midpoint
of $qr$, at incidence distance four from $p$, contradicting its verified
point distance four. Thus the two first lines alternate. Diametrically
opposite targets use the same first line, in both directions.

## The four remaining geodesics

Choose diametrically opposite $p,p'$ in $P$ and $q,q'$ in $Q$. Denote
the common first lines at the four points by $L,L',M,M'$ respectively.
Each cross pair of lines is at incidence distance four along the
six-edge geodesic of the corresponding points.

Suppose the routes from $M$ toward $L,L'$ attach at different points
on $M$. Then the eight-edge route from $L$ to $L'$ through $M$ is
nonbacktracking, so $L,L'$ are opposite. The route through $M'$ also
has length eight and is geodesic. The two arms through $M,M'$ are
distinct and use the two points other than $p$ on $L$ and the two
other than $p'$ on $L'$. An arm attaching at $p$ on $L$, or at $p'$
on $L'$, would backtrack on the corresponding cross geodesic. The third arm must
then match $p$ to $p'$ at point distance three, contradicting opposition.

Interchanging the two sides gives a common attachment point $u$ on $L$
toward $M,M'$ and a common point $u'$ on $L'$. The routes from $q$ and
$q'$ to $u$ have incidence length four. Together they form a geodesic
since $q,q'$ are opposite. The same holds for $u'$. The two geodesics
start on $M$ at $q$, so their remaining seven edges are the unique
geodesic from $M$ to $q'$. Their midpoints coincide. This point would
be collinear with both $p$ and $p'$, again contradicting opposition.
This excludes the assignment respecting all four moves.

## Counting long cycles

In a component with $n$ opposite pairs, each of the 97 fixed words
permutes the $12n$ frames. From every starting frame some word must
land on an indexed cycle of length at least twelve. The union of the
97 inverse images therefore gives at least $12n/97$ long-cycle frames.

For an unordered index $\{X,C\}$ each edge gives two frames, one with
$X$ first and one with $C$ first. The same fraction bound holds for
indexed edge incidences. In the whole octagon there are
$N=2925\cdot2048/2=2,995,200$ opposite pairs. The sum of the lengths
of the long cycles is at least $6N/97$, and each index accounts for at
most 128 edges. This yields
\[
\left\lceil\frac{6N}{97\cdot128}\right\rceil=1448.
\]

## Strengthening the count with rooted surfaces

A frame $(X,Y;C,D,E)$ determines a rooted face at $X$ with partner $Y$.
Call the frame large when this face belongs to a rooted piece with more
than eight faces. In each opposite-pair component at least **one frame in
twenty is large**.

To prove this, use the same reference graph with a different assignment
tree. Choose the twenty vertex labels other than $0,5,19,22$ in the
verifier's coordinate order. For each chosen $x$ first take a spanning
tree in the residue fixing $x$ (moves 1, 2 and 3). Extend this forest to
a spanning tree of the whole flag graph. Use the tree word to one fixed
flag in each of the twenty residues.

Suppose these twenty words all take the starting frame to eight-face
pieces. Such a piece has the eight triples of the octahedron and the
three permitted moves on its 48 flags. There is an isomorphism from the
reference residue to the rooted piece taking the designated flag to its
assigned flag by matching its ordered face and then the opposite label in
each of the three colors. Thus **every internal edge** in each chosen residue
is consistent with the tree assignment. Starting with these edges and
the tree, the verifier fills all remaining edges using only unconditional
two-color relations. This would give the frame map already excluded.
The twenty fixed tree words permute the frames in a component, proving
the fraction bound as before.

Consider one rooted piece with $f>8$ faces and $v$ split vertices. The
degree of a split vertex is half its indexed cycle length. The links are
simple even cycles, every degree is at least four, and $v-f/2\le2$ by
Euler characteristic. Write $h$ for the number of vertices of degree
greater than four and $L$ for their degree sum. Then
\[
 L=3f-4(v-h)\ge f-8+4h.
\]
There must be at least two such vertices. If $h=0$ the displayed bound
gives $f\le8$. If $h=1$ and the unique degree is $d$, each of the other
two colors has $f/4$ vertices. The link at that vertex uses $d/2$
distinct vertices of each of those colors, so $d\le f/2$. Euler gives
$d\ge f-4$, again yielding $f\le8$. Thus $L\ge f$ for every large
piece. At least a third of its $3f$ face–vertex incidences correspond to
long cycles. [Rooted surfaces](opposite_link_surfaces.md)

All faces in a rooted piece stay in the same opposite-pair component, since
crossing a face edge changes partners via their common neighbor. Each
face–vertex incidence at $X$ occurs in two frames with first endpoint
$X$ and the given middle line first. Hence the two bounds multiply. The
long-cycle fraction in each component is at least $1/(20\cdot3)=1/60$.
Globally this gives at least **299,520 indexed edge incidences** on long
cycles, hence at least **2,340 indices** of capacity 128.

There are also at least $2N/20=299,520$ rooted faces in large pieces,
spread across at least 147 root lines (2,048 faces per root). Writing
$k_X$ for the number of pieces at $X$, we obtain
\[
 \sum_X(256-k_X)=\sum_{\text{pieces}}(f-8)/8\ge12,480.
\]
Here each large piece has $f\ge12$, hence $(f-8)/8\ge f/24$.
This bound also applies when every failed pair has the 18-cycle shape.

Run `python thin/indexed_eight_obstruction.py` from the repository root.
The [verifier](indexed_eight_obstruction.py) uses the Python standard
library. The earlier eight-face calculation can be replayed with
`python thin/opposite_link_surfaces.py`. Compatibility of the systems
containing longer cycles remains open.

## Dependence on the 16-cycle shape and the equality case

In a cubic component let $p$ be the fraction of its opposite pairs of
16-cycle shape. Let $b$ be the fraction of its frames for which the rooted
piece is large at **both** endpoints of the underlying pair. Both endpoints
of a 16-shape pair are large, because all faces of an eight-face piece have
18-cycle shape. Thus $b\ge p$. [Eight-face calculation](opposite_link_surfaces.md)

The verifier chooses the assignment tree so the twenty tested residues
are joined by nineteen tree edges of move 0. Contracting their internal
trees gives a tree on the twenty tests. For each starting frame mark the
tests landing in large pieces. There is at least one mark. The subgraph
induced by the marks is a forest, so the number of marks is at least one
plus the number of its edges. Choose the starting frame uniformly in
the component. Each individual tree word acts bijectively. On a join,
the two flags at its ends are assigned by exchanging the endpoints of
the same opposite pair. Within either residue its internal tree path
preserves the rooted piece, even when that piece is large. Thus each
join has two marked ends with probability exactly $b$.

If $a$ is the fraction of large frames and $t$ the long-cycle fraction,
taking expectations and applying the rooted degree bound gives
\[
 20a\ge 1+19b,\qquad t\ge a/3\ge(1+19b)/60.
\]
In particular, if $N_{16}$ is the total number of 16-shape pairs in a
grid-free octagon, at least
\[
 \left\lceil 2340+19N_{16}/1280\right\rceil
\]
distinct indices contain long cycles. This follows by summing $6nt$ over
components of size $n$ and dividing by the capacity 128 of an index.

Suppose exactly 2,340 indices contain long cycles. Every inequality
giving this bound is then an equality. Each of those indices consists
entirely of long cycles, $b=0$ in every component, and every large rooted
piece satisfies $L=f$. The Euler argument above gives exactly two vertices
of degree greater than four in each such piece.

In fact the piece is a suspension of an even cycle. Take a high-degree
vertex $u$. All but at most one vertex on its link have degree four.
For any such degree-four neighbor $x$, the link at $x$ is a square.
There is a unique second neighbor of $x$ in the color of $u$. Two
successive degree-four neighbors on the link of $u$ have the same second
neighbor: the other triangle on their common edge supplies it. Removing
at most one vertex from a cycle leaves it connected. Hence these second
neighbors are one vertex $u'$. Every edge of the link has a degree-four
endpoint, so its two incident faces use $u,u'$. These triangles already
form a closed suspension. Closure under face adjacency exhausts the piece.

Every face of this suspension has exactly one high-degree vertex, its
pole. Thus at a cubic vertex $\{A,B\}$ a large piece at $A$ supplies
exactly one long index through $A$, and $b=0$ makes the piece at $B$
small. If both pieces are small there is no long index there. Hence
**the long indexed cycles are vertex-disjoint in the cubic graph**.
Also at a fixed root all its long indices must attach at the same point
on the root line. Two middle lines in different attachment colors occur
together in a face at that root. If both indices were long, purity of
the indices would give two high-degree vertices on that face.

These are necessary conditions for equality. Grid-free systems with
more long indices, and the existence of the extremal configuration
itself, remain open.
