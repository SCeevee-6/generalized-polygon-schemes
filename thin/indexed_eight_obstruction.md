# Every opposite-pair component needs a longer indexed cycle

**Theorem.** In a grid-free generalized octagon of order $(2,4)$, every
component of the cubic opposite-pair graph contains an indexed cycle of
length at least twelve. In each component at least $1/97$ of its indexed
edge incidences belong to such cycles. Globally at least **1,448 distinct
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
other than $p'$ on $L'$. An arm using one of those endpoints there
would backtrack on the corresponding cross geodesic. The third arm must
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

Run `python thin/indexed_eight_obstruction.py` from the repository root.
The [verifier](indexed_eight_obstruction.py) uses the Python standard
library. The earlier eight-face calculation can be replayed with
`python thin/opposite_link_surfaces.py`. Compatibility of the systems
containing longer cycles remains open.
