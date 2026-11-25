# Compatibility of indexed systems at a fixed line

Fix a line $A$ in a grid-free generalized octagon of order $(2,4)$.
The 96 lines at line distance two from $A$ are divided by their
attachment point on $A$ into three classes of size 32. For each of
the 2048 opposite partners $B$ of $A$, make a triangle whose vertices
are the three middle lines of $\{A,B\}$. Its three vertices have
different classes.

These triangles are distinct. Two identical triangles would give
a four-cycle in the opposite-pair graph. Every pair of vertices
from different classes occurs in exactly two triangles. Indeed
the two lines are opposite by the length-eight path through $A$.
Their other two middle lines are precisely the two possible
partners of $A$. These facts are also used in
[the grid criterion](opposite_seed.md).

Glue the two triangles at each edge. The link at a line $C$
is the graph whose vertices are the other middle lines occurring
with $C$ and whose edges come from opposite partners $B$ using $C$.
It is a union of cycles. Subdivide every edge of this link once,
label its new vertex $\{A,B\}$ and label each old vertex $D$ by
$\{C,D\}$. The result is exactly $J_{A,C}$ from
[the indexed cycle construction](opposite_cycles.md).
The 96 indexed systems through $A$ are therefore links of one
joint triangle complex.

Split each vertex into one vertex for each component of its link.
The result is a union of closed triangulated surfaces. Each edge
still has two incident triangles, and each vertex now has a single
cycle as its link. The connected components are precisely the
components of the graph on triangles joined across their edges.
In fact a vertex link is connected after splitting, so all triangles
at one vertex belong to one component of that graph.

## An Euler restriction on joint cycle counts

For a connected piece with $f$ triangles and $v$ split vertices,
the edge count is $3f/2$. Each split vertex belongs to a link cycle
of length at least four, since an indexed cycle has length at least
eight before removing its subdivision. Its Euler characteristic
satisfies
\[
v-f/2\le 2.
\]
The piece contains at least two vertices of each of the three
classes. Starting from any face, cross its edge opposite a given
class. The new face has a different line of that class because
there are no repeated triangles. This gives $v\ge6$ and $f\ge8$.
In particular the 2048 opposite partners split into at most
256 components.

Let $c(A,C)$ be the number of cycles of $J_{A,C}$ and $k_A$
the number of connected pieces. The exact sum of Euler
characteristics of the surfaces is
\[
\sum_{C:\,d_{\rm line}(A,C)=2} c(A,C)-1024.
\]
In particular $\sum_C c(A,C)\le 1024+2k_A$. Combining this with
the [local matrix rank formula](opposite_cycles.md) gives
\[
\sum_C \operatorname{rank} M_{A,C}\ge 5120-2k_A
\]
over every field. For example, if the triangle adjacency graph
at $A$ is connected, some indexed system through $A$ has rank
at least 54. If all indexed systems through $A$ consist of
eight-cycles in the grid-free case, then all 256 pieces have
exactly eight faces.

## The case of eight faces

An eight-face piece has exactly six split vertices, two in each
class. The argument using adjacent faces shows that the two labels
in a class are distinct even before splitting. The faces must be
all eight triples on these three pairs of labels. The resulting
surface is the boundary of an octahedron, and its triangle
adjacency graph is a cube.

Write its six middle lines as $C_{i,s}$ and its eight opposite
partners as $B_{s_0,s_1,s_2}$. For any fixed $C_{i,s}$ the link is
a square, which gives an indexed eight-cycle. Choose a color
different from $i$. Its two middle lines occur on the same side
and in the same branch of this indexed eight-cycle. The
[different-ancestor restriction](indexed_cycle_paths.md)
shows that the two middle lines in that class have distinct
ancestor points on the paths from $A$. Thus they have either
different intermediate lines or the same intermediate line
and different points on it. This applies to all three classes.

The full-line replay leaves at most one class in which the two
middle lines have the same intermediate line. There are eight
labelled choices of the set of classes sharing an intermediate
line. Four pass the replay, giving two types up to permuting
classes. The order of attachments is forced up to relabelling:
at each $C_{i,s}$ the attachment toward $A$ uses one point.
The four adjacent partners form a square under changes of a
single remaining coordinate. Two such partners are opposite,
since they are both middle lines of the pair of middle lines
they share. Hence the other two attachment points on $C_{i,s}$
alternate on the square. At $B_{s_0,s_1,s_2}$ its three middle
lines use distinct points. The verifier imposes these incidences,
the three-line distance constraints and the forced oppositions.
The diagonal and failed-pair checks below exclude the shared case.

## Surface components with at most fourteen faces

Every component has an even number of faces. Write its color sizes after
splitting as $v_0,v_1,v_2$. The graph of edges between any two colors is
simple. Splitting a vertex cannot create parallel edges because there was
only one edge between the two original labels. The links are simple even
cycles. In a given color their lengths sum to $f$, hence $2 \le v_i \le
\lfloor f/4\rfloor$. Ten faces would give just two vertices in each color,
permitting at most eight distinct triples.

If one color has size two, every edge between the other two colors supports
both possible faces. Its graph is the link of either pole and is a single
even cycle. The component is consequently the suspension of that cycle and
has a multiple of four faces. With twelve faces there must be a two-vertex
color, since nine split vertices would give Euler characteristic three.
Thus the only twelve-face surface is the suspension of a hexagon.

With fourteen faces all three colors have size three. In each color the
three link lengths are $4,4,6$. Take the links of the two degree-four
vertices in one color. In the bipartite graph on the other colors these
are two rectangles. Their symmetric difference is the link of the third
vertex, since each edge supports two faces. The rectangles share exactly
one edge, as their symmetric difference has length six. Labelling that
edge by $(0,0)$ gives the following face set, with the degree-six vertex
of each color labelled zero:
\[
  (\{0,1\}^3\setminus\{000\})\;\cup\;
  (\{0,2\}^3\setminus\{000\}).
\]
It consists of two octahedra joined along the boundary of a removed face.

In each of these three surfaces every two vertices of the same color
appear on opposite sides of some face edge. This follows directly from
the displayed face set and the square or hexagonal ring of the suspensions.
Their original line labels must be distinct, since the two faces at that
edge are distinct triples. Their ancestor points must also be distinct.
Indeed the [eight-output repeat test](opposite_leaf_maps.md) applies to
the other two colors on that edge. The two faces give its two places and
a common ancestor would put them in the same block. Thus in a class
of three, at most one pair of labels can use the same intermediate line,
which has only two available ancestor points.

## Replaying the attachments

Two adjacent faces have opposite partner lines. Their common middle lines
are an opposite pair whose three middle lines are the root and the two
partners. At each vertex of the surface, the face partners therefore use
the two outward points of its middle line alternately around the link.
Choosing their names fixes the assignments up to swapping those two
points. At each partner the three middle lines use distinct points.

For a face with partner $B$, call the other partners across its three
edges $B_0,B_1,B_2$. These are precisely the third middle lines recorded
for $\{A,B\}$ in [the boundary calculation](opposite_boundary.md). Their
mutual line distances are either $(3,3,3)$ or $(2,3,3)$. The verifier
tries every ordered line square of distance two or three for each such
pair and applies [three-line propagation](grid_path_closure.md). A square
is removed only when that choice gives an empty domain. Unspecified slots
may coincide.

There is a second check on a chosen face. The full failed-pair metric on
its 36 slots is already named. In fact the geodesic between the two tips
of colors $j,k$ runs through the partner across the edge retaining $j,k$.
Its two internal points are on that partner and its three line spares are
also named. The verifier tries all 32 [labelled failed metrics](failed_pair_types.md)
on these slots. The following counts refer to fixed labels and, in the
last column, to a single designated face.

| Faces | Attachment choices | After diagonal test | After first-face test |
|---:|---:|---:|---:|
| 8 | 8 | 1 | 1 |
| 12 | 32 | 2 | 1 |
| 14 | 64 | 8 | 8 (all eight pass simultaneous face checks) |

For eight and twelve faces the surviving type uses distinct intermediate
lines in every color. In the octahedron all twelve diagonals have line
distance three. Consequently **all eight faces have the 18-cycle shape**.
For the hexagonal suspension the extra survivor after the diagonal test
has just the two poles sharing an intermediate line. None of the 32
first-face metrics fits that case. In the all-distinct case eight
18-cycle and six 16-cycle metrics pass the single-face test.

In a fourteen-face piece only labels 1 and 2 of a color may share an
intermediate line. In particular the only possible edges across which the
changed label keeps its intermediate line are the three edges with both
retained labels zero. Permuting colors reduces the eight possibilities to
four, distinguished by how many of these three edges have sharing.
For each of them, and for both smaller surviving surfaces, an explicit
assignment of failed metrics to **all** faces passes the diagonal and
three-line checks. These assignments establish local compatibility.
Whether any occurs in an octagon is open.

Run `python thin/opposite_link_surfaces.py` from the repository root.
[Verifier and joint assignments](opposite_link_surfaces.py)

## A global count of shared intermediate lines

In a grid-free octagon let $N_{16}$ be the number of unordered opposite
line pairs with the 16-cycle shape. Over all root lines and all components,
the number of face edges across which the changed color keeps its
intermediate line is **exactly** $N_{16}$.

Indeed an edge at root $A$ retains two middle lines $D,E$. Put
$V=\{D,E\}$. Its three middle lines are $A,B,B'$, where the two faces
at the edge have partners $B,B'$. Let $u$ be the unused point on $A$,
and $b,b'$ the unused points on $B,B'$. The geodesics from $u$ to $b,b'$
use the third arms of $\{A,B\}$ and $\{A,B'\}$, respectively. Their first
lines are precisely the intermediate lines of the two changed labels.
This is the same pair of paths used in the [repeat test](opposite_leaf_maps.md).
Conversely each of the three middle lines of $V$ gives one such rooted edge.
By the [middle-tip shape classification](opposite_edge_labels.md), a 16-cycle
has a shared first edge at exactly one of its three tips, and an 18-cycle
has none. This proves the count. At a shared edge both incident faces
also have the 16-cycle shape. The two edges of the opposite-pair graph
from $V$ to $\{A,B\}$ and $\{A,B'\}$ have label 4, since both retain the
tip where sharing occurs. A face has at most two edges with label 4.

In particular **one 16-cycle pair forces a rooted surface component with
at least fourteen faces**. Its root is the unique middle line carrying the
shared initial edge. If the witness has exactly fourteen faces, the shared
edge must be one of the three identified above.

There is a stronger test at the two endpoints of a 16-cycle pair $U$.
The two label-4 neighbors of $U$ each supply a shared edge incident to
its face at exactly one endpoint of $U$, namely the distinguished middle
line of that neighbor. The numbers of shared edges incident to the two
faces of $U$ therefore sum to two. In every piece of at most fourteen
faces the possible shared edges form a matching in the face adjacency
graph: for fourteen faces their endpoints are the six faces with exactly
one nonzero coordinate. Thus **if a 16-cycle pair belongs to a twelve-face
piece at one endpoint, its piece at the other endpoint has at least sixteen
faces**. The same conclusion holds for a fourteen-face piece when its face
misses all shared edges.

If every piece everywhere has at most fourteen faces, each 16-cycle
pair must consequently touch a shared edge at both endpoints. Under that
hypothesis the fourteen-face assignments in the verifier still pass:
their 16-cycle faces are exactly the endpoints of the shared edges.
This joint check leaves the small-component alternative open.

There is also a bound on all the component counts $k_A$. In a piece
with $f$ faces, let $s$ count its 16-cycle faces and $h$ its shared edges.
Then $h\le s\le f$, by counting the two incidences at each shared edge.
On an eight-face piece $s=h=0$, on a twelve-face piece $h=0$, and on a
fourteen-face piece $h\le3$. These bounds give $6(f-8)\ge2s+h$ in each
case, and the inequality also holds for $f\ge16$ because $2s+h\le3f$.
Summing over all pieces and roots gives
\[
  48\sum_A (256-k_A) \ge 5N_{16}.
\]
Here every opposite pair of 16-cycle shape is a face at both of its
endpoints and contributes one shared edge at its distinguished middle
line. The bound permits $N_{16}=0$ and does not force a grid.
