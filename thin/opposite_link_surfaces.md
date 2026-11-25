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
Run [the replay](opposite_link_surfaces.py) with
`python thin/opposite_link_surfaces.py`.

Both attachment types survive this test. The surface constraint
couples all indexed systems at a fixed line. Compatibility of
these surfaces across different roots remains open.
