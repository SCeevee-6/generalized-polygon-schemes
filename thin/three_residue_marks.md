# Three marked residues among the twenty-four

In each opposite-pair component of a grid-free generalized octagon of order
$(2,4)$, at least **one eighth** of the rooted faces belong to pieces with
more than eight faces. At least **one twenty-fourth** of its indexed edge
incidences belong to cycles of length at least twelve. Globally this
requires **748,800 incidences** on long cycles and **5,850 distinct unordered
indices** containing them.

The previous count gave 3,900 indices. Its only remaining way to have two
marked residues was to mark two adjacent vertices of the reference 24-cell.
The argument below excludes all 96 choices by a finite computation.
[Two-mark argument and assignment tree](all_residue_marks.md)

## What can be identified with two exceptional residues

Assign actual frames to the 1,152 reference frames by the fixed tree words.
Omit the two putative marked residues. All internal edges of every other
vertex residue respect the frame moves, since its assigned piece is an
octahedron. Tree edges also respect their moves. Repeatedly fill the last
missing edge of an unconditional two-color polygon. Record which edges have
been filled. At the omitted residues each tree word still gives an actual
frame, even if some loops there fail to close.

Start with five line slots per frame $(A,B;C,D,E)$. Identify entries preserved
by a known move, exchanging entries for the three swap moves. Move 1 retains
the first, third and fourth entries. Record the full middle triple of each
pair $A,B$. A known move 1 also identifies the full middle triple of $C,D$
as $A,B$ and the new partner. Whenever two triples for the same pair
already share two entries, identify their third entries. Iterate. Each
identification follows from the three arms of an opposite pair.

Each end of a recorded pair has line distance two from each middle.
The ends are opposite and the three middles are pairwise opposite.
At every assumed small root the two middles of each attachment color
have distinct intermediate lines, giving line distance three via the
common attachment point on the root. Two partners on a face diagonal of
its cube also have line distance three. The reference test for the latter
is that their face triples share exactly one facet.
[Eight-face attachment and diagonal calculation](opposite_link_surfaces.md)

On a line $L$, a known neighbor at distance two has a unique attachment
point. Opposite neighbors must attach at different points. Enumerate
the proper three-colorings of this opposition graph, up to color names.
Two neighbors assigned the same color in every coloring can share one
point slot. Use a clique of up to three neighbors to name distinct
points on $L$, completing the line with spare slots as necessary.
Any remaining attachment is at point distance zero or one from each
of its three points. Insert the full connector line for every recorded
distance-two pair. For opposite lines require every cross point distance
to be at least three. For lines at distance three require at least two.
Other unspecified cross distances include zero. Thus separate names for
lines or points can still have identical images.

Propagate these domains with the existing 39 ordered full-line distance
patterns and their three-line composition table. Retain singleton point
distances for the next test. [Seed construction](partial_frame_metrics.py)
and [line propagation](grid_path_closure.md)

## First lines of short geodesics

Fix a point $p$. For every target $x$ at known point distance one, two or
three, the incidence geodesic to $x$ has a unique first line. Put
$i=d(p,x), j=d(p,y), k=d(x,y)$.

- If $i+j+k<8$, the three incidence geodesics form a tree by girth sixteen.
  For $i+j>k$ the first lines to $x,y$ coincide.
- If $k=i+j$, the geodesic through $p$ forces different first lines.
- If $i+j+k=8$ and the first lines coincide, the union again has fewer
  than sixteen distinct incidence edges and is a tree. For $k=2,3$ any
  named interior point $z$ on the geodesic from $x$ to $y$ must satisfy
  $2d(p,z)=i+j-k+|i-j+k-2d(x,z)|$. This is the distance to the branch
  point plus the distance along the $xy$ path, in incidence units. Any
  violation forces different first lines.

Identify the forced equalities and retain the known inequalities between
the resulting classes. Two classes that both avoid the same four pairwise
different first lines must use the fifth line through $p$. Iterate these
identifications. For each omitted pair the certificate specifies four
point slots $p,r,q,s$. The pairs $p,r$ and $q,s$ are opposite and the
four cross distances are three. At each corner the two cross geodesics
must start on the same line, unless the first-line constraints already
give a contradiction there. This is exactly the impossible configuration
proved in [The four remaining geodesics](indexed_eight_obstruction.md#the-four-remaining-geodesics).
Its obstruction uses the three points per line as well as girth.

## Coverage and count

The checker regenerates the reference graph and independently matches the
certificate keys against all adjacent omissions. For the 180 nonadjacent
pairs, single-missing-edge filling gives the entire forbidden frame map.
The same test covers the 24 single omissions and, by monotonicity, zero
omissions. The rectangles cover the other 96 pairs. Thus **every starting
frame marks at least three residues** on this assignment tree.

For a uniform starting frame in a component, each of the 24 fixed words
acts bijectively. Each residue is marked with probability $a$, the
large-piece fraction. Taking expectations gives $24a\ge3$. In every
large piece the sum of the degrees above four is at least its number of
faces, and there are three face–vertex incidences per face. Thus the long
fraction satisfies $t\ge a/3\ge1/24$. There are $N=2925\cdot2048/2$
opposite pairs and six indexed edge incidences per pair. The count follows
from $6N/24=748800$ and the capacity of 128 incidences per unordered index.
[Rooted degree and index counts](indexed_eight_obstruction.md#counting-long-cycles)

Run `python thin/frame_direction_obstruction.py` for the exhaustive replay,
or append `--pair 0 6` for one certificate. The optional `--model-control`
checks inferred first lines against actual geodesics at three roots of
the stored Ree–Tits model. The optional `--twelve` applies the same
first-line and rectangle search to each of the thirteen stored joint
twelve-cycle witnesses. All thirteen pass this local test. Realizability
of those witnesses and existence of a grid in an arbitrary octagon remain
open. [Checker](frame_direction_obstruction.py) and
[certificate slots](frame_direction_rectangles.txt)
