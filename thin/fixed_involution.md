# Fixed elements of a grid-free octagon

In a grid-free generalized octagon of order $(2,4)$, the fixed incidence
graph of every nonidentity involution has the following form. It contains
the entire radius-three ball about one point, and it is contained in the
radius-four ball about that point. It has exactly 45 lines and one of
$31,51,71,91$ points. The value 91 is the point-centered central
involution case.

The proof uses the failed-pair stabilizer parity and the fixed-flag
identity in [failed_pair_types.md](failed_pair_types.md).

## The fixed incidence graph is a tree

The fixed graph is nonempty since there are an odd number of points.
Two fixed elements at incidence distance below eight have a unique
geodesic, which is fixed pointwise. Two opposite fixed points have five
geodesics. The involution fixes at least one of the five as a set and
fixes it pointwise because it fixes its endpoints. Two opposite fixed
lines constitute a failed opposite pair. The stabilizer parity forces
the involution to act evenly on the three arms. An even permutation of
order dividing two on three objects is the identity. Thus all three
geodesics are fixed pointwise. In particular, **two opposite fixed lines
are both fixed pointwise, including all their incident points**.

Consequently the fixed graph is connected and its distances agree with
the ambient distances. If it contains a cycle, its shortest cycle has
length 16. Indeed a longer shortest cycle would have length at least 18.
Two vertices separated by nine edges along it have a fixed path of odd
length at most seven. The short path and the nine-edge arc would give
a cycle of length at most 16.

Let $C$ be this fixed 16-cycle. Any fixed line $L$ has an opposite line
on $C$. Take a shortest path from $L$ to $C$ and continue along $C$ to
length eight. There is no backtrack at the join by minimality of the
initial path. The resulting path is geodesic by girth, and its endpoint
is a line. The preceding parity argument now shows that every fixed
line has all three points fixed. For any fixed point on $L$, the
geodesic from $L$ to its opposite line on $C$ starting at that point
supplies another fixed line through the point. Every fixed point has
odd fixed degree, since the involution acts on its five incident lines.
Thus the fixed graph is thick, has diameter eight and girth 16, and
contains all three points of each of its lines.

A finite thick generalized octagon is regular. The fixed suboctagon
would have order $(2,2)$ or $(2,4)$, since its point degrees are three
or five. Feit–Higman excludes $(2,2)$ because $2st=8$ is not a square.
Order $(2,4)$ gives the whole ambient octagon, contradicting nonidentity.
Hence the fixed graph is a tree.
[Regularity and Feit–Higman](https://cage.ugent.be/~hvm/artikels/225.pdf)

## The number of points and lines

Let $P,F,I$ count fixed points, lines and flags. The tree gives
$I=P+F-1$. The fixed-flag identity with $B=0$ gives
$64F+45=2925$, so $F=45$. Write $T=(3F-I)/2=(91-P)/2$
for the number of point transpositions on fixed lines and
$S=(5P-I)/2=160-4T$ for the number of line transpositions at fixed points.

Let $d_i$ count points mapped to point-graph distance $i$. For $i=1,2,3$
the geodesic has a fixed midpoint. A swapped pair on a fixed line
yields one pair at point distance one and eight at distance three.
For distance three, extend the path on one specified side by a choice
of four lines and two points. A swapped pair of lines at a fixed point
yields two pairs at point distance two, from the two choices of point on
one specified line. Each unordered pair contributes twice to $d_i$.
The complement gives $d_4$. Therefore
[
(d_0,d_1,d_2,d_3,d_4)=(91-2T,,2T,,640-16T,,16T,,1024).
]

The point-graph adjacency eigenvalues are $10,5,1,-3,-5$. The distance
polynomials computed from its intersection array give the following
columns of distance-matrix eigenvalues in that order.
[
M=egin{pmatrix}
1&1&1&1&1\
10&5&1&-3&-5\
80&10&-10&2&20\
640&0&-8&16&-80\
1024&-16&16&-16&64
end{pmatrix}.
]
For the permutation matrix of the involution, the trace of its
product with the $i$th distance matrix is $d_i$. Thus the traces on
the five adjacency eigenspaces are $M^{-1}d$, namely
$(1,31-3T/5,10,35-T,14-2T/5)$. Applying the same formula to the
identity gives dimensions $(1,351,650,675,78)$. Each trace is integral
and congruent to its dimension modulo two. The eigenvalue-five
eigenspace alone gives $T\equiv0\pmod{10}$.

## Location of the fixed tree

The tree has diameter at most eight. If its center is a line and its
diameter is at most six, there are at most $1+3\cdot4=13$ fixed lines.
If its center is a line and its diameter is eight, the endpoints of
a diameter are opposite fixed lines. Both have every incident point
fixed, extending the diameter in the tree. This is impossible.

If the center is an edge, delete it and consider the two rooted sides.
Each side has depth at most three. On the point side there are at most
$4+4\cdot2\cdot4=36$ lines. On the line side there are at most
$1+2\cdot4=9$ lines. Equality throughout is necessary for $F=45$.
All points at depths zero and two on the point side and depth one on
the line side are then fixed. On the line side the eight depth-two
lines each have zero or two fixed children, so $P=11+2a$ with
$0\le a\le8$. This gives $32\le T\le40$. The congruence forces $T=40$
and $a=0$. The two rooted depths are then three and two, which rules
out an edge center.

The center is therefore a point. The tree is contained in its
radius-four ball, which contains exactly $5+5\cdot2\cdot4=45$ lines.
It contains all those lines and the radius-three ball. The 40 outer
lines each contribute either both of their remaining points or neither.
Together with $T\equiv0\pmod{10}$ this gives $P\in\{11,31,51,71,91\}$.

Finally take a point $y$ opposite the center. It is moved. If its
displacement were less than eight in the incidence graph, the unique
geodesic to its image would have a fixed midpoint at distance at
most three from $y$. Every fixed element is within four of the center,
contrary to $y$ being opposite. Thus $y$ and its image are opposite.
Their five geodesics are permuted by the involution, and at least one
is preserved as a set. Its midpoint is a fixed point at distance four
from $y$. It must be at distance four from the center by the triangle
inequality. This excludes $P=11$, when the fixed graph would have no
point in the outer shell.

An involution with any other fixed structure forces a grid and hence
the Ree–Tits octagon by [grid_path_closure.md](grid_path_closure.md).
The four residual point counts are necessary conditions for a
grid-free octagon with an involution. This argument leaves open the
existence of a grid in an arbitrary octagon.
