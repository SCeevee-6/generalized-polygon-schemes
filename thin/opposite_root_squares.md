# Squares indexed by opposite root lines

If an octagon is grid-free, then **every line has at least 32 opposite
partners whose auxiliary graph contains two vertex-disjoint squares**.
A uniform upper bound of one square per opposite pair would therefore force a grid.

## The auxiliary graph

Use line-graph distances. For opposite lines $X,Y$, take the two sets
\[
  P=\{A:d(X,A)=2,\ d(Y,A)=4\},\qquad
  Q=\{D:d(Y,D)=2,\ d(X,D)=4\}.
\]
Join $A\in P$ to $D\in Q$ when $d(A,D)=2$. Both sides have 66 vertices.
For example, the line intersection array gives
$L_2=L_1^2-3L_1-12I$. From a line in shell four about $Y$, three
neighbors lie in shell three and nine in shell four. The corresponding
entry of the shell quotient for $L_2$ is
$3\cdot8+9^2-3\cdot9-12=66$.

Every vertex has degree two. At $A\in P$ the three middle lines of
$\{A,Y\}$ attach at the three different points of $A$. The two that
attach away from the point towards $X$ are opposite $X$, by the
nonbacktracking eight-edge path through $A$. The third has an incidence
walk of length at most six to $X$ through the common attachment point.
Interchanging the roots proves the degree on the other side. Write
$c(X,Y)$ for the number of four-cycle components. In particular
$0\le c(X,Y)\le33$.

A square with sides $\{A,B\}\subset P$ and $\{D,E\}\subset Q$ gives
an edge of the cubic opposite-pair graph. Indeed $A,B$ are two middle
lines of $\{X,D\}$ and hence are opposite. Likewise $D,E$ are two
middle lines of $\{Y,A\}$. The middle triples of $\{A,B\}$ and
$\{D,E\}$ are exactly $\{D,E,X\}$ and $\{A,B,Y\}$. Conversely every
cubic edge with omitted lines $X,Y$ supplies exactly this square.
The [omitted-line calculation](opposite_boundary.md) puts those two
lines at distance four on every failed edge, and at distance two on a
grid edge.

## Counting at a single line

Fix $X$. There are 3072 pairs $U$ with $X$ among their middle lines.
Choose the two attachment points on $X$ in three ways, then one line
at distance two through each attachment in $32^2$ ways. The two chosen
lines are opposite via $X$. Each $U$ gives one cubic edge omitting
$X$, by selecting its other two middle lines as its neighbor.

Let $g_X$ count grid components of the cubic graph whose six lines
include $X$. In such a component there are three edges omitting $X$,
one for each line of the other family. All remaining edges just counted
have an opposite omitted partner $Y$. The bijection with squares gives
the exact identity
\[
   \sum_{Y:\,d(X,Y)=4} c(X,Y)=3072-3g_X .
\]
In the grid-free case the 2048 summands average $3/2$. If only $k$ of
them are at least two, their sum is at most $2048+32k$, proving
$k\ge32$. Globally at least $2925\cdot32/2=46{,}800$ unordered opposite
pairs must have two such squares. Counting pairs of squares with
multiplicity gives a further bound
\[
  \sum_{\{X,Y\}:\,d(X,Y)=4}\binom{c(X,Y)}2
      \ge \tfrac12\,2925(3072-2048)=1{,}497{,}600,
\]
using $\binom c2\ge c-1$. All squares at fixed roots are
vertex-disjoint because the graph is 2-regular.

## Attachment counts

Label the points on $X,Y$ so that $x_i,y_i$ are at incidence distance
six. All other cross distances are eight. Color each vertex on either
side by its attachment point on its root line. The $3\times3$ edge
count by color is 18 on the diagonal and 13 off the diagonal.

Indeed an auxiliary edge $AD$ specifies a nonbacktracking incidence
walk of length twelve from $X$ via $A,D$ to $Y$, by joining the three
four-edge geodesics. Its two eight-edge segments are geodesic by
opposition. Conversely the lines at positions four and eight on any
such walk give an auxiliary edge. Prescribing the attachment colors
leaves a ten-edge walk from $x_i$ to $y_j$ whose first edge avoids $X$
and whose last edge avoids $Y$.

For the count write $b_r(u,v)$ for the number of nonbacktracking walks
of length $r$. Set $(u_0,u_1)=(x_i,X)$ and $(v_0,v_1)=(y_j,Y)$.
Removing walks starting with the forbidden edge replaces the count by
the unrestricted count minus a shorter count with the first edge
reversed. The same recurrence holds at the last edge. The two endpoint
sets are disjoint, so the length-zero boundary introduces no extra term.
Expanding gives
\[
 \sum_{a+b\le10}(-1)^{a+b}
 b_{10-a-b}(u_{a\bmod2},v_{b\bmod2}),\qquad a,b\ge0.
\]
The [verifier](opposite_root_squares.py) evaluates this using the nine
distance shells about the target. At length two it subtracts all immediate
reversals. At subsequent lengths it subtracts $(\deg(u)-1)b_{r-2}(u,v)$
from the sum over first neighbors. This gives 18 for $i=j$ and 13 otherwise.

Each square uses two colors at each root. If the omitted colors agree,
the square uses two diagonal edges. Otherwise it uses one. Agreement
is exactly [edge label 3](opposite_edge_labels.md), since the omitted
points are the tips that define the label. Thus if $c_3,c_4$ count
squares of the two labels at fixed roots, then
\[
 2c_3+c_4\le54,\qquad 2c_3+3c_4\le78 .
\]
For example, if all failed pairs have the 18-cycle shape in a grid-free
octagon, all edges have label 3 and $c(X,Y)\le27$. The preceding count
then forces at least $\lceil1024/26\rceil=40$ partners with two squares
at every line.

## Testing two squares together

Take the first lexicographically ordered failed metric in
[the surface verifier](opposite_link_surfaces.py). Its omitted roots for
the edge retaining arms zero and one give a square. Take two copies
and identify the root lines, trying all six permutations of the points
on the first root. Projection fixes the identification on the second
root. Require all four corner lines of one square to differ from all
four of the other, and impose the same failed metric on the shared
opposite pair. All six choices pass three-line propagation. This checks
one choice of local metrics and supplies no ambient completion.
[Overlap replay](opposite_square_overlap.py)

## Check in the Ree–Tits model

For all 2048 partners opposite line 0, the component lengths are
$6,12,12,14,14,14,30,30$. This checks the construction in the known
octagon. The argument above leaves open whether the required double
squares can coexist throughout a grid-free octagon.

From the repository root, run `python thin/opposite_square_overlap.py`
for the local test and `python thin/opposite_root_squares.py` for the
[walk counts and model check](opposite_root_squares.py).
