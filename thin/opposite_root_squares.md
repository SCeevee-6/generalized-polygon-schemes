# Squares indexed by opposite root lines

If an octagon is grid-free, then **every line has at least 45 opposite
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
## Separating the spare points on the root arms

A vertex of color $i$ has neighbors of two distinct colors, since they
are middle lines of an opposite pair with the other root as endpoint.
Call the vertex special if its neighbors omit color $i$. Its color class
has $(18+13+13)/2=22$ vertices. Exactly four are special, since 18
vertices have a neighbor of color $i$.

These four lines have a common point. On arm $i$ from $X$ to $Y$,
let $u_i,v_i$ be the unused points of the first and last connector,
as in [the opposite seed](opposite_seed.md). The four lines through
$u_i$ other than its connector have color $i$ and are opposite $Y$.
Indeed the seven-edge geodesic from $u_i$ to $Y$ starts through the
connector and ends through $y_i$. For any of the four lines, its arm
towards $Y$ through $u_i$ therefore ends at $y_i$. The two arms supplying
auxiliary neighbors end at the other two points of $Y$. This identifies
all four special lines. On the other side the common point is $v_i$.

For $i\ne j$ there is **exactly one edge** between special vertices
of left color $i$ and right color $j$. The seed gives $d(u_i,v_j)=3$
in point distance, with both other points of either connector at distance
four from the opposite spare. The unique six-edge geodesic between
the spares thus starts and ends on special lines, which are at line
distance two. Conversely a special-to-special edge gives this geodesic
by using its four-edge connection between lines. The connection starts
and ends away from the spares by the auxiliary neighbor rule.
In this ordered color block four edges touch a special left vertex and
four touch a special right vertex. Exactly $13-4-4+1=6$ edges have
two ordinary endpoints.

In a square every left vertex misses the omitted right color, and
conversely. If the omitted colors differ, precisely one vertex on each
side is special. Their edge determines the square because every vertex
has degree two. There are only six such edges, so $c_4\le6$.
If the omitted color is the same on both sides, all four vertices
are ordinary. For each fixed omitted color the square uses an edge in
each of the two off-diagonal blocks on the remaining colors. Either
block has only six ordinary-to-ordinary edges, so $c_3\le3\cdot6=18$.

Thus $c(X,Y)\le24$. In the grid-free case the 3072 squares at $X$
force at least $\lceil1024/23\rceil=45$ opposite partners with two or
more squares, giving at least 65,813 unordered root pairs globally.
If all pairs have the 18-cycle shape then $c_4=0$, and the bound rises
to $\lceil1024/17\rceil=61$ partners at every root.

## Six collinearity tests for label 4

Fix opposite roots $X,Y$. Each ordered pair $i\ne j$ gives a test with
no choices. Write the unique geodesic between the spares $u_i,v_j$ as
three successive lines $A,K,D$, and let $p,q$ be the unused points on
$A,D$ along this path. Put $k=3-i-j$ with colors numbered $0,1,2$.
Let $E$ be the middle line of the geodesic from $p$ to $y_k$, and
$B$ the middle line of the geodesic from $q$ to $x_k$. Write $s,t$
for the unused points on $E,B$ along these paths. Then **the special
edge $AD$ belongs to a square exactly when $s,t$ are collinear**.

Indeed the three points on $A$ project bijectively to the three points
on its opposite line $Y$. The point $u_i$ projects to $y_i$. The
attachment toward $D$ projects to $y_j$ along the arm through $D$.
Thus $p$ projects to $y_k$ and the six-edge geodesic used above exists.
Its middle line $E$ is the other auxiliary neighbor of $A$. Symmetrically
$B$ is the other neighbor of $D$. A square closes exactly when $B,E$
are at line distance two. In that case the connecting arm at $E$ must
use its point away from both $A$ and $Y$, namely $s$: the two auxiliary
neighbors attach at the two points away from the root attachment.
The same argument gives $t$ at $B$, proving necessity. Collinearity
gives distance at most two between $B,E$, and their distances from $X$
differ by two, proving sufficiency. Each label-4 square has just one
special-to-special edge, so the six indicators sum to $c_4(X,Y)$.

These local indicators have an exact global total. Sum over **unordered**
opposite roots, choosing either orientation once. Every successful test
counts the unique cubic edge with those omitted lines and that square.
The label-4 subgraph has degree two at every 16-cycle pair and zero
elsewhere. Hence the sum is exactly the number of 16-cycle pairs.
For example, if all pairs are grid-free and some have the 18-cycle shape,
at least as many root pairs have zero successful tests as there are
18-cycle pairs. A root pair with multiple successes increases that deficit
by its number of successes beyond the first.

The [replay](opposite_square_geodesics.py) tries to force a contradiction
by requiring all six tests to succeed at once. For each of the five
failed root metric orbits it names the connecting lines, requires all
24 corner lines to be distinct, and adds the three geodesics between
the unused middle-line tips of each of the twelve cubic endpoints.
All five pass three-line propagation on 116 line slots. This gives
surviving distance domains only. The compatibility of these tests
between different root pairs remains open.

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
