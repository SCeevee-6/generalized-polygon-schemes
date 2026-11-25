# Distance blocks in the local rank test

The local rank test associates a $64\times64$ matrix $M$ with two
lines $A,C$ at line distance two. Its rows are lines $B$ opposite
$A$ at distance two from $C$. Its columns are lines $D$ opposite
$C$ at distance two from $A$. Entry $M_{B,D}$ is one when $B,D$
are at line distance two. A rank below 48 forces a grid.
[The opposite-pair cycle proof](opposite_cycles.md)

Here are universal distance constraints on the matrix and its
surrounding incidence distances. An abstract rank-48 example satisfies
these constraints, shell bounds and projector positivity, so those
tests together do not force a grid.

## The branch partitions

Continue the unique incidence path $A$ to $C$ by four edges. The
branches through the two available exit points of $C$ have 32 lines
$B$ each. Within a branch there are four choices of intermediate line
and then two choices of its point $b$ away from $C$. There are four
outward lines $B$ on each such $b$. This partitions the row side
into groups of 32, 8 and 4. The column side partitions in the same way
via paths from $C$ through $A$, with ancestor point $d$.

Distances within either side are the rooted tree distances, because
the tree path between two leaves has incidence length at most eight.
For fixed $B$, each of the two column branches supplies exactly one
$D$ at line distance two. They are the two middle lines of the
opposite pair $\{A,B\}$ other than $C$. These middle lines use the
two other attachment points on $A$.

## The four by four blocks

Fix ancestor points $b,d$. Their path through $C,A$ has incidence
length ten and is nonbacktracking. The incidence distance between
them is therefore six or eight. A shorter geodesic together with
the ten-edge path would contain a cycle of length below 16.

If the distance is six, the short geodesic cannot share an edge with
the ten-edge path. It therefore uses an outward line $B_*$ at $b$
and an outward line $D_*$ at $d$. The block has one distance-two
entry $(B_*,D_*)$. The rest of that row and column is at distance
three. The remaining nine entries are at distance four. In each
case, keep the corresponding end edges of the six-edge path or
extend them to the specified lines. The resulting paths have lengths
four, six and eight and are geodesic by girth.

If the distance is eight, there are five geodesics between $b,d$,
one for each incident line at each end. The two inward lines cannot
belong to the same geodesic. That geodesic would share two edges
with the ten-edge path and give a cycle of length at most 14.
Exactly three geodesics use outward lines at both ends. Each gives
one entry at line distance three and those entries lie in distinct
rows and columns. Conversely, a line geodesic of incidence length
six between a row and column, extended to $b,d$, must give a
geodesic of length eight. Thus the other thirteen entries are four.

## Exact margins on the groups of eight

Use the normalized vectors $e_L$ in the line graph eigenspace for
eigenvalue $-3$. Their inner products are $(-1/4)^{d(L,L')}$,
by the line-graph intersection array
$\{12,8,8,8;1,1,1,3\}$. This is the kernel of point-line incidence.
Indeed the transpose-times-incidence matrix on lines is the
line adjacency matrix plus $3I$. In particular the sum of the five
vectors on lines through a point is zero.

For an intermediate line $T$ in one column branch, sum the vectors
of all eight descendants $D$. The four outward lines at either of
the two ancestor points have vector sum $-e_T$. Hence the sum over
the group is $-2e_T$.

For a fixed row $B$, the line $T$ has distance three from $B$ if
it occurs on the corresponding $A$ to $B$ geodesic, and distance
four otherwise. Each branch contains one such $T$ at distance three.
That group is precisely the one containing the entry of distance two.
Take inner products with $e_B$ and scale by $256$. An entry of distance
two, three or four contributes $16,-4,1$, respectively. The target sum
is $8$ in the group with a distance-two entry and $-2$ in each of the
other groups. The resulting counts are:

| Group of eight | Distance two | Distance three | Distance four |
|---|---:|---:|---:|
| Contains the matched column | 1 | 3 | 4 |
| Each other group | 0 | 2 | 6 |

By symmetry the same rule holds for columns. Each row and column has
counts $(2,18,44)$ over the full block.

## An abstract example at the threshold

Encode a row or column by $(i,u,s,p)$ with branch $i\in\{0,1\}$,
intermediate line $u\in\mathbb Z/4$, ancestor child $s\in\{0,1\}$,
and outward slot $p\in\mathbb Z/4$. For each row and target branch
$j$, put a distance-two entry at
\[
(j,u+p,s\mathbin{\oplus}(ij),p).
\]
The four matchings around the two row branches and the two column
branches return to the same $u,p$ and flip $s$. They give sixteen
eight-cycles, so $M$ has rank 48 over every field.

Fill each four by four block having a distance-two entry with the
required row and column of threes. In each remaining block, omit
the row and column whose respective match is in the other ancestor
of the same group of eight. Match the remaining three rows and
columns in order, putting a three on each of those positions and
four elsewhere. This realizes both the block rules and the margins.

The verifier also adjoins the ancestors, intermediate lines and
root lines to this table. The 179 named elements consist of 147
lines and 32 points. It checks triangle inequalities, the short
four-point test, all two-anchor shell upper bounds on named elements,
and the absence of a grid among the named lines. For the spectral
test, pair each normalized point projector with its corresponding
line projector. In each point eigenspace, represent a line by the sum
of its three point vectors. The mixed entry at incidence distance
$2h+1$ is then $k_h+2k_{h+1}$ for the point cosine sequence $k$.
The line block is scaled by the squared norm $3(1+2k_1)$.
The remaining line-only eigenspace is checked separately. All the
Gram tests use integer arithmetic.

The example certifies compatibility with these tests. Completing
its named lines subsequently excludes the matching table itself.
[Full-line obstruction](indexed_cycle_geometry.md)

Run the [verifier](indexed_cycle_blocks.py) using
`python thin/indexed_cycle_blocks.py`.
