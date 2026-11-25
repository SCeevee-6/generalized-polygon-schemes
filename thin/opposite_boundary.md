# A global complex and spectral constraints on opposite pairs

The [surfaces at individual root lines](opposite_link_surfaces.md)
can be tied together without assuming automorphisms. This gives
an exact boundary matrix over $\mathbb F_2$ and identifies which
moments of the same matrix are fixed by the ambient line spectrum.
The constraints below leave the existence of a grid open.

## One complex for all root lines

Suppose the octagon is grid-free. Make a simplicial complex with
the ambient lines as vertices and the unordered opposite pairs as
edges. For each opposite pair $U$, insert the triangle consisting
of its three middle lines. The triangle is unique to $U$. If two
distinct pairs had the same triangle, those two vertices of the
opposite-pair graph would share three neighbors. This gives a
four-cycle and hence a grid.

Index both the edges and the triangles by the opposite pairs. Over
$\mathbb F_2$ the triangle boundary matrix is exactly the adjacency
matrix $J$ of the cubic opposite-pair graph. In fact the three edges
of the triangle for $U$ are exactly the three neighbors of $U$.
If $D$ is the line-by-opposite-pair endpoint incidence matrix, then
$DJ=0$ over $\mathbb F_2$. The line opposition graph is connected
because any two lines have a common opposite line, by the size 2048
of each opposite shell in a set of 2925 lines. Thus
\[
\operatorname{nullity}_{\mathbb F_2} J\ge 2924.
\]

The link at a line $A$ has vertices $B$ opposite $A$ and an edge
$BB'$ if $\{A,B,B'\}$ is one of the middle-line triangles. This is
exactly the triangle adjacency graph of the fixed-root construction.
Indeed $B,B'$ are adjacent in that construction when they share
two middle lines $C,D$ with $A$. The three middle lines of $\{C,D\}$
are then $A,B,B'$.

There is a stronger rank constraint that uses all the components of
all these links. Put $k_A$ for the number of components of the link at
$A$. Split the vertex $A$ by those components. Each edge $\{A,B\}$
goes to the copy of $A$ indexed by the component containing $B$, and
likewise at $B$. The three edges of an original triangle still form a
triangle, since the two neighbors at each vertex were adjacent in its
link. If the resulting graph on the split vertices has $k$ components,
its vertex-edge incidence matrix has rank $\sum_A k_A-k$. The
triangle boundary still vanishes after applying that incidence matrix,
so
\[
\operatorname{nullity}_{\mathbb F_2} J\ge \sum_A k_A-k.
\]
In particular the link component counts constrain a single global
matrix. This bound alone allows grid-free systems.

## A real-valued version

In the following expressions all matrices have real entries. Write
$L_i$ for the matrix on ambient lines indicating line distance $i$.
Let $R$ be the line-by-opposite-pair matrix indicating membership in
the triple of middle lines. There are 2048 opposite partners of a
line. Each pair of lines at distance two has 64 continuations to an
opposite partner via the second line. Each pair of opposite lines
shares three middle lines. The resulting identities are
\[
DD^{\mathsf T}=2048I+L_4,\qquad
DJ=2R,\qquad DJD^{\mathsf T}=128L_2,\qquad
RR^{\mathsf T}=3072I+3L_4.
\]
For the last identity, two distinct lines can occur together among
the middle lines of $U$ exactly when their opposite pair is a
neighbor of $U$, which gives three choices. A given line occurs in
3072 triples, since its 96 lines at distance two each have 64
opposite continuations and each pair using the designated middle
line is counted twice.

Let $S$ have rows indexed by opposite pairs $U$ and columns by lines.
For each neighbor of $U=\{A,B\}$ record its third middle line beyond
$A,B$, with multiplicity. In a grid-free octagon these three lines
are distinct. A repeated third line would give a second common
neighbor of two neighbors of $U$, forming a four-cycle. Thus
\[
J^2D^{\mathsf T}=6D^{\mathsf T}+2S,\qquad DS=3L_4.
\]
For the second identity, a third line $F$ recorded for $U$ is opposite
both endpoints of $U$. Conversely, for fixed opposite lines $A,F$,
each neighbor of $\{A,F\}$ has three middle lines $A,F,B$ and records
$F$ once for $\{A,B\}$. This accounts for the three contributions.

## The omitted middle lines

Let $t$ be the fraction of opposite pairs with middle-tip shape
$(3,4,4)$ in the [edge labels](opposite_edge_labels.md). The three lines
selected by a row of $S$ have mutual line distances $(2,3,3)$ for this
shape and $(3,3,3)$ for the $(3,3,3)$ shape.

Consider adjacent pairs $U,V$. Write their middle triples as
$V\cup\{X\}$ and $U\cup\{Y\}$. Eight of the nine cross distances
between the triples are two. **The ninth distance $d(X,Y)$ is four
on every failed edge and two on every grid edge.** The
[verifier](opposite_boundary.py) finds the original middle line $X$
by its unused tip. It checks all three choices of omitted line in
each of the 24 labelled 16-cycle metrics, eight 18-cycle metrics and
one centered metric. An earlier calculation selected the last connector
of each arm in place of its middle line. This corrects the
distance-three frequency to zero.

In a grid-free octagon let $H_{XY}$ count the undirected edges with
omitted pair $\{X,Y\}$. Set $H_{XX}=0$. Then
\[
RJR^{\mathsf T}=256 L_2+H.
\]
Indeed an entry counts directed edges $(U,V)$ with the first specified
line in the middle triple of $U$ and the second in that of $V$.
For lines at distance two, requiring the second line to be an endpoint
of $U$ gives $64\cdot3$ choices. Requiring the first line to be an
endpoint of $V$ gives the same number, with $64\cdot2$ choices in both
counts. The remaining contributions are precisely the omitted pairs,
one in each orientation of an undirected edge. Their entries count
[squares indexed by opposite roots](opposite_root_squares.md).

## Exact projected moments

For a line adjacency eigenvalue $e$, let $E$ be its orthogonal
projector and let $a_i$ be its cosine at line distance $i$. These
cosines come from the line intersection array
$\{12,8,8,8;1,1,1,3\}$. In particular
$\operatorname{tr}(EL_i)/\operatorname{rank}E$ equals the size of
the $i$th shell times $a_i$. Define
\[
m_j=\frac{\operatorname{tr}(E D J^jD^{\mathsf T})}{\operatorname{rank}E}.
\]
The incidence and distance counts give
\[
\begin{aligned}
m_0 &=2048(1+a_4),\\
m_1 &=12288a_2,\\
m_2 &=12288(1+2a_4),\\
m_3 &=12288(8a_2+a_4),\\
m_4 &=36m_0+24\cdot6144a_4
      +4096[3+6a_3+2t(a_2-a_3)].
\end{aligned}
\]
For $m_3$, use $DJ^3D^{\mathsf T}=4RJR^{\mathsf T}$ and the nine
cross distances on each directed edge. For $m_4$ use the squared
norm of $6D^{\mathsf T}+2S$ on the eigenspace. Each row of $S$
contributes three diagonal terms and six off-diagonal terms to
that last Gram calculation. The ratio of opposite pairs to lines
is 1024.

The cubic spectrum of $J$ lies in $[-3,3]$. Therefore the moment
matrices $(m_{i+j})_{0\le i,j\le2}$ and
$(9m_{i+j}-m_{i+j+2})_{0\le i,j\le1}$ are positive semidefinite.
For example the second matrix amounts to inserting the positive
semidefinite matrix $9I-J^2$ into the squared norm. Exact rational
arithmetic verifies these two inequalities for every $0\le t\le1$.
By convexity it suffices to check the two endpoints for each
nonconstant line eigenspace. The constant
eigenspace has $m_j=4096\cdot3^j$. Thus these low moments do not yet
give an incompatibility.

Run the [verifier](opposite_boundary.py) with
`python thin/opposite_boundary.py`.
