# Exact component overlaps at collinear roots

In any generalized octagon of order $(2,4)$ let $X(x)$ be the graph on
points opposite $x$, with collinearity as adjacency. If $x,y$ are
collinear and components $C,D$ of $X(x),X(y)$ have sizes $256i,256j$,
then
\[
                         |C\cap D|=32ij.
\]
In particular all four intersections for two splits of type $512+512$
have size 128. The component contrasts at collinear roots are orthogonal
in the 351-dimensional point eigenspace of eigenvalue five. This settles
the [cross-cut question](../thin/tip_coupling_audit.md) without a symmetry
or grid assumption.

## Counting across a line

Write the root line as $F=\{x,y,z\}$. In its far point shell let $S_a$
be the 512 points projecting to $a$ at point distance three. A line
opposite $F$ has one point in each class, since the three arms between
opposite lines give a bijection on their points. Each point of a class
is on four opposite lines: its unique inward line toward $F$ accounts
for the fifth. The edge between the $S_y,S_z$ points of any such line
belongs to a single component at $x$. Call that the line's label at $x$,
and define the labels at $y,z$ cyclically.

We use the [boundary count](README.md): at every point $b$ of point
distance three from $x$, exactly $i$ of the four outward lines have
their two outer points in $C$. In particular this applies at every
$b\in S_x$. Also $|D\cap S_x|=128j$: the outer shell at $y$ splits
as $S_x,S_z$, and the four cross edges per vertex form a regular
bipartite graph even upon restriction to the component $D$.

Count lines opposite $F$ with labels $C,D$ at $x,y$. From $S_x$ there
are $128j$ possible points in $D$, and exactly $i$ choices of line at
each. From $S_z$ there are exactly four choices at each point in
$C\cap D$. Equating the counts proves the formula. For the contrasts
$1_C-(i/4)1_{X(x)}$ and $1_D-(j/4)1_{X(y)}$, expansion of the inner
product gives $|C\cap D|-32ij=0$.

## Three independent component labels

Take also a component $E$ of size $256k$ at $z$. The number of lines
opposite $F$ with all three prescribed labels $C,D,E$ is **$32ijk$**.
Indeed their point in $S_x$ belongs to $D\cap E$, which has $32jk$
points by the pair formula. Exactly $i$ outward lines there have label
$C$. Thus for a uniformly chosen line opposite $F$ the three labels
are mutually independent with respective probabilities $i/4,j/4,k/4$.
The counts sum to 2048. This couples the component decompositions
simultaneously at all incident root triples. It permits every previously
allowed component-size partition at the level of these cell counts.

Run `python octagon/check_collinear_components.py` from the repository
root. The model checker reconstructs every outer component and counts
the intersections and line labels directly on all 2925 lines.

## A rank constraint coupling all roots

Let $Q_x$ be the orthogonal projector onto the component contrasts at $x$.
Its rank is the number of components at $x$ minus one, since the
component indicators have disjoint supports and one weighted mean is
removed. Write $R$ for the sum of these ranks. Every projector has image
in the same 351-dimensional eigenspace. The pair formula gives $Q_xQ_y=0$
at distance one. Applying Cauchy–Schwarz to the eigenvalues of
$\sum_x Q_x$ and expanding the trace of its square therefore gives
\[
  \sum_{d(x,y)\ge2}\operatorname{tr}(Q_xQ_y)\ \ge\ R^2/351-R .
\]
The sum uses ordered pairs. Equality for positive $R$ requires the sum
of projectors to be a scalar on the entire eigenspace.

**Model calculation.** In Ree–Tits every root has two equal components.
Let $F_{xy}$ be $+1$ or $-1$ according to the component of $y$ at $x$,
and zero off the outer shell. The ordering of components at each root
is arbitrary. The verifier checks all entries of the two identities
\[
 FF^T=1024I+64H,\qquad H^2=1024I+48H,\qquad H_{xy}=F_{xy}F_{yx}.
\]
Thus distinct contrasts are orthogonal even at distances two and three
in this model. At opposite points their normalized inner product has
magnitude $1/16$. The eigenvalues of the signed opposition matrix $H$
are 64 and $-16$ with multiplicities 351 and 1404 by its zero trace.
The normalized contrasts attain the global trace bound with $R=1755$.
This exact model identity suggests a stronger recognition target. Its
extension to arbitrary octagons, including those with connected outer
shells, remains open.
