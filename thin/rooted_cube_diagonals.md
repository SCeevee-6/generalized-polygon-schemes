# Long indexed cycles forced by diagonals of rooted cubes

**Theorem.** In each opposite-pair component of a grid-free generalized
octagon of order $(2,4)$, let $t$ be the fraction of indexed edge
incidences on cycles of length at least twelve. Let $b$ be the fraction
of opposite pairs whose rooted pieces at both endpoints have more than
eight faces. Then
\[
                         43t\ \ge\ 1+2b .
\]
Consequently at least **3,266 distinct unordered indices** contain a long
cycle, and their long cycles account for at least **417,936 indexed edge
incidences**. The earlier bound was 2,340 indices. The argument uses no
automorphism assumption.

Counting marks over all 24 residues improves this to $t\ge1/36$ and
3,900 indices, and shows that $1/43$ is the limit of the present
argument. [Two marked residues](all_residue_marks.md)

## A count on one rooted surface

Use the [rooted surfaces](opposite_link_surfaces.md), splitting each
vertex into its separate link cycles. Call a degree-four vertex low and
a vertex of higher degree high. A high vertex corresponds to a long
indexed cycle. A small piece has eight faces and is an octahedron.
Its eight faces, viewed as partner lines, form a cube under adjacency
across a face edge.

In a large piece, take a triangle $xyz$ of low vertices. Across its three
edges write $x',y',z'$ for the replacements of $x,y,z$. The primed and
unprimed vertices in each color are distinct. The square links at
$x,y,z$ force $xy'z',x'yz',x'y'z$, in addition to the starting triangle
and its three neighbors. If any primed vertex were low, its square link
would force the eighth triangle $x'y'z'$. These eight triangles already
use both sides of every one of their edges. Connectivity across edges
would make them the entire piece. Hence all three primed vertices are
high. In particular each low vertex in a large piece has at least two
high neighbors. Otherwise it lies on an all-low triangle, which already
supplies two.

At a low vertex $B$ let $u(B)$ count its high neighbors whose opposite
vertex in the square link of $B$ is low. Write $F_j(B)$ for the number
of triangles at $B$ with exactly $j$ high vertices. Going around the
square gives the complete table.

| High neighbors of $B$ | $u(B)$ | $F_0(B)$ | $F_2(B)$ |
|---|---:|---:|---:|
| Two adjacent | 2 | 1 | 1 |
| Two opposite | 0 | 0 | 0 |
| Three | 1 | 0 | 2 |
| Four | 0 | 0 | 4 |

Thus $u(B)\le 3F_2(B)-F_0(B)$. Let $n_j$ count all faces of the
piece with $j$ high vertices, $f$ its face count, and $L$ the sum of its
high degrees. Summing over low vertices gives
\[
 U:=\sum_B u(B)\le 3n_2-3n_0\le 3(L-f),
 \qquad L-f=n_2+2n_3-n_0.
\]
In particular $L\ge f$ in every large piece. This also follows from the
[earlier Euler bound](indexed_eight_obstruction.md).

## From a surface link to a diagonal at another root

Choose a uniform frame $(R,Y;A,B,D)$ in the cubic component. At root
$R$ the face with partner $Y$ has the ordered vertices $A,B,D$. Let $a$
be the probability that this face lies in a large piece. The probability
that $A$ is high is $t$, since its link is the indexed cycle at $R,A$
with its subdivisions removed. Each rooted face supplies six frames.
In particular $a\le3t$ by the preceding count of $L$.

Exchange roles by passing to the frame $(B,A;R,Y,F)$, where $F$ is
the third middle line of $\{B,A\}$. This is a bijection on frames of the
component. The inverse recovers $D$ as the third middle line of
$\{R,Y\}$. If $A$ was high at $R$, the rooted piece at $A$ with partner
$B$ is large. Indeed the same long cycle of $J_{R,A}$ passes through
the adjacent cubic vertices $\{R,Y\}$ and $\{A,B\}$.

Define $w$ as the following event probability for the uniform frame
$(B,A;R,Y,F)$. If its piece at $B$
is small, take the unique square of its partner cube through $A$
retaining the middle line $R$. Let $A'$ be its diagonal partner.
The event counted by $w$ requires that the pieces at $A$ and $A'$ with
partner $B$ are both large. Given the small piece, this construction
is uniform over its 24 directed cube diagonals. A cube face used as a
partner has three choices of retained label, and either ordering of
the remaining two gives the same diagonal.

Start from an occurrence where $A$ is high at $R$. If the piece at $B$
with partner $A$ is large, both endpoints of $\{A,B\}$ have large
pieces, an event of probability at most $b$ after restricting to these
occurrences. Otherwise $B$ has square link at $R$. The four vertices of
that link are exactly the four partner lines in the square at $B$
retaining $R$. In both descriptions consecutive partners share an edge
of a rooted face containing the retained vertex. If the vertex opposite
$A$ in this link is high, the same indexed-cycle argument makes its
reciprocal piece large and supplies the event counted by $w$.

The remaining possibilities are bounded by $U$. Each ordered choice of
a high vertex $A$ and its low neighbor $B$ occurs in two frames at the
root, one for each incident triangle. Summing over large pieces and
dividing by the total number of frames bounds the probability lost
here by $3(t-a/3)$. A whole large piece contributes $2L$ frames with
$A$ high and $6f$ frames altogether. Consequently
\[
                       w\ge t-b-3(t-a/3)=a-2t-b .
\]

## Seventeen comparisons in the twenty-root test

Use the 1,152 reference frames and the obstruction to a map of all
four moves from the [24-cell argument](indexed_eight_obstruction.md).
The verifier chooses a new assignment tree. The twenty tested root
residues are joined by nineteen move-0 edges forming a tree $T$.
As before, each tested residue is first spanned internally. If all
twenty assigned pieces were small, each internal assignment would
respect every residue edge. Starting with these edges and the assignment
tree, the verifier fills **all 2,304 reference edges using only
unconditional relations**. Thus at least one test is large for each
starting frame.

For each vertex $x$ of $T$, consider its neighbors in $T$, viewed as
partner faces of the reference cube at $x$. Join two if they are a cube
diagonal apart. Their two face triples then share precisely one middle
label. Choose a spanning forest of this graph at each $x$. The
[verifier](indexed_eight_obstruction.py) checks that the forests have
**seventeen edges in total**. This count and the unconditional filling
are checked for the very same assignment tree.

For a uniform starting frame, mark a vertex of $T$ when its assigned
piece is large. Each vertex is marked with probability $a$. Both ends
of an edge are marked with probability $b$, since a move-0 join
exchanges the two endpoints and the internal tree words stay in their
rooted piece. For any chosen diagonal at $x$, the probability that $x$
is unmarked and its two designated neighbors are marked is $w$.
Indeed, at a small $x$ the reference and assigned cubes agree throughout.
A uniform frame at the reference entry position gives a uniform
identification of that cube with its assigned cube. The six orderings
at each of its eight faces give all 48 cube symmetries, which act
transitively on directed diagonals.

On the marked vertices keep the original edges with both ends marked.
At each unmarked center also keep every chosen diagonal with both ends
marked. This graph is a forest. One way to see this is to identify the
marked neighbors of an unmarked center using its original star in $T$.
Before processing that star the neighbors are in distinct components,
since the unique tree path between any two runs through that center.
Replacing the star by a forest on those neighbors preserves acyclicity.
The same argument works for all centers because shortcuts only traverse
unmarked internal vertices and have marked endpoints.

There is at least one mark, so the marked count is at least one plus
the retained edge count. Taking expectations yields
\[
 20a\ge1+19b+17w\ge 1+2b+17a-34t.
\]
It follows that $3a+34t\ge1+2b$. Substituting $a\le3t$ proves
the theorem.

There are $N=2{,}995{,}200$ unordered opposite pairs and $6N$
indexed edge incidences. Long cycles have lengths divisible by four.
The least multiple of four at least $6N/43$ is 417,936. At most 128
incidences belong to one unordered index. Dividing and rounding up
gives 3,266. Every pair of the 16-cycle shape has both pieces large, so
one can also retain $b$ and sum the stronger componentwise inequality.
An upper bound incompatible with these lower bounds would force a grid.
Whether a grid-free octagon can satisfy them is still open.

## When no opposite pair is large at both endpoints

There is an additional structural consequence if $b=0$ in a component.
Its long indexed cycles group into circuits of at least eight distinct
cycles, all of the same length within each circuit.

In any rooted piece, two high vertices $A,B$ cannot share an edge.
They would make both rooted pieces at the opposite pair $\{A,B\}$
large by the exchange of roles above. Take a high vertex $A$. Every
neighbor on its link is low. Each of these neighbors has a square link,
so it has a unique second neighbor in the color of $A$. At two successive
neighbors these second neighbors agree, by taking the other triangle
on their common edge. This produces a single vertex $A'$ across the
entire link of $A$. The triangles using $A$ or $A'$ on every edge of
the link form a closed suspension and exhaust the connected piece.
The two poles have the same high degree and the other vertices are low.
Every large piece has a high vertex by $L\ge f$.

Make a graph whose vertices are the large rooted pieces in the component.
For each long cycle of $J_{R,A}$ join the piece at $R$ containing its
split vertex $A$ to the piece at $A$ containing its split vertex $R$.
Each piece has exactly two poles and hence degree two. Its two poles
have distinct ambient line labels, since the two faces at any equatorial
edge have distinct triples of labels. This excludes parallel edges,
and the two roots of one index are distinct. The graph is therefore
a union of simple circuits. The indexed cycle length is twice the pole
degree, so it is constant along each circuit.

For a piece rooted at $R$, mark the point of $R$ at which its two poles
attach. They have the same attachment color. Along a circuit, the marks
of consecutive pieces are the two points on the intermediate line between
their roots, so they are distinct and collinear. The two intermediate
lines at an individual mark are distinct as well. To see this, take an
equatorial edge $DE$ of the suspension. The two incident faces have
partners $Y,Y'$ and poles $A,A'$. At the cubic vertex $\{D,E\}$ the
neighbors $\{R,Y\}$ and $\{R,Y'\}$ have third middle lines $A,A'$
beyond $D,E$. A pair of the 16-cycle shape would be large at both
ends, so every pair in this component has the 18-cycle shape.
Its three third middle lines are pairwise at line distance three by
the [failed-pair metric](opposite_boundary.md). If $A,A'$ used the same
intermediate line from $R$, their line distance would be at most two.
The marks on a circuit of length $s$ thus give a closed nonbacktracking
incidence walk of length $2s$. Girth sixteen implies $s\ge8$.

Run `python thin/indexed_eight_obstruction.py` from the repository root.
