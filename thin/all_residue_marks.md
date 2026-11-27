# Two marked residues among all twenty-four

The [three-mark argument](three_residue_marks.md) excludes the adjacent two-mark
case below and raises the index bound to 5,850.

**Theorem.** In each opposite-pair component of a grid-free generalized
octagon of order $(2,4)$, let $a$ be the fraction of rooted faces lying in
pieces with more than eight faces and $t$ the fraction of indexed edge
incidences on cycles of length at least twelve. Then
\[
  a\ \ge\ \tfrac1{12},\qquad t\ \ge\ \tfrac1{36}.
\]
Globally at least **499,200 indexed edge incidences** lie on long cycles,
in at least **3,900 distinct unordered indices**. The previous bounds were
$1/43$ and 3,266 from the [cube-diagonal count](rooted_cube_diagonals.md).
The new bound does not involve $b$, the fraction of pairs that are large
at both endpoints.

## Testing every residue

Use the 1,152 reference frames of the
[24-cell argument](indexed_eight_obstruction.md). Build the assignment
tree from a spanning tree inside **each** of the 24 vertex residues, followed
by the 23 move-0 joins listed in the verifier. For a uniform starting
frame, the tree words in a residue stay inside one rooted piece, since moves
1, 2 and 3 fix the root. Mark the residue when that piece is large. Each
tree word acts bijectively on the frames of the component, so every residue
is marked with probability $a$.

Suppose a single residue $x$ is marked. The other 23 assigned pieces are
octahedra, so every internal edge of those residues is consistent with the
tree, exactly as in the twenty-root test. For each of the 24 choices of $x$
the verifier fills all 2,304 reference edges from these edges and the tree,
using only unconditional two-color relations. This produces the frame map
respecting all four moves, which is impossible. Hence every starting frame
has at least two marked residues, and
\[
  24a\ \ge\ 2 .
\]
The rooted degree bound $L\ge f$ gives $a\le 3t$. Summing $6nt$ over the
components and rounding up to a multiple of four gives 499,200 incidences.
At most 128 belong to one index.

The same tree also fills after omitting any two residues at non-adjacent
vertices of the 24-cell (all 180 such pairs). So if a starting frame has
exactly two marks, they sit at adjacent vertices. All 96 adjacent
omissions fail, already over GF(2) (see below). Equality in the theorem
would need exactly two marks, at adjacent vertices, for every starting
frame, together with $L=f$ in every large piece.

## An upper bound in the other direction

The fraction of frames whose first two middle lines are both high is at
most $b$. If $A,B$ are both high at the root $R$, the reciprocal pieces at
$A$ (partner $B$) and at $B$ (partner $A$) are both large, by the exchange
of roles in the [cube-diagonal note](rooted_cube_diagonals.md), and that
exchange is a bijection on frames. In a large piece write $n_j$ for the
number of faces with $j$ high vertices. Then $2n_2+6n_3$ counts those
frames, and $L-f=n_2+2n_3-n_0\le(2n_2+6n_3)/2$. Summing over the pieces,
\[
  \tfrac a3\ \le\ t\ \le\ \tfrac a3+b,\qquad a\le\tfrac{1+b}2,
\]
where the last bound holds because a pair has two rooted faces and at most
a fraction $b$ of pairs have both large. In particular $t\le(1+7b)/6$, and
$t=a/3\le1/6$ when $b=0$.

Adding the diagonal events does not improve $1/36$. The point $a=1/12$,
$b=t=1/36$, $w=0$ satisfies every inequality in this note and in
[rooted_cube_diagonals.md](rooted_cube_diagonals.md), including
$w\ge a-2t-b$, $20a\ge1+19b+17w$ and $43t\ge1+2b$.

## Why the twenty-residue test stops at 43

Work over GF(2) in the cycle space of the frame graph, of dimension
$2304-1152+1=1153$. A successful filling makes every fundamental cycle of
the tree a sum of filled relations. The unconditional two-color faces span
a subspace of dimension 1,056. The 144 color-$\{1,2\}$ octagons, six in each
vertex residue, bring the span to the full 1,153; the 97 tests of the
original argument are this difference. Leaving out the residues in a set
$\Omega$ removes their octagons. Exactly **72** omitted sets of size four keep
full rank, and **no** omitted set of size five does. Rank only drops when
more residues are omitted. So no 19 vertex residues can force the frame map,
whatever tree is used. The 72 good sets all consist of two orthogonal pairs
in two different 16-cells, with the four cross pairs at inner product $-1$.

At a small root, two partners are a cube diagonal apart exactly when the
corresponding 24-cell neighbors lie in the same one of the two other
16-cells. A join tree all of whose centers keep their neighbors in one
half-cube therefore stays inside two 16-cells. Any 20 tested residues meet
all three 16-cells, so some center has neighbors in both halves and the
diagonal forests have at most $r-3$ edges. The resulting bound
$(3r-D)\,t\ge1$ then has $3r-D\ge2r+3\ge43$. The value $1/43$ is thus
the best this form of the argument can give; counting all 24 residues at
once removes the obstruction.

Run `python thin/all_residue_marks.py` from the repository root.
The [verifier](all_residue_marks.py) uses the Python standard library.
Whether a grid-free octagon can meet these bounds is still open.
