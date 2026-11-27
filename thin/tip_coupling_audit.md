# Rooting opposite pairs at their unused tips

This audits the November 2026 draft *Structural Couplings, Bipartite Cycle
Factorizations, and Spectral Obstructions* against branch revision 393f255.
Write $J$ for the cubic opposite-pair adjacency matrix and use the
[three-arm notation](opposite_seed.md). A tip is the unused point on a middle line.

## A simple graph at every point, including in the Ree–Tits octagon

Fix a point $x$. Choosing its tip line and the extensions on its two other
points gives $5(4\cdot2\cdot4)^2=5120$ unordered opposite pairs.
Each end line has point profile $(2,3,3)$ from $x$. The two other middle
lines have profile $(3,4,4)$ and each determines an outer-shell edge.

The map to unordered pairs of outer edges is **always injective**.
The two edges recover their supporting lines $C,D$. Each possible
preimage consists of two of the three middle lines of $C,D$. Two
preimages would therefore share an endpoint $A$. The unique incidence
geodesic from $x$ to $A$ has length five, and its first line is the
tip line in both constructions. Call it $F$. The three middle lines
of $C,D$ attach at distinct points on $C$, and also at distinct
points on $F$ if the two preimages have the same triple $C,D,F$:
indeed the three endpoints in their union are pairwise opposite
and each is at distance two from $F$. Their attachment points on $F$
must therefore be distinct. Exactly one endpoint pair leaves $x$
unused on $F$. This proves injectivity in the remaining case too.

For a fixed $F=\{x,y,z\}$, the far shell from $F$ has three projection
classes of 512 points. Each vertex has one neighbor in its own class and
four in each other class. Its unique inward line toward $F$ supplies the
same-class neighbor. Each of the other four lines is opposite $F$ and
has one point projecting to each of $x,y,z$. The $y,z$ cross edges thus
give **2048** supporting lines. Each supports one edge of the matching
construction in the following sense. If $C$ is a cross line, it is
opposite $F$. The arms from $F$ to $C$ starting at $y,z$ uniquely
determine their middle lines $A,B$. The pair $A,B$ has $F,C$ and one
other middle line $D$; the arms show $D$ is also a cross line, since its
attachments on the $A,B$ sides are at distance three from $y,z$
respectively. Repeating the construction at $D$ returns $C$.
This gives **one matching of 1024 pairs**. Swapping $y,z$ reverses the
names $A,B$ and leaves their unordered pair unchanged.

An outer line with boundary point $b$ participates for exactly four of
the five choices of $F$ through $x$. Exclude the first line toward $b$.
For each remaining $F$ the path from $F$ through $x,b$ to the outer line
is geodesic of incidence length eight and the other two points project
to the other two points of $F$. Consequently the five matchings form
a simple 4-regular graph on the **2560** outer edges with 5120 links.
All these counts hold without a grid hypothesis.

## The tip-incidence spectrum

Let $P$ indicate the three tips of each opposite pair. For points $x,y$
at distance three let $N$ be the middle line of their geodesic. Choose
a line at either tip other than its first geodesic line, four choices
at each end. The chosen lines $C,D$ are opposite, with one arm through
$x,N,y$. Take the middle lines $A,B$ of their other two arms. Then
$C,D$ are middles of $A,B$, and the unused points on $C,D$ are $x,y$
since the three arms have distinct attachments. Conversely these choices
are recovered from any such pair by its two tip lines. Other distinct
tip distances cannot occur. Hence $PP^T=5120I+16A_3$. Evaluation of
the point distance polynomial gives eigenvalues 15360, 5120, 4992,
5376, 3840 with multiplicities 1, 351, 650, 675, 78 respectively.
The model check enumerates all 5120 pairs at one root and checks every
off-diagonal entry there. Run `python thin/check_tip_gram.py`.

## Audit of the other claims

- **Half-shell balance.** The matching of 256 internal edges per half
  and the $128k$ vertices per half in a component of size $256k$ are
  correct. Count the four cross neighbors in each component on both
  sides. Trace bits describe which half contains a vertex. Each component
  covers one parity class of the five-bit Clebsch model, giving $16k$
  points per fiber there. A common global labeling by even words needs
  the components to use the same parity or separate relabelings.
  [Component cover](../octagon/regular256/PROOF.md)
- **Tripartite slices.** Opposition on the 96 neighbors at distance two
  gives $K_{32,32,32}$ by the eight-edge paths through the root line.
  The degree-two slices sum to $2K_{32,32}$ and repeated triples detect
  grids. These are the existing two-coordinate margins. The 128-vertex
  indexed system is the **subdivision** of the 64-vertex slice, one
  additional vertex per edge. Thus a slice cycle of length four already
  gives indexed length eight. Long-cycle fraction bounds are established
  per cubic component and globally; localization at each root needs
  another proof. [Root arrays](opposite_seed.md) and [cycles](opposite_cycles.md)
- **Exits.** At a slice vertex the two opposite partners $X,X'$ and
  the fixed root line are the three middle lines of the same opposite
  pair. Thus $X,X'$ are opposite and use different exits at both retained
  lines. Both exit classes are perfect matchings of size 32. Merged
  vertices cannot occur. Binary flips give ordinary even cycle length.
  The proposed distance-three labels and twin-branch forcing do not
  follow. The twin-branch premise itself is excluded by the four-step
  ancestor test. [Indexed paths](indexed_cycle_paths.md)
- **Middle-line matrices.** With $M$ the repository's $R$, the claimed
  Gram spectrum, trace 27,657,676,800 and grid-free one-middle degree
  9201 follow from $MM^T=3072I+3A_4$. Of the other 3071 pairs using a
  fixed middle, four share a second middle with the starting pair.
  The kernel identity holds for the sum of the operators. The endpoint
  compression table, its two unconditional min–max bounds and joint
  dimension 5849 check out. The four nonconstant block Gram determinants
  are 4951040, 5999616, 5737472 and 6217920. The omitted-line distance
  has already been proved four on failed edges, giving $q=1$.
  Eigenspace-averaged diagonal blocks do not supply a min–max bound for
  every vector of the joint space. [Boundary identities](opposite_boundary.md)
- **A replacement for the joint interlacing proof.** Every grid-free
  cubic component has girth at least eight. In the whole graph greedily
  choose vertices at mutual distance at least eight, deleting at most
  $1+3(2^7-1)=382$ vertices per choice. The induced radius-three balls
  are trees with no edges between balls. Their radial adjacency on
  normalized layers is tridiagonal with off-diagonals
  $\sqrt3,\sqrt2,\sqrt2$ and characteristic polynomial
  $t^4-7t^2+6$. Principal-submatrix interlacing gives **at least 7841
  eigenvalues on each side** of $\pm\sqrt6$ (at least $\sqrt6$ and at
  most $-\sqrt6$ respectively). This proves all four numerical tail
  claims by girth alone. It supplies no additional octagon-specific
  obstruction.
- **Component integrality and binary rank.** $DJ=2M$ gives
  $2M1_C=3D1_C$ and hence the stated endpoint and middle multiples.
  The alternating binary rank is even and at most 2,992,276 by $DJ=0$.
  The proposed deduction of a further drop by two fails: $J1=1$
  over the binary field, each column of $J$ has odd weight, and a
  four-cycle has even weight. The cited nullity bound is not a rank
  lower bound. A lower bound of 2924 is nevertheless easy: greedily
  select edges and discard all edges incident with either endpoint or
  a neighbor of an endpoint (at most 18 edges per selection). The
  induced matching gives a nonsingular principal submatrix of twice
  its size. [Existing binary bounds](opposite_boundary.md)
- **Automorphisms.** The fixed point for a 2-group is valid (also by
  odd point count). Commuting involutions have centers at point distance
  at most two. Neighbors on three different lines through one point
  have these distances and span a tree. Thus the argument for collinear
  centers, rank at most three and recognition at order 16 fails.
  A weaker bound is rank at most five. The union of the unique incidence
  geodesics between centers is a tree of diameter at most four:
  rooting at one center, every triangle of paths has total length at
  most twelve, hence is a tree. The union of the root paths is a tree
  by uniqueness below distance eight, and contains every pair path.
  Its midpoint puts the centers on a line or in a point ball of radius
  one in the point graph, at most eleven centers. For each center the
  identity and its involutions form the subgroup fixing its radius-four
  incidence ball pointwise (all nonidentity fixed balls have equal size
  and unique center). It acts freely on the four outward lines at a
  boundary point of that ball, so has size at most four. Therefore
  $|E|\le1+11\cdot3=34$, hence $|E|\le32$.
  [Fixed balls](central_involution.md)
- **Multiple roots.** The norm, dimension bound 351, independence for
  adjacent roots and inner product $|C\cap C'|-32kk'$ are correct.
  The shared half-shell contains exactly half of each squared norm.
  Cauchy–Schwarz restricted to it gives the non-strict bound
  $32\sqrt{kk'(4-k)(4-k')}$. At $k=k'=2$ this allows overlaps from
  0 through 256 in even steps. Independence on the full space does
  not imply strictness on the shared half-shell. Also a strict bound
  with irrational right side cannot in general be tightened by
  subtracting two. A subsequent [boundary count](../octagon/collinear_components.md) gives exact overlap 128.
